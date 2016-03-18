#include <libcloudph++/lgrngn/factory.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <stdio.h>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>
#include <iostream>
#if defined(USE_MPI)
  #include "mpi.h"
#endif


using namespace std;
using namespace libcloudphxx::lgrngn;
  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, double>
    mean_rd1 = double(30e-6) * si::metres,
    mean_rd2 = double(40e-6) * si::metres;
  const quantity<si::dimensionless, double>
    sdev_rd1 = double(1.4),
    sdev_rd2 = double(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, double>
    n1_stp = double(60e6) / si::cubic_metres,
    n2_stp = double(40e6) / si::cubic_metres;



// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, double>(lnrd)) +
        lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, double>(lnrd)) 
      ) * si::cubic_metres
    );  
  }   

  log_dry_radii *do_clone() const 
  { return new log_dry_radii( *this ); }
};  


void two_step(particles_proto_t<double> *prtcls, 
             arrinfo_t<double> th,
             arrinfo_t<double> rhod,
             arrinfo_t<double> rv,
             opts_t<double> opts)
{
    prtcls->step_sync(opts,th,rv,rhod);
    prtcls->step_async(opts);
}


int main(int argc, char *argv[]){
  opts_init_t<double> opts_init;

  int ndims;
  sscanf(argv[1], "%d", &ndims);
  printf("ndims %d\n", ndims);

  int rank = 0;
#if defined(USE_MPI)
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  opts_init.dt=100.;
  opts_init.sstp_coal = 1; 
  opts_init.kernel = kernel_t::geometric;
  opts_init.terminal_velocity = vt_t::beard76;
  opts_init.dx = 1;
  opts_init.nx = 1*(rank+1); 
  opts_init.x1 = 1*(rank+1);
  opts_init.sd_conc = 64;
  opts_init.n_sd_max = 2000;
  opts_init.rng_seed = 4444 + rank;
  if(ndims>1)
  {
    opts_init.dz = 1; 
    opts_init.nz = 1; 
    opts_init.z1 = 1; 
  }
//  opts_init.sd_const_multi = 1;

  boost::assign::ptr_map_insert<
    log_dry_radii<double> // value type
  >(  
    opts_init.dry_distros // map
  )(  
    0.001 // key
  ); 

  particles_proto_t<double> *prtcls;
     prtcls = factory<double>(
        (backend_t)serial, 
        opts_init
      );
  double pth[] = {300., 300., 300.};
  double prhod[] = {1., 1., 1.};
  double prv[] = {.01, 0.01, 0.01};
  double pCx[] = {1, 1, 1, 1};
  double pCz[] = {0., 0., 0.};
  //long int strides[] = {sizeof(double)};
  long int strides[] = {1, 3};
  long int xstrides[] = {1, 4};
  long int ystrides[] = {1, 1};

  arrinfo_t<double> th(pth, strides);
  arrinfo_t<double> rhod(prhod, strides);
  arrinfo_t<double> rv(prv, strides);
  arrinfo_t<double> Cx(pCx, xstrides);
  arrinfo_t<double> Cz(pCz, ystrides);

  if(ndims==1)
    prtcls->init(th,rhod,rv, Cx);
  else if(ndims==2)
    prtcls->init(th,rhod,rv, Cx, arrinfo_t<double>(), Cz);

  opts_t<double> opts;
  opts.adve = 0;
  opts.sedi = 0;
  opts.cond = 0;
  opts.coal = 1;
//  opts.chem = 0;


  prtcls->diag_all();
  prtcls->diag_sd_conc();
  double *out = prtcls->outbuf();
  printf("---sd_conc init---\n");
  printf("%d: %lf %lf %lf\n",rank, out[0], out[1], out[2]);
#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  

  for(int i=0;i<70;++i)
  {
//    if(rank==0)
      two_step(prtcls,th,rhod,rv,opts);
   //   MPI_Barrier(MPI_COMM_WORLD);
  //  if(rank==1)
   //   two_step(prtcls,th,rhod,rv,opts);
   // MPI_Barrier(MPI_COMM_WORLD);
  }
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  out = prtcls->outbuf();
  printf("---sd_conc po coal---\n");
  printf("%d: %lf %lf %lf\n",rank, out[0], out[1], out[2]);
  opts.coal = 0;
  opts.adve = 1;
  two_step(prtcls,th,rhod,rv,opts);
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  out = prtcls->outbuf();

#if defined(USE_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  printf("---sd_conc po adve---\n");
  printf("%d: %lf %lf %lf\n",rank, out[0], out[1], out[2]);
#if defined(USE_MPI)
  MPI_Finalize();
#endif
}
