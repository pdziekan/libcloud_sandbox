#!/bin/sh
#module add el6/openmpi/1.10.0-gnu-4.9.2-ib-cuda
#module add libs/boost/1.59.0
module add el6/openmpi/1.10.0-gnu-4.9.2-ib-cuda
#mpic++ -std=c++11 libcloud_mpi_test.cpp -lcloudphxx_lgrngn

USE_MPI=1 VERBOSE=1 mpic++ -DSTD_FUTURE_WORKS -std=c++11  -fopenmp -pthread -g -I/people/plgpdziekan/code/libcloudphxx/include -I/software/local/el6/COMMON/libs/boost/1.59.0/include -L/people/plgpdziekan/code/libcloudphxx/build/src/  libcloud_mpi_test.cpp -lcloudphxx_lgrngn

#-I/software/local/el6/COMMON/libs/boost/1.54.0-python2.7.5/include
#LD_LIBRARY_PATH="/people/plgpdziekan/code/libcloudphxx/build/src/" 
