#!/bin/sh
mpic++ -std=c++11 libcloud_mpi_test.cpp -lcloudphxx_lgrngn

#VERBOSE=1 LD_LIBRARY_PATH="/home/piotr/praca/libcloudphxx/build/src/" mpic++ -DSTD_FUTURE_WORKS -std=c++11  -fopenmp -pthread -g -I/home/piotr/praca/libcloudphxx/include profiling_test.cpp -lcloudphxx_lgrngn 


