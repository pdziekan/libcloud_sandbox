#!/bin/sh
module add libs/boost/1.59.0

VERBOSE=1 g++ -DSTD_FUTURE_WORKS -std=c++11  -fopenmp -pthread -g -I/people/plgpdziekan/code/libcloudphxx/include -I/software/local/el6/COMMON/libs/boost/1.59.0/include -L/people/plgpdziekan/code/libcloudphxx/build/src/ libcloud_mpi_test.cpp -lcloudphxx_lgrngn

#-I/software/local/el6/COMMON/libs/boost/1.54.0-python2.7.5/include
#LD_LIBRARY_PATH="/people/plgpdziekan/code/libcloudphxx/build/src/" 
