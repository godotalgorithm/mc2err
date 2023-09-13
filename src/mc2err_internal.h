// mc2err internal header file
#ifndef MC2ERR_INTERNAL_H
#define MC2ERR_INTERNAL_H

// standard C headers
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// external function prototypes for BLAS & LAPACK
// NOTE: switch to dsyevr for better performance when its non-orthogonal eigenvector bug is fixed
#define MC2ERR_LAPACK_DSYEV dsyev_
void MC2ERR_LAPACK_DSYEV(char*, char*, int*, double*, int*, double*, double*, int*, int*);
#define MC2ERR_BLAS_DGEMV dgemv_
void MC2ERR_BLAS_DGEMV(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

// acronyms used in comments:
//  ACF = autocorrelation function
//  EQP = equilibration point

// pointer comment format:
//  square brackets denote the memory footprint for each pointer
//  for multiple pointers to arrays of non-uniform size, numbers refer to which index value (0-based) is used for sizing
//  e.g. "***data; // [length][length2[0]][length3[0][1]]" denotes that data[i] is defined for i in [0,length-1],
//       data[i][j] is defined for j in [0,length2[i]-1], & data[i][j][k] is defined for k in [0,length3[i][j]-1]

// mc2err data accumulator
struct mc2err_data
{
    // fixed parameters (cannot change after creation, must be equal to merge mc2err_data structures)
    int width; // vector dimension of each data point
    int length; // number of data points retained at each level of coarse graining

    // global parameters
    int num_level; // number of active coarse-graining levels
    int num_chain; // number of Markov chains
    long int num_data; // total number of data points accumulated

    // local data for each Markov chain
    int *chain_level; // number of active coarse-graining levels in each chain [num_chain]
    long int *chain_count; // number of data points accumulated in each chain [num_chain]
    double **chain_sum; // buffer of partial sums for block data accumulation [num_chain][2*chain_level[0]*length*width]

    // global data for each ACF & EQP choice
    long int *data_count; // number of data points in summation [(num_level+1)*length]
    double *data_sum; // partial sums of data points [(num_level+1)*length*width]
    long int *pair_count; // number of lagged data point pairs in summation [(num_level+1)^2*length^2]
    double *pair_sum; // partial sums of lagged data point pairs [(num_level+1)^2*length^2*width*(width+1)/2]
    double *pair_tail; // tail sums of lagged data point pairs [num_level*(num_level+1)*length*width*(width+1)/2]

    // statistical analysis results (from the most recent call to mc2err_output)
    int eqp_cut; // index of the first retained interval of equilibrated data in the mean vector
    int acf_cut; // index of the last retained interval of the autocorrelation function in the covariance matrix
    double *eqp_p_value; // P-values for the coarse-grained EQP cutoff decisions [2*num_level*length]
    double *acf_p_value; // P-values for the coarse-grained ACF cutoff decisions [2*num_level*length]

    // other information
    int lapack_info; // 'info' output of the last call to a LAPACK function (LAPACK error codes)
};

#endif
