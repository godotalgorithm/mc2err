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

// malloc wrapper w/ error handling
#define MC2ERR_MALLOC(PTR, TYPE, NUM) {\
    if(NUM != 0)\
    {\
        PTR = (TYPE*)malloc(sizeof(TYPE)*NUM);\
        if(PTR == NULL) { return 5; }\
    }\
    else\
    { PTR = NULL; }\
}

// realloc wrapper w/ error handling
#define MC2ERR_REALLOC(PTR, TYPE, NUM) {\
    if(NUM != 0)\
    {\
        PTR = (TYPE*)realloc(PTR, sizeof(TYPE)*NUM);\
        if(PTR == NULL) { return 5; }\
    }\
    else\
    {\
        free(PTR);\
        PTR = NULL;\
    }\
}

// free wrapper w/ NULL pointer hygiene
#define MC2ERR_FREE(PTR) {\
    free(PTR);\
    PTR = NULL;\
}

// fill wrapper
#define MC2ERR_FILL(PTR, TYPE, NUM, VALUE) {\
    TYPE* _ptr = PTR;\
    size_t _num = NUM;\
    TYPE _value = VALUE;\
    for(size_t _i=0 ; _i<_num ; _i++)\
    { _ptr[_i] = _value; }\
}

// pointer comment format:
//  square brackets denote the memory footprint for each pointer
//  for multiple pointers to arrays of non-uniform size, numbers refer to which index value (0-based) is used for sizing
//  e.g. "***data; // [length][length2[0]][length3[0][1]]" denotes that data[i] is defined for i in [0,length-1],
//       data[i][j] is defined for j in [0,length2[i]-1], & data[i][j][k] is defined for k in [0,length3[i][j]-1]

// mc2err data accumulator
struct mc2err_data
{
    // fixed parameters (cannot change after creation, must be equal to merge mc2err_data structures)
    int width; // number of observables for which data is being gathered
    int length; // number of observable vectors retained at each level of coarse graining

    // active parameters
    int num_chain; // number of Markov chains
    int max_level; // maximum number of coarse-graining levels
    long max_step; // maximum number of steps in a Markov chain
    long *max_count; // total number of data points accumulated for each observable [width]
    long long *max_pair; // maximum number of data pairs for each observable [width]

    // local data for each Markov chain
    int *num_level; // number of coarse-graining levels in each chain [num_chain]
    long *num_step; // number of steps in each chain [num_chain]
    long **local_count; // number of data points in each local buffer [num_chain][2*num_level[0]*length*width]
    double **local_sum; // local buffer of partial sums for each chain [num_chain][2*num_level[0]*length*width]

    // global data for each choice of equilibration point (EQP)
    long *global_count; // global number of data points [2*max_level*length*width]
    double *global_sum; // partial sums of data points [2*max_level*length*width]

    // global pair data for each choice of equilibration point (EQP) & autocorrelation cutoff (ACC)
    long long **pair_count; // global number of data pairs [2*max_level*length][2*max_level*length*width^2]
    double **pair_sum; // partial sums of data pairs [2*max_level*length][2*max_level*length*width^2]
};

#endif
