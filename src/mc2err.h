// mc2err: standardized error bars for Markov-chain data streams
// C99 standard compliant (C89 + permissive variable declaration + fixed-width integer types)
// permissive open-source MIT license (see LICENSE file)
// <WIP> OpenMP ?.? compatible & thread-safe
// <WIP> optional MPI ?.? support 
#ifndef MC2ERR_H
#define MC2ERR_H

// standard C99 headers
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// external prototype for LAPACK real-symmetric eigensolver & BLAS matrix-vector product
// NOTE: switch to dsyevr for better performance when its non-orthogonal eigenvector bug is fixed
#define MC2ERR_LAPACK_DSYEV dsyev_
void MC2ERR_LAPACK_DSYEV(char*, char*, int*, double*, int*, double*, double*, int*, int*);
#define MC2ERR_BLAS_DGEMV dgemv_
void MC2ERR_BLAS_DGEMV(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

// the primary data structure of mc2err (defined below)
struct mc2err_data;

// initialize parameters & pointers for an mc2err data structure
void mc2err_begin(uint16_t width, struct mc2err_data *m2e);

// deallocate all memory used by an mc2err data structure
void mc2err_end(struct mc2err_data *m2e);

// input a new data vector from a Markov chain into an mc2err data structure
uint8_t mc2err_input(uint64_t chain_index, double *data, struct mc2err_data *m2e);

// output an estimated mean vector & covariance matrix from all data inside an mc2err data structure
uint8_t mc2err_output(double *mean, double *covariance, struct mc2err_data *m2e);

// save an mc2err data structure to disk in a portable binary format
uint8_t mc2err_save(char *filename, struct mc2err_data *m2e);

// load an mc2err data structure from a portable binary format on disk
uint8_t mc2err_load(char *filename, struct mc2err_data *m2e);

// append accumulated statistical data from the source to the target mc2err data structure
uint8_t mc2err_append(struct mc2err_data *source, struct mc2err_data *target);

// all-to-all append operation over an MPI communicator (ignores NULL pointers)
//uint8_t mc2err_append_mpi(struct mc2err_data *source, struct mc2err_data *target);

// returned error codes:
//  0 = successful return
//  1 = memory allocation failure (malloc or realloc)
//  2 = global data overflow (num_input > UINT64_MAX)
//  3 = data width mismatch between source & target
//  4 = file I/O error
//  5 = LAPACK error
//  6 = no data to analyze

// acronyms used in comments:
//  ACF = autocorrelation function
//  AIC = Akaike information criterion
//  EQP = equilibration point
//  MLE = maximum likelihood estimation

// the primary data structure of mc2err (square brackets denote the memory footprint for each pointer)
struct mc2err_data
{
    // basic parameters
    uint16_t width; // vector dimension of each data point (constant)
    uint8_t num_block; // number of power-of-2 blocks, max_i{chain_block[i]} (increased by mc2err_input)
    uint64_t num_chain; // number of Markov chains (increased by mc2err_input)
    uint64_t num_input; // number of data points accumulated (increased by mc2err_input)
    uint64_t num_output; // number of data points analyzed (increased by mc2err_output)

    // local data for each Markov chain (i == chain index)
    uint8_t *chain_block; // number of power-of-2 blocks in each Markov chain [num_chain]
    uint64_t *chain_length; // length of each Markov chains [num_chain]
    double **chain_buffer; // block data buffer for each Markov chain [num_chain][2*chain_block[i]*width]

    // global data for each block size & EQP interval
    uint64_t *data_num; // data count [num_block*(num_block+1)/2]
    double **data_sum1; // data sum [num_block*(num_block+1)/2][width]
    double **data_sum2; // data squared sum [num_block*(num_block+1)/2][width*(width+1)/2]
    uint64_t *lag_num; // lagged pair count [num_block*(num_block+1)/2]
    double **lag_sum2; // lagged pair sum [num_block*(num_block+1)/2][width*(width+1)/2]
    uint64_t *haar_num; // Haar wavelet count [num_block*(num_block+1)/2]
    double **haar_sum1; // Haar wavelet sum [num_block*(num_block+1)/2][width]
    double **haar_sum2; // Haar wavelet squared sum [num_block*(num_block+1)/2][width*(width+1)/2]

    // statistical analysis results (updated by mc2err_output)
    uint8_t acf_cut; // AIC-based ACF cutoff estimate
    uint8_t eqp_cut; // AIC-based EQP cutoff estimate
    double *likelihood; // -2*log(likelihood) for each ACF & EQP truncation decision [num_block^2]
    double *penalty; // AIC bias penalty for each ACF & EQP truncation decision [num_block^2]
    double **mean; // mean vector for each EQP truncation [num_block][width]
    double **covariance; // covariance matrix for each ACF & EQP truncation [num_block^2][width*(width+1)/2]

    // <WIP> OpenMP locks for thread safety ...
};

// NOTE: symmetric covariance matrices and upper-triangular block/EQP-indexed matrices
//       are stored using the LAPACK-style packed upper-triangle matrix format where
//       row i & column j correspond to matrix index i+j*(j+1)/2 for 0-based indexing

#endif