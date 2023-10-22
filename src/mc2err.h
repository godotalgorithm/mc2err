// mc2err: error bar standard for Markov chain data
// C99 standard compliant (C89 + permissive variable declaration + '//' comment delimiter)
// MIT license (see LICENSE file)
#ifndef MC2ERR_H
#define MC2ERR_H

// structure and function prototypes for the main C API for the mc2err library:

// mc2err data accumulator, the primary data structure of mc2err
struct mc2err_data;

// Initialize the mc2err data structure 'm2e' for data points of dimension 'width'
// and a statistical analysis that coarse grains data in pairs every 'length' number of data points.
int mc2err_initialize(struct mc2err_data *m2e, int width, int length);

// Deallocate the memory of the mc2err data structure 'm2e'.
int mc2err_deallocate(struct mc2err_data *m2e);

// Input the data point 'data' from the Markov chain with index 'chain_index' into the mc2err data structure 'm2e'.
// (The first Markov chain index is 0.)
int mc2err_input(struct mc2err_data *m2e, int chain_index, double *data);

// Output the estimated mean vector 'mean', the Markov chain covariance matrix 'covariance', and the
// uncorrelated covariance matrix 'covariance0' for a false-positive error rate less than or equal to 'error_bound'
// by analyzing all of the data contained in the mc2err data structure 'm2e'. (Both matrices are symmetric with matrix
// elements stored contiguously in memory, corresponding to both row-major and column-major dense matrix formats.)
int mc2err_output(struct mc2err_data *m2e, double error_bound, double *mean, double *covariance, double *covariance0);

// Save the mc2err data structure 'm2e' to disk in a file named 'filename'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_save(struct mc2err_data *m2e, char *filename);

// Load the mc2err data structure 'm2e' from disk in a file named 'filename'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_load(struct mc2err_data *m2e, char *filename);

// Append the mc2err data structure, 'm2e_source', to another mc2err structure, 'm2e'.
// (Chain indices from 'm2e_source' are shifted by the number of Markov chains in 'm2e'.)
int mc2err_append(struct mc2err_data *m2e, struct mc2err_data *m2e_source);

// returned error codes:
//  0 = successful return
//  1 = invalid function argument
//  2 = no data to analyze
//  3 = width or length mismatch between source and target
//  4 = file I/O error
//  5 = memory allocation failure (malloc or realloc)
//  6 = LAPACK error
//  7 = integer overflow (INT_MAX or LONG_MAX)

// OLD returned error codes:
//  0 = successful return
//  1 = no data to analyze
//  2 = width or length mismatch between source and target
//  3 = file I/O error
//  4 = memory allocation failure (malloc or realloc)
//  5 = LAPACK error (see 'lapack_info' in 'm2e' for more details)
//  6 = integer overflow (INT_MAX or LONG_MAX)

#endif
