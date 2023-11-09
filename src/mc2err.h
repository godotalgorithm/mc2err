// mc2err: error bar standard for Markov chain data
// C99 standard compliant (C89 + permissive variable declaration + '//' comment delimiter)
// MIT license (see LICENSE file)
#ifndef MC2ERR_H
#define MC2ERR_H

// structure and function prototypes for the main C API of the mc2err library:

// mc2err data accumulator
struct mc2err_data;

// mc2err analysis results
struct mc2err_analysis
{
    // analysis parameters
    int width; // vector dimension of each data point
    int length; // number of data points retained at each level of coarse graining
    int num_level; // number of coarse-graining levels in the statistical analysis
    double error; // target upper bound on false-positive error rate

    // main outputs of the statistical analysis
    long int num_data; // total number of data points used for sampled observables
    double *mean; // sample mean vector of dimension width
    double *variance; // width-by-width covariance matrix of the sample mean in row-major format
    double *variance0; // width-by-width covariance matrix of the observables in row-major format

    // intermediate quantities from the statistical analysis
    int eqp_cut; // index in eqp_p of the boundary between accepted and rejected hypothesis tests
    int acf_cut; // index in acf_p of the boundary between accepted and rejected hypothesis tests
    double *eqp_p; // num_level-by-2*length matrix of equilibration-point P-values in row-major format
    double *acf_p; // num_level-by-2*length matrix of autocorrelation-cut P-values in row-major format
};

// Initialize the data accumulator 'data' for data points of dimension 'width' and buffer size 'length'.
int mc2err_initialize(struct mc2err_data *data, int width, int length);

// Merge the array of data accumulators, 'source', of size 'num' to form the new data accumulator 'data'.
// (The order of chains in 'source' is preserved in 'data', and chain indices are offset accordingly.)
int mc2err_merge(struct mc2err_data *data, struct mc2err_data *source, int num);

// Input the data point 'point' from the Markov chain with 0-based index 'chain' into the data accumulator 'data'.
int mc2err_input(struct mc2err_data *data, int chain, double *point);

// Output the statistical analysis of the data accumulator 'data' to the analysis results 'analysis'
// for a false-positive error rate less than or equal to 'error'.
int mc2err_output(struct mc2err_data *data, struct mc2err_analysis *analysis, double error);

// Save the data accumulator 'data' to the file on disk named 'file'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_save(struct mc2err_data *data, char *file);

// Load the data accumulator 'data' from the file on disk named 'file'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_load(struct mc2err_data *data, char *file);

// Trim the data vectors in the data accumulator 'source' to form the new data accumulator 'data' with data
// points of dimension 'width' and buffer size 'length' that is less than or equal to the buffer size of 'source'.
// The array 'map' of size 'width' contains the 0-based indices of the data vector elements that are kept in 'data'.
int mc2err_trim(struct mc2err_data *data, struct mc2err_data *source, int width, int length, int *map);

// Delete the memory of the data accumulator 'data' and/or the analysis results 'analysis'.
// (Use a NULL pointer for one of the data structures to delete only an instance of the other data structure.)
int mc2err_delete(struct mc2err_data *data, struct mc2err_analysis *analysis);

// returned error codes:
//  0 = successful return
//  1 = invalid function argument
//  2 = size mismatch between data structures
//  3 = file I/O error
//  4 = memory allocation failure (malloc or realloc)
//  5 = LAPACK error
//  6 = integer overflow (INT_MAX or LONG_MAX)

#endif
