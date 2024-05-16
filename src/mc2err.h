// mc2err: data accumulator for Markov chains
// C99 standard compliant
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
    int width; // number of observables for which data is being gathered
    int length; // number of observable vectors retained at each level of coarse graining
    int num_level; // number of coarse-graining levels in statistical analysis
    double eqp_error; // target upper bound on false-positive error rate for equilibration point (EQP)
    double acc_error; // target upper bound on false-positive error rate for autocorrelation cutoff (ACC)

    // main outputs of the statistical analysis
    long *count; // total number of data points in the expected value of each observable
    double *mean; // width-dimensional vector of sample means
    double *variance; // width-by-width covariance matrix of the sample mean in row-major format
    double *variance0; // width-by-width covariance matrix of the observables in row-major format

    // workspace for the statistical analysis of EQP & ACC decisions
    int eqp_level; // coarse-graining level of EQP
    int acc_level; // coarse-graining level of ACC
    int eqp_index; // position of EQP
    int acc_index; // position of ACC
    double *eqp_p; // num_level-by-(2*length) matrix of P values for EQP hypothesis tests in row-major format
    double *acc_p; // num_level-by-(2*length) matrix of P values for ACC hypothesis tests in row-major format
};

// Begin the sampling process by initializing the new data accumulator 'data' for
// observable vectors of dimension 'width' and for accumulation buffers of size 'length'.
int mc2err_begin(struct mc2err_data *data, int width, int length);

// End the sampling process and deallocate the memory of the data accumulator 'data'.
int mc2err_end(struct mc2err_data *data);

// Input the observable vector 'observable' from the Markov chain with index 'chain' into the data
// accumulator 'data'. Any missing elements of the observable vector should be recorded as NaN, and
// a completely empty observable vector can be input as a NULL pointer.
int mc2err_input(struct mc2err_data *data, int chain, double *observable);

// Output the statistical analysis of the data accumulator 'data' to the analysis results 'analysis'
// for a false-positive error rate less than or equal to 'eqp_error' for the equilibration point decision
// and a false-positive error rate less than or equal to 'acc_error' for the autocorrelation cutoff decision.
int mc2err_output(struct mc2err_data *data, struct mc2err_analysis *analysis, double eqp_error, double acc_error);

// Clear and deallocate the memory of the analysis results 'analysis' after it is no longer needed
// or before it is reused in another call to 'mc2err_output'.
int mc2err_clear(struct mc2err_analysis *analysis);

// Save the data accumulator 'data' to the file on disk named 'file' in a non-portable binary format.
int mc2err_save(struct mc2err_data *data, char *file);

// Load the data accumulator 'data' from the file on disk named 'file' in a non-portable binary format.
int mc2err_load(struct mc2err_data *data, char *file);

// Map the data accumulator 'source' to form the new data accumulator 'data' for observable vectors of
// dimension 'width' and the smaller or equal buffer size 'length'. The vector 'index' of dimension
// 'width' contains the indices of the observable vectors from 'source' that are kept in 'data', and
// any out-of-bounds indices correspond to new observables with no previously recorded data.
int mc2err_map(struct mc2err_data *data, struct mc2err_data *source, int width, int length, int *index);

// Append all data from the data accumulator 'source' to the data accumulator 'data'.
// The chain indices from 'source' are offset by the number of Markov chains already in 'data'.
int mc2err_append(struct mc2err_data *data, struct mc2err_data *source);

// returned error codes:
//  0 = successful return
//  1 = invalid function argument
//  2 = invalid data point (+/- infinity)
//  3 = size mismatch between data structures
//  4 = file I/O error
//  5 = memory allocation failure (malloc or realloc)
//  6 = LAPACK error
//  7 = integer overflow (INT_MAX, LONG_MAX, or LLONG_MAX)

#endif
