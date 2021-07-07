// mc2err library: online averages, error bars, & equilibration analysis for Markov chain Monte Carlo data streams
// C99 standard compliant (C89 + permissive variable declaration + fixed-width integer types)
#ifndef MC2ERR_H
#define MC2ERR_H

// common acronyms:
//  ACF = autocorrelation function
//  EQP = equilibration point

// error codes:
//  0 = successful return
//  1 = ACF precision is invalid (<2.4e-10 or >1.0)
//  2 = EQP precision is invalid (<2.4e-10 or >1.0)
//  3 = memory allocation failure
//  4 = maximum number of data points exceeded

// mc2err uses C99 fixed-width integer types
#include <stdint.h>

// the primary data structure of mc2err (square brackets denote the maximum memory footprint for each pointer)
struct mc2err
{
    // basic parameters
    uint8_t data_dim; // vector dimension of each data point
    uint64_t data_max; // maximum number of allowed data points (cannot exceed UINT64_MAX)
    uint64_t data_num; // number of data points that have been collected
    double acf_precision; // relative precision for truncating the ACF
    double eqp_precision; // relative precision for assigning the EQP

    // derived parameters
    uint8_t acf_pow; // number of power-of-2 block grids for the ACF
    uint32_t acf_num; // number of points per block grid for the ACF (base grid has 2x points)
    uint8_t eqp_pow; // number of power-of-2 block grids for the EQP
    uint32_t eqp_num; // number of points per block grid for the EQP (base grid has 2x points)
    uint8_t data_eqp; // active EQP block size of the data

    // power-of-2 block data buffers
    uint8_t *buffer_eqp; // active EQP block size of the buffer [acf_pow]
    uint64_t *buffer_offset; // block data offset for the buffer [acf_pow]
    uint64_t *buffer_transfer; // buffer index to insert new block data [acf_pow]
    uint64_t *buffer_head; // head buffer index to process block data [acf_pow]
    uint64_t *buffer_tail; // head buffer index to process block data [acf_pow]
    double **buffer; // data buffers to store block data points [acf_pow][4*acf_num*data_dim]

    // partial sums of data & autocorrelation functions
    double **data_partial; // partial data sum for each EQP [eqp_pow][2*eqp_num*data_dim]
    double ****acf_partial; // partial ACF sum for each EQP [eqp_pow][2*eqp_num][acf_pow][2*acf_num*data_dim*data_dim]

    // analysis buffers
    double *data_sum; // accumulated data sum [data_dim]
    double **data_mean; // accumulated average for each EQP [eqp_pow][2*eqp_num*data_dim]
    double **data_error; // estimated error bar for each EQP [eqp_pow][2*eqp_num*data_dim*data_dim]
    double **acf_sum; // accumulated ACF sums [acf_pow][2*acf_num*data_dim*data_dim]
    double **acf_mean; // accumulated ACF averages [acf_pow][2*acf_num*data_dim*data_dim]
    // likelihood & ACF truncation point for each EQP

    // FFT data structures

    // analysis results
    // false-positive error tolerance??? (maximum fraction of data to ignore in equilibration analysis)
    // error bar from each equilibrium point
    // equilibration point
    // average, error bar, & number of independent samples
};

// initializes parameters & allocates memory for an mc2err structure
uint8_t mc2err_initialize(uint8_t data_dim, uint64_t data_max, double acf_precision, double eqp_precision, struct mc2err *m2e);

// input new MCMC data into the mc2err data structure
uint8_t mc2err_input(uint64_t new_data_num, double *new_data, struct mc2err *m2e);

// analyze all available data in mc2err
//uint8_t mc2err_output(double acf_bias, double eqp_bias, struct mc2err *m2e);
uint8_t mc2err_output(uint64_t acf_cut, uint64_t eqp_cut, struct mc2err *m2e);

// deallocate all of the memory used by mc2err
void mc2err_finalize(struct mc2err *m2e);

#endif