// mc2err library header file
#ifndef MC2ERR_HEADER
#define MC2ERR_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

struct mc2err
{
    // basic parameters
    uint8_t data_dim; // vector dimension of data points
    uint64_t data_max; // maximum number of allowed data points (cannot exceed the maximum representable integer of uint64_t)
    double acf_precision; // relative precision for truncating the autocorrelation function
    double eqp_precision; // relative precision for assigning the equilibration point

    // derived parameters
    uint8_t acf_pow; // number of power-of-2 grids for the autocorrelation function
    uint8_t eqp_pow; // number of power-of-2 grids for the equilibration point
    uint32_t acf_num; // number of points per grid for the autocorrelation function
    uint32_t eqp_num; // number of points per grid for the equilibration point
    uint64_t data_num; // number of data points that have been collected

    // power-of-2 block data buffers
    uint64_t *buffer_offset; // block data offset for the buffer [acf_pow-1]
    uint64_t *buffer_transfer; // buffer index to insert new block data [acf_pow-1]
    uint64_t *buffer_process; // buffer index to process block data [acf_pow-1]
    double **buffer; // data buffers to store block data points [acf_pow-1][4*acf_num*data_dim]

    // partial sums of data & autocorrelation functions
    double **data_partial; // partial data sum for each equilibration point [eqp_pow][eqp_num*data_dim]
    double ****acf_partial; // partial autocorrelation sum for each equilibration point [eqp_pow][eqp_num][acf_pow][acf_num*data_dim*data_dim]

    // analysis buffers
    double **data_mean; // accumulated average for each equilibration point [eqp_pow][eqp_num*data_dim]
    double **data_error; // accumulated error bar for each equilibration point [eqp_pow][eqp_num*data_dim]
    double **acf_mean; // accumulated autocorrelation sums [acf_pow][acf_num*data_dim*data_dim]

    // FFT data structures

    // analysis results
    // false-positive error tolerance??? (maximum fraction of data to ignore in equilibration analysis)
    // error bar from each equilibrium point
    // equilibration point
    // average, error bar, & number of independent samples
};

void mc2err_initialize(uint8_t data_dim, uint64_t data_max, double acf_precision, double eqp_precision, struct mc2err *m2e);
void mc2err_finalize(struct mc2err *m2e);

uint64_t mc2err_input(uint64_t new_data_num, double *new_data, struct mc2err *m2e);
uint8_t mc2err_output(double *mean, double *error, double *size, struct mc2err *m2e);

#endif