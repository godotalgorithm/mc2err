// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// prototype of internal helper function
uint8_t pow2grid_initialize(uint64_t max, double precision, uint8_t *pow, uint32_t *num);

// initializes parameters & allocates memory for an mc2err structure
uint8_t mc2err_initialize(uint8_t data_dim, uint64_t data_max, double acf_precision, double eqp_precision, struct mc2err *m2e)
{
    // initialize the basic parameters
    m2e->data_num = 0;
    m2e->data_dim = data_dim;
    m2e->data_max = data_max;
    m2e->data_eqp = 0;
    m2e->acf_precision = acf_precision;
    m2e->eqp_precision = eqp_precision;

    // initialize power-of-2 logarithmic grids
    if(pow2grid_initialize(m2e->data_max, m2e->acf_precision, &(m2e->acf_pow), &(m2e->acf_num)))
    { return 1; }
    if(pow2grid_initialize(m2e->data_max, m2e->eqp_precision, &(m2e->eqp_pow), &(m2e->eqp_num)))
    { return 2; }

    // allocate memory for buffer pointers & indices
    m2e->buffer_eqp = NULL; MEM_CHECK_SET(m2e->buffer_eqp, uint8_t, m2e->acf_pow, 0);
    m2e->buffer_offset = NULL; MEM_CHECK_SET(m2e->buffer_offset, uint64_t, m2e->acf_pow, 0);
    m2e->buffer_transfer = NULL; MEM_CHECK_SET(m2e->buffer_transfer, uint64_t, m2e->acf_pow, 0);
    m2e->buffer_head = NULL; MEM_CHECK_SET(m2e->buffer_head, uint64_t, m2e->acf_pow, 0);
    m2e->buffer_tail = NULL; MEM_CHECK_SET(m2e->buffer_tail, uint64_t, m2e->acf_pow, 0);
    m2e->buffer = NULL; MEM_CHECK_SET(m2e->buffer, double*, m2e->acf_pow, NULL);

    // allocate the 1st data buffer
    MEM_CHECK(m2e->buffer[0], double, 4*m2e->acf_num*m2e->data_dim);

    // allocate memory for data & acf storage
    m2e->data_sum = NULL; MEM_CHECK(m2e->data_sum, double, m2e->data_dim);
    m2e->data_mean = NULL; MEM_CHECK_SET(m2e->data_mean, double*, m2e->eqp_pow, NULL);
    m2e->data_partial = NULL; MEM_CHECK_SET(m2e->data_partial, double*, m2e->eqp_pow, NULL);
    m2e->data_error = NULL; MEM_CHECK_SET(m2e->data_error, double*, m2e->eqp_pow, NULL);
    m2e->acf_sum = NULL; MEM_CHECK_SET(m2e->acf_sum, double*, m2e->acf_pow, NULL);
    m2e->acf_mean = NULL; MEM_CHECK_SET(m2e->acf_mean, double*, m2e->acf_pow, NULL);
    m2e->acf_partial = NULL; MEM_CHECK_SET(m2e->acf_partial, double***, m2e->eqp_pow, NULL);

    // TODO: setup FFT data structures for convolution
}

// initialize a power-of-2 logarithmic grid
uint8_t pow2grid_initialize(uint64_t max, double precision, uint8_t *pow, uint32_t *num)
{
    // check if there are too few grid points
    if(precision > 1.0)
    { return 1; }

    // check if there are too many grid points
    if(precision < 2.4e-10) // overflow-free proxy of: if(2*ceil(0.5 + 0.5/precision) > UINT32_MAX)
    { return 1; }

    *num = ceil(0.5 + 0.5/precision);
    *pow = 1;
    uint64_t exp2pow = max/(2*(*num)) + (max%(2*(*num)) != 0);
    while((1<<(*pow-1)) < exp2pow)
    { (*pow)++; }

    return 0;
}