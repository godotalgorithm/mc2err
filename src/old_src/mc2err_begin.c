// mc2err header file
#include "mc2err.h"

// initialize parameters & pointers for an mc2err data structure
void mc2err_begin(uint16_t width, struct mc2err_data *m2e)
{
    // initialize all parameters
    m2e->width = width;
    m2e->num_block = 0;
    m2e->num_chain = 0;
    m2e->num_input = 0;
    m2e->num_output = 0;
    m2e->eqp_cut = 0;
    m2e->acf_cut = 0;

    // initialize all pointers to NULL
    m2e->chain_block = NULL;
    m2e->chain_length = NULL;
    m2e->chain_buffer = NULL;
    m2e->data_num = NULL;
    m2e->data_sum1 = NULL;
    m2e->data_sum2 = NULL;
    m2e->lag_num = NULL;
    m2e->lag_sum2 = NULL;
    m2e->haar_num = NULL;
    m2e->haar_sum1 = NULL;
    m2e->haar_sum2 = NULL;
    m2e->likelihood = NULL;
    m2e->penalty = NULL;
    m2e->mean = NULL;
    m2e->covariance = NULL;
}
