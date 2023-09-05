// mc2err header file
#include "mc2err.h"

// deallocate all memory used by an mc2err data structure
void mc2err_end(struct mc2err_data *m2e)
{
    // local loop size variables
    uint8_t num_block = m2e->num_block;
    uint16_t num_block2 = m2e->num_block*m2e->num_block;
    uint16_t num_block2t = m2e->num_block*(m2e->num_block+1)/2;
    uint64_t num_chain = m2e->num_chain;

    // loop over double-pointer indices to deallocate inner pointers
    for(uint64_t i=0 ; i<num_chain ; i++)
    { free(m2e->chain_buffer[i]); }
    for(uint8_t i=0 ; i<num_block ; i++)
    { free(m2e->mean[i]); }
    for(uint16_t i=0 ; i<num_block2t ; i++)
    {
        free(m2e->data_sum1[i]);
        free(m2e->data_sum2[i]);
        free(m2e->lag_sum2[i]);
        free(m2e->haar_sum1[i]);
        free(m2e->haar_sum2[i]);
    }
    for(uint16_t i=0 ; i<num_block2 ; i++)
    { free(m2e->covariance[i]); }

    // deallocate outer pointers
    free(m2e->chain_block);
    free(m2e->chain_length);
    free(m2e->chain_buffer);
    free(m2e->data_num);
    free(m2e->data_sum1);
    free(m2e->data_sum2);
    free(m2e->lag_num);
    free(m2e->lag_sum2);
    free(m2e->haar_num);
    free(m2e->haar_sum1);
    free(m2e->haar_sum2);
    free(m2e->likelihood);
    free(m2e->penalty);
    free(m2e->mean);
    free(m2e->covariance);
}
