// mc2err header file
#include "mc2err.h"

// expand memory footprint of mc2err data structure to accommodate index maxima
uint8_t mc2err_expand(uint64_t max_chain, uint8_t max_block, struct mc2err_data *m2e)
{
    // expand local memory workspaces for new Markov chains
    if(max_chain > m2e->num_chain)
    {
        m2e->chain_block = realloc(m2e->chain_block, sizeof(uint8_t)*max_chain);
        if(m2e->chain_block == NULL) { return 1; }
        m2e->chain_length = realloc(m2e->chain_length, sizeof(uint64_t)*max_chain);
        if(m2e->chain_length == NULL) { return 1; }
        m2e->chain_buffer = realloc(m2e->chain_buffer, sizeof(double*)*max_chain);
        if(m2e->chain_buffer == NULL) { return 1; }

        for(uint64_t i=m2e->num_chain ; i<max_chain ; i++)
        {
            m2e->chain_block[i] = 0;
            m2e->chain_length[i] = 0;
            m2e->chain_buffer[i] = NULL;
        }

        m2e->num_chain = max_chain;
    }

    // expand local memory workspaces for new local power-of-2 blocks
    if(max_chain > 0 && max_block > m2e->chain_block[max_chain-1])
    {
        m2e->chain_buffer[max_chain-1] = realloc(m2e->chain_buffer[max_chain-1], 
            sizeof(double)*2*max_block*m2e->width);
        if(m2e->chain_buffer[max_chain-1] == NULL) { return 1; }

        m2e->chain_block[max_chain-1] = max_block;
    }

    // expand global memory workspace for new power-of-2 blocks
    if(max_block > m2e->num_block)
    {
        // local loop size variables
        uint8_t num_block = m2e->num_block;
        uint16_t num_block2 = m2e->num_block*m2e->num_block;
        uint16_t num_block2t = m2e->num_block*(m2e->num_block+1)/2;
        uint16_t max_block2 = max_block*max_block;
        uint16_t max_block2t = max_block*(max_block+1)/2;
        uint16_t width = m2e->width;
        uint32_t width2t = m2e->width*(m2e->width+1)/2;

        // expand outer pointers for raw data
        m2e->data_num = realloc(m2e->data_num, sizeof(uint64_t)*max_block2t);
        if(m2e->data_num == NULL) { return 1; }
        m2e->data_sum1 = realloc(m2e->data_sum1, sizeof(double*)*max_block2t);
        if(m2e->data_sum1 == NULL) { return 1; }
        m2e->data_sum2 = realloc(m2e->data_sum2, sizeof(double*)*max_block2t);
        if(m2e->data_sum2 == NULL) { return 1; }
        m2e->lag_num = realloc(m2e->lag_num, sizeof(uint64_t)*max_block2t);
        if(m2e->lag_num == NULL) { return 1; }
        m2e->lag_sum2 = realloc(m2e->lag_sum2, sizeof(double*)*max_block2t);
        if(m2e->lag_sum2 == NULL) { return 1; }
        m2e->haar_num = realloc(m2e->haar_num, sizeof(uint64_t)*max_block2t);
        if(m2e->haar_num == NULL) { return 1; }
        m2e->haar_sum1 = realloc(m2e->haar_sum1, sizeof(double*)*max_block2t);
        if(m2e->haar_sum1 == NULL) { return 1; }
        m2e->haar_sum2 = realloc(m2e->haar_sum2, sizeof(double*)*max_block2t);
        if(m2e->haar_sum2 == NULL) { return 1; }

        // expand outer pointers for analysis buffers
        m2e->likelihood = realloc(m2e->likelihood, sizeof(double)*max_block2);
        if(m2e->likelihood == NULL) { return 1; }
        m2e->penalty = realloc(m2e->penalty, sizeof(double)*max_block2);
        if(m2e->penalty == NULL) { return 1; }
        m2e->mean = realloc(m2e->mean, sizeof(double*)*max_block);
        if(m2e->mean == NULL) { return 1; }
        m2e->covariance = realloc(m2e->covariance, sizeof(double*)*max_block2);
        if(m2e->covariance == NULL) { return 1; }

        // allocate new inner pointers & initialize new accumulators for raw data
        for(uint16_t i=num_block2t; i<max_block2t ; i++)
        {
            m2e->data_sum1[i] = malloc(sizeof(double)*width);
            if(m2e->data_sum1[i] == NULL) { return 1; }
            m2e->data_sum2[i] = malloc(sizeof(double)*width2t);
            if(m2e->data_sum2[i] == NULL) { return 1; }
            m2e->lag_sum2[i] = malloc(sizeof(double)*width2t);
            if(m2e->lag_sum2[i] == NULL) { return 1; }
            m2e->haar_sum1[i] = malloc(sizeof(double)*width);
            if(m2e->haar_sum1[i] == NULL) { return 1; }
            m2e->haar_sum2[i] = malloc(sizeof(double)*width2t);
            if(m2e->haar_sum2[i] == NULL) { return 1; }

            for(uint16_t j=0 ; j<width ; j++)
            { 
                m2e->data_sum1[i][j] = 0.0;
                m2e->haar_sum1[i][j] = 0.0;
            }
            for(uint32_t j=0 ; j<width2t ; j++)
            {
                m2e->data_sum2[i][j] = 0.0;
                m2e->lag_sum2[i][j] = 0.0;
                m2e->haar_sum2[i][j] = 0.0;
            }
        }

        // allocate new inner pointers for analysis buffers
        for(uint8_t i=num_block ; i<max_block ; i++)
        {
            m2e->mean[i] = malloc(sizeof(double)*width);
            if(m2e->mean[i] == NULL) { return 1; }
        }
        for(uint16_t i=num_block2 ; i<max_block2 ; i++)
        {
            m2e->covariance[i] = malloc(sizeof(double)*width2t);
            if(m2e->covariance[i] == NULL) { return 1; }
        }

        // save larger num_block value & reset analysis
        m2e->num_block = max_block;
        m2e->num_output = 0;
    }

    return 0;
}
