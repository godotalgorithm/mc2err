// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: Vector addition is performed here with simple for loops, assuming that compiler optimizations
//       will unroll them if that is more efficient, and no attempt is made to use level-1 BLAS,
//       which might be more efficient on certain machines in certain size regimes.

// Append the mc2err data structure, 'm2e_source', to another mc2err structure, 'm2e'.
// (Chain indices from 'm2e_source' are shifted by the number of Markov chains in 'm2e'.)
int mc2err_append(struct mc2err_data *m2e, struct mc2err_data *m2e_source)
{
    // check for invalid arguments
    if(m2e == NULL || m2e_source == NULL)
    { return 1; }

    // check for consistency
    if(m2e->width != m2e_source->width || m2e->length != m2e_source->length)
    { return 3; }

    // check for overflow in the total number of chains
    if(INT_MAX - m2e->num_chain < m2e_source->num_chain)
    { return 7; }

    // check for overflow in the total number of accumulated data points
    if(LONG_MAX - m2e->num_data < m2e_source->num_data)
    { return 7; }

    // append local data
    m2e->chain_level = (int*)realloc(m2e->chain_level, sizeof(int)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_level == NULL) { return 5; }
    memcpy(m2e->num_level+m2e->num_chain, m2e_source->num_level, sizeof(int)*m2e_source->num_chain);

    m2e->chain_count = (long int*)realloc(m2e->chain_count, sizeof(long int)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_count == NULL) { return 5; }
    memcpy(m2e->chain_count+m2e->num_chain, m2e_source->chain_count, sizeof(long int)*m2e_source->num_chain);

    m2e->chain_sum = (double**)realloc(m2e->chain_sum, sizeof(double*)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_sum == NULL) { return 5; }
    for(int i=0 ; i<m2e_source->num_chain ; i++)
    {
        size_t chain_sum_size = sizeof(double)*2*m2e_source->chain_level[i]*m2e->length*m2e->width;
        m2e->chain_sum[m2e->num_chain+i] = (double*)malloc(chain_sum_size);
        if(m2e->chain_sum[m2e->num_chain+i] == NULL) { return 5; }
        memcpy(m2e->chain_sum[m2e->num_chain+i], m2e_source->chain_sum[i], chain_sum_size);
    }

    // update num_chain & num_data
    m2e->num_chain += m2e_source->num_chain;
    m2e->num_data += m2e_source->num_data;

    // expand global data if num_level increases
    if(m2e->num_level < m2e_source->num_level)
    {
        int status = mc2err_expand(m2e, m2e_source->num_level);
        if(status != 0) { return status; }
    }

    // append global data
    for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length ; i++)
    { m2e->data_count[i] += m2e_source->data_count[i]; }
    for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length*m2e->width ; i++)
    { m2e->data_sum[i] += m2e_source->data_sum[i]; }

    int stride1 = (m2e_source->num_level+1)*m2e->length;
    int stride2 = (m2e->num_level+1)*m2e->length;
    for(int i=0 ; i<stride1 ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_count[j+i*stride2] += m2e_source->pair_count[j+i*stride1]; }

    stride1 = (m2e_source->num_level+1)*m2e->length*(m2e->width+1)/2;
    stride2 = (m2e->num_level+1)*m2e->length*(m2e->width+1)/2;
    for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_sum[j+i*stride2] += m2e_source->pair_sum[j+i*stride1]; }

    stride1 = m2e_source->num_level*(m2e->width+1)/2;
    stride2 = m2e->num_level*(m2e->width+1)/2;
    for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_tail[j+i*stride2] += m2e_source->pair_tail[j+i*stride1]; }

    // return without errors
    return 0;
}
