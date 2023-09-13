// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: Vector addition is performed here with simple for loops, assuming that compiler optimizations
//       will unroll them if that is more efficient, and no attempt is made to use level-1 BLAS,
//       which might be more efficient on certain machines in certain size regimes.

// allocate all memory using malloc with a loop for the double pointer
int mc2err_append(struct mc2err_data *m2e, struct mc2err_data *m2e_source)
{
    // check for consistency
    if(m2e->width != m2e_source->width || m2e->length != m2e_source->length)
    { return 2; }

    // check for overflow in the total number of chains
    if(INT_MAX - m2e->num_chain < m2e_source->num_chain)
    { return 6; }

    // check for overflow in the total number of accumulated data points
    if(LONG_MAX - m2e->num_data < m2e_source->num_data)
    { return 6; }

    // append local data
    m2e->chain_level = (int*)realloc(m2e->chain_level, sizeof(int)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_level == NULL) { return 4; }
    memcpy(m2e->num_level+m2e->num_chain, m2e_source->num_level, sizeof(int)*m2e_source->num_chain);

    m2e->chain_count = (long int*)realloc(m2e->chain_count, sizeof(long int)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_count == NULL) { return 4; }
    memcpy(m2e->chain_count+m2e->num_chain, m2e_source->chain_count, sizeof(long int)*m2e_source->num_chain);

    m2e->chain_sum = (double**)realloc(m2e->chain_sum, sizeof(double*)*(m2e->num_chain + m2e_source->num_chain));
    if(m2e->chain_sum == NULL) { return 4; }
    for(int i=0 ; i<m2e_source->num_chain ; i++)
    {
        size_t chain_sum_size = sizeof(double)*2*m2e_source->chain_level[i]*m2e->length*m2e->width;
        m2e->chain_sum[m2e->num_chain+i] = (double*)malloc(chain_sum_size);
        if(m2e->chain_sum[m2e->num_chain+i] == NULL) { return 4; }
        memcpy(m2e->chain_sum[m2e->num_chain+i], m2e_source->chain_sum[i], chain_sum_size);
    }

    // update num_chain & num_data
    m2e->num_chain += m2e_source->num_chain;
    m2e->num_data += m2e_source->num_data;

    // append global data in place if num_level doesn't change
    if(m2e->num_level >= m2e_source->num_level)
    {
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
    }
    else // reallocate memory to append global data if num_level increases
    {
        // expand memory for data
        size_t new_size = sizeof(long int)*(m2e_source->num_level+1)*m2e->length;
        m2e->data_count = (long int*)realloc(m2e->data_count, new_size);
        if(m2e->data_count == NULL) { return 4; }

        new_size = sizeof(double)*(m2e_source->num_level+1)*m2e->length*m2e->width;
        m2e->data_sum = (double*)realloc(m2e->data_sum, new_size);
        if(m2e->data_sum == NULL) { return 4; }

        // add data from m2e_source
        for(int i=0 ; i<(m2e->num_level+1)*m2e->length ; i++)
        { m2e->data_count[i] += m2e_source->data_count[i]; }
        for(int i=(m2e->num_level+1)*m2e->length ; i<(m2e_source->num_level+1)*m2e->length ; i++)
        { m2e->data_count[i] = m2e_source->data_count[i]; }
        for(int i=0 ; i<(m2e->num_level+1)*m2e->length*m2e->width ; i++)
        { m2e->data_sum[i] += m2e_source->data_sum[i]; }
        for(int i=(m2e->num_level+1)*m2e->length*m2e->width ; i<(m2e_source->num_level+1)*m2e->length*m2e->width ; i++)
        { m2e->data_sum[i] = m2e_source->data_sum[i]; }

        // save pointers to old pair data
        long int *old_pair_count = m2e->pair_count;
        double *old_pair_sum = m2e->pair_sum;
        double *old_pair_tail = m2e->pair_tail;

        // allocate new memory & move pair data from m2e_source
        new_size = sizeof(long int)*(m2e_source->num_level+1)*(m2e_source->num_level+1)*m2e->length*m2e->length;
        m2e->pair_count = (long int*)malloc(new_size);
        if(m2e->pair_count == NULL) { return 4; }
        memcpy(m2e->pair_count, m2e_source->pair_count, new_size);

        new_size = sizeof(double)*(m2e_source->num_level+1)*(m2e_source->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2;
        m2e->pair_sum = (double*)malloc(new_size);
        if(m2e->pair_sum == NULL) { return 4; }
        memcpy(m2e->pair_sum, m2e_source->pair_sum, new_size);

        new_size = sizeof(double)*m2e_source->num_level*(m2e_source->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
        m2e->pair_tail = (double*)malloc(new_size);
        if(m2e->pair_tail == NULL) { return 4; }
        memcpy(m2e->pair_tail, m2e_source->pair_tail, new_size);

        // add pair data from m2e
        int stride1 = (m2e_source->num_level+1)*m2e->length;
        int stride2 = (m2e->num_level+1)*m2e->length;
        for(int i=0 ; i<stride1 ; i++)
        for(int j=0 ; j<stride1 ; j++)
        { m2e->pair_count[j+i*stride1] += old_pair_count[j+i*stride2]; }

        stride1 = (m2e_source->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
        stride2 = (m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
        for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length ; i++)
        for(int j=0 ; j<stride1 ; j++)
        { m2e->pair_sum[j+i*stride1] += old_pair_sum[j+i*stride2]; }

        stride1 = m2e_source->num_level*m2e->width*(m2e->width+1)/2;
        stride2 = m2e->num_level*m2e->width*(m2e->width+1)/2;
        for(int i=0 ; i<(m2e_source->num_level+1)*m2e->length ; i++)
        for(int j=0 ; j<stride1 ; j++)
        { m2e->pair_tail[j+i*stride1] += old_pair_tail[j+i*stride2]; }

        // deallocate memory of old pair data
        free(old_pair_count);
        free(old_pair_sum);
        free(old_pair_tail);

        // expand statistical analysis buffers
        m2e->eqp_p_value = (double*)realloc(m2e->eqp_p_value, sizeof(double)*2*m2e_source->num_level*m2e->length);
        if(m2e->eqp_p_value == NULL) { return 4; }
        m2e->acf_p_value = (double*)realloc(m2e->acf_p_value, sizeof(double)*2*m2e_source->num_level*m2e->length);
        if(m2e->acf_p_value == NULL) { return 4; }

        // update num_level
        m2e->num_level = m2e_source->num_level;
    }

    // return without errors
    return 0;
}
