// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Expand the memory for global data in the mc2err data structure 'm2e' to match a given larger value of 'num_level'.
int mc2err_expand(struct mc2err_data *m2e, int num_level)
{
    // check for invalid arguments
    if(m2e == NULL || num_level <= m2e->num_level)
    { return 1; }

    // expand memory for data
    size_t new_size = sizeof(long int)*(num_level+1)*m2e->length;
    m2e->data_count = (long int*)realloc(m2e->data_count, new_size);
    if(m2e->data_count == NULL) { return 5; }

    new_size = sizeof(double)*(num_level+1)*m2e->length*m2e->width;
    m2e->data_sum = (double*)realloc(m2e->data_sum, new_size);
    if(m2e->data_sum == NULL) { return 5; }

    // pad new data with zero
    for(int i=(m2e->num_level+1)*m2e->length ; i<(num_level+1)*m2e->length ; i++)
    { m2e->data_count[i] = 0; }
    for(int i=(m2e->num_level+1)*m2e->length*m2e->width ; i<(num_level+1)*m2e->length*m2e->width ; i++)
    { m2e->data_sum[i] = 0.0; }

    // save pointers to old pair data
    long int *old_pair_count = m2e->pair_count;
    double *old_pair_sum = m2e->pair_sum;
    double *old_pair_tail = m2e->pair_tail;

    // allocate new memory for pair data
    new_size = sizeof(long int)*(num_level+1)*(num_level+1)*m2e->length*m2e->length;
    m2e->pair_count = (long int*)malloc(new_size);
    if(m2e->pair_count == NULL) { return 5; }

    new_size = sizeof(double)*(num_level+1)*(num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2;
    m2e->pair_sum = (double*)malloc(new_size);
    if(m2e->pair_sum == NULL) { return 5; }

    new_size = sizeof(double)*num_level*(num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    m2e->pair_tail = (double*)malloc(new_size);
    if(m2e->pair_tail == NULL) { return 5; }

    // copy old pair data
    int stride1 = (num_level+1)*m2e->length;
    int stride2 = num_level*m2e->length;
    for(int i=0 ; i<num_level*m2e->length ; i++)
    { memcpy(m2e->pair_count+i*stride1, old_pair_count+i*stride2, sizeof(long int)*stride2); }

    stride1 = (num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    stride2 = num_level*m2e->length*m2e->width*(m2e->width+1)/2;
    for(int i=0 ; i<num_level*m2e->length ; i++)
    { memcpy(m2e->pair_sum+i*stride1, old_pair_sum+i*stride2, sizeof(double)*stride2); }

    stride1 = (num_level+1)*m2e->width*(m2e->width+1)/2;
    stride2 = num_level*m2e->width*(m2e->width+1)/2;
    for(int i=0 ; i<num_level*m2e->length ; i++)
    { memcpy(m2e->pair_tail+i*stride1, old_pair_tail+i*stride2, sizeof(double)*stride2); }

    // pad new pair data with zero
    stride1 = (num_level+1)*m2e->length;
    stride2 = num_level*m2e->length;
    for(int i=0 ; i<stride2 ; i++)
    for(int j=stride2 ; j<stride1 ; j++)
    { m2e->pair_count[j+i*stride1] = 0.0; }
    for(int i=stride2 ; i<stride1 ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_count[j+i*stride1] = 0.0; }

    stride1 = (num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    stride2 = num_level*m2e->length*m2e->width*(m2e->width+1)/2;
    for(int i=0 ; i<num_level*m2e->length ; i++)
    for(int j=stride2 ; j<stride1 ; j++)
    { m2e->pair_sum[j+i*stride1] = 0.0; }
    for(int i=num_level*m2e->length ; i<(num_level+1)*m2e->length ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_sum[j+i*stride1] = 0.0; }

    stride1 = num_level*m2e->width*(m2e->width+1)/2;
    stride2 = num_level*m2e->width*(m2e->width+1)/2;
    for(int i=0 ; i<num_level*m2e->length ; i++)
    for(int j=stride2 ; j<stride1 ; j++)
    { m2e->pair_tail[j+i*stride1] = 0.0; }
    for(int i=m2e->num_level*m2e->length ; i<(num_level+1)*m2e->length ; i++)
    for(int j=0 ; j<stride1 ; j++)
    { m2e->pair_tail[j+i*stride1] = 0.0; }

    // deallocate memory of old pair data
    free(old_pair_count);
    free(old_pair_sum);
    free(old_pair_tail);

    // expand statistical analysis buffers
    m2e->eqp_p_value = (double*)realloc(m2e->eqp_p_value, sizeof(double)*2*num_level*m2e->length);
    if(m2e->eqp_p_value == NULL) { return 5; }
    m2e->acf_p_value = (double*)realloc(m2e->acf_p_value, sizeof(double)*2*num_level*m2e->length);
    if(m2e->acf_p_value == NULL) { return 5; }

    // update num_level
    m2e->num_level = num_level;

    // return without errors
    return 0;
}
