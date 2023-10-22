// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Initialize the mc2err data structure 'm2e' for data points of dimension 'width'
// and a statistical analysis that coarse grains data in pairs every 'length' number of data points.
int mc2err_initialize(struct mc2err_data *m2e, int width, int length)
{
    // check for invalid arguments
    if(m2e == NULL || width < 1 || length < 1)
    { return 1; }

    // pass through width & length
    m2e->width = width;
    m2e->length = length;

    // initialize sizes & info to 0, except for num_level
    m2e->num_level = 1;
    m2e->num_chain = 0;
    m2e->num_data = 0;
    m2e->eqp_cut = 0;
    m2e->acf_cut = 0;

    // initialize all chain pointers to NULL
    m2e->chain_level = NULL;
    m2e->chain_count = NULL;
    m2e->chain_sum = NULL;

    // initialize all global data
    m2e->data_count = (long int*)malloc(sizeof(long int)*(m2e->num_level+1)*m2e->length);
    if(m2e->data_count == NULL) { return 5; }
    m2e->data_sum = (double*)malloc(sizeof(double)*(m2e->num_level+1)*m2e->length*m2e->width);
    if(m2e->data_sum == NULL) { return 5; }
    m2e->pair_count = (long int*)malloc(sizeof(long int)*(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length);
    if(m2e->pair_count == NULL) { return 5; }
    m2e->pair_sum = (double*)malloc(sizeof(double)*(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2);
    if(m2e->pair_sum == NULL) { return 5; }
    m2e->pair_tail = (double*)malloc(sizeof(double)*m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2);
    if(m2e->pair_tail == NULL) { return 5; }
    m2e->eqp_p_value = (double*)malloc(sizeof(double)*m2e->num_level*m2e->length);
    if(m2e->eqp_p_value == NULL) { return 5; }
    m2e->acf_p_value = (double*)malloc(sizeof(double)*m2e->num_level*m2e->length);
    if(m2e->acf_p_value == NULL) { return 5; }

    // set the data accumulators to zero
    for(int i=0 ; i<(m2e->num_level+1)*m2e->length ; i++)
    { m2e->data_count[i] = 0; }
    for(int i=0 ; i<(m2e->num_level+1)*m2e->length*m2e->width ; i++)
    { m2e->data_sum[i] = 0.0; }
    for(int i=0 ; i<(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length ; i++)
    { m2e->pair_count[i] = 0; }
    for(int i=0 ; i<(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2 ; i++)
    { m2e->pair_sum[i] = 0.0; }
    for(int i=0 ; i<m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2 ; i++)
    { m2e->pair_tail[i] = 0.0; }

    // return without errors
    return 0;
}
