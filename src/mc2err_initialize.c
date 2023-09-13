// include details of the mc2err_data structure
#include "mc2err_internal.h"

// sizes are all set to 0, except for num_level=1, which avoids any bookkeeping associated with num_level=0
int mc2err_initialize(struct mc2err_data *m2e, int width, int length)
{
    // pass through width & length
    m2e->width = width;
    m2e->length = length;

    // initialize sizes & info to 0, except for num_level
    m2e->num_level = 1;
    m2e->num_chain = 0;
    m2e->num_data = 0;
    m2e->eqp_cut = 0;
    m2e->acf_cut = 0;
    m2e->lapack_info = 0;

    // initialize all chain pointers to NULL
    m2e->chain_level = NULL;
    m2e->chain_count = NULL;
    m2e->chain_sum = NULL;

    // initialize all global data
    m2e->data_count = (long int*)malloc(sizeof(long int)*(m2e->num_level+1)*m2e->length);
    if(m2e->data_count == NULL) { return 4; }
    m2e->data_sum = (double*)malloc(sizeof(double)*(m2e->num_level+1)*m2e->length*m2e->width);
    if(m2e->data_sum == NULL) { return 4; }
    m2e->pair_count = (long int*)malloc(sizeof(long int)*(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length);
    if(m2e->pair_count == NULL) { return 4; }
    m2e->pair_sum = (double*)malloc(sizeof(double)*(m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2);
    if(m2e->pair_sum == NULL) { return 4; }
    m2e->pair_tail = (double*)malloc(sizeof(double)*m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2);
    if(m2e->pair_tail == NULL) { return 4; }
    m2e->eqp_p_value = (double*)malloc(sizeof(double)*m2e->num_level*m2e->length);
    if(m2e->eqp_p_value == NULL) { return 4; }
    m2e->acf_p_value = (double*)malloc(sizeof(double)*m2e->num_level*m2e->length);
    if(m2e->acf_p_value == NULL) { return 4; }

    // return without errors
    return 0;
}
