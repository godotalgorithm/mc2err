// include details of the mc2err_data structure
#include "mc2err_internal.h"

// sizes are all set to 0 and pointers are all set to NULL so that realloc can be used for all memory allocation
int mc2err_initialize(struct mc2err_data *m2e, int width, int length)
{
    // pass through width & length
    m2e->width = width;
    m2e->length = length;

    // initialize sizes & info to 0
    m2e->num_chain = 0;
    m2e->num_level = 0;
    m2e->num_data = 0;
    m2e->eqp_cut = 0;
    m2e->acf_cut = 0;
    m2e->lapack_info = 0;

    // initialize all pointers to NULL
    m2e->chain_level = NULL;
    m2e->chain_count = NULL;
    m2e->chain_sum = NULL;
    m2e->data_count = NULL;
    m2e->data_sum = NULL;
    m2e->pair_count = NULL;
    m2e->pair_sum = NULL;
    m2e->pair_tail = NULL;
    m2e->eqp_p_value = NULL;
    m2e->acf_p_value = NULL;

    // no reportable errors are possible in this function
    return 0;
}
