// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Delete the memory of the mc2err data structure 'm2e'.
int mc2err_delete(struct mc2err_data *m2e)
{
    // check for invalid arguments
    if(m2e == NULL)
    { return 1; }

    // free inner pointers of the double pointer
    for(int i=0 ; i<m2e->num_chain ; i++)
    { free(m2e->chain_sum[i]); }

    // free all remaining pointers
    free(m2e->chain_level);
    free(m2e->chain_count);
    free(m2e->chain_sum);
    free(m2e->data_count);
    free(m2e->data_sum);
    free(m2e->pair_count);
    free(m2e->pair_sum);
    free(m2e->pair_tail);
    free(m2e->eqp_p_value);
    free(m2e->acf_p_value);

    // return without errors
    return 0;
}
