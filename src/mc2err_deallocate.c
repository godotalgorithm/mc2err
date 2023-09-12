// include details of the mc2err_data structure
#include "mc2err_internal.h"

// deallocate all memory using free with a loop for the double pointer
int mc2err_deallocate(struct mc2err_data *m2e)
{
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

    // no reportable errors are possible in this function
    return 0;
}
