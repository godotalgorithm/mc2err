// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Delete the memory of the data accumulator 'data' and/or the analysis results 'analysis'.
// (Use a NULL pointer for one of the data structures to delete only an instance of the other data structure.)
int mc2err_delete(struct mc2err_data *data, struct mc2err_analysis *analysis)
{
    // free memory in data if the pointer is not NULL
    if(data != NULL)
    {
        // free inner pointers of the double pointer
        for(int i=0 ; i<data->num_chain ; i++)
        { free(data->chain_sum[i]); }

        // free all remaining pointers
        free(data->chain_level);
        free(data->chain_count);
        free(data->chain_sum);
        free(data->data_count);
        free(data->data_sum);
        free(data->pair_count);
        free(data->pair_sum);
    }

    // free memory in analysis if the pointer is not NULL
    if(analysis != NULL)
    {
        free(analysis->mean);
        free(analysis->variance);
        free(analysis->variance0);
        free(analysis->eqp_p);
        free(analysis->acf_p);
    }

    // return without errors
    return 0;
}
