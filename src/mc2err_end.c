// include details of the mc2err_data structure
#include "mc2err_internal.h"

// End the sampling process and deallocate the memory of the data accumulator 'data'.
int mc2err_end(struct mc2err_data *data)
{
    // check for invalid arguments
    if(data == NULL)
    { return 1; }

    // free inner pointers of the double pointers
    for(int i=0 ; i<data->num_chain ; i++)
    {
        MC2ERR_FREE(data->local_count[i]);
        MC2ERR_FREE(data->local_sum[i]);
    }
    for(size_t i=0, i_max=2*data->max_level*data->length ; i<i_max ; i++)
    {
        MC2ERR_FREE(data->pair_count[i]);
        MC2ERR_FREE(data->pair_sum[i]);
    }

    // free all remaining pointers
    MC2ERR_FREE(data->max_count);
    MC2ERR_FREE(data->max_pair);
    MC2ERR_FREE(data->num_level);
    MC2ERR_FREE(data->num_step);
    MC2ERR_FREE(data->local_count);
    MC2ERR_FREE(data->local_sum);
    MC2ERR_FREE(data->global_count);
    MC2ERR_FREE(data->global_sum);
    MC2ERR_FREE(data->pair_count);
    MC2ERR_FREE(data->pair_sum);

    // set sizes to zero for hygiene
    data->width = 0;
    data->length = 0;
    data->num_chain = 0;
    data->max_level = 0;

    // return without errors
    return 0;
}
