// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Begin the sampling process by initializing the new data accumulator 'data' for
// observable vectors of dimension 'width' and for accumulation buffers of size 'length'.
int mc2err_begin(struct mc2err_data *data, int width, int length)
{
    // check for invalid arguments
    if(data == NULL || width < 1 || length < 1)
    { return 1; }

    // pass through width & length
    data->width = width;
    data->length = length;

    // initial memory allocation
    MC2ERR_MALLOC(data->max_count, long, width);
    MC2ERR_MALLOC(data->max_pair, long long, width);

    // initialize sizes to 0
    data->num_chain = 0;
    data->max_level = 0;
    data->max_step = 0;
    MC2ERR_FILL(data->max_count, long, width, 0);
    MC2ERR_FILL(data->max_pair, long long, width, 0);

    // initialize remaining pointers to NULL
    data->num_level = NULL;
    data->num_step = NULL;
    data->local_count = NULL;
    data->local_sum = NULL;
    data->global_count = NULL;
    data->global_sum = NULL;
    data->pair_count = NULL;
    data->pair_sum = NULL;

    // return without errors
    return 0;
}
