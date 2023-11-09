// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Initialize the mc2err data accumulator 'data' for data points of dimension 'width' and buffer size 'length'.
int mc2err_initialize(struct mc2err_data *data, int width, int length)
{
    // check for invalid arguments
    if(data == NULL || width < 1 || length < 1)
    { return 1; }

    // pass through width & length
    data->width = width;
    data->length = length;

    // initialize sizes & info to 0
    data->num_level = 0;
    data->num_chain = 0;
    data->num_data = 0;

    // initialize all pointers to NULL
    data->chain_level = NULL;
    data->chain_count = NULL;
    data->chain_sum = NULL;
    data->data_count = NULL;
    data->data_sum = NULL;
    data->pair_count = NULL;
    data->pair_sum = NULL;

    // return without errors
    return 0;
}
