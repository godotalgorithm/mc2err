// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Append all data from the data accumulator 'source' to the data accumulator 'data'.
// The chain indices from 'source' are offset by the number of Markov chains already in 'data'.
int mc2err_append(struct mc2err_data *data, struct mc2err_data *source)
{
    // check for invalid arguments
    if(data == NULL || source == NULL || data == source)
    { return 1; }

    // local copies of width & length for convenience
    int const width = data->width;
    int const length = data->length;

    // check for size consistency
    if(width != source->width || length != source->length)
    { return 3; }

    // check for overflow in the total number of chains
    if(data->num_chain > INT_MAX - source->num_chain)
    { return 7; }

    // check for overflow in the total number of accumulated data points
    for(int i=0 ; i<width ; i++)
    {
        if(data->max_count[i] > LONG_MAX - source->max_count[i])
        { return 7; }
    }

    // check for overflow in the total number of data pairs
    for(int i=0 ; i<width ; i++)
    {
        if(data->max_pair[i] > LLONG_MAX - source->max_pair[i])
        { return 7; }
    }

    // update max_level, reallocate & initialize global buffer as needed
    if(data->max_level < source->max_level)
    {
        // allocate more global buffers
        size_t old_size = 2*data->max_level*length;
        size_t new_size = 2*source->max_level*length;
        MC2ERR_REALLOC(data->global_count, long, new_size*width);
        MC2ERR_REALLOC(data->global_sum, double, new_size*width);
        MC2ERR_REALLOC(data->pair_count, long long*, new_size);
        MC2ERR_REALLOC(data->pair_sum, double*, new_size);
        for(size_t i=0 ; i<old_size ; i++)
        {
            MC2ERR_REALLOC(data->pair_count[i], long long, new_size*width*width);
            MC2ERR_REALLOC(data->pair_sum[i], double, new_size*width*width);
        }
        for(size_t i=old_size ; i<new_size ; i++)
        {
            MC2ERR_MALLOC(data->pair_count[i], long long, new_size*width*width);
            MC2ERR_MALLOC(data->pair_sum[i], double, new_size*width*width);
        }

        // initialize new global buffers to zero
        MC2ERR_FILL(data->global_count+old_size*width, long, (new_size-old_size)*width, 0);
        MC2ERR_FILL(data->global_sum+old_size*width, double, (new_size-old_size)*width, 0.0);
        for(size_t i=0 ; i<old_size ; i++)
        {
            MC2ERR_FILL(data->pair_count[i]+old_size*width*width, long long, (new_size-old_size)*width*width, 0);
            MC2ERR_FILL(data->pair_sum[i]+old_size*width*width, double, (new_size-old_size)*width*width, 0.0);
        }
        for(size_t i=old_size ; i<new_size ; i++)
        {
            MC2ERR_FILL(data->pair_count[i], long long, new_size*width*width, 0);
            MC2ERR_FILL(data->pair_sum[i], double, new_size*width*width, 0.0);
        }

        // update max_level
        data->max_level = source->max_level;
    }

    // update other size information
    if(data->max_step < source->max_step)
    { data->max_step = source->max_step; }
    for(int i=0 ; i<width ; i++)
    {
        data->max_count[i] += source->max_count[i];
        data->max_pair[i] += source->max_pair[i];
    }

    // append local data
    MC2ERR_REALLOC(data->num_level, int, data->num_chain+source->num_chain);
    MC2ERR_REALLOC(data->num_step, long, data->num_chain+source->num_chain);
    MC2ERR_REALLOC(data->local_count, long*, data->num_chain+source->num_chain);
    MC2ERR_REALLOC(data->local_sum, double*, data->num_chain+source->num_chain);
    memcpy(data->num_level+data->num_chain, source->num_level, sizeof(int)*source->num_chain);
    memcpy(data->num_step+data->num_chain, source->num_step, sizeof(long)*source->num_chain);
    for(int i=0 ; i<source->num_chain ; i++)
    {
        size_t size = 2*source->num_level[i]*length*width;
        MC2ERR_MALLOC(data->local_count[data->num_chain+i], long, size);
        MC2ERR_MALLOC(data->local_sum[data->num_chain+i], double, size);
        memcpy(data->local_count[data->num_chain+i], source->local_count[i], sizeof(long)*size);
        memcpy(data->local_sum[data->num_chain+i], source->local_sum[i], sizeof(double)*size);
    }
    data->num_chain += source->num_chain;

    // merge global data
    for(size_t i=0, i_max=2*source->max_level*length*width ; i<i_max ; i++)
    {
        data->global_count[i] += source->global_count[i];
        data->global_sum[i] += source->global_sum[i];
    }
    for(size_t i=0, i_max=2*source->max_level*length ; i<i_max ; i++)
    for(size_t j=0, j_max=2*source->max_level*length*width*width ; j<j_max ; j++)
    {
        data->pair_count[i][j] += source->pair_count[i][j];
        data->pair_sum[i][j] += source->pair_sum[i][j];
    }

    // return without errors
    return 0;
}
