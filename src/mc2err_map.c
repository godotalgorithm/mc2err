// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Map the data accumulator 'source' to form the new data accumulator 'data' for observable vectors of
// dimension 'width' and the smaller or equal buffer size 'length'. The vector 'index' of dimension
// 'width' contains the indices of the observable vectors from 'source' that are kept in 'data', and
// any out-of-bounds indices correspond to new observables with no previously recorded data.
int mc2err_map(struct mc2err_data *data, struct mc2err_data *source, int width, int length, int *index)
{
    // check for invalid arguments
    if(data == NULL || source == NULL || data == source || index == NULL || width < 1 || length < 1)
    { return 1; }

    // check for consistency of sizes
    if(length > source->length)
    { return 2; }

    // copy size information
    data->width = width;
    data->length = length;
    data->num_chain = source->num_chain;
    data->max_level = source->max_level;
    data->max_step = source->max_step;
    MC2ERR_MALLOC(data->max_count, long, width);
    MC2ERR_MALLOC(data->max_pair, long long, width);
    memcpy(data->max_step, source->max_step, sizeof(long)*width);
    memcpy(data->max_pair, source->max_pair, sizeof(long long)*width);

    // allocate local memory
    MC2ERR_MALLOC(data->num_level, int, data->num_chain);
    MC2ERR_MALLOC(data->num_step, long, data->num_chain);
    MC2ERR_MALLOC(data->local_count, long*, data->num_chain);
    MC2ERR_MALLOC(data->local_sum, double*, data->num_chain);
    for(int i=0 ; i<data->num_chain ; i++)
    {
        MC2ERR_MALLOC(data->local_count[i], long, 2*data->num_level[i]*length*width);
        MC2ERR_MALLOC(data->local_sum[i], double, 2*data->num_level[i]*length*width);
    }

    // allocate global buffer
    MC2ERR_MALLOC(data->global_count, long, 2*data->max_level*length*width);
    MC2ERR_MALLOC(data->global_sum, double, 2*data->max_level*length*width);
    MC2ERR_MALLOC(data->pair_count, long long*, 2*data->max_level*length);
    MC2ERR_MALLOC(data->pair_sum, double*, 2*data->max_level*length);
    for(size_t i=0, i_max=2*data->max_level*length ; i<i_max ; i++)
    {
        MC2ERR_MALLOC(data->pair_count[i], long long, i_max*width*width);
        MC2ERR_MALLOC(data->pair_sum[i], double, i_max*width*width);
    }

    // transfer local data
    memcpy(data->num_level, source->num_level, sizeof(int)*data->num_chain);
    memcpy(data->num_step, source->num_step, sizeof(long)*data->num_chain);

    // map data in the local chain buffers
    for(int i=0 ; i<data->num_chain ; i++)
    for(int j=0 ; j<data->num_level[i] ; j++)
    for(int k=0 ; k<2*length ; k++)
    {
        size_t data_offset = (j*2*length + k)*width;
        size_t source_offset = (j*2*source->length + k)*source->width;
        for(int l=0 ; l<width ; l++)
        {
            if(index[l] >= 0 && index[l] < source->width)
            {
                data->local_count[i][data_offset+l] = source->local_count[i][source_offset+index[l]];
                data->local_sum[i][data_offset+l] = source->local_sum[i][source_offset+index[l]];
            }
            else
            {
                data->local_count[i][data_offset+l] = 0;
                data->local_sum[i][data_offset+l] = 0.0;
            }
        }
    }

    // map data in global buffer
    for(int i=0 ; i<data->max_level ; i++)
    for(int j=0 ; j<2*length ; j++)
    {
        size_t data_offset = (i*2*length + j)*width;
        size_t source_offset = (i*2*source->length + j)*source->width;
        for(int k=0 ; k<width; k++)
        {
            if(index[k] >= 0 && index[k] < source->width)
            {
                data->global_count[data_offset+k] = source->global_count[source_offset+index[k]];
                data->global_sum[data_offset+k] = source->global_sum[source_offset+index[k]];
            }
            else
            {
                data->global_count[data_offset+k] = 0;
                data->global_sum[data_offset+k] = 0.0;
            }
        }
    }

    // map data in pair buffer
    for(int i=0 ; i<data->max_level ; i++)
    for(int j=0 ; j<2*length ; j++)
    {
        long long *data_count_ptr = data->pair_count[(i*2*length + j)*width];
        double *data_sum_ptr = data->pair_sum[(i*2*length + j)*width];
        long long *source_count_ptr = source->pair_count[(i*2*source->length + j)*source->width];
        double *source_sum_ptr = source->pair_sum[(i*2*source->length + j)*source->width];
        for(int k=0 ; k<data->max_level ; k++)
        for(int l=0 ; l<2*length ; l++)
        {
            size_t data_offset = (k*2*length + l)*width;
            size_t source_offset = (k*2*source->length + l)*source->width;
            for(int m=0 ; m<width; m++)
            {
                if(index[m] >= 0 && index[m] < source->width)
                {
                    data_count_ptr[data_offset+m] = source_count_ptr[source_offset+index[m]];
                    data_sum_ptr[data_offset+m] = source_sum_ptr[source_offset+index[m]];
                }
                else
                {
                    data_count_ptr[data_offset+m] = 0;
                    data_sum_ptr[data_offset+m] = 0.0;
                }
            }
        }
    }

    // return without errors
    return 0;
}
