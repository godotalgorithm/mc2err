// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Trim the data vectors in the data accumulator 'source' to form the new data accumulator 'data' with data
// points of dimension 'width' and buffer size 'length' that is less than or equal to the buffer size of 'source'.
// The array 'map' of size 'width' contains the 0-based indices of the data vector elements that are kept in 'data'.
int mc2err_trim(struct mc2err_data *data, struct mc2err_data *source, int width, int length, int *map)
{
    // check for invalid arguments
    if(data == NULL || source == NULL || map == NULL)
    { return 1; }

    // check for consistency of sizes
    if(length > source->length)
    { return 2; }

    for(int i=0 ; i<width ; i++)
    {
        if(map[i] < 0 || map[i] > source->width-1)
        { return 2; }
    }

    // copy size information
    data->width = width;
    data->length = length;
    data->num_chain = source->num_chain;
    data->num_data = source->num_data;
    data->num_level = source->num_level;

    // allocate memory
    MC2ERR_MALLOC(data->chain_level, int, data->num_chain);
    MC2ERR_MALLOC(data->chain_count, long int, data->num_chain);
    MC2ERR_MALLOC(data->chain_sum, double*, data->num_chain);
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_MALLOC(data->chain_sum[i], double, 2*source->chain_level[i]*length*width); }
    size_t global_size = data->num_level*(data->num_level+3)*length/2;
    MC2ERR_MALLOC(data->data_count, long int, global_size);
    MC2ERR_MALLOC(data->data_sum, double, global_size*width);
    MC2ERR_MALLOC(data->pair_count, long int, global_size*2*length);
    MC2ERR_MALLOC(data->pair_sum, double, global_size*length*width*(width+1));

    // copy local data
    memcpy(data->chain_level, source->chain_level, sizeof(int)*data->num_chain);
    memcpy(data->chain_count, source->chain_count, sizeof(long int)*data->num_chain);
    for(int i=0 ; i<data->num_chain ; i++)
    for(int j=0 ; j<data->chain_level[i] ; j++)
    for(int k=0 ; k<2*length ; k++)
    for(int l=0 ; l<width ; l++)
    {
        data->chain_sum[i][(j*2*length+k)*width+l] =
        source->chain_sum[i][(j*2*source->length+k)*source->width+map[l]];
    }

    // copy global data
    for(int i=0 ; i<data->num_level ; i++)
    for(int j=0 ; j<data->num_level+1-i ; j++)
    {
        size_t k1 = i*(2*data->num_level+3-i)*data->length/2 + j*data->length;
        size_t k2 = i*(2*data->num_level+3-i)*source->length/2 + j*source->length;
        for(int k=0 ; k<data->length ; k++)
        {
            data->data_count[k1+k] = source->data_count[k2+k];

            for(int l=0 ; l<data->width ; l++)
            { data->data_sum[(k1+k)*data->width+l] = source->data_sum[(k2+k)*source->width+map[l]]; }

            for(int l=0 ; l<2*data->length ; l++)
            {
                data->pair_count[(k1+k)*2*data->length+l] = source->pair_count[(k2+k)*2*source->length+l];

                size_t m1 = ((k1+k)*2*data->length+l)*data->width*(data->width+1)/2;
                size_t m2 = ((k2+k)*2*source->length+l)*source->width*(source->width+1)/2;
                for(int m=0 ; m<data->width ; m++)
                for(int n=m ; n<data->width ; n++)
                {
                    int map_min = map[m], map_max = map[n];
                    if(map_min > map_max)
                    { int swap = map_min; map_min = map_max; map_max = swap; }

                    data->pair_sum[m1+m*(2*data->width-1-m)/2+n] =
                    source->pair_sum[m2+map_min*(2*source->width-1-map_min)/2+map_max];
                }
            }
        }
    }

    // return without errors
    return 0;
}
