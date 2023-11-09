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

    // allocate memory for local data & copy it from source
    data->chain_level = (int*)malloc(sizeof(int)*data->num_chain);
    if(data->chain_level == NULL) { return 4; }
    memcpy(data->chain_level, source->chain_level, sizeof(int)*data->num_chain);

    data->chain_count = (long int*)malloc(sizeof(long int)*data->num_chain);
    if(data->chain_count == NULL) { return 4; }
    memcpy(data->chain_count, source->chain_count, sizeof(long int)*data->num_chain);

    data->chain_sum = (double**)malloc(sizeof(double*)*data->num_chain);
    if(data->chain_sum == NULL) { return 4; }
    for(int i=0 ; i<data->num_chain ; i++)
    {
        size_t chain_sum_size = 2*data->chain_level[i]*data->length*data->width;
        data->chain_sum[i] = (double*)malloc(sizeof(double)*chain_sum_size);
        if(data->chain_sum[i] == NULL) { return 4; }
        for(int j=0 ; j<data->chain_level[i] ; j++)
        for(int k=0 ; k<data->length ; k++)
        for(int l=0 ; l<data->width ; l++)
        {
            data->chain_sum[i][(j*data->length+k)*data->width+l] =
            source->chain_sum[i][(j*source->length+k)*source->width+map[l]];
        }
    }

    // allocate memory for global data
    size_t num_global = data->num_level*(data->num_level+3)*data->length/2;
    data->data_count = (long int*)malloc(sizeof(long int)*num_global);
    if(data->data_count == NULL) { return 4; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->width/2;
    data->data_sum = (double*)malloc(sizeof(double)*num_global);
    if(data->data_sum == NULL) { return 4; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->length;
    data->pair_count = (long int*)malloc(sizeof(long int)*num_global);
    if(data->pair_count == NULL) { return 4; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->length*data->width*(data->width+1)/2;
    data->pair_sum = (double*)malloc(sizeof(double)*num_global);
    if(data->pair_sum == NULL) { return 4; }

    // copy global data
    size_t width2 = data->width*(data->width+1)/2;
    size_t width3 = source->width*(source->width+1)/2;
    for(int i=0 ; i<data->num_level ; i++)
    {
        size_t index1 = i*(2*data->num_level+3-i)*data->length/2;
        size_t index2 = i*(2*data->num_level+3-i)*source->length/2;
        for(int j=0 ; j<data->num_level+1-i ; j++)
        for(int k=0 ; k<data->length ; k++)
        {
            data->data_count[index1+j*data->length+k] = source->data_count[index2+j*source->length+k];
            for(int l=0 ; l<data->width ; l++)
            {
                data->data_sum[(index1+j*data->length+k)*data->width+l] =
                source->data_sum[(index2+j*source->length+k)*source->width+map[l]];
            }

            size_t index3 = (index1+j*data->length+k)*2*data->length;
            size_t index4 = (index2+j*source->length+k)*2*source->length;
            for(int l=0 ; l<2*data->length ; l++)
            {
                data->pair_count[index3+l] = source->pair_count[index4+l];
                for(int m=0 ; m<data->width ; m++)
                for(int n=m ; n<data->width ; n++)
                {
                    int map_min = map[m];
                    int map_max = map[n];
                    if(map_min > map_max)
                    { int swap = map_min; map_min = map_max; map_max = swap; }
                    data->pair_sum[(index3+k)*width2+m*(2*data->width-1-m)/2+n] =
                    source->pair_sum[(index4+k)*width3+map_min*(2*source->width-1-map_min)/2+map_max];
                }
            }
        }
    }

    // return without errors
    return 0;
}
