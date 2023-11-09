// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: Vector addition is performed here with simple for loops, assuming that compiler optimizations
//       will unroll them if that is more efficient, and no attempt is made to use level-1 BLAS,
//       which might be more efficient on certain machines in certain size regimes.

// Merge the array of data accumulators, 'source', of size 'num' to form the new data accumulator 'data'.
// (The order of chains in 'source' is preserved in 'data', and chain indices are offset accordingly.)
int mc2err_merge(struct mc2err_data *data, struct mc2err_data *source, int num)
{
    // check for invalid arguments
    if(data == NULL || source == NULL || num < 1)
    { return 1; }

    // check for consistency
    for(int i=1 ; i<num ; i++)
    {
        if(source[i-1].width != source[i].width || source[i-1].length != source[i].length)
        { return 2; }
    }

    // check for overflow in the total number of chains
    int data->num_chain = 0;
    for(int i=0 ; i<num ; i++)
    {
        if(INT_MAX - source[i].num_chain < data->num_chain)
        { return 6; }

        data->num_chain += source[i].num_chain;
    }

    // check for overflow in the total number of accumulated data points
    long int data->num_data = 0;
    for(int i=0 ; i<num ; i++)
    {
        if(LONG_MAX - source[i].num_data < data->num_data)
        { return 6; }

        data->num_data += source[i].num_data;
    }

    // assign num_level
    data->num_level = source[0].num_level;
    for(int i=1 ; i<num ; i++)
    {
        if(source[i].num_level > data->num_level)
        { data->num_level = source[i].num_level; }
    }

    // allocate memory & copy local data
    data->chain_level = (int*)malloc(sizeof(int)*data->num_chain);
    if(data->chain_level == NULL) { return 4; }
    size_t offset = 0;
    for(int i=0 ; i<num ; i++)
    {
        memcpy(data->chain_level+offset, source[i].chain_level, sizeof(int)*source[i].num_chain);
        offset += source[i].num_chain;
    }

    data->chain_count = (long int*)malloc(sizeof(long int)*data->num_chain);
    if(data->chain_count == NULL) { return 4; }
    offset = 0;
    for(int i=0 ; i<num ; i++)
    {
        memcpy(data->chain_count+offset, source[i].chain_count, sizeof(long int)*source[i].num_chain);
        offset += source[i].num_chain;
    }

    data->chain_sum = (double**)malloc(sizeof(double*)*data->num_chain);
    if(data->chain_sum == NULL) { return 4; }
    offset = 0;
    for(int i=0 ; i<num ; i++)
    {
        for(int j=0 ; j<data->num_chain ; j++)
        {
            size_t chain_sum_size = 2*source[i].chain_level[j]*data->length*data->width;
            data->chain_sum[offset+j] = (double*)malloc(sizeof(double)*chain_sum_size);
            if(data->chain_sum[offset+j] == NULL) { return 4; }
            memcpy(data->chain_sum[offset+j], source[i].chain_sum[j], sizeof(double)*chain_sum_size);
        }
        offset += source[i].num_chain;
    }

    // allocate memory & initialize global data
    size_t num_global = data->num_level*(data->num_level+3)*data->length/2;
    data->data_count = (long int*)malloc(sizeof(long int)*num_global);
    if(data->data_count == NULL) { return 4; }
    for(int i=0 ; i<num_global ; i++)
    { data->data_count[i] = 0; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->width/2;
    data->data_sum = (double*)malloc(sizeof(double)*num_global);
    if(data->data_sum == NULL) { return 4; }
    for(int i=0 ; i<num_global ; i++)
    { data->data_sum[i] = 0.0; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->length;
    data->pair_count = (long int*)malloc(sizeof(long int)*num_global);
    if(data->pair_count == NULL) { return 4; }
    for(int i=0 ; i<num_global ; i++)
    { data->pair_count[i] = 0; }

    num_global = data->num_level*(data->num_level+3)*data->length*data->length*data->width*(data->width+1)/2;
    data->pair_sum = (double*)malloc(sizeof(double)*num_global);
    if(data->pair_sum == NULL) { return 4; }
    for(int i=0 ; i<num_global ; i++)
    { data->data_sum[i] = 0.0; }

    // merge global data
    size_t width2 = data->width*(data->width+1)/2;
    for(int i=0 ; i<num ; i++)
    {
        for(int j=0 ; j<source[i].num_level ; j++)
        {
            size_t index1 = j*(2*data->num_level+3-j)*data->length/2;
            size_t index2 = j*(2*source[i].num_level+3-j)*data->length/2;
            for(int k=0 ; k<(source[i].num_level+1-j)*data->length ; k++)
            {
                data->data_count[index1+k] += source[i].data_count[index2+k];
                for(int l=0 ; l<data->width ; l++)
                { data->data_sum[(index1+k)*data->width+l] += source[i].data_sum[(index2+k)*data->width+l]; }

                size_t index3 = (index1+k)*2*data->length;
                size_t index4 = (index2+k)*2*data->length;
                for(int l=0 ; l<2*data->length ; l++)
                {
                    data->pair_count[index3+l] += source[i].pair_count[index4+l];
                    for(int m=0 ; m<width2 ; m++)
                    { data->pair_sum[(index3+l)*width2+m] += source[i].pair_sum[(index4+l)*width2+m]; }
                }
            }
        }
    }

    // return without errors
    return 0;
}
