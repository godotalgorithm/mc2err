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
    data->width = source->width;
    data->length = source->length;
    for(int i=1 ; i<num ; i++)
    {
        if(data->width != source[i].width || data->length != source[i].length)
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

    // allocate memory
    MC2ERR_MALLOC(data->chain_level, int, data->num_chain);
    MC2ERR_MALLOC(data->chain_count, long int, data->num_chain);
    MC2ERR_MALLOC(data->chain_sum, double*, data->num_chain);
    for(int i=0, j=0 ; i<num ; j+=source[i++].num_chain)
    {
        for(int k=0 ; k<source[i].num_chain ; k++)
        { MC2ERR_MALLOC(data->chain_sum[j+k], double, 2*source[i].chain_level[k]*data->length*data->width); }
    }
    size_t global_size = data->num_level*(data->num_level+3)*data->length/2;
    MC2ERR_MALLOC(data->data_count, long int, global_size);
    MC2ERR_MALLOC(data->data_sum, double, global_size*data->width);
    MC2ERR_MALLOC(data->pair_count, long int, global_size*2*data->length);
    MC2ERR_MALLOC(data->pair_sum, double, global_size*data->length*data->width*(data->width+1));

    // copy local data
    for(int i=0, j=0 ; i<num ; j+=source[i++].num_chain)
    {
        memcpy(data->chain_level+j, source[i].chain_level, sizeof(int)*source[i].num_chain);
        memcpy(data->chain_count+j, source[i].chain_count, sizeof(long int)*source[i].num_chain);
        for(int k=0 ; k<source[i].num_chain ; k++)
        {
            size_t chain_sum_size = sizeof(double)*2*source[i].chain_level[k]*data->length*data->width;
            memcpy(data->chain_sum[j+k], source[i].chain_sum[k], chain_sum_size);
        }
    }

    // initialize global data
    for(int i=0 ; i<global_size ; i++)
    { data->data_count[i] = 0; }

    for(int i=0 ; i<global_size*data->width ; i++)
    { data->data_sum[i] = 0.0; }

    for(int i=0 ; i<global_size*2*data->length ; i++)
    { data->pair_count[i] = 0; }

    for(int i=0 ; i<global_size*data->length*data->width*(data->width+1) ; i++)
    { data->pair_sum[i] = 0.0; }

    // merge global data
    for(int i=0 ; i<num ; i++)
    for(int j=0 ; j<source[i].num_level ; j++)
    {
        size_t k1 = j*(2*data->num_level+3-j)*data->length/2;
        size_t k2 = j*(2*source[i].num_level+3-j)*data->length/2;
        for(int k=0 ; k<(source[i].num_level+1-j)*data->length ; k++)
        {
            data->data_count[k1+k] += source[i].data_count[k2+k];

            for(int l=0 ; l<data->width ; l++)
            { data->data_sum[(k1+k)*data->width+l] += source[i].data_sum[(k2+k)*data->width+l]; }

            size_t l1 = (k1+k)*2*data->length;
            size_t l2 = (k2+k)*2*data->length;
            for(int l=0 ; l<2*data->length ; l++)
            {
                data->pair_count[l1+l] += source[i].pair_count[l2+l];

                size_t m1 = (l1+l)*data->width*(data->width+1)/2;
                size_t m2 = (l2+l)*data->width*(data->width+1)/2;
                for(int m=0 ; m<data->width*(data->width+1)/2 ; m++)
                { data->pair_sum[m1+m] += source[i].pair_sum[m2+m]; }
            }
        }
    }

    // return without errors
    return 0;
}
