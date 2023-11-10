// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: All possible data accumulation is performed whenever a data point is inputted,
//       rather than collecting data points in a buffer and processing the buffer when
//       it is full. Buffered data processing can be more efficient overall if it uses
//       level-3 BLAS operations instead of individual vector outer products, but then
//       the data accumulation cost is much less uniform and extra logic is needed to
//       track the buffer and process an incomplete buffer before performing statistical
//       analysis. Future implementations might benefit from the option of buffered input.

// Input the data point 'point' from the Markov chain with 0-based index 'chain' into the data accumulator 'data'.
int mc2err_input(struct mc2err_data *data, int chain, double *point)
{
    // check for invalid arguments
    if(data == NULL || chain < 0 || point == NULL)
    { return 1; }

    // check for data overflow
    if(data->num_data == LONG_MAX || chain == INT_MAX)
    { return 6; }

    // expand number of chains as needed
    if(chain >= data->num_chain)
    {
        MC2ERR_REALLOC(data->chain_level, int, chain+1);
        MC2ERR_REALLOC(data->chain_count, long int, chain+1);
        MC2ERR_REALLOC(data->chain_sum, double*, chain+1);

        for(int i=data->num_chain ; i<=chain ; i++)
        {
            data->chain_level[i] = 0;
            data->chain_count[i] = 0;
            data->chain_sum[i] = NULL;
        }

        data->num_chain = chain+1;
    }

    // expand local memory of a chain as needed
    if(data->chain_count[chain]+1 >= 1<<data->chain_level[chain])
    {
        MC2ERR_REALLOC(data->chain_sum[chain], double, 2*(data->chain_level[chain]+1)*data->length*data->width);
        data->chain_level[chain]++;
    }

    // expand global memory as needed
    if(data->chain_count[chain]+1 >= 1<<data->num_level)
    {
        size_t global_size = (data->num_level+1)*(data->num_level+4)*data->length/2;
        MC2ERR_REALLOC(data->data_count, long int, global_size);
        MC2ERR_REALLOC(data->data_sum, double, global_size*data->width);
        MC2ERR_REALLOC(data->data_count, long int, global_size*2*data->length);
        MC2ERR_REALLOC(data->data_count, long int, global_size*data->length*data->width*(data->width+1));

        // initialize new rows to zero
        size_t global_new = data->num_level*(data->num_level+5)*data->length/2;
        for(int i=0 ; i<global_size-global_new ; i++)
        { data->data_count[global_new+i] = 0; }
        for(int i=0 ; i<(global_size-global_new)*data->width ; i++)
        { data->data_sum[global_new*data->width+i] = 0.0; }
        for(int i=0 ; i<(global_size-global_new)*2*data->length ; i++)
        { data->pair_count[global_new*2*data->length+i] = 0; }
        for(int i=0 ; i<(global_size-global_new)*data->length*data->width*(data->width+1) ; i++)
        { data->pair_sum[global_new*data->length*data->width*(data->width+1)+i] = 0.0; }

        // adjust stride of old data & initialize new data in old rows to zero
        for(int i=data->num_level-1 ; i>=0 ; i--)
        {
            size_t i1 = i*(2*data->num_level+3-i)/2;
            size_t i2 = data->num_level+1-i;
            size_t num = data->length;
            memmove(data->data_count+(i1+i)*num, data->data_count+i1*num, sizeof(long int)*i2*num);
            for(int j=0 ; j<num ; j++)
            { data->data_count[(i1+i2+i)*num+j] = 0; }

            num = data->length*data->width;
            memmove(data->data_sum+(i1+i)*num, data->data_sum+i1*num, sizeof(double)*i2*num);
            for(int j=0 ; j<num ; j++)
            { data->data_sum[(i1+i2+i)*num+j] = 0.0; }

            num = 2*data->length*data->length;
            memmove(data->pair_count+(i1+i)*num, data->pair_count+i1*num, sizeof(long int)*i2*num);
            for(int j=0 ; j<num ; j++)
            { data-pair_count[(i1+i2+i)*num+j] = 0; }

            num = data->length*data->length*data->width*(data->width+1);
            memmove(data->pair_sum+(i1+i)*num, data->pair_sum+i1*num, sizeof(double)*i2*num);
            for(int j=0 ; j<num ; j++)
            { data-pair_sum[(i1+i2+i)*num+j] = 0.0; }
        }

        // update num_level
        data->num_level++;
    }

    // move new data point to the chain buffer
    int index = data->chain_count[chain]%(2*data->length);
    memcpy(data->chain_sum+index*data->width, point, sizeof(double)*data->width);

    // process new data in the chain buffer
    for(int i=0 ; i<data->chain_level[chain] ; i++)
    {
        // add data to mean vector accumulator

        // add pairwise products to pair matrix accumulator

        // merge data if a block is complete
        if(index%2 == 0)
        {
        }
        // otherwise halt the loop
        else { break; }
    }

    // increment the data counts
    data->chain_count[chain]++;
    data->num_data++;

    // return without errors
    return 0;
}
