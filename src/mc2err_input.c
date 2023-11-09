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
        data->chain_level = (int*)realloc(data->chain_level, sizeof(int)*(chain+1));
        if(data->chain_level == NULL) { return 4; }
        for(int i=data->num_chain ; i<=chain ; i++)
        { data->chain_level[i] = 0; }

        data->chain_count = (long int*)realloc(data->chain_count, sizeof(long int)*(chain+1));
        if(data->chain_count == NULL) { return 4; }
        for(int i=data->num_chain ; i<=chain ; i++)
        { data->chain_count[i] = 0; }

        data->chain_sum = (double**)realloc(data->chain_sum, sizeof(double*)*(chain+1));
        if(data->chain_sum == NULL) { return 4; }
        for(int i=data->num_chain ; i<=chain ; i++)
        { data->chain_sum[i] = NULL; }

        data->num_chain = chain+1;
    }

    // expand local memory of a chain as needed
    if(data->chain_count[chain]+1 >= 1<<data->chain_level[chain])
    {
        size_t new_num = 2*(data->chain_level[chain]+1)*data->length*data->width;
        data->chain_sum[chain] = (double*)realloc(data->chain_sum[chain], sizeof(double)*new_num);
        if(data->chain_sum[chain] == NULL) { return 4; }

        data->chain_level[chain]++;
    }

    // expand global memory as needed
    if(data->chain_count[chain]+1 >= 1<<data->num_level)
    {
        // expand memory & initialize new rows to zero
        size_t new_row = data->num_level*(data->num_level+5)*data->length/2;
        size_t new_num = (data->num_level+1)*(data->num_level+4)*data->length/2;
        data->data_count = (long int*)realloc(data->data_count, sizeof(long int)*new_num);
        if(data->data_count == NULL) { return 4; }
        for(int i=new_row ; i<new_num ; i++)
        { data->data_count[i] = 0; }

        new_row = data->num_level*(data->num_level+5)*data->length*data->width/2;
        new_num = (data->num_level+1)*(data->num_level+4)*data->length*data->width/2;
        data->data_sum = (double*)realloc(data->data_sum, sizeof(double)*new_num);
        if(data->data_sum == NULL) { return 4; }
        for(int i=new_row ; i<new_num ; i++)
        { data->data_sum[i] = 0.0; }

        new_row = data->num_level*(data->num_level+5)*data->length*data->length;
        new_num = (data->num_level+1)*(data->num_level+4)*data->length*data->length;
        data->pair_count = (long int*)realloc(data->pair_count, sizeof(long int)*new_num);
        if(data->pair_count == NULL) { return 4; }
        for(int i=new_row ; i<new_num ; i++)
        { data->pair_count[i] = 0; }

        new_row = data->num_level*(data->num_level+5)*data->length*data->length*data->width*(data->width+1)/2;
        new_num = (data->num_level+1)*(data->num_level+4)*data->length*data->length*data->width*(data->width+1)/2;
        data->pair_sum = (double*)realloc(data->pair_sum, sizeof(double)*new_num);
        if(data->pair_sum == NULL) { return 4; }
        for(int i=new_row ; i<new_num ; i++)
        { data->pair_sum[i] = 0.0; }

        // adjust stride of old data & initialize new data in old rows to zero
        for(int i=data->num_level-1 ; i>=0 ; i--)
        {
            size_t index1 = i*(2*data->num_level+3-i)/2;
            size_t index2 = data->num_level+1-i;
            size_t index3 = data->length;
            memmove(data->data_count+(index1+i)*index3, data->data_count+index1*index3, sizeof(long int)*index2*index3);
            for(int j=0 ; j<index3 ; j++)
            { data->data_count[(index1+index2+i)*index3+j] = 0; }

            index3 = data->length*data->width;
            memmove(data->data_sum+(index1+i)*index3, data->data_sum+index1*index3, sizeof(double)*index2*index3);
            for(int j=0 ; j<index3 ; j++)
            { data->data_sum[(index1+index2+i)*index3+j] = 0.0; }

            index3 = 2*data->length*data->length;
            memmove(data->pair_count+(index1+i)*index3, data->pair_count+index1*index3, sizeof(long int)*index2*index3);
            for(int j=0 ; j<index3 ; j++)
            { data-pair_count[(index1+index2+i)*index3+j] = 0; }

            index3 = 2*data->length*data->length*data->width*(data->width+1)/2;
            memmove(data->pair_sum+(index1+i)*index3, data->pair_sum+index1*index3, sizeof(double)*index2*index3);
            for(int j=0 ; j<index3 ; j++)
            { data-pair_sum[(index1+index2+i)*index3+j] = 0.0; }

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
