// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: This implementation was designed for simplicity over efficiency. A more complicated but more
//       efficient implementation would have an input buffer to delay the processing of data so that
//       data can be processed more efficiently in batches. Also, the local chain buffer could reduce
//       data movement by switching to a cyclic buffer that moves the head of the buffer as data is
//       added instead of shifting data and always adding data to the front of the buffer. Although
//       both changes could improve performance modestly, their implementations are too complicated
//       to be justified right now.

// Input the observable vector 'observable' from the Markov chain with index 'chain' into the data
// accumulator 'data'. Any missing elements of the observable vector should be recorded as NaN, and
// a completely empty observable vector can be input as a NULL pointer.
int mc2err_input(struct mc2err_data *data, int chain, double *observable)
{
    // check for invalid arguments
    if(data == NULL || chain < 0)
    { return 1; }

    // local copies of width & length for convenience
    int const width = data->width;
    int const length = data->length;

    // check for invalid data
    if(observable != NULL)
    {
        for(int i=0 ; i<width ; i++)
        {
            if(isinf(observable[i])) { return 2; }
        }
    }

    // check for data overflows
    if(chain == INT_MAX)
    { return 7; }
    if(chain < data->num_chain && data->num_step[chain] == LONG_MAX)
    { return 7; }
    if(observable != NULL)
    {
        for(int i=0 ; i<width ; i++)
        {
            if(isnan(observable[i])) { continue; }
            if(data->max_count[i] == LONG_MAX) { return 7; }
            if(data->max_pair[i] >= LLONG_MAX - data->max_count[i]) { return 7; }
        }
    }

    // expand number of chains as needed
    if(chain >= data->num_chain)
    {
        // expand memory footprint of chain list
        MC2ERR_REALLOC(data->num_level, int, chain+1);
        MC2ERR_REALLOC(data->num_step, long, chain+1);
        MC2ERR_REALLOC(data->local_count, long*, chain+1);
        MC2ERR_REALLOC(data->local_sum, double*, chain+1);

        // initialize empty chains
        MC2ERR_FILL(data->num_level+data->num_chain, int, chain-data->num_chain+1, 0);
        MC2ERR_FILL(data->num_step+data->num_chain, long, chain-data->num_chain+1, 0);
        MC2ERR_FILL(data->local_count+data->num_chain, long*, chain-data->num_chain+1, NULL);
        MC2ERR_FILL(data->local_sum+data->num_chain, double*, chain-data->num_chain+1, NULL);

        // initialize new chain
        data->num_level[chain] = 1;
        MC2ERR_MALLOC(data->local_count[chain], long, 2*length*width);
        MC2ERR_MALLOC(data->local_sum[chain], double, 2*length*width);
        MC2ERR_FILL(data->local_count[chain], long, 2*length*width, 0);
        MC2ERR_FILL(data->local_sum[chain], double, 2*length*width, 0.0);

        // update num_chain
        data->num_chain = chain+1;
    }

    // expand local memory of a chain as needed
    if(data->num_step[chain]<<1 == 1<<data->num_level[chain])
    {
        size_t new_size = 2*(data->num_level[chain]+1)*length*width;
        MC2ERR_REALLOC(data->local_count[chain], long, new_size);
        MC2ERR_REALLOC(data->local_sum[chain], double, new_size);

        // initialize expanded local buffer
        size_t old_size = 2*data->num_level[chain]*length*width;
        MC2ERR_FILL(data->local_count[chain]+old_size, long, new_size-old_size, 0);
        MC2ERR_FILL(data->local_sum[chain]+old_size, double, new_size-old_size, 0.0);

        // fill front of new local buffer with data from previous coarse-graining level
        size_t offset = 2*(data->num_level[chain]-1)*length*width;
        memcpy(data->local_count[chain]+old_size, data->local_count[chain]+offset, sizeof(long)*width);
        memcpy(data->local_sum[chain]+old_size, data->local_sum[chain]+offset, sizeof(double)*width);

        // update num_level
        data->num_level[chain]++;
    }

    // expand global & pair buffers as needed
    if(data->max_level < data->num_level[chain])
    {
        // expand global buffer
        size_t new_size = 2*(data->max_level+1)*length;
        MC2ERR_REALLOC(data->global_count, long, new_size*width);
        MC2ERR_REALLOC(data->global_sum, double, new_size*width);

        // initialize new global buffer to zero
        size_t old_size = 2*data->max_level*length;
        MC2ERR_FILL(data->global_count+old_size*width, long, (new_size-old_size)*width, 0);
        MC2ERR_FILL(data->global_sum+old_size*width, double, (new_size-old_size)*width, 0.0);

        // expand & initialize pair buffer
        MC2ERR_REALLOC(data->pair_count, long long*, new_size);
        MC2ERR_REALLOC(data->pair_sum, double*, new_size);
        for(int i=0 ; i<data->max_level ; i++)
        for(int j=0 ; j<2*length ; j++)
        {
            MC2ERR_REALLOC(data->pair_count[2*length*i+j], long long, new_size*width*width);
            MC2ERR_REALLOC(data->pair_sum[2*length*i+j], double, new_size*width*width);
            MC2ERR_FILL(data->pair_count[2*length*i+j]+2*(data->max_level-i)*length*width*width,
                long long, (new_size-old_size)*width*width, 0);
            MC2ERR_FILL(data->pair_sum[2*length*i+j]+2*(data->max_level-i)*length*width*width,
                double, (new_size-old_size)*width*width, 0.0);
        }
        for(int i=0 ; i<2*length ; i++)
        {
            MC2ERR_MALLOC(data->pair_count[2*length*data->max_level+i], long long, 2*length*width*width);
            MC2ERR_MALLOC(data->pair_sum[2*length*data->max_level+i], double, 2*length*width*width);
            MC2ERR_FILL(data->pair_count[2*length*data->max_level+i], long long, 2*length*width*width, 0);
            MC2ERR_FILL(data->pair_sum[2*length*data->max_level+i], double, 2*length*width*width, 0.0);
        }

        // fill front of new global buffers with data from previous coarse-graining level
        if(data->max_level > 0)
        {
            size_t offset = 2*(data->max_level-1)*length;
            memcpy(data->global_count+old_size*width, data->global_count+offset*width, sizeof(long)*width);
            memcpy(data->global_sum+old_size*width, data->global_sum+offset*width, sizeof(double)*width);
            memcpy(data->pair_count[offset]+2*length*width*width, data->pair_count[offset], sizeof(long long)*width*width);
            memcpy(data->pair_sum[offset]+2*length*width*width, data->pair_sum[offset], sizeof(double)*width*width);
            memcpy(data->pair_count[old_size], data->pair_count[offset], sizeof(long long)*width*width);
            memcpy(data->pair_sum[old_size], data->pair_sum[offset], sizeof(double)*width*width);
        }

        // update max_level
        data->max_level++;
    }
    const int max_level = data->max_level;

    // shift data in local buffer
    for(int i=0 ; i<data->num_level[i] ; i++)
    {
        // shift data in local buffer by one block
        size_t offset = 2*i*length*width;
        memmove(data->local_count[chain]+offset+width, data->local_count[chain]+offset, sizeof(long)*(2*length-1)*width);
        memmove(data->local_sum[chain]+offset+width, data->local_sum[chain]+offset, sizeof(double)*(2*length-1)*width);

        // fill front of local buffer
        MC2ERR_FILL(data->local_count[chain]+offset, long, width, 0);
        MC2ERR_FILL(data->local_sum[chain]+offset, double, width, 0.0);

        // criteria to stop shifting
        if((data->num_step[chain]>>i)&1)
        { break; }
    }

    // add new data to all buffers if there is any
    if(observable != NULL)
    {
        // add data to local buffer
        for(int i=0 ; i<data->num_level[chain] ; i++)
        {
            size_t offset = 2*i*length*width;
            for(int j=0 ; j<width ; j++)
            {
                if(isnan(observable[j])) { continue; }
                data->local_count[chain][offset+j]++;
                data->local_sum[chain][offset+j] += observable[j];
            }
        }

        // add data to global buffer
        for(int i=max_level-1 ; i>=0 ; i--)
        {
            // offset & shift for the coarse-graining level
            size_t offset = 2*i*length;
            long shift = data->num_step[chain]/(1<<i);
            if(shift >= 2*length)
            { break; }

            // accumulate the average
            for(int j=0 ; j<width ; j++)
            {
                if(isnan(observable[j])) { continue; }
                data->global_count[(offset+shift)*width+j]++;
                data->global_sum[(offset+shift)*width+j] += observable[j];
            }
        }

        // add data to pair buffer
        for(int i=0 ; i<max_level ; i++) // loop over ACC level
        {
            int local_level = (i < data->num_level[chain]) ? i : data->num_level[chain];
            long local_max = ((data->num_step[chain]/(1<<i)) < 2*length-1) ? data->num_step[chain]/(1<<i) : 2*length-1;
            for(int j=0 ; j<local_max ; j++) // loop over ACC offset
            for(int k=max_level-1 ; k>=i ; k--) // loop over EQP level
            {
                // offset & shift for the coarse-graining level
                size_t offset = 2*k*length;
                long shift = (data->num_step[chain] - j*(1<<i))/(1<<k);
                if(shift >= 2*length)
                { break; }

                // accumulate the covariance
                for(int l=0 ; l<width ; l++)
                for(int m=0 ; m<width ; m++)
                {
                    if(isnan(observable[l])) { continue; }
                    data->pair_count[2*length*i+j][(offset+shift)*width*width+width*l+m] +=
                        data->local_count[chain][(2*length*local_level+j)*width+m];
                    data->pair_sum[2*length*i+j][(offset+shift)*width*width+width*l+m] +=
                        observable[l]*data->local_sum[chain][(2*length*local_level+j)*width+m];
                }
            }
        }

        // update total number of data points
        for(int i=0 ; i<width ; i++)
        {
            if(isnan(observable[i])) { continue; }
            data->max_count[i]++;
            data->max_pair[i] += data->max_count[i];
        }
    }

    // update number of steps
    data->num_step[chain]++;
    if(data->num_step[chain] > data->max_step)
    { data->max_step = data->num_step[chain]; }

    // return without errors
    return 0;
}
