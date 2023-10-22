// include details of the mc2err_data structure
#include "mc2err_internal.h"

// NOTE: All possible data accumulation is performed whenever a data point is inputted,
//       rather than collecting data points in a buffer and processing the buffer when
//       it is full. Buffered data processing can be more efficient overall if it uses
//       level-3 BLAS operations instead of individual vector outer products, but then
//       the data accumulation cost is much less uniform and extra logic is needed to
//       track the buffer and process an incomplete buffer before performing statistical
//       analysis. Future implementations might benefit from the option of buffered input.

// Input the data point 'data' from the Markov chain with index 'chain_index' into the mc2err data structure 'm2e'.
// (The first Markov chain index is 0.)
int mc2err_input(struct mc2err_data *m2e, int chain_index, double *data)
{
    // check for invalid arguments
    if(m2e == NULL || chain_index < 0)
    { return 1; }

    // check for data overflow
    if(m2e->num_data == LONG_MAX || chain_index == INT_MAX)
    { return 7; }

    // expand number of chains as needed
    if(chain_index >= m2e->num_chain)
    {
        m2e->chain_level = (int*)realloc(m2e->chain_level, sizeof(int)*(chain_index+1));
        if(m2e->chain_level == NULL) { return 5; }
        for(int i=m2e->num_chain ; i<=chain_index ; i++)
        { m2e->chain_level[i] = 0; }

        m2e->chain_count = (long int*)realloc(m2e->chain_count, sizeof(long int)*(chain_index+1));
        if(m2e->chain_count == NULL) { return 5; }
        for(int i=m2e->num_chain ; i<=chain_index ; i++)
        { m2e->chain_count[i] = 0; }

        m2e->chain_sum = (double**)realloc(m2e->chain_sum, sizeof(double*)*(chain_index+1));
        if(m2e->chain_sum == NULL) { return 5; }
        for(int i=m2e->num_chain ; i<=chain_index ; i++)
        { m2e->chain_sum[i] = NULL; }

        m2e->num_chain = chain_index+1;
    }

    // expand local memory of a chain as needed
    if(m2e->chain_count[chain_index]+1 >= 1<<m2e->chain_level[chain_index])
    {
        size_t size = sizeof(double)*2*(m2e->chain_level[chain_index]+1)*m2e->length*m2e->width;
        m2e->chain_sum[chain_index] = (double*)realloc(m2e->chain_sum[chain_index], size);
        if(m2e->chain_sum[chain_index] == NULL) { return 5; }

        m2e->chain_level[chain_index]++;
    }

    // expand global memory as needed
    if(m2e->chain_count[chain_index]+1 >= 1<<m2e->num_level)
    {
        int status = mc2err_expand(m2e, m2e->num_level+1);
        if(status != 0) { return status; }
    }

    // process data point
    {
        // add data to chain buffer
        // add data to mean vector accumulator
        // add pairwise products to pair accumulators (check for tail contributions)
        // merge data & move to the next block size if the next power of 2 is reached (loop over previous steps)
    }

    // return without errors
    return 0;
}

// mc2err header file
#include "mc2err.h"
#include "mc2err_internal.h"

// input a new data vector from a Markov chain into an mc2err data structure
uint8_t mc2err_input(uint64_t chain_index, double *data, struct mc2err_data *m2e)
{
    // identify the active EQP index
    uint8_t eqp_index = 0;
    if(chain_index < m2e->num_chain)
    {
        eqp_index = m2e->chain_block[chain_index]-1;
        if(m2e->chain_length[chain_index]+1 == 1<<m2e->chain_block[chain_index])
        { eqp_index++; }
    }

    // allocate memory if necessary
    int status = mc2err_expand(chain_index+1, eqp_index+1, m2e);
    if(status) { return status; }

    // propagate block data through the buffer & accumulate statistics
    uint8_t block_index = 0;
    uint8_t block_offset;
    uint16_t width = m2e->width;
    uint64_t data_index = m2e->chain_length[chain_index];
    do
    {
        // calculate the buffer offset of the active block data point
        block_offset = data_index&1;

        // adjust effective EQP index for the packed storage layout
        if(block_index > eqp_index)
        { eqp_index = block_index; }
        uint16_t new_index = block_index + eqp_index*(eqp_index+1)/2;

        // aliased pointers to old & new block data
        double *new_data = m2e->chain_buffer[chain_index] + (2*block_index + block_offset)*width;
        double *old_data = m2e->chain_buffer[chain_index] + (2*block_index + !block_offset)*width;

        // move data to the active block in the buffer
        if(block_index == 0)
        {
            for(uint16_t i=0 ; i<width ; i++)
            { new_data[i] = data[i]; }
        }
        else
        {
            for(uint16_t i=0 ; i<width ; i++)
            { new_data[i] += 0.5*(m2e->chain_buffer[chain_index][(2*block_index-1)*width]
                + m2e->chain_buffer[chain_index][(2*block_index-2)*width]); }
        }

        // accumulate self statistics
        for(uint16_t i=0 ; i<width ; i++)
        { m2e->data_sum1[new_index][i] += new_data[i]; }
        for(uint16_t i=0 ; i<width ; i++)
        for(uint16_t j=0 ; j<=i ; j++)
        { m2e->data_sum2[new_index][j+i*(i+1)/2] += new_data[i]*new_data[j]; }
        m2e->data_num[new_index]++;

        // accumulate pair statistics
        if(eqp_index > block_index)
        {
            // adjust the index of the old block data point
            uint16_t old_index = new_index;
            if(data_index == 1<<(eqp_index-block_index-1))
            { old_index -= eqp_index; }

            // accumulate lagged pair statistics
            for(uint16_t i=0 ; i<width ; i++)
            for(uint16_t j=0 ; j<=i ; j++)
            {
                m2e->lag_sum2[old_index][j+i*(i+1)/2] += 0.5*
                    (new_data[i]*old_data[j] + old_data[i]*new_data[j]);
            }
            m2e->lag_num[old_index]++;

            // accumulate Haar wavelet statistics
            if(block_offset == 1)
            {
                for(uint16_t i=0 ; i<width ; i++)
                { m2e->haar_sum1[old_index][i] += 0.5*(old_data[i] - new_data[i]); }
                for(uint16_t i=0 ; i<width ; i++)
                for(uint16_t j=0 ; j<=i ; j++)
                {
                    m2e->haar_sum2[old_index][j+i*(i+1)/2] += 0.25*
                        (old_data[i] - new_data[i])*(old_data[j] - new_data[j]);
                }
                m2e->haar_num[old_index]++;
            }
        }

        // update indices
        block_index++;
        data_index >>= 1;
    } while(block_offset == 1);

    // extend the active Markov chain
    m2e->chain_length[chain_index]++;
    m2e->num_input++;
    return 0;
}
