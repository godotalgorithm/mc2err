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
