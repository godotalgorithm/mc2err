// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// process data up to a new head value, shift buffer when it is half-processed
uint8_t mc2err_process(uint8_t index, uint64_t new_head, struct mc2err *m2e)
{
    // update for new tail w/ old head
    if(m2e->buffer_tail[index] < m2e->buffer_transfer[index])
    {

    }

    // update for new head
    if(new_head < m2e->buffer_head[index])
    {

    }

    // only process buffers with work to do
    if(m2e->buffer_transfer[index] < 2*m2e->acf_num)
    { return 0; }

    // loop over the data in the 1st half of the buffer
    while(m2e->buffer_process[index] < m2e->buffer_transfer[index] - m2e->acf_num)
    {
        // expand storage space for ACF partial sums as needed
        uint8_t status = mc2err_expand(index, m2e);
        if(status) { return status; }

        // identify the amount of data to be processed for the active EQP
        uint32_t process_num = MIN(m2e->acf_num - m2e->buffer_process[index],
                                   1 + offset_max - m2e->buffer_offset[index]);

        // update acf_partial
        convolution(m2e->data_dim,
                    process_num,
                    process_num + m2e->acf_num,
                    m2e->acf_num,
                    m2e->buffer[0] + m2e->buffer_process[0]*m2e->data_dim,
                    m2e->buffer[0] + m2e->buffer_process[0]*m2e->data_dim,
                    m2e->acf_partial[eqp_grid][eqp_sub][0]);

        // update the processing indices
        m2e->buffer_process[index] += process_num;
    }

    // shift a buffer when half of the data has been processed
    if(m2e->buffer_head[index] >= 2*m2e->acf_num && m2e->buffer_tail[index] == 4*m2e->acf_num)
    {
        // copy 2nd half of the buffer to the 1st half
        memcpy(m2e->buffer[index], m2e->buffer[index]+2*m2e->acf_num*m2e->data_dim,
               2*m2e->acf_num*m2e->data_dim*sizeof(double));

        // shift buffer indices to coincide w/ shifted buffer
        m2e->buffer_transfer[index] -= 2*m2e->acf_num;
        m2e->buffer_head[index] -= 2*m2e->acf_num;
        m2e->buffer_tail[index] -= 2*m2e->acf_num;
        m2e->buffer_offset[index] += 2*m2e->acf_num;
        while(m2e->buffer_offset[index])

        // shift EQP indices
        uint32_t shift_num = 2*m2e->acf_num;
        while(shift_num != 0)
        {
            m2e->buffer_eqp1[index] +=;
            m2e->buffer_eqp2[index] +=;
        }
    }

    return 0;
}

// expand memory allocations to accommodate new ACF partial sums
uint8_t mc2err_expand(uint8_t index, struct mc2err *m2e)
{
    // 1st-level allocation of acf_partial
    if(m2e->acf_partial[m2e->buffer_eqp1[index]] == NULL)
    {
        m2e->acf_partial[m2e->buffer_eqp1[index]] =
        (double***)malloc(sizeof(double**)*m2e->eqp_num*(1 + !m2e->buffer_eqp1[index]));
        if(m2e->acf_partial[m2e->buffer_eqp1[index]] == NULL) { return 3; }
    }

    // 2nd-level allocation of acf_partial
    if(m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]] == NULL)
    {
        m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]] = (double**)malloc(sizeof(double*)*m2e->acf_pow);
        if(m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]] == NULL) { return 3; }
    }

    // 3rd-level allocations of acf_partial
    if(m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]][index] == NULL)
    {
        uint64_t size = (1 + !index)*m2e->acf_num*m2e->data_dim*m2e->data_dim;
        m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]][index] = (double*)malloc(sizeof(double)*size);
        if(m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]][index] == NULL) { return 3; }

        for(uint64_t i=0 ; i<size ; i++)
        { m2e->acf_partial[m2e->buffer_eqp1[index]][m2e->buffer_eqp2[index]][index][i] = 0.0; }
    }

    return 0;
}