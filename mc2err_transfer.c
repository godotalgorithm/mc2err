// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// transfer data from one buffer to the next buffer
uint8_t mc2err_transfer(uint8_t index, struct mc2err *m2e)
{
    // nothing to do if we are in the last buffer
    if(index >= m2e->acf_pow-1)
    { return; }

    // identify the beginning & ending locations for the data transfer between buffers
    uint64_t offset = (m2e->buffer_transfer[index+1] > 2*m2e->acf_num) ? 2*m2e->acf_num : 0;
    uint64_t begin = m2e->buffer_transfer[index+1] - offset;
    uint64_t end = m2e->buffer_transfer[index]/2;

    // stop if there is no new data to propagate to the next buffer
    if(begin == end)
    { return; }

    // allocate memory for new buffers as needed
    MEM_CHECK(m2e->buffer[index+1], double, 4*m2e->acf_num*m2e->data_dim);

    // write the block data to the next buffer
    for(uint64_t i=begin ; i<end ; i++)
    for(uint8_t j=0 ; j<m2e->data_dim ; j++)
    {
        m2e->buffer[index+1][j+(offset+i)*m2e->data_dim] =
        0.5*(m2e->buffer[index][j+2*i*m2e->data_dim] + m2e->buffer[index][j+(2*i+1)*m2e->data_dim]);
    }

    // update the write position
    m2e->buffer_transfer[index+1] = offset + end;
    return 0;
}