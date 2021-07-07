// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// input new MCMC data into the mc2err data structure
uint8_t mc2err_input(uint64_t new_data_num, double *new_data, struct mc2err *m2e)
{
    // check if there are too many data points
    if(new_data_num > m2e->data_max - m2e->data_num)
    { return 4; }

    // primary aggregation of partial sums for the mean
    uint8_t status = mc2err_aggregate(new_data_num, new_data, m2e);
    if(status) { return status; }

    // propagate data through the buffers
    while(new_data_num != 0)
    {
        // identify the number of data points to copy,
        // either copy all remaining data or fill up the primary data buffer
        uint32_t copy_num = MIN(new_data_num, 4*m2e->acf_num - m2e->buffer_transfer[0]);

        // use memcpy to copy data efficiently to the 1st buffer
        memcpy(m2e->buffer[0] + m2e->data_dim*m2e->buffer_transfer[0], new_data, copy_num*m2e->data_dim*sizeof(double));

        // update indices & pointers
        m2e->buffer_transfer[0] += copy_num;
        new_data_num -= copy_num;
        new_data += copy_num*m2e->data_dim;

        // loop over filled data buffers
        for(uint8_t i=0 ; i<m2e->acf_pow && m2e->buffer_transfer[i]==4*m2e->acf_num ; i++)
        {
            // transfer data between buffers
            status = mc2err_transfer(i, m2e);
            if(status) { return status; }

            // process all data from the 1st half of the buffer
            status = mc2err_process(i, 2*m2e->acf_num, m2e);
            if(status) { return status; }
        }
    }
    return 0;
}