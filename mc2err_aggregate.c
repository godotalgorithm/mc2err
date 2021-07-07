// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// aggregate data into partial sums
uint8_t mc2err_aggregate(uint64_t new_data_num, double *new_data, struct mc2err *m2e)
{
    while(new_data_num != 0)
    {
        // identify the number of data points to aggregate,
        // either all available data or until the EQP block is filled up
        uint32_t aggregate_num = MIN(new_data_num,
                                     1<<m2e->data_eqp - (m2e->data_num&((1<<m2e->data_eqp)-1)) );

        // allocation & initialization of data_partial
        MEM_CHECK_SET(m2e->data_partial[m2e->data_eqp], double, (1 + !m2e->data_eqp)*m2e->eqp_num*m2e->data_dim, 0.0);

        // perform the data aggregation
        uint32_t index = m2e->data_num>>m2e->data_eqp - (1 - !m2e->data_eqp)*m2e->eqp_num;
        for(uint32_t i=0 ; i<aggregate_num ; i++)
        for(uint8_t j=0 ; j<m2e->data_dim ; j++)
        { m2e->data_partial[m2e->data_eqp][j+index*m2e->data_dim] += new_data[j+i*m2e->data_dim]; }

        // update the aggregation indices
        new_data += aggregate_num*m2e->data_dim;
        new_data_num -= aggregate_num;
        m2e->data_num += aggregate_num;
        while( (m2e->data_num>>m2e->data_eqp) >= 2*m2e->eqp_num)
        { m2e->data_eqp++; } 
    }

    return 0;
}