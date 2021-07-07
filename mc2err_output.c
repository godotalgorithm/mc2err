// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// analyze all available data in mc2err
//uint8_t mc2err_output(double acf_bias, double eqp_bias, struct mc2err *m2e)
uint8_t mc2err_output(uint64_t acf_cut, uint64_t eqp_cut, struct mc2err *m2e)
{
    // transfer & process all available data in block data buffers
    uint8_t status;
    for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
    {
        // transfer complete blocked data to the next buffer
        status = mc2err_transfer(i, m2e);
        if(status) { return status; }

        // process all available data in the buffer (head = tail)
        status = mc2err_process(i, m2e->buffer_tail[i], m2e);
        if(status) { return status; }
    }

    // expand data analysis buffers as necessary
    for(uint8_t i=0 ; i<=m2e->data_eqp ; i++)
    {
        MEM_CHECK(m2e->data_mean[i], double, (1 + !i)*m2e->eqp_num*m2e->data_dim);
        MEM_CHECK(m2e->data_error[i], double, (1 + !i)*m2e->eqp_num*m2e->data_dim*m2e->data_dim);
    }

    // clear data accumulation buffers
    for(uint8_t i=0 ; i<m2e->data_dim ; i++)
    { m2e->data_sum[i] = 0.0; }
 
    // clear ACF accumulation buffers
    for(uint8_t i=0 ; i<m2e->acf_pow && m2e->buffer[i]!=NULL ; i++)
    {
        MEM_CHECK(m2e->acf_sum[i], double, (1 + !i)*m2e->acf_num*m2e->data_dim*m2e->data_dim);
        MEM_CHECK(m2e->acf_mean[i], double, (1 + !i)*m2e->acf_num*m2e->data_dim*m2e->data_dim);

        for(uint64_t j=0 ; j<(1 + !i)*m2e->acf_num*m2e->data_dim*m2e->data_dim ; j++)
        { m2e->acf_sum[i][j] = 0.0; }
    }

    // countdown through active EQP
    uint8_t block = m2e->data_eqp;
    uint32_t index = m2e->data_num>>m2e->data_eqp - ((m2e->data_eqp == 0) ? 0 : m2e->eqp_num);
    do
    {
        // number of data points at this EQP
        uint64_t data_num0 = m2e->data_num - (index + !block)<<block;

        // update data sum & mean
        for(uint8_t i=0 ; i<m2e->data_dim ; i++)
        {
            m2e->data_sum[i] += m2e->data_partial[block][i+index*m2e->data_dim];
            m2e->data_mean[block][i+index*m2e->data_dim] = m2e->data_sum[i]/(double)data_num0;
        }

        // update ACF sum & mean
        uint16_t offset = m2e->data_dim*m2e->data_dim;
        for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
        for(uint32_t j=0 ; j<(1 + !i)*m2e->acf_num ; j++)
        for(uint16_t k=0 ; k<offset ; k++)
        {
            m2e->acf_sum[i][k+j*offset] += m2e->acf_partial[block][index][i][k+j*offset];
            m2e->acf_mean[i][k+j*offset] = m2e->acf_sum[i][k+j*offset]/(double)((data_num0>>i)-j);
        }

        // statistical estimate of the ACF truncation & error bar (presently neutered by acf_cut)
        mc2err_acf(acf_cut, block, index, m2e);

        // update indices
        if(index == 0)
        {
            if(block != 0)
            {
                block--;
                index = (1 + !block)*m2e->eqp_num;
            }
        }
        else
        { index--; }
    } while(block != 0 && index !=0);

    // statistical estimate of the EQP (presently neutered by eqp_cut)
    mc2err_eqp(eqp_cut, m2e);

    return 0;
}
