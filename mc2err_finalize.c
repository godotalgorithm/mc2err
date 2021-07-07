// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// deallocate all of the memory used by mc2err
void mc2err_finalize(struct mc2err *m2e)
{
    free(m2e->buffer_eqp);
    free(m2e->buffer_offset);
    free(m2e->buffer_transfer);
    free(m2e->buffer_head);
    free(m2e->buffer_tail);

    for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
    { free(m2e->buffer[i]); }
    free(m2e->buffer);

    free(m2e->data_sum);
    for(uint8_t i=0 ; i<m2e->eqp_pow ; i++)
    {
        free(m2e->data_mean[i]);
        free(m2e->data_partial[i]);
        free(m2e->data_error[i]);
    }
    free(m2e->data_mean);
    free(m2e->data_partial);
    free(m2e->data_error);

    for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
    {
        free(m2e->acf_sum[i]);
        free(m2e->acf_mean[i]);
    }
    free(m2e->acf_sum);
    free(m2e->acf_mean);

    for(uint8_t i=0 ; i<m2e->eqp_pow ; i++)
    {
        if(m2e->acf_partial[i] != NULL)
        {
            for(uint64_t j=0 ; j<(1 + !m2e->eqp_pow)*m2e->eqp_num ; j++)
            {
                if(m2e->acf_partial[i][j] != NULL)
                {
                    for(uint8_t k=0 ; k<m2e->acf_pow ; k++)
                    { free(m2e->acf_partial[i][j][k]); }
                }
                free(m2e->acf_partial[i][j]);
            }
        }
        free(m2e->acf_partial[i]);
    }
    free(m2e->acf_partial);

    // TODO: destroy FFT data structures
}