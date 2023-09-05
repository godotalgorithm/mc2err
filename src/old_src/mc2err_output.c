// mc2err header file
#include "mc2err.h"
#include "mc2err_internal.h"

// output an estimated mean vector & covariance matrix from all data inside an mc2err data structure
uint8_t mc2err_output(double *mean, double *covariance, struct mc2err_data *m2e)
{
    // perform the statistical analysis
    uint8_t status = mc2err_analyze(m2e);
    if(status) { return status; }

    // transfer output data
    uint16_t width = m2e->width;
    for(uint16_t i=0 ; i<width ; i++)
    {
        mean[i] = m2e->mean[m2e->eqp_cut][i];
        for(uint16_t j=i ; j<width ; j++)
        { covariance[j+i*width] = m2e->covariance[m2e->acf_cut+m2e->eqp_cut*m2e->num_block][i+j*(j+1)/2]; }
        for(uint16_t j=0 ; j<i ; j++)
        { covariance[j+i*width] = covariance[i+j*width]; }
    }
    return 0;
}
