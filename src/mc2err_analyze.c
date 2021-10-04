// mc2err header files
#include "mc2err.h"
#include "mc2err_internal.h"

// perform a statistical analysis of all data inside an mc2err data structure
uint8_t mc2err_analyze(struct mc2err_data *m2e)
{
    // report an error if there was no data to analyze
    if(m2e->num_input == 0)
    { return 6; }

    // exit if there is no new data to analyze
    if(m2e->num_input == m2e->num_output)
    { return 0; }

    // local loop size variables
    uint8_t num_block = m2e->num_block;
    uint16_t num_block2t = m2e->num_block*(m2e->num_block+1)/2;
    uint16_t width = m2e->width;
    uint32_t width2t = m2e->width*(m2e->width+1)/2;

    // allocate temporary workspaces
    struct mc2err_likelihood_work work;
    uint8_t status = mc2err_likelihood_begin(width, &work);
    if(status) { return status; }
    uint64_t accumulate_num;
    double *accumulate1 = (double*)malloc(sizeof(double)*width);
    double *accumulate2 = (double*)malloc(sizeof(double)*width2t);

    // initialize AIC likelihood values & assign AIC penalty values
    for(uint8_t i=0 ; i<num_block ; i++)
    for(uint8_t j=0 ; j<num_block ; j++)
    {
        m2e->likelihood[j+i*num_block] = 0.0;
        m2e->penalty[j+i*num_block] = (double)(width*(1+(i+1)*j) + width2t*(i+1)*(j+1));
    }

    // calculate Haar wavelet likelihoods in pre-EQP intervals
    for(uint8_t i=1 ; i+2<num_block ; i++)
    for(uint8_t j=0 ; j<i ; j++)
    {
        double mle;
        status = mc2err_likelihood(width, m2e->haar_num[j+i*(i+1)/2],
            m2e->haar_sum1[j+i*(i+1)/2], m2e->haar_sum2[j+i*(i+1)/2], &mle, &work);
        if(status) { return status; }

        m2e->likelihood[j+1+(i+2)*num_block] += mle;
    }
    for(uint8_t i=0 ; i+1<num_block ; i++)
    {
        double mle;
        status = mc2err_likelihood(width, m2e->haar_num[i+i*(i+1)/2],
            m2e->haar_sum1[i+i*(i+1)/2], m2e->haar_sum2[i+i*(i+1)/2], &mle, &work);
        if(status) { return status; }

        m2e->likelihood[i+1+num_block] += mle;
    }

    // accumulate Haar wavelet likelihoods in pre-EQP intervals over EQP
    for(uint8_t i=2 ; i<num_block ; i++)
    for(uint8_t j=1 ; j<num_block ; j++)
    { m2e->likelihood[j+i*num_block] += m2e->likelihood[j+(i-1)*num_block]; }

    // calculate Haar wavelet likelihoods of accumulated post-EQP data
    for(uint8_t i=0 ; i+1<num_block ; i++)
    {
        accumulate_num = 0;
        for(uint16_t j=0 ; j<width ; j++)
        { accumulate1[j] = 0.0; }
        for(uint32_t j=0 ; j<width2t ; j++)
        { accumulate2[j] = 0.0; }

        double mle;
        for(uint8_t j=num_block ; j>i+2 ; j--)
        {
            accumulate_num += m2e->haar_num[i+(j-2)*(j-1)/2];
            for(uint32_t k=0 ; k<width2t ; k++)
            { accumulate2[k] += m2e->haar_sum2[i+(j-2)*(j-1)/2][k]; }

            // zero-mean model
            status = mc2err_likelihood(width, accumulate_num, accumulate1, accumulate2, &mle, &work);
            if(status) { return status; }

            m2e->likelihood[i+1+(j-1)*num_block] += mle;
        }
        for(uint8_t j=1 ; j<=i+1 ; j++)
        { m2e->likelihood[i+1+j*num_block] += mle; }

        accumulate_num += m2e->haar_num[i];
        for(uint32_t k=0 ; k<width2t ; k++)
        { accumulate2[k] += m2e->haar_sum2[i][k]; }

        // zero-mean model
        status = mc2err_likelihood(width, accumulate_num, accumulate1, accumulate2, &mle, &work);
        if(status) { return status; }

        m2e->likelihood[i+1] += mle;
    }

    // accumulate Haar wavelet likelihoods over ACF
    for(uint8_t i=0 ; i<num_block ; i++)
    for(uint8_t j=2 ; j<num_block ; j++)
    { m2e->likelihood[j+i*num_block] += m2e->likelihood[j-1+i*num_block]; }

    // calculate Haar scaling likelihoods in pre-EQP intervals
    for(uint8_t i=1 ; i+1<num_block ; i++)
    for(uint8_t j=0 ; j<i ; j++)
    {
        double mle;
        status = mc2err_likelihood(width, m2e->data_num[j+i*(i+1)/2],
            m2e->data_sum1[j+i*(i+1)/2], m2e->data_sum2[j+i*(i+1)/2], &mle, &work);
        if(status) { return status; }

        m2e->likelihood[j+(i+1)*num_block] += mle;
    }
    for(uint8_t i=0 ; i<num_block ; i++)
    {
        double mle;
        status = mc2err_likelihood(width, m2e->data_num[i+i*(i+1)/2],
            m2e->data_sum1[i+i*(i+1)/2], m2e->data_sum2[i+i*(i+1)/2], &mle, &work);
        if(status) { return status; }

        m2e->likelihood[i+num_block] += mle;
    }

    // calculate Haar scaling likelihoods of accumulated post-EQP data
    for(uint8_t i=0 ; i<num_block ; i++)
    {
        accumulate_num = 0;
        for(uint16_t j=0 ; j<width ; j++)
        { accumulate1[j] = 0.0; }
        for(uint32_t j=0 ; j<width2t ; j++)
        { accumulate2[j] = 0.0; }

        double mle;
        for(uint8_t j=num_block ; j>i+1 ; j--)
        {
            accumulate_num += m2e->data_num[i+(j-1)*j/2];
            for(uint16_t k=0 ; k<width ; k++)
            { accumulate1[k] += m2e->data_sum1[i+(j-1)*j/2][k]; }
            for(uint32_t k=0 ; k<width2t ; k++)
            { accumulate2[k] += m2e->data_sum2[i+(j-1)*j/2][k]; }

            status = mc2err_likelihood(width, accumulate_num, accumulate1, accumulate2, &mle, &work);
            if(status) { return status; }

            m2e->likelihood[i+(j-1)*num_block] += mle;
        }
        for(uint8_t j=1 ; j<=i ; j++)
        { m2e->likelihood[i+j*num_block] += mle; }

        accumulate_num += m2e->data_num[i];
        for(uint16_t k=0 ; k<width ; k++)
        { accumulate1[k] += m2e->data_sum1[i][k]; }
        for(uint32_t k=0 ; k<width2t ; k++)
        { accumulate2[k] += m2e->data_sum2[i][k]; }

        status = mc2err_likelihood(width, accumulate_num, accumulate1, accumulate2, &mle, &work);
        if(status) { return status; }

        m2e->likelihood[i] += mle;
    }

    // accumulate data into mean & non-lagged covariance estimates over EQP
    if(num_block > 0)
    {
        for(uint16_t j=0 ; j<width ; j++)
        { m2e->mean[num_block-1][j] = m2e->data_sum1[(num_block-1)*num_block/2][j]; }
        for(uint32_t j=0 ; j<width2t ; j++)
        { m2e->covariance[(num_block-1)*num_block][j] = m2e->data_sum2[(num_block-1)*num_block/2][j]; }
    }
    for(uint8_t i=num_block ; i>1 ; i--)
    {
        for(uint16_t j=0 ; j<width ; j++)
        { m2e->mean[i-2][j] = m2e->mean[i-1][j] + m2e->data_sum1[(i-2)*(i-1)/2][j]; }
        for(uint32_t j=0 ; j<width2t ; j++)
        { m2e->covariance[(i-2)*num_block][j] = m2e->covariance[(i-1)*num_block][j]
            + m2e->data_sum2[(i-2)*(i-1)/2][j]; }
    }

    // weight the mean & non-lagged covariance estimates
    accumulate_num = 0;
    for(uint8_t i=num_block ; i>1 ; i--)
    {
        accumulate_num += m2e->data_num[(i-1)*i/2];
        if(accumulate_num)
        {
            double wt = 1.0/(double)accumulate_num;
            for(uint16_t j=0 ; j<width ; j++)
            { m2e->mean[i-1][j] *= wt; }
            for(uint32_t j=0 ; j<width2t ; j++)
            { m2e->covariance[(i-1)*num_block][j] *= wt; }
        }
    }

    // accumulate data into lagged covariance estimates over EQP
    for(uint8_t i=1 ; i<num_block ; i++)
    {
        for(uint32_t j=0 ; j<width2t ; j++)
        { m2e->covariance[i+(num_block-1)*num_block][j] = m2e->lag_sum2[i-1+(num_block-1)*num_block/2][j]; }

        for(uint8_t j=num_block ; j>1 ; j--)
        {
            for(uint32_t k=0 ; k<width2t ; k++)
            { m2e->covariance[i+(j-2)*num_block][k] = m2e->covariance[i+(j-1)*num_block][k]
                + m2e->lag_sum2[i-1+(j-2)*(j-1)/2][k]; }
        }
    }

    // weight the lagged covariance estimates
    for(uint8_t i=1 ; i<num_block ; i++)
    {
        accumulate_num = 0;
        for(uint8_t j=num_block ; j>1 ; j--)
        {
            accumulate_num += m2e->lag_num[i-1+(j-1)*j/2];
            if(accumulate_num)
            {
                double wt = 1.0/(double)accumulate_num;
                for(uint32_t k=0 ; k<width2t ; k++)
                { m2e->covariance[i+(j-1)*num_block][k] *= wt; }
            }
        }
    }

    // accumulate data into covariance estimates over ACF
    for(uint8_t i=0 ; i<num_block ; i++)
    {
        for(uint32_t j=0 ; j<width2t ; j++)
        { accumulate2[j] = m2e->covariance[i*num_block][j]; }

        for(uint8_t j=1 ; j<num_block ; j++)
        {
            // accumulate lagged block covariances
            for(uint32_t k=0 ; k<width2t ; k++)
            { accumulate2[k] += m2e->covariance[j+i*num_block][k]; }

            // double-count the final lagged block
            for(uint32_t k=0 ; k<width2t ; k++)
            { m2e->covariance[j+i*num_block][k] += accumulate2[k]; }
        }
    }

    // include the mean-squared terms in the covariance matrices
    for(uint8_t i=0 ; i<num_block ; i++)
    for(uint8_t j=0 ; j<num_block ; j++)
    for(uint16_t k=0 ; k<width ; k++)
    for(uint16_t l=0 ; l<=k ; l++)
    { m2e->covariance[i+j*num_block][l+k*(k+1)/2] -= m2e->mean[j][k]*m2e->mean[j][l]; }

    // determine the most likely statistical estimate
    m2e->acf_cut = 0;
    m2e->eqp_cut = 0;
    double min_aic = m2e->likelihood[0] + m2e->penalty[0];
    for(uint8_t i=0 ; i<num_block ; i++)
    for(uint8_t j=0 ; j<num_block ; j++)
    {
        double aic = m2e->likelihood[j+i*num_block] + m2e->penalty[j+i*num_block];
        if(aic < min_aic)
        {
            m2e->acf_cut = j;
            m2e->eqp_cut = i;
            min_aic = aic;
        }
    }

    // free temporary workspaces
    mc2err_likelihood_end(&work);
    free(accumulate1);
    free(accumulate2);

    // mark all data as analyzed
    m2e->num_output = m2e->num_input;
    return 0;
}

void dpotrf_(char*, int*, double*, int*, int*);

// old, unregularized likelihood implementation (REMOVE AFTER USING IT FOR VERIFICATION TESTS)
uint8_t mc2err_likelihood0(uint16_t width, uint64_t num, double *covariance, double *likelihood)
{
    char uplo = 'U';
    int size = width;
    int info;
    dpotrf_(&uplo, &size, covariance, &size, &info);
    if(info != 0) { return 5; }

    *likelihood = (double)width*(1.0 + log(2.0*M_PI));
    for(uint16_t i=0 ; i<width ; i++)
    { *likelihood += 2.0*log(covariance[i+i*size]); }
    *likelihood *= (double)num;

    return 0;
}