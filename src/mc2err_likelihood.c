// mc2err header file
#include "mc2err.h"
#include "mc2err_internal.h"

// local definition of pi to guarantee availability
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// numerical regularization parameters for MLE
#define MC2ERR_RELATIVE_MIN 1e-8
#define MC2ERR_ABSOLUTE_MIN 1e-308

uint8_t mc2err_likelihood_begin(uint16_t width, struct mc2err_likelihood_work *work)
{
    work->sum1 = (double*)malloc(sizeof(double)*width);
    if(work->sum1 == NULL) { return 1; }
    work->sum2 = (double*)malloc(sizeof(double)*width*width);
    if(work->sum2 == NULL) { return 1; }

    char jobz = 'V', uplo = 'U';
    int size = width, info;
    double work0;
    work->w = (double*)malloc(sizeof(double)*size);
    if(work->w == NULL) { return 1; }
    work->lwork = -1;
    MC2ERR_LAPACK_DSYEV(&jobz, &uplo, &size, NULL, &size, work->w, &work0, &work->lwork, &info);
    if(info != 0) { return 5; }
    work->lwork = (int)work0;
    work->work = (double*)malloc(sizeof(double)*work->lwork);
    if(work->work == NULL) { return 1; }
    return 0;
}

uint8_t mc2err_likelihood(uint16_t width, uint64_t num, double *sum1, double *sum2, double *likelihood,
                          struct mc2err_likelihood_work *work)
{
    // trivial result without data
    if(num == 0)
    {
        *likelihood = 0.0;
        return 0;
    }

    // MLE offset
    *likelihood = (double)width*log(2.0*M_PI);

    // rescale local copies of sum1 & sum2, account for scale factor contributions to MLE
    for(uint16_t i=0 ; i<width ; i++)
    {
        double wt = sum2[i+i*(i+1)/2];
        if(wt < MC2ERR_ABSOLUTE_MIN*(double)num)
        { wt = MC2ERR_ABSOLUTE_MIN*(double)num; }

        *likelihood += log(wt/(double)num);

        wt = 1.0/sqrt(wt);
        work->sum1[i] = wt*sum1[i];
        for(uint16_t j=0 ; j<i ; j++)
        { work->sum2[j+i*width] = wt*sum2[j+i*(i+1)/2]; }
        work->sum2[i+i*width] = 1.0;
        for(uint16_t j=i+1 ; j<width ; j++)
        { work->sum2[i+j*width] = wt*sum2[i+j*(j+1)/2]; }
    }

    // diagonalize sum2 copy
    char jobz = 'N', uplo = 'U';
    if(sum1 != NULL) { jobz = 'V'; }
    int size = width, info;
    MC2ERR_LAPACK_DSYEV(&jobz, &uplo, &size, work->sum2, &size, work->w, work->work, &work->lwork, &info);
    if(info != 0) { return 5; }

    // calculate minimum variance
    double min_variance = work->w[width-1];
    for(uint16_t i=0 ; i<width ; i++)
    { min_variance += work->sum1[i]*work->sum1[i]; }
    min_variance *= MC2ERR_RELATIVE_MIN;

    // transform sum1 copy
    char trans = 'T';
    int inc = 1;
    double one = 1.0, zero = 0.0;
    MC2ERR_BLAS_DGEMV(&trans, &size, &size, &one, work->sum2, &size, work->sum1, &inc, &zero, work->work, &inc);

    // add normalized & regularized variance contributions to MLE
    for(uint16_t i=0 ; i<width ; i++)
    {
        double variance = work->w[i] - work->work[i]*work->work[i];
        if(variance > min_variance)
        { *likelihood += log(variance) + 1.0; }
        else
        { *likelihood += log(min_variance) + variance/min_variance; }
    }

    // rescale likelihood by the number of data points
    *likelihood *= (double)num;
    return 0;
}

void mc2err_likelihood_end(struct mc2err_likelihood_work *work)
{
    free(work->w);
    free(work->work);
}