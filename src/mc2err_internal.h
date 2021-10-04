// mc2err internal header file
#ifndef MC2ERR_INTERNAL_H
#define MC2ERR_INTERNAL_H

// mc2err header file
#include "mc2err.h"

// expand memory footprint of mc2err data structure to accommodate index maxima
uint8_t mc2err_expand(uint64_t max_chain, uint8_t max_block, struct mc2err_data *m2e);

// perform a statistical analysis of all data inside an mc2err data structure
uint8_t mc2err_analyze(struct mc2err_data *m2e);

// subsystem for MLE
struct mc2err_likelihood_work
{
    int lwork;
    double *sum1, *sum2, *w, *work;
};
uint8_t mc2err_likelihood_begin(uint16_t width, struct mc2err_likelihood_work *work);
uint8_t mc2err_likelihood(uint16_t width, uint64_t num, double *sum1, double *sum2, double *likelihood,
                   struct mc2err_likelihood_work *work);
void mc2err_likelihood_end(struct mc2err_likelihood_work *work);

#endif