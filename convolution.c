// mc2err library for streaming MCMC error analysis (C99 standard, MIT license)
#include "mc2err.h"
#include "mc2err_internal.h"

// vectorial convolution of 2 input arrays added to an output array
void convolution(uint8_t dim, uint64_t num_input1, uint64_t num_input2, uint64_t num_output, double *input1, double *input2, double *output)
{
    // simple, slow convolution
    for(uint64_t i=0 ; i<num_output ; i++)
    {
        uint64_t max_sum = MIN(num_input1, num_input2-i);
        for(uint64_t j=0 ; j<max_sum ; j++)
        for(uint8_t k=0 ; k<dim ; k++)
        for(uint8_t l=0 ; l<dim ; l++)
        { output[(i*dim+k)*dim+l] += input1[j*dim+k]*input2[(i+j)*dim+l]; }
    }

    // TODO: fast vectorial convolution using radix-2 Cooley-Tukey FFT
}