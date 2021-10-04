#include "mc2err.h"

// Example 1: trivial outputs for trivial inputs (constant & normal data)

#define NUM_DATA 10000

int main(void)
{
    struct mc2err_data data;

    // constant example
    mc2err_begin(1, &data);
    double one = 1;
    for(uint64_t i=0 ; i<NUM_DATA ; i++)
    { mc2err_input(0, &one, &data); }
    double mean, covariance;
    mc2err_output(&mean, &covariance, &data);
    printf("constant data (value = 1): %e +/- %e\n", mean, covariance);
    mc2err_end(&data);

    // normal example
    srand(1);
    mc2err_begin(1, &data);
    for(uint64_t i=0 ; i<NUM_DATA ; i++)
    {
        double normal = sqrt(fabs(2.0*log((double)rand()/(double)RAND_MAX))) * cos(2.0*M_PI*rand());
        mc2err_input(0, &normal, &data);
    }
    mc2err_output(&mean, &covariance, &data);
    printf("normal data (mean = 0, variance = 1): %e +/- %e\n", mean, covariance);
    mc2err_end(&data);

    return 0;
}
