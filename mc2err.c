// mc2err: online averages, error bars, & equilibration analysis for Markov chain Monte Carlo data streams
// C99 standard compliant (C89 + permissive variable declaration + fixed-width integer types)

#include "mc2err.h"

// convenient macros for the minimum and maximum of two comparable values
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

// power-of-2 size associated with a given grid index
#define GRID_POW(grid) (((grid) == 0) ? 0 : ((grid)-1))

// floor(log2(x)) function (slow implementation)
uint8_t log2floor(uint64_t x)
{
    uint8_t num = 0;
    while(x != 0)
    { x >>= 1; num++; }
    return num;
}

// grid & subgrid indices of a power-of-2 grid w/ "num" points per subgrid for a given data "offset"
void pow2index(uint64_t offset, uint32_t num, uint8_t *grid, uint32_t *sub)
{
    *grid = log2floor(offset/num);
    *sub = (offset>>GRID_POW(*grid)) - ((*grid) ? num : 0);
}

// initialize a power-of-2 logarithmic grid
uint8_t pow2grid_initialize(uint64_t max, double precision, uint8_t *pow, uint32_t *num)
{
    if(ceil(0.5*(1.0 + 1.0/precision)) > UINT32_MAX)
    { return 1; }

    *num = (uint32_t)ceil(0.5*(1.0 + 1.0/precision));
    *pow = 1;
    while( (((uint64_t)(*num)) << *pow) < max &&
           ((((uint64_t)(*num)) << (*pow-1)) >> 63) == 0 )
    { (*pow)++; }

    return 0;
}

// initializes parameters & allocates memory for an mc2err structure
void mc2err_initialize(uint8_t data_dim, uint64_t data_max, double acf_precision, double eqp_precision, struct mc2err *m2e)
{
    // initialize the basic parameters
    m2e->data_num = 0;
    m2e->data_dim = data_dim;
    m2e->data_max = data_max;
    m2e->acf_precision = acf_precision;
    m2e->eqp_precision = eqp_precision;

    // initialize power-of-2 logarithmic grids
    if(pow2grid_initialize(m2e->data_max, m2e->acf_precision, &(m2e->acf_pow), &(m2e->acf_num)))
    { printf("ERROR: acf_precision is too small in mc2err\n"); exit(1); }
    if(pow2grid_initialize(m2e->data_max, m2e->eqp_precision, &(m2e->eqp_pow), &(m2e->eqp_num)))
    { printf("ERROR: eqp_precision is too small in mc2err\n"); exit(1); }

    // allocate memory for buffer pointers & indices
    m2e->buffer_offset = (uint64_t*)malloc(sizeof(uint64_t)*(m2e->acf_pow-1));
    m2e->buffer_transfer = (uint64_t*)malloc(sizeof(uint64_t)*(m2e->acf_pow-1));
    m2e->buffer_process = (uint64_t*)malloc(sizeof(uint64_t)*(m2e->acf_pow-1));
    m2e->buffer = (double**)malloc(sizeof(double*)*(m2e->acf_pow-1));
    for(uint8_t i=0 ; i<m2e->acf_pow-1 ; i++)
    { m2e->buffer[i] = NULL; m2e->buffer_offset[i] = m2e->buffer_transfer[i] = m2e->buffer_process[i] = 0; }

    // allocate memory for data & acf storage
    m2e->data_partial = (double**)malloc(sizeof(double*)*m2e->eqp_pow);
    m2e->data_mean = (double**)malloc(sizeof(double*)*m2e->eqp_pow);
    m2e->data_error = (double**)malloc(sizeof(double*)*m2e->eqp_pow);
    m2e->acf_partial = (double****)malloc(sizeof(double***)*m2e->eqp_pow);
    m2e->acf_mean = (double**)malloc(sizeof(double*)*m2e->acf_pow);
    for(uint8_t i=0 ; i<m2e->eqp_pow ; i++)
    { m2e->data_partial[i] = m2e->data_mean[i] = m2e->data_error[i] = m2e->acf_mean[i] = NULL; m2e->acf_partial[i] = NULL; }
    for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
    { m2e->acf_mean[i] = NULL; }

    // allocate the 1st data buffer
    m2e->buffer[0] = (double*)malloc(sizeof(double)*4*m2e->acf_num*m2e->data_dim);

    // TEMPORARY INITIALIZATION TO SUPPORT TEMPORARY INPUT/OUTPUT STUB
    for(uint8_t i=0 ; i<2*m2e->data_dim ; i++)
    { m2e->buffer[0][i] = 0.0; }

    // TODO: setup FFT data structures
}

// deallocate all of the memory used by mc2err
void mc2err_finalize(struct mc2err *m2e)
{
    free(m2e->buffer_offset);
    free(m2e->buffer_transfer);
    free(m2e->buffer_process);

    for(uint8_t i=0 ; i<m2e->acf_pow-1 ; i++)
    { free(m2e->buffer[i]); }
    free(m2e->buffer);

    for(uint8_t i=0 ; i<m2e->eqp_pow ; i++)
    {
        free(m2e->data_partial[i]);
        free(m2e->data_mean[i]);
        free(m2e->data_error[i]);
    }
    free(m2e->data_partial);
    free(m2e->data_mean);
    free(m2e->data_error);

    for(uint8_t i=0 ; i<m2e->acf_pow ; i++)
    { free(m2e->acf_mean[i]); }
    free(m2e->acf_mean);

    for(uint8_t i=0 ; i<m2e->eqp_pow ; i++)
    {
        if(m2e->acf_partial[i] != NULL)
        {
            for(uint32_t j=0 ; j<m2e->eqp_num ; j++)
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

// check if memory is allocated for processing a buffer & allocate memory as necessary
void mc2err_touch(uint8_t buffer_index, uint8_t eqp_grid, uint32_t eqp_sub, struct mc2err *m2e)
{
    // 1st-level allocation
    if(m2e->acf_partial[eqp_grid] == NULL)
    { m2e->acf_partial[eqp_grid] = (double***)malloc(sizeof(double**)*m2e->eqp_num); }

    // 2nd-level allocation
    if(m2e->acf_partial[eqp_grid][eqp_sub] == NULL)
    { m2e->acf_partial[eqp_grid][eqp_sub] = (double**)malloc(sizeof(double*)*m2e->acf_pow); }

    // memory that is only relevant to the 1st buffer
    if(buffer_index == 0)
    {
        if(m2e->data_partial[eqp_grid] == NULL)
        {
            m2e->data_partial[eqp_grid] = (double*)malloc(sizeof(double)*m2e->eqp_num*m2e->data_dim);
            for(uint64_t i=0 ; i<m2e->eqp_num*m2e->data_dim ; i++)
            { m2e->data_partial[eqp_grid][i] = 0.0; }
        }

        if(m2e->acf_partial[eqp_grid][eqp_sub][0] == NULL)
        {
            m2e->acf_partial[eqp_grid][eqp_sub][0] = (double*)malloc(sizeof(double)*m2e->acf_num*m2e->data_dim*m2e->data_dim);
            for(uint64_t i=0 ; i<m2e->acf_num*m2e->data_dim*m2e->data_dim ; i++)
            { m2e->acf_partial[eqp_grid][eqp_sub][0][i] = 0.0; }
        }
    }

    // memory that is relevant to all buffers
    if(m2e->acf_partial[eqp_grid][eqp_sub][buffer_index+1] == NULL)
    {
        m2e->acf_partial[eqp_grid][eqp_sub][buffer_index+1] = (double*)malloc(sizeof(double)*m2e->acf_num*m2e->data_dim*m2e->data_dim);

        for(uint64_t i=0 ; i<m2e->acf_num*m2e->data_dim*m2e->data_dim ; i++)
        { m2e->acf_partial[eqp_grid][eqp_sub][buffer_index+1][i] = 0.0; }
    }
}

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

// transfer data from one buffer to the next buffer
void mc2err_transfer(uint8_t buffer_index, struct mc2err *m2e)
{
    // nothing to do if we are in the last buffer
    if(buffer_index >= m2e->acf_pow-2)
    { return; }

    // identify the beginning & ending locations for the data transfer between buffers
    uint64_t offset = (m2e->buffer_transfer[buffer_index+1] > 2*m2e->acf_num) ? 2*m2e->acf_num : 0;
    uint64_t begin = m2e->buffer_transfer[buffer_index+1] - offset;
    uint64_t end = m2e->buffer_transfer[buffer_index]/2;

    // stop if there is no new data to propagate to the next buffer
    if(begin == end)
    { return; }

    // allocate the next buffer if necessary
    if(m2e->buffer[buffer_index+1] == NULL)
    { m2e->buffer[buffer_index+1] = (double*)malloc(sizeof(double)*4*m2e->acf_num*m2e->data_dim); }

    // write the block data to the next buffer
    for(uint64_t i=begin ; i<end ; i++)
    for(uint8_t j=0 ; j<m2e->data_dim ; j++)
    {
        m2e->buffer[buffer_index+1][j+(offset+i)*m2e->data_dim] = 0.5*(m2e->buffer[buffer_index][j+2*i*m2e->data_dim] +
                                                                       m2e->buffer[buffer_index][j+(2*i+1)*m2e->data_dim]);
    }

    // update the write position
    m2e->buffer_transfer[buffer_index+1] = offset + end;
}

// shift filled buffer & process data that will be erased by the shift
void mc2err_process(uint8_t buffer_index, struct mc2err *m2e)
{
    // only process buffers with complete autocorrelation function pairs
    if(m2e->buffer_transfer[buffer_index] <= 2*m2e->acf_num)
    { return; }

    // loop over the data in the 1st half of the buffer
    while(m2e->buffer_process[buffer_index] < m2e->buffer_transfer[buffer_index] - 2*m2e->acf_num)
    {
        // identify the active equilibration point & block offset of its last data point
        uint8_t eqp_grid;
        uint32_t eqp_sub;
        pow2index(m2e->buffer_offset[buffer_index]<<buffer_index, m2e->eqp_num, &eqp_grid, &eqp_sub);
        uint64_t offset_max = (m2e->eqp_num*(1 + eqp_sub + (1<<GRID_POW(eqp_grid))) - 1)>>buffer_index;

        // identify the amount of data to be processed for the active equilibration point
        uint32_t process_num = MIN(2*m2e->acf_num - m2e->buffer_process[buffer_index], 1+offset_max-m2e->buffer_offset[buffer_index]);

        // on-the-fly memory allocation of partial analysis buffers
        mc2err_touch(buffer_index, eqp_grid, eqp_sub, m2e);

        // update data_partial & acf_partial[0] for the 1st buffer
        if(buffer_index == 0)
        {
            // update data_partial
            for(uint32_t i=0 ; i<process_num ; i++)
            for(uint8_t j=0 ; j<m2e->data_dim ; j++)
            { m2e->data_partial[eqp_grid][j+(i+eqp_sub)*m2e->data_dim] += m2e->buffer[0][j+(i+m2e->buffer_process[0])*m2e->data_dim]; }

            // update acf_partial[0]
            convolution(m2e->data_dim,
                        process_num,
                        process_num + m2e->acf_num,
                        m2e->acf_num,
                        m2e->buffer[0] + m2e->buffer_process[0]*m2e->data_dim,
                        m2e->buffer[0] + m2e->buffer_process[0]*m2e->data_dim,
                        m2e->acf_partial[eqp_grid][eqp_sub][0]);
        }

        // update acf_partial[buffer_index+1]
        convolution(m2e->data_dim,
                    process_num,
                    process_num + m2e->acf_num,
                    m2e->acf_num,
                    m2e->buffer[buffer_index] + m2e->buffer_process[buffer_index]*m2e->data_dim,
                    m2e->buffer[buffer_index] + (m2e->buffer_process[buffer_index] + m2e->acf_num)*m2e->data_dim,
                    m2e->acf_partial[eqp_grid][eqp_sub][buffer_index+1]);

        // update the processing indices
        m2e->buffer_offset[buffer_index] += process_num;
        m2e->buffer_process[buffer_index] += process_num;
    }

    // shift a buffer when half of the data has been processed
    if(m2e->buffer_process[buffer_index] == 2*m2e->acf_num)
    {
        // copy 2nd half of the buffer to the 1st half
        memcpy(m2e->buffer[buffer_index], m2e->buffer[buffer_index]+2*m2e->acf_num*m2e->data_dim, 2*m2e->acf_num*m2e->data_dim*sizeof(double));

        // shift buffer indices to coincide w/ shifted buffer
        m2e->buffer_transfer[buffer_index] -= 2*m2e->acf_num;
        m2e->buffer_process[buffer_index] -= 2*m2e->acf_num;
    }
}

// input new data into the mc2err structure & return # of data points that couldn't be processed
uint64_t mc2err_input(uint64_t new_data_num, double *new_data, struct mc2err *m2e)
{
    // FAKE INPUT STUB FOR PRE-ALPHA API TESTING (simple sample mean & variance)
    for(uint8_t i=0 ; i<m2e->data_dim ; i++)
    {
        for(uint64_t j=0 ; j<new_data_num ; j++)
        {
            m2e->buffer[0][i] += new_data[i+j*m2e->data_dim];
            m2e->buffer[0][i+m2e->data_dim] += new_data[i+j*m2e->data_dim]*new_data[i+j*m2e->data_dim];
        }
    }
    m2e->data_num += new_data_num;
    return 0;

    // check if the maximum number of samples will be exceeded
    uint64_t excess_data_num = 0;
    if(new_data_num > m2e->data_max - m2e->data_num)
    {
        excess_data_num = new_data_num - (m2e->data_max - m2e->data_num);
        new_data_num = m2e->data_max - m2e->data_num;
    }
    m2e->data_num += new_data_num;

    // copy data into the main buffer
    while(new_data_num != 0)
    {
        // identify the number of data points to copy,
        // either copy all remaining data or fill up the primary data buffer
        uint32_t copy_num = MIN(new_data_num, 4*m2e->acf_num-m2e->buffer_transfer[0]);

        // use memcpy to copy data efficiently
        memcpy(m2e->buffer[0]+m2e->data_dim*m2e->buffer_transfer[0], new_data, copy_num*m2e->data_dim*sizeof(double));

        // update indices & pointers
        m2e->buffer_transfer[0] += copy_num;
        new_data_num -= copy_num;
        new_data += copy_num;

        // loop over filled data buffers
        for(uint8_t i=0 ; i<m2e->acf_pow-2 && m2e->buffer_transfer[i]==4*m2e->acf_num ; i++)
        {
            // transfer data between buffers
            mc2err_transfer(i, m2e);

            // process all data from the 1st half of the buffer
            mc2err_process(i, m2e);
        }
    }
    return excess_data_num;
}

// analyze all available data in mc2err
uint8_t mc2err_output(double *mean, double *error, double *size, struct mc2err *m2e)
{
    // transfer all available data between block data buffers
    for(uint8_t i=0 ; i<m2e->acf_pow-1 ; i++)
    {
        // transfer complete blocked data to the next buffer
        mc2err_transfer(i, m2e);

        // process all data with complete pair information for the autocorrelation function
        mc2err_process(i, m2e);
    }

    // FAKE OUTPUT STUB FOR PRE-ALPHA API TESTING (simple sample mean & variance)
    for(uint8_t i=0 ; i<m2e->data_dim ; i++)
    {
        size[i] = (double)m2e->data_num;
        mean[i] = m2e->buffer[0][i]/m2e->data_num;
        error[i] = sqrt(m2e->buffer[0][i+m2e->data_dim]/m2e->data_num - mean[i]*mean[i])/sqrt(m2e->data_num);
    }
    return 0;

    // calculate maximum acf & eqp indices

    // expand online analysis buffers as necessary

    // clear data analysis buffers
 
    // clear autocorrelation analysis buffers

    // countdown through active equilibration points
    {
        // add in residual data

        // add in partial summations

        // statistical estimate of the autocorrelation truncation & error bar
    }

    // statistical estimate of the equilibration point

    // calculate the autocorrelation function for the estimated equilibration point
}
