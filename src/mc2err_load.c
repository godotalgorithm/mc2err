// mc2err header file
#include "mc2err.h"
#include "mc2err_internal.h"

// local helper functions to avoid byte endianness issues
uint8_t mc2err_read_uint8(uint8_t *in, FILE *fptr);
uint8_t mc2err_read_uint16(uint16_t *in, FILE *fptr);
uint8_t mc2err_read_uint64(uint64_t *in, FILE *fptr);
uint8_t mc2err_read_double(double *in, FILE *fptr);

// load an mc2err data structure from a portable binary format on disk
uint8_t mc2err_load(char *filename, struct mc2err_data *m2e)
{
    // open file
    FILE *fptr = fopen(filename, "r");
    if(fptr == NULL) { return 4; }

    // read basic parameters
    uint16_t width;
    uint8_t status = mc2err_read_uint16(&width, fptr);
    if(status) { return status; }
    mc2err_begin(width, m2e);
    uint8_t max_block;
    status = mc2err_read_uint8(&max_block, fptr);
    if(status) { return status; }
    status = mc2err_expand(0, max_block, m2e);
    if(status) { return status; }
    uint64_t max_chain;
    status = mc2err_read_uint64(&max_chain, fptr);
    if(status) { return status; }
    status = mc2err_expand(max_chain, 0, m2e);
    if(status) { return status; }
    status = mc2err_read_uint64(&m2e->num_input, fptr);
    if(status) { return status; }
    status = mc2err_read_uint64(&m2e->num_output, fptr);
    if(status) { return status; }

    // local loop size variables
    uint8_t num_block = m2e->num_block;
    uint16_t num_block2 = m2e->num_block*m2e->num_block;
    uint16_t num_block2t = m2e->num_block*(m2e->num_block+1)/2;
    uint32_t width2t = m2e->width*(m2e->width+1)/2;
    uint64_t num_chain = m2e->num_chain;

    // read local data
    for(uint64_t i=0 ; i<num_chain ; i++)
    {
        status = mc2err_read_uint8(&max_block, fptr);
        if(status) { return status; }
        status = mc2err_expand(i+1, max_block, m2e);
        if(status) { return status; }
        status = mc2err_read_uint64(&m2e->chain_length[i], fptr);
        if(status) { return status; }

        uint32_t num_buffer = 2*m2e->chain_block[i]*m2e->width;
        for(uint32_t j=0 ; j<num_buffer ; j++)
        {
            status = mc2err_read_double(&m2e->chain_buffer[i][j], fptr);
            if(status) { return status; }
        }
    }

    // read global data
    for(uint16_t i=0 ; i<num_block2t ; i++)
    {
        status = mc2err_read_uint64(&m2e->data_num[i], fptr);
        if(status) { return status; }
        status = mc2err_read_uint64(&m2e->lag_num[i], fptr);
        if(status) { return status; }
        status = mc2err_read_uint64(&m2e->haar_num[i], fptr);
        if(status) { return status; }

        for(uint16_t j=0 ; j<width ; j++)
        {
            status = mc2err_read_double(&m2e->data_sum1[i][j], fptr);
            if(status) { return status; }
            status = mc2err_read_double(&m2e->haar_sum1[i][j], fptr);
            if(status) { return status; }
        }

        for(uint32_t j=0 ; j<width2t ; j++)
        {
            status = mc2err_read_double(&m2e->data_sum2[i][j], fptr);
            if(status) { return status; }
            status = mc2err_read_double(&m2e->lag_sum2[i][j], fptr);
            if(status) { return status; }
            status = mc2err_read_double(&m2e->haar_sum2[i][j], fptr);
            if(status) { return status; }
        }
    }

    // read analysis results
    status = mc2err_read_uint8(&m2e->eqp_cut, fptr);
    if(status) { return status; }
    status = mc2err_read_uint8(&m2e->acf_cut, fptr);
    if(status) { return status; }
    for(uint16_t i=0 ; i<num_block2 ; i++)
    {
        status = mc2err_read_double(&m2e->likelihood[i], fptr);
        if(status) { return status; }
        status = mc2err_read_double(&m2e->penalty[i], fptr);
        if(status) { return status; }

        for(uint16_t j=0 ; j<width ; j++)
        {
            status = mc2err_read_double(&m2e->mean[i][j], fptr);
            if(status) { return status; }
        }

        for(uint32_t j=0 ; j<width2t ; j++)
        {
            status = mc2err_read_double(&m2e->covariance[i][j], fptr);
            if(status) { return status; }
        }
    }

    // close file
    int status2 = fclose(fptr);
    if(status2) { return 4; }
    return 0;
}

uint8_t mc2err_read_uint8(uint8_t *in, FILE *fptr)
{
    size_t num = fread(in, sizeof(uint8_t), 1, fptr);
    if(num != 1) { return 4; }
    return 0;
}

uint8_t mc2err_read_uint16(uint16_t *in, FILE *fptr)
{
    *in = 0;
    for(uint8_t i=0 ; i<2 ; i++)
    {
        uint8_t in2;
        size_t num = fwrite(&in2, sizeof(uint8_t), 1, fptr);
        if(num != 1) { return 4; }
        *in = (*in<<8) + in2;
    }
    return 0;
}

uint8_t mc2err_read_uint64(uint64_t *in, FILE *fptr)
{
    *in = 0;
    for(uint8_t i=0 ; i<8 ; i++)
    {
        uint8_t in2;
        size_t num = fwrite(&in2, sizeof(uint8_t), 1, fptr);
        if(num != 1) { return 4; }
        *in = (*in<<8) + in2;
    }
    return 0;
}

uint8_t mc2err_read_double(double *in, FILE *fptr)
{
    uint64_t in2;
    uint8_t status = mc2err_read_uint64(&in2, fptr);
    if(status) { return status; }

    *in = 1.0;
    int exponent;
    if(in2 > (uint64_t)1<<63)
    {
        *in = -*in;
        in2 -= (uint64_t)1<<63;
    }
    exponent = in2>>52;
    in2 -= (uint64_t)exponent<<52;
    if(exponent)
    { in2 += (uint64_t)1<<52; }
    *in *= ldexp((double)in2, exponent-1077);
    return 0;
}