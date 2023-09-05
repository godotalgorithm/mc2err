// mc2err header file
#include "mc2err.h"

// local helper functions to avoid byte endianness issues
uint8_t mc2err_write_uint8(uint8_t *out, FILE *fptr);
uint8_t mc2err_write_uint16(uint16_t *out, FILE *fptr);
uint8_t mc2err_write_uint64(uint64_t *out, FILE *fptr);
uint8_t mc2err_write_double(double *out, FILE *fptr);

// save an mc2err data structure to disk in a portable binary format
uint8_t mc2err_save(char *filename, struct mc2err_data *m2e)
{
    // open file
    FILE *fptr = fopen(filename, "w");
    if(fptr == NULL) { return 4; }

    // write basic parameters
    uint8_t status = mc2err_write_uint16(&m2e->width, fptr);
    if(status) { return status; }
    status = mc2err_write_uint8(&m2e->num_block, fptr);
    if(status) { return status; }
    status = mc2err_write_uint64(&m2e->num_chain, fptr);
    if(status) { return status; }
    status = mc2err_write_uint64(&m2e->num_input, fptr);
    if(status) { return status; }
    status = mc2err_write_uint64(&m2e->num_output, fptr);
    if(status) { return status; }

    // local loop size variables
    uint8_t num_block = m2e->num_block;
    uint16_t num_block2 = m2e->num_block*m2e->num_block;
    uint16_t num_block2t = m2e->num_block*(m2e->num_block+1)/2;
    uint16_t width = m2e->width;
    uint32_t width2t = m2e->width*(m2e->width+1)/2;
    uint64_t num_chain = m2e->num_chain;

    // write local data
    for(uint64_t i=0 ; i<num_chain ; i++)
    {
        status = mc2err_write_uint8(&m2e->chain_block[i], fptr);
        if(status) { return status; }
        status = mc2err_write_uint64(&m2e->chain_length[i], fptr);
        if(status) { return status; }

        uint32_t num_buffer = 2*m2e->chain_block[i]*m2e->width;
        for(uint32_t j=0 ; j<num_buffer ; j++)
        {
            status = mc2err_write_double(&m2e->chain_buffer[i][j], fptr);
            if(status) { return status; }
        }
    }

    // write global data
    for(uint16_t i=0 ; i<num_block2t ; i++)
    {
        status = mc2err_write_uint64(&m2e->data_num[i], fptr);
        if(status) { return status; }
        status = mc2err_write_uint64(&m2e->lag_num[i], fptr);
        if(status) { return status; }
        status = mc2err_write_uint64(&m2e->haar_num[i], fptr);
        if(status) { return status; }

        for(uint16_t j=0 ; j<width ; j++)
        {
            status = mc2err_write_double(&m2e->data_sum1[i][j], fptr);
            if(status) { return status; }
            status = mc2err_write_double(&m2e->haar_sum1[i][j], fptr);
            if(status) { return status; }
        }

        for(uint32_t j=0 ; j<width2t ; j++)
        {
            status = mc2err_write_double(&m2e->data_sum2[i][j], fptr);
            if(status) { return status; }
            status = mc2err_write_double(&m2e->lag_sum2[i][j], fptr);
            if(status) { return status; }
            status = mc2err_write_double(&m2e->haar_sum2[i][j], fptr);
            if(status) { return status; }
        }
    }

    // write analysis results
    status = mc2err_write_uint8(&m2e->eqp_cut, fptr);
    if(status) { return status; }
    status = mc2err_write_uint8(&m2e->acf_cut, fptr);
    if(status) { return status; }
    for(uint16_t i=0 ; i<num_block2 ; i++)
    {
        status = mc2err_write_double(&m2e->likelihood[i], fptr);
        if(status) { return status; }
        status = mc2err_write_double(&m2e->penalty[i], fptr);
        if(status) { return status; }

        for(uint16_t j=0 ; j<width ; j++)
        {
            status = mc2err_write_double(&m2e->mean[i][j], fptr);
            if(status) { return status; }
        }

        for(uint32_t j=0 ; j<width2t ; j++)
        {
            status = mc2err_write_double(&m2e->covariance[i][j], fptr);
            if(status) { return status; }
        }
    }

    // close file
    int status2 = fclose(fptr);
    if(status2) { return 4; }
    return 0;
}

uint8_t mc2err_write_uint8(uint8_t *out, FILE *fptr)
{
    size_t num = fwrite(out, sizeof(uint8_t), 1, fptr);
    if(num != 1) { return 4; }
    return 0;
}

uint8_t mc2err_write_uint16(uint16_t *out, FILE *fptr)
{
    for(uint8_t i=0 ; i<2 ; i++)
    {
        uint8_t out2 = *out & 255;
        size_t num = fwrite(&out2, sizeof(uint8_t), 1, fptr);
        if(num != 1) { return 4; }
        *out >>= 8;
    }
    return 0;
}

uint8_t mc2err_write_uint64(uint64_t *out, FILE *fptr)
{
    for(uint8_t i=0 ; i<8 ; i++)
    {
        uint8_t out2 = *out & 255;
        size_t num = fwrite(&out2, sizeof(uint8_t), 1, fptr);
        if(num != 1) { return 4; }
        *out >>= 8;
    }
    return 0;
}

uint8_t mc2err_write_double(double *out, FILE *fptr)
{
    double significand;
    int exponent, shift = 53;
    significand = frexp(*out, &exponent);
    exponent += 1024;
    if(exponent < 0)
    {
        shift += exponent;
        exponent = 0;
    }
    uint64_t out2 = 0;
    if(significand < 0.0)
    {
        significand = -significand;
        out2 = (uint64_t)1<<63;
    }
    out2 += (uint64_t)ldexp(significand, shift);
    if(exponent)
    { out2 -= (uint64_t)1<<52; }
    out2 += (uint64_t)exponent<<52;

    return mc2err_write_uint64(&out2, fptr);
}