// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Load the mc2err data structure 'm2e' from disk in a file named 'filename'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_load(struct mc2err_data *m2e, char *filename)
{
    // check for invalid arguments
    if(m2e == NULL || filename == NULL || *filename == '\0')
    { return 1; }

    // open the file
    FILE *fptr = fopen(filename, "w");
    if(fptr == NULL) { return 4; }

    // read data in order, allocate memory as necessary
    size_t num_in = 1;
    size_t num_out = fread(&m2e->width, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fread(&m2e->length, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fread(&m2e->num_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fread(&m2e->num_chain, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fread(&m2e->num_data, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = m2e->num_chain;
    m2e->chain_level = (int*)malloc(sizeof(int)*num_in);
    if(m2e->chain_level == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->chain_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    m2e->chain_count = (long int*)malloc(sizeof(long int)*num_in);
    if(m2e->chain_count == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->chain_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    m2e->chain_sum = (double**)malloc(sizeof(double*)*num_in);
    if(m2e->chain_sum == NULL) { fclose(fptr); return 5; }
    for(int i=0 ; i<m2e->num_chain ; i++)
    {
        num_in = 2*m2e->chain_level[i]*m2e->length*m2e->width;
        m2e->chain_sum[i] = (double*)malloc(sizeof(double)*num_in);
        if(m2e->chain_sum[i] == NULL) { fclose(fptr); return 5; }
        num_out = fread(m2e->chain_sum[i], sizeof(double), num_in, fptr);
        if(num_in != num_out) { fclose(fptr); return 4; }
    }

    num_in = (m2e->num_level+1)*m2e->length;
    m2e->data_count = (long int*)malloc(sizeof(long int)*num_in);
    if(m2e->data_count == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->data_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*m2e->length*m2e->width;
    m2e->data_sum = (double*)malloc(sizeof(double)*num_in);
    if(m2e->data_sum == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->data_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length;
    m2e->pair_count = (long int*)malloc(sizeof(long int)*num_in);
    if(m2e->pair_count == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->pair_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2;
    m2e->pair_sum = (double*)malloc(sizeof(double)*num_in);
    if(m2e->pair_sum == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->pair_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    m2e->pair_tail = (double*)malloc(sizeof(double)*num_in);
    if(m2e->pair_tail == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->pair_tail, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = 1;
    num_out = fread(&m2e->eqp_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fread(&m2e->acf_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = 2*m2e->num_level*m2e->length;
    m2e->eqp_p_value = (double*)malloc(sizeof(double)*num_in);
    if(m2e->eqp_p_value == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->eqp_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    m2e->acf_p_value = (double*)malloc(sizeof(double)*num_in);
    if(m2e->acf_p_value == NULL) { fclose(fptr); return 5; }
    num_out = fread(m2e->acf_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    // close the file
    int status = fclose(fptr);
    if(status) { return 4; }

    // return without errors
    return 0;
}
