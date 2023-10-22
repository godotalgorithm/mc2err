// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Save the mc2err data structure 'm2e' to disk in a file named 'filename'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_save(struct mc2err_data *m2e, char *filename)
{
    // check for invalid arguments
    if(m2e == NULL || filename == NULL || *filename == '\0')
    { return 1; }

    // open the file
    FILE *fptr = fopen(filename, "w");
    if(fptr == NULL) { return 4; }

    // write data in order
    size_t num_in = 1;
    size_t num_out = fwrite(&m2e->width, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(&m2e->length, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(&m2e->num_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(&m2e->num_chain, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(&m2e->num_data, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = m2e->num_chain;
    num_out = fwrite(m2e->chain_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(m2e->chain_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    for(int i=0 ; i<m2e->num_chain ; i++)
    {
        num_in = 2*m2e->chain_level[i]*m2e->length*m2e->width;
        num_out = fwrite(m2e->chain_sum[i], sizeof(double), num_in, fptr);
        if(num_in != num_out) { fclose(fptr); return 4; }
    }

    num_in = (m2e->num_level+1)*m2e->length;
    num_out = fwrite(m2e->data_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*m2e->length*m2e->width;
    num_out = fwrite(m2e->data_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length;
    num_out = fwrite(m2e->pair_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2;
    num_out = fwrite(m2e->pair_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    num_out = fwrite(m2e->pair_tail, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = 1;
    num_out = fwrite(&m2e->eqp_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(&m2e->acf_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_in = 2*m2e->num_level*m2e->length;
    num_out = fwrite(m2e->eqp_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    num_out = fwrite(m2e->acf_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 4; }

    // close the file
    int status = fclose(fptr);
    if(status) { return 4; }

    // return without errors
    return 0;
}
