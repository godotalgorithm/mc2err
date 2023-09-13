// include details of the mc2err_data structure
#include "mc2err_internal.h"

// save all data to disk using fwrite commands
int mc2err_save(struct mc2err_data *m2e, char *filename)
{
    // open the file
    FILE *fptr = fopen(filename, "w");
    if(fptr == NULL) { return 3; }

    // write data in order
    size_t num_in = 1;
    size_t num_out = fwrite(&m2e->width, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&m2e->length, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&m2e->num_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&m2e->num_chain, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&m2e->num_data, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = m2e->num_chain;
    num_out = fwrite(m2e->chain_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(m2e->chain_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    for(int i=0 ; i<m2e->num_chain ; i++)
    {
        num_in = 2*m2e->chain_level[i]*m2e->length*m2e->width;
        num_out = fwrite(m2e->chain_sum[i], sizeof(double), num_in, fptr);
        if(num_in != num_out) { fclose(fptr); return 3; }
    }

    num_in = (m2e->num_level+1)*m2e->length;
    num_out = fwrite(m2e->data_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = (m2e->num_level+1)*m2e->length*m2e->width;
    num_out = fwrite(m2e->data_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length;
    num_out = fwrite(m2e->pair_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = (m2e->num_level+1)*(m2e->num_level+1)*m2e->length*m2e->length*m2e->width*(m2e->width+1)/2;
    num_out = fwrite(m2e->pair_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = m2e->num_level*(m2e->num_level+1)*m2e->length*m2e->width*(m2e->width+1)/2;
    num_out = fwrite(m2e->pair_tail, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = 1;
    num_out = fwrite(&m2e->eqp_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&m2e->acf_cut, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = 2*m2e->num_level*m2e->length;
    num_out = fwrite(m2e->eqp_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(m2e->acf_p_value, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = 1;
    num_out = fwrite(&m2e->lapack_info, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    // close the file
    int status = fclose(fptr);
    if(status) { return 3; }

    // return without errors
    return 0;
}
