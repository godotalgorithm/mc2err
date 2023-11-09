// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Save the mc2err data accumulator 'data' to the file on disk named 'file'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_save(struct mc2err_data *data, char *file)
{
    // check for invalid arguments
    if(data == NULL || file == NULL || *file == '\0')
    { return 1; }

    // open the file
    FILE *fptr = fopen(file, "w");
    if(fptr == NULL) { return 3; }

    // write data in order
    size_t num_in = 1;
    size_t num_out = fwrite(&data->width, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&data->length, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&data->num_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&data->num_chain, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(&data->num_data, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_chain;
    num_out = fwrite(data->chain_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fwrite(data->chain_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    for(int i=0 ; i<data->num_chain ; i++)
    {
        num_in = 2*data->chain_level[i]*data->length*data->width;
        num_out = fwrite(data->chain_sum[i], sizeof(double), num_in, fptr);
        if(num_in != num_out) { fclose(fptr); return 3; }
    }

    num_in = data->num_level*(data->num_level+3)*data->length/2;
    num_out = fwrite(data->data_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->width/2;
    num_out = fwrite(data->data_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->length;
    num_out = fwrite(data->pair_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->length*data->width*(data->width+1)/2;
    num_out = fwrite(data->pair_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    // close the file
    int status = fclose(fptr);
    if(status) { return 3; }

    // return without errors
    return 0;
}
