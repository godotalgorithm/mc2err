// include details of the mc2err_data structure
#include "mc2err_internal.h"

// Load the mc2err data accumulator 'data' from the file on disk named 'file'.
// (This file is not portable between computing environments with different endianness or integer sizes.)
int mc2err_load(struct mc2err_data *data, char *file)
{
    // check for invalid arguments
    if(data == NULL || file == NULL || *file == '\0')
    { return 1; }

    // open the file
    FILE *fptr = fopen(file, "w");
    if(fptr == NULL) { return 3; }

    // read data in order, allocate memory as necessary
    size_t num_in = 1;
    size_t num_out = fread(&data->width, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fread(&data->length, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fread(&data->num_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fread(&data->num_chain, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_out = fread(&data->num_data, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_chain;
    data->chain_level = (int*)malloc(sizeof(int)*num_in);
    if(data->chain_level == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->chain_level, sizeof(int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    data->chain_count = (long int*)malloc(sizeof(long int)*num_in);
    if(data->chain_count == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->chain_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    data->chain_sum = (double**)malloc(sizeof(double*)*num_in);
    if(data->chain_sum == NULL) { fclose(fptr); return 4; }
    for(int i=0 ; i<data->num_chain ; i++)
    {
        num_in = 2*data->chain_level[i]*data->length*data->width;
        data->chain_sum[i] = (double*)malloc(sizeof(double)*num_in);
        if(data->chain_sum[i] == NULL) { fclose(fptr); return 4; }
        num_out = fread(data->chain_sum[i], sizeof(double), num_in, fptr);
        if(num_in != num_out) { fclose(fptr); return 3; }
    }

    num_in = data->num_level*(data->num_level+3)*data->length/2;
    data->data_count = (long int*)malloc(sizeof(long int)*num_in);
    if(data->data_count == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->data_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->width/2;
    data->data_sum = (double*)malloc(sizeof(double)*num_in);
    if(data->data_sum == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->data_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->length;
    data->pair_count = (long int*)malloc(sizeof(long int)*num_in);
    if(data->pair_count == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->pair_count, sizeof(long int), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    num_in = data->num_level*(data->num_level+3)*data->length*data->length*data->width*(data->width+1)/2;
    data->pair_sum = (double*)malloc(sizeof(double)*num_in);
    if(data->pair_sum == NULL) { fclose(fptr); return 4; }
    num_out = fread(data->pair_sum, sizeof(double), num_in, fptr);
    if(num_in != num_out) { fclose(fptr); return 3; }

    // close the file
    int status = fclose(fptr);
    if(status) { return 3; }

    // return without errors
    return 0;
}
