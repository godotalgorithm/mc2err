// include details of the mc2err_data structure
#include "mc2err_internal.h"

// local macro for reading from a file & allocating memory as necessary
#define MC2ERR_FREAD(PTR, TYPE, NUM, FILE) {\
    size_t _mc2err_fread_num = fread(PTR, sizeof(TYPE), NUM, FILE);\
    if(NUM > 1)\
    {\
        PTR = (TYPE*)malloc(sizeof(TYPE)*NUM);\
        if(PTR == NULL) { fclose(FILE); return 4; }\
    }\
    if(NUM != _mc2err_fread_num) { fclose(FILE); return 3; }\
}

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

    // read data in order
    MC2ERR_FREAD(&data->width, int, 1, fptr);
    MC2ERR_FREAD(&data->length, int, 1, fptr);
    MC2ERR_FREAD(&data->num_level, int, 1, fptr);
    MC2ERR_FREAD(&data->num_chain, int, 1, fptr);
    MC2ERR_FREAD(&data->num_data, long int, 1, fptr);
    MC2ERR_FREAD(data->chain_level, int, data->num_chain, fptr);
    MC2ERR_FREAD(data->chain_count, long int, data->num_chain, fptr);
    MC2ERR_MALLOC(data->chain_sum, double*, data->num_chain);
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_FREAD(data->chain_sum[i], double, 2*data->chain_level[i]*data->length*data->width, fptr); }
    size_t global_size = data->num_level*(data->num_level+3)*data->length/2;
    MC2ERR_FREAD(data->data_count, long int, global_size, fptr);
    MC2ERR_FREAD(data->data_sum, double, global_size*data->width, fptr);
    MC2ERR_FREAD(data->pair_count, long int, global_size*2*data->length, fptr);
    MC2ERR_FREAD(data->pair_sum, double, global_size*data->length*data->width*(data->width+1), fptr);

    // close the file
    int status = fclose(fptr);
    if(status) { return 3; }

    // return without errors
    return 0;
}
