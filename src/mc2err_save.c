// include details of the mc2err_data structure
#include "mc2err_internal.h"

// local macro for writing to a file
#define MC2ERR_FWRITE(PTR, TYPE, NUM, FILE) {\
    size_t _mc2err_fwrite_num = fwrite(PTR, sizeof(TYPE), NUM, FILE);\
    if(NUM != _mc2err_fwrite_num) { fclose(FILE); return 4; }\
}

// Save the mc2err data accumulator 'data' to the file on disk named 'file' in a non-portable binary format.
int mc2err_save(struct mc2err_data *data, char *file)
{
    // check for invalid arguments
    if(data == NULL || file == NULL || *file == '\0')
    { return 1; }

    // local copies of width, length, & max_level for convenience
    const int width = data->width;
    const int length = data->length;
    const int max_level = data->max_level;

    // open the file
    FILE *fptr = fopen(file, "wb");
    if(fptr == NULL) { return 4; }

    // write size info
    MC2ERR_FWRITE(&width, int, 1, fptr);
    MC2ERR_FWRITE(&length, int, 1, fptr);
    MC2ERR_FWRITE(&data->num_chain, int, 1, fptr);
    MC2ERR_FWRITE(&max_level, int, 1, fptr);
    MC2ERR_FWRITE(&data->max_step, long, 1, fptr);
    MC2ERR_FWRITE(data->max_count, long, width, fptr);
    MC2ERR_FWRITE(data->max_pair, long long, width, fptr);

    // write local data
    MC2ERR_FWRITE(data->num_level, int, data->num_chain, fptr);
    MC2ERR_FWRITE(data->num_step, long, data->num_chain, fptr);
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_FWRITE(data->local_count[i], long, 2*data->num_level[i]*length*width, fptr); }
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_FWRITE(data->local_sum[i], double, 2*data->num_level[i]*length*width, fptr); }

    // write global data
    MC2ERR_FWRITE(data->global_count, long, 2*max_level*length*width, fptr);
    MC2ERR_FWRITE(data->global_sum, double, 2*max_level*length*width, fptr);
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_FWRITE(data->pair_count[i], long long, 2*max_level*length*width*width, fptr); }
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_FWRITE(data->pair_sum[i], double, 2*max_level*length*width*width, fptr); }

    // close the file
    int status = fclose(fptr);
    if(status) { return 4; }

    // return without errors
    return 0;
}
