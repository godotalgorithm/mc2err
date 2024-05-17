// include details of the mc2err_data structure
#include "mc2err_internal.h"

// local macro for reading from a file
#define MC2ERR_FREAD(PTR, TYPE, NUM, FILE) {\
    size_t _mc2err_fread_num = fread(PTR, sizeof(TYPE), NUM, FILE);\
    if(NUM != _mc2err_fread_num) { fclose(FILE); return 4; }\
}

// Load the mc2err data accumulator 'data' from the file on disk named 'file' in a non-portable binary format.
int mc2err_load(struct mc2err_data *data, char *file)
{
    // check for invalid arguments
    if(data == NULL || file == NULL || *file == '\0')
    { return 1; }

    // open the file
    FILE *fptr = fopen(file, "rb");
    if(fptr == NULL) { return 4; }

    // read main size info
    int width, length;
    MC2ERR_FREAD(&data->width, int, 1, fptr);
    MC2ERR_FREAD(&data->length, int, 1, fptr);
    MC2ERR_FREAD(&data->num_chain, int, 1, fptr);
    MC2ERR_FREAD(&data->max_level, int, 1, fptr);

    // local copies of width & length for convenience
    const int width = data->width;
    const int length = data->length;
    const int max_level = data->max_level;

    // initialize outer pointers
    MC2ERR_MALLOC(data->max_count, long, width);
    MC2ERR_MALLOC(data->max_pair, long long, width);
    MC2ERR_MALLOC(data->num_level, int, data->num_chain);
    MC2ERR_MALLOC(data->num_step, long, data->num_chain);
    MC2ERR_MALLOC(data->local_count, long*, data->num_chain);
    MC2ERR_MALLOC(data->local_sum, double*, data->num_chain);
    MC2ERR_MALLOC(data->global_count, long, 2*max_level*length*width);
    MC2ERR_MALLOC(data->global_sum, double, 2*max_level*length*width);
    MC2ERR_MALLOC(data->pair_count, long long*, 2*max_level*length);
    MC2ERR_MALLOC(data->pair_sum, double*, 2*max_level*length);

    // read remaining size info
    MC2ERR_FREAD(&data->max_step, long, 1, fptr);
    MC2ERR_FREAD(data->max_count, long, width, fptr);
    MC2ERR_FREAD(data->max_pair, long long, width, fptr);

    // read local data for secondary size info
    MC2ERR_FREAD(data->num_level, int, data->num_chain, fptr);
    MC2ERR_FREAD(data->num_step, long, data->num_chain, fptr);

    // initialize inner pointers
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_MALLOC(data->local_count[i], long, 2*data->num_level[i]*length*width); }
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_MALLOC(data->local_sum[i], double, 2*data->num_level[i]*length*width); }
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_MALLOC(data->pair_count[i], long long, 2*max_level*length*width*width); }
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_MALLOC(data->pair_sum[i], double, 2*max_level*length*width*width); }

    // read remaining local data
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_FREAD(data->local_count[i], long, 2*data->num_level[i]*length*width, fptr); }
    for(int i=0 ; i<data->num_chain ; i++)
    { MC2ERR_FREAD(data->local_sum[i], double, 2*data->num_level[i]*length*width, fptr); }

    // read global data
    MC2ERR_FREAD(data->global_count, long, 2*max_level*length*width, fptr);
    MC2ERR_FREAD(data->global_sum, double, 2*max_level*length*width, fptr);
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_FREAD(data->pair_count[i], long long, 2*max_level*length*width*width, fptr); }
    for(size_t i=0 ; i<2*max_level*length ; i++)
    { MC2ERR_FREAD(data->pair_sum[i], double, 2*max_level*length*width*width, fptr); }

    // close the file
    int status = fclose(fptr);
    if(status) { return 4; }

    // return without errors
    return 0;
}
