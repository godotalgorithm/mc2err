// mc2err header file
#include "mc2err.h"
#include "mc2err_internal.h"

// append accumulated statistical data from the source to the target mc2err data structure
uint8_t mc2err_append(struct mc2err_data *source, struct mc2err_data *target)
{
    // return if there is nothing to do
    if(source->num_input == 0)
    { return 0; }

    // check for a global data overflow
    if(source->num_input > UINT64_MAX-target->num_input)
    { return 2; }

    // check for a data width mismatch
    if(source->width != target->width)
    { return 3; }

    // resize workspaces in the target as necessary
    uint8_t status = mc2err_expand(0, source->num_block, target);
    if(status) { return status; }

    // local loop size variables
    uint16_t num_block2t = source->num_block*(source->num_block+1)/2;
    uint16_t width = source->width;
    uint32_t width2t = source->width*(source->width+1)/2;

    // aggregate block data
    for(uint16_t i=0 ; i<num_block2t ; i++)
    {
        target->data_num[i] += source->data_num[i];
        target->lag_num[i] += source->lag_num[i];
        target->haar_num[i] += source->haar_num[i];

        for(uint16_t j=0 ; j<width ; j++)
        {
            target->data_sum1[i][j] += source->data_sum1[i][j];
            target->haar_sum1[i][j] += source->haar_sum1[i][j];
        }

        for(uint32_t j=0 ; j<width2t ; j++)
        {
            target->data_sum2[i][j] += source->data_sum2[i][j];
            target->lag_sum2[i][j] += source->lag_sum2[i][j];
            target->haar_sum2[i][j] += source->haar_sum2[i][j];
        }
    }

    // update the total amount of data & return
    target->num_input += source->num_input;
    return 0;
}
