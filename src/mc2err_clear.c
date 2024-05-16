// include definition of NULL
#include <stdlib.h>

// include details of the mc2err_analysis structure
#include "mc2err.h"

// Clear and deallocate the memory of the analysis results 'analysis' after it is no longer needed
// or before it is reused in another call to 'mc2err_output'.
int mc2err_clear(struct mc2err_analysis *analysis)
{
    // check for invalid arguments
    if(analysis == NULL)
    { return 1; }

    // free all pointers
    MC2ERR_FREE(analysis->count);
    MC2ERR_FREE(analysis->mean);
    MC2ERR_FREE(analysis->variance);
    MC2ERR_FREE(analysis->variance0);
    MC2ERR_FREE(analysis->eqp_p);
    MC2ERR_FREE(analysis->acc_p);

    // set sizes to zero for hygiene
    analysis->width = 0;
    analysis->length = 0;
    analysis->num_level = 0;

    // return without errors
    return 0;
}
