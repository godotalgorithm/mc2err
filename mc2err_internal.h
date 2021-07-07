// internal functions & dependencies of the mc2err library
#ifndef MC2ERR_INTERNAL_H
#define MC2ERR_INTERNAL_H

// standard C library dependencies
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// convenient macros for the minimum and maximum of two comparable values
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

// convenient macros for memory checking, allocation, & assignment
#define MEM_CHECK(ptr, type, size)  if((ptr) == NULL) {\
(ptr) = malloc(sizeof(type)*(size)); \
if((ptr) == NULL) { return 3; } \
}
#define MEM_CHECK_SET(ptr, type, size, val)  if((ptr) == NULL) {\
(ptr) = malloc(sizeof(type)*(size)); \
if((ptr) == NULL) { return 3; } \
for(size_t set_counter=0 ; set_counter<(size) ; set_counter++) \
{ (ptr)[set_counter] = (val); } \
}

// aggregate data into partial sums
uint8_t mc2err_aggregate(uint64_t new_data_num, double *new_data, struct mc2err *m2e);

// transfer data from one buffer to the next buffer
uint8_t mc2err_transfer(uint8_t index, struct mc2err *m2e);

// process data up to a new head value, shift buffer when it is half-processed
uint8_t mc2err_process(uint8_t index, uint64_t new_head, struct mc2err *m2e);

#endif