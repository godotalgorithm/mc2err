# C FFI builder for mc2err library

from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("struct mc2err { uint64_t data_num; ...; };")
ffibuilder.cdef("void mc2err_initialize(uint8_t data_dim, uint64_t data_max, double acf_precision, double eqp_precision, struct mc2err *m2e);")
ffibuilder.cdef("void mc2err_finalize(struct mc2err *m2e);")
ffibuilder.cdef("uint64_t mc2err_input(uint64_t new_data_num, double *new_data, struct mc2err *m2e);")
ffibuilder.cdef("uint8_t mc2err_output(double *mean, double *error, double *size, struct mc2err *m2e);")

ffibuilder.set_source("_mc2err_cffi",'#include "mc2err.h"', sources=['mc2err.c'], libraries=['m'])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
