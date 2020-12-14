# mc2err Python wrapper
from _mc2err_cffi import ffi, lib

class ExpValue:
    """object wrapper for mc2err functionality"""

    def __init__(self, acf_precision = 0.01, eqp_precision = 0.01, data_max = 18446744073709551615):
        self._workspace = ffi.new("struct mc2err *")
        lib.mc2err_initialize(1, data_max, acf_precision, eqp_precision, self._workspace)
        self._mean = ffi.new("double *")
        self._error = ffi.new("double *")
        self._size = ffi.new("double *")
        self.unprocessed_data = False

    def __del__(self):
        lib.mc2err_finalize(self._workspace)

    def add(self, data):
        lib.mc2err_input(len(data), data, self._workspace)
        self.unprocessed_data = True

    def _update(self):
        if self.unprocessed_data:
            lib.mc2err_output(self._mean, self._error, self._size, self._workspace)
            self.unprocessed_data = False

    def mean(self):
        self._update()
        return self._mean[0]

    def error(self):
        self._update()
        return self._error[0]

    def size(self):
        """effective sample size"""
        self._update()
        return self._size[0]

# basic testing of wrapped structure & functions
if __name__ == "__main__":

    workspace = ffi.new("struct mc2err *")
    mean = ffi.new("double *")
    error = ffi.new("double *")
    size = ffi.new("double *")

    lib.mc2err_initialize(1, 18446744073709551615, 1e-3, 1e-3, workspace)

    data = [0,1]*500

    lib.mc2err_input(1000, data, workspace)

    lib.mc2err_output(mean, error, size, workspace)
    print(mean[0],"+/-",error[0])

    data = [0,1]*1000

    lib.mc2err_input(2000, data, workspace)

    lib.mc2err_output(mean, error, size, workspace)
    print(mean[0],"+/-",error[0])

    lib.mc2err_finalize(workspace)

    print(workspace.data_num)
    print(dir(lib))
