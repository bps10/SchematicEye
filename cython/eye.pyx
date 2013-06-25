# distutils: language = c++
# distutils: sources = ../src/eye_eye.cpp ../src/eye_analysis.cpp

cdef extern from "../src/eye_analysis.hh" namespace "SchematicEye":
    cdef cppclass Analysis:
        Analysis()
        float *IntensityAnalysis(float, float, float, float)

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
 
# Import the Python-level symbols of numpy
import numpy as np
 
# Import the C-level symbols of numpy
cimport numpy as np
 
# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()
 
# We need to build an array-wrapper class to deallocate our array when
# the Python object is deleted.

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size
 
    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array
 
        This cannot be done in the constructor as it must recieve C-level
        arguments.
 
        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data            
 
        """
        self.data_ptr = data_ptr
        self.size = size
 
    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_FLOAT, self.data_ptr)
        return ndarray
 
    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
 
 
def py_eye(float object_distance, float off_axis, float pupil_size, float diopters):
    """ Python binding of the 'compute' function in 'c_code.c' that does
        not copy the data allocated in C.
    """
    cdef float *array
    cdef np.ndarray ndarray
    # Call the C function
    array = Analysis().IntensityAnalysis(object_distance, off_axis, pupil_size, diopters)
 
    size = 399
    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(size, <void*> array) 
    ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)
 
 
    return ndarray