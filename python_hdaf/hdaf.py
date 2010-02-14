
from ctypes import *

#from numpy import *
import numpy


c_hdaf = cdll.LoadLibrary('/usr/lib/libhdaf.so.1')

hdaf_data_dir = '/workspace/hdaf_data'


c_hdaf.sigma_from_cutoff_frequency.argtyprs = [ c_double, c_int ]
c_hdaf.sigma_from_cutoff_frequency.restype = c_double

def freq_to_sigma( f, m ):
    return c_hdaf.sigma_from_cutoff_frequency( c_double(f), c_int(m) )



c_hdaf.hdaf_equate_arrays.argtypes = [ c_void_p, c_void_p, c_ulong ]

c_hdaf.hdaf_free_array.argtypes = [ c_void_p ]




c_hdaf.get_hdaf_kernel.argtypes=[ c_void_p, c_void_p, c_double, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel.restype = c_int


def hdaf_kernel( sampling_period, order, sigma ):
    
    ptr = c_void_p()
    size = c_ulong()

    error = c_hdaf.get_hdaf_kernel( byref(ptr), byref(size), c_double(sampling_period), c_int(order), c_double(sigma), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel'; code ", error
    
    
    ar = numpy.zeros(size.value)
    c_hdaf.hdaf_equate_arrays( ar.ctypes.data_as(c_void_p), ptr, size )
    c_hdaf.hdaf_free_array( ptr )
    ptr = c_void_p()
    
    return ar




c_hdaf.get_hdaf_kernel_arbitrary_points.argtypes=[ c_void_p, c_void_p, c_int, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel_arbitrary_points.restype = c_int


def hdaf_kernel_bypts( eval_points, order, sigma ):
    
    eval_points = numpy.array(eval_points)
    ar = numpy.zeros(len(eval_points))
    length = c_int(len(eval_points))
    
    error = c_hdaf.get_hdaf_kernel_arbitrary_points( eval_points.ctypes.data_as(c_void_p), ar.ctypes.data_as(c_void_p), length, c_int(order), c_double(sigma), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel_arbitrary_points'; code ", error
    
    return ar






c_hdaf.get_hdaf_kernel_lp.argtypes=[ c_void_p, c_void_p, c_double, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel_lp.restype = c_int


def lp_hdaf_kernel( sampling_period, order, cutoff_frequency ):
    
    ptr = c_void_p()
    size = c_ulong()

    error = c_hdaf.get_hdaf_kernel_lp( byref(ptr), byref(size), c_double(sampling_period), c_int(order), c_double(cutoff_frequency), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel_lp'; code ", error
    
    
    ar = numpy.zeros(size.value)
    c_hdaf.hdaf_equate_arrays( ar.ctypes.data_as(c_void_p), ptr, size )
    c_hdaf.hdaf_free_array( ptr )
    ptr = c_void_p()
    
    return ar










"""

int get_hdaf_kernel(double *&kernel, int & kernel_size, double sampling_period, int order, double sigma, const char *hdaf_data_dir );

int get_hdaf_kernel_lp(double *&kernel, int & kernel_size, double sampling_period, int order, double cutoff_frequency, const char *hdaf_data_dir );

int get_hdaf_kernel_bp(double *&kernel, int & kernel_size, double sampling_period, int low_pass_order, double low_pass_frequency, int high_pass_order, double high_pass_frequency, const char *hdaf_data_dir );
"""



        





