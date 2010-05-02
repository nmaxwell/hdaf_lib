

from ctypes import *
from ctypes.util import *

import numpy


c_hdaf = cdll.LoadLibrary(find_library('hdaf'))

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








c_hdaf.render_png_scalar.argtypes = [ c_char_p, c_int, c_int, c_void_p, c_int,  c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double  ]

c_hdaf.render_png_scalar_resample.argtypes = [ c_char_p, c_void_p, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double,  c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double,  c_double, c_double, c_double ]



def write_png( data, fname, major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ), ordering='rm' ):
    
    nx,ny = numpy.shape(data)
    ordering_int = 0
    if ordering not in [ 'rm', 'row_major', 'row' ]:
        ordering_int = 1
    
    try:
        err = c_hdaf.render_png_scalar( c_char_p(fname), c_int(nx), c_int(ny), data.ctypes.data_as(c_void_p), c_int(ordering_int), center, major_scale, c_double(red_params[0]),c_double(red_params[1]),c_double(red_params[2]),c_double(red_params[3]), c_double(green_params[0]),c_double(green_params[1]),c_double(green_params[2]),c_double(green_params[3]), c_double(blue_params[0]),c_double(blue_params[1]),c_double(blue_params[2]),c_double(blue_params[3]) )
    except:
        print "write_png error, c_waveprop.write_png, exception."


def write_png_resample( data, fname, grid_new, grid_old,  major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ), default_color=(0.0,0.0,0.0) ):
    
    nx,ny = numpy.shape(data)
    
    try:
        err = c_hdaf.render_png_scalar_resample( c_char_p(fname), data.ctypes.data_as(c_void_p), grid_old.a1, grid_old.b2, grid_old.a2, grid_old.b2, grid_old.n1, grid_old.n2, grid_new.a1, grid_new.b2, grid_new.a2, grid_new.b2, grid_new.n1, grid_new.n2, center, major_scale, c_double(red_params[0]),c_double(red_params[1]),c_double(red_params[2]),c_double(red_params[3]), c_double(green_params[0]),c_double(green_params[1]),c_double(green_params[2]),c_double(green_params[3]), c_double(blue_params[0]),c_double(blue_params[1]),c_double(blue_params[2]),c_double(blue_params[3]), c_double(default_color[0]), c_double(default_color[1]), c_double(default_color[2]) )
    except:
        print "write_png_resample error, c_waveprop.render_png_scalar, exception."









class grid2d:
    
    def __init__(self, a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, dx2=None, ):
        
        if b1 != None and b2 != None:
            if n1 != None and n2 != None:
                dx1 = float(b1-a1)/n1
                dx2 = float(b2-a2)/n2
            elif dx1 != None and dx2 != None:
                n1 = float(b1-a1)/dx1
                n2 = float(b2-a2)/dx2
        elif dx1 != None and dx2 != None:
            if n1 != None and n2 != None:
                b1 = dx1*n1 + a1
                b2 = dx2*n2 + a2
            elif b1 != None and b2 != None:
                n1 = float(b1-a1)/dx1
                n2 = float(b2-a2)/dx2
        elif n1 != None and n2 != None:
            if dx1 != None and dx2 != None:
                b1 = float(dx1)*n1 + a1
                b2 = float(dx2)*n2 + a2
            elif b1 != None and b2 != None:
                dx1 = float(b1-a1)/n1
                dx2 = float(b2-a2)/n2
        
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.b1 = float(n1*dx1+a1)
        self.b2 = float(n2*dx2+a2)
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.dx1 = float(dx1)
        self.dx2 = float(dx2)
        
        [self.X1, self.X2] = self.coord_mesh()
    
    def coord_mesh(self ):
        return numpy.mgrid[ self.a1: self.a1+self.n1*self.dx1: self.dx1, self.a2: self.a2+self.n2*self.dx2: self.dx2 ]
    
    def zeros(self ):
        return numpy.zeros((self.n1,self.n2))
    
    def ones(self ):
        return numpy.ones((self.n1,self.n2))
    
    def evaluate(self, f):
        g = numpy.vectorize(f)
        try:
            return g(self.X1, self.X2)
        except:
            print "grid2d error: evaluate"
            return self.zeros()
    
    def shape(self ):
        return (self.n1, self.n2 )
    
    def index1(self, x ):
        return int((x-self.a1)/self.dx1)
    
    def index2(self, y ):
        return int((y-self.a2)/self.dx2)
    
    def index(self, x):
        return self.index1(x[0]), self.index2(x[1])
    
    def discretize(self, x):
        i = self.index1(x[0])
        j = self.index2(x[1])
        return self.X1[i][j], self.X2[i][j]



c_hdaf.laplacian2d_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_int, c_int, c_double, c_double ]
c_hdaf.laplacian2d_free.argtyprs = [ c_void_p ]
c_hdaf.laplacian2d_execute.argtyprs = [ c_void_p, c_void_p, c_void_p ]


class hdafLaplacian2D:
    
    def __init__(self, ):
        self.data=c_void_p(0)
    
    def __del__(self ):
        self.free()
    
    def free(self ):
        try:
            if self.data.value is not None:
                err = c_hdaf.laplacian2d_free( byref(self.data) )
                if err != 0:
                    print "error: hdafLaplacian2D.free; error", err
        except:
            print "error: hdafLaplacian2D.free"
    
    
    def set(self, grid=None, a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, dx2=None, m1=8, m2=8, gamma1=0.8, gamma2=0.8):
        
        if grid is not None:
            self.grid = grid
        else:
            self.grid = grid2d( a1=a1, b1=b1, a2=a2, b2=b2, n1=n1, n2=n2, dx1=dx1, dx2=dx2 )
        
        self.free()
        
        try:
            self.data = c_void_p(0)
            err = c_hdaf.laplacian2d_init( byref(self.data), c_int(self.grid.n1), c_int(self.grid.n2), c_double(self.grid.b1-self.grid.a1), c_double(self.grid.b2-self.grid.a2), c_int(m1), c_int(m2), c_double(gamma1), c_double(gamma2) )
            
            if err != 0 or self.data.value==0:
                print "error: hdafLaplacian2D initialization error: ", err
        except:
            print "error: hdafLaplacian2D initialization error, exception."
        
    
    def __call__(self, input, output):
        try:
            err = c_hdaf.laplacian2d_execute( self.data, input.ctypes.data_as(c_void_p), output.ctypes.data_as(c_void_p) )
            if err != 0:
                print "error: hdafLaplacian2D call error, error: ", err
        except:
            print "error: hdafLaplacian2D call error, exception."
            quit()










    



