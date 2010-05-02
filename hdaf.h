#ifndef HDAF_H
#define HDAF_H

//#define HDAF_STANDALONE


/*

documentation:

    nicholas.maxwell@gmail.com for questions.

    Google: 'Hermite distributed approximating functionals as almost ideal low pass filters'


    
get_hdaf_kernel:

    Gets the convolution kernel for the (m,sigma) hdaf, see
        Google: 'Hermite distributed approximating functionals as almost ideal low pass filters'
    for details.
    
    *kernel will be free'ed and re malloc'ed. You may pass it as a null pointer.
    
    kernel_size will be the array size of kernel.
    
    kernel[0] is associated to time=0; kernel[k] is associated to time = h*sampling_period.
    sampling_period is the sampling period of kernel.

    All hdaf kernels are even functions, so are symmetric in time around t=0; kernel is only half of the actual kernel.
    
    so for example, to apply this to a set of data, double *y, 
    
    y[m] = 0.0;
    for (int k=0; k<kernel_size; k++)
        y[m] += kernel[k]*y[m+k];
    for (int k=1; k<kernel_size; k++)
        y[m] += kernel[k]*y[m-k];
    
    assuming that y is sampled with a period of sampling_period.
    
    get_hdaf_kernel reads data files containing the raw kernel data.
    
    get_hdaf_kernel will return 0 upon successfull completion, !0 if some error occurs.
        

*/

/*
    read_std_hdaf_kernel_file error codes:
        0 no error
        1 error opening file
        2 file too short
        3 buffer malloc error
        4 file read error
        5 corrupt header data
        6 not enough data in file
        7 kernel malloc error
        8 corrupt kernel data
*/



#include <stdio.h>
#include <math.h>

#define hdaf_max_file_name_length 200
#define hdaf_sqrt2 1.414213562373095
#define hdaf_pi 3.14159265358979323846264338
#define hdaf_2pi    6.2831853071795865
#define hdaf_sqrtpi  1.772453850905516



 #ifdef __cplusplus
 extern "C" {
 #endif

void hdaf_equate_arrays( double *x, double *y, unsigned int size );

void hdaf_free_array( double *p );


void std_kernel_file_naming_convention( char * file_name, const char *hdaf_data_dir, int hdaf_order );

double sigma_from_cutoff_frequency( double cutoff_frequency, int hdaf_order );


int write_std_hdaf_kernel_file( const char * file_name, int hdaf_order, double step_size, int n_points, double *std_kernel );

int read_std_hdaf_kernel_file( const char * file_name, int *hdaf_order, double *step_size, int *n_points, double **std_kernel );

double hdaf_interpolate( double * std_kernel, int n_points, double step_size, double x );


int get_hdaf_kernel(double **kernel, int * kernel_size, double sampling_period, int order, double sigma, const char *hdaf_data_dir );


int get_hdaf_kernel_arbitrary_points(double *eval_points, double *kernel, int n_points, int order, double sigma, const char *hdaf_data_dir );


int get_hdaf_kernel_lp(double **kernel, int * kernel_size, double sampling_period, int order, double cutoff_frequency, const char *hdaf_data_dir );

int get_hdaf_kernel_bp(double **kernel, int * kernel_size, double sampling_period, int low_pass_order, double low_pass_frequency, int high_pass_order, double high_pass_frequency, const char *hdaf_data_dir );






#ifndef HDAF_STANDALONE


// tools:


int render_png_scalar(
    char *fname, int nx, int ny, double *data, int ordering,
    double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue );


int render_png_scalar_resample(
    char *fname, double *data,
    double a1, double b1, double a2, double b2, int nx, int ny,
    double a1_new, double b1_new, double a2_new, double b2_new, int nx_new, int ny_new,
    double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue,
    double red_default, double green_default, double blue_default );




//differential operators:

int laplacian2d_init( void **data, int n1, int n2, double L1, double L2, int m1, int m2, double gamma1, double gamma2 );
int laplacian2d_free( void **data );
int laplacian2d_execute( void *data, double *in, double *out );




#endif



// internal

double hdaf_data_access( int k, double *kernel, int n_points );

double hdaf_neville_Pij(double x, int i, int j, int n_data, double * data, double step_size );

double hdaf_interpolate( double * kernel, int n_points, double step_size, double x );



//double hdaf_fourier_transform(double omega, int m, double sigma);

 #ifdef __cplusplus
 }
 #endif









#endif
















/*


compute_hdaf_kernel_truncation:

    The hdafs have an infinite impulse response (are IIR filters).
    
    This function returns where the impulse response, which is the convolution kernel, can be truncated.
    
    The error that this introduces is computed exactly, and ensured to be within some tolerance range (eps_min, eps_max).
    

compute_hdaf_kernel_width_double:
    
    The hdafs have an infinite impulse response (are IIR filters).
    
    This function computes the size kernel to use when using a sampling period of sampling_period.
    
    The error introduced can be arbitraryily small, but is specifie in the cpp file to be (1E-17, 5E-17), which is below machine precision for double percision floats.
    
    This calls compute_hdaf_kernel_truncation



double compute_hdaf_kernel_truncation(int m, double sigma, double eps_min, double eps_max);

int compute_hdaf_kernel_width_double(int m, double sigma, double sampling_period);

*/
