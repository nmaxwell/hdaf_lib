#ifndef HDAF_CPP
#define HDAF_CPP

#include "hdaf.h"

//#include <iostream>
#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

 #ifdef __cplusplus
 extern "C" {
 #endif

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );





void hdaf_equate_arrays( double *x, double *y, unsigned int size )
{
    if ( size != 0 && x != NULL && y != NULL )
    {
        for (unsigned int k=0; k<size; k++)
            x[k] = y[k];
    }
}

void hdaf_free_array( double *p )
{
    if (p != NULL)
        free( p );
    
    p = NULL;
}



void std_kernel_file_naming_convention( char * file_name, const char *hdaf_data_dir, int hdaf_order )
{
    sprintf( file_name, "%s/%05d.hdaf", hdaf_data_dir, hdaf_order );
}

double sigma_from_cutoff_frequency( double cutoff_frequency, int hdaf_order )
{
    return hdaf_2pi*sqrt((double)(2*hdaf_order+1))/(cutoff_frequency);
}


int write_std_hdaf_kernel_file( const char * file_name, int hdaf_order, double step_size, int n_points, double *kernel )
{
/*
    
    file contents:
        int hdaf_order (4 bytes)
        int n_points (4 bytes)
        doulbe step_size (8 bytes)
        immediatly after that is the kernel data, which is an array of doubles (8 bytes) of length n_points.
    
    error codes:
        0 no error
        1 invalid arguements
        2 corrupt kernel data
        3 error opening file
        4 error writing header
        5 error writing kernel data
    
*/
    
    if ( kernel==NULL || hdaf_order<0 || n_points<=0 || step_size<=0 || isinf(step_size) || isnan(step_size) || !(step_size==step_size) )
        return 1;
    
    for (int k=0; k<n_points; k++)
        if ( isinf(kernel[k]) || isnan(kernel[k]) || !(kernel[k]==kernel[k]) )
            return 2;
    
    ofstream out;
    out.open(file_name, ios::out | ios::binary );
    if ( !out.good() || !out.is_open() || (out.rdstate() & ifstream::failbit ) != 0 )
        return 3;
    
    const int size_hdr = 16;
    
    char hdr[size_hdr];
    *( (int*)(&(hdr[0])+0) )    = hdaf_order;
    *( (int*)(&(hdr[0])+4) )    = n_points;
    *( (double*)(&(hdr[0])+8) ) = step_size;
    
    out.write( &(hdr[0]), size_hdr );
    if ( (out.rdstate() & ofstream::failbit ) != 0  )
        return 4;
    
    out.write( (char *)kernel, n_points*8 );
    if ( (out.rdstate() & ofstream::failbit ) != 0  )
        return 5;
    
    out.close();
    
    return 0;
}





int read_std_hdaf_kernel_file( const char * file_name, int *hdaf_order, double *step_size, int *n_points, double **kernel )
{
/*
    
    file contents:
        int hdaf_order (4 bytes)
        int n_points (4 bytes)
        doulbe step_size (8 bytes)
        immediatly after that is the kernel data, which is an array of doubles (8 bytes) of length n_points.
    
    error codes:
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
    
    if (kernel==NULL || hdaf_order==NULL || n_points==NULL || step_size==NULL )
        return -1;
    
    ifstream in;
    in.open(file_name, ios::in | ios::binary );
    if ( !in.good() || !in.is_open() || (in.rdstate() & ifstream::failbit ) != 0 )
        return 1;
    
    const int size_hdr = 16;
    
    in.seekg(0, ios::end);
    int file_size = in.tellg();
    in.seekg(0, ios::beg);
    
    if ( file_size < size_hdr )
        return 2;
    
    char * buff = (char *)malloc(file_size);
    if (buff==NULL)
        return 3;
    
    in.read(buff, file_size);
    if ( (in.rdstate() & ifstream::failbit ) != 0  )
        return 4;
    in.close();
    
    *hdaf_order = *( (int*)(buff+0) );
    *n_points   = *( (int*)(buff+4) );
    *step_size  = *( (double*)(buff+8) );
    
    if ( *hdaf_order<0 || *n_points<=0 || *step_size<=0 || isinf(*step_size) || isnan(*step_size) || !(*step_size==*step_size) )
        return 5;
    
    int n_points_read = (file_size-size_hdr)/8;
    if ( n_points_read != *n_points )
        return 6;
    
    if (*kernel != NULL)
    {
        free( *kernel );
        kernel = 0;
    }
    
    *kernel = (double *)malloc( (*n_points)*8 );
    
    if (*kernel == NULL)
    {
        *n_points = 0;
        return 7;
    }
    
    for (int k=0; k<(*n_points); k++)
    {
        (*kernel)[k] = *( (double*)(buff+size_hdr+k*8) );
        if ( isinf((*kernel)[k]) || isnan((*kernel)[k]) || !((*kernel)[k]==(*kernel)[k]) )
            return 8;
    }
    
    free(buff);
    buff=0;
    
    return 0;
}


double hdaf_data_access( int k, double *kernel, int n_points )
{
    k = abs(k);
    if (k>=n_points) return 0.0;
    return kernel[k];
}

double hdaf_neville_Pij(double x, int i, int j, int n_data, double * data, double step_size )
{
	if (i == j)
		return hdaf_data_access( i, data, n_data );
    
    double DX = step_size*(i-j);
    if (fabs(DX) < 1E-15) DX = 1E-15;
    return ((x-step_size*j)*hdaf_neville_Pij(x,i,j-1,n_data,data,step_size) +(step_size*i-x)*hdaf_neville_Pij(x,i+1,j,n_data,data,step_size))/DX;
}

double hdaf_interpolate( double * kernel, int n_points, double step_size, double x )
{
    // linear interpolation for now
    
    x = fabs(x);
    
    int k=floor(x/step_size);
    if (k+1>=n_points)
        return 0.0;
    
    int L=5;
    
    return hdaf_neville_Pij( x, k-L, k+L+1, n_points, kernel, step_size );
    
    //if (k<interp_order) return (kernel[k+1]-kernel[k])*(x-k*step_size)/step_size + kernel[k];
}


int get_hdaf_kernel(double **kernel, int *kernel_size, double sampling_period, int order, double sigma, const char *hdaf_data_dir )
{
    if ( kernel == NULL || kernel_size==NULL )
        return -1;
    
    if ( sampling_period<= 1E-15 || order < 0 || sigma <=0.0 )
        return 1;
    
    char file_name[hdaf_max_file_name_length];
    std_kernel_file_naming_convention( file_name, hdaf_data_dir, order );
    
    int error=0;
    int read_order=0;
    double step_size=0;
    int n_points=0;
    double * std_kernel=0;
    
    error = read_std_hdaf_kernel_file( file_name, &read_order, &step_size,  &n_points, &std_kernel );
    
    if ( error )
        return 100+error;
    
    if ( read_order != order )
        return 2;
    
    if ( std_kernel==NULL )
        return 3;
    
    double s = 1.0/(hdaf_sqrt2*sigma);
    
    *kernel_size = (int)ceil((step_size*n_points)/(s*sampling_period));
    
    if (*kernel != NULL)
    {
        free( *kernel );
        *kernel = 0;
    }
    
    *kernel = (double *)malloc( (*kernel_size)*8 );
    
    if (*kernel == NULL)
    {
        *kernel_size = 0;
        return 4;
    }
    
    for (int k=0; k<*kernel_size; k++)
        (*kernel)[k] = s*hdaf_interpolate( std_kernel, n_points, step_size, (sampling_period*k)*s );
    
    // std kernel is not scaled by it's step size, doing so for returned kernel.
    
    free( std_kernel );
    std_kernel=0;
    
    return 0;
}

int get_hdaf_kernel_arbitrary_points(double *eval_points, double *kernel, int n_points, int order, double sigma, const char *hdaf_data_dir )
{
    if (n_points<=0 || kernel==NULL || eval_points==NULL )
        return 1;
    
    char file_name[hdaf_max_file_name_length];
    std_kernel_file_naming_convention( file_name, hdaf_data_dir, order );
    
    int error=0;
    int read_order=0;
    double step_size=0;
    int n_stdker_points=0;
    double * std_kernel=0;
    
    error = read_std_hdaf_kernel_file( file_name, &read_order, &step_size,  &n_stdker_points, &std_kernel );
    
    if (error)
        return 100+error;
    
    if ( read_order != order )
        return 2;
    
    if ( std_kernel==NULL )
        return 3;
        
    double s = 1.0/(hdaf_sqrt2*sigma);
    
    for (int k=0; k<n_points; k++)
        (kernel)[k] = s*hdaf_interpolate( std_kernel, n_stdker_points, step_size, eval_points[k]*s );
    
    free( std_kernel );
    std_kernel=0;
    
    return 0;
}






int get_hdaf_kernel_lp(double **kernel, int *kernel_size, double sampling_period, int order, double cutoff_frequency, const char *hdaf_data_dir )
{
    double sigma = sigma_from_cutoff_frequency( cutoff_frequency, order );
    int error=0;
    error = get_hdaf_kernel( kernel, kernel_size, sampling_period, order, sigma, hdaf_data_dir );
    
    if ( (*kernel != NULL) && *kernel_size >= 0 && error==0 )
    {
        for (int k=0; k<*kernel_size; k++ )
            (*kernel)[k] *= sampling_period;
    }
    
    return error;
}

int get_hdaf_kernel_bp(double **kernel, int *kernel_size, double sampling_period, int low_pass_order, double low_pass_frequency, int high_pass_order, double high_pass_frequency, const char *hdaf_data_dir )
{
    if ( kernel == NULL || kernel_size==NULL )
        return -1;
    
    double *low_pass_ker=0;
    double *high_pass_ker=0;
    int low_pass_ker_size=0;
    int high_pass_ker_size=0;
    
    int low_pass_err=0;
    int high_pass_err=0;
    
    low_pass_err = get_hdaf_kernel_lp( &low_pass_ker, &low_pass_ker_size, sampling_period, low_pass_order, low_pass_frequency, hdaf_data_dir );
    if (low_pass_err)
        return 100+low_pass_err;
    
    high_pass_err = get_hdaf_kernel_lp( &high_pass_ker, &high_pass_ker_size, sampling_period, high_pass_order, high_pass_frequency, hdaf_data_dir );
    if (high_pass_err)
        return 200+high_pass_err;
    
    // always: band_pass = low_pass - high_pass
    
    if (low_pass_ker_size >= high_pass_ker_size)
    {
        for (int k=0; k<high_pass_ker_size; k++)
            low_pass_ker[k] -= high_pass_ker[k];
        
        *kernel = low_pass_ker;
        *kernel_size = low_pass_ker_size;
        
        free(high_pass_ker);
        high_pass_ker = 0;
    }
    
    else // high_pass_ker_size > low_pass_ker_size
    {
        for (int k=0; k<low_pass_ker_size; k++)
            high_pass_ker[k] = low_pass_ker[k] - high_pass_ker[k];
        
        for (int k=low_pass_ker_size; k<high_pass_ker_size; k++)
            high_pass_ker[k] = 0.0 - high_pass_ker[k];
        
        *kernel = high_pass_ker;
        *kernel_size = high_pass_ker_size;
        
        free(low_pass_ker);
        low_pass_ker = 0;
    }
    
    return 0;
}













 #ifdef __cplusplus
 }
 #endif

















/*
double neville_Pij(double x, int i, int j, int n_data, double * data, double step_size )
{
	assert(i>=0); assert(i<n_data); assert(j<n_data);assert(i<=j);
	if (i == j)
		return data[i];
	else // i!=j;
	{
		double DX = tstep_size*(i-j);
		if (fabs(DX) < 1E-14) DX = 1E-14;
		return (x-step_size*j)*neville_Pij(x,i,j-1,n_data,data,step_size) +(step_size*i-x)*neville_Pij(x,i+1,j,n_data,data,step_size))/DX;
	}
}*/
/*

double neville_interpolate( double x, int interp_order, int n_data, double * data, double step_size )
D signal1<T,D>::nev_interp(D t,int d) // neville's algorithm
{
	if (n_data <= interp_order) return neville_Pij(x,0,n_data-1,n_data,data,step_size);
    
	for (k = 0;k<N-1 && this->t(k) < t ;k++)
    
    
	if (k==N) {return Pij(t,N-d-1,N-1); } //extrap right
	if (k==0) {return Pij(t,0,d); } //extrap left
	int i,j;
	i = k-1;
	j = k;
	while (j-i<d)
	{
		if ( t-this->t(i-1) < this->t(j+1)-t  ) i--;
		else j++;
		if (i<0) { i=0;j++;}
		if (j>=N) {j = N-1;i--;}
	}
	return Pij(t,i,j);
}
*/







//#include <mathlib/math/std_math.h>
//#include <mathlib/math/hdaf/hdaf.h>

//#include <mathlib/link.cpp>
//#include <mathlib/non-lib_things.h>


/*

double compute_hdaf_kernel_truncation(int m, double sigma, double eps_min, double eps_max)
{
    return hdaf_truncate_point(eps_min, eps_max, m, sigma, 0 );
}


int compute_hdaf_kernel_size(int m, double sigma, double sampling_period)
{
    return ceil( hdaf_truncate_point( m, sigma, 1E-17, 5E-17 )/sampling_period )+1;
}


*/



/*
double hdaf_fourier_transform(double omega, int m, double sigma)
{
    if (fabs(omega)<1E-12) return 1.0;
    double s = 0.0;
	double x = k*k*sigma*sigma/2.0;
	double r = log(fabs(omega))*2.0+log(sigma)*2.0;
	for (int n = 0; n<=m; n++)
		s += exp(-x+r*n-ml_log2*n-lognfact(n));
	return s;
}
*/







/*

class hdaf_delta_hat
{
public:
	double sigma;
	int m;
	dlognfact lognfact;
    
	hdaf_deltaHat(int m_, double sig, int mmax_ = 2000)
	:m(m_),sigma(sig),lognfact(mmax_) {};
	
	double operator() (double const & k) const 
	{
		if (m >= lognfact.nmax)
			*(const_cast<dlognfact *> (&lognfact)) = dlognfact(lognfact.nmax*2);
		
		double s = 0.0;
		double x = k*k*sigma*sigma/2.0;
		double r = log(fabs(k))*2.0+log(sigma)*2.0;
		if (fabs(k)<1E-12) return 1.0;
		for (int n = 0;n<=m;n++)
			s += exp(-x+r*n-ml_log2*n-lognfact(n));
		return s;
	}
};

*/






#endif
