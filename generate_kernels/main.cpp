
#include <mathlib/math/std_math.h>

#include <mathlib/math/hdaf/hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>


#include "../hdaf.h"


#define mp_0 ((mp_real)(0.0))
#define mp_1 ((mp_real)(1.0))
#define mp_2 ((mp_real)(2.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)
#define mp_sqrt2 (sqrt(mp_2))
#define mp_iu ( complex<mp_real > (mp_0, mp_1 ) )


int main()
{
    std_setup();
    
    
    
    char hdaf_data_dir[] = "/workspace/hdaf_data";
    int max_order = 100;
    double eps_min = 1E-17;
    double eps_max = 1E-18;
    double over_sample = 4.0;
    
    
    for (int order=0; order <= max_order; order++ )
    {
        mp::mp_init(70);
        
        mp_real sigma = mp_1/mp_sqrt2;
        
        double step_size = dble( mp_1/( (mp_2pi*sqrt((mp_real)(2*order+1))/sigma) *2.0*over_sample) );
        
        ml_poly<mp_real > P;
        make_hdaf_ml_poly(P, order );
        
        int n_points = (int)ceil(hdaf_truncate_point (eps_min, eps_max, order, dble(sigma), 0 )/step_size)+2;
        
        double * kernel = ml_alloc<double > ( n_points );
        
        for (int k=0; k<n_points; k++ )
        {
            kernel[k] = dble( exp(-step_size*step_size*(k*k))*P(dble(step_size*k))  ); 
            
            //kernel[k] = dble((s*k)*order);
        }
        
        // not scaling by step size here.
        
        char fname[200];
        std_kernel_file_naming_convention( fname, hdaf_data_dir, order );
        
        cout << order << "\t" << step_size << "\t" << n_points << "\t" << fname << endl;
        
        int werr=0;
        werr = write_std_hdaf_kernel_file( fname, order, step_size, n_points, kernel );
        if ( werr )
            cerr << "write error:\t" << werr << endl;
        
        double * kernel2=0;
        int order2;
        double h2;
        int n_points2;
        
        int rerr = 0;
        rerr = read_std_hdaf_kernel_file( fname, order2, h2, n_points2, kernel2 );
        if ( rerr )
            cerr << "read error:\t" << rerr << endl;
        
        if ( order != order2 )
            cout << "order error\n";
        if ( step_size != h2 )
            cout << "step_size error\n";
        if ( n_points != n_points2 )
            cout << "n_points error\n";
        
        int ker_err=0;
        for (int k=0; k<n_points; k++ )
            if (kernel[k] != kernel2[k])
            {
                ker_err = 1;
                cout << k << "\t" << kernel[k] << "\t" << kernel2[k] << endl;
                }    
        
        if (ker_err != 0)
            cout << "kernel data error\n";
        
        ml_free( kernel );
        free( kernel2);
        kernel = 0;
    }
    
    
    
    std_exit();
}




