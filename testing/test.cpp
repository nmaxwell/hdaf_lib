




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
    
    
    int order = 13;
    double sampling_period=0.2;
    double sigma= dble(1.0/mp_sqrt2);
    
    double * kernel=0;
    int kernel_size=0;
    
    int error=0;
    error = get_hdaf_kernel( &kernel, &kernel_size, sampling_period, order, sigma, "/workspace/hdaf_data" );
    
    if (error)
        cout << "get_hdaf_kernel error: " << error << endl;
    
    ml_poly<mp_real > P;
    make_hdaf_ml_poly(P, order );
    
    mp_real s = (mp_real)sampling_period/(mp_sqrt2*sigma);
    mp_real ss = s*s;
    mp_real f = pow(mp_sqrt2*sigma,-1)*sampling_period;
    
    for (int k=0; k<kernel_size; k++)
        cout << k << "\t" <<  log10( fabs( (kernel[k]-dble(exp(-ss*(k*k))*P(dble(s*k))*f ))   /kernel[k] ))  << "\t" << log10( fabs(kernel[k]) ) << endl;
    
    
    
    return 0;
    
}






