
#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "/workspace/hdaf_lib/hdaf.h"


int main()
{
    int n1 = 256;
    int n2 = 256;
    double L1 = 1.0;
    double L2 = 1.0;
    
    void *data=0;
    
    double *in= (double*)malloc(8*n1*n2);
    double *out= (double*)malloc(8*n1*n2);
    
    
    laplacian2d_init(&data, n1, n2, L1, L2, 8,8, 0.8,0.8 );
    
    laplacian2d_execute( data, in, out );
    
    laplacian2d_free( &data );
    
}




