


from hdaf import *
from pylab import *
from math import *

cat = numpy.concatenate

def nz(T, f):
    sum = 0.0
    for k in range(len(T)-1):
        sum += (T[k+1]-T[k])*f[k]
    return f/sum



f = 1
m = 32
h = 0.1
delta = lp_hdaf_kernel( h, m, f )

delta = cat((delta[::-1][0:len(delta)-1],delta))

X = [ h*(x-len(delta)/2)  for x in range(len(delta)) ]




delta = nz(X, delta)

var = 0.0

for k,x in enumerate(X):
    
    var += x**2*delta[k]*h

print sqrt(var) , 1.0/(freq_to_sigma(f , m)*pi)




plot(X,delta)



show()








