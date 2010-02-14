


from hdaf import *
from pylab import *

cat = numpy.concatenate


for f in [ float(x) for x in range(100)  ]:



    m = 24
    h = 0.005
    delta = lp_hdaf_kernel( h, 2*m, f )

    delta = cat((delta[::-1][0:len(delta)-1],delta))

    X = [ h*(x-len(delta)/2)  for x in range(len(delta)) ]


    plot(X,delta)



show()








