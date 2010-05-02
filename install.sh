
#icpc -Wall -fPIC -c hdaf.cpp -o libhdaf.o $std_link -lpng -lpngwriter -lz -lfreetype -lrt
#icpc -shared -Wl,-soname,libhdaf.so.1 -o libhdaf.so.1.0 libhdaf.o -lpng -lpngwriter -lz -lfreetype -lrt



icpc -Wall -fPIC -c hdaf.cpp -o libhdaf.o  `freetype-config --cflags` -lpng -lpngwriter -lz -lfreetype -lrt -pthread -lgsl -lgslcblas -larprec -I/usr/local/include -L/usr/local/lib -lz -lfftw3_threads -lfftw3  
icpc -shared -Wl,-soname,libhdaf.so.1 -o libhdaf.so.1.0 libhdaf.o `freetype-config --cflags` -lpng -lpngwriter -lz -lfreetype -lrt -pthread -lgsl -lgslcblas -larprec  -I/usr/local/include -L/usr/local/lib -lz -lfftw3_threads -lfftw3


cp libhdaf.so.1.0 /usr/local/lib/
cp libhdaf.so.1.0 /usr/lib/
cp libhdaf.so.1.0 /lib64/
rm libhdaf.o
rm libhdaf.so.1.0

ln -sf /usr/local/lib/libhdaf.so.1.0 /usr/local/lib/libhdaf.so
ln -sf /usr/local/lib/libhdaf.so.1.0 /usr/local/lib/libhdaf.so.1
ln -sf /usr/lib/libhdaf.so.1.0 /usr/lib/libhdaf.so
ln -sf /usr/lib/libhdaf.so.1.0 /usr/lib/libhdaf.so.1
ln -sf /lib64/libhdaf.so.1.0 /lib64/libhdaf.so
ln -sf /lib64/libhdaf.so.1.0 /lib64/libhdaf.so.1


ln -svf  /workspace/hdaf_lib/python_hdaf/hdaf.py  /usr/lib/python2.6/site-packages/hdaf.py


