
g++ -Wall -fPIC -c hdaf.cc -o libhdaf.o
g++ -shared -Wl,-soname,libhdaf.so.1 -o libhdaf.so.1.0   libhdaf.o
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


