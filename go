cmake .
make
make install DESTDIR=/data/apps/armadillo/svn/
rm -rf /data/apps/armadillo/svn/lib
rm -rf /data/apps/armadillo/svn/include
mv /data/apps/armadillo/svn/usr/* /data/apps/armadillo/svn/
rm -r /data/apps/armadillo/svn/usr
