cmake .
make
TAG=mmf
make install DESTDIR=/data/apps/armadillo/${TAG}/
rm -rf /data/apps/armadillo/${TAG}/lib
rm -rf /data/apps/armadillo/${TAG}/include
rm -rf /data/apps/armadillo/${TAG}/share
mv /data/apps/armadillo/${TAG}/usr/local/* /data/apps/armadillo/${TAG}/
rm -r /data/apps/armadillo/${TAG}/usr
ln -shf /data/apps/armadillo/${TAG}/include/* /usr/local/include/
ln -shf /data/apps/armadillo/${TAG}/lib/* /usr/local/lib/
