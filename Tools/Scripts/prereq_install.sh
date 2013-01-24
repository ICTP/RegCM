#!/bin/bash
#
# INSTALL NETCDF and MPI2 LIBRARY IN DEST.
#
# Requires wget program
#
# CHEK HERE BELOW THE COMPILERS
#
# Working CC Compiler
#CC=icc
# Working Fortran Compiler
#FC=ifort
# Destination directory
#DEST=$PWD

if [ -z "$DEST" ]
then
  echo "SCRIPT TO INSTALL NETCDF V4 and MPICH LIBRARIES."
  echo "EDIT ME TO DEFINE DEST, CC AND FC VARIABLE"
  exit 1
fi

UNIDATA=http://www.unidata.ucar.edu/downloads/netcdf/ftp
MPICH=http://www.mpich.org/static/tarballs/3.0.1

WGET=`which wget 2> /dev/null`
if [ -z "$WGET" ]
then
  echo "wget programs must be installed to download netCDF lib."
  exit 1
fi

echo "This script installs the netCDF/mpi librares in the"
echo
echo -e "\t $DEST"
echo
echo "directory. If something goes wrong, logs are saved in"
echo
echo -e "\t$DEST/logs"
echo

cd $DEST
mkdir $DEST/logs
echo "Downloading ZLIB library..."
$WGET -c http://zlib.net/zlib-1.2.7.tar.gz -o $DEST/logs/download_Z.log
if [ $? -ne 0 ]
then
  echo "Error downloading ZLIB library from zlib.net"
  exit 1
fi
echo "Downloading HDF5 library..."
$WGET -c http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.9.tar.gz \
      -o $DEST/logs/download_H.log
if [ $? -ne 0 ]
then
  echo "Error downloading HDF5 library from www.hdfgroup.org"
  exit 1
fi
echo "Downloading netCDF Library..."
$WGET -c $UNIDATA/netcdf-4.2.1.1.tar.gz -o $DEST/logs/download_C.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF C library from www.unidata.ucar.edu"
  exit 1
fi
$WGET -c $UNIDATA/netcdf-fortran-4.2.tar.gz -o $DEST/logs/download_F.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF Fortran library from www.unidata.ucar.edu"
  exit 1
fi
echo "Downloading MPICH Library..."
wget -c $MPICH/mpich-3.0.1.tar.gz -o $DEST/logs/download_M.log
if [ $? -ne 0 ]
then
  echo "Error downloading MPICH from MPICH website"
  exit 1
fi

echo "Compiling MPI library."
tar zxvf mpich-3.0.1.tar.gz > /dev/null
if [ $? -ne 0 ]
then
  echo "Error uncompressing mpich library"
  exit 1
fi
cd mpich-3.0.1
CC="$CC" FC="$FC" ./configure --prefix=$DEST > $DEST/logs/configure.log 2>&1
make > $DEST/logs/compile.log 2>&1 && \
  make install > $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling mpich library"
  exit 1
fi
cd $DEST
rm -fr mpich-3.0.1
echo "Compiled MPI library."
echo "Compiling zlib Library."
tar zxvf zlib-1.2.7.tar.gz > /dev/null
if [ $? -ne 0 ]
then
  echo "Error uncompressing zlib library"
  exit 1
fi
cd zlib-1.2.7
CC="$CC" FC="$FC" ./configure --prefix=$DEST --static >> \
             $DEST/logs/configure.log 2>&1
make >> $DEST/logs/compile.log 2>&1 && \
  make install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling zlib library"
  exit 1
fi
cd $DEST
rm -fr zlib-1.2.7
echo "Compiled zlib library."
echo "Compiling HDF5 library."
tar zxvf hdf5-1.8.9.tar.gz > $DEST/logs/extract.log
cd hdf5-1.8.9
./configure CC="$CC" --prefix=$DEST --with-zlib=$DEST --disable-shared \
        --disable-cxx --disable-fortran >> $DEST/logs/configure.log 2>&1
make > $DEST/logs/compile.log 2>&1 && \
  make install > $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling HDF5 library"
  exit 1
fi
cd $DEST
rm -fr hdf5-1.8.9
echo "Compiled HDF5 library."
echo "Compiling netCDF Library."
tar zxvf netcdf-4.2.1.1.tar.gz > $DEST/logs/extract.log
cd netcdf-4.2.1.1
./configure CC="$CC" FC="$FC" --prefix=$DEST --enable-netcdf-4 \
  CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib LIBS="-lhdf5_hl -lhdf5 -lz" \
  --disable-shared --disable-dap >> $DEST/logs/configure.log 2>&1
make > $DEST/logs/compile.log 2>&1 && \
  make install > $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF C library"
  exit 1
fi
cd $DEST
rm -fr netcdf-4.2.1.1
echo "Compiled netCDF C library."
tar zxvf netcdf-fortran-4.2.tar.gz >> $DEST/logs/extract.log
cd netcdf-fortran-4.2
PATH=$DEST/bin:$PATH ./configure CC="$CC" FC="$FC" \
     CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib --prefix=$DEST \
     --disable-shared >> $DEST/logs/configure.log 2>&1
make >> $DEST/logs/compile.log 2>&1 && \
  make install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF Fortran library"
  exit 1
fi
cd $DEST
rm -fr netcdf-fortran-4.2
echo "Compiled netCDF Fortran library."

# Done
echo "Done!"
echo "To link RegCM with this librares use:"
echo
echo  PATH=$DEST/bin:$PATH ./configure CC=$CC FC=$FC '\'
echo         CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib '\'
echo         LIBS=\"-lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz\"
echo

echo "Cleanup..."
mkdir src && mv *gz src || exit 1
exit 0
