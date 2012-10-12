#!/bin/bash
#
# INSTALL NETCDF V3 LIBRARY IN DEST.
#
# Requires wget program
#
# CHEK HERE BELOW THE COMPILERS
#
# Working CC Compiler
CC=icc
# Working Fortran Compiler
FC=ifort
# Destination directory
DEST=$PWD

if [ -z "$DEST" ]
then
  echo "SCRIPT TO INSTALL NETCDF V3 LIBRARY"
  echo "EDIT ME TO DEFINE DEST, CC AND FC VARIABLE"
  exit 1
fi

UNIDATA=http://www.unidata.ucar.edu/downloads/netcdf/ftp

WGET=`which wget 2> /dev/null`
if [ -z "$WGET" ]
then
  echo "wget programs must be installed to download netCDF lib."
  exit 1
fi

echo "This script installs the netCDF library in the $DEST directory."
echo "If something goes wrong, logs are saved in $DEST/logs"

cd $DEST
mkdir $DEST/logs
echo "Downloading netCDF Library..."
$WGET $UNIDATA/netcdf-4.2.1.1.tar.gz -o $DEST/logs/download_C.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF C library from www.unidata.ucar.edu"
  exit 1
fi
$WGET $UNIDATA/netcdf-fortran-4.2.tar.gz -o $DEST/logs/download_F.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF Fortran library from www.unidata.ucar.edu"
  exit 1
fi

echo "Compiling netCDF Library."
tar zxvf netcdf-4.2.1.1.tar.gz > $DEST/logs/extract.log
cd netcdf-4.2.1.1
./configure CC="$CC" FC="$FC" --prefix=$DEST --disable-netcdf-4 \
  --disable-shared --disable-dap > $DEST/logs/configure.log 2>&1
make > $DEST/logs/compile.log 2>&1 && make install > $DEST/logs/install.log
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF C library"
  exit 1
fi
cd $DEST
rm -fr netcdf-4.2.1.1*
echo "Compiled netCDF C library. Half way..."

tar zxvf netcdf-fortran-4.2.tar.gz >> $DEST/logs/extract.log
cd netcdf-fortran-4.2
sed -i configure -e 's/ nc_def_opaque//'
PATH=$DEST/bin:$PATH ./configure CC="$CC" FC="$FC" \
     CPPFLAGS=-I$DEST/include LD_FLAGS=-L$DEST/lib --prefix=$DEST \
     --disable-shared >> $DEST/logs/configure.log 2>&1
make >> $DEST/logs/compile.log 2>&1 && make install >> $DEST/logs/install.log
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF Fortran library"
  exit 1
fi
cd $DEST
rm -fr netcdf-fortran-4.2*

# Done
echo "Done!"
echo "To link RegCM with this library use:"
echo
echo  --with-netcdf=$DEST
echo

echo "Cleanup..."
rm -fr $DEST/share $DEST/logs $DEST/lib/pkgconfig $DEST/lib/*.la
exit 0
