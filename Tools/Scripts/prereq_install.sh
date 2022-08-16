#!/bin/bash
#
# INSTALL NETCDF and MPI2 LIBRARY IN DEST.
#
# Requires wget program
#
# CHEK HERE BELOW THE COMPILERS
#
# Working CC Compiler
#CC="gcc -fPIC"
#CC="icc -fPIC"
#CC="pgcc -fpic -DpgiFortran"
#CC="xlc_r -qpic"
# Working C++ Compiler
#CXX="g++ -fPIC"
#CXX="icpc -fPIC"
#CXX="pgCC -fpic"
#CXX="xlc++_r -qpic"
# Working Fortran Compiler
#FC="gfortran -fPIC"
#FC="ifort -fPIC"
#FC="pgf90 -fpic"
#FC="xlf2003_r -qpic"
# Destination directory
#DEST=$PWD

netcdf_c_ver=4.9.0
netcdf_f_ver=4.5.4
hdf5_ver=1.13.1
zlib_ver=1.2.12
ompi_ver=4.1.4
ompi_major=`echo $ompi_ver | cut -d "." -f 1-2`
hdf5_major=`echo $hdf5_ver | cut -d "." -f 1-2`

UNIDATA=https://downloads.unidata.ucar.edu/
OPENMPI=http://www.open-mpi.org/software/ompi/v${ompi_major}/downloads
HDFGROUP=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_major}/hdf5-${hdf5_ver}/src
ZLIB=http://zlib.net

export LD_LIBRARY_PATH=$DEST/lib:$LD_LIBRARY_PATH

if [ -z "$DEST" ]
then
  echo "SCRIPT TO INSTALL NETCDF V5 and MPICH LIBRARIES."
  echo "EDIT ME TO DEFINE DEST, CC AND FC VARIABLE"
  exit 1
fi

MAKE=`which gmake 2> /dev/null`
if [ -z "$MAKE" ]
then
  echo "Assuming make program is GNU make program"
  MAKE=make
fi

WGET=`which wget 2> /dev/null`
if [ -z "$WGET" ]
then
  echo "wget programs must be installed to download netCDF lib."
  exit 1
fi

WGET="$WGET --no-check-certificate"

echo "This script installs the netCDF/mpi librares in the"
echo
echo -e "\t $DEST"
echo
echo "directory. If something goes wrong, logs are saved in"
echo
echo -e "\t$DEST/logs"
echo

cd $DEST
mkdir -p $DEST/logs
echo "Downloading ZLIB library..."
$WGET -c $ZLIB/zlib-${zlib_ver}.tar.gz -o $DEST/logs/download_Z.log
if [ $? -ne 0 ]
then
  echo "Error downloading ZLIB library from zlib.net"
  exit 1
fi
echo "Downloading HDF5 library..."
$WGET -c $HDFGROUP/hdf5-${hdf5_ver}.tar.gz \
      -o $DEST/logs/download_H.log
if [ $? -ne 0 ]
then
  echo "Error downloading HDF5 library from www.hdfgroup.org"
  exit 1
fi
echo "Downloading netCDF Library..."
$WGET -c $UNIDATA/netcdf-c/${netcdf_c_ver}/netcdf-c-${netcdf_c_ver}.tar.gz -o $DEST/logs/download_C.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF C library from www.unidata.ucar.edu"
  exit 1
fi
$WGET -c $UNIDATA/netcdf-fortran/${netcdf_f_ver}/netcdf-fortran-${netcdf_f_ver}.tar.gz -o $DEST/logs/download_F.log
if [ $? -ne 0 ]
then
  echo "Error downloading netCDF Fortran library from www.unidata.ucar.edu"
  exit 1
fi
echo "Downloading OPENMPI Library..."
$WGET -c $OPENMPI/openmpi-${ompi_ver}.tar.bz2 -o $DEST/logs/download_M.log
if [ $? -ne 0 ]
then
  echo "Error downloading OPENMPI from OPENMPI website"
  exit 1
fi

rm -f logs/*.log

echo "Compiling zlib Library."
tar zxvf zlib-${zlib_ver}.tar.gz >> $DEST/logs/extract.log
if [ $? -ne 0 ]
then
  echo "Error uncompressing zlib library"
  exit 1
fi
cd zlib-${zlib_ver}
echo CC="$CC" FC="$FC" ./configure --prefix=$DEST --shared >> \
	$DEST/logs/configure.log
CC="$CC" FC="$FC" ./configure --prefix=$DEST --shared >> \
             $DEST/logs/configure.log 2>&1
$MAKE >> $DEST/logs/compile.log 2>&1 && \
  $MAKE install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling zlib library"
  exit 1
fi
cd $DEST
rm -fr zlib-${zlib_ver}
echo "Compiled zlib library."
echo "Compiling MPI library."
tar jxvf openmpi-${ompi_ver}.tar.bz2 >> $DEST/logs/extract.log
if [ $? -ne 0 ]
then
  echo "Error uncompressing openmpi library"
  exit 1
fi
cd openmpi-${ompi_ver}
echo ./configure CC="$CC" FC="$FC" F77="$FC" CXX="$CXX" \
      --prefix=$DEST --disable-cxx >> $DEST/logs/configure.log
./configure CC="$CC" FC="$FC" F77="$FC" CXX="$CXX" \
	--prefix=$DEST >> $DEST/logs/configure.log 2>&1
$MAKE >> $DEST/logs/compile.log 2>&1 && \
  $MAKE install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling openmpi library"
  exit 1
fi
cd $DEST
rm -fr openmpi-${ompi_ver}
echo "Compiled MPI library."
echo "Compiling HDF5 library."
tar zxvf hdf5-${hdf5_ver}.tar.gz >> $DEST/logs/extract.log
cd hdf5-${hdf5_ver}
echo ./configure CC="$CC" CXX="$CXX" FC="$FC" \
        --prefix=$DEST --with-zlib=$DEST --enable-shared \
        --disable-cxx --disable-fortran >> $DEST/logs/configure.log
./configure CC="$CC" CXX="$CXX" FC="$FC" \
	--prefix=$DEST --with-zlib=$DEST --enable-shared \
        --disable-cxx --disable-fortran >> $DEST/logs/configure.log 2>&1
$MAKE >> $DEST/logs/compile.log 2>&1 && \
  $MAKE install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling HDF5 library"
  exit 1
fi
cd $DEST
rm -fr hdf5-${hdf5_ver}
echo "Compiled HDF5 library."
echo "Compiling netCDF Library."
tar zxvf netcdf-c-${netcdf_c_ver}.tar.gz >> $DEST/logs/extract.log
cd netcdf-c-${netcdf_c_ver}
H5LIBS="-lhdf5_hl -lhdf5 -lz"
if [ "X$FC" == "Xgfortran" ]
then
  H5LIBS="$H5LIBS -lm -ldl"
fi
echo ./configure CC="$CC" FC="$FC" --prefix=$DEST --enable-netcdf-4 \
  CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib LIBS="$H5LIBS" \
  --enable-dap --enable-shared >> $DEST/logs/configure.log
./configure CC="$CC" FC="$FC" --prefix=$DEST --enable-netcdf-4 \
  CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib LIBS="$H5LIBS" \
  --enable-dap --enable-shared >> $DEST/logs/configure.log 2>&1
$MAKE >> $DEST/logs/compile.log 2>&1 && \
  $MAKE install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF C library"
  exit 1
fi
cd $DEST
rm -fr netcdf-c-${netcdf_c_ver}
echo "Compiled netCDF C library."
tar zxvf netcdf-fortran-${netcdf_f_ver}.tar.gz >> $DEST/logs/extract.log
cd netcdf-fortran-${netcdf_f_ver}
echo ./configure PATH=$DEST/bin:$PATH CC="$CC" FC="$FC" F77="$FC" \
     CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib --prefix=$DEST \
     --enable-shared >> $DEST/logs/configure.log
./configure PATH=$DEST/bin:$PATH CC="$CC" FC="$FC" F77="$FC" \
     CPPFLAGS=-I$DEST/include LDFLAGS=-L$DEST/lib --prefix=$DEST \
     --enable-shared >> $DEST/logs/configure.log 2>&1
$MAKE >> $DEST/logs/compile.log 2>&1 && \
  $MAKE install >> $DEST/logs/install.log 2>&1
if [ $? -ne 0 ]
then
  echo "Error compiling netCDF Fortran library"
  exit 1
fi
cd $DEST
rm -fr netcdf-fortran-${netcdf_f_ver}
echo "Compiled netCDF Fortran library."

# Done
CC=`echo $CC | cut -d " " -f 1`
FC=`echo $FC | cut -d " " -f 1`
echo
echo                 "Done!"
echo
echo "To link RegCM with this librares use:"
echo
echo  "./configure PATH=$DEST/bin:\$PATH \\"
echo  "            CC=\"$CC\" \\"
echo  "            FC=\"$FC\" \\"
echo  "            MPIFC=\"$DEST/bin/mpif90\" \\"
echo  "            CPPFLAGS=-I$DEST/include \\"
echo  "            LDFLAGS=-L$DEST/lib \\"
echo  "            LIBS=\"-lnetcdff -lnetcdf $H5LIBS\""
echo
echo "To run the model use these PATH and LD_LIBRARY_PATH variable:"
echo
echo "export LD_LIBRARY_PATH=$DEST/lib:\$LD_LIBRARY_PATH"
echo "export PATH=$DEST/bin:\$PATH"
echo
echo or
echo
echo "setenv LD_LIBRARY_PATH $DEST/lib:\$LD_LIBRARY_PATH"
echo "setenv PATH $DEST/bin:\$PATH"
echo "rehash"
echo

echo "Cleanup..."
mkdir -p src && mv -f *gz *.bz2 src || exit 1
exit 0
