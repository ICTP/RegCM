#!/usr/bin/python
# Configuration script for the RegCM distribution

import os,sys

def main():
    print "Welcome to the RegCM configuration!"
    print "BETA VERSION"

    print "Please enter the path to your RegCM distribution [",os.getcwd(),"]"
    regcm_root=raw_input("regcm_root = ")
    # will need a more sophisticated test here!
    if not os.path.isdir(regcm_root):
        regcm_root = os.getcwd()

    print "Please enter the path where the RegCM binaries will be stored [",regcm_root+"/bin","]"
    bin_dir=raw_input("bin_dir = ")
    if not os.path.isdir(bin_dir):
        bin_dir=regcm_root+"/bin"

    ncpath=netcdf_search()

    print "Do you want a debug - 0 or a production - 1  binary [1]?"
    dbg=raw_input("dbg = ")
    if len(dbg) == 0 :
        dbg = 1

    print "Do you want a serial - 0 or MPI-parallel - 1 binary [1]?"
    mpi=raw_input("mpi = ")
    if len(mpi) == 0 :
        mpi = 1

    print "The following compiler/architecture combinations are tested and known to work,\nplease choose one:\n\n\t1. GNU Fortran v. 4.4 (Linux x86-64) [default]\n\t2. Intel Fortran v. 10 or 11 (Linux x86-64)\n\t3. PGI Fortran v. 9 (Linux x86-64)\n\t4. IBM Xlf Compiler (AIX on SPx)\n\t5. Other"
    compiler=raw_input("\ncompiler = ")
    if compiler == 0 :
        compiler = 1
        
    print regcm_root
    print bin_dir
    print ncpath
    print dbg
    print mpi
    print compiler

def netcdf_search():

    ncok=True

    # check for path in env variables
    try :
        path_a=os.environ['NETCDF']
        if (os.path.isfile(path_a+"/lib/libnetcdf.a") or os.path.isfile(path_a+"/lib/libnetcdf.so")) and (os.path.isfile(path_a+"/include/netcdf.mod")):
            ncpath=path_a
        else :
            ncok=False
    except KeyError :
        ncok=False

    if not ncok:
        try :
            path_b=os.environ['NETCDF_LIB']
            #print path_b[:-3]+"include"
            if (os.path.isfile(path_b+"/libnetcdf.a") or os.path.isfile(path_b+"/libnetcdf.so")) and (os.path.isfile(path_b[:-3]+"/include/netcdf.mod")):
                ncpath=path_b[:-3]
                ncok=True
                #print ncpath
            else :
                ncok=False
        except KeyError :
            ncok=False
    #print ncpath
    if ncok:
        return ncpath
    else:
        print "Unable to find a working NetCDF... please input a valid path."
        ncpath=raw_input("ncpath = ")
        if (os.path.isfile(ncpath+"/lib/libnetcdf.a") or os.path.isfile(ncpath+"/lib/libnetcdf.so")) and (os.path.isfile(ncpath+"/include/netcdf.mod")):
            print "NetCDF found in ",ncpath
            return ncpath
        else :
            print "NetCDF not found! Cannot continue. Please provide a working NetCDF library."
            os.sys.exit(1)

if __name__ == "__main__":
    main()
