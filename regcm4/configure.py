#!/usr/bin/python
# Configuration script for the RegCM distribution by M. Scarcia

import os,sys,shutil,fileinput

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
        bin_dir=regcm_root+"/Bin"

    ncpath=netcdf_search()

    print "Do you want a debug - 0 or a production - 1  binary [1]?"
    dbg=raw_input("dbg = ")
    if len(dbg) == 0 :
        dbg = 1

    print "Do you want a serial - 0 or MPI-parallel - 1 binary [1]?"
    mpi=raw_input("mpi = ")
    if len(mpi) == 0 :
        mpi = 1

    # Let's assume that the compiler is always mpif90
    # We can change that later
    if mpi == 1 :
        mpi_compiler="mpif90"

    print "Do you want to enable the DCSST  - 1 or not - 0 [0]?"
    dcsst=raw_input("DCSST = ")
    if len(dcsst) == 0 :
        dcsst = 0

    print "Do you want to enable the SeaIce  - 1 or not - 0 [0]?"
    seaice=raw_input("SeaIce = ")
    if len(seaice) == 0 :
        seaice = 0

    print "Do you want to enable CLM  - 1 or not - 0 [0]?"
    clm=raw_input("CLM = ")
    if len(clm) == 0 :
        clm = 0
    
    isnotint = True

    print "The following compiler/architecture combinations are tested and known to work,\nplease choose one:\n\n\t1. GNU Fortran v. 4.4 (Linux x86-64)\n\t2. Intel Fortran v. 10 or 11 (Linux x86-64)\n\t3. PGI Fortran v. 9 (Linux x86-64)\n\t4. IBM Xlf Compiler (AIX on SPx)\n\t5. Other"
    while isnotint :
        try :
            compiler=int(raw_input("\ncompiler = "))
            if (1 <= compiler <= 5):
                isnotint=False
            else :
                print "Please choose one of the available!"
        except:
            print "Please choose one of the available!" 

    # Add a non-blocking check to see if the compiler binaries actually exist

    print "\nChosen configuration :"
    print regcm_root
    print bin_dir
    print ncpath
    print dbg
    print mpi
    print dcsst
    print seaice
    print clm
    print compiler

    # Start real stuff

    choose_template(regcm_root,compiler,dbg)
    makefile_edit(regcm_root,bin_dir,ncpath,mpi,mpi_compiler)

    print "Configuration complete! You should be able to compile RegCM"

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

def choose_template(regcm_root,compiler,dbg):
    if dbg == 1 and compiler < 5 :
        if compiler == 1 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_gnu4.4",regcm_root+"/Makefile.inc")
        elif compiler == 2 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_intel10+",regcm_root+"/Makefile.inc")
        elif compiler == 3 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_pgi9",regcm_root+"/Makefile.inc")
        elif compiler == 4 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_xlf",regcm_root+"/Makefile.inc")
    elif dbg == 0 and compiler < 5 :
        if compiler == 1 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_gnu4.4_debug",regcm_root+"/Makefile.inc")
        elif compiler == 2 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_intel10+_debug",regcm_root+"/Makefile.inc")
        elif compiler == 3 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_pgi9_debug",regcm_root+"/Makefile.inc")
        elif compiler == 4 :
            shutil.copyfile(regcm_root+"/Arch/Makefile.inc_xlf_debug",regcm_root+"/Makefile.inc")
    else :
        shutil.copyfile(regcm_root+"/Arch/Makefile.inc_other",regcm_root+"/Makefile.inc")
        
    return


def makefile_edit(regcm_root,bin_dir,ncpath,mpi,mpi_compiler) :
    #fp = open(regcm_root+"/Makefile.inc","rU")

    for line in fileinput.FileInput(regcm_root+"/Makefile.inc",inplace=1):
        line=line.replace("!REGCM_ROOT",regcm_root)
        line=line.replace("!BIN_DIR",bin_dir)

        line=line.replace("!NETCDFINC","-I"+ncpath+"include")
        line=line.replace("!NETCDFLIB","-L"+ncpath+"lib -lnetcdf")
        line=line.replace("!NETCDFC++","-L"+ncpath+"lib -lnetcdf_c++ -lnetcdf")

        if mpi == 1 :
            line=line.replace("# PARALLEL = MPP1","PARALLEL = MPP1")
            line=line.replace("!MPIF90",mpi_compiler)
    
        print line.rstrip()

    fileinput.close()

if __name__ == "__main__":
    main()
