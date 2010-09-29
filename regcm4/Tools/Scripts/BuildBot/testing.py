#!/usr/bin/python

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Script for running RegCM regression tests by M. Scarcia
#

import os,sys,shutil,fileinput,subprocess

def edit_namelist(namelist,datadir,simdir):

    for line in fileinput.FileInput(namelist,inplace=1):

        line = line.replace("/set/this/to/where/your/surface/dataset/is",datadir+"/")
        line = line.replace("/set/this/to/where/your/input/global/data/is",datadir+"/")

        line = line.replace("/set/this/to/where/your/domain/file/is",simdir+"/input")
        line = line.replace("/set/this/to/where/your/icbc/for/model/is",simdir+"/input")

        line = line.replace("/set/this/to/where/your/output/files/will/be/written",simdir+"/output")

        line = line.replace("/set/this/to/where/your/input/clm/data/are",simdir+"/input")

        print line.rstrip()
        
    fileinput.close()

COMMENT_CHAR = '#'
OPTION_CHAR =  '='
 
def parse_config(filename):
    options = {}
    try:
        f = open(filename)
    except :
        print "File "+filename+" does not exist or is not accessible!"
        sys.exit(1)
        
    for line in f:
        # First, remove comments:
        if COMMENT_CHAR in line:
            # split on comment char, keep only the part before
            line, comment = line.split(COMMENT_CHAR, 1)
        # Second, find lines with an option=value:
        if OPTION_CHAR in line:
            # split on option char:
            option, value = line.split(OPTION_CHAR, 1)
            # strip spaces:
            option = option.strip()
            value = value.strip()
            # store in dictionary:
            options[option] = value
    f.close()
    return options

def compare_nc_file(filename,refname,varname):

    #print filename
    #print refname
    
    try :
        p_1 = subprocess.Popen("ncdiff -v "+varname+" "+filename+" "+refname+" temp.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
        #print "ncdiff ok"
    except OSError :
        print "Could not run ncdiff!"
        output,error = p_1.communicate()
        return output

    if p_1.wait() == 0 :
        try :
            p_2 = subprocess.Popen("ncwa -y rms temp.nc rms.nc && rm temp.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
        except OSError :
           print "Could not run ncwa!"
           output,error = p_2.communicate()
           return output  
    else :
        print "Step 1 failed!"
        output,error = p_1.communicate()
        return output

    if p_2.wait() == 0 :
        try :
            p_3 = subprocess.Popen('ncks -H -s "%g\n" -v '+varname+' rms.nc && rm rms.nc',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        except OSError :
            print "Could not run ncks!"
            output,error = p_3.communicate()
            return output
    else :
        print "Step 2 failed!"
        output,error = p_2.communicate()
        return output    
    
    if p_3.wait() != 0:
        print "Step 3 failed!"
        output,error = p_3.communicate()
        return output+error
    else:
        output,error = p_3.communicate()
        
    return output
                
def main(argv):

    if (len(sys.argv) < 2):
        print "Please specify a configuration file!"
        sys.exit(1)
        
    cfg=sys.argv[1]

    options = parse_config(cfg)

    datadir = options["DATADIR"]
    bindir = options["BINDIR"]
    testdir = options["TESTDIR"]
    namelistdir = options["NLDIR"]
    referencedir = options["REFDIR"]
    teststodo = options["TESTSTODO"]

    datadir = os.path.abspath(datadir)
    testdir = os.path.abspath(testdir)
    bindir = os.path.abspath(bindir)
    referencedir = os.path.abspath(referencedir)
    
    if teststodo.rfind(",") > -1 :
        tests=teststodo.split(",")
        listtype=0
    elif teststodo.rfind("-") > -1:
        tests=teststodo.split("-")
        imin=int(tests[0])
        imax=int(tests[1])
        listtype=1
    else :
        tests=int(teststodo)
        listtype=2

    # main loop over tests
    # number of total tests present
    TOT_TESTS = 10

    if not os.path.isdir(testdir):
        os.mkdir(testdir)

    if listtype == 2 :
        if tests == 0 :
            imin = 1
            imax = TOT_TESTS
        else :
            imin = int(tests)
            imax = int(tests)
    elif listtype == 1 :
        imin = int(tests[0])
        imax = int(tests[1])
    else :
        imin = 0
        imax = len(tests)-1


    #print "imin =",imin
    #print "imax =",imax
        
    for i in range(imin,imax+1):

        if listtype == 0 :
            testname="test_00"+str(tests[i])
        else :
            # make this better!
            testname="test_00"+str(i)

        simdir=testdir+"/"+testname
        testrefdir=referencedir+"/"+testname

        if not os.path.isdir(simdir):
            os.mkdir(simdir)

        # create simulation directory tree
        if not os.path.isdir(simdir+"/input"):
            os.mkdir(simdir+"/input")
        if not os.path.isdir(simdir+"/output"):
            os.mkdir(simdir+"/output")

        namelist = simdir+"/regcm.in"

        shutil.copy(namelistdir+"/"+testname+".in",namelist)

        # parse for idate0 and idate2

        infile = open(namelist,"r")
        file_content = infile.readlines()
        infile.close()
        
        for line in file_content:
            if line.find("globidate1") > -1 :
                linea=line.rsplit("=")
                globidate1=linea[1].rstrip(" ,")
            if line.find(" idate0 ") > -1 :
                linea=line.rsplit("=")
                idate0=linea[1]
                idate0=filter(lambda x:x.isdigit(),idate0)
            if line.find(" idate1 ") > -1 :
                linea=line.rsplit("=")
                idate1=linea[1].rstrip(" ,")
            if line.find(" idate2 ") > -1 :
                linea=line.rsplit("=")
                idate2=linea[1].rstrip(" ,")

        #print "a"+idate0+"a"
                
        #edit the namelist here

        edit_namelist(namelist,datadir,simdir)

        mpistring=options["MPISTRING"]
        run_preproc=int(options["PREPROC"])
        run_clm=int(options["USECLM"])
        run_band=int(options["USEBAND"])

        # open log file

        writelog=True

        try:
            log = open(testname+".log","w")
        except :
            print "Unable to write log!"
            writelog=False

        exit_status = 0

        runMain = True # won't run Main if PreProc crashes...
        
        # run preproc

        if (run_preproc == 1):
        
            p_terrain = subprocess.Popen(bindir+"/terrain "+namelist,stdout=log,stderr=log,shell=True)
            if p_terrain.wait() != 0:
                print "\nError: Terrain in",testname,"crashed!!\n"
                runMain = False
                exit_status = 1
            else:
                print "Terrain in",testname,"passed."
    
            p_sst=subprocess.Popen(bindir+"/sst "+namelist,stdout=log,stderr=log,shell=True)
            if p_sst.wait() != 0:
                print "\nError: SST in",testname,"crashed!!\n"
                runMain = False
                exit_status = 1
            else :
                print "SST in",testname,"passed."
            
            p_icbc=subprocess.Popen(bindir+"/icbc "+namelist,stdout=log,stderr=log,shell=True)
            if p_icbc.wait() != 0:
                print "\nError: ICBC in",testname,"crashed!!\n"
                runMain = False
                exit_status = 1
            else :
                print "ICBC in",testname,"passed."

            if run_clm == 1:
                print bindir+"/clm2rcm "+namelist
                p_clmpre=subprocess.Popen(bindir+"/clm2rcm "+namelist)
                print bindir+"/clm2rcm "+namelist
                if p_clmpre.wait() != 0:
                    print "\nError: clm2rcm in",testname,"crashed!!\n"
                    runMain = False
                    exit_status = 1
                else :
                    print "clm2rcm in",testname,"passed."

            sys.stdout.flush()

            # compare preproc output only if everything went ok
            if runMain :

                dom_diff={}
                icbc_diff={}

                domain_file = "/input/"+testname+"_DOMAIN000.nc"
                icbc_file = "/input/"+testname+"_ICBC."+idate0+".nc"

                domain_vars = ["topo","landuse"]
                icbc_vars = ["u","v","t","ts"]

                # domain
                for var in domain_vars :    
                    dom_diff[var] = compare_nc_file(simdir+domain_file,testrefdir+domain_file,var).rstrip("\n")
                    print var+" =",dom_diff[var]

                # icbc

                for var in icbc_vars :
                    icbc_diff[var] = compare_nc_file(simdir+icbc_file,testrefdir+icbc_file,var).rstrip("\n")
                    print var+" =",icbc_diff[var]
            
        #runMain = False

        # if preproc is ok, run main
        if runMain :
            if (int(options["SERIAL"]) == 1):
                p_regcm=subprocess.Popen(bindir+"/regcmSerial "+namelist,stdout=log,stderr=log,shell=True)
            elif int(run_clm == 1):
                p_regcm=subprocess.Popen(mpistring+" "+bindir+"/regcm_clM "+namelist,stdout=log,stderr=log,shell=True)
            elif int(run_band == 1):
                p_regcm=subprocess.Popen(mpistring+" "+bindir+"/regcm_band "+namelist,stdout=log,stderr=log,shell=True)
            else :
                p_regcm=subprocess.Popen(mpistring+" "+bindir+"/regcmMPI "+namelist,stdout=log,stderr=log,shell=True)
                
            if p_regcm.wait() != 0:
                print "\nError: RegCM",testname,"crashed!!\n"
                exit_status = 1
            else :
                print "RegCM",testname,"passed."
        else :
            print "PreProc did not complete correctly, RegCM main skipped..."

        log.close()

        if exit_status == 1:
            outlog = open(testname+".log","r")
            stdouterr = outlog.read()
            print stdouterr

        # if everything ok compare output
        
        if exit_status == 0:

            srf_diff={}
            srf_file="/output/"+testname+"_SRF."+idate0+".nc"

            srf_vars = ["t2m"]

            for var in srf_vars :
                srf_diff[var] = compare_nc_file(simdir+srf_file,testrefdir+srf_file,var).rstrip("\n")
                print var+" =",srf_diff[var]
            

	sys.stdout.flush()

        

    # end of the big loop
    if exit_status == 1:
        print "Warning! Some tests failed!"

    print "\n****  Test script terminated.  ****"

    sys.exit(exit_status)
    
if __name__ == "__main__":
    main(sys.argv[1:])

