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

import os,sys,shutil,fileinput
from optparse import OptionParser

class MyOptions(object):
    
    def __init__(self):
        regcm_root=os.getcwd()+"/../../../"
        usage = "usage: %s [-t test directory] [-b binaries directory] [-d data directory]" % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-t","--test-dir",dest="testdir",default=regcm_root+"/sandbox",
                      help="Directory where the tests will be run", metavar="DIR")
#        self.parser.add_option("-q", "--quiet",
#                      action="store_false", dest="verbose", default=True,
#                      help="Don't print status messages to stdout")
        self.parser.add_option("-b","--bin-dir", dest="bindir",default=regcm_root+"/Bin",
                  help="Directory where RegCM binaries are", metavar="DIR")
	self.parser.add_option("-d","--data-dir", dest="datadir",default=regcm_root+"/data",
                  help="Directory where RegCM preprocessing data is", metavar="DIR")

    def parse(self,args):
        
        (self.options, self.args) = self.parser.parse_args(args=args)
        self.datadir = self.options.datadir
        self.testdir = self.options.testdir
        self.bindir = self.options.bindir

#	if not os.path.isdir(self.testdir):
#            print "The directory",self.datadir,"does not exist or is not accessible!"
#            sys.exit(1)

        if not os.path.isdir(self.datadir):
            print "The directory",self.datadir,"does not exist or is not accessible!"
            sys.exit(1)

	if not os.path.isdir(self.bindir):
            print "The directory",self.datadir,"does not exist or is not accessible!"
            sys.exit(1)

        if not os.path.isdir(self.testdir):
            try :
                os.mkdir(self.testdir)
            except OSError:
                print "Cannot create",self.testdir,"directory!"
                sys.exit(1)

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
    
            
def main():

    regcm_root=os.getcwd()+"/../../../"

    options = MyOptions()
    options.parse(sys.argv[1:])

    datadir = options.datadir
    bindir = options.bindir
    testdir = options.testdir

    datadir = os.path.abspath(datadir)
    testdir = os.path.abspath(testdir)
    bindir = os.path.abspath(bindir)

    # main loop over tests
    # number of tests temporarily fixed
    total_tests=2

    if not os.path.isdir(testdir):
        os.mkdir(testdir)

    for i in range(1,total_tests+1):

        # make this better!
        testname="test_00"+str(i)

        simdir=testdir+"/"+testname

        if not os.path.isdir(simdir):
            os.mkdir(simdir)

        # create simulation directory tree
        if not os.path.isdir(simdir+"/input"):
            os.mkdir(simdir+"/input")
        if not os.path.isdir(simdir+"/output"):
            os.mkdir(simdir+"/output")

        namelist = simdir+"/regcm.in"
        
        shutil.copy(regcm_root+"Testing/"+testname+".in",namelist)
        
        #os.symlink(bindir,simdir+"/Bin")
        
        #edit the namelist here

        edit_namelist(namelist,datadir,simdir)

        # hardcoded machine params

        mpistring="mpirun -np 4"
        run_clm=0
        run_band=0

        # run selected test

        log=""

        err_terrain=os.system(bindir+"/terrain "+namelist)
        if err_terrain != 0:
            print "Terrain crashed!!",err_terrain
        err_sst=os.system(bindir+"/sst "+namelist)
        if err_terrain != 0:
            print "SST crashed!!"
        err_icbc=os.system(bindir+"/icbc "+namelist)
        if err_terrain != 0:
            print "ICBC crashed!!"

        #print datadir
        #print simdir
        #print namelist

if __name__ == "__main__":
    main()

