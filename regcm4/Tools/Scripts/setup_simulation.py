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

# Helper script for setting up a RegCM simulation by M. Scarcia
#

import os,sys,shutil,fileinput
from optparse import OptionParser

class MyOptions(object):
    
    def __init__(self):
        usage = "usage: %s [-d data directory] [-s simulation directory] [-b binaries directory] [-n namelist]" % sys.argv[0]
        self.parser = OptionParser(usage)
        
        self.parser.add_option("-d","--data-dir",dest="datadir",default="./data",
                  help="Directory with the data needed by preproc", metavar="DIR")
        
        self.parser.add_option("-s","--sim-dir", dest="simdir",default ="./",
                  help="Directory where to run the simulation", metavar="DIR")
        
        self.parser.add_option("-b","--bin-dir", dest="bindir",default="../../Bin",
                  help="Directory where RegCM binaries are", metavar="DIR")
        
        self.parser.add_option("-n","--namelist", dest="namelist",default="./regcm.in",
                  help="RegCM input namelist file", metavar="FILE")

    def parse(self,args):
        
        (self.options, self.args) = self.parser.parse_args(args=args)
        self.datadir = self.options.datadir
        self.simdir = self.options.simdir
        self.bindir = self.options.bindir
        self.namelist = self.options.namelist

        if not self.simdir:
              print self.parser.format_help()
              self.parser.error("not enough arguments.")
        if not os.path.isdir(self.datadir):
            print "The directory",self.datadir,"does not exist or is not accessible!"
            sys.exit(1)

        if not os.path.isdir(self.bindir):
            print "The directory",self.bindir,"does not exist or is not accessible!"
            sys.exit(1)

        if not os.path.isdir(self.simdir):
            try :
                os.mkdir(self.simdir)
            except OSError:
                print "Cannot create",self.simdir,"directory!"
                sys.exit(1)

        if not os.path.isfile(self.namelist):
            print "The namelist provided (",self.namelist.strip(),") does not exist or is not accessible!"
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

    options = MyOptions()
    options.parse(sys.argv[1:])

    datadir = options.datadir
    simdir = options.simdir
    bindir = options.bindir
    namelist = options.namelist

    datadir = os.path.abspath(datadir)
    simdir = os.path.abspath(simdir)
    bindir = os.path.abspath(bindir)

    # create simulation directory tree
    if not os.path.isdir(simdir+"/input"):
        os.mkdir(simdir+"/input")
    if not os.path.isdir(simdir+"/output"):
        os.mkdir(simdir+"/output")

    # copy namelist file
    if os.path.isfile(simdir+"/regcm.in"):
        shutil.copy(simdir+"/regcm.in",simdir+"/regcm.bak")
        shutil.copy(namelist,simdir+"/regcm.in")
    else:
        shutil.copy(namelist,simdir+"/regcm.in")
    
    namelist = simdir+"/regcm.in"

    # link binaries directory
    if os.path.islink(simdir+"/Bin"):
        os.rename(simdir+"/Bin",simdir+"/Bin-old")
        os.symlink(bindir,simdir+"/Bin")
    else :
        os.symlink(bindir,simdir+"/Bin")

    #edit the namelist
    edit_namelist(namelist,datadir,simdir)

    # print some info
    print "Preprocessing data will be read from",datadir
    print "Simulation will be run in",simdir,"using binaries found in",bindir
    print "Please note that no checks are performed on the actual presence of data and/or binaries."

if __name__ == "__main__":
    main()

