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

import os,sys,shutil,fileinput
from optparse import OptionParser

class MyOptions(object):
    
    def __init__(self):
        regcm_root=os.getcwd()+"/../../"
        usage = "usage: %s [-d data directory] [-s simulation directory] [-n namelist]" % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-d","--data-dir",dest="datadir",default=regcm_root+"/data",
                      help="Directory with the data needed by preproc", metavar="DIR")
#        self.parser.add_option("-q", "--quiet",
#                      action="store_false", dest="verbose", default=True,
#                      help="Don't print status messages to stdout")
        self.parser.add_option("-s","--sim-dir", dest="simdir",
                  help="Directory where to run the simulation", metavar="DIR")
        self.parser.add_option("-n","--namelist", dest="namelist",default=regcm_root+"/regcm.in_template",
                  help="RegCM input namelist file", metavar="FILE")

    def parse(self,args):
        
        (self.options, self.args) = self.parser.parse_args(args=args)
        self.datadir = self.options.datadir
        self.simdir = self.options.simdir
        self.namelist = self.options.namelist

        if not self.simdir:
              print self.parser.format_help()
              self.parser.error("not enough arguments.")
        if not os.path.isdir(self.datadir):
            print "The directory",self.datadir,"does not exist or is not accessible!"
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

        line = line.replace("/set/this/to/where/your/surface/dataset/is",datadir)
        line = line.replace("/set/this/to/where/your/input/global/data/is",datadir)

        line = line.replace("/set/this/to/where/your/domain/file/is",simdir+"/input")
        line = line.replace("/set/this/to/where/your/icbc/for/model/is",simdir+"/input")

        line = line.replace("/set/this/to/where/your/output/files/will/be/written",simdir+"/output")

        print line.rstrip()
        
    fileinput.close()
    
            
def main():

    regcm_root=os.getcwd()+"/../../"

    options = MyOptions()
    options.parse(sys.argv[1:])

    datadir = options.datadir
    simdir = options.simdir
    namelist = options.namelist

    # create simulation directory tree
    if not os.path.isdir(simdir+"/input"):
        os.mkdir(simdir+"/input")
    if not os.path.isdir(simdir+"/output"):
        os.mkdir(simdir+"/output")

    if os.path.isfile(simdir+"/regcm.in"):
        shutil.copy(simdir+"/regcm.in",simdir+"/regcm.bak")
        shutil.copy(namelist,simdir+"/regcm.in")
    else:
        shutil.copy(namelist,simdir+"/regcm.in")
    
    namelist = simdir+"/regcm.in"

    if os.path.islink(simdir+"/Bin"):
        os.rename(simdir+"/Bin",simdir+"/Bin-old")
        os.symlink(regcm_root+"/Bin",simdir+"/Bin")
    else :
        os.symlink(regcm_root+"/Bin",simdir+"/Bin")

    #edit the namelist here

    edit_namelist(namelist,datadir,simdir)

    print datadir
    print simdir
    print namelist

if __name__ == "__main__":
    main()

