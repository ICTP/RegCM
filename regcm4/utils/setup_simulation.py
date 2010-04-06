#!/usr/bin/env python
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of RegCM model.
#
#    RegCM model is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RegCM model is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
#
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# python script to read domani.param file (regcm v3) and generate all input files
# needed by all regcm4 steps:
#  mod_regcm_param
#  mod_param 
#  
# written by S.Cozzini
#
   
import os
import sys
import shutil
from optparse import OptionParser
import re
import fileinput
import string



def write_modulef90(parameters,filename,verbose):
   """This function writes filename with new values"""
   #backup original file
   shutil.copyfile(filename,filename+".bak")	
   text_file=open(filename).read()
   for key, value in parameters.iteritems():
       pattern=key + "\s*=\s*\S+"
       valuereplace= key + " = " + value
       text_file=re.sub(pattern,valuereplace,text_file)
   if verbose:
      print text_file 
   f= open(filename,'w')
   f.write(text_file)
   f.close() 
   return 


def readinput(filename,verbose):
   """This function reads the input file"""
   parameters={}
   try:
     f = open(filename,'r')
     data = f.readlines()
     f.close()
     for line in data:
       g=re.search("parameter\((\w+)\s*=\s*(\S+)\)",line)
       if g:
       # keys are variable in the code so let us consider them all lowercase 
         parameters[string.lower(g.group(1))] = g.group(2)
   except IOError, error:
     print  filename,' is not available (%s)' % error
   return parameters

def cleandir(verbose=False):
    """This function remove all directories in the current directory"""
    for filename in os.listdir(os.getcwd()):
        if os.path.isdir(filename):
            if verbose: sys.stdout.write("Removing directory %s\n" % filename)
            shutil.rmtree(filename)


def write_ICBC(parameters,filename,verbose=False):
    """This functions writes ICBC script to generate ICBC input files"""
    
    
    fhandler = open(filename,"w")
    fhandler.write('#!/bin/csh \n')
    fhandler.write('#cleaning files.. \n')
    fhandler.write("foreach FILE (RCM_SST.ctl RCM_SST.dat SST.RCM)\n")
    fhandler.write("if ( -f $FILE ) /bin/rm $FILE\n")
    fhandler.write("end\n")

    # identify the variable we are interested in..  
    DATTYP=parameters['dattyp'].strip("'")
    SSTTYP=parameters['ssttyp'].strip("'")
    AERTYP=parameters['aertyp'].strip("'")
 
    default="./sst_1deg \n" 
    # various conditions to check 
    if ((DATTYP=='FVGCM') or  (DATTYP == 'FNEST' and (SSTTYP== 'FV_RF' or SSTTYP == 'FV_A2') ) ):
        default="./sst_fvgcm\n"
    if ((DATTYP=='EH5OM') or (DATTYP=='FNEST' and (SSTTYP=='EH5RF' or SSTTYP=='EH5A2' or SSTTYP=='EH5B1' or SSTTYP=='EHA1B'))):
        default="./sst_eh50m\n"
    if ((DATTYP=='EIN15' or DATTYP=='ERAIN') and (SSTTYP=='ERSST' or SSTTYP=='ERSKT')):
        default="./sst_ersst\n" 
 
    fhandler.write(default) 

    # chemistry enabled ? 
    if (AERTYP[3:5]!='00'):
          fhandler.write('./aereosol \n') 
      
    if (DATTYP=='ERAHI'): 
          fhandler.writelines('echo: not yet ported ERAHI \nexit \n') 

#         write(23,'(a)') 'cp ../../Commons/tools/srcERAHI/*.f .'
#         write(23,'(a)') 'make ERAHI_HT'
#         write(23,'(a)') './ERAHI_HT'
#         write(23,'(a)') 'make ERAHI_PS'
#         write(23,'(a)') './ERAHI_PS'
#         write(23,'(a)') 'make ERAHI_T'
#         write(23,'(a)') './ERAHI_T'
#         write(23,'(a)') 'make ERAHI_Q'
#         write(23,'(a)') './ERAHI_Q'
#         write(23,'(a)') 'make ERAHI_U'
#         write(23,'(a)') './ERAHI_U'
#         write(23,'(a)') 'make ERAHI_V'
#         write(23,'(a)') './ERAHI_V'
#      ENDIF
    
    fhandler.writelines('./icbc \n') 
    fhandler.close() 


class MyOptions(object):
    def __init__(self):
        usage = "usage: %s [-q] [-c|-d directory] " % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-f", "--file", dest="filename",
                      help="input file to read ", metavar="FILE")
        self.parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
        self.parser.add_option("-c", "--clean",
                      action="store_true", dest="clean", default=False,
                      help="remove all input files")
    
    def parse(self,args):
      (self.options, self.args) = self.parser.parse_args(args=args) 
      self.verbose=self.options.verbose
      self.clean = self.options.clean
      self.filename = self.options.filename

      #import pdb; pdb.set_trace()
      if not self.filename:
          if not self.clean:
              print self.parser.format_help()
              self.parser.error("not enough arguments.")
      else:
          if self.clean:
              print self.parser.format_help()
              self.parser.error("Only one options among '-c' and '-d' has to be supplied")

def check():
   filename="../Makefile.inc"
   try:
      f=open(filename,"r")
      makeinc=f.readlines()
      f.close()

      datadir=''
      basedir=''
      ncdfinc=''
      ncdflib=''
      
      for line in makeinc :
         if line.find('REGCM_DATA_DIR') > -1 : 
            tmp = line.rsplit('=')
            datadir = tmp[1].rstrip()
            datadir = datadir.replace(' ','')
            if os.path.isdir(datadir) :
               print 'The data directory',datadir,'exists.'
            else :
               print 'WARNING! The data directory',datadir,'does not exist!'
         if line.find('REGCM_BASE_DIR') > -1 : 
            tmp = line.rsplit('=')
            basedir = tmp[1].rstrip()
            basedir = basedir.replace(' ','')
            if os.path.isdir(basedir) :
               print 'The base directory',basedir,'exists.'
            else :
               print 'WARNING! The data directory',basedir,'does not exist!'
               #sys.exit(1)
         if line.find('NETCDFINC =') > -1 : 
            tmp = line.rsplit('=')
            ncdfinc = tmp[1].rstrip()
            if ncdfinc ==  '':
               print 'WARNING! The NetCDF include directory',ncdfinc,'was not defined!'
               #sys.exit(1)
            else :
               print ncdfinc,'defined as NetCDF include directory. Please be sure is the correct one.'
         if line.find('NETCDFLIB =') > -1 : 
            tmp = line.rsplit('=')
            ncdflib = tmp[1].rstrip()
            if ncdflib ==  '':
               print 'WARNING! The NetCDF library',ncdflib,'was not defined!'
               #sys.exit(1)
            else :
               print ncdflib,'defined as NetCDF library. Please be sure is the correct one.'  
   except IOError, error:
      print  filename,' is not available!'
      sys.exit(1)
	
def main():
    options = MyOptions()
    options.parse(sys.argv[1:])
    parameters = {}

    check()
# check if Makefile inc exists and all the variables are defined:
# variable to be there:
#   netcdf       error 
#   datadir ( ?) warning 
#   regcmbasedir error 
#   Martin  

    if options.clean:
        cleandir(options.verbose)
    elif options.filename:
        # read parameters from doman.parama files.. 
        parameters = readinput(options.filename, options.verbose)
	for key, value in parameters.iteritems():
		print "%s has value %s" % (key, value)
        # write F90 modules 
        # define three standard places: this should be improved at a later stage 
        default_mod_regcm_param="../Config/mod_regcm_param.F90"
        default_mod_paramT="../Config/mod_preproc_param.f90" 
        write_modulef90(parameters,default_mod_regcm_param,options.verbose)
        write_modulef90(parameters,default_mod_paramT,options.verbose)
        default_ICBC="../PreProc/ICBC/icbc_v4.x"
        write_ICBC(parameters,default_ICBC,options.verbose)
    else:
        sys.stderr.write("I don't know what to do!!!\n")
            
    sys.stdout.write("Completed\n")
    sys.exit(0)


if __name__ == "__main__":
	main()

