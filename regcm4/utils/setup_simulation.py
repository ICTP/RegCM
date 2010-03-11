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
   open(filename,'w').write(text_file)
   return 


def readinput(filename,verbose, parameters):
   """This function reads the input file"""
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



class MyOptions(object):
    def __init__(self):
        usage = "usage: %s [-q] [-c|-d directory] " % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-f", "--file", dest="filename",
                      help="input fle to read ", metavar="FILE")
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

def main():
    options = MyOptions()
    options.parse(sys.argv[1:])
    parameters = {}
    if options.clean:
        cleandir(options.verbose)
    elif options.filename:
        readinput(options.filename, options.verbose, parameters)
	for key, value in parameters.iteritems():
		print "%s has value %s" % (key, value)
        # write F90 modules 
        write_modulef90(parameters,'mod_regcm_param.F90',options.verbose)
        write_modulef90(parameters,'mod_param.f90',options.verbose)
    else:
        sys.stderr.write("I don't know what to do!!!\n")
            
    sys.stdout.write("Completed\n")
    sys.exit(0)


if __name__ == "__main__":
	main()

