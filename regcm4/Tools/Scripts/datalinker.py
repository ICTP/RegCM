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
#  simple (?)  python script to link DATA directories
# it performs  some checks on existence of DATA directories
# improvements: check against a well defined list of directories

# written by S.Cozzini, improved by A.Messina 
##
   
import os
import sys
import shutil
from optparse import OptionParser

def cleandir(verbose=False):
    """This function remove all directories in the current directory"""
    for filename in os.listdir(os.getcwd()):
        if os.path.isdir(filename) and filename != ".svn":
            if verbose: sys.stdout.write("Removing directory %s\n" % filename)
            shutil.rmtree(filename)

def createlinks(dirname, verbose):
    """Create needed symlinks from directory 'dirname'"""
    wrong=False
    if verbose:
       print "start processing %s " % dirname
    #check if DATA exists and what it contains.. 
    if not os.path.isdir(dirname):
        sys.stderr.write("'%s is not a directory\n" % dirname)
        raise SystemExit, 1
    #check content.. 
    for filename in os.listdir(dirname):
        absolute_name=os.path.join(dirname,filename)
        sys.stdout.write("filename to process is %s \n" % absolute_name)
        if not os.path.isdir(absolute_name):
            sys.stdout.write("%s is not a  directory: will skip link\n" % filename )
            continue
        sys.stdout.write("%s is a directory: try to create link\n" % filename)
        if os.path.isdir(os.path.join(os.getcwd(),filename)):
           sys.stdout.write("%s already present, skipping..\n " % filename)
           wrong=True
           continue
        os.mkdir(filename)
        cmd='ln -sf %s/* %s/.' % (absolute_name, filename)
        print 'executing', cmd
        failure=os.system(cmd)
        if failure:
           sys.stdout.write("something wrong for %s  \n" % filename )
           wrong=True
    # post-opt check
    if wrong:
        print ' script completed but something went wrong: please check !'
    else:
        print 'successfully created links to DATA directories'

class MyOptions(object):
    def __init__(self):
        usage = "usage: %s [-q] [-c|-d directory] " % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-d", "--dir", dest="dirname",
                      help="specify DATA directory ", metavar="FILE")
        self.parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
        self.parser.add_option("-c", "--clean",
                      action="store_true", dest="clean", default=False,
                      help="remove all directories in DATA")
    
    def parse(self,args):
      (self.options, self.args) = self.parser.parse_args(args=args) 
      self.verbose=self.options.verbose
      self.clean = self.options.clean
      self.dirname = self.options.dirname

      #import pdb; pdb.set_trace()
      if not self.dirname:
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

    if options.clean:
        cleandir(options.verbose)
    elif options.dirname:
        createlinks(options.dirname, options.verbose) 
    else:
        sys.stderr.write("I don't know what to do!!!\n")
            
    sys.stdout.write("Completed\n")
    sys.exit(0)


if __name__ == "__main__":
	main()

