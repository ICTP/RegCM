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
# simple (?)  python script to generate the icbc.x 
# reading from domain.param file.
# written by S.Cozzini

import os
import sys
from optparse import OptionParser
import re
import shutil

def cleandir(verbose=False):
    """This function remove all files created by this script  """
    if verbose: sys.stdout.write("Removing files...  \n" )
    shutil.rmtree("icbc.x")

def readparam(param, verbose=False):
    """This function reads domain.param """
    #define regular expression: only interested in three lines here..
    f = file(param)
    data = f.readlines()
    f.close()
    for line in data:
    #searching the three matches I need 
           dp=re.search("parameter\(DATTYP='(\w+)'",line)
           sp=re.search("parameter\(SSTTYP='(\w+)'",line)
           ap=re.search("parameter\(AERTYP='\D+(\d+)\w+'",line)
           if dp:
	      dattyp=dp.group(1)
           if sp:
	      ssttyp=sp.group(1)
           if ap:
	      aertyp=ap.group(1)

def execute_stuff()
    """ This function execute the correct executable in the iCBC dir""" 
    # remove files if there.. 
        




class MyOptions(object):
    def __init__(self):
        usage = "usage: %s [-q] [-f regcm.param] " % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-f", "--file", dest="param",
                      help="specify regcm.param file  ", metavar="FILE")
        self.parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
        self.parser.add_option("-c", "--clean",
                      action="store_true", dest="clean", default=False,
                      help="remove all previuous generated files")
    
    def parse(self,args):
      (self.options, self.args) = self.parser.parse_args(args=args) 
      self.verbose=self.options.verbose
      self.clean = self.options.clean
      self.param = self.options.param

      #import pdb; pdb.set_trace()
      if not self.param:
          if not self.clean:
              print self.parser.format_help()
              self.parser.error("param file not supplied ")
      else:
          if self.clean:
              print self.parser.format_help()
              self.parser.error("Only one options among '-c' and '-f' has to be supplied")

def main():
    options = MyOptions()
    options.parse(sys.argv[1:])

    if options.clean:
        cleandir(options.verbose)
    elif options.param:
        readparam(options.param, options.verbose) 
        execute_stuff 
    else:
        sys.stderr.write("I don't know what to do!!!\n")
            
    sys.stdout.write("Completed\n")
    sys.exit(0)

if __name__ == "__main__":
    main()

#open(23,file='../ICBC/icbc.x')
#      write(23,'(a)') '#!/bin/csh -f'
#      write(23,'(a)') 'make clean'
#      write(23,'(a)') 'foreach FILE (RCM_SST.ctl RCM_SST.dat SST.RCM)'
#      write(23,'(a)') 'if ( -f $FILE ) /bin/rm $FILE'
#      write(23,'(a)') 'end'
#      IF(DATTYP.EQ.'FVGCM'.or.
#     &        (DATTYP.EQ.'FNEST'.and.
#     &         (SSTTYP.EQ.'FV_RF'.or.SSTTYP.EQ.'FV_A2'))) THEN
#         write(23,'(a)') 'make SST_FVGCM'
#         write(23,'(a)') './SST_FVGCM'
#         write(23,'(a)') '/bin/rm -f SST_FVGCM*.o SST_FVGCM'
#      ELSE IF( DATTYP.EQ.'EH5OM'.or.
#     &        (DATTYP.EQ.'FNEST'.and.
#     &         (SSTTYP.EQ.'EH5RF'.or.SSTTYP.EQ.'EH5A2'.or.
#     &          SSTTYP.EQ.'EH5B1'.or.SSTTYP.EQ.'EHA1B'))) THEN
#         write(23,'(a)') 'make SST_EH5OM'
#         write(23,'(a)') './SST_EH5OM'
#         write(23,'(a)') '/bin/rm -f SST_EH5OM*.o SST_EH5OM'
#      ELSE IF((DATTYP.EQ.'EIN15'.or.DATTYP.EQ.'ERAIN').and.
#     &        (SSTTYP.EQ.'ERSST'.or.SSTTYP.EQ.'ERSKT')) THEN
#         write(23,'(a)') 'make SST_ERSST'
#         write(23,'(a)') './SST_ERSST'
#         write(23,'(a)') '/bin/rm -f SST_ERSST*.o SST_ERSST'
#      ELSE
#         write(23,'(a)') 'make SST'
#         write(23,'(a)') './SST'
#         write(23,'(a)') '/bin/rm -f SST_1DEG*.o SST'
#      ENDIF
#      IF(.not.(AERTYP(4:5).eq.'00')) THEN
#         write(23,'(a)') 'make AEROSOL'
#         write(23,'(a)') './AEROSOL'
#         write(23,'(a)') '/bin/rm -f AEROSOL.o AEROSOL'
#      ENDIF
#      IF(DATTYP.EQ.'ERAHI') THEN
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
#      write(23,'(a)') 'make ICBC'
#      write(23,'(a)') './ICBC'
#      write(23,'(a)') '/bin/rm -f ICBC.o ICBC SST.RCM'
##
