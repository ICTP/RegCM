#!/usr/bin/env python


import os
import sys
from optparse import OptionParser
import re


class MyOptions(object):
    def __init__(self):
        usage = "usage: %s [-q] [-c|-d directory] " % sys.argv[0]
        self.parser = OptionParser(usage)
        self.parser.add_option("-f", "--file", dest="filename",default="regcm.in",
                      help="file to read ( default regcm.in in current directory", metavar="FILE")
        self.parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")
	self.parser.add_option("-d", "--dir", dest="dirname",default="../Input",
                  help="specify ICBC files directory (default ../Input directoty) ", metavar="FILE")
	self.parser.add_option("-s","--start", dest="startdate",
		  help="starting date (i.e. 1961) ")
	self.parser.add_option("-e","--end", dest="finaldate",
		  help=" final date (i. e. 1963) ", )
        self.parser.add_option("-u","--unit", dest="fortranunit",default="101",
		  help=" fortran unit to start (default: 101) ", )
       
   
    def parse(self,args):
      (self.options, self.args) = self.parser.parse_args(args=args) 
      self.verbose=self.options.verbose
      self.filename = self.options.filename
      self.dirname = self.options.dirname
      self.finaldate=self.options.finaldate
      self.startdate=self.options.startdate
      self.fortranunit=self.options.fortranunit
      #import pdb; pdb.set_trace()
      # check if dates are provided..
      if self.startdate:
	if not self.finaldate:
           print "final date not specified"
           self.parser.error("final date not specified")
      if self.options.finaldate:
	if not self.options.startdate:
           print "final date not specified"
           self.parser.error("starting date not specified")
      # error checks on inputs:
      if not os.path.isdir(self.dirname): 
          self.parser.error(" is not a directory.. ")

def parse_date_from_file(filename,verbose):
    """Parse from file starting and ending date"""
    if verbose: sys.stdout.write("Reading Input file: %s\n" % filename)
    f = open(filename,'r')
    data = f.readlines()
    f.close()
    for line in data:
      start=re.search("\idate1\s*=\s*(\S+),",line)
      end=re.search("\idate2\s*=\s*(\S+),",line)
      if start:
        startdate=start.group(1)
      if end:
	finaldate=end.group(1)
    return startdate,finaldate

def check_dates(startdate,finaldate):
    "simple function to check if final date is greater than starting date" 
    if (startdate > finaldate):
      print " error: initial date is greater than final date: please check" 
      sys.exit(1)  
    

def check_and_create_link(suffix, directory, unit):
    """ checks if ICBC files are available in the directory specified and create symlink"""
    string="%s" % (unit)
    namefile=directory+ "ICBC" + suffix 
    if not os.path.isfile(namefile):
       print "%s does not exist" % namefile 
       sys.exit(1)
    os.symlink(namefile,"fort."+string)
	

def main():
    options = MyOptions()
    options.parse(sys.argv[1:])
  
    if options.verbose:
       print "set ICBC directory to  %s..." % options.dirname
       print "starting date is %s " % options.startdate
       print "final date  is %s " % options.finaldate
    startdate=options.startdate
    finaldate=options.finaldate

    # if dates are not given as input let us read from file..    
    if not options.finaldate:
       (startdate,finaldate)=parse_date_from_file(options.filename,options.verbose)    
    if (startdate > finaldate):
      print " error: initial date is greater than final date: please check" 
      sys.exit(1)
    

    # identifying month and years.. 
    startingyear=startdate[0:4]
    finalyear=finaldate[0:4]
    startingmonth="01"
    finalmonth="12"
    if len(startdate) > 4: 
      startingmonth=startdate[4:6]      
    if len(finaldate) >4:
      finalmonth=finaldate[4:6]

    startingstring=startingyear + startingmonth 
    finalstring=finalyear +finalmonth
    
    # start looping over month and year:
    if options.verbose: sys.stdout.write("starting loop from: Year: %s Month: %s to Year: %s Month: %s \n" % (startingyear, startingmonth,finalyear, finalmonth))
    string=startingstring
    month=int(startingmonth)
    year=int(startingyear)
    unit=int(options.fortranunit) 
    while (string <= finalstring) :
        suffix=string + "0100"
        check_and_create_link(suffix,options.dirname,unit)
        month= month + 1
        unit=unit + 1 
        if (month == 12):
           year=year+1
           month=1 
	string="%s%02d" % (year, month)
        
        #define string in the format 0digit
    print "done for range",finalyear, finalmonth ,startingyear,startingmonth 


if __name__ == "__main__":
    main()
