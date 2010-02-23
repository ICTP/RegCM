#!/usr/bin/env python

import os
import sys
from optparse import OptionParser



def main():

#print os.getcwd()

	usage = "usage: %prog [ -d -q ] -s starting date  -f final date "

	parser = OptionParser(usage)
	parser.add_option("-d", "--dir", dest="dirname",
                  help="specify ICBC files directory (default Input) ", metavar="FILE")
	parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
	parser.add_option("-s","--start", dest="startdate",
		  help="starting date (1961) ")
	parser.add_option("-f","--final", dest="finaldate",
		  help=" final date ( 1963) ", )
	# only one default here: the ICBC directory:
	parser.set_defaults(dirname="../Input")

	(options, args) = parser.parse_args()

	if len(args) == 0:
            print "ciao"
	    sys.exit() 
	if options.verbose:
       	  print "set ICBC directory to  %s..." % options.dirname
       	  print "starting date is %s " % options.startdate
       	  print "final date  is %s " % options.finaldate
       


	sys.exit()

def altro():
 f = file("gal2o3.ph.in")
 data = f.read()
 f.close()
 f1=file("qsubmit-template.in")
 data1=f1.read()
 f1.close()

 for i in range(1, 12):
        # define string:
        dd="q_%02d" % i
        # create directory 
        res=os.system("mkdir %s" % dd )
        if res!=0 : 
		print "directory exist" 
        # copy data directory 
        res=os.system("cp -r data %s" % dd) 
        # Write the input file
        #define string in the format 0digit
        f=file("%s/gal2o3.in.%s" % (dd,dd) , "w")
        f.write(data % (i, i))
        f.close()
        #write the submission file:
        nome="submit-q_%02d.sh" % (i)
	f1=file("submit-q_%02d.sh" % (i), "w")
        f1.write(data1 % (i, i))
        f1.close()
        # Run the code
        for coda in ["blade"]:
            res=os.system("qsub -q %s %s" % (coda, nome))

if __name__ == "__main__":
    main()
