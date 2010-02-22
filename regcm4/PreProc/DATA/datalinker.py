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

# written by S.Cozzini   
import os
import os.path
import sys
from optparse import OptionParser


def main():
#print os.getcwd()
        usage = "usage: %s [-q] -d directory " % sys.argv[0]
        parser = OptionParser(usage)
        parser.add_option("-d", "--dir", dest="dirname",
                  help="specify DATA directory ", metavar="FILE")
        parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

        (options, args) = parser.parse_args()

        if options.dirname == None:
            print parser.format_help()
            parser.error("not enough arguments, a directory must be specified")
        if options.verbose:
          print "start processing %s " % options.dirname
        wrong=False
        data_location=options.dirname
        #check if DATA exists and what it contains.. 
	if os.path.isdir(data_location):
        	#check content.. 
           	for filename in os.listdir(data_location):
                    absolute_name=os.path.join(data_location,filename)
                    sys.stdout.write("filename to process is %s \n" % absolute_name)
                    if os.path.isdir(absolute_name):
     	     	        sys.stdout.write("%s is a directory: try to create link\n" % filename)
			if os.path.isdir(os.path.join(os.getcwd(),filename)):
			   sys.stdout.write("%s already present, skipping..\n " % filename)
                           wrong=True
			   continue
                        cmd='ln -sf %s/* %s/.' % (absolute_name, filename)
                        print 'executing', cmd
                        failure=os.system(cmd)
			if failure:
			   sys.stdout.write("something wrong for %s  \n" % filename )
		           wrong=True
                    else:
		        sys.stdout.write("%s is not a  directory: will skip link\n" % filename ) 	
	else:
		sys.stderr.write("'%s is not a directory\n" % data_location)
		raise SystemExit, 1 
        if wrong:
	    print ' script completed but something went wrong: please check !' 
	else: 
            print 'successfully created links to DATA directories' 

if __name__ == "__main__":
	main()

