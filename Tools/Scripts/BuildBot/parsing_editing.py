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

# Module that contains various parsing/editing functions used in other scripts

import fileinput,calendar

# replaces paths in RegCM namelist file
# requirements: fileinput

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

# parses namelist for all the simulation dates (tbi)
# requirements: fileinput, calendar

def parse_dates(namelist,simdays):

    infile = open(namelist,"r")
    file_content = infile.readlines()
    infile.close()
        
    for line in file_content:
            #if line.find("globidate1") > -1 :
            #    linea=line.rsplit("=")
            #    globidate1=linea[1]
            #    globidate1=filter(lambda x:x.isdigit(),globidate1)
            if line.find(" idate0 ") > -1 :
                linea=line.rsplit("=")
                idate0=linea[1]
                idate0=filter(lambda x:x.isdigit(),idate0)
            #if line.find(" idate1 ") > -1 :
            #   linea=line.rsplit("=")
            #   idate1=linea[1]
            #   idate1=filter(lambda x:x.isdigit(),idate1)
            if line.find(" idate2 ") > -1 :
                linea=line.rsplit("=")
                idate2=linea[1]
                idate2=filter(lambda x:x.isdigit(),idate2)
                
    year = int(idate0[:4])
    month = int(idate0[4:6])
    day_start = int(idate0[6:8])
    day_end = int(idate2[6:8])

    maxdate = calendar.monthrange(year,month)[1]

    if (day_start + simdays) <= maxdate :
    	day_end = day_start + simdays
    else :
	day_end = maxdate

    #print idate2
    #print str(year)+str(month).zfill(2)+str(day_end).zfill(2)+"00"
        
    for line in fileinput.FileInput(namelist,inplace=1):
        line = line.replace(idate2,str(year)+str(month).zfill(2)+str(day_end).zfill(2)+"00")
        print line.rstrip()

    fileinput.close()

    return idate0 # if others needed will put a list as output

# parses config file and returns options; expects a list like
#
#          OPTION1 = value
#          OPTION2 = value
#          #comment
#
# requirements: sys
def parse_config(filename):

    COMMENT_CHAR = '#'
    OPTION_CHAR =  '='
    
    options = {}
    try:
        f = open(filename)
    except :
        print "File "+filename+" does not exist or is not accessible!"
        sys.exit(1)
        
    for line in f:
        # First, remove comments:
        if COMMENT_CHAR in line:
            # split on comment char, keep only the part before
            line, comment = line.split(COMMENT_CHAR, 1)
        # Second, find lines with an option=value:
        if OPTION_CHAR in line:
            # split on option char:
            option, value = line.split(OPTION_CHAR, 1)
            # strip spaces:
            option = option.strip()
            value = value.strip()
            # store in dictionary:
            options[option] = value
    f.close()
    return options

