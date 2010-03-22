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
# python script to perform daily and monthly average in one single shot
# written by S.Cozzini
#

from average import average_day
from average import month_average
#import shutils 
import os 
import re
#

def get_size(filename):
  ix=0
  iy=0
  try:
    f = open(filename,'r')
    data = f.readlines()
    f.close()
    for line in data:
       g=re.search("pdef\s+(\d+)\s+(\d+)\s*lccr",line)
       if g:
       # keys are variable in the code so let us consider them all lowercase 
         iy=g.group(1)
         ix=g.group(2)
  except IOError, error:
    print  filename,' is not available (%s)' % error
  return (ix,iy)

# days for months 

days_in_months= {'01':31,'02':28,'03':31,'04':30,'05':31,'06':30,'07':31,'08':31,'09':30,'10':31,'11':30,'12':31 }
listafile=[]
# check SRF files in the dir.. 
dirname=os.getcwd()
listafile=dict()
listadays=dict()
for namefile in os.listdir(dirname):
    if not os.path.isfile(namefile):
       print "%s not a file: skipping  " % namefile 
    else:   
      g=re.search("SRF\.(\d\d\d\d)(\d\d)(\d\d)\d\d$",namefile) 
      if g:
        year=g.group(1)
	month=g.group(2)
	day=g.group(3)
        print "  processing SRF file for year: %s" %  namefile 
        (ix,iy)=get_size(namefile+".ctl") 
        day_to_process= int(days_in_months[month]) +1 - int(day)
        print " SRF file for year: %s , month: %s, day for months: %s; size is %d,%d" % (year ,month, day_to_process,int(ix),int(iy)) 
        fileoutput= namefile+"mon" 
        listafile.setdefault(year, []).append(fileoutput)  	
        listadays.setdefault(year, []).append(day_to_process)  	
#       
#        average_day(int(ix),int(iy),namefile,fileoutput,day_to_process)	
#
for year,value  in listafile.iteritems():
#  print year,value,listadays[year]
  print " processing year  %s" % year 
  outputfile= "out" + year
  files= "".join(value) 
  month_average(int(iy),int(ix),files,outputfile,listadays[year])
  print " completed: monthly average written in %s " % outputfile 
   
#  
