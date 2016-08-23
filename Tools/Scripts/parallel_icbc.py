#!/usr/bin/env python
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#  This file is part of RegCM model.
#
#  RegCM model is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RegCM model is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#
# Simple python script to parallelize ICBC jobs.
# Requires two external libraries in python, f90nml and mpi4py
#

# Basic python libraries. Should be available on any decent install
import sys
import os
import copy
import tempfile
import calendar
import subprocess

# Mission specific libraries.
try:
    import f90nml
except:
    print("Please install f90nml: https://pypi.python.org/pypi/f90nml")
    sys.exit(1)
try:
    from mpi4py import MPI
except:
    print("Please install mpi4py: https://pypi.python.org/pypi/mpi4py/2.0.0")
    sys.exit(1)

# Classes and functions

class mydatetime(object):
    def __init__(self,ival):
        self.year = ival/1000000
        self.month = (ival-self.year*1000000)/10000
        self.day = (ival-self.year*1000000-self.month*10000)/100
        self.hour = ival-self.year*1000000-self.month*10000-self.day*100
    def __str__(self):
        return "%04d%02d%02d%02d" % (self.year, self.month, self.day, self.hour)
    def __repr__(self):
        return "%04d%02d%02d%02d" % (self.year, self.month, self.day, self.hour)
    def __int__(self):
        return self.year*1000000+self.month*10000+self.day*100+self.hour

def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month

def usage(message):
    if ( rank == 0 ):
        print("Python helper to run an icbc job in parallel.")
        print(message)
        print("USAGE :")
        print("   "+sys.argv[0]+" icbc_executable namelist_file")
        print("")
    sys.exit(1)

def monthadd(date, nmon):
    xdate = copy.deepcopy(date)
    xdate.year = xdate.year + nmon/12
    xdate.month = xdate.month + (nmon % 12)
    if xdate.month > 12:
        xdate.month = xdate.month - 12
        xdate.year = xdate.year + 1
    return xdate

def eomonth(date, cal):
    dayomon = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 )
    xdate = copy.deepcopy(date)
    if cal == "gregorian":
        if calendar.isleap(date.year):
            if date.month == 2:
                xdate.day = 29
            else:
                xdate.day = dayomon[date.month-1]
        else:
            xdate.day = dayomon[date.month-1]
    elif cal == "noleap":
        xdate.day = dayomon[date.month-1]
    elif cal == "360_days":
        xdate.day = 30
    else:
        raise RuntimeError("Unknown calendar")
    xdate.hour = 18
    return xdate

def pprint(message):
    if ( rank == 0 ):
        print(message)

#
# Execution start here
#

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

try:
    icbcexec = sys.argv[1]
    namelist = sys.argv[2]
except:
    usage("Missing Argument")

if not os.path.isfile(icbcexec):
    usage("Not a regular file : "+icbcexec)
if not os.path.isfile(namelist):
    usage("Not a regular file : "+namelist)

pprint("Parsing namelist file...")

nml = f90nml.read(namelist)
try:
    start_icbc_date = nml["globdatparam"]["gdate1"]
    stop_icbc_date = nml["globdatparam"]["gdate2"]
except:
    usage("gdate1,gdate2 not present in namelist file "+namelist)

try:
    xcal = nml["globdatparam"]["calendar"]
except:
    xcal = "gregorian"

gdate1 = mydatetime(start_icbc_date)
gdate2 = mydatetime(stop_icbc_date)
ntasks = diff_month(gdate2, gdate1)

if ntasks < size:
    pprint("Too many processors on job. Maximum required is "+repr(ntasks))
    sys.exit(1)

mytasks = ntasks/size
pp0tasks = ntasks - mytasks*size + mytasks

if rank == 0:
    mystart = 0
    myend = pp0tasks
else:
    mystart = pp0tasks + (rank-1)*mytasks
    myend = mystart + mytasks

basen = os.path.splitext(namelist)[0]
ofile = open(basen+".icbc.out"+("%04d" % (rank)),'w')
efile = open(basen+".icbc.err"+("%04d" % (rank)),'w')
for imon in xrange(mystart,myend):
    datest = monthadd(gdate1, imon)
    d1 = datest
    if imon == ntasks-1:
        d2 = gdate2
    else:
        d2 = eomonth(datest, xcal)
    print("Processor "+repr(rank)+" working on "+repr(d1)+" - "+repr(d2))
    nml["globdatparam"]["gdate1"] = int(d1)
    nml["globdatparam"]["gdate2"] = int(d2)
    nmlfile = tempfile.NamedTemporaryFile(delete=False)
    nml.write(nmlfile)
    nmlfile.close()
    subprocess.call([icbcexec, nmlfile.name], stdout=ofile, stderr=efile)
    os.unlink(nmlfile.name)
ofile.close()
efile.close()
