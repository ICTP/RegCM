#
# Copyright (C) 2010 by Graziano Giuliani
#
#  graziano.giuliani at aquila.infn.it
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. (see COPYING)
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# LIC: GPL
#

include ../Makefile.inc

all: dom2netcdf icbc2netcdf rcm2netcdf

CXX = icpc
CXX = g++
CXXFLAGS = -O2 -g
CPPFLAGS = -I. $(NETCDFINC) $(SVNDEF)
CPPFLAGS = -I. $(SVNDEF)
LIBS = $(NETCDFLIB) -lnetcdf_c++ -lnetcdf
LIBS = -lnetcdf_c++ -lnetcdf

OBJS = rcminp.o rcmio.o rcmNc.o calc.o

calc.o: calc.cpp calc.h

rcminp.o: rcminp.cpp rcminp.h

rcmio.o: rcmio.cpp rcmio.h

rcmNc.o: rcmNc.cpp rcmNc.h rcmio.h calc.h

rcm2netcdf: rcm2netcdf.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(OBJS) $(LIBS)

dom2netcdf: dom2netcdf.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(OBJS) $(LIBS)

icbc2netcdf: icbc2netcdf.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(OBJS) $(LIBS)

clean:
	rm -f rcm2netcdf dom2netcdf icbc2netcdf *.o
