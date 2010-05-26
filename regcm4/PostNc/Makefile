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

all:: dom2netcdf icbc2netcdf rcm2netcdf ncsigtop
	@cp dom2netcdf icbc2netcdf rcm2netcdf ncsigtop $(BINDIR)
	@echo Programs all installed in $(BINDIR)

LIBS = $(CXXNETCDF)

ifeq ($(AIXSP6),1)
OBJS = rcminp.o  getopt_long.o rcmio.o rcmNc.o calc.o gradsctl.o
CPPFLAGS = -I. $(NETCDFINC) $(CXXSVN) -Ilibgetopt
else
OBJS = rcminp.o rcmio.o rcmNc.o calc.o gradsctl.o
CPPFLAGS = -I. $(NETCDFINC) $(SVNDEF)
endif

getopt_long.o: libgetopt/getopt_long.c
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DREPLACE_GETOPT -c $<

gradsctl.o: gradsctl.cpp gradsctl.h

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

ncsigtop: ncsigtop.cpp calc.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< calc.o $(LIBS)

clean:
	rm -f rcm2netcdf dom2netcdf icbc2netcdf ncsigtop *.o
