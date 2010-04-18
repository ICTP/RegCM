
all: rcm2netcdf

CXX = g++
CXXFLAGS = -O2 -g
CPPFLAGS = -I. -I/usr/include
LIBS = -L/usr/lib64 -lnetcdf_c++ -lnetcdf

OBJS = rcminp.o rcmio.o rcmNc.o

rcminp.o: rcminp.cpp rcminp.h

rcmio.o: rcmio.cpp rcmio.h

rcmNc.o: rcmNc.cpp rcmNc.h rcmio.h

rcm2netcdf: rcm2netcdf.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(OBJS) $(LIBS)

clean:
	rm -f rcm2netcdf *.o
