# Use the c++4 library at
#    http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-cxx4-4.2.tar.gz

CXX = g++
CXXFLAGS = -O2 -mtune=native -g
CPPFLAGS = `ncxx4-config --cflags`
LD_FLAGS = `ncxx4-config --libs`

srfread: srfread.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LD_FLAGS)

clean:
	rm -f srfread
