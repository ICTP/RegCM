
all: hgt.so vertint.so humid.so

%.so : %.F90
	f2py2 -c -m $* $<

hgt.so: hgt.F90
vertint.so: vertint.F90
humid.so: humid.F90

clean:
	rm -f hgt.so vertint.so humid.so
