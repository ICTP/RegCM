all: hgt.so vertint.so humid.so

%.so : %.F90
	$(F2PYC) $(F2PYOPT) -c -m $* $<

clean:
	rm -f hgt.so vertint.so humid.so hgt.*.so vertint.*.so humid.*.so
