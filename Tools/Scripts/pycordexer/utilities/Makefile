all: hgt.so vertint.so capecin.so

%.so : %.F90
	$(F2PYC) $(F2PYOPT) -c -m $* $<

clean:
	rm -f hgt.so vertint.so hgt.*.so vertint.*.so capecin.so capecin.*.so
