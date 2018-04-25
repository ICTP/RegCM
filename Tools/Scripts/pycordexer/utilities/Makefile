all: hgt.so vertint.so humid.so

%.so : %.F90
	$(F2PYC) -c -m $* $<

clean:
	rm -f hgt.so vertint.so humid.so hgt.*.so vertint.*.so humid.*.so
