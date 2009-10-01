# makefile for foo2_to_v5d.f conversion program

# The input to the conversion program is your 5-D grid format.
# The ouput is a v5d file.

PROGRAM = postpv5d

# If DEC or Linux, (Little-endian), add -DLITTLE to CFLAGS
# If AIX, remove the -DUNDERSCORE
CFLAGS = -c -g -DUNDERSCORE 
FFLAGS = -byteswapio -fast -Ktrap=fp -g 
CC = cc
F77 =pgf90
LIBS = -lm

OBJECTS  =  $(PROGRAM).o 


$(PROGRAM): $(OBJECTS)
	$(F77) $(OBJECTS) $(LIBS) -o $@ libv5dnew.a 


$(PROGRAM).o: $(PROGRAM).f
	$(F77) $(FFLAGS) $(LIBS)  $(PROGRAM).f

clean:
	rm -f *.o $(PROGRAM) 
