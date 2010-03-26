
CC = gcc
CFLAGS = -O2 -g

ENDIANESS = LITTLE
USCORING = UNDERSCORE
# ENDIANESS = BIG
# USCORING = NOUNDE$RSCORE

all: libv5d.a

libv5d.a: v5d.o binio.o
	$(AR) $(ARFLAGS) $@ v5d.o binio.o

v5d.o: binio.h v5d.h vis5d.h
	$(CC) $(CFLAGS) -D$(ENDIANESS) -D$(USCORING) -o $@ $<

binio.o: binio.h binio.c
	$(CC) $(CFLAGS) -D$(ENDIANESS) -D$(USCORING) -o $@ $<

clean:
	rm -f libv5d.a *.o
