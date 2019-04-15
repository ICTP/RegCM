# ----------------------------------------------
export
SHELL     = sh

F90       = ifort
F90FLAGS  = -cpp -O3 
#F90FLAGS = -cpp -O0 -g \
             -check arg_temp_created,bounds,format,output_conversion,uninit \
             -traceback -warn all -fp-model precise -fp-model source -ftz \
             -fimf-arch-consistency=true
LIBS      = -L/opt/netcdf/4.1.2/intel/2011/lib -lnetcdf -lnetcdff
INCLUDES  = -I/opt/netcdf/4.1.2/intel/2011/include

# ----------------------------------------------

rootdir     = $(PWD)
srcdir     = ./src
bindir     = ./bin
includedir = ./include
PROG       = emcre

# ----------------------------------------------

all:
	cd $(srcdir) ; \
	$(MAKE)

install:
	if test ! -d $(bindir) ; then mkdir $(bindir) ; fi
	cp -f $(srcdir)/$(PROG) $(bindir)/.
	if test ! -d $(includedir) ; then mkdir $(includedir) ; fi
	cp -f $(srcdir)/*.mod $(includedir)/.

clean:
	rm -f $(srcdir)/*.o $(srcdir)/*~ $(srcdir)/*.mod $(srcdir)/$(PROG)

distclean: clean
	rm -fr $(bindir) $(includedir) 

zip: clean 
	zip -r emcre.zip .
        

# ----------------------------------------------
