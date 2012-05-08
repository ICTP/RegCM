# ----------------------------------------------
SRCS =  mo_f2kcli.f90 \
	emcre_tools.f90 \
	emcre_netcdf.f90 \
	emcre.f90 

OBJS := $(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) $(INCLUDES) $(LIBS) -c $<

# ----------------------------------------------
mo_f2kcli.o : mo_f2kcli.f90
emcre_tools.o : emcre_tools.f90
emcre_netcdf.o : emcre_netcdf.f90 \
	emcre_tools.o 
emcre.o : emcre.f90 \
	emcre_tools.o \
	emcre_netcdf.o \
	mo_f2kcli.o
# ----------------------------------------------
