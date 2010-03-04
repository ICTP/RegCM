# To build this by itself, use the make target esmf_time_f90_only 
# from the top-level WRF Makefile.  
# > cd ../..
# configure
# make esmf_time_f90_only

cpp_dirs := . ..
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line
# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH    := $(subst $(space),:,$(VPATH))


ifeq ($(WRFTEST),$(null))
CPPDEF := -DNO_LEAP_CALENDAR
endif
ifneq ($(ESMF_BLD),$(null))
CPPDEF += -DREAL_ESMF -DWRFU_=ESMF_
endif

include ../../../csm_share/unit_testers/make.Macros

#.SUFFIXES: .F90 .o .f

# a hack, should pass this in from top-level Makefile
SED_FTN := sed -f standard.sed

AR  := ar
CPP := cpp

OBJS := ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        Meat.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
        ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF_Mod.o \
	ESMF_AlarmClockMod.o

debug: $(OBJS)
      echo "VPATH : $VPATH"
      echo "OBJS  : $OBJ"
      echo "FFLAGS: $FFLAGS"

default: libesmf_time.a

tests: Test1.exe

libesmf_time.a : $(OBJS)
	\rm -f libesmf_time.a
	$(AR) $(ARFLAGS) libesmf_time.a $(OBJS)
	$(RANLIB) libesmf_time.a

ifeq ($(ESMF_BLD),$(null))
Test1.exe : libesmf_time.a Test1.o
	$(FC) $(LDFLAGS) -o Test1.exe Test1.o -L. -lesmf_time
else
Test1.exe : Test1.o
	$(FC) $(LDFLAGS) -o Test1.exe Test1.o
endif

#.F90.o :
#	$(RM) $@
#	$(SED_FTN) $*.F90 > $*.b
#	$(CPP) -C -P -I. $*.b > $*.f
#	$(RM) $*.b
#	$(FC) -c $*.f
#
#.F90.f :
#	$(RM) $@
#	$(SED_FTN) $*.F90 > $*.b
#	$(CPP) -C -P -I. $*.b > $*.f
#	$(RM) $*.b
#
#.f.o :
#	$(RM) $@
#	$(RM) $*.b
#	$(FC) -c $*.f

superclean: 
	/bin/rm -f *.b *.f *.o libesmf_time.a *.mod Test1.exe Test1.out *.f90

# DEPENDENCIES : only dependencies after this line 

#$$$ update dependencies!  

ifeq ($(ESMF_BLD),$(null))

include Depends

else

endif

Test1.o :  Test1.F90


