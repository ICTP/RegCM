#-----------------------------------------------------------------------
# This Makefile is for doing csm_share unit testing
#------------------------------------------------------------------------
cpp_dirs := . ../shr ../eshr ../../utils/esmf_wrf_timemgr
ifneq ($(SPMD),TRUE)
cpp_dirs += ../../utils/mct/mpi-serial
endif
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Platform specific macros
include make.Macros

space := $(null) $(null)

ifneq ($(ESMF_BLD),$(null))
cpp_dirs += $(ESMF_LIB)
endif

# Expand any tildes in directory names.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
# Change spaces to colons.
VPATH    := $(subst $(space),:,$(VPATH))

all: test_eshr_timemgr


ifneq ($(SPMD),TRUE)
   OBJS_NOMPI   := fort.o group.o collective.o comm.o list.o handles.o mpi.o recv.o req.o \
	      		send.o time.o
else
   OBJS_NOMPI   := $(null)
endif

OBJS       := test_shr_sys.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o \
	      $(OBJS_NOMPI)
OBJS_FILE  := test_shr_file.o shr_sys_mod.o shr_kind_mod.o shr_file_mod.o shr_mpi_mod.o \
	      $(OBJS_NOMPI)
OBJS_IN    := test_shr_inputinfo.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o \
		shr_cal_mod.o shr_file_mod.o shr_string_mod.o \
		shr_inputinfo_mod.o shr_ncio_mod.o $(OBJS_NOMPI)
OBJS_TIME  := test_eshr_timemgr.o eshr_timemgr_mod.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o \
		shr_cal_mod.o shr_file_mod.o shr_string_mod.o \
		shr_inputinfo_mod.o shr_ncio_mod.o shr_orb_mod.o eshr_rc_mod.o \
		eshr_estate_mod.o $(OBJS_NOMPI)
OBJS_ORB := test_shr_orb.o shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o shr_mpi_mod.o \
	    shr_const_mod.o $(OBJS_NOMPI)
OBJS_ESTATE := test_eshr_estate.o eshr_inputinfo_mod.o eshr_timemgr_mod.o \
		eshr_estate_mod.o eshr_rc_mod.o shr_sys_mod.o shr_kind_mod.o \
                shr_mpi_mod.o shr_const_mod.o shr_string_mod.o \
                shr_inputinfo_mod.o shr_ncio_mod.o shr_orb_mod.o shr_file_mod.o \
                shr_cal_mod.o gridndecomp_mod.o $(OBJS_NOMPI)
OBJS_MPI := test_shr_mpi.o shr_mpi_mod.o shr_kind_mod.o shr_sys_mod.o shr_const_mod.o $(OBJS_NOMPI)

WRFESMF_OBJS := ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
                Meat.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
                ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF_Mod.o \
                ESMF_AlarmClockMod.o wrf_error_fatal.o wrf_message.o

ifeq ($(ESMF_BLD),$(null))
   OBJS_TIME   += $(WRFESMF_OBJS)
   OBJS_ESTATE += $(WRFESMF_OBJS)
endif

#
# Executables:
#
debug:
	@echo "VPATH:      " $(VPATH)
	@echo "ESMF_MOD:   " $(ESMF_MOD)
	@echo "ESMF_ARCH:  " $(ESMF_ARCH)
	@echo "FC:         " $(FC)
	@echo "INC_NETCDF: " $(INC_NETCDF)
	@echo "LIB_MPI:    " $(LIB_MPI)
test_shr_sys: $(OBJS)
	$(LD) -o test_shr_sys     $(OBJS) $(LDFLAGS)
test_shr_file: $(OBJS_FILE)
	$(LD) -o test_shr_file    $(OBJS_FILE) $(LDFLAGS)
test_shr_orb: $(OBJS_ORB)
	$(LD) -o test_shr_orb     $(OBJS_ORB) $(LDFLAGS)
test_shr_inputinfo: $(OBJS_IN)
	$(LD) -o test_shr_inputinfo   $(OBJS_IN) $(LDFLAGS)
test_shr_mpi: $(OBJS_MPI)
	$(LD) -o test_shr_mpi $(OBJS_MPI) $(LDFLAGS)
test_eshr_estate: $(OBJS_ESTATE)
	$(LD) -o test_eshr_estate $(OBJS_ESTATE) $(LDFLAGS)
test_estate: test_estate.o gridndecomp_mod.o shr_sys_mod.o shr_kind_mod.o eshr_rc_mod.o \
             shr_mpi_mod.o
	$(LD) -o test_estate test_estate.o gridndecomp_mod.o \
		   shr_sys_mod.o shr_kind_mod.o eshr_rc_mod.o shr_mpi_mod.o $(LDFLAGS)
test_eshr_timemgr: $(OBJS_TIME)
	$(LD) -o test_eshr_timemgr $(OBJS_TIME) $(LDFLAGS)

clean:
	$(RM) *.mod *.o *.f *.f90 test_shr_inputinfo test_shr_sys test_shr_orb test_eshr_timemgr \
		libesmf.a test_eshr_estate

#
# Dependencies
#

# Test programs
test_shr_file.o: test_shr_file.F90 shr_file_mod.o shr_sys_mod.o shr_kind_mod.o
test_eshr_timemgr.o: test_eshr_timemgr.F90 eshr_timemgr_mod.o shr_inputinfo_mod.o \
		      shr_kind_mod.o shr_sys_mod.o shr_ncio_mod.o eshr_estate_mod.o
test_eshr_estate.o: test_eshr_estate.F90 eshr_timemgr_mod.o eshr_inputinfo_mod.o \
		    eshr_estate_mod.o shr_inputinfo_mod.o shr_kind_mod.o shr_sys_mod.o \
		    shr_ncio_mod.o gridndecomp_mod.o
test_estate.o: test_estate.F90 gridndecomp_mod.o shr_sys_mod.o eshr_rc_mod.o shr_kind_mod.o
test_shr_inputinfo.o: shr_inputinfo_mod.o test_shr_inputinfo.F90
test_shr_sys.o: shr_sys_mod.o shr_sys_mod.o test_shr_sys.F90
test_shr_orb.o: shr_orb_mod.o test_shr_orb.F90
test_shr_mpi.o: test_shr_mpi.F90 shr_mpi_mod.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o

# Test modules
gridndecomp_mod.o: gridndecomp_mod.F90 shr_sys_mod.o shr_kind_mod.o eshr_rc_mod.o

# share codes
shr_alarm_mod.o : shr_alarm_mod.F90 shr_cal_mod.o shr_date_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_cal_mod.o : shr_cal_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_const_mod.o : shr_const_mod.F90 shr_kind_mod.o 
shr_date_mod.o : shr_date_mod.F90 shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_file_mod.o : shr_file_mod.F90 shr_sys_mod.o shr_kind_mod.o 
shr_inputinfo_mod.o : shr_inputinfo_mod.F90 shr_file_mod.o shr_string_mod.o shr_sys_mod.o shr_ncio_mod.o shr_kind_mod.o 
shr_kind_mod.o : shr_kind_mod.F90  
shr_map_mod.o : shr_map_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_mpi_mod.o : shr_mpi_mod.F90 shr_kind_mod.o 
shr_msg_mod.o : shr_msg_mod.F90 shr_file_mod.o
shr_ncio_mod.o : shr_ncio_mod.F90 shr_file_mod.o shr_string_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_orb_mod.o : shr_orb_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_string_mod.o : shr_string_mod.F90 shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_sys_mod.o : shr_sys_mod.F90 shr_mpi_mod.o shr_kind_mod.o 
shr_tInterp_mod.o : shr_tInterp_mod.F90 shr_const_mod.o shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o 
shr_timer_mod.o : shr_timer_mod.F90 shr_sys_mod.o shr_kind_mod.o 
shr_mpi_mod.o : shr_kind_mod.F90
ifeq ($(ESMF_BLD),$(null))
  eshr_timemgr_mod.o : eshr_timemgr_mod.F90 shr_file_mod.o shr_string_mod.o \
			shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o shr_ncio_mod.o \
			eshr_rc_mod.o libesmf.a eshr_estate_mod.o
  eshr_inputinfo_mod.o : eshr_inputinfo_mod.F90 shr_inputinfo_mod.o shr_sys_mod.o \
			libesmf.a eshr_estate_mod.o
  eshr_estate_mod.o : eshr_estate_mod.F90 libesmf.a
  eshr_rc_mod.o : shr_sys_mod.o libesmf.a
else
  eshr_timemgr_mod.o : eshr_timemgr_mod.F90 shr_file_mod.o shr_string_mod.o \
			shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o shr_ncio_mod.o \
			eshr_rc_mod.o eshr_estate_mod.o
  eshr_inputinfo_mod.o : eshr_inputinfo_mod.F90 shr_inputinfo_mod.o shr_sys_mod.o \
			eshr_estate_mod.o
  eshr_estate_mod.o : eshr_estate_mod.F90
  eshr_rc_mod.o : shr_sys_mod.o
endif


# ESMF code...
ifeq ($(ESMF_BLD),$(null))

include ../../utils/esmf_wrf_timemgr/unittests/Depends
AR  := ar
CPP := cpp

libesmf.a : $(WRFESMF_OBJS)
	$(RM) -f libesmf.a
	$(AR) $(ARFLAGS) libesmf.a $(WRFESMF_OBJS)
	$(RANLIB) libesmf.a

endif
wrf_error_fatal.o : wrf_error_fatal.F90 shr_sys_mod.o
wrf_message.o : wrf_message.F90
