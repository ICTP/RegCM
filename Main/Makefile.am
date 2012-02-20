
SUBDIRS = netlib batslib chemlib clmlib cumlib radlib pbllib

include $(REGCM_ROOT)/makeinc

CLMLIBPTH = ./clmlib
BATLIBPTH = ./batslib
NETLIBPTH = ./netlib
CHELIBPTH = ./chemlib
CUMLIBPTH = ./cumlib
RADLIBPTH = ./radlib
PBLLIBPTH = ./pbllib
CBMZLIBPTH = chemlib/GAS_CBMZ
RRTMSWLIBPTH = radlib/RRTMG_SW


CLMLIB = $(CLMLIBPTH)/libclm.a
BATLIB = $(BATLIBPTH)/libbats.a
NETLIB = $(NETLIBPTH)/libnetlib.a
CHELIB = $(CHELIBPTH)/libchem.a
CUMLIB = $(CUMLIBPTH)/libcum.a
RADLIB = $(RADLIBPTH)/librad.a
PBLLIB = $(PBLLIBPTH)/libpbl.a
CBMZLIB = $(CBMZLIBPTH)/libcbmz.a
RRTMSWLIB = $(RRTMSWLIBPTH)/librrtmgsw.a


if REGCM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif  

if DO_COMPILE_LIBCLM
rcaddsuffix = _clm
endif

if DO_COMPILE_BAND
rcaddsuffix = _band
endif

program_transform_name = s&$$&$(rcsuffix)$(rcaddsuffix)&

LBATLIB = -L$(BATLIBPTH) -lbats
LCLMLIB = -L$(CLMLIBPTH) -lclm
LNETLIB = -L$(NETLIBPTH) -lnetlib
LCHELIB = -L$(CHELIBPTH) -lchem
LCUMLIB = -L$(CUMLIBPTH) -lcum
LRADLIB = -L$(RADLIBPTH) -lrad
LPBLLIB = -L$(PBLLIBPTH) -lpbl
LCBMZLIB = -L$(CBMZLIBPTH) -lcbmz
LRRTMSWLIB = -L$(RRTMSWLIBPTH) -lrrtmgsw

if COMPILER_SUN
  CPPFLAGS += -M$(CLMLIBPTH) -M$(NETLIBPTH) -M$(CBMZLIBPTH) -M$(CHELIBPTH) \
              -M$(BATLIBPTH) -M$(CUMLIBPTH) -M$(RADLIBPTH) -M$(PBLLIBPTH) -M$(RRTMSWLIBPTH)
else
  CPPFLAGS += -I$(CLMLIBPTH) -I$(NETLIBPTH) -I$(CBMZLIBPTH) -I$(CHELIBPTH) \
              -I$(BATLIBPTH) -I$(CUMLIBPTH) -I$(RADLIBPTH) -I$(PBLLIBPTH) -I$(RRTMSWLIBPTH)
endif

LIBS += $(LCLMLIB) $(LNETLIB) $(LCHELIB) $(LCBMZLIB) $(LBATLIB) $(LCUMLIB) \
        $(LRADLIB) $(LPBLLIB) $(LRRTMSWLIB)

bin_PROGRAMS = regcm

regcm_SOURCES = mod_diagnosis.F90 mod_header.F90 mod_cu_interface.f90 \
                mod_advection.F90 mod_bdycod.F90 mod_diffusion.F90 \
                mod_che_interface.F90 mod_rad_interface.f90 \
                mod_pbl_interface.F90 mod_init.F90 mod_atm_interface.f90 \
                mod_precip.f90 mod_mppio.F90 mod_ncio.F90 mod_output.F90 \
                mod_params.F90 mod_runparams.F90 mod_savefile.F90 \
                mod_slice.F90 mod_sun.F90 mod_tendency.F90 mod_tstep.f90 \
                mod_regcm_interface.F90 mod_vmodes.F90 mod_split.F90  \
                mod_lm_interface.F90 regcm.F90

FCLD = $(MPIFC)

#if DO_COMPILE_LIBESMF
#LIBS += $(ESMFLIBPTH) 
#CPPFLAGS += $(ESMFINCPTH) -I$(ROMS_PREFIX)
#regcm_SOURCES += mod_esmf_atm.F90 mod_esmf_ocn.F90 \
#                 mod_esmf_cpl.F90 mod_couplerr.F90
#skip_files  = $(ROMS_PREFIX)/master.o
#regcm_LDADD = $(filter-out $(skip_files), $(wildcard $(ROMS_PREFIX)/*.o))
#endif

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

#
# Modules
#
mod_advection.o: mod_advection.F90 mod_atm_interface.o mod_runparams.o
mod_atm_interface.o: mod_atm_interface.f90 mod_runparams.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_bdycod.o: mod_bdycod.F90 mod_runparams.o mod_atm_interface.o \
  mod_lm_interface.o mod_che_interface.o mod_ncio.o mod_mppio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_che_interface.o: mod_che_interface.F90 $(CHELIB)
mod_rad_interface.o: mod_rad_interface.f90 mod_lm_interface.o \
  mod_atm_interface.o mod_che_interface.o $(RADLIB)
mod_cu_interface.o: mod_cu_interface.f90 mod_atm_interface.o $(CUMLIB)
mod_diagnosis.o: mod_diagnosis.F90 mod_runparams.o mod_atm_interface.o \
  mod_che_interface.o mod_mppio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_diffusion.o: mod_diffusion.F90 mod_runparams.o
mod_gridfuncs.o: mod_gridfuncs.f90 mod_atm_interface.o mod_runparams.o
mod_header.o: mod_header.F90
mod_init.o: mod_init.F90 mod_bdycod.o mod_pbl_interface.o \
  mod_diagnosis.o mod_atm_interface.o mod_precip.o \
  mod_runparams.o mod_cu_interface.o mod_rad_interface.o \
  mod_mppio.o mod_ncio.o mod_che_interface.o \
  mod_lm_interface.o mod_sun.o mod_savefile.o $(RADLIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_lm_interface.o: mod_lm_interface.F90 mod_atm_interface.o $(BATLIB) $(CLMLIB)
mod_mppio.o: mod_mppio.F90 mod_runparams.o mod_pbl_interface.o \
  mod_lm_interface.o mod_che_interface.o mod_atm_interface.o \
  mod_rad_interface.o mod_cu_interface.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_ncio.o: mod_ncio.F90 mod_che_interface.o mod_cu_interface.o \
  mod_lm_interface.o mod_runparams.o mod_rad_interface.o
mod_output.o: mod_output.F90 mod_bdycod.o mod_cu_interface.o \
  mod_rad_interface.o mod_che_interface.o mod_atm_interface.o mod_precip.o \
  mod_runparams.o mod_diagnosis.o mod_lm_interface.o \
  mod_split.o mod_mppio.o mod_ncio.o mod_savefile.o $(CHELIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_params.o: mod_params.F90 mod_lm_interface.o mod_cu_interface.o \
  mod_atm_interface.o mod_runparams.o mod_mppio.o mod_precip.o \
  mod_che_interface.o mod_split.o mod_pbl_interface.o mod_ncio.o \
  mod_rad_interface.o mod_bdycod.o mod_tendency.o mod_che_interface.o mod_ncio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_precip.o: mod_precip.f90 mod_slice.o mod_atm_interface.o mod_runparams.o
mod_runparams.o: mod_runparams.F90
mod_savefile.o: mod_savefile.F90 mod_runparams.o mod_che_interface.o \
  mod_atm_interface.o mod_bdycod.o mod_diagnosis.o mod_mppio.o \
  mod_lm_interface.o mod_rad_interface.o
mod_slice.o: mod_slice.F90 mod_atm_interface.o mod_che_interface.o \
  mod_runparams.o mod_pbl_interface.o
mod_split.o: mod_split.F90 mod_atm_interface.o mod_runparams.o \
  mod_vmodes.o mod_savefile.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_sun.o: mod_sun.F90 mod_atm_interface.o
mod_pbl_interface.o: mod_pbl_interface.F90 mod_atm_interface.o $(PBLLIB)
mod_tendency.o: mod_tendency.F90 mod_bdycod.o mod_cu_interface.o \
  mod_atm_interface.o mod_precip.o mod_runparams.o mod_rad_interface.o \
  mod_slice.o mod_sun.o mod_precip.o mod_pbl_interface.o \
  mod_diagnosis.o mod_diffusion.o mod_advection.o mod_che_interface.o \
  mod_lm_interface.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_tstep.o: mod_tstep.f90 mod_runparams.o
mod_vmodes.o: mod_vmodes.F90 mod_runparams.o
#
# ESMF modules 
#
#if DO_COMPILE_LIBESMF
#mod_couplerr.o: mod_couplerr.F90
#mod_esmf_atm.o: mod_esmf_atm.F90 mod_regcm_interface.o mod_couplerr.o
#mod_esmf_ocn.o: mod_esmf_ocn.F90 mod_regcm_interface.o mod_couplerr.o
#mod_esmf_cpl.o: mod_esmf_cpl.F90 mod_regcm_interface.o mod_couplerr.o
#endif
#
# RegCM Main Program
#
mod_regcm_interface.o: mod_regcm_interface.F90 mod_runparams.o mod_init.o \
  mod_ncio.o mod_output.o mod_split.o mod_header.o mod_params.o mod_tendency.o \
  mod_tstep.o mod_mppio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
#if DO_COMPILE_LIBESMF
#regcm.o: regcm.F90 regcm_nonesmf.F90 regcm_esmf.F90 \
#         mod_couplerr.o mod_esmf_atm.o mod_esmf_ocn.o mod_esmf_cpl.o
#	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
#else
regcm.o: regcm.F90 regcm_nonesmf.F90 regcm_esmf.F90
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
#endif
