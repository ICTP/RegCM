#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#    
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include $(REGCM_ROOT)/makeinc

vpath %.f90 $(REGCM_ROOT)/external/aerosols
vpath %.F90 $(REGCM_ROOT)/external/cosp2

noinst_LIBRARIES = libaerosol.a libbmi.a libcosp.a

libaerosol_a_SOURCES = mo_simple_plumes.f90
libbmi_a_SOURCES = bmi.f90

# COSP2 satellite simulators (CFMIP, BSD license, see cosp2/), ISCCP and
# MODIS used by the model (icosp = 1); the other simulators only to link.
libcosp_a_SOURCES = cosp_kinds.F90 cosp_constants.F90 \
  cosp_errorHandling.F90 cosp_config.F90 cosp_rttov_utilSTUB.F90 \
  cosp_stats.F90 cosp_rttovSTUB.F90 cosp_rttov_interfaceSTUB.F90 \
  icarus.F90 modis_simulator.F90 MISR_simulator.F90 parasol.F90 \
  lidar_simulator.F90 quickbeam.F90 cosp_isccp_interface.F90 \
  cosp_modis_interface.F90 cosp_misr_interface.F90 \
  cosp_parasol_interface.F90 cosp_calipso_interface.F90 \
  cosp_cloudsat_interface.F90 cosp_grLidar532_interface.F90 \
  cosp_atlid_interface.F90 cosp.F90 mo_rng.F90 scops.F90 cosp_optics.F90

# COSP2 has free form lines longer than 132 characters
if COMPILER_GNU
COSPFLAGS = -ffree-line-length-none
else
COSPFLAGS =
endif

%.o: %.c
	$(CC) $(CFLAGS) -c $<

%.o: %.f90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) $(COSPFLAGS) -c $<

cosp_constants.o : cosp_kinds.o
cosp_errorHandling.o : cosp_kinds.o
cosp_config.o : cosp_kinds.o
cosp_rttov_utilSTUB.o : cosp_kinds.o
cosp_stats.o : cosp_kinds.o cosp_constants.o cosp_config.o \
  cosp_rttov_utilSTUB.o
cosp_rttovSTUB.o : cosp_kinds.o cosp_constants.o cosp_config.o
cosp_rttov_interfaceSTUB.o : cosp_kinds.o cosp_rttovSTUB.o \
  cosp_rttov_utilSTUB.o
icarus.o : cosp_kinds.o cosp_constants.o cosp_config.o cosp_stats.o
modis_simulator.o : cosp_kinds.o cosp_config.o cosp_stats.o
MISR_simulator.o : cosp_kinds.o cosp_config.o cosp_stats.o
parasol.o : cosp_kinds.o cosp_constants.o cosp_config.o
lidar_simulator.o : cosp_kinds.o cosp_config.o cosp_stats.o
quickbeam.o : cosp_kinds.o cosp_config.o cosp_stats.o
cosp_isccp_interface.o : cosp_kinds.o cosp_stats.o icarus.o
cosp_modis_interface.o : cosp_kinds.o cosp_config.o cosp_stats.o \
  modis_simulator.o
cosp_misr_interface.o : cosp_kinds.o cosp_stats.o
cosp_parasol_interface.o : cosp_kinds.o cosp_stats.o
cosp_calipso_interface.o : cosp_kinds.o cosp_stats.o lidar_simulator.o
cosp_cloudsat_interface.o : cosp_kinds.o cosp_stats.o quickbeam.o
cosp_grLidar532_interface.o : cosp_kinds.o
cosp_atlid_interface.o : cosp_kinds.o cosp_stats.o
cosp.o : cosp_kinds.o cosp_constants.o cosp_config.o cosp_stats.o \
  cosp_rttovSTUB.o cosp_rttov_interfaceSTUB.o cosp_rttov_utilSTUB.o \
  icarus.o modis_simulator.o MISR_simulator.o parasol.o \
  lidar_simulator.o quickbeam.o cosp_isccp_interface.o \
  cosp_modis_interface.o cosp_misr_interface.o cosp_parasol_interface.o \
  cosp_calipso_interface.o cosp_cloudsat_interface.o \
  cosp_grLidar532_interface.o cosp_atlid_interface.o
mo_rng.o : cosp_kinds.o
scops.o : cosp_kinds.o cosp_errorHandling.o mo_rng.o
cosp_optics.o : cosp_kinds.o cosp_constants.o modis_simulator.o

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt mpif.h

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
