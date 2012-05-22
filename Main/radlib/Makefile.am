#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


SUBDIRS = RRTMG_SW
RRTMGSW = RRTMG_SW/librrtmgsw.a

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

noinst_LIBRARIES = librad.a
librad_a_SOURCES = mod_rad_colmod3.F90 mod_rad_o3blk.F90 mod_rad_outrad.F90 \
  mod_rad_common.f90 mod_rad_radiation.F90 mod_rad_tracer.f90 \
  mod_rad_aerosol.F90 mod_rad_scenarios.F90 mod_atm_msis.f90 \
  mod_rrtmg_driver.f90 rrtmg_sw_cldprmc.f90  rrtmg_sw_read_nc.f90 \
  rrtmg_sw_spcvmc.f90  rrtmg_sw_vrtqdr.f90 rrtmg_sw_cldprop.f90 \
  rrtmg_sw_rad.f90 rrtmg_sw_reftra.f90 rrtmg_sw_spcvrt.f90 \
  rrtmg_sw_init.f90 rrtmg_sw_rad.nomcica.f90  rrtmg_sw_setcoef.f90 \
  rrtmg_sw_taumol.f90 mcica_random_numbers.f90  mcica_subcol_gen_sw.f90\
  rrtmg_lw_cldprop.f90  rrtmg_lw_rad.f90     rrtmg_lw_rtrn.f90 \
  rrtmg_lw_setcoef.f90 mcica_subcol_gen_lw.f90   rrtmg_lw_init.f90 \
  rrtmg_lw_rad.nomcica.f90 rrtmg_lw_rtrnmc.f90  rrtmg_lw_taumol.f90  \
  rrtmg_lw_cldprmc.f90  rrtmg_lw_read_nc.f90    rrtmg_lw_rtrnmr.f90 \
	mod_rad_atmosphere.f90

if COMPILER_SUN
  CPPFLAGS += -MRRTMG_SW
else
  CPPFLAGS += -IRRTMG_SW
endif

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.f90

clean-local:
	rm -f *.mod *__genmod.f90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_atm_msis.o mod_atm_msis.mod: mod_atm_msis.f90 \
  $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_constants.mod \
  $(REGCM_ROOT)/Share/mod_spline.mod
mod_rad_atmosphere.o mod_rad_atmosphere.mod: mod_rad_atmosphere.f90 \
  $(REGCM_ROOT)/Share/mod_realkinds.mod
mod_rad_aerosol.o mod_rad_aerosol.mod: mod_rad_aerosol.F90 \
  $(REGCM_ROOT)/Share/mod_realkinds.mod $(MPP)/mod_mpmessage.mod \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod mod_rad_common.o
mod_rad_common.o mod_rad_common.mod: mod_rad_common.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod
mod_rad_o3blk.o mod_rad_o3blk.mod: mod_rad_o3blk.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(REGCM_ROOT)/Share/mod_memutil.mod $(MPP)/mod_mppparam.mod mod_rad_common.o
mod_rad_colmod3.o mod_rad_colmod3.mod: mod_rad_colmod3.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(MPP)/mod_service.mod \
  mod_rad_common.o mod_rad_radiation.o mod_rad_outrad.o mod_rrtmg_driver.o
mod_rad_outrad.o mod_rad_outrad.mod: mod_rad_outrad.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_rad_common.o
mod_rad_scenarios.o mod_rad_scenarios.mod: mod_rad_scenarios.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(MPP)/mod_mpmessage.mod
mod_rad_tracer.o mod_rad_tracer.mod: mod_rad_tracer.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod
mod_rad_radiation.o mod_rad_radiation.mod: mod_rad_radiation.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(MPP)/mod_mpmessage.mod $(MPP)/mod_service.mod \
  $(REGCM_ROOT)/Share/mod_memutil.mod mod_rad_tracer.o mod_rad_scenarios.o \
  mod_rad_aerosol.o
mod_rrtmg_driver.o mod_rrtmg_driver.mod: mod_rrtmg_driver.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(RRTMGSW) rrtmg_sw_rad.o rrtmg_lw_rad.o mcica_subcol_gen_sw.o mcica_subcol_gen_lw.o mod_rad_common.o \
  mod_rad_scenarios.o mod_rad_tracer.o mod_rad_outrad.o 
rrtmg_sw_cldprmc.o rrtmg_sw_cldprmc.mod : rrtmg_sw_cldprmc.f90 $(RRTMGSW)
rrtmg_sw_read_nc.o rrtmg_sw_read_nc.mod : rrtmg_sw_read_nc.f90  $(RRTMGSW)
rrtmg_sw_spcvmc.o  rrtmg_sw_spcvmc.mod : rrtmg_sw_spcvmc.f90   $(RRTMGSW) rrtmg_sw_reftra.o \
	rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o
rrtmg_sw_vrtqdr.o  rrtmg_sw_vrtqdr.mod : rrtmg_sw_vrtqdr.f90  $(RRTMGSW)
rrtmg_sw_cldprop.o  rrtmg_sw_cldprop.mod :rrtmg_sw_cldprop.f90  $(RRTMGSW)
rrtmg_sw_rad.o rrtmg_sw_rad.mod : rrtmg_sw_rad.f90  $(RRTMGSW) mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o \
	rrtmg_sw_setcoef.o rrtmg_sw_spcvmc.o
rrtmg_sw_reftra.o  rrtmg_sw_reftra.mod : rrtmg_sw_reftra.f90   $(RRTMGSW)
rrtmg_sw_spcvrt.o rrtmg_sw_spcvrt.mod : rrtmg_sw_spcvrt.f90  $(RRTMGSW)
rrtmg_sw_init.o rrtmg_sw_init.mod : rrtmg_sw_init.f90  $(RRTMGSW)
rrtmg_sw_rad.nomcica.o  rrtmg_sw_rad.nomcica.mod : rrtmg_sw_rad.nomcica.f90  $(RRTMGSW)
rrtmg_sw_setcoef.o rrtmg_sw_setcoef.mod  : rrtmg_sw_setcoef.f90   $(RRTMGSW)
rrtmg_sw_taumol.o rrtmg_sw_taumol.mod : rrtmg_sw_taumol.f90   $(RRTMGSW)
mcica_random_numbers.o mcica_random_numbers.mod : mcica_random_numbers.f90   $(RRTMGSW)
mcica_subcol_gen_sw.o mcica_subcol_gen_sw.mod : mcica_subcol_gen_sw.f90  $(RRTMGSW) mcica_random_numbers.o
rrtmg_lw_cldprop.o rrtmg_lw_cldprop.mod:  rrtmg_lw_cldprop.f90  $(RRTMGSW)
rrtmg_lw_rad.o  rrtmg_lw_rad.mod: rrtmg_lw_rad.f90  $(RRTMGSW)  mcica_subcol_gen_lw.o rrtmg_lw_cldprmc.o \
	rrtmg_lw_setcoef.o  rrtmg_lw_rtrnmc.o rrtmg_lw_taumol.o
rrtmg_lw_rtrn.o rrtmg_lw_rtrn.mod:  rrtmg_lw_rtrn.f90  $(RRTMGSW)
rrtmg_lw_setcoef.o rrtmg_lw_setcoef.mod:  rrtmg_lw_setcoef.f90  $(RRTMGSW)
mcica_subcol_gen_lw.o mcica_subcol_gen_lw.mod: mcica_subcol_gen_lw.f90  $(RRTMGSW)  mcica_random_numbers.o
rrtmg_lw_init.o   rrtmg_lw_init.mod:  rrtmg_lw_init.f90  $(RRTMGSW) rrtmg_lw_setcoef.o
rrtmg_lw_rad.nomcica.o rrtmg_lw_rad.nomcica.mod: rrtmg_lw_rad.nomcica.f90  $(RRTMGSW) rrtmg_lw_rtrnmr.o rrtmg_lw_rtrn.o rrtmg_lw_setcoef.o \
	rrtmg_lw_taumol.o
rrtmg_lw_rtrnmc.o  rrtmg_lw_rtrnmc.mod:  rrtmg_lw_rtrnmc.f90   $(RRTMGSW) 
rrtmg_lw_taumol.o  rrtmg_lw_taumol.mod:  rrtmg_lw_taumol.f90   $(RRTMGSW)
rrtmg_lw_cldprmc.o rrtmg_lw_cldprmc.mod: rrtmg_lw_cldprmc.f90 $(RRTMGSW)
rrtmg_lw_read_nc.o    rrtmg_lw_read_nc.mod:  rrtmg_lw_read_nc.f90  $(RRTMGSW)
rrtmg_lw_rtrnmr.o   rrtmg_lw_rtrnmr.mod:  rrtmg_lw_rtrnmr.f90  $(RRTMGSW)
