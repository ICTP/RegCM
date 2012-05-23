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

vpath %.mod $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:RRTMG_SW/:.

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

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

distclean-local:
	rm -f *.mod *__genmod.f90

clean-local:
	rm -f *.mod *__genmod.f90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mcica_random_numbers.o mcica_random_numbers.mod: mcica_random_numbers.f90 parkind.mod parkind.mod
mcica_subcol_gen_lw.o mcica_subcol_gen_lw.mod: mcica_subcol_gen_lw.f90 parkind.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod rrlw_vsn.mod mcica_random_numbers.mod
mcica_subcol_gen_sw.o mcica_subcol_gen_sw.mod: mcica_subcol_gen_sw.f90 parkind.mod parrrsw.mod rrsw_con.mod rrsw_wvn.mod rrsw_vsn.mod mcica_random_numbers.mod
mod_atm_msis.o mod_atm_msis.mod: mod_atm_msis.f90 mod_constants.mod mod_spline.mod mod_realkinds.mod
mod_rad_aerosol.o mod_rad_aerosol.mod: mod_rad_aerosol.F90 mod_realkinds.mod mod_dynparam.mod mod_constants.mod mod_memutil.mod mod_mpmessage.mod mod_rad_common.mod
mod_rad_atmosphere.o mod_rad_atmosphere.mod: mod_rad_atmosphere.f90 mod_realkinds.mod
mod_rad_colmod3.o mod_rad_colmod3.mod: mod_rad_colmod3.F90 mod_realkinds.mod mod_dynparam.mod mod_service.mod mod_rad_radiation.mod mod_rad_common.mod mod_rad_outrad.mod mod_rrtmg_driver.mod mod_rad_aerosol.mod mod_runparams.mod
mod_rad_common.o mod_rad_common.mod: mod_rad_common.f90 mod_constants.mod mod_dynparam.mod mod_realkinds.mod mod_memutil.mod
mod_rad_o3blk.o mod_rad_o3blk.mod: mod_rad_o3blk.F90 mod_constants.mod mod_date.mod mod_interp.mod mod_vertint.mod mod_mppparam.mod mod_mpmessage.mod mod_dynparam.mod mod_memutil.mod mod_rad_common.mod mod_stdio.mod
mod_rad_outrad.o mod_rad_outrad.mod: mod_rad_outrad.F90 mod_dynparam.mod mod_mpmessage.mod mod_rad_common.mod
mod_rad_radiation.o mod_rad_radiation.mod: mod_rad_radiation.F90 mod_dynparam.mod mod_mpmessage.mod mod_service.mod mod_memutil.mod mod_rad_common.mod mod_rad_tracer.mod mod_rad_scenarios.mod mod_rad_aerosol.mod mod_runparams.mod
mod_rad_scenarios.o mod_rad_scenarios.mod: mod_rad_scenarios.F90 mod_dynparam.mod mod_mpmessage.mod mod_realkinds.mod
mod_rad_tracer.o mod_rad_tracer.mod: mod_rad_tracer.f90 mod_constants.mod mod_dynparam.mod mod_realkinds.mod
mod_rrtmg_driver.o mod_rrtmg_driver.mod: mod_rrtmg_driver.f90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod mod_rad_common.mod mod_rad_scenarios.mod mod_rad_tracer.mod mod_rad_aerosol.mod rrtmg_sw_rad.mod mcica_subcol_gen_sw.mod parrrsw.mod rrsw_wvn.mod parrrtm.mod rrtmg_lw_rad.mod mod_rad_outrad.mod
rrtmg_lw_cldprmc.o rrtmg_lw_cldprmc.mod: rrtmg_lw_cldprmc.f90 parkind.mod parrrtm.mod rrlw_cld.mod rrlw_wvn.mod rrlw_vsn.mod
rrtmg_lw_cldprop.o rrtmg_lw_cldprop.mod: rrtmg_lw_cldprop.f90 parkind.mod parrrtm.mod rrlw_cld.mod rrlw_vsn.mod
rrtmg_lw_init.o rrtmg_lw_init.mod: rrtmg_lw_init.f90 parkind.mod rrlw_wvn.mod rrtmg_lw_setcoef.mod parrrtm.mod rrlw_tbl.mod rrlw_vsn.mod parrrtm.mod rrlw_con.mod rrlw_vsn.mod parrrtm.mod rrlw_kg01.mod parrrtm.mod rrlw_kg02.mod parrrtm.mod rrlw_kg03.mod parrrtm.mod rrlw_kg04.mod parrrtm.mod rrlw_kg05.mod parrrtm.mod rrlw_kg06.mod parrrtm.mod rrlw_kg07.mod parrrtm.mod rrlw_kg08.mod parrrtm.mod rrlw_kg09.mod parrrtm.mod rrlw_kg10.mod parrrtm.mod rrlw_kg11.mod parrrtm.mod rrlw_kg12.mod parrrtm.mod rrlw_kg13.mod parrrtm.mod rrlw_kg14.mod parrrtm.mod rrlw_kg15.mod parrrtm.mod rrlw_kg16.mod rrlw_cld.mod
rrtmg_lw_k_g.o rrtmg_lw_k_g.mod: rrtmg_lw_k_g.f90 parkind.mod rrlw_kg01.mod rrlw_vsn.mod parkind.mod rrlw_kg02.mod parkind.mod rrlw_kg03.mod parkind.mod rrlw_kg04.mod parkind.mod rrlw_kg05.mod parkind.mod rrlw_kg06.mod parkind.mod rrlw_kg07.mod parkind.mod rrlw_kg08.mod parkind.mod rrlw_kg09.mod parkind.mod rrlw_kg10.mod parkind.mod rrlw_kg11.mod parkind.mod rrlw_kg12.mod parkind.mod rrlw_kg13.mod parkind.mod rrlw_kg14.mod parkind.mod rrlw_kg15.mod parkind.mod rrlw_kg16.mod
rrtmg_lw_rad.o rrtmg_lw_rad.mod: rrtmg_lw_rad.f90 parkind.mod rrlw_vsn.mod mcica_subcol_gen_lw.mod rrtmg_lw_cldprmc.mod rrtmg_lw_rtrnmc.mod rrtmg_lw_setcoef.mod rrtmg_lw_taumol.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod
rrtmg_lw_rad.nomcica.o rrtmg_lw_rad.nomcica.mod: rrtmg_lw_rad.nomcica.f90 parkind.mod rrlw_vsn.mod rrtmg_lw_cldprop.mod rrtmg_lw_rtrn.mod rrtmg_lw_rtrnmr.mod rrtmg_lw_setcoef.mod rrtmg_lw_taumol.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod
rrtmg_lw_read_nc.o rrtmg_lw_read_nc.mod: rrtmg_lw_read_nc.f90 rrlw_kg01.mod rrlw_ncpar.mod rrlw_kg02.mod rrlw_ncpar.mod rrlw_kg03.mod rrlw_ncpar.mod rrlw_kg04.mod rrlw_ncpar.mod rrlw_kg05.mod rrlw_ncpar.mod rrlw_kg06.mod rrlw_ncpar.mod rrlw_kg07.mod rrlw_ncpar.mod rrlw_kg08.mod rrlw_ncpar.mod rrlw_kg09.mod rrlw_ncpar.mod rrlw_kg10.mod rrlw_ncpar.mod rrlw_kg11.mod rrlw_ncpar.mod rrlw_kg12.mod rrlw_ncpar.mod rrlw_kg13.mod rrlw_ncpar.mod rrlw_kg14.mod rrlw_ncpar.mod rrlw_kg15.mod rrlw_ncpar.mod rrlw_kg16.mod rrlw_ncpar.mod
rrtmg_lw_rtrn.o rrtmg_lw_rtrn.mod: rrtmg_lw_rtrn.f90 parkind.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod rrlw_tbl.mod rrlw_vsn.mod
rrtmg_lw_rtrnmc.o rrtmg_lw_rtrnmc.mod: rrtmg_lw_rtrnmc.f90 parkind.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod rrlw_tbl.mod rrlw_vsn.mod
rrtmg_lw_rtrnmr.o rrtmg_lw_rtrnmr.mod: rrtmg_lw_rtrnmr.f90 parkind.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod rrlw_tbl.mod rrlw_vsn.mod
rrtmg_lw_setcoef.o rrtmg_lw_setcoef.mod: rrtmg_lw_setcoef.f90 parkind.mod parrrtm.mod rrlw_wvn.mod rrlw_ref.mod rrlw_vsn.mod
rrtmg_lw_taumol.o rrtmg_lw_taumol.mod: rrtmg_lw_taumol.f90 parkind.mod parrrtm.mod rrlw_con.mod rrlw_wvn.mod rrlw_vsn.mod parrrtm.mod rrlw_kg01.mod parrrtm.mod rrlw_kg02.mod parrrtm.mod rrlw_ref.mod rrlw_kg03.mod parrrtm.mod rrlw_ref.mod rrlw_kg04.mod parrrtm.mod rrlw_ref.mod rrlw_kg05.mod parrrtm.mod rrlw_ref.mod rrlw_kg06.mod parrrtm.mod rrlw_ref.mod rrlw_kg07.mod parrrtm.mod rrlw_ref.mod rrlw_kg08.mod parrrtm.mod rrlw_ref.mod rrlw_kg09.mod parrrtm.mod rrlw_kg10.mod parrrtm.mod rrlw_kg11.mod parrrtm.mod rrlw_ref.mod rrlw_kg12.mod parrrtm.mod rrlw_ref.mod rrlw_kg13.mod parrrtm.mod rrlw_kg14.mod parrrtm.mod rrlw_ref.mod rrlw_kg15.mod parrrtm.mod rrlw_ref.mod rrlw_kg16.mod
rrtmg_sw_cldprmc.o rrtmg_sw_cldprmc.mod: rrtmg_sw_cldprmc.f90 parkind.mod parrrsw.mod rrsw_cld.mod rrsw_wvn.mod rrsw_vsn.mod
rrtmg_sw_cldprop.o rrtmg_sw_cldprop.mod: rrtmg_sw_cldprop.f90 parkind.mod parrrsw.mod rrsw_cld.mod rrsw_wvn.mod rrsw_vsn.mod
rrtmg_sw_init.o rrtmg_sw_init.mod: rrtmg_sw_init.f90 parkind.mod rrsw_wvn.mod rrtmg_sw_setcoef.mod parrrsw.mod rrsw_tbl.mod rrsw_vsn.mod rrsw_con.mod rrsw_vsn.mod rrsw_aer.mod rrsw_kg16.mod rrsw_kg17.mod rrsw_kg18.mod rrsw_kg19.mod rrsw_kg20.mod rrsw_kg21.mod rrsw_kg22.mod rrsw_kg23.mod rrsw_kg24.mod rrsw_kg25.mod rrsw_kg26.mod rrsw_kg27.mod rrsw_kg28.mod rrsw_kg29.mod rrsw_cld.mod
rrtmg_sw_k_g.o rrtmg_sw_k_g.mod: rrtmg_sw_k_g.f90 parkind.mod rrsw_kg16.mod rrsw_vsn.mod parkind.mod rrsw_kg17.mod parkind.mod rrsw_kg18.mod parkind.mod rrsw_kg19.mod parkind.mod rrsw_kg20.mod parkind.mod rrsw_kg21.mod parkind.mod rrsw_kg22.mod parkind.mod rrsw_kg23.mod parkind.mod rrsw_kg24.mod parkind.mod rrsw_kg25.mod parkind.mod rrsw_kg26.mod parkind.mod rrsw_kg27.mod parkind.mod rrsw_kg28.mod parkind.mod rrsw_kg29.mod
rrtmg_sw_rad.o rrtmg_sw_rad.mod: rrtmg_sw_rad.f90 parkind.mod rrsw_vsn.mod mcica_subcol_gen_sw.mod rrtmg_sw_cldprmc.mod rrtmg_sw_setcoef.mod rrtmg_sw_spcvmc.mod parrrsw.mod rrsw_aer.mod rrsw_con.mod rrsw_wvn.mod rrsw_con.mod parrrsw.mod rrsw_con.mod rrsw_wvn.mod
rrtmg_sw_rad.nomcica.o rrtmg_sw_rad.nomcica.mod: rrtmg_sw_rad.nomcica.f90 parkind.mod rrsw_vsn.mod rrtmg_sw_cldprop.mod rrtmg_sw_setcoef.mod rrtmg_sw_spcvrt.mod parrrsw.mod rrsw_aer.mod rrsw_con.mod rrsw_wvn.mod rrsw_con.mod parrrsw.mod rrsw_con.mod rrsw_wvn.mod
rrtmg_sw_read_nc.o rrtmg_sw_read_nc.mod: rrtmg_sw_read_nc.f90 rrsw_kg16.mod rrsw_ncpar.mod rrsw_kg17.mod rrsw_ncpar.mod rrsw_kg18.mod rrsw_ncpar.mod rrsw_kg19.mod rrsw_ncpar.mod rrsw_kg20.mod rrsw_ncpar.mod rrsw_kg21.mod rrsw_ncpar.mod rrsw_kg22.mod rrsw_ncpar.mod rrsw_kg23.mod rrsw_ncpar.mod rrsw_kg24.mod rrsw_ncpar.mod rrsw_kg25.mod rrsw_ncpar.mod rrsw_kg26.mod rrsw_ncpar.mod rrsw_kg27.mod rrsw_ncpar.mod rrsw_kg28.mod rrsw_ncpar.mod rrsw_kg29.mod rrsw_ncpar.mod
rrtmg_sw_reftra.o rrtmg_sw_reftra.mod: rrtmg_sw_reftra.f90 parkind.mod rrsw_tbl.mod rrsw_vsn.mod
rrtmg_sw_setcoef.o rrtmg_sw_setcoef.mod: rrtmg_sw_setcoef.f90 parkind.mod parrrsw.mod rrsw_ref.mod rrsw_vsn.mod
rrtmg_sw_spcvmc.o rrtmg_sw_spcvmc.mod: rrtmg_sw_spcvmc.f90 parkind.mod parrrsw.mod rrsw_tbl.mod rrsw_vsn.mod rrsw_wvn.mod rrtmg_sw_reftra.mod rrtmg_sw_taumol.mod rrtmg_sw_vrtqdr.mod
rrtmg_sw_spcvrt.o rrtmg_sw_spcvrt.mod: rrtmg_sw_spcvrt.f90 parkind.mod parrrsw.mod rrsw_tbl.mod rrsw_vsn.mod rrsw_wvn.mod rrtmg_sw_reftra.mod rrtmg_sw_taumol.mod rrtmg_sw_vrtqdr.mod
rrtmg_sw_taumol.o rrtmg_sw_taumol.mod: rrtmg_sw_taumol.f90 parkind.mod rrsw_con.mod rrsw_wvn.mod rrsw_vsn.mod parrrsw.mod rrsw_kg16.mod parrrsw.mod rrsw_kg17.mod parrrsw.mod rrsw_kg18.mod parrrsw.mod rrsw_kg19.mod parrrsw.mod rrsw_kg20.mod parrrsw.mod rrsw_kg21.mod parrrsw.mod rrsw_kg22.mod parrrsw.mod rrsw_kg23.mod parrrsw.mod rrsw_kg24.mod parrrsw.mod rrsw_kg25.mod parrrsw.mod rrsw_kg26.mod parrrsw.mod rrsw_kg27.mod parrrsw.mod rrsw_kg28.mod parrrsw.mod rrsw_kg29.mod
rrtmg_sw_vrtqdr.o rrtmg_sw_vrtqdr.mod: rrtmg_sw_vrtqdr.f90 parkind.mod
