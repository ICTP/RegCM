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

include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/external:$(REGCM_ROOT)/Main/mpplib:RRTMG/:.

SUBDIRS = RRTMG
RRTMG = RRTMG/librrtmg.a

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) -M$(REGCM_ROOT)/external
NOOPTFC += -M$(MPP_INCLUDE) -M$(REGCM_ROOT)/external
else
FCFLAGS += -I$(MPP_INCLUDE) -I$(REGCM_ROOT)/external
NOOPTFC += -I$(MPP_INCLUDE) -I$(REGCM_ROOT)/external
endif

noinst_LIBRARIES = librad.a
librad_a_SOURCES = mod_rad_colmod3.F90 mod_rad_o3blk.F90 mod_rad_outrad.F90 \
  mod_rad_common.F90 mod_rad_radiation.F90 mod_rad_tracer.F90 \
  mod_rad_aerosol.F90 mod_atm_msis.F90 mod_rrtmg_driver.F90 \
  rrtmg_sw_cldprmc.F90 rrtmg_sw_k_g.F90 rrtmg_sw_spcvmc.F90 \
  rrtmg_sw_vrtqdr.F90 rrtmg_sw_cldprop.F90 rrtmg_sw_rad.F90 \
  rrtmg_sw_reftra.F90 rrtmg_sw_spcvrt.F90 rrtmg_sw_init.F90 \
  rrtmg_sw_rad.nomcica.F90 rrtmg_sw_setcoef.F90 rrtmg_sw_taumol.F90 \
  mcica_random_numbers.F90 mcica_subcol_gen_sw.F90 \
  rrtmg_lw_cldprop.F90 rrtmg_lw_rad.F90 rrtmg_lw_rtrn.F90 \
  rrtmg_lw_setcoef.F90 mcica_subcol_gen_lw.F90 rrtmg_lw_init.F90 \
  rrtmg_lw_rad.nomcica.F90 rrtmg_lw_rtrnmc.F90 rrtmg_lw_taumol.F90  \
  rrtmg_lw_cldprmc.F90  rrtmg_lw_k_g.F90 rrtmg_lw_rtrnmr.F90

if COMPILER_SUN
  CPPFLAGS += -MRRTMG
else
  CPPFLAGS += -IRRTMG
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.F90

clean-local:
	rm -f *.mod *__genmod.F90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mcica_random_numbers.o : mcica_random_numbers.F90 parkind.o
mcica_subcol_gen_lw.o : mcica_subcol_gen_lw.F90 parkind.o parrrtm.o rrlw_con.o rrlw_wvn.o rrlw_vsn.o mcica_random_numbers.o
mcica_subcol_gen_sw.o : mcica_subcol_gen_sw.F90 parkind.o parrrsw.o rrsw_con.o rrsw_wvn.o rrsw_vsn.o mcica_random_numbers.o
mod_atm_msis.o : mod_atm_msis.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_spline.o
mod_rad_aerosol.o : mod_rad_aerosol.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_regcm_types.o mod_runparams.o mod_constants.o mod_memutil.o mod_mpmessage.o mod_rad_common.o mod_mppparam.o mod_nhinterp.o mod_zita.o mo_simple_plumes.o mod_date.o mod_stdio.o
mod_rad_colmod3.o : mod_rad_colmod3.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_service.o mod_runparams.o mod_rad_radiation.o mod_rad_common.o mod_rad_outrad.o mod_rrtmg_driver.o mod_rad_aerosol.o mod_regcm_types.o
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(NOOPTFC) -c $<
mod_rad_common.o : mod_rad_common.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o
mod_rad_o3blk.o : mod_rad_o3blk.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_date.o mod_kdinterp.o mod_vertint.o mod_mppparam.o mod_mpmessage.o mod_dynparam.o mod_memutil.o mod_rad_common.o mod_stdio.o mod_runparams.o mod_regcm_types.o
mod_rad_outrad.o : mod_rad_outrad.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_rad_common.o mod_runparams.o mod_outvars.o mod_regcm_types.o mod_constants.o
mod_rad_radiation.o : mod_rad_radiation.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_service.o mod_runparams.o mod_memutil.o mod_rad_common.o mod_rad_tracer.o mod_ipcc_scenario.o mod_rad_aerosol.o mod_mppparam.o mod_constants.o
mod_rad_tracer.o : mod_rad_tracer.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o
mod_rrtmg_driver.o : mod_rrtmg_driver.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_memutil.o mod_rad_common.o mod_ipcc_scenario.o mod_rad_tracer.o mod_rad_aerosol.o mod_stdatm.o rrtmg_sw_rad.o mcica_subcol_gen_sw.o parrrsw.o rrsw_wvn.o parrrtm.o rrtmg_lw_rad.o mod_rad_outrad.o mod_mpmessage.o mod_runparams.o mod_regcm_types.o rrtmg_lw_rad.nomcica.o rrtmg_sw_rad.nomcica.o mod_mppparam.o
rrtmg_lw_cldprmc.o : rrtmg_lw_cldprmc.F90 parkind.o parrrtm.o rrlw_cld.o rrlw_wvn.o rrlw_vsn.o
rrtmg_lw_cldprop.o : rrtmg_lw_cldprop.F90 parkind.o parrrtm.o rrlw_cld.o rrlw_vsn.o
rrtmg_lw_init.o : rrtmg_lw_init.F90 parkind.o parrrtm.o rrlw_wvn.o rrtmg_lw_setcoef.o rrlw_tbl.o rrlw_vsn.o rrlw_con.o rrlw_vsn.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o rrlw_cld.o rrtmg_lw_k_g.o
rrtmg_lw_rad.o : rrtmg_lw_rad.F90 parkind.o parrrtm.o rrlw_vsn.o mcica_subcol_gen_lw.o rrtmg_lw_cldprmc.o rrtmg_lw_rtrnmc.o rrtmg_lw_setcoef.o rrtmg_lw_taumol.o rrlw_con.o rrlw_wvn.o
rrtmg_lw_rad.nomcica.o : rrtmg_lw_rad.nomcica.F90 parkind.o parrrtm.o rrlw_vsn.o rrtmg_lw_cldprop.o rrtmg_lw_rtrn.o rrtmg_lw_rtrnmr.o rrtmg_lw_setcoef.o rrtmg_lw_taumol.o rrlw_con.o rrlw_wvn.o
rrtmg_lw_k_g.o : rrtmg_lw_k_g.F90 parkind.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o rrlw_vsn.o
	$(FC) $(CPPFLAGS) $(NOOPTFC) -c $<
rrtmg_lw_read_nc.o : rrtmg_lw_read_nc.F90 rrlw_kg01.o rrlw_ncpar.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o
rrtmg_lw_rtrn.o : rrtmg_lw_rtrn.F90 parkind.o parrrtm.o rrlw_con.o rrlw_wvn.o rrlw_tbl.o rrlw_vsn.o
rrtmg_lw_rtrnmc.o : rrtmg_lw_rtrnmc.F90 parkind.o parrrtm.o rrlw_con.o rrlw_wvn.o rrlw_tbl.o rrlw_vsn.o
rrtmg_lw_rtrnmr.o : rrtmg_lw_rtrnmr.F90 parkind.o parrrtm.o rrlw_con.o rrlw_wvn.o rrlw_tbl.o rrlw_vsn.o
rrtmg_lw_setcoef.o : rrtmg_lw_setcoef.F90 parkind.o parrrtm.o rrlw_wvn.o rrlw_ref.o rrlw_vsn.o
rrtmg_lw_taumol.o : rrtmg_lw_taumol.F90 parkind.o parrrtm.o rrlw_con.o rrlw_wvn.o rrlw_vsn.o rrlw_ref.o rrlw_kg01.o rrlw_kg02.o rrlw_kg03.o rrlw_kg04.o rrlw_kg05.o rrlw_kg06.o rrlw_kg07.o rrlw_kg08.o rrlw_kg09.o rrlw_kg10.o rrlw_kg11.o rrlw_kg12.o rrlw_kg13.o rrlw_kg14.o rrlw_kg15.o rrlw_kg16.o
rrtmg_sw_cldprmc.o : rrtmg_sw_cldprmc.F90 parkind.o parrrsw.o rrsw_cld.o rrsw_wvn.o rrsw_vsn.o
rrtmg_sw_cldprop.o : rrtmg_sw_cldprop.F90 parkind.o parrrsw.o rrsw_cld.o rrsw_wvn.o rrsw_vsn.o
rrtmg_sw_init.o : rrtmg_sw_init.F90 parkind.o rrsw_wvn.o rrtmg_sw_setcoef.o parrrsw.o rrsw_tbl.o rrsw_vsn.o rrsw_con.o rrsw_aer.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_cld.o rrtmg_sw_k_g.o
rrtmg_sw_rad.o : rrtmg_sw_rad.F90 parkind.o rrsw_vsn.o mcica_subcol_gen_sw.o rrtmg_sw_cldprmc.o rrtmg_sw_setcoef.o rrtmg_sw_spcvmc.o parrrsw.o rrsw_aer.o rrsw_con.o rrsw_wvn.o
rrtmg_sw_rad.nomcica.o : rrtmg_sw_rad.nomcica.F90 parkind.o rrsw_vsn.o rrtmg_sw_cldprop.o rrtmg_sw_setcoef.o rrtmg_sw_spcvrt.o parrrsw.o rrsw_aer.o rrsw_con.o rrsw_wvn.o
rrtmg_sw_k_g.o : rrtmg_sw_k_g.F90 parkind.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o rrsw_vsn.o
	$(FC) $(CPPFLAGS) $(NOOPTFC) -c $<
rrtmg_sw_read_nc.o : rrtmg_sw_read_nc.F90 rrsw_kg16.o rrsw_ncpar.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o
rrtmg_sw_reftra.o : rrtmg_sw_reftra.F90 parkind.o rrsw_tbl.o rrsw_vsn.o
rrtmg_sw_setcoef.o : rrtmg_sw_setcoef.F90 parkind.o parrrsw.o rrsw_ref.o rrsw_vsn.o
rrtmg_sw_spcvmc.o : rrtmg_sw_spcvmc.F90 parkind.o parrrsw.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o rrtmg_sw_reftra.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o
rrtmg_sw_spcvrt.o : rrtmg_sw_spcvrt.F90 parkind.o parrrsw.o rrsw_tbl.o rrsw_vsn.o rrsw_wvn.o rrtmg_sw_reftra.o rrtmg_sw_taumol.o rrtmg_sw_vrtqdr.o
rrtmg_sw_taumol.o : rrtmg_sw_taumol.F90 parkind.o rrsw_con.o rrsw_wvn.o rrsw_vsn.o parrrsw.o rrsw_kg16.o rrsw_kg17.o rrsw_kg18.o rrsw_kg19.o rrsw_kg20.o rrsw_kg21.o rrsw_kg22.o rrsw_kg23.o rrsw_kg24.o rrsw_kg25.o rrsw_kg26.o rrsw_kg27.o rrsw_kg28.o rrsw_kg29.o
rrtmg_sw_vrtqdr.o : rrtmg_sw_vrtqdr.F90 parkind.o
