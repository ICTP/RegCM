
vpath %.o $(REGCM_ROOT)/Share:cmip6:.
vpath %.inc $(REGCM_ROOT)/Share:cmip6:.

include $(REGCM_ROOT)/makeinc

SUBDIRS = cmip6

if DO_COMPILE_LIBCLM45
bin_PROGRAMS = sst icbc chem_icbc clmbc pgwbc pgw_icbc
else
bin_PROGRAMS = sst icbc chem_icbc pgwbc pgw_icbc
endif

CMIP6PTH = cmip6
CMIP6LIB = cmip6/libcmip6.a

if COMPILER_SUN
  CPPFLAGS += -M$(CMIP6PTH)
else
  CPPFLAGS += -I$(CMIP6PTH)
endif

LCMIP6LIB = -L$(CMIP6PTH) -lcmip6

FCLD = $(MPIFC)

LIBS = $(LCMIP6LIB) $(SYSLIBS)

sst_SOURCES = mod_sst_ersst.F90 \
              mod_sst_1deg.F90 mod_sst_gnmnc.F90 mod_sst_grid.F90 \
              mod_sst_gnhnc.F90 mod_header.F90 mod_ccsm3_helper.F90 \
              mod_hadgem_helper.F90 mod_csiro_helper.F90 mod_canesm_helper.F90 \
              mod_miroc_helper.F90 mod_ipsl_helper.F90 mod_gfdl_helper.F90 \
              mod_cnrm_helper.F90 mod_mpiesm_helper.F90 mod_ecearth_helper.F90 \
              mod_noresm_helper.F90 mod_sst_gndnc.F90 mod_ccsm4_helper.F90 \
              mod_lgm_helper.F90 mod_sst_cmip6.F90 sst.F90
icbc_SOURCES = mod_vectutil.F90 mod_ein.F90 mod_grid.F90  \
               mod_ncep.F90 mod_nest.F90 mod_write.F90 mod_header.F90 \
               mod_mksst.F90 mod_gn6hnc.F90 mod_ifs.F90 \
               mod_ccsm3_helper.F90 mod_hadgem_helper.F90 mod_csiro_helper.F90 \
               mod_canesm_helper.F90 mod_miroc_helper.F90 mod_ipsl_helper.F90 \
               mod_gfdl_helper.F90 mod_cnrm_helper.F90 mod_mpiesm_helper.F90 \
               mod_ecearth_helper.F90 mod_noresm_helper.F90 mod_era5.F90 \
               mod_lgm_helper.F90 mod_ccsm4_helper.F90 mod_cmip6.F90 icbc.F90
clmbc_SOURCES = mod_vectutil.F90 mod_grid.F90 mod_write.F90 mod_header.F90 \
               mod_mksst.F90 mod_era5.F90 clmbc.F90
chem_icbc_SOURCES = mod_ch_param.F90 mod_ch_icbc.F90 mod_ch_icbc_clim.F90 \
                    mod_ox_icbc.F90 mod_ae_icbc.F90 mod_wrtoxd.F90 \
                    mod_ch_fnest.F90 mod_grid.F90 mod_header.F90 chem_icbc.F90
pgwbc_SOURCES = mod_vectutil.F90 mod_grid.F90 mod_write.F90 mod_header.F90 \
                mod_pgw.F90 pgwbc.F90
pgw_icbc_SOURCES = mod_header.F90 pgw_icbc.F90 mod_vectutil.F90

LDADD = $(OMPFLAGS)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
chem_icbc.o : chem_icbc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o mod_date.o mod_dynparam.o mod_space.o mod_memutil.o mod_nchelper.o mod_domain.o mod_grid.o mod_ncstream_types.o mod_ncstream.o mod_wrtoxd.o mod_header.o mod_ch_icbc.o mod_ch_icbc_clim.o mod_ox_icbc.o mod_ae_icbc.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
pgw_icbc.o : pgw_icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_header.o mod_stdio.o mod_memutil.o mod_date.o mod_vertint.o mod_date.o mod_vectutil.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
pgwbc.o : pgwbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_header.o mod_stdio.o mod_memutil.o mod_date.o mod_grid.o mod_date.o mod_pgw.o mod_write.o mod_projections.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
icbc.o : icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_header.o mod_stdio.o mod_memutil.o mod_mksst.o mod_date.o mod_grid.o mod_date.o mod_ein.o mod_ifs.o mod_era5.o mod_ncep.o mod_nest.o mod_gn6hnc.o mod_write.o mod_projections.o $(CMIP6LIB)
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
mod_ch_icbc.o : mod_ch_icbc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_grid.o mod_wrtoxd.o mod_kdinterp.o mod_date.o mod_memutil.o mod_message.o mod_nchelper.o mod_ch_param.o mod_constants.o
mod_ch_icbc_clim.o : mod_ch_icbc_clim.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_grid.o mod_wrtoxd.o mod_kdinterp.o mod_date.o mod_memutil.o mod_message.o mod_nchelper.o mod_ch_param.o mod_constants.o
mod_ox_icbc.o : mod_ox_icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_grid.o mod_wrtoxd.o mod_kdinterp.o mod_date.o mod_nchelper.o mod_constants.o mod_stdio.o
mod_ae_icbc.o : mod_ae_icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_grid.o mod_wrtoxd.o mod_kdinterp.o mod_date.o mod_nchelper.o mod_stdio.o mod_constants.o
mod_sst_cmip6.o : mod_sst_cmip6.F90 mod_cmip6_helper.o mod_intkinds.o mod_realkinds.o mod_date.o mod_message.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_kdinterp.o mod_date.o mod_stdio.o
mod_cmip6.o : mod_cmip6.F90 mod_cmip6_helper.o mod_cmip6_cnrm.o mod_cmip6_cesm.o mod_cmip6_normm.o mod_cmip6_miroc6.o mod_cmip6_hadmm.o mod_cmip6_mpihr.o mod_cmip6_ecea.o mod_intkinds.o mod_realkinds.o mod_date.o mod_message.o mod_dynparam.o mod_memutil.o mod_grid.o mod_kdinterp.o mod_write.o mod_vectutil.o mod_mksst.o mod_humid.o mod_hgt.o mod_vertint.o mod_date.o mod_stdio.o
mod_ccsm_helper.o : mod_ccsm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o
mod_ccsm3_helper.o : mod_ccsm3_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_stdio.o
mod_hadgem_helper.o : mod_hadgem_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_csiro_helper.o : mod_csiro_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_canesm_helper.o : mod_canesm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_ecearth_helper.o : mod_ecearth_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_miroc_helper.o : mod_miroc_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_ipsl_helper.o : mod_ipsl_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_gfdl_helper.o : mod_gfdl_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_cnrm_helper.o : mod_cnrm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_mpiesm_helper.o : mod_mpiesm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o mod_message.o mod_stdio.o
mod_lgm_helper.o : mod_lgm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o mod_message.o
mod_noresm_helper.o : mod_noresm_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_ccsm4_helper.o : mod_ccsm4_helper.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o
mod_ch_fnest.o : mod_ch_fnest.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_grid.o mod_wrtoxd.o mod_kdinterp.o mod_date.o mod_nchelper.o mod_posix.o mod_message.o mod_vertint.o mod_stdio.o mod_constants.o
mod_ein.o : mod_ein.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_projections.o mod_vectutil.o mod_message.o mod_nchelper.o mod_kdinterp.o mod_earth.o mod_date.o mod_constants.o
mod_ifs.o : mod_ifs.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_vertint.o mod_hgt.o mod_humid.o mod_projections.o mod_vectutil.o mod_message.o mod_nchelper.o mod_kdinterp.o mod_earth.o mod_constants.o mod_date.o
mod_era5.o : mod_era5.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_projections.o mod_vectutil.o mod_message.o mod_nchelper.o mod_kdinterp.o mod_earth.o mod_date.o mod_constants.o
mod_pgw.o : mod_pgw.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_vertint.o mod_hgt.o mod_humid.o mod_projections.o mod_vectutil.o mod_message.o mod_nchelper.o mod_kdinterp.o mod_earth.o mod_constants.o
mod_gn6hnc.o : mod_gn6hnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_constants.o mod_grid.o mod_write.o mod_kdinterp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_projections.o mod_vectutil.o mod_message.o mod_nchelper.o mod_ccsm3_helper.o mod_hadgem_helper.o mod_canesm_helper.o mod_miroc_helper.o mod_ipsl_helper.o mod_gfdl_helper.o mod_csiro_helper.o mod_cnrm_helper.o mod_mpiesm_helper.o mod_ecearth_helper.o mod_noresm_helper.o mod_ccsm4_helper.o mod_lgm_helper.o mod_date.o
mod_grid.o : mod_grid.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_memutil.o mod_stdio.o mod_message.o mod_nchelper.o mod_domain.o mod_dynparam.o mod_nhinterp.o mod_projections.o mod_constants.o
mod_header.o : mod_header.F90 mod_intkinds.o mod_stdio.o $(REGCM_ROOT)/version
mod_mksst.o : mod_mksst.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_grid.o mod_memutil.o mod_constants.o mod_dynparam.o mod_message.o mod_nchelper.o mod_earth.o mod_stdio.o mod_date.o
mod_ncep.o : mod_ncep.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_grid.o mod_write.o mod_kdinterp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_projections.o mod_vectutil.o mod_message.o mod_memutil.o mod_nchelper.o mod_grid.o mod_date.o mod_constants.o
mod_nest.o : mod_nest.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_constants.o mod_date.o mod_grid.o mod_write.o mod_kdinterp.o mod_vertint.o mod_hgt.o mod_humid.o mod_projections.o mod_vectutil.o mod_message.o mod_memutil.o mod_nchelper.o mod_stdatm.o
mod_sst_1deg.o : mod_sst_1deg.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_sst_grid.o mod_kdinterp.o mod_nchelper.o mod_message.o mod_date.o
mod_sst_ersst.o : mod_sst_ersst.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_sst_grid.o mod_kdinterp.o mod_message.o mod_memutil.o mod_nchelper.o mod_date.o
mod_sst_gnmnc.o : mod_sst_gnmnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_kdinterp.o mod_message.o mod_nchelper.o mod_hadgem_helper.o mod_canesm_helper.o mod_miroc_helper.o mod_ipsl_helper.o mod_gfdl_helper.o mod_csiro_helper.o mod_cnrm_helper.o mod_ecearth_helper.o mod_ccsm4_helper.o mod_date.o
mod_sst_gnhnc.o : mod_sst_gnhnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_message.o mod_nchelper.o mod_ccsm3_helper.o mod_mpiesm_helper.o mod_ecearth_helper.o mod_earth.o mod_kdinterp.o mod_lgm_helper.o mod_date.o
mod_sst_gndnc.o : mod_sst_gndnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_message.o mod_nchelper.o mod_noresm_helper.o mod_earth.o mod_kdinterp.o mod_date.o
mod_sst_grid.o : mod_sst_grid.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_message.o mod_nchelper.o mod_ncstream_types.o mod_ncstream.o mod_domain.o
mod_vectutil.o : mod_vectutil.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_message.o
mod_write.o : mod_write.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_grid.o mod_memutil.o mod_message.o mod_ncstream_types.o mod_ncstream.o mod_nhinterp.o mod_vectutil.o mod_zita.o mod_stdio.o
mod_wrtoxd.o : mod_wrtoxd.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o mod_grid.o mod_memutil.o mod_message.o mod_stdio.o mod_nchelper.o mod_ncstream.o
sst.o : sst.F90 mod_intkinds.o mod_realkinds.o mod_header.o mod_dynparam.o mod_memutil.o mod_stdio.o mod_message.o mod_sst_grid.o mod_sst_1deg.o mod_sst_ersst.o mod_sst_gnmnc.o mod_sst_gndnc.o mod_sst_gnhnc.o mod_cmip6.o mod_date.o $(CMIP6LIB)
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
mod_ch_param.o : mod_ch_param.F90 mod_intkinds.o mod_realkinds.o
clmbc.o : clmbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_header.o mod_stdio.o mod_memutil.o mod_date.o mod_grid.o mod_date.o mod_write.o mod_era5.o mod_projections.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
