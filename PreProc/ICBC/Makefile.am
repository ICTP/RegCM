
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

bin_PROGRAMS = sst icbc chem_icbc

sst_SOURCES = mod_sst_eh5om.F90 mod_sst_ersst.F90 mod_sst_fvgcm.F90 \
              mod_sst_1deg.F90 mod_sst_gnmnc.F90 mod_sst_grid.F90 \
              mod_sst_gnhnc.F90 mod_header.F90 sst.F90 mod_ccsm3_helper.F90
icbc_SOURCES = mod_vectutil.F90 mod_ecwcp.F90 mod_eh5om.F90 mod_ein.F90 \
               mod_era40.F90 mod_erahi.F90 mod_fvgcm.F90 mod_grid.F90  \
               mod_ncep.F90 mod_nest.F90 mod_write.F90 mod_header.F90 \
               mod_mksst.F90 mod_uvrot.F90 mod_ecens.F90 mod_gn6hnc.F90 \
               mod_ccsm3_helper.F90 icbc.F90
chem_icbc_SOURCES = mod_ch_param.F90 mod_ch_icbc.F90 mod_ch_icbc_clim.F90 \
                    mod_ox_icbc.F90 mod_ae_icbc.F90 mod_wrtoxd.F90 \
                    mod_grid.F90 mod_header.F90 chem_icbc.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
chem_icbc.o : chem_icbc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o mod_date.o mod_dynparam.o mod_space.o mod_memutil.o mod_nchelper.o mod_domain.o mod_grid.o mod_ncstream_types.o mod_ncstream.o mod_wrtoxd.o mod_header.o mod_interp.o mod_ch_icbc.o mod_ch_icbc_clim.o mod_ox_icbc.o mod_ae_icbc.o
icbc.o : icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_header.o mod_stdio.o mod_memutil.o mod_mksst.o mod_date.o mod_grid.o mod_date.o mod_ecwcp.o mod_eh5om.o mod_ein.o mod_era40.o mod_erahi.o mod_ecens.o mod_fvgcm.o mod_ncep.o mod_nest.o mod_gn6hnc.o mod_write.o
mod_ch_icbc.o : mod_ch_icbc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_grid.o mod_wrtoxd.o mod_interp.o mod_date.o mod_memutil.o mod_message.o mod_nchelper.o  mod_ch_param.o
mod_ch_icbc_clim.o : mod_ch_icbc_clim.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_grid.o mod_wrtoxd.o mod_interp.o mod_date.o mod_memutil.o mod_message.o mod_nchelper.o  mod_ch_param.o
mod_ox_icbc.o : mod_ox_icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_grid.o mod_wrtoxd.o mod_interp.o mod_date.o mod_nchelper.o
mod_ae_icbc.o : mod_ae_icbc.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_grid.o mod_wrtoxd.o mod_interp.o mod_date.o mod_nchelper.o
mod_ecens.o : mod_ecens.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_nchelper.o
mod_ecwcp.o : mod_ecwcp.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_stdio.o mod_message.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_eh5om.o : mod_eh5om.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_grid.o mod_memutil.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o
mod_ein.o : mod_ein.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_nchelper.o
mod_era40.o : mod_era40.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_memutil.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_nchelper.o
mod_erahi.o : mod_erahi.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_message.o mod_stdio.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_fvgcm.o : mod_fvgcm.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_message.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_gn6hnc.o : mod_gn6hnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_constants.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_nchelper.o mod_ccsm3_helper.o
mod_grid.o : mod_grid.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_stdio.o mod_message.o mod_nchelper.o mod_domain.o mod_dynparam.o mod_nhinterp.o
mod_header.o : mod_header.F90 mod_intkinds.o mod_stdio.o
mod_mksst.o : mod_mksst.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_grid.o mod_memutil.o mod_constants.o mod_dynparam.o mod_message.o mod_nchelper.o
mod_ncep.o : mod_ncep.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_memutil.o mod_nchelper.o mod_grid.o
mod_nest.o : mod_nest.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_constants.o mod_date.o mod_grid.o mod_write.o mod_interp.o mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o mod_vectutil.o mod_message.o mod_memutil.o mod_nchelper.o
mod_sst_1deg.o : mod_sst_1deg.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_sst_grid.o mod_interp.o mod_nchelper.o mod_message.o
mod_sst_eh5om.o : mod_sst_eh5om.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o mod_dynparam.o mod_sst_grid.o mod_interp.o
mod_sst_ersst.o : mod_sst_ersst.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_sst_grid.o mod_interp.o mod_message.o mod_memutil.o mod_nchelper.o
mod_sst_fvgcm.o : mod_sst_fvgcm.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o mod_dynparam.o mod_interp.o mod_sst_grid.o
mod_sst_gnmnc.o : mod_sst_gnmnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_interp.o mod_message.o mod_nchelper.o
mod_sst_gnhnc.o : mod_sst_gnhnc.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_sst_grid.o mod_interp.o mod_message.o mod_nchelper.o mod_ccsm3_helper.o
mod_sst_grid.o : mod_sst_grid.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_message.o mod_nchelper.o mod_ncstream_types.o mod_ncstream.o mod_domain.o
mod_uvrot.o : mod_uvrot.F90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_vectutil.o : mod_vectutil.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_message.o
mod_write.o : mod_write.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_grid.o mod_memutil.o mod_message.o mod_ncstream_types.o mod_ncstream.o mod_nhinterp.o mod_vectutil.o 
mod_wrtoxd.o : mod_wrtoxd.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_date.o mod_grid.o mod_memutil.o mod_message.o mod_stdio.o mod_nchelper.o mod_ncstream.o
sst.o : sst.F90 mod_intkinds.o mod_realkinds.o mod_header.o mod_dynparam.o mod_memutil.o mod_stdio.o mod_message.o mod_sst_grid.o mod_sst_1deg.o mod_sst_eh5om.o mod_sst_ersst.o mod_sst_fvgcm.o mod_sst_gnmnc.o mod_sst_gnhnc.o
mod_ch_param.o : mod_ch_param.F90 mod_intkinds.o mod_realkinds.o
