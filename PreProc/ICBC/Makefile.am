
include $(REGCM_ROOT)/makeinc

vpath %.mod $(REGCM_ROOT)/Share:.

bin_PROGRAMS = sst icbc chem_icbc

sst_SOURCES = mod_sst_eh5om.f90 mod_sst_ersst.F90 mod_sst_fvgcm.f90 \
              mod_sst_1deg.F90 mod_sst_gnmnc.F90 mod_sst_grid.F90 \
              mod_header.F90 sst.f90
icbc_SOURCES = mod_vectutil.f90 mod_ecwcp.f90 mod_eh5om.F90 mod_ein.F90 \
               mod_era40.F90 mod_erahi.f90 mod_fvgcm.f90 mod_grid.F90  \
               mod_ncep.F90 mod_nest.F90 mod_write.F90 mod_header.F90 \
               mod_hgt.F90 mod_humid.f90 mod_mksst.F90 mod_uvrot.f90 \
               mod_ecens.F90 mod_gn6hnc.F90 icbc.f90
chem_icbc_SOURCES = mod_ch_icbc.F90 mod_ch_oxcl.F90 mod_wrtoxd.F90 \
                    mod_grid.F90 mod_header.F90 chem_icbc.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
chem_icbc.o chem_icbc.mod: chem_icbc.f90 mod_dynparam.mod mod_date.mod mod_grid.mod mod_wrtoxd.mod mod_header.mod mod_ch_icbc.mod mod_ch_oxcl.mod mod_memutil.mod
icbc.o icbc.mod: icbc.f90 mod_dynparam.mod mod_message.mod mod_header.mod mod_stdio.mod mod_memutil.mod mod_mksst.mod mod_date.mod mod_grid.mod mod_date.mod mod_ecwcp.mod mod_eh5om.mod mod_ein.mod mod_era40.mod mod_erahi.mod mod_ecens.mod mod_fvgcm.mod mod_ncep.mod mod_nest.mod mod_gn6hnc.mod mod_write.mod
mod_ch_icbc.o mod_ch_icbc.mod: mod_ch_icbc.F90 mod_dynparam.mod mod_grid.mod mod_wrtoxd.mod mod_interp.mod mod_date.mod mod_realkinds.mod mod_memutil.mod mod_message.mod mod_nchelper.mod
mod_ch_oxcl.o mod_ch_oxcl.mod: mod_ch_oxcl.F90 mod_dynparam.mod mod_memutil.mod mod_grid.mod mod_wrtoxd.mod mod_interp.mod mod_date.mod mod_realkinds.mod mod_nchelper.mod
mod_ecens.o mod_ecens.mod: mod_ecens.F90 mod_dynparam.mod mod_realkinds.mod mod_stdio.mod mod_memutil.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_nchelper.mod
mod_ecwcp.o mod_ecwcp.mod: mod_ecwcp.f90 mod_dynparam.mod mod_memutil.mod mod_stdio.mod mod_realkinds.mod mod_message.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod
mod_eh5om.o mod_eh5om.mod: mod_eh5om.F90 mod_dynparam.mod mod_realkinds.mod mod_stdio.mod mod_grid.mod mod_memutil.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod
mod_ein.o mod_ein.mod: mod_ein.F90 mod_dynparam.mod mod_realkinds.mod mod_stdio.mod mod_memutil.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_nchelper.mod
mod_era40.o mod_era40.mod: mod_era40.F90 mod_dynparam.mod mod_realkinds.mod mod_stdio.mod mod_memutil.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_nchelper.mod
mod_erahi.o mod_erahi.mod: mod_erahi.f90 mod_dynparam.mod mod_memutil.mod mod_realkinds.mod mod_message.mod mod_stdio.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod
mod_fvgcm.o mod_fvgcm.mod: mod_fvgcm.f90 mod_dynparam.mod mod_memutil.mod mod_realkinds.mod mod_message.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod
mod_gn6hnc.o mod_gn6hnc.mod: mod_gn6hnc.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_memutil.mod mod_constants.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_nchelper.mod
mod_grid.o mod_grid.mod: mod_grid.F90 mod_memutil.mod mod_realkinds.mod mod_stdio.mod mod_message.mod mod_nchelper.mod mod_domain.mod mod_dynparam.mod
mod_header.o mod_header.mod: mod_header.F90 mod_stdio.mod
mod_hgt.o mod_hgt.mod: mod_hgt.F90 mod_constants.mod mod_message.mod mod_realkinds.mod
mod_humid.o mod_humid.mod: mod_humid.f90 mod_constants.mod mod_realkinds.mod
mod_mksst.o mod_mksst.mod: mod_mksst.F90 mod_realkinds.mod mod_stdio.mod mod_grid.mod mod_memutil.mod mod_constants.mod mod_dynparam.mod mod_message.mod mod_nchelper.mod
mod_ncep.o mod_ncep.mod: mod_ncep.F90 mod_dynparam.mod mod_realkinds.mod mod_stdio.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_memutil.mod mod_nchelper.mod mod_grid.mod
mod_nest.o mod_nest.mod: mod_nest.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_constants.mod mod_date.mod mod_grid.mod mod_write.mod mod_interp.mod mod_vertint.mod mod_hgt.mod mod_humid.mod mod_mksst.mod mod_uvrot.mod mod_vectutil.mod mod_message.mod mod_memutil.mod mod_nchelper.mod
mod_sst_1deg.o mod_sst_1deg.mod: mod_sst_1deg.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_sst_grid.mod mod_interp.mod mod_nchelper.mod mod_message.mod
mod_sst_eh5om.o mod_sst_eh5om.mod: mod_sst_eh5om.f90 mod_realkinds.mod mod_stdio.mod mod_message.mod mod_dynparam.mod mod_sst_grid.mod mod_interp.mod
mod_sst_ersst.o mod_sst_ersst.mod: mod_sst_ersst.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_sst_grid.mod mod_interp.mod mod_message.mod mod_nchelper.mod
mod_sst_fvgcm.o mod_sst_fvgcm.mod: mod_sst_fvgcm.f90 mod_realkinds.mod mod_stdio.mod mod_message.mod mod_dynparam.mod mod_interp.mod mod_sst_grid.mod
mod_sst_gnmnc.o mod_sst_gnmnc.mod: mod_sst_gnmnc.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_memutil.mod mod_sst_grid.mod mod_interp.mod mod_message.mod mod_nchelper.mod
mod_sst_grid.o mod_sst_grid.mod: mod_sst_grid.F90 mod_realkinds.mod mod_stdio.mod mod_dynparam.mod mod_memutil.mod mod_message.mod mod_nchelper.mod mod_domain.mod
mod_uvrot.o mod_uvrot.mod: mod_uvrot.f90 mod_constants.mod mod_realkinds.mod
mod_vectutil.o mod_vectutil.mod: mod_vectutil.f90 mod_realkinds.mod mod_stdio.mod
mod_write.o mod_write.mod: mod_write.F90 mod_realkinds.mod mod_dynparam.mod mod_grid.mod mod_memutil.mod mod_message.mod mod_nchelper.mod
mod_wrtoxd.o mod_wrtoxd.mod: mod_wrtoxd.F90 mod_dynparam.mod mod_date.mod mod_grid.mod mod_realkinds.mod mod_memutil.mod mod_message.mod mod_stdio.mod mod_nchelper.mod
sst.o sst.mod: sst.f90 mod_header.mod mod_dynparam.mod mod_memutil.mod mod_stdio.mod mod_message.mod mod_sst_grid.mod mod_sst_1deg.mod mod_sst_eh5om.mod mod_sst_ersst.mod mod_sst_fvgcm.mod mod_sst_gnmnc.mod
