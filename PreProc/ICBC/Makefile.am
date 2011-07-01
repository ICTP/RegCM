
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = aerosol sst icbc chem_icbc

aerosol_SOURCES = mod_header.F90 aerosol.F90
sst_SOURCES = mod_sst_eh5om.f90 mod_sst_ersst.f90 mod_sst_fvgcm.f90 \
              mod_sst_1deg.f90 mod_sst_ccsm.f90 mod_sst_cam2.f90 mod_interp.f90 \
              mod_sst_grid.F90 mod_header.F90 sst.f90
icbc_SOURCES = mod_ccsm.f90 mod_vectutil.f90 mod_ecwcp.f90      \
               mod_eh5om.f90 mod_ein.f90 mod_era40.f90 mod_erahi.f90 \
               mod_fvgcm.f90 mod_gfs11.f90 mod_grid.f90         \
               mod_ncep.f90 mod_nest.f90 mod_write.F90          \
               mod_interp.f90 mod_vertint.f90 mod_header.F90    \
               mod_hgt.f90 mod_humid.f90 mod_mksst.f90          \
               mod_uvrot.f90 mod_cam2.f90 icbc.f90
chem_icbc_SOURCES = mod_ch_icbc.f90 mod_ch_oxcl.f90 mod_wrtoxd.F90 mod_grid.f90 \
                    mod_header.F90 mod_interp.f90 chem_icbc.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod

distclean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

aerosol.o: aerosol.F90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod
chem_icbc.o: chem_icbc.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_date.mod mod_grid.o mod_wrtoxd.mod mod_header.mod \
 mod_ch_icbc.o mod_ch_oxcl.o $(REGCM_ROOT)/Share/mod_memutil.mod
icbc.o: icbc.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_mksst.o \
 $(REGCM_ROOT)/Share/mod_date.mod mod_grid.o $(REGCM_ROOT)/Share/mod_date.mod \
 mod_ecwcp.mod mod_eh5om.mod mod_ein.mod mod_era40.mod mod_erahi.mod \
 mod_fvgcm.mod mod_gfs11.mod mod_ncep.mod mod_nest.mod mod_ccsm.mod \
 mod_cam2.mod mod_write.o mod_header.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod
mod_cam2.o mod_cam2.mod: mod_cam2.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_zeit.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/Share/mod_constants.mod mod_grid.o \
 mod_write.o mod_interp.mod mod_vertint.o mod_hgt.o mod_humid.o \
 mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_ccsm.o mod_ccsm.mod: mod_ccsm.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_mall.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_constants.mod mod_grid.o \
 mod_write.o mod_interp.mod mod_vertint.o mod_hgt.o mod_humid.o \
 mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_ch_icbc.o mod_ch_icbc.mod: mod_ch_icbc.f90 \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_grid.o mod_wrtoxd.o mod_interp.o \
 $(REGCM_ROOT)/Share/mod_date.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/Share/mod_memutil.mod
mod_ch_oxcl.o mod_ch_oxcl.mod: mod_ch_oxcl.f90 \
 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_grid.o \
 mod_wrtoxd.o mod_interp.o $(REGCM_ROOT)/Share/mod_date.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_realkinds.mod
mod_ecwcp.o mod_ecwcp.mod: mod_ecwcp.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_zeit.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o
mod_eh5om.o mod_eh5om.mod: mod_eh5om.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod \
 mod_grid.o $(REGCM_ROOT)/Share/mod_memutil.mod \
 mod_write.o mod_interp.mod mod_vertint.o mod_hgt.o mod_humid.o \
 mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_ein.o mod_ein.mod: mod_ein.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o 
mod_era40.o mod_era40.mod: mod_era40.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o 
mod_erahi.o mod_erahi.mod: mod_erahi.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_zeit.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o
mod_fvgcm.o mod_fvgcm.mod: mod_fvgcm.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod mod_grid.o \
 mod_write.o mod_interp.mod mod_vertint.o mod_hgt.o mod_humid.o \
 mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_gfs11.o mod_gfs11.mod: mod_gfs11.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod mod_grid.o \
 mod_write.o mod_interp.mod mod_vertint.o mod_hgt.o mod_humid.o \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_mksst.o mod_uvrot.o mod_vectutil.o
mod_grid.o mod_grid.mod: mod_grid.f90 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod
mod_header.o mod_header.mod: mod_header.F90 $(REGCM_ROOT)/external/m_stdio.mod
mod_hgt.o mod_hgt.mod: mod_hgt.f90 $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod
mod_humid.o mod_humid.mod: mod_humid.f90 $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod
mod_interp.o mod_interp.mod: mod_interp.f90 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/Share/mod_constants.mod
mod_mksst.o mod_mksst.mod: mod_mksst.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod
mod_ncep.o mod_ncep.mod: mod_ncep.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_mall.mod \
 $(REGCM_ROOT)/external/m_zeit.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o mod_grid.o 
mod_nest.o mod_nest.mod: mod_nest.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/Share/mod_date.mod mod_grid.o mod_write.o mod_interp.mod \
 mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
 mod_vectutil.o
mod_sst_1deg.o mod_sst_1deg.mod: mod_sst_1deg.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_sst_grid.o \
 mod_interp.o
mod_sst_cam2.o mod_sst_cam2.mod: mod_sst_cam2.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_mall.mod \
 $(REGCM_ROOT)/external/m_zeit.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_sst_grid.o mod_interp.o
mod_sst_ccsm.o mod_sst_ccsm.mod: mod_sst_ccsm.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
 mod_sst_grid.o mod_interp.o 
mod_sst_eh5om.o mod_sst_eh5om.mod: mod_sst_eh5om.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_sst_grid.o mod_interp.o
mod_sst_ersst.o mod_sst_ersst.mod: mod_sst_ersst.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_sst_grid.o mod_interp.o 
mod_sst_fvgcm.o mod_sst_fvgcm.mod: mod_sst_fvgcm.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_sst_grid.o mod_interp.o
mod_sst_grid.o mod_sst_grid.mod: mod_sst_grid.F90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_zeit.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod 
mod_uvrot.o mod_uvrot.mod: mod_uvrot.f90 $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod
mod_vectutil.o mod_vectutil.mod: mod_vectutil.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod
mod_vertint.o mod_vertint.mod: mod_vertint.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/Share/mod_constants.mod
mod_write.o mod_write.mod: mod_write.F90 $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_zeit.mod $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_grid.o $(REGCM_ROOT)/Share/mod_memutil.mod
mod_wrtoxd.o mod_wrtoxd.mod: mod_wrtoxd.F90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_date.mod mod_grid.o $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/Share/mod_memutil.mod
sst.o: sst.f90 mod_header.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
 mod_sst_grid.mod mod_sst_1deg.mod mod_sst_eh5om.mod mod_sst_ersst.mod \
 mod_sst_fvgcm.mod mod_sst_ccsm.mod mod_sst_cam2.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/external/m_die.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_zeit.mod
