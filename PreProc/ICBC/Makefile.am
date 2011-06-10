
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:.
vpath %.f90 $(REGCM_ROOT)/Config:.

bin_PROGRAMS = aerosol sst icbc oxidant

aerosol_SOURCES = mod_header.F90 aerosol.F90
sst_SOURCES = mod_sst_eh5om.f90 mod_sst_ersst.f90 mod_sst_fvgcm.f90 \
							mod_sst_1deg.f90 mod_sst_ccsm.f90 mod_interp.f90      \
							mod_sst_grid.F90 mod_header.F90 sst.f90
icbc_SOURCES = mod_ccsm.f90 mod_vectutil.f90 mod_ecwcp.f90      \
							 mod_eh5om.f90 mod_ein.f90 mod_era40.f90 mod_erahi.f90 \
               mod_fvgcm.f90 mod_gfs11.f90 mod_grid.f90         \
							 mod_ncep.f90 mod_nest.f90 mod_write.F90          \
							 mod_interp.f90 mod_vertint.f90 mod_header.F90    \
							 mod_hgt.f90 mod_humid.f90 mod_mksst.f90          \
							 mod_uvrot.f90 icbc.f90 mod_cam2.f90
oxidant_SOURCES = mod_oxidant.f90 mod_wrtoxd.F90 mod_grid.f90      \
									mod_header.F90 mod_interp.f90 oxidant.f90

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

mod_ccsm.o:  mod_ccsm.f90  mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_cam2.o:  mod_cam2.f90  mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ecwcp.o: mod_ecwcp.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_eh5om.o: mod_eh5om.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ein.o: mod_ein.f90 mod_grid.o mod_write.o mod_interp.o \
           mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
           mod_vectutil.o
mod_era40.o: mod_era40.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_erahi.o: mod_erahi.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_fvgcm.o: mod_fvgcm.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_gfs11.o: mod_gfs11.f90 mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ncep.o:  mod_ncep.f90  mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_grid.o: mod_grid.f90
mod_hgt.o: mod_hgt.f90
mod_humid.o: mod_humid.f90
mod_interp.o: mod_interp.f90
mod_mksst.o: mod_mksst.f90
mod_nest.o: mod_nest.f90 mod_grid.o mod_write.o mod_interp.o
mod_uvrot.o: mod_uvrot.f90
mod_vectutil.o: mod_vectutil.f90
mod_vertint.o: mod_vertint.f90
mod_write.o: mod_write.F90 mod_grid.o
mod_wrtoxd.o : mod_wrtoxd.F90 mod_grid.o
mod_oxidant.o : mod_oxidant.f90 mod_wrtoxd.o \
	              mod_grid.o mod_interp.o
oxidant.o : oxidant.f90 mod_grid.o mod_wrtoxd.o \
            mod_header.o mod_oxidant.o
aerosol.o: aerosol.F90
mod_sst_grid.o: mod_sst_grid.F90
mod_sst_1deg.o: mod_sst_1deg.f90 mod_interp.o mod_sst_grid.o
mod_sst_ccsm.o: mod_sst_ccsm.f90 mod_interp.o mod_sst_grid.o
mod_sst_eh5om.o: mod_sst_eh5om.f90 mod_interp.o \
	               mod_sst_grid.o
mod_sst_ersst.o: mod_sst_ersst.f90 mod_interp.o \
	               mod_sst_grid.o
mod_sst_fvgcm.o: mod_sst_fvgcm.f90 mod_interp.o
sst.o: sst.f90 mod_sst_1deg.o mod_sst_eh5om.o mod_sst_ersst.o \
       mod_sst_fvgcm.o mod_sst_ccsm.o mod_sst_grid.o mod_header.o
icbc.o: icbc.f90 mod_grid.o mod_write.o \
        mod_ein.o mod_ccsm.o mod_cam2.o mod_ncep.o mod_interp.o mod_header.o \
				mod_mksst.o
