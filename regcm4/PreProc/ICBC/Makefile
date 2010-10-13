# Makefile for RegCM4
# for ICBC code.
  
include ../../Makefile.inc

all: aerosol sst icbc oxidant install

OXIDANTOBJ =  mod_oxidant.o oxidant.o mod_wrtoxd.o mod_date.o mod_grid.o \
							mod_header.o mod_interp.o
AEROSOLOBJ = aerosol.o mod_date.o
SSTOBJ = sst.o mod_sst_eh5om.o mod_sst_ersst.o mod_sst_fvgcm.o mod_sst_1deg.o \
         mod_sst_ccsm.o mod_interp.o mod_date.o mod_sst_grid.o mod_header.o

SHRMOD = mod_dynparam.o mod_constants.o

MODOBJ = mod_ccsm.o mod_date.o mod_vectutil.o mod_ecwcp.o mod_eh5om.o \
         mod_ein15.o mod_ein25.o mod_ein75.o mod_era40.o mod_erahi.o \
         mod_fvgcm.o mod_gfs11.o mod_grid.o mod_ncep.o mod_nest.o \
				 mod_write.o mod_interp.o mod_vertint.o \
         mod_header.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o

ICBCOBJ = icbc.o

oxidant : $(SHRMOD) $(OXIDANTOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMOD) $(OXIDANTOBJ) $(NETCDFLIB)

aerosol: $(AEROSOLOBJ) $(SHRMOD)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(AEROSOLOBJ) $(SHRMOD) \
               $(NETCDFLIB)

sst: $(SHRMOD) $(SSTOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SSTOBJ) $(SHRMOD) \
               $(NETCDFLIB)

icbc: $(SHRMOD) $(MODOBJ) $(ICBCOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMOD) $(MODOBJ) $(ICBCOBJ) \
               $(NETCDFLIB)

install: aerosol icbc sst oxidant
	@cp aerosol icbc sst oxidant $(BINDIR)
	@chmod 755 $(BINDIR)/aerosol $(BINDIR)/icbc $(BINDIR)/sst $(BINDIR)/oxidant
	@echo "Programs all installed in $(BINDIR)" 

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_dynparam.o: ../../Config/mod_dynparam.F90 mod_constants.o
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_header.o: mod_header.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) -c $(CPPFLAGS) $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o aerosol sst icbc oxidant

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_date.o: mod_date.f90
mod_ccsm.o:  mod_ccsm.f90  mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ecwcp.o: mod_ecwcp.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_eh5om.o: mod_eh5om.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ein15.o: mod_ein15.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ein25.o: mod_ein25.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ein75.o: mod_ein75.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_era40.o: mod_era40.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_erahi.o: mod_erahi.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_fvgcm.o: mod_fvgcm.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_gfs11.o: mod_gfs11.f90 mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_ncep.o:  mod_ncep.f90  mod_dynparam.o mod_grid.o mod_write.o mod_interp.o \
             mod_vertint.o mod_hgt.o mod_humid.o mod_mksst.o mod_uvrot.o \
             mod_vectutil.o
mod_grid.o: mod_grid.f90 mod_dynparam.o
mod_hgt.o: mod_hgt.f90 mod_constants.o
mod_humid.o: mod_humid.f90 mod_constants.o
mod_interp.o: mod_interp.f90 mod_constants.o
mod_mksst.o: mod_mksst.f90 mod_date.o
mod_nest.o: mod_nest.f90 mod_dynparam.o mod_grid.o mod_write.o \
            mod_constants.o mod_date.o mod_interp.o
mod_uvrot.o: mod_uvrot.f90 mod_constants.o
mod_vectutil.o: mod_vectutil.f90
mod_vertint.o: mod_vertint.f90 mod_constants.o
mod_write.o: mod_write.F90 mod_dynparam.o mod_date.o
mod_wrtoxd.o : mod_wrtoxd.F90 mod_grid.o mod_date.o mod_dynparam.o
mod_oxidant.o : mod_oxidant.f90 mod_dynparam.o mod_date.o mod_wrtoxd.o \
	              mod_grid.o mod_interp.o
oxidant.o : oxidant.f90 mod_dynparam.o mod_date.o mod_grid.o mod_wrtoxd.o \
            mod_header.o mod_oxidant.o
aerosol.o: aerosol.f90 mod_dynparam.o mod_date.o
mod_sst_1deg.o: mod_sst_1deg.f90 mod_dynparam.o mod_interp.o mod_sst_grid.o
mod_sst_ccsm.o: mod_sst_ccsm.f90 mod_dynparam.o mod_interp.o mod_sst_grid.o
mod_sst_eh5om.o: mod_sst_eh5om.f90 mod_dynparam.o mod_date.o mod_interp.o \
	               mod_sst_grid.o
mod_sst_ersst.o: mod_sst_ersst.f90 mod_dynparam.o mod_date.o mod_interp.o \
	               mod_sst_grid.o
mod_sst_fvgcm.o: mod_sst_fvgcm.f90 mod_dynparam.o mod_interp.o
sst.o: mod_dynparam.o mod_sst_1deg.o mod_sst_eh5om.o mod_sst_ersst.o \
       mod_sst_fvgcm.o mod_sst_ccsm.o mod_sst_grid.o mod_header.o
icbc.o: icbc.f90 mod_dynparam.o mod_date.o mod_grid.o mod_write.o \
        mod_ein15.o mod_interp.o mod_header.o mod_mksst.o
