# Makefile for RegCM4
# for ICBC code.
  
include ../../Makefile.inc


all: aerosol sst_1deg sst_ccsm sst_eh50m sst_ersst sst_fvgcm icbc

AEROSOLOBJ = aerosol.o bilinx.o
SST1DEGOBJ = sst_1deg.o bilinx.o printl.o
SSTCCSMOBJ = sst_ccsm.o bilinx.o printl.o
SSTEH50OBJ = sst_ccsm.o bilinx.o printl.o
SSTERAOBJ = sst_ersst.o bilinx.o printl.o mod_datenum.o
SSTFVGCOBJ = sst_fvgcm.o bilinx.o printl.o mod_datenum.o

SHRMOD = mod_regcm_param.o mod_preproc_param.o mod_constants.o

MODOBJ = mod_ccsm.o mod_datenum.o mod_datewk.o mod_ecwcp.o mod_eh5om.o \
         mod_ein15.o mod_ein25.o mod_ein75.o mod_era40.o mod_erahi.o \
         mod_fvgcm.o mod_gfs11.o mod_grid.o mod_ingrid.o mod_mxncom.o \
         mod_ncep.o mod_nest.o mod_write.o
ICBCOBJ = icbc.o bilinx.o cressmcr.o cressmdt.o fexist.o gradsctlb.o \
          gradsctl.o height.o height_o.o htsig.o htsig_o.o humid1.o \
          humid1fv.o humid1_o.o humid2.o humid2fv.o hydrost.o intgtb.o \
          intlin.o intlin_o.o intlog.o intlog_o.o intpsn.o intv1.o \
          intv2.o intv3.o julian.o mksst2.o mksst2a.o mksst3.o mksst.o \
          mkssta.o mxmn3d.o p1p2.o printl.o top2btm.o uvrot4.o uvrot4nx.o \
          header_ICBC.o

aerosol: $(SHRMOD) $(AEROSOLOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(AEROSOLOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sst_1deg: $(SHRMOD) $(SST1DEGOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SST1DEGOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sst_ccsm: $(SHRMOD) $(SSTCCSMOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SSTCCSMOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sst_eh50m: $(SHRMOD) $(SSTEH50OBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SSTEH50OBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sst_ersst: $(SHRMOD) $(SSTERAOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SSTERAOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sst_fvgcm: $(SHRMOD) $(SSTFVGCOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SSTFVGCOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

icbc: $(SHRMOD) $(MODOBJ) $(ICBCOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMOD) $(MODOBJ) $(ICBCOBJ) \
               $(LIB) $(NETCDFLIB)

mod_regcm_param.o: ../../Config/mod_regcm_param.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_preproc_param.o: ../../Config/mod_preproc_param.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

header_ICBC.o: header_ICBC.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o aerosol sst_1deg sst_ccsm \
                        sst_eh50m sst_ersst sst_fvgcm icbc

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_ccsm.o: mod_ccsm.f90 mod_regcm_param.o mod_preproc_param.o \
            mod_grid.o mod_write.o
mod_datenum.o: mod_datenum.f90
mod_datewk.o: mod_datewk.f90
mod_ecwcp.o: mod_ecwcp.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_eh5om.o: mod_eh5om.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_ein15.o: mod_ein15.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_ein25.o: mod_ein25.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_ein75.o: mod_ein75.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_era40.o: mod_era40.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_erahi.o: mod_erahi.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_fvgcm.o: mod_fvgcm.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_gfs11.o: mod_gfs11.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_grid.o: mod_grid.f90 mod_regcm_param.o
mod_ingrid.o: mod_ingrid.f90 mod_regcm_param.o mod_grid.o
mod_mxncom.o: mod_mxncom.f90
mod_ncep.o: mod_ncep.f90 mod_preproc_param.o mod_regcm_param.o \
             mod_grid.o mod_write.o
mod_nest.o: mod_nest.f90 mod_preproc_param.o mod_regcm_param.o \
            mod_grid.o mod_write.o mod_datenum.o mod_mxncom.o
mod_write.o: mod_write.f90 mod_regcm_param.o
bilinx.o: bilinx.f90
cressmcr.o: mod_mxncom.o
cressmdt.o: mod_mxncom.o
fexist.o: fexist.f90
gradsctlb.o: gradsctlb.f90 mod_preproc_param.o mod_regcm_param.o mod_grid.o
gradsctl.o: gradsctl.f90 mod_preproc_param.o mod_regcm_param.o mod_grid.o
height.o: height.f90 mod_constants.o
height_o.o: height_o.f90 mod_constants.o
htsig.o: htsig.f90 mod_constants.o
htsig_o.o: htsig_o.f90 mod_constants.o
humid1.o: humid1.f90 mod_constants.o
humid1fv.o: humid1fv.f90 mod_constants.o
humid1_o.o: humid1_o.f90 mod_constants.o
humid2.o: humid2.f90 mod_constants.o
humid2fv.o: humid2fv.f90 mod_constants.o
hydrost.o: hydrost.f90 mod_constants.o
intgtb.o: intgtb.f90
intlin.o: intlin.f90
intlin_o.o: intlin_o.f90
intlog.o: intlog.f90 mod_constants.o
intlog_o.o: intlog_o.f90 mod_constants.o
intpsn.o: intpsn.f90 mod_constants.o
intv1.o: intv1.f90
intv2.o: intv2.f90 mod_constants.o
intv3.o: intv3.f90 mod_constants.o
julian.o: julian.f90
mksst2.o: mksst2.f90 mod_datewk.o mod_datenum.o
mksst2a.o: mksst2a.f90
mksst3.o: mksst3.f90
mksst.o: mksst.f90
mkssta.o: mkssta.f90
mxmn3d.o: mxmn3d.f90
p1p2.o: p1p2.f90
printl.o: printl.f90
top2btm.o: top2btm.f90
uvrot4.o: uvrot4.f90
uvrot4nx.o: uvrot4nx.f90
aerosol.o: aerosol.f90 mod_preproc_param.o mod_regcm_param.o
sst_1deg.o: sst_1deg.f90 mod_preproc_param.o mod_regcm_param.o
sst_ccsm.o: sst_ccsm.f90 mod_preproc_param.o mod_regcm_param.o
sst_eh50m.o: sst_eh50m.f90 mod_preproc_param.o mod_regcm_param.o mod_datenum.o
sst_ersst.o: sst_ersst.f90 mod_datenum.o mod_preproc_param.o mod_regcm_param.o
sst_fvgcm.o: sst_fvgcm.f90 mod_preproc_param.o mod_regcm_param.o
icbc.o: icbc.f90 mod_regcm_param.o mod_preproc_param.o mod_datenum.o \
        mod_grid.o mod_write.o mod_ingrid.o mod_datewk.o mod_ein15.o
