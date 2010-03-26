# Makefile for RegCM4
# for PostProc code.
  
include ../Makefile.inc

all: postproc test_read

SHRMODS =  mod_postproc_param.o mod_regcm_param.o mod_constants.o \
           mod_bats_param.o

MODOBJS =  mod_batflds.o   mod_bcflds.o   mod_cheflds.o   mod_outflds.o  \
           mod_point.o     mod_radflds.o  mod_subflds.o

OBJS =     average_day.o   average_month.o  avgdata2d.o  avgdata3d.o \
           avgrain.o       calchgt.o        calcmse2d.o  calcmse.o   \
           calcrh2d.o      calcrh.o         calcslp.o    calcvd.o    \
           fexist.o        fexistnew.o      getminmax.o  height.o    \
           htsig.o         intlin.o         intlog.o     julian.o    \
           nclib.o         one2two.o        param.o      postproc.o  \
           rdhead.o        rdheadicbc.o     setconst.o   two2three.o \
           writebathead.o  writegrads.o     writehead.o  xtdot.o

mod_regcm_param.o: ../Config/mod_regcm_param.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_postproc_param.o: ../Config/mod_postproc_param.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_bats_param.o: ../Config/mod_bats_param.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

postproc: $(SHRMODS) $(MODOBJS) $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMODS) $(MODOBJS) $(OBJS) \
               $(LIB) $(NETCDFLIB)

test_read: test_read.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $< $(LIB) $(NETCDFLIB)

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o postproc test_read

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

average_day.o: average_day.f90
average_month.o: average_month.f90
avgdata2d.o: avgdata2d.f90
avgdata3d.o: avgdata3d.f90
avgrain.o: avgrain.f90
calchgt.o: calchgt.f90 mod_constants.o
calcmse2d.o: calcmse2d.f90 mod_constants.o
calcmse.o: calcmse.f90 mod_constants.o
calcrh2d.o: calcrh2d.f90 mod_constants.o
calcrh.o: calcrh.f90 mod_constants.o
calcslp.o: calcslp.f90 mod_constants.o
calcvd.o: calcvd.f90
fexist.o: fexist.f90
fexistnew.o: fexistnew.f90
getminmax.o: getminmax.f90
height.o: height.f90 mod_constants.o
htsig.o: htsig.f90
intlin.o: intlin.f90
intlog.o: intlog.f90 mod_constants.o
julian.o: julian.f90
mod_batflds.o: mod_batflds.f90 mod_regcm_param.o mod_postproc_param.o \
               mod_point.o
mod_bcflds.o: mod_bcflds.f90 mod_regcm_param.o mod_postproc_param.o \
              mod_point.o
mod_cheflds.o: mod_cheflds.f90 mod_regcm_param.o mod_postproc_param.o \
               mod_point.o
mod_outflds.o: mod_outflds.f90 mod_regcm_param.o mod_postproc_param.o \
               mod_point.o
mod_point.o: mod_point.f90
mod_radflds.o: mod_radflds.f90 mod_regcm_param.o mod_postproc_param.o \
               mod_point.o
mod_subflds.o: mod_subflds.f90 mod_regcm_param.o mod_postproc_param.o \
               mod_point.o
nclib.o: nclib.f90
one2two.o: one2two.f90
param.o: param.f90
postproc.o: postproc.f90 mod_regcm_param.o mod_batflds.o mod_bcflds.o \
            mod_cheflds.o mod_outflds.o mod_radflds.o mod_subflds.o \
            mod_point.o mod_postproc_param.o
rdhead.o: rdhead.f90 mod_regcm_param.o mod_postproc_param.o
rdheadicbc.o: rdheadicbc.f90
setconst.o: setconst.f90
test_read.o: test_read.f90
two2three.o: two2three.f90
writebathead.o: writebathead.f90 mod_regcm_param.o mod_bats_param.o
writegrads.o: writegrads.f90
writehead.o: writehead.f90 mod_regcm_param.o
xtdot.o: xtdot.f90
