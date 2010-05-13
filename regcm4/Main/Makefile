# Makefile for Main regcm Model

###########################################################################
# CONFIGURABLE PATAMETERS
###########################################################################

# Set to 1 to enable CLM coupling
USECLM = 0
# Set to 1 to enable daily cycle SST
USEDCSST = 0
# Set to 1 to enable Sea Ice effects
USESEAICE = 0

SCENARIO = A1B

###########################################################################
# END CONFIGURABLE PARAMETERS
###########################################################################

include ../Makefile.inc

destdir = ../Run
prgsuffix = Serial

# check if we need the MPI wrapper 
ifneq ($(origin PARALLEL), undefined)
  F90 = $(MPIF90)
  prgsuffix = MPI
endif  

ifeq ($(USECLM),1)
CLMPTH = ./clmlib
CLMLIB = $(CLMPTH)/libclm.a
all:: $(CLMLIB) program
prgsuffix = _clM
else
CLMLIB =
all:: program
endif

SHROBJS = mod_constants.o mod_dynparam.o mod_bats_param.o

MODOBJS = \
mod_mppio.o      mod_crdcae.o     mod_mainchem.o   mod_radbuf.o      \
mod_aerosol.o    mod_convect.o    mod_main.o       mod_rad.o         \
mod_bats.o       mod_culookup.o   mod_message.o    mod_slice.o       \
mod_bdycod.o     mod_cvaria.o     mod_o3blk.o      mod_split.o       \
mod_blh_tmp.o    mod_date.o       mod_outrad.o     mod_tiedparam0.o  \
mod_bmparam.o    mod_diagnosis.o  mod_param1.o     mod_tmpsav.o      \
mod_bxq.o        mod_dust.o       mod_param2.o     mod_tracer0.o     \
mod_comctl.o     mod_param3.o     mod_tracer.o     mod_comozp.o      \
mod_when.o       mod_pmoist.o     mod_trachem.o    mod_comtim.o      \
mod_ictp01.o     mod_iunits.o     mod_pbldim.o     mod_lake.o        \
mod_interfaces.o mod_aero_param.o mod_isrchf.o

CLMOBJS = \
initializeMod.o  shr_orb_mod.o    shr_kind_mod.o  clm_varpar.o       \
clm_varsur.o     atmdrvMod.o      spmdMod.o       program_offMod.o   \
clm_comp.o       clmtype.o        perf_mod.o      clm_time_manager.o \
restFileMod.o

OBJS =            cup.o         printsp.o   radclw.o    tractend2.o  \
intmax.o          depth.o       interf.o    radcsw.o    trcab.o      \
radabs.o          deriv.o                   radctl.o    trcabn.o     \
albedov.o         diffu.o       cuparan.o   radded.o    trcems.o     \
radclr.o          diffut.o      lenchr.o    radems.o    trcmix.o     \
bdyin.o           dragdn.o      lfdrag.o    radini.o    trcplk.o     \
bdyuv.o           drag.o        lftemp.o    radinp.o    trcpth.o     \
bdyval.o          drip.o        linpack.o   radout.o    tseaice.o    \
blhnew.o          dsflx.o       ljchrn.o    radoz2.o    tstep.o      \
bmpara.o          zenitm.o      lshfch.o    radtap.o    vadv.o       \
bndry.o           eispack.o     lutbl.o     radtpl.o    vcheke.o     \
carbon.o          eomb.o        mapsmp.o    regcm.o     vcheki.o     \
chdrydep.o        finddate.o    maximi.o    reorder.o   vchekt.o     \
chemtap.o         for_next.o    tracbud.o   root.o      vcover.o     \
chrsetc.o         frawat.o      minimi.o    satur.o     vecbats.o    \
chsrfem.o         getdat.o      tgrund.o    slice1D.o   vmodes.o     \
cldefr.o          grads_stuf.o  mkfile.o    slice3D.o   pcp.o        \
cldems.o          hadv.o        nconvp.o    snow.o      vnorml.o     \
cldfrac.o         header.o      nudge.o     soilbc.o    vorder.o     \
co2.o             holtbl.o      o3data.o    solar1.o    vtlaps.o     \
colmod3.o         htdiff.o      outprt.o    spinit.o    water.o      \
conadv.o                        output.o    spline.o    cupeman.o    \
condch.o          tracdiag.o    outsav.o    splitf.o    invmtrx.o    \
condcq.o          inidust.o     outsrf.o    sponge.o    cupemandrv.o \
condtq.o          inirad.o      outsub.o    spstep.o    tpfc.o       \
conmas.o          initb.o       outtap0.o   stomat.o    zengocndrv.o \
convect_tables.o  initdate.o    outtap.o                zengocn.o    \
cumtran.o         init.o        param.o     tend.o      zenith.o     \
cupara.o

ifeq ($(USECLM),1)
  CPPFLAGS += -DCLM
  MODOBJS += mod_clm.o
  OBJS += initclm.o interfclm.o mod_clm.o mtrxclm.o solar1clm.o albedoclm.o
  LIB += -L$(CLMPTH) -lclm
else
ifeq ($(USEDCSST),1)
  CPPFLAGS += -DDCSST
  dcsstsuffix = _ds
endif
ifeq ($(USESEAICE),1)
  CPPFLAGS += -DSEAICE
  seaicesuffix = _si
endif
endif

program: regcm
	@cp regcm \
          $(destdir)/regcm$(prgsuffix)$(dcsstsuffix)$(seaicesuffix)
	@chmod 755  $(destdir)/regcm$(prgsuffix)$(dcsstsuffix)$(seaicesuffix) 
	@echo "Executable is installed now in $(destdir)"

$(CLMPTH)/libclm.a: mod_dynparam.o mod_clm.o
ifeq ($(origin PARALLEL), undefined)
	@echo "ERROR: CLM wants PARALLEL activated"
	@exit 1
endif
	make -C $(CLMPTH)

regcm: $(SHROBJS) $(CLMLIB) $(MODOBJS) $(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(MODOBJS) $(OBJS) $(SHROBJS)\
               $(LIB) $(NETCDFLIB)

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c -o $@ $<

mod_dynparam.o: ../Config/mod_dynparam.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c -o $@ $<

mod_bats_param.o: ../Config/mod_bats_param.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c -o $@ $<

%.o: %.F90
	$(F90) -c $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c $(F90FLAGS) $(USCORING) -c $<

ifeq ($(USECLM),1)
clean:
	rm -f *.mod *.o regcm
	cd clmlib && make clean
else
clean:
	rm -f *.mod *.o regcm
endif

mod_aero_param.o: mod_aero_param.F90 mod_dynparam.o
mod_aerosol.o: mod_aerosol.F90 mod_constants.o mod_main.o mod_mainchem.o \
               mod_message.o mod_param2.o mod_dynparam.o mod_trachem.o \
               mod_slice.o
mod_bats.o: mod_bats.F90 mod_dynparam.o mod_bats_param.o
mod_bdycod.o: mod_bdycod.F90 mod_dynparam.o
mod_blh_tmp.o: mod_blh_tmp.F90 mod_dynparam.o
mod_bmparam.o: mod_bmparam.f90
mod_bxq.o: mod_bxq.F90 mod_dynparam.o
mod_clm.o: mod_clm.F90 mod_dynparam.o
mod_comctl.o: mod_comctl.f90
mod_comozp.o: mod_comozp.f90
mod_comtim.o: mod_comtim.f90
mod_convect.o: mod_convect.f90
mod_crdcae.o: mod_crdcae.f90
mod_culookup.o: mod_culookup.f90
mod_cvaria.o: mod_cvaria.F90 mod_dynparam.o
mod_date.o: mod_date.f90
mod_diagnosis.o: mod_diagnosis.f90 mod_dynparam.o
mod_dust.o: mod_dust.F90 mod_dynparam.o
mod_ictp01.o: mod_ictp01.F90 mod_dynparam.o
mod_iunits.o: mod_iunits.f90
mod_lake.o: mod_lake.f90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
            mod_dynparam.o
mod_main.o: mod_main.F90 mod_dynparam.o
mod_mainchem.o: mod_mainchem.F90 mod_dynparam.o
mod_message.o: mod_message.F90 mod_dynparam.o
mod_mppio.o: mod_mppio.F90 mod_dynparam.o mod_dust.o
mod_o3blk.o: mod_o3blk.F90 mod_dynparam.o
mod_outrad.o: mod_outrad.F90 mod_dynparam.o
mod_param1.o: mod_param1.f90 mod_dynparam.o
mod_param2.o: mod_param2.F90
mod_param3.o: mod_param3.f90 mod_dynparam.o
mod_pbldim.o: mod_pbldim.F90 mod_dynparam.o
mod_pmoist.o: mod_pmoist.F90 mod_dynparam.o
mod_rad.o: mod_rad.F90 mod_dynparam.o
mod_radbuf.o: mod_radbuf.F90
mod_slice.o: mod_slice.F90 mod_dynparam.o
mod_split.o: mod_split.F90 mod_dynparam.o
mod_tiedparam0.o: mod_tiedparam0.f90
mod_tmpsav.o: mod_tmpsav.F90 mod_dynparam.o
mod_tracer0.o: mod_tracer0.F90
	$(F90) -c $(CPPFLAGS) -D$(SCENARIO) $(F90FLAGS) $(USCORING) -c $<
mod_tracer.o: mod_tracer.f90
mod_trachem.o: mod_trachem.F90 mod_dynparam.o
mod_when.o: mod_when.f90
#
albedoclm.o: albedoclm.F90 mod_bats.o mod_dynparam.o
albedov.o: albedov.F90 mod_constants.o mod_bats.o mod_dynparam.o
bdyin.o: bdyin.F90 mod_bats.o mod_bdycod.o mod_date.o mod_iunits.o mod_main.o \
         mod_param1.o mod_param3.o mod_dynparam.o mod_tmpsav.o
bdyuv.o: bdyuv.F90 mod_bdycod.o mod_cvaria.o mod_main.o mod_dynparam.o
bdyval.o: bdyval.F90 mod_bdycod.o mod_date.o mod_main.o mod_mainchem.o \
          mod_param1.o mod_param2.o mod_dynparam.o
blhnew.o: blhnew.F90 mod_blh_tmp.o mod_constants.o mod_main.o mod_param2.o \
          mod_param3.o mod_pbldim.o mod_dynparam.o mod_slice.o
bmpara.o: bmpara.F90 mod_bats.o mod_bmparam.o mod_constants.o mod_date.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o \
          mod_rad.o mod_dynparam.o mod_trachem.o
bndry.o: bndry.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_param1.o \
         mod_dynparam.o
carbon.o: carbon.f90
chdrydep.o: chdrydep.f90 mod_constants.o mod_aero_param.o
chemtap.o: chemtap.F90 mod_date.o mod_iunits.o mod_main.o mod_mainchem.o \
           mod_param2.o mod_dynparam.o mod_trachem.o
chrsetc.o: chrsetc.f90
ifeq ($(USECLM),1)
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_iunits.o mod_mainchem.o \
           mod_message.o mod_dynparam.o mod_trachem.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_iunits.o mod_mainchem.o \
           mod_message.o mod_dynparam.o mod_trachem.o
endif
cldefr.o: cldefr.f90 mod_dynparam.o
cldems.o: cldems.f90 mod_dynparam.o
cldfrac.o: cldfrac.f90 mod_pmoist.o mod_rad.o mod_dynparam.o mod_slice.o
co2.o: co2.f90 mod_bats.o mod_ictp01.o mod_param1.o mod_dynparam.o
colmod3.o: colmod3.F90 mod_constants.o mod_bats.o mod_comctl.o mod_comtim.o \
           mod_dynparam.o
conadv.o: conadv.F90 mod_constants.o mod_diagnosis.o mod_main.o mod_param1.o \
          mod_param3.o mod_dynparam.o
condch.o: condch.f90 mod_constants.o mod_bats.o mod_dynparam.o
condcq.o: condcq.f90 mod_bats.o mod_dynparam.o
condtq.o: condtq.f90 mod_constants.o mod_cvaria.o mod_main.o mod_param1.o \
          mod_param3.o mod_pmoist.o mod_dynparam.o mod_slice.o
conmas.o: conmas.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
          mod_param1.o mod_param2.o mod_param3.o mod_dynparam.o
convect_tables.o: convect_tables.f90 mod_constants.o mod_culookup.o \
                  mod_tiedparam0.o
cumtran.o: cumtran.F90 mod_mainchem.o mod_param3.o mod_dynparam.o \
           mod_trachem.o
cup.o: cup.f90 mod_param2.o mod_pmoist.o mod_rad.o mod_dynparam.o \
       mod_trachem.o
cupara.o: cupara.f90 mod_bats.o mod_constants.o mod_cvaria.o mod_date.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o \
          mod_rad.o mod_dynparam.o mod_trachem.o
cuparan.o: cuparan.f90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
           mod_param1.o mod_param3.o mod_pmoist.o mod_rad.o mod_dynparam.o
cupeman.o: cupeman.f90 mod_convect.o
cupemandrv.o: cupemandrv.f90 mod_bats.o mod_cvaria.o mod_date.o mod_main.o \
              mod_param1.o mod_param3.o mod_pmoist.o mod_rad.o \
              mod_dynparam.o mod_slice.o
depth.o: depth.f90 mod_bats.o mod_dynparam.o
deriv.o: deriv.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_dynparam.o
diffu.o: diffu.F90 mod_main.o mod_dynparam.o mod_slice.o
diffut.o: diffut.F90 mod_main.o mod_dynparam.o mod_slice.o
drag.o: drag.f90 mod_constants.o mod_bats.o mod_dynparam.o
dragdn.o: dragdn.F90 mod_constants.o mod_bats.o mod_dynparam.o
drip.o: drip.f90 mod_constants.o mod_bats.o mod_dynparam.o
dsflx.o: dsflx.F90 mod_constants.o mod_dust.o mod_message.o mod_dynparam.o
eispack.o: eispack.f90
eomb.o: eomb.f90 mod_constants.o
finddate.o: finddate.F90 mod_message.o mod_dynparam.o
for_next.o: for_next.F90 mod_convect.o mod_date.o mod_param1.o mod_param2.o \
            mod_param3.o mod_pmoist.o mod_dynparam.o mod_trachem.o
frawat.o: frawat.f90 mod_constants.o mod_bats.o mod_dynparam.o
getdat.o: getdat.f90 mod_comtim.o mod_constants.o mod_date.o \
          mod_main.o mod_param1.o mod_param3.o mod_rad.o mod_dynparam.o
grads_stuf.o: grads_stuf.F90 mod_date.o mod_iunits.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_dynparam.o
hadv.o: hadv.F90 mod_cvaria.o mod_main.o mod_dynparam.o
header.o: header.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
holtbl.o: holtbl.F90 mod_bats.o mod_blh_tmp.o mod_constants.o mod_cvaria.o \
          mod_diagnosis.o mod_main.o mod_mainchem.o mod_param1.o mod_param2.o \
          mod_param3.o mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o \
          mod_trachem.o
htdiff.o: htdiff.F90 mod_pmoist.o mod_dynparam.o
ifport.o: ifport.f90
inidust.o: inidust.F90 mod_constants.o mod_dust.o mod_dynparam.o
inirad.o: inirad.F90 mod_date.o mod_rad.o mod_dynparam.o
ifeq ($(USECLM),1)
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_iunits.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_radbuf.o mod_dynparam.o mod_tmpsav.o mod_trachem.o \
        $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_iunits.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_radbuf.o mod_dynparam.o mod_tmpsav.o mod_trachem.o
endif
initb.o: initb.F90 mod_bats.o mod_main.o mod_param2.o mod_dynparam.o
initdate.o: initdate.F90 mod_message.o mod_param1.o mod_dynparam.o
interf.o: interf.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
          mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_interfaces.o \
          mod_slice.o
intmax.o: intmax.f90
invmtrx.o: invmtrx.F90 mod_message.o
mod_isrchf.o: mod_isrchf.f90
lenchr.o: lenchr.f90
lfdrag.o: lfdrag.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_dynparam.o
lftemp.o: lftemp.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_param1.o \
          mod_dynparam.o
linpack.o: linpack.f90
ljchrn.o: ljchrn.f90
lshfch.o: lshfch.F90
lutbl.o: lutbl.f90 mod_constants.o mod_bmparam.o
mapsmp.o: mapsmp.f90
maximi.o: maximi.f90
minimi.o: minimi.f90
mkfile.o: mkfile.F90 mod_bats.o mod_date.o mod_iunits.o mod_mainchem.o \
          mod_outrad.o mod_param1.o mod_param2.o mod_dynparam.o
nconvp.o: nconvp.F90 mod_constants.o mod_param3.o mod_pmoist.o mod_dynparam.o
nudge.o: nudge.F90 mod_bdycod.o mod_main.o mod_param3.o mod_dynparam.o
o3data.o: o3data.F90 mod_main.o mod_o3blk.o mod_param3.o mod_rad.o \
          mod_dynparam.o
outprt.o: outprt.F90 mod_bdycod.o mod_constants.o mod_date.o mod_main.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_dynparam.o
ifeq ($(USECLM),1)
output.o: output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_rad.o \
          mod_radbuf.o mod_dynparam.o mod_split.o mod_tmpsav.o \
          mod_trachem.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
output.o: output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_rad.o \
          mod_radbuf.o mod_dynparam.o mod_split.o mod_tmpsav.o mod_trachem.o
endif
outsav.o: outsav.F90 mod_bats.o mod_bdycod.o mod_date.o mod_diagnosis.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_param2.o mod_pmoist.o \
          mod_rad.o mod_radbuf.o mod_dynparam.o mod_split.o mod_trachem.o
outsrf.o: outsrf.F90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
          mod_param2.o mod_dynparam.o
outsub.o: outsub.F90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
          mod_param2.o mod_dynparam.o
outtap.o: outtap.F90 mod_bats.o mod_cvaria.o mod_date.o mod_iunits.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_dynparam.o
outtap0.o: outtap0.F90 mod_constants.o mod_bats.o mod_date.o \
           mod_main.o mod_param2.o mod_param3.o mod_dynparam.o
param.o: param.F90 mod_constants.o mod_bats.o mod_convect.o mod_date.o \
         mod_iunits.o mod_main.o mod_message.o mod_param1.o \
         mod_param2.o mod_param3.o mod_pmoist.o mod_dynparam.o mod_trachem.o
pcp.o: pcp.f90 mod_bats.o mod_constants.o mod_cvaria.o mod_date.o mod_main.o \
       mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_dynparam.o \
       mod_slice.o mod_trachem.o
printsp.o: printsp.f90
radabs.o: radabs.f90 mod_crdcae.o mod_constants.o mod_radbuf.o \
          mod_dynparam.o
radclr.o: radclr.f90 mod_aerosol.o mod_dynparam.o mod_when.o
radclw.o: radclw.f90 mod_constants.o mod_date.o mod_param1.o mod_param2.o \
          mod_radbuf.o mod_dynparam.o mod_when.o
radcsw.o: radcsw.f90 mod_aerosol.o mod_bats.o mod_constants.o mod_param2.o \
          mod_dynparam.o mod_tracer.o
radctl.o: radctl.f90 mod_aerosol.o mod_comtim.o mod_dynparam.o \
          mod_interfaces.o
radded.o: radded.f90 mod_dynparam.o mod_when.o
radems.o: radems.f90 mod_crdcae.o mod_constants.o mod_radbuf.o \
          mod_dynparam.o
radini.o: radini.F90 mod_comozp.o mod_crdcae.o mod_constants.o mod_date.o \
          mod_message.o mod_dynparam.o mod_tracer.o mod_tracer0.o
ifeq ($(USECLM),1)
radinp.o: radinp.F90 mod_comtim.o mod_constants.o mod_dynparam.o mod_clm.o
else
radinp.o: radinp.F90 mod_comtim.o mod_constants.o mod_dynparam.o
endif
radout.o: radout.F90 mod_bats.o mod_date.o mod_main.o mod_outrad.o \
          mod_param1.o mod_param2.o mod_param3.o mod_rad.o mod_dynparam.o
radoz2.o: radoz2.f90 mod_comozp.o mod_dynparam.o
radtap.o: radtap.F90 mod_date.o mod_iunits.o mod_outrad.o mod_param1.o \
          mod_param2.o mod_radbuf.o mod_dynparam.o
radtpl.o: radtpl.f90 mod_dynparam.o mod_constants.o
ifeq ($(USECLM),1)
regcm.o: regcm.F90 mod_date.o mod_message.o mod_param1.o mod_param2.o \
         mod_param3.o mod_dynparam.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
regcm.o: regcm.F90 mod_date.o mod_message.o mod_param1.o mod_param2.o \
         mod_param3.o mod_dynparam.o
endif
reorder.o: reorder.f90
root.o: root.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_dynparam.o
satur.o: satur.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_dynparam.o
slice1D.o: slice1D.f90 mod_bats.o mod_constants.o mod_main.o mod_pbldim.o \
           mod_dynparam.o
slice3D.o: slice3D.F90 mod_main.o mod_constants.o mod_mainchem.o mod_param2.o \
           mod_param3.o mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o
snow.o: snow.f90 mod_constants.o mod_bats.o mod_param1.o mod_dynparam.o
soilbc.o: soilbc.f90 mod_bats.o mod_dynparam.o
solar1.o: solar1.f90 mod_constants.o mod_date.o mod_message.o
spinit.o: spinit.F90 mod_bdycod.o mod_bxq.o mod_date.o mod_iunits.o mod_main.o \
          mod_param1.o mod_param2.o mod_dynparam.o mod_split.o mod_tmpsav.o \
          mod_constants.o
spline.o: spline.f90
splitf.o: splitf.F90 mod_bxq.o mod_main.o mod_param1.o mod_dynparam.o \
          mod_split.o
sponge.o: sponge.F90 mod_bdycod.o mod_dynparam.o
spstep.o: spstep.F90 mod_bxq.o mod_main.o mod_dynparam.o
stomat.o: stomat.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_dynparam.o
ifeq ($(USECLM),1)
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_dynparam.o mod_slice.o mod_trachem.o mod_interfaces.o mod_clm.o \
        $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_dynparam.o mod_slice.o mod_trachem.o mod_interfaces.o
endif
tgrund.o: tgrund.F90 mod_constants.o mod_bats.o mod_param1.o mod_dynparam.o
tpfc.o: tpfc.f90 mod_constants.o
tracbud.o: tracbud.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
           mod_mainchem.o mod_param1.o mod_param2.o mod_param3.o \
           mod_dynparam.o mod_trachem.o
tracdiag.o: tracdiag.F90 mod_constants.o mod_diagnosis.o mod_main.o \
            mod_mainchem.o mod_param1.o mod_param3.o mod_dynparam.o
tractend2.o: tractend2.f90 mod_constants.o mod_bats.o mod_cvaria.o mod_date.o \
             mod_dust.o mod_main.o mod_mainchem.o mod_message.o mod_param1.o \
             mod_param3.o mod_pbldim.o mod_pmoist.o mod_rad.o \
             mod_dynparam.o mod_trachem.o
trcab.o: trcab.f90 mod_constants.o mod_dynparam.o
trcabn.o: trcabn.f90 mod_constants.o mod_dynparam.o
trcems.o: trcems.f90 mod_constants.o mod_dynparam.o
trcmix.o: trcmix.f90 mod_dynparam.o mod_tracer.o
trcplk.o: trcplk.f90 mod_constants.o mod_dynparam.o
trcpth.o: trcpth.f90 mod_constants.o mod_dynparam.o
tseaice.o: tseaice.F90 mod_constants.o mod_bats.o mod_param1.o mod_dynparam.o
tstep.o: tstep.f90 mod_date.o mod_param1.o mod_dynparam.o
vadv.o: vadv.f90 mod_cvaria.o mod_main.o mod_param3.o mod_dynparam.o
vcheke.o: vcheke.f90
vcheki.o: vcheki.f90
vchekt.o: vchekt.f90
vcover.o: vcover.f90 mod_bats.o mod_dynparam.o
vecbats.o: vecbats.F90 mod_date.o mod_param1.o mod_param2.o mod_interfaces.o
vmodes.o: vmodes.F90 mod_message.o mod_param1.o mod_dynparam.o mod_split.o
vmultm.o: vmultm.f90
vnorml.o: vnorml.f90
vorder.o: vorder.f90
vtlaps.o: vtlaps.f90 mod_constants.o
water.o: water.f90 mod_constants.o mod_bats.o mod_param1.o mod_dynparam.o
zengocn.o: zengocn.f90 mod_constants.o
ifeq ($(USECLM),1)
zengocndrv.o: zengocndrv.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_slice.o \
              mod_bdycod.o mod_interfaces.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
zengocndrv.o: zengocndrv.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_slice.o \
              mod_interfaces.o mod_dynparam.o mod_bdycod.o
endif
zenith.o: zenith.f90 mod_constants.o
ifeq ($(USECLM),1)
zenitm.o: zenitm.F90 mod_dynparam.o mod_clm.o mod_date.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
zenitm.o: zenitm.F90 mod_constants.o mod_date.o mod_main.o mod_dynparam.o
endif
initclm.o: initclm.F90 $(CLMLIB) mod_constants.o mod_dynparam.o mod_param1.o \
           mod_mppio.o mod_date.o mod_param2.o mod_param3.o mod_clm.o \
           mod_bats.o mod_mppio.o mod_main.o mod_pbldim.o mod_slice.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
interfclm.o: interfclm.F90 $(CLMPTH)/clm_varsur.o $(CLMPTH)/clmtype.o \
             mod_clm.o mod_dynparam.o mod_date.o mod_main.o mod_param1.o \
             mod_param2.o mod_param3.o mod_slice.o mod_pbldim.o mod_bats.o \
             mod_constants.o mod_mppio.o mod_interfaces.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
mtrxclm.o: mtrxclm.F90 $(CLMLIB) mod_message.o mod_clm.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
solar1clm.o: solar1clm.f90 $(CLMLIB) mod_clm.o mod_message.o mod_date.o \
             mod_constants.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
