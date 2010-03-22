include ../Makefile.inc

USECLM = 1

SCENARIO = A1B
SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

all:: regcm

MODOBJS = \
mod_mppio.o      mod_crdcae.o     mod_mainchem.o   mod_radbuf.o      \
mod_aerosol.o    mod_convect.o    mod_main.o       mod_rad.o         \
mod_bats.o       mod_culookup.o   mod_message.o    mod_slice.o       \
mod_bdycod.o     mod_cvaria.o     mod_o3blk.o      mod_split.o       \
mod_blh_tmp.o    mod_date.o       mod_outrad.o     mod_tiedparam0.o  \
mod_bmparam.o    mod_diagnosis.o  mod_param1.o     mod_tmpsav.o      \
mod_bxq.o        mod_dust.o       mod_param2.o     mod_tracer0.o     \
mod_comctl.o     mod_param3.o     mod_tracer.o     mod_comozp.o      \
mod_grads.o      mod_pmoist.o     mod_trachem.o    mod_comtim.o      \
mod_ictp01.o     mod_constants.o  mod_iunits.o     mod_regcm_param.o \
mod_pbldim.o     mod_lake.o       mod_interfaces.o

CLMOBJS = \
initializeMod.o  shr_orb_mod.o    shr_kind_mod.o  clm_varpar.o       \
clm_varsur.o     atmdrvMod.o      spmdMod.o       program_offMod.o   \
clm_comp.o       clmtype.o        perf_mod.o      clm_time_manager.o \
restFileMod.o

OBJS = \
aermix.o          cup.o         isrchfgt.o  radclw.o    tractend2.o  \
aeroppt.o         depth.o       isrchfle.o  radcsw.o    trcab.o      \
aerout.o          deriv.o       lakedrv.o   radctl.o    trcabn.o     \
albedov.o         diffu.o       lake.o      radded.o    trcems.o     \
radclr.o          diffut.o      lenchr.o    radems.o    trcmix.o     \
bdyin.o           dragdn.o      lfdrag.o    radini.o    trcplk.o     \
bdyuv.o           drag.o        lftemp.o    radinp.o    trcpth.o     \
bdyval.o          drip.o        linpack.o   radout.o    tseaice.o    \
blhnew.o          dsflx.o       ljchrn.o    radoz2.o    tstep.o      \
bmpara.o          eddy.o        lshfch.o    radtap.o    vadv.o       \
bndry.o           eispack.o     lutbl.o     radtpl.o    vcheke.o     \
carbon.o          eomb.o        mapsmp.o    regcm.o     vcheki.o     \
chdrydep.o        finddate.o    maximi.o    reorder.o   vchekt.o     \
chemtap.o         for_next.o    message.o   root.o      vcover.o     \
chrsetc.o         frawat.o      minimi.o    satur.o     vecbats.o    \
chsrfem.o         getdat.o      mixer.o     slice1D.o   vmodes.o     \
cldefr.o          grads_stuf.o  mkfile.o    slice3D.o   vmultm.o     \
cldems.o          hadv.o        nconvp.o    snow.o      vnorml.o     \
cldfrac.o         header.o      nudge.o     soilbc.o    vorder.o     \
co2.o             holtbl.o      o3data.o    solar1.o    vtlaps.o     \
colmod3.o         htdiff.o      outprt.o    spinit.o    water.o      \
conadv.o          ice.o         output.o    spline.o    wheneq.o     \
condch.o          tracdiag.o    outsav.o    splitf.o    whenfgt.o    \
condcq.o          inidust.o     outsrf.o    sponge.o    whenflt.o    \
condtq.o          inirad.o      outsub.o    spstep.o    whenne.o     \
conmas.o          initb.o       outtap0.o   stomat.o    zengocndrv.o \
convect_tables.o  initdate.o    outtap.o    temp.o      zengocn.o    \
cumtran.o         init.o        param.o     tend.o      zenith.o     \
cupara.o          initlk.o      pcp.o       tgrund.o    zenitm.o     \
cuparan.o         interf.o      printsp.o   tpfc.o      cupemandrv.o \
intmax.o          radabs.o      tracbud.o   cupeman.o   invmtrx.o

ifeq ($(USECLM),1)
  CLMPTH = ./clmlib
  CPPFLAGS += -DCLM
  MODOBJS += mod_clm.o
  OBJS += initclm.o interfclm.o mod_clm.o mtrxclm.o solar1clm.o albedoclm.o
  LIB += -L$(CLMPTH) -lclm
  CLMLIB = $(CLMPTH)/libclm.a
else
  CLMLIB =
endif

$(CLMPTH)/libclm.a: mod_clm.o mod_regcm_param.o
	cd $(CLMPTH) && touch patchset/* && make patchsrc && make

regcm: $(CLMLIB) $(MODOBJS) $(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(MODOBJS) $(OBJS) \
               $(LIB) $(NETCDFLIB)

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

mod_aerosol.o: mod_aerosol.F90 mod_regcm_param.o
mod_bats.o: mod_bats.F90 mod_regcm_param.o
mod_bdycod.o: mod_bdycod.F90 mod_regcm_param.o
mod_blh_tmp.o: mod_blh_tmp.F90 mod_regcm_param.o
mod_bmparam.o: mod_bmparam.f90
mod_bxq.o: mod_bxq.F90 mod_regcm_param.o
mod_clm.o: mod_clm.f90 mod_regcm_param.o
mod_comctl.o: mod_comctl.f90
mod_comozp.o: mod_comozp.f90
mod_comtim.o: mod_comtim.f90
mod_constants.o: mod_constants.F90
mod_convect.o: mod_convect.f90
mod_crdcae.o: mod_crdcae.f90
mod_culookup.o: mod_culookup.f90
mod_cvaria.o: mod_cvaria.F90 mod_regcm_param.o
mod_date.o: mod_date.f90
mod_diagnosis.o: mod_diagnosis.f90 mod_regcm_param.o
mod_dust.o: mod_dust.F90 mod_regcm_param.o
mod_grads.o: mod_grads.f90
mod_ictp01.o: mod_ictp01.F90 mod_regcm_param.o
mod_iunits.o: mod_iunits.f90
mod_lake.o: mod_lake.f90
mod_main.o: mod_main.F90 mod_regcm_param.o
mod_mainchem.o: mod_mainchem.F90 mod_regcm_param.o
mod_message.o: mod_message.f90
mod_mppio.o: mod_mppio.F90 mod_regcm_param.o
mod_o3blk.o: mod_o3blk.F90 mod_regcm_param.o
mod_outrad.o: mod_outrad.F90 mod_regcm_param.o
mod_param1.o: mod_param1.f90 mod_regcm_param.o
mod_param2.o: mod_param2.F90
mod_param3.o: mod_param3.f90 mod_regcm_param.o
mod_regcm_param.o: mod_regcm_param.F90
mod_pbldim.o: mod_pbldim.F90 mod_regcm_param.o
mod_pmoist.o: mod_pmoist.F90 mod_regcm_param.o
mod_rad.o: mod_rad.F90 mod_regcm_param.o
mod_radbuf.o: mod_radbuf.F90
mod_slice.o: mod_slice.F90 mod_regcm_param.o
mod_split.o: mod_split.F90 mod_regcm_param.o
mod_tiedparam0.o: mod_tiedparam0.f90
mod_tmpsav.o: mod_tmpsav.F90 mod_regcm_param.o
mod_tracer0.o: mod_tracer0.F90
	$(F90) -c $(CPPFLAGS) -D$(SCENARIO) $(F90FLAGS) $(USCORING) -c $<
mod_tracer.o: mod_tracer.f90
mod_trachem.o: mod_trachem.F90 mod_regcm_param.o
aermix.o: aermix.f90 mod_aerosol.o mod_main.o mod_mainchem.o \
          mod_param2.o mod_slice.o mod_constants.o mod_interfaces.o
aeroppt.o: aeroppt.F90 mod_aerosol.o mod_constants.o \
           mod_message.o mod_param2.o mod_regcm_param.o mod_trachem.o
aerout.o: aerout.F90 mod_aerosol.o mod_param2.o mod_regcm_param.o \
          mod_trachem.o
albedoclm.o: albedoclm.F90 mod_bats.o mod_regcm_param.o
albedov.o: albedov.f90 mod_constants.o mod_bats.o mod_regcm_param.o
bdyin.o: bdyin.F90 mod_bats.o mod_bdycod.o mod_date.o mod_iunits.o mod_main.o \
         mod_param1.o mod_param3.o mod_regcm_param.o mod_tmpsav.o
bdyuv.o: bdyuv.F90 mod_bdycod.o mod_cvaria.o mod_main.o mod_regcm_param.o
bdyval.o: bdyval.F90 mod_bdycod.o mod_date.o mod_main.o mod_mainchem.o \
          mod_param1.o mod_param2.o mod_regcm_param.o
blhnew.o: blhnew.F90 mod_blh_tmp.o mod_constants.o mod_main.o mod_param2.o \
          mod_param3.o mod_pbldim.o mod_regcm_param.o mod_slice.o
bmpara.o: bmpara.F90 mod_bats.o mod_bmparam.o mod_constants.o mod_date.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o \
          mod_rad.o mod_regcm_param.o mod_trachem.o
bndry.o: bndry.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_param1.o \
         mod_regcm_param.o
carbon.o: carbon.f90
chdrydep.o: chdrydep.f90 mod_constants.o
chemtap.o: chemtap.F90 mod_date.o mod_iunits.o mod_main.o mod_mainchem.o \
           mod_param2.o mod_regcm_param.o mod_trachem.o
chrsetc.o: chrsetc.f90
ifeq ($(USECLM),1)
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_iunits.o mod_mainchem.o \
           mod_message.o mod_regcm_param.o mod_trachem.o \
        $(CLMPTH)/clm_varsur.o $(CLMPTH)/surfrdMod.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_iunits.o mod_mainchem.o \
           mod_message.o mod_regcm_param.o mod_trachem.o
endif
cldefr.o: cldefr.f90 mod_regcm_param.o
cldems.o: cldems.f90 mod_regcm_param.o
cldfrac.o: cldfrac.f90 mod_pmoist.o mod_rad.o mod_regcm_param.o mod_slice.o
co2.o: co2.f90 mod_bats.o mod_ictp01.o mod_param1.o mod_regcm_param.o
colmod3.o: colmod3.F90 mod_constants.o mod_bats.o mod_comctl.o mod_comtim.o \
           mod_regcm_param.o
conadv.o: conadv.F90 mod_constants.o mod_diagnosis.o mod_main.o mod_param1.o \
          mod_param3.o mod_regcm_param.o
condch.o: condch.f90 mod_constants.o mod_bats.o mod_regcm_param.o
condcq.o: condcq.f90 mod_bats.o mod_regcm_param.o
condtq.o: condtq.f90 mod_constants.o mod_cvaria.o mod_main.o mod_param1.o \
          mod_param3.o mod_pmoist.o mod_regcm_param.o mod_slice.o
conmas.o: conmas.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
          mod_param1.o mod_param2.o mod_param3.o mod_regcm_param.o
convect_tables.o: convect_tables.f90 mod_constants.o mod_culookup.o \
                  mod_tiedparam0.o
cumtran.o: cumtran.F90 mod_mainchem.o mod_param3.o mod_regcm_param.o \
           mod_trachem.o
cup.o: cup.f90 mod_param2.o mod_pmoist.o mod_rad.o mod_regcm_param.o \
       mod_trachem.o
cupara.o: cupara.f90 mod_bats.o mod_constants.o mod_cvaria.o mod_date.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o \
          mod_rad.o mod_regcm_param.o mod_trachem.o
cuparan.o: cuparan.f90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
           mod_param1.o mod_param3.o mod_pmoist.o mod_rad.o mod_regcm_param.o
cupeman.o: cupeman.f90 mod_convect.o
cupemandrv.o: cupemandrv.f90 mod_bats.o mod_cvaria.o mod_date.o mod_main.o \
              mod_param1.o mod_param3.o mod_pmoist.o mod_rad.o \
              mod_regcm_param.o mod_slice.o
depth.o: depth.f90 mod_bats.o mod_regcm_param.o
deriv.o: deriv.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_regcm_param.o
diffu.o: diffu.F90 mod_main.o mod_regcm_param.o mod_slice.o
diffut.o: diffut.F90 mod_main.o mod_regcm_param.o mod_slice.o
drag.o: drag.f90 mod_constants.o mod_bats.o mod_regcm_param.o
dragdn.o: dragdn.f90 mod_constants.o mod_bats.o mod_regcm_param.o
drip.o: drip.f90 mod_constants.o mod_bats.o mod_regcm_param.o
dsflx.o: dsflx.F90 mod_constants.o mod_dust.o mod_message.o mod_regcm_param.o
eddy.o: eddy.f90 mod_constants.o
eispack.o: eispack.f90
eomb.o: eomb.f90 mod_constants.o
finddate.o: finddate.F90 mod_message.o mod_regcm_param.o
for_next.o: for_next.F90 mod_convect.o mod_date.o mod_param1.o mod_param2.o \
            mod_param3.o mod_pmoist.o mod_regcm_param.o mod_trachem.o
frawat.o: frawat.f90 mod_constants.o mod_bats.o mod_regcm_param.o
getdat.o: getdat.f90 mod_comtim.o mod_constants.o mod_date.o \
          mod_main.o mod_param1.o mod_param3.o mod_rad.o mod_regcm_param.o
grads_stuf.o: grads_stuf.F90 mod_date.o mod_grads.o mod_iunits.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_regcm_param.o
hadv.o: hadv.F90 mod_cvaria.o mod_main.o mod_regcm_param.o
header.o: header.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
holtbl.o: holtbl.F90 mod_bats.o mod_blh_tmp.o mod_constants.o mod_cvaria.o \
          mod_diagnosis.o mod_main.o mod_mainchem.o mod_param1.o mod_param2.o \
          mod_param3.o mod_pbldim.o mod_pmoist.o mod_regcm_param.o mod_slice.o \
          mod_trachem.o
htdiff.o: htdiff.F90 mod_pmoist.o mod_regcm_param.o
ice.o: ice.f90 mod_constants.o
ifport.o: ifport.f90
inidust.o: inidust.F90 mod_constants.o mod_dust.o mod_regcm_param.o
inirad.o: inirad.F90 mod_date.o mod_rad.o mod_regcm_param.o
ifeq ($(USECLM),1)
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_iunits.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_radbuf.o mod_regcm_param.o mod_tmpsav.o mod_trachem.o \
        $(CLMPTH)/clm_varsur.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_iunits.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_radbuf.o mod_regcm_param.o mod_tmpsav.o mod_trachem.o
endif
initb.o: initb.F90 mod_bats.o mod_main.o mod_param2.o mod_regcm_param.o
initdate.o: initdate.F90 mod_message.o mod_param1.o mod_regcm_param.o
initlk.o: initlk.f90 mod_lake.o mod_regcm_param.o
interf.o: interf.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
          mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_interfaces.o \
          mod_slice.o
intmax.o: intmax.f90
invmtrx.o: invmtrx.F90 mod_message.o
isrchfgt.o: isrchfgt.f90
isrchfle.o: isrchfle.f90
lake.o: lake.f90 mod_constants.o mod_lake.o
lakedrv.o: lakedrv.f90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
           mod_regcm_param.o
lenchr.o: lenchr.f90
lfdrag.o: lfdrag.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_regcm_param.o
lftemp.o: lftemp.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_param1.o \
          mod_regcm_param.o
linpack.o: linpack.f90
ljchrn.o: ljchrn.f90
lshfch.o: lshfch.F90
lutbl.o: lutbl.f90 mod_constants.o mod_bmparam.o
mapsmp.o: mapsmp.f90
maximi.o: maximi.f90
message.o: message.F90 mod_message.o mod_regcm_param.o
minimi.o: minimi.f90
mixer.o: mixer.f90 mod_constants.o
mkfile.o: mkfile.F90 mod_bats.o mod_date.o mod_iunits.o mod_mainchem.o \
          mod_outrad.o mod_param1.o mod_param2.o mod_regcm_param.o
nconvp.o: nconvp.F90 mod_constants.o mod_param3.o mod_pmoist.o mod_regcm_param.o
nudge.o: nudge.F90 mod_bdycod.o mod_main.o mod_param3.o mod_regcm_param.o
o3data.o: o3data.F90 mod_main.o mod_o3blk.o mod_param3.o mod_rad.o \
          mod_regcm_param.o
outprt.o: outprt.F90 mod_bdycod.o mod_constants.o mod_date.o mod_main.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_regcm_param.o
ifeq ($(USECLM),1)
output.o: output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_rad.o \
          mod_radbuf.o mod_regcm_param.o mod_split.o mod_tmpsav.o \
          mod_trachem.o $(CLMPTH)/clm_varctl.o $(CLMPTH)/restFileMod.o \
          $(CLMPTH)/clm_time_manager.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
output.o: output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_param1.o mod_param2.o mod_pmoist.o mod_rad.o \
          mod_radbuf.o mod_regcm_param.o mod_split.o mod_tmpsav.o mod_trachem.o
endif
outsav.o: outsav.F90 mod_bats.o mod_bdycod.o mod_date.o mod_diagnosis.o \
          mod_iunits.o mod_main.o mod_mainchem.o mod_param2.o mod_pmoist.o \
          mod_rad.o mod_radbuf.o mod_regcm_param.o mod_split.o mod_trachem.o
outsrf.o: outsrf.F90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
          mod_param2.o mod_regcm_param.o
outsub.o: outsub.F90 mod_bats.o mod_date.o mod_iunits.o mod_param1.o \
          mod_param2.o mod_regcm_param.o
outtap.o: outtap.F90 mod_bats.o mod_cvaria.o mod_date.o mod_iunits.o \
          mod_main.o mod_param1.o mod_param2.o mod_param3.o mod_regcm_param.o
outtap0.o: outtap0.F90 mod_constants.o mod_bats.o mod_date.o mod_grads.o \
           mod_main.o mod_param2.o mod_param3.o mod_regcm_param.o
param.o: param.F90 mod_constants.o mod_bats.o mod_convect.o mod_date.o \
         mod_grads.o mod_iunits.o mod_main.o mod_message.o mod_param1.o \
         mod_param2.o mod_param3.o mod_pmoist.o mod_regcm_param.o mod_trachem.o
pcp.o: pcp.f90 mod_bats.o mod_constants.o mod_cvaria.o mod_date.o mod_main.o \
       mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_regcm_param.o \
       mod_slice.o mod_trachem.o
printsp.o: printsp.f90
radabs.o: radabs.f90 mod_crdcae.o mod_constants.o mod_radbuf.o \
          mod_regcm_param.o
radclr.o: radclr.f90 mod_aerosol.o mod_regcm_param.o
radclw.o: radclw.f90 mod_constants.o mod_date.o mod_param1.o mod_param2.o \
          mod_radbuf.o mod_regcm_param.o
radcsw.o: radcsw.f90 mod_aerosol.o mod_bats.o mod_constants.o mod_param2.o \
          mod_regcm_param.o mod_tracer.o
radctl.o: radctl.f90 mod_aerosol.o mod_comtim.o mod_regcm_param.o \
          mod_interfaces.o
radded.o: radded.f90 mod_regcm_param.o
radems.o: radems.f90 mod_crdcae.o mod_constants.o mod_radbuf.o \
          mod_regcm_param.o
radini.o: radini.F90 mod_comozp.o mod_crdcae.o mod_constants.o mod_date.o \
          mod_message.o mod_regcm_param.o mod_tracer.o mod_tracer0.o
ifeq ($(USECLM),1)
radinp.o: radinp.F90 mod_comtim.o mod_constants.o mod_regcm_param.o mod_clm.o
else
radinp.o: radinp.F90 mod_comtim.o mod_constants.o mod_regcm_param.o
endif
radout.o: radout.F90 mod_bats.o mod_date.o mod_outrad.o mod_param1.o \
          mod_param2.o mod_rad.o mod_regcm_param.o
radoz2.o: radoz2.f90 mod_comozp.o mod_regcm_param.o
radtap.o: radtap.F90 mod_date.o mod_iunits.o mod_outrad.o mod_param1.o \
          mod_param2.o mod_radbuf.o mod_regcm_param.o
radtpl.o: radtpl.f90 mod_regcm_param.o mod_constants.o
ifeq ($(USECLM),1)
regcm.o: regcm.F90 mod_date.o mod_message.o mod_param1.o mod_param2.o \
         mod_regcm_param.o $(CLMPTH)/perf_mod.o $(CLMPTH)/spmdMod.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
regcm.o: regcm.F90 mod_date.o mod_message.o mod_param1.o mod_param2.o \
         mod_regcm_param.o
endif
reorder.o: reorder.f90
root.o: root.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_regcm_param.o
satur.o: satur.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_regcm_param.o
slice1D.o: slice1D.f90 mod_bats.o mod_constants.o mod_main.o mod_pbldim.o \
           mod_regcm_param.o
slice3D.o: slice3D.F90 mod_main.o mod_constants.o mod_mainchem.o mod_param2.o \
           mod_param3.o mod_pbldim.o mod_pmoist.o mod_regcm_param.o mod_slice.o
snow.o: snow.f90 mod_constants.o mod_bats.o mod_param1.o mod_regcm_param.o
soilbc.o: soilbc.f90 mod_bats.o mod_regcm_param.o
solar1.o: solar1.f90 mod_constants.o mod_date.o mod_message.o
spinit.o: spinit.F90 mod_bdycod.o mod_bxq.o mod_date.o mod_iunits.o mod_main.o \
          mod_param1.o mod_param2.o mod_regcm_param.o mod_split.o mod_tmpsav.o
spline.o: spline.f90
splitf.o: splitf.F90 mod_bxq.o mod_main.o mod_param1.o mod_regcm_param.o \
          mod_split.o
sponge.o: sponge.F90 mod_bdycod.o mod_regcm_param.o
spstep.o: spstep.F90 mod_bxq.o mod_main.o mod_regcm_param.o
stomat.o: stomat.f90 mod_constants.o mod_bats.o mod_ictp01.o mod_regcm_param.o
temp.o: temp.f90 mod_constants.o
ifeq ($(USECLM),1)
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_regcm_param.o mod_slice.o mod_trachem.o mod_interfaces.o mod_clm.o \
        $(CLMPTH)/clm_varsur.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_param1.o mod_param2.o mod_param3.o mod_pmoist.o mod_rad.o \
        mod_regcm_param.o mod_slice.o mod_trachem.o mod_interfaces.o
endif
tgrund.o: tgrund.F90 mod_constants.o mod_bats.o mod_param1.o mod_regcm_param.o
tpfc.o: tpfc.f90 mod_constants.o
tracbud.o: tracbud.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
           mod_mainchem.o mod_param1.o mod_param2.o mod_param3.o \
           mod_regcm_param.o mod_trachem.o
tracdiag.o: tracdiag.F90 mod_constants.o mod_diagnosis.o mod_main.o \
            mod_mainchem.o mod_param1.o mod_param3.o mod_regcm_param.o
tractend2.o: tractend2.f90 mod_constants.o mod_bats.o mod_cvaria.o mod_date.o \
             mod_dust.o mod_main.o mod_mainchem.o mod_message.o mod_param1.o \
             mod_param3.o mod_pbldim.o mod_pmoist.o mod_rad.o \
             mod_regcm_param.o mod_trachem.o
trcab.o: trcab.f90 mod_constants.o mod_regcm_param.o
trcabn.o: trcabn.f90 mod_constants.o mod_regcm_param.o
trcems.o: trcems.f90 mod_constants.o mod_regcm_param.o
trcmix.o: trcmix.f90 mod_regcm_param.o mod_tracer.o
trcplk.o: trcplk.f90 mod_constants.o mod_regcm_param.o
trcpth.o: trcpth.f90 mod_constants.o mod_regcm_param.o
tseaice.o: tseaice.f90 mod_constants.o mod_bats.o mod_param1.o mod_regcm_param.o
tstep.o: tstep.f90 mod_date.o mod_param1.o mod_regcm_param.o
vadv.o: vadv.f90 mod_cvaria.o mod_main.o mod_param3.o mod_regcm_param.o
vcheke.o: vcheke.f90
vcheki.o: vcheki.f90
vchekt.o: vchekt.f90
vcover.o: vcover.f90 mod_bats.o mod_regcm_param.o
vecbats.o: vecbats.F90 mod_date.o mod_param1.o mod_param2.o mod_interfaces.o
vmodes.o: vmodes.F90 mod_message.o mod_param1.o mod_regcm_param.o mod_split.o
vmultm.o: vmultm.f90
vnorml.o: vnorml.f90
vorder.o: vorder.f90
vtlaps.o: vtlaps.f90 mod_constants.o
water.o: water.f90 mod_constants.o mod_bats.o mod_param1.o mod_regcm_param.o
wheneq.o: wheneq.f90
whenfgt.o: whenfgt.f90
whenflt.o: whenflt.f90
whenne.o: whenne.f90
zengocn.o: zengocn.f90 mod_constants.o
ifeq ($(USECLM),1)
zengocndrv.o: zengocndrv.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_slice.o \
              mod_interfaces.o $(CLMPTH)/clm_varsur.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
zengocndrv.o: zengocndrv.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
              mod_param1.o mod_param2.o mod_param3.o mod_pbldim.o mod_slice.o \
              mod_interfaces.o mod_regcm_param.o
endif
zenith.o: zenith.f90
ifeq ($(USECLM),1)
zenitm.o: zenitm.F90 mod_regcm_param.o mod_clm.o $(CLMPTH)/shr_orb_mod.o \
          $(CLMPTH)/clm_time_manager.o mod_date.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
else
zenitm.o: zenitm.F90 mod_constants.o mod_date.o mod_main.o mod_regcm_param.o
endif
initclm.o: initclm.F90 $(CLMPTH)/initializeMod.o $(CLMPTH)/shr_orb_mod.o \
           $(CLMPTH)/shr_kind_mod.o $(CLMPTH)/clm_varpar.o \
           $(CLMPTH)/clm_varsur.o $(CLMPTH)/atmdrvMod.o $(CLMPTH)/spmdMod.o \
           $(CLMPTH)/program_offMod.o $(CLMPTH)/clm_comp.o \
           $(CLMPTH)/clmtype.o $(CLMPTH)/perf_mod.o \
           mod_constants.o mod_regcm_param.o mod_param1.o \
           mod_date.o mod_param2.o mod_param3.o mod_clm.o mod_bats.o \
           mod_mppio.o mod_main.o mod_pbldim.o mod_slice.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
interfclm.o: interfclm.F90 $(CLMPTH)/clm_varsur.o $(CLMPTH)/clmtype.o \
             mod_clm.o mod_regcm_param.o mod_date.o mod_main.o mod_param1.o \
             mod_param2.o mod_param3.o mod_slice.o mod_pbldim.o mod_bats.o \
             mod_constants.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
mtrxclm.o: mtrxclm.F90 $(CLMPTH)/atmdrvMod.o $(CLMPTH)/clm_comp.o \
           $(CLMPTH)/shr_kind_mod.o mod_message.o mod_clm.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
solar1clm.o: solar1clm.f90 $(CLMPTH)/shr_orb_mod.o $(CLMPTH)/clm_varsur.o \
             $(CLMPTH)/shr_kind_mod.o mod_clm.o mod_message.o mod_date.o \
             mod_constants.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
