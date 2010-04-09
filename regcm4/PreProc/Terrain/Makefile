# Makefile for RegCM4
# for Terrain code.
  
include ../../Makefile.inc


all: terrain

SHRMODS =  mod_preproc_param.o mod_regcm_param.o

MODOBJS =  mod_block.o     mod_icontr.o     mod_const.o      mod_maps.o   \
           mod_aa.o        mod_maps_s.o     mod_addstack2.o  mod_verts.o  \
           mod_addstack.o  mod_a.o

OBJS =     anal2.o         lndfudge.o       oned.o           smthtr.o     \
           bint.o          output.o         surf.o           fexist.o     \
           mappol.o        rdldtr.o         rdldtr_s.o       rdldtr_nc.o  \
           smther.o        rot2nrot.o       texfudge.o       interp.o     \
           mxmnll.o        rotmer.o         xyobsll.o        interp_s.o   \
           normer.o        setup.o          lakeadj.o        nrot2rot.o   \
           smth121.o       lambrt.o         rdldtr_s_nc.o    terrain.o    \
           header_terrain.o

mod_regcm_param.o: ../../Config/mod_regcm_param.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_preproc_param.o: ../../Config/mod_preproc_param.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

terrain: $(SHRMODS) $(MODOBJS) $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMODS) $(MODOBJS) $(OBJS) \
               $(LIB) $(NETCDFLIB)

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

header_terrain.o: header_terrain.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
  
clean:
	rm -f *.mod *.o terrain

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_aa.o: mod_aa.f90
mod_addstack2.o: mod_addstack2.f90 mod_regcm_param.o
mod_addstack.o: mod_addstack.f90
mod_a.o: mod_a.f90
mod_block.o: mod_block.f90
mod_const.o: mod_const.f90
mod_icontr.o: mod_icontr.f90
mod_maps.o: mod_maps.f90 mod_regcm_param.o mod_preproc_param.o
mod_maps_s.o: mod_maps_s.f90 mod_regcm_param.o mod_preproc_param.o
mod_verts.o: mod_verts.f90 mod_regcm_param.o
anal2.o: anal2.f90 mod_addstack.o mod_block.o mod_const.o
bint.o: bint.f90
fexist.o: fexist.f90
interp.o: interp.f90 mod_preproc_param.o mod_aa.o mod_block.o mod_const.o \
          mod_maps.o
interp_s.o: interp_s.f90 mod_preproc_param.o mod_aa.o mod_block.o mod_const.o \
            mod_maps_s.o
lakeadj.o: lakeadj.f90
lambrt.o: lambrt.f90
lndfudge.o: lndfudge.f90
mappol.o: mappol.f90
mxmnll.o: mxmnll.f90 mod_aa.o
normer.o: normer.f90
nrot2rot.o: nrot2rot.f90
oned.o: oned.f90
output.o: output.f90
rdldtr.o: rdldtr.f90 mod_preproc_param.o mod_a.o mod_aa.o mod_block.o \
          mod_const.o
rdldtr_nc.o: rdldtr_nc.f90 mod_preproc_param.o mod_a.o mod_aa.o mod_block.o \
             mod_const.o
rdldtr_s.o: rdldtr_s.f90 mod_preproc_param.o mod_a.o mod_aa.o mod_block.o \
            mod_const.o
rdldtr_s_nc.o: rdldtr_s_nc.f90 mod_preproc_param.o mod_a.o mod_aa.o \
               mod_block.o mod_const.o
rot2nrot.o: rot2nrot.f90
rotmer.o: rotmer.f90
setup.o: setup.f90 mod_const.o
smth121.o: smth121.f90
smther.o: smther.f90
smthtr.o: smthtr.f90
surf.o: surf.f90 mod_a.o
xyobsll.o: xyobsll.f90 mod_block.o mod_const.o
texfudge.o: texfudge.f90
terrain.o: terrain.f90 mod_regcm_param.o mod_preproc_param.o mod_aa.o \
           mod_addstack.o mod_addstack2.o mod_block.o mod_const.o \
           mod_icontr.o mod_maps.o mod_maps_s.o mod_verts.o
