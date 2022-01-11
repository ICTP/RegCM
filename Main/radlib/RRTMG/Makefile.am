#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegC.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later versio.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOS.  See the
#    GNU General Public License for more detail.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegC.  If not, see <http://ww.gn.org/licenses/.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = librrtmg.a
librrtmg_a_SOURCES = \
parkind.F90    rrlw_kg06.F90  rrlw_kg16.F90   rrsw_kg17.F90  rrsw_kg27.F90  \
parrrsw.F90    rrlw_kg07.F90  rrlw_ncpar.F90  rrsw_kg18.F90  rrsw_kg28.F90  \
parrrtm.F90    rrlw_kg08.F90  rrlw_ref.F90    rrsw_kg19.F90  rrsw_kg29.F90  \
rrlw_cld.F90   rrlw_kg09.F90  rrlw_tbl.F90    rrsw_kg20.F90  rrsw_ncpar.F90 \
rrlw_con.F90   rrlw_kg10.F90  rrlw_vsn.F90    rrsw_kg21.F90  rrsw_ref.F90   \
rrlw_kg01.F90  rrlw_kg11.F90  rrlw_wvn.F90    rrsw_kg22.F90  rrsw_tbl.F90   \
rrlw_kg02.F90  rrlw_kg12.F90  rrsw_aer.F90    rrsw_kg23.F90  rrsw_vsn.F90   \
rrlw_kg03.F90  rrlw_kg13.F90  rrsw_cld.F90    rrsw_kg24.F90  rrsw_wvn.F90   \
rrlw_kg04.F90  rrlw_kg14.F90  rrsw_con.F90    rrsw_kg25.F90  rrlw_kg05.F90  \
rrlw_kg15.F90  rrsw_kg16.F90  rrsw_kg26.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.F90

clean-local:
	rm -f *.mod *__genmod.F90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on i.
parkind.o : parkind.F90 mod_intkinds.o mod_realkinds.o
parrrsw.o : parrrsw.F90 parkind.o
parrrtm.o : parrrtm.F90 parkind.o
rrlw_cld.o : rrlw_cld.F90 parkind.o
rrlw_con.o : rrlw_con.F90 parkind.o
rrlw_kg01.o : rrlw_kg01.F90 parkind.o
rrlw_kg02.o : rrlw_kg02.F90 parkind.o
rrlw_kg03.o : rrlw_kg03.F90 parkind.o
rrlw_kg04.o : rrlw_kg04.F90 parkind.o
rrlw_kg05.o : rrlw_kg05.F90 parkind.o
rrlw_kg06.o : rrlw_kg06.F90 parkind.o
rrlw_kg07.o : rrlw_kg07.F90 parkind.o
rrlw_kg08.o : rrlw_kg08.F90 parkind.o
rrlw_kg09.o : rrlw_kg09.F90 parkind.o
rrlw_kg10.o : rrlw_kg10.F90 parkind.o
rrlw_kg11.o : rrlw_kg11.F90 parkind.o
rrlw_kg12.o : rrlw_kg12.F90 parkind.o
rrlw_kg13.o : rrlw_kg13.F90 parkind.o
rrlw_kg14.o : rrlw_kg14.F90 parkind.o
rrlw_kg15.o : rrlw_kg15.F90 parkind.o
rrlw_kg16.o : rrlw_kg16.F90 parkind.o
rrlw_ncpar.o : rrlw_ncpar.F90 parkind.o
rrlw_ref.o : rrlw_ref.F90 parkind.o
rrlw_tbl.o : rrlw_tbl.F90 parkind.o
rrlw_vsn.o : rrlw_vsn.F90
rrlw_wvn.o : rrlw_wvn.F90 parkind.o parrrtm.o
rrsw_aer.o : rrsw_aer.F90 parkind.o parrrsw.o
rrsw_cld.o : rrsw_cld.F90 parkind.o
rrsw_con.o : rrsw_con.F90 parkind.o
rrsw_kg16.o : rrsw_kg16.F90 parkind.o parrrsw.o
rrsw_kg17.o : rrsw_kg17.F90 parkind.o parrrsw.o
rrsw_kg18.o : rrsw_kg18.F90 parkind.o parrrsw.o
rrsw_kg19.o : rrsw_kg19.F90 parkind.o parrrsw.o
rrsw_kg20.o : rrsw_kg20.F90 parkind.o parrrsw.o
rrsw_kg21.o : rrsw_kg21.F90 parkind.o parrrsw.o
rrsw_kg22.o : rrsw_kg22.F90 parkind.o parrrsw.o
rrsw_kg23.o : rrsw_kg23.F90 parkind.o parrrsw.o
rrsw_kg24.o : rrsw_kg24.F90 parkind.o parrrsw.o
rrsw_kg25.o : rrsw_kg25.F90 parkind.o parrrsw.o
rrsw_kg26.o : rrsw_kg26.F90 parkind.o parrrsw.o
rrsw_kg27.o : rrsw_kg27.F90 parkind.o parrrsw.o
rrsw_kg28.o : rrsw_kg28.F90 parkind.o parrrsw.o
rrsw_kg29.o : rrsw_kg29.F90 parkind.o parrrsw.o
rrsw_ncpar.o : rrsw_ncpar.F90 parkind.o
rrsw_ref.o : rrsw_ref.F90 parkind.o
rrsw_tbl.o : rrsw_tbl.F90 parkind.o
rrsw_vsn.o : rrsw_vsn.F90
rrsw_wvn.o : rrsw_wvn.F90 parkind.o parrrsw.o
