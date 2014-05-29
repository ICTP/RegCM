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

noinst_LIBRARIES = librrtmgsw.a
librrtmgsw_a_SOURCES = \
parkind.f90    rrlw_kg06.f90  rrlw_kg16.f90   rrsw_kg17.f90  rrsw_kg27.f90  \
parrrsw.f90    rrlw_kg07.f90  rrlw_ncpar.f90  rrsw_kg18.f90  rrsw_kg28.f90  \
parrrtm.f90    rrlw_kg08.f90  rrlw_ref.f90    rrsw_kg19.f90  rrsw_kg29.f90  \
rrlw_cld.f90   rrlw_kg09.f90  rrlw_tbl.f90    rrsw_kg20.f90  rrsw_ncpar.f90 \
rrlw_con.f90   rrlw_kg10.f90  rrlw_vsn.f90    rrsw_kg21.f90  rrsw_ref.f90   \
rrlw_kg01.f90  rrlw_kg11.f90  rrlw_wvn.f90    rrsw_kg22.f90  rrsw_tbl.f90   \
rrlw_kg02.f90  rrlw_kg12.f90  rrsw_aer.f90    rrsw_kg23.f90  rrsw_vsn.f90   \
rrlw_kg03.f90  rrlw_kg13.f90  rrsw_cld.f90    rrsw_kg24.f90  rrsw_wvn.f90   \
rrlw_kg04.f90  rrlw_kg14.f90  rrsw_con.f90    rrsw_kg25.f90  rrlw_kg05.f90  \
rrlw_kg15.f90  rrsw_kg16.f90  rrsw_kg26.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.f90

clean-local:
	rm -f *.mod *__genmod.f90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on i.
parkind.o : parkind.f90
parrrsw.o : parrrsw.f90 parkind.o
parrrtm.o : parrrtm.f90 parkind.o
rrlw_cld.o : rrlw_cld.f90 parkind.o
rrlw_con.o : rrlw_con.f90 parkind.o
rrlw_kg01.o : rrlw_kg01.f90 parkind.o
rrlw_kg02.o : rrlw_kg02.f90 parkind.o
rrlw_kg03.o : rrlw_kg03.f90 parkind.o
rrlw_kg04.o : rrlw_kg04.f90 parkind.o
rrlw_kg05.o : rrlw_kg05.f90 parkind.o
rrlw_kg06.o : rrlw_kg06.f90 parkind.o
rrlw_kg07.o : rrlw_kg07.f90 parkind.o
rrlw_kg08.o : rrlw_kg08.f90 parkind.o
rrlw_kg09.o : rrlw_kg09.f90 parkind.o
rrlw_kg10.o : rrlw_kg10.f90 parkind.o
rrlw_kg11.o : rrlw_kg11.f90 parkind.o
rrlw_kg12.o : rrlw_kg12.f90 parkind.o
rrlw_kg13.o : rrlw_kg13.f90 parkind.o
rrlw_kg14.o : rrlw_kg14.f90 parkind.o
rrlw_kg15.o : rrlw_kg15.f90 parkind.o
rrlw_kg16.o : rrlw_kg16.f90 parkind.o
rrlw_ncpar.o : rrlw_ncpar.f90 parkind.o
rrlw_ref.o : rrlw_ref.f90 parkind.o
rrlw_tbl.o : rrlw_tbl.f90 parkind.o
rrlw_vsn.o : rrlw_vsn.f90
rrlw_wvn.o : rrlw_wvn.f90 parkind.o parrrtm.o
rrsw_aer.o : rrsw_aer.f90 parkind.o parrrsw.o
rrsw_cld.o : rrsw_cld.f90 parkind.o
rrsw_con.o : rrsw_con.f90 parkind.o
rrsw_kg16.o : rrsw_kg16.f90 parkind.o parrrsw.o
rrsw_kg17.o : rrsw_kg17.f90 parkind.o parrrsw.o
rrsw_kg18.o : rrsw_kg18.f90 parkind.o parrrsw.o
rrsw_kg19.o : rrsw_kg19.f90 parkind.o parrrsw.o
rrsw_kg20.o : rrsw_kg20.f90 parkind.o parrrsw.o
rrsw_kg21.o : rrsw_kg21.f90 parkind.o parrrsw.o
rrsw_kg22.o : rrsw_kg22.f90 parkind.o parrrsw.o
rrsw_kg23.o : rrsw_kg23.f90 parkind.o parrrsw.o
rrsw_kg24.o : rrsw_kg24.f90 parkind.o parrrsw.o
rrsw_kg25.o : rrsw_kg25.f90 parkind.o parrrsw.o
rrsw_kg26.o : rrsw_kg26.f90 parkind.o parrrsw.o
rrsw_kg27.o : rrsw_kg27.f90 parkind.o parrrsw.o
rrsw_kg28.o : rrsw_kg28.f90 parkind.o parrrsw.o
rrsw_kg29.o : rrsw_kg29.f90 parkind.o parrrsw.o
rrsw_ncpar.o : rrsw_ncpar.f90 parkind.o
rrsw_ref.o : rrsw_ref.f90 parkind.o
rrsw_tbl.o : rrsw_tbl.f90 parkind.o
rrsw_vsn.o : rrsw_vsn.f90
rrsw_wvn.o : rrsw_wvn.f90 parkind.o parrrsw.o
