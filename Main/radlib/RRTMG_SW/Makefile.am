
noinst_LIBRARIES = librrtmgsw.a
librrtmgsw_a_SOURCES =    parkind.f90  parrrsw.f90  rrsw_cld.f90   \
  rrsw_kg17.f90  rrsw_kg20.f90  rrsw_kg23.f90  rrsw_kg26.f90 \
  rrsw_kg29.f90   rrsw_tbl.f90  parrrsw.f90   rrsw_con.f90  \
  rrsw_kg18.f90  rrsw_kg21.f90  rrsw_kg24.f90  rrsw_kg27.f90 \
  rrsw_ncpar.f90  rrsw_vsn.f90 rrsw_aer.f90  rrsw_kg16.f90  \
  rrsw_kg19.f90  rrsw_kg22.f90  rrsw_kg25.f90  rrsw_kg28.f90 \
  rrsw_ref.f90    rrsw_wvn.f90


%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.


parkind.o parkind.mod: parkind.f90 
parrrsw.o parrrsw.mod: parrrsw.f90
rrsw_cld.o  rrsw_cld.mod:  rrsw_cld.f90
rrsw_kg17.o rrsw_kg17.mod: rrsw_kg17.f90 parrrsw.o
rrsw_kg20.o rrsw_kg20.mod: rrsw_kg20.f90 
rrsw_kg23.o rrsw_kg23.mod :rrsw_kg23.f90
rrsw_kg26.o rrsw_kg26.mod: rrsw_kg26.f90 
rrsw_kg29.o  rrsw_kg29.mod:rrsw_kg29.f90 
rrsw_tbl.o  rrsw_tbl.mod: rrsw_tbl.f90 
rrsw_con.o rrsw_con.mod: rrsw_con.f90
rrsw_kg18.o  rrsw_kg18.mod: rrsw_kg18.f90
rrsw_kg21.o  rrsw_kg21.mod: rrsw_kg21.f90 
rrsw_kg24.o  rrsw_kg24.mod: rrsw_kg24.f90
rrsw_kg27.o  rrsw_kg27.mod: rrsw_kg27.f90 
rrsw_ncpar.o rrsw_ncpar.mod: rrsw_ncpar.f90 
rrsw_vsn.o  rrsw_vsn.mod: rrsw_vsn.f90
rrsw_aer.o rrsw_aer.mod: rrsw_aer.f90
rrsw_kg16.o rrsw_kg16.mod: rrsw_kg16.f90
rrsw_kg19.o rrsw_kg19.mod: rrsw_kg19.f90
rrsw_kg22.o rrsw_kg22.mod: rrsw_kg22.f90 
rrsw_kg25.o  rrsw_kg25.mod: rrsw_kg25.f90
rrsw_kg28.o rrsw_kg28.mod: rrsw_kg28.f90
rrsw_ref.o  rrsw_ref.mod: rrsw_ref.f90 
rrsw_wvn.o  rrsw_wvn.mod: rrsw_wvn.f90 

