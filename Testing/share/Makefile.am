#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#    
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Unit tests of the Share library, built and run only by "make check".
# Each program prints PASS or FAIL for every check and exits with a non
# zero status if any check failed.

vpath %.o $(REGCM_ROOT)/Share:.

include $(REGCM_ROOT)/makeinc

check_PROGRAMS = test_date test_spbarcoord test_projections test_hgt \
                 test_vertint test_humid test_interp test_spline \
                 test_memutil test_kdtree2 test_kdinterp test_nchelper \
                 test_ncstream test_async_netcdf test_domain test_stdatm \
                 test_wavelet test_ensemble test_sunorbit test_intldtr \
                 test_rdldtr
test_date_SOURCES = share_test_utils.F90 test_date.F90
test_spbarcoord_SOURCES = share_test_utils.F90 test_spbarcoord.F90
test_projections_SOURCES = share_test_utils.F90 test_projections.F90
test_hgt_SOURCES = share_test_utils.F90 test_hgt.F90
test_vertint_SOURCES = share_test_utils.F90 test_vertint.F90
test_humid_SOURCES = share_test_utils.F90 test_humid.F90
test_interp_SOURCES = share_test_utils.F90 test_interp.F90
test_spline_SOURCES = share_test_utils.F90 test_spline.F90
test_memutil_SOURCES = share_test_utils.F90 test_memutil.F90
test_kdtree2_SOURCES = share_test_utils.F90 test_kdtree2.F90
test_kdinterp_SOURCES = share_test_utils.F90 test_kdinterp.F90
test_nchelper_SOURCES = share_test_utils.F90 test_nchelper.F90
test_ncstream_SOURCES = share_test_utils.F90 test_ncstream.F90
test_async_netcdf_SOURCES = share_test_utils.F90 test_async_netcdf.F90
test_domain_SOURCES = share_test_utils.F90 test_domain.F90
test_stdatm_SOURCES = share_test_utils.F90 test_stdatm.F90
test_wavelet_SOURCES = share_test_utils.F90 test_wavelet.F90
test_ensemble_SOURCES = share_test_utils.F90 test_ensemble.F90
test_sunorbit_SOURCES = share_test_utils.F90 test_sunorbit.F90
test_intldtr_SOURCES = share_test_utils.F90 test_intldtr.F90
test_rdldtr_SOURCES = share_test_utils.F90 test_rdldtr.F90
LDADD = $(OMPFLAGS)

FCLD = $(MPIFC)

TESTS = $(check_PROGRAMS)

# A small memory cap: with --enable-async-netcdf the writes of
# test_async_netcdf go through the worker thread and wait for space.
AM_TESTS_ENVIRONMENT = RCM_ASYNC_OUTPUT_GB=0.00005; export RCM_ASYNC_OUTPUT_GB;

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt test_nchelper.nc test_ncstream*.nc \
	      test_async_netcdf.nc test_domain*.nc test_rdldtr.nc

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt test_nchelper.nc test_ncstream*.nc \
	      test_async_netcdf.nc test_domain*.nc test_rdldtr.nc

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
share_test_utils.o : share_test_utils.F90 mod_intkinds.o mod_realkinds.o
test_date.o : test_date.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_date.o
test_spbarcoord.o : test_spbarcoord.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_spbarcoord.o
test_projections.o : test_projections.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_projections.o mod_memutil.o
test_hgt.o : test_hgt.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_constants.o mod_hgt.o
test_vertint.o : test_vertint.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_constants.o mod_vertint.o mod_memutil.o
test_humid.o : test_humid.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_constants.o mod_humid.o
test_interp.o : test_interp.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_interp.o
test_spline.o : test_spline.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_spline.o
test_memutil.o : test_memutil.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_memutil.o mod_date.o
test_kdtree2.o : test_kdtree2.F90 share_test_utils.o mod_kdtree2.o
test_kdinterp.o : test_kdinterp.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_kdinterp.o
test_nchelper.o : test_nchelper.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_nchelper.o
test_ncstream.o : test_ncstream.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_date.o mod_dynparam.o mod_ncstream.o
test_async_netcdf.o : test_async_netcdf.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_async_netcdf.o
test_domain.o : test_domain.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_ncstream.o mod_nchelper.o mod_domain.o
test_stdatm.o : test_stdatm.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_stdatm.o
test_wavelet.o : test_wavelet.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_wavelet.o
test_ensemble.o : test_ensemble.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_ensemble.o
test_sunorbit.o : test_sunorbit.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_sunorbit.o
test_intldtr.o : test_intldtr.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_intldtr.o
test_rdldtr.o : test_rdldtr.F90 share_test_utils.o mod_intkinds.o mod_realkinds.o mod_memutil.o mod_kdinterp.o mod_rdldtr.o
