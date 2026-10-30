!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Unit test of Share/mod_ncstream.F90 without MPI: an output stream with
! static, time dependent and time averaged variables is written in the
! current directory and read back with the netCDF library (data, time,
! time bounds, coordinates, attributes). A file with a single record is
! also written.

program test_ncstream
  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_dynparam
  use mod_ncstream
  use share_test_utils
  implicit none
  character(len=*), parameter :: fname = 'test_ncstream.nc'
  character(len=*), parameter :: fname1 = 'test_ncstream1.nc'
  integer(ik4), parameter :: nj = 16, ni = 12, nk = 4
  type(nc_output_stream) :: ncout
  type(ncoutstream_params) :: opar
  type(ncvariable2d_mixed) :: vtopo
  type(ncvariable3d_mixed) :: vta
  type(ncvariable2d_mixed) :: vpr
  type(ncvariable2d_integer) :: vmask
  real(rkx), pointer, contiguous, dimension(:) :: sigma => null()
  real(rkx), target, dimension(nj,ni) :: topo, pr, r2
  real(rkx), target, dimension(nj,ni,nk) :: ta, r3
  integer(ik4), target, dimension(nj,ni) :: mask, im
  integer(ik4) :: ifill
  real(rk8), dimension(2,2) :: tb
  real(rk8), dimension(2) :: tv
  type(rcm_time_and_date) :: idate, d0, d1
  integer(ik4) :: ncid, ivar, istat, i, j, k, n
  character(len=64) :: ca

  jx = nj
  iy = ni
  kz = nk
  ds = 50.0_rkx
  nsg = 1
  domname = 'test'
  calendar = 'gregorian'
  ical = gregorian
  iproj = 'LAMCON'
  idynamic = 1
  xcone = 0.71_rkx
  clat = 43.0_rkx
  clon = 15.0_rkx
  cntrj = real(nj+1,rkx)*0.5_rkx
  cntri = real(ni+1,rkx)*0.5_rkx
  ptop = 5.0_rkx
  truelatl = 30.0_rkx
  truelath = 60.0_rkx
  allocate(sigma(nk))
  sigma = [0.2_rkx, 0.5_rkx, 0.8_rkx, 0.95_rkx]
  do i = 1, ni
    do j = 1, nj
      topo(j,i) = real(j+100*i,rkx)
      do k = 1, nk
        ta(j,i,k) = real(j+100*i+10000*k,rkx)
      end do
    end do
  end do

  vtopo%vname = 'topo'
  vtopo%vunit = 'm'
  vtopo%long_name = 'Surface height'
  vtopo%standard_name = 'surface_altitude'
  vtopo%lrecords = .false.
  vtopo%rval => topo
  vta%vname = 'ta'
  vta%vunit = 'K'
  vta%long_name = 'Temperature'
  vta%standard_name = 'air_temperature'
  vta%lrecords = .true.
  vta%rval => ta
  vpr%vname = 'pr'
  vpr%vunit = 'kg m-2 s-1'
  vpr%long_name = 'Precipitation'
  vpr%standard_name = 'precipitation_flux'
  vpr%lrecords = .true.
  vpr%lfillvalue = .true.
  vpr%cell_method = 'time: mean'
  vpr%rval => pr

  vmask%vname = 'mask'
  vmask%vunit = '1'
  vmask%long_name = 'Land mask'
  vmask%standard_name = 'land_binary_mask'
  vmask%lrecords = .false.
  vmask%lfillvalue = .true.
  mask = 1
  mask(1,:) = 0
  vmask%ival => mask

  opar%fname = fname
  opar%pname = 'test_ncstream'
  opar%zero_date = 1979022200
  opar%l_bound = .true.
  call outstream_setup(ncout,opar)
  call outstream_addvar(ncout,vtopo)
  call outstream_addvar(ncout,vta)
  call outstream_addvar(ncout,vpr)
  call outstream_addvar(ncout,vmask)
  call outstream_addatt(ncout,ncattribute_string('lake_fudging','No'))
  call outstream_addatt(ncout,ncattribute_logical('flag',.true.))
  call outstream_addatt(ncout,ncattribute_integer('count',7))
  call outstream_addvaratt(ncout,vtopo,ncattribute_string('notes','test'))
  call outstream_enable(ncout,sigma)
  call outstream_writevar(ncout,vtopo)
  call outstream_writevar(ncout,vmask)
  do n = 1, 3
    idate = 1979022200
    idate = idate + rcm_time_interval(6*n,uhrs)
    call outstream_addrec(ncout,idate)
    ta = ta + 1.0_rkx
    pr = real(n,rkx)
    call outstream_writevar(ncout,vta)
    call outstream_writevar(ncout,vpr)
  end do
  call outstream_dispose(ncout)
  call outstream_async_finalize( )

  ! Check the file with the netCDF library.
  d0 = 1950010100
  istat = nf90_open(fname,nf90_nowrite,ncid)
  istat = nf90_inq_varid(ncid,'time',ivar)
  istat = nf90_get_var(ncid,ivar,tv,[2],[2])
  d1 = 1979022212
  call check(abs(tv(1)-hourdiff(d1,d0)) < 1.0e-9_rk8 &
    .and. abs(tv(2)-tv(1)-6.0_rk8) < 1.0e-9_rk8, &
    'time in hours since 1950-01-01')
  istat = nf90_get_att(ncid,ivar,'units',ca)
  call check(ca == 'hours since 1950-01-01 00:00:00 UTC','time units')
  istat = nf90_inq_varid(ncid,'time_bnds',ivar)
  call check(istat == nf90_noerr,'time bounds for a time mean')
  istat = nf90_get_var(ncid,ivar,tb,[1,1],[2,2])
  d1 = 1979022200
  call check(abs(tb(1,1)-hourdiff(d1,d0)) < 1.0e-9_rk8 &
    .and. abs(tb(2,1)-tb(1,2)) < 1.0e-9_rk8 .and. &
    abs(tb(2,2)-tb(1,2)-6.0_rk8) < 1.0e-9_rk8,'time bounds values')
  istat = nf90_inq_varid(ncid,'pr',ivar)
  istat = nf90_get_att(ncid,ivar,'cell_methods',ca)
  call check(ca == 'time: mean','cell_methods')
  istat = nf90_get_att(ncid,ivar,'grid_mapping',ca)
  call check(ca == 'crs','grid_mapping')
  istat = nf90_inquire_attribute(ncid,ivar,'_FillValue')
  call check(istat == nf90_noerr,'_FillValue')
  istat = nf90_inq_varid(ncid,'crs',ivar)
  istat = nf90_get_att(ncid,ivar,'grid_mapping_name',ca)
  call check(ca == 'lambert_conformal_conic','grid mapping of LAMCON')
  istat = nf90_inq_varid(ncid,'jx',ivar)
  istat = nf90_get_var(ncid,ivar,tv,[1],[2])
  call check(abs(tv(1)+7.5_rk8*ds*1000.0_rk8) < 1.0e-6_rk8 .and. &
             abs(tv(2)-tv(1)-ds*1000.0_rk8) < 1.0e-6_rk8, &
    'x coordinate in m, centred')
  istat = nf90_inq_varid(ncid,'ptop',ivar)
  istat = nf90_get_var(ncid,ivar,tv(1))
  call check(abs(tv(1)-50.0_rk8) < 1.0e-9_rk8,'ptop in hPa')
  istat = nf90_get_att(ncid,nf90_global,'flag',ca)
  call check(ca == 'Yes','logical global attribute')
  istat = nf90_inq_varid(ncid,'topo',ivar)
  istat = nf90_get_att(ncid,ivar,'notes',ca)
  call check(ca == 'test','variable attribute')
  istat = nf90_get_var(ncid,ivar,r2)
  call check(all(r2 == topo),'static 2-D variable')
  istat = nf90_inq_varid(ncid,'ta',ivar)
  istat = nf90_get_var(ncid,ivar,r3,[1,1,1,2],[nj,ni,nk,1])
  call check(all(abs(r3-(ta-1.0_rkx)) < 1.0e-3_rkx), &
    'record of a 3-D variable')
  istat = nf90_inq_varid(ncid,'pr',ivar)
  istat = nf90_get_var(ncid,ivar,r2,[1,1,3],[nj,ni,1])
  call check(all(r2 == 3.0_rkx),'last record of a 2-D variable')
  istat = nf90_inq_varid(ncid,'mask',ivar)
  istat = nf90_get_var(ncid,ivar,im)
  call check(all(im == mask),'integer variable')
  istat = nf90_get_att(ncid,ivar,'_FillValue',ifill)
  call check(istat == nf90_noerr .and. ifill == -9999, &
    'integer _FillValue')
  istat = nf90_close(ncid)

  ! A file with a single record.
  opar%fname = fname1
  call outstream_setup(ncout,opar)
  call outstream_addvar(ncout,vta)
  call outstream_enable(ncout,sigma)
  idate = 1979022206
  call outstream_addrec(ncout,idate)
  call outstream_writevar(ncout,vta)
  call outstream_dispose(ncout)
  call outstream_async_finalize( )
  istat = nf90_open(fname1,nf90_nowrite,ncid)
  istat = nf90_inq_dimid(ncid,'time',ivar)
  istat = nf90_inquire_dimension(ncid,ivar,len=n)
  call check(n == 1,'one record file')
  istat = nf90_close(ncid)

  call test_end('test_ncstream')
end program test_ncstream
