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

subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort

program sigma2z
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_message
  use mod_vertint
  use mod_nchelper
  use mod_dynparam, only : iomode
#ifdef NETCDF4_HDF5
  use mod_dynparam, only : ncfilter, ncfilter_nparams, ncfilter_params
#endif
  use mod_hgt
  use mod_humid
  use mod_stdio
  use mod_memutil
  use mod_sigma, only : init_sigma, half_sigma_coordinate
  use netcdf

  implicit none

  character(256) :: prgname, ncsfile, ncpfile
  character(128) :: attname, dimname, varname
  integer(ik4) :: numarg, istatus, ncid, ncout

  integer(ik4), allocatable, dimension(:) :: dimids, dimlen
  real(rk4), allocatable, dimension(:) :: sigma, ak, bk
  real(rk4), allocatable, dimension(:,:,:) :: tazvar, hzvar, qazvar
  real(rk4), allocatable, dimension(:,:,:) :: pp, press, pai
  real(rk4), allocatable, save, dimension(:,:,:) :: zvar, xvar
  real(rk4), allocatable, dimension(:,:) :: ps, topo, mslpr, ps0
  real(rk4), allocatable, dimension(:) :: avar
  character, allocatable, dimension(:) :: txtvar
  real(rk4), allocatable, dimension(:) :: azvar
  real(rkx), allocatable, dimension(:) :: times
  real(rkx), allocatable, dimension(:) :: sigfix
  logical, allocatable, dimension(:) :: lkvarflag, ltvarflag, lchnameflag
  integer(ik4), allocatable, dimension(:) :: varsize
  integer(ik4), allocatable, dimension(:) :: intscheme
  integer(ik4), allocatable, dimension(:) :: nvdims
  integer(ik4), allocatable, dimension(:,:) :: dimsize
  integer(ik4), allocatable, dimension(:) :: istart, icount
  integer(ik4), allocatable, dimension(:) :: invarid
  integer(ik4), allocatable, dimension(:) :: outvarid
  integer(ik4) :: ndims, nvars, natts, udimid, nvatts
  integer(ik4) :: ivarid, idimid, xtype
  integer(ik4) :: jxdimid, iydimid, kzdimid, itdimid, itvarid, ikvarid
  integer(ik4) :: ipsvarid, ishvarid, ppvarid, ip0varid
  integer(ik4) :: avarid, bvarid, paivarid
  integer(ik4) :: jx, iy, kz, nt
  real(rkx) :: ptop
  integer(ik4), dimension(4) :: tdimids
  integer(ik4), dimension(3) :: psdimids
  integer(ik4) :: i, j, k, it, iv, iid1, iid2, ii, i3d, p3d, ich
  integer(ik4) :: tvarid, qvarid, irhvar, imslzvar, ircm_map
  logical :: has_t, has_q, has_rh, has_sph, is_icbc
  logical :: make_rh, make_mslp
  integer(ik4) :: n3d, iz3d, iodyn, nz, ipunit, iresult
  real(rk4), allocatable, dimension(:) :: zlevs

  data has_t /.false./
  data has_q /.false./
  data has_rh /.false./
  data make_rh /.false./
  data make_mslp /.false./

!$OMP THREADPRIVATE(xvar,zvar)

#ifdef PNETCDF
  iomode = ior(nf90_clobber, nf90_64bit_offset)
#endif

  call get_command_argument(0,value=prgname)
  numarg = command_argument_count()
  if (numarg < 1) then
    write (6,*) 'Not enough arguments.'
    write (6,*) ' '
    write (6,*) 'Usage : ', trim(prgname), '[namelist_file.in] Rcmfile.nc'
    write (6,*) ' '
    stop
  end if

  call get_command_argument(1,value=ncsfile)
  istatus = nf90_open(ncsfile, nf90_nowrite, ncid)
  if ( istatus /= nf90_noerr ) then
    ! Assume we have been provided a namelist with pressure levels.
    open(newunit=ipunit, file=ncsfile, status='old', &
         action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(ncsfile)
      stop
    end if
#ifdef NETCDF4_HDF5
    call get_ncfilter(ipunit)
#endif
    call get_nz(ipunit)
    allocate(zlevs(nz))
    call get_zlevs(ipunit)
    call get_command_argument(2,value=ncsfile)
    istatus = nf90_open(ncsfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
           'Error Opening Input file '//trim(ncsfile))
  else
    call checkncerr(istatus,__FILE__,__LINE__, &
           'Error Opening Input file '//trim(ncsfile))
    nz = 14
    allocate(zlevs(nz))
    zlevs = [ 20.0,50.0,80.0,100.0,150.0,200.0,500.0,750.0,1000.0, &
              1500.0,2000.0,5000.0,7000.0,10000.0 ]
  end if

  iid1 = scan(ncsfile, '/', .true.)
  iid2 = scan(ncsfile, '.', .true.)
  ncpfile = trim(ncsfile(iid1+1:iid2-1))//'_hgt.nc'

  jxdimid = -1
  iydimid = -1
  kzdimid = -1
  itdimid = -1
  itvarid = -1
  ikvarid = -1
  avarid = -1
  bvarid = -1
  iy = -1
  jx = -1
  istatus = nf90_create(ncpfile, iomode, ncout)
  call checkncerr(istatus,__FILE__,__LINE__, &
          'Error Opening Output file '//trim(ncpfile))

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  call checkncerr(istatus,__FILE__,__LINE__, &
          'Error Reading Netcdf file '//trim(ncsfile))

  do i = 1, natts
    istatus = nf90_inq_attname(ncid, nf90_global, i, attname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read global attribute')
    istatus = nf90_copy_att(ncid, nf90_global, attname, ncout, nf90_global)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Error copy attribute '//trim(attname))
  end do

  istatus = nf90_get_att(ncid, nf90_global, 'dynamical_core', iodyn)
  if ( istatus /= nf90_noerr ) then
    iodyn = 1
  end if

  allocate(dimlen(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimlen')

  kz = 0
  nt = 0
  do i = 1, ndims
    istatus = nf90_inquire_dimension(ncid, i, dimname, dimlen(i))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading dimension info')
    if (dimname == 'iy' .or. dimname == 'y') then
      iy = dimlen(i)
      iydimid = i
    else if (dimname == 'jx' .or. dimname == 'x') then
      jx = dimlen(i)
      jxdimid = i
    else if (dimname == 'kz' .or. dimname == 'z' .or. dimname == 'lev') then
      kz = dimlen(i)
      kzdimid = i
    end if
    if (dimname == 'time') then
      itdimid = i
      nt = dimlen(i)
      istatus = nf90_def_dim(ncout, dimname, nf90_unlimited, idimid)
    else if (dimname == 'kz' .or. dimname == 'z' .or. dimname == 'lev') then
      istatus = nf90_def_dim(ncout, 'zlev', nz, idimid)
    else
      istatus = nf90_def_dim(ncout, dimname, dimlen(i), idimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Error define dimension '//trim(dimname))
  end do

  i3d = jx*iy*kz
  p3d = jx*iy*nz

  if (kz == 0 .or. nt == 0) then
    write (6,*) 'Nothing to do: no vertical dimension or no time'
    istatus = nf90_close(ncout)
    istatus = nf90_close(ncid)
    stop
  end if

  allocate(dimids(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimids')
  allocate(lkvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'lkvarflag')
  allocate(ltvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'ltvarflag')
  allocate(lchnameflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'lchnameflag')
  allocate(varsize(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'varsize')
  allocate(invarid(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'invarid')
  allocate(outvarid(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'outvarid')
  allocate(intscheme(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'intscheme')
  allocate(nvdims(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'nvdims')
  allocate(dimsize(ndims,nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimsize')
  allocate(istart(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'istart')
  allocate(icount(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'icount')
  allocate(ps(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'ps')
  allocate(xvar(jx,iy,kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'xvar')
  allocate(tazvar(jx,iy,kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'tazvar')
  allocate(hzvar(jx,iy,kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'hzvar')
  allocate(zvar(jx,iy,nz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'zvar')
  allocate(times(nt), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'times')
  allocate(sigma(kz), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'sigma')
  allocate(topo(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'topo')

  if ( iodyn == 2 ) then
    allocate(ps0(jx,iy), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'ps0')
    allocate(pp(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'pp')
    allocate(press(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'press')
  else if ( iodyn == 3 ) then
    allocate(pai(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'pai')
    allocate(press(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'press')
    allocate(ak(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'ak')
    allocate(bk(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'bk')
  end if

  paivarid = -1
  ppvarid = -1
  tvarid = -1
  ip0varid = -1
  ircm_map = -1
  invarid = -1
  outvarid = -1

  has_sph = .false.
  is_icbc = .false.

  do i = 1, nvars
    lkvarflag(i) = .false.
    ltvarflag(i) = .false.
    lchnameflag(i) = .false.
    varsize(i) = 1
    intscheme(i) = 1
    dimsize(:,i) = 1
    istatus = nf90_inquire_variable(ncid,i,varname,xtype,nvdims(i), &
                                    dimids,nvatts)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire variable')
    if (varname == 'a') then
      avarid = i
      cycle
    else if (varname == 'b') then
      bvarid = i
      cycle
    else if (varname == 'time') then
      itvarid = i
    else if (varname == 'kz' .or. &
             varname == 'sigma' .or. &
             varname == 'lev') then
      varname = 'zlev'
      ikvarid = i
    else if (varname == 'ptop') then
      istatus = nf90_get_var(ncid,i,ptop)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read variable ptop')
      ptop = ptop * d_100
    else if (varname == 'crs' .or. varname == 'rcm_map' ) then
      ircm_map = i
    else if (varname == 'ps') then
      ipsvarid = i
      psdimids = dimids(1:3)
    else if (varname == 'ppa' .or. varname == 'pp') then
      ppvarid = i
      cycle
    else if (varname == 'pai') then
      paivarid = i
      cycle
    else if (varname == 'p0') then
      ip0varid = i
    else if (varname == 'topo') then
      ishvarid = i
    else if (varname == 'rh') then
      has_rh = .true.
    else if (varname == 'ta' .or. varname == 't') then
      has_t = .true.
      intscheme(i) = 2
      tvarid = i
      tdimids = dimids(1:4)
      if ( varname == 't' ) is_icbc = .true.
    else if (varname == 'qas' .or. &
             varname == 'hus' .or. &
             varname == 'qv') then
      if ( varname == 'qas' .or. varname == 'hus' ) has_sph = .true.
      has_q = .true.
      qvarid = i
    else if (varname == 'chtrname') then
      lchnameflag(i) = .true.
    end if
    iv = 1
    do j = 1, nvdims(i)
      if (dimids(j) == kzdimid) then
        lkvarflag(i) = .true.
      else if (dimids(j) == itdimid) then
        ltvarflag(i) = .true.
        iv = iv + 1
        cycle
      end if
      dimsize(iv,i) = dimlen(dimids(j))
      iv = iv + 1
      varsize(i) = varsize(i) * dimlen(dimids(j))
    end do
    istatus = nf90_def_var(ncout, varname, xtype, dimids(1:nvdims(i)), ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Error define variable '//trim(varname))
    invarid(i) = i
    outvarid(i) = ivarid
#ifdef NETCDF4_HDF5
#ifdef NCFILTERS_AVAIL
    if (nvdims(i) > 2) then
      istatus = nf90_def_var_filter(ncout, ivarid, &
                  ncfilter,ncfilter_nparams,ncfilter_params)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error set filter for '//trim(varname))
    end if
#endif
#endif
    if (varname == 'zlev') then
      istatus = nf90_put_att(ncout, ivarid, 'standard_name', 'height')
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error adding hgt standard name')
      istatus = nf90_put_att(ncout, ivarid, 'long_name', &
                             'elevation')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt long name')
      istatus = nf90_put_att(ncout, ivarid, 'units', 'm')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt units')
      istatus = nf90_put_att(ncout, ivarid, 'axis', 'Z')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt axis')
      istatus = nf90_put_att(ncout, ivarid, 'positive', 'up')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding hgt positive')
      cycle
    end if
    do j = 1, nvatts
      istatus = nf90_inq_attname(ncid, i, j, attname)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error read att for '//trim(varname))
      istatus = nf90_copy_att(ncid, i, attname, ncout, ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error copy att '//trim(attname)//' for '//trim(varname))
    end do
  end do

  if ( has_t ) then
    make_mslp = .true.
    istart(1) = 1
    istart(2) = 1
    icount(1) = jx
    icount(2) = iy
    allocate(mslpr(jx,iy), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'mslpr')
    istatus = nf90_def_var(ncout, 'mslp', nf90_float, psdimids, imslzvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable mslp')
#ifdef NETCDF4_HDF5
#ifdef NCFILTERS_AVAIL
    istatus = nf90_def_var_filter(ncout, imslzvar, &
                  ncfilter,ncfilter_nparams,ncfilter_params)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set filter for mslp')
#endif
#endif
    istatus = nf90_put_att(ncout, imslzvar, 'standard_name', &
                     'air_pressure_at_sea_level')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, imslzvar, 'long_name', &
                     'Sea Level pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, imslzvar, 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
    istatus = nf90_put_att(ncout, imslzvar, 'coordinates', 'xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding coordinates')
    istatus = nf90_put_att(ncout, imslzvar, 'grid_mapping', 'crs')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding grid_mapping')
  end if
  if ( has_t .and. has_q .and. .not. has_rh ) then
    make_rh = .true.
    allocate(qazvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'qazvar')
    istatus = nf90_def_var(ncout, 'rh', nf90_float, tdimids, irhvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable rh')
#ifdef NETCDF4_HDF5
#ifdef NCFILTERS_AVAIL
    istatus = nf90_def_var_filter(ncout, irhvar, &
                  ncfilter,ncfilter_nparams,ncfilter_params)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set filter for rh')
#endif
#endif
    istatus = nf90_put_att(ncout, irhvar, 'standard_name', 'relative_humidity')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, irhvar, 'long_name', 'Relative Humidity')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, irhvar, 'units', '%')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
    istatus = nf90_put_att(ncout, irhvar, 'coordinates', 'xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding coordinates')
    istatus = nf90_put_att(ncout, irhvar, 'grid_mapping', 'crs')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding grid_mapping')
  end if

  istatus = nf90_inq_varid(ncid, "kz", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "sigma", ivarid)
  end if
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lev", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error searching variable kz, sigma or lev.')
    call memory_init
    call init_sigma(kz,0.05_rkx,0.01_rkx)
    do k = 1, kz
      sigma(k) = real(half_sigma_coordinate(k),rk4)
    end do
    call memory_destroy
  else
    istatus = nf90_get_var(ncid, ivarid, sigma)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error reading variable kz, sigma or lev.')
    if ( sigma(1) < dlowval ) then
      ! Fix for a buggy RegCM 4.3.x revision
      allocate(sigfix(kz+1))
      sigfix(1:kz) = sigma
      sigfix(kz+1) = d_one
      do k = 1, kz
        sigma(k) = 0.5*real(sigfix(k)+sigfix(k+1))
      end do
      deallocate(sigfix)
    end if
  end if

  istart(1) = 1
  istart(2) = 1
  icount(1) = jx
  icount(2) = iy
  istatus = nf90_get_var(ncid, ishvarid, topo, istart(1:2), icount(1:2))
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading topo.')

  if ( iodyn == 3 ) then
    istatus = nf90_inq_varid(ncid, "a", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable a.')
    istatus = nf90_get_var(ncid, ivarid, ak)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable a.')
    istatus = nf90_inq_varid(ncid, "b", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable b.')
    istatus = nf90_get_var(ncid, ivarid, bk)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable b.')
    do k = 1, kz
      hzvar(:,:,k) = ak(k) + bk(k) * topo(:,:)
    end do
  end if

  if ( iodyn == 2 ) then
    istatus = nf90_get_var(ncid, ip0varid, ps0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable p0.')
    ps0 = ps0 - real(ptop)
  end if

  istatus = nf90_inq_varid(ncid, "time", ivarid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable time.')
  istatus = nf90_get_var(ncid, ivarid, times)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable time.')

  istatus = nf90_enddef(ncout)
  call checkncerr(istatus,__FILE__,__LINE__, &
          'Error preparing for write on output')

  ! Write time independent variables

  do i = 1, nvars
    if ( ltvarflag(i) ) cycle
    if (i == avarid) cycle
    if (i == bvarid) cycle
    if (i == itvarid) cycle
    if (i == ircm_map) cycle
    if (i == paivarid) cycle
    if (i == ppvarid) cycle
    if (i == ikvarid) then
      istatus = nf90_put_var(ncout, outvarid(i), zlevs)
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable zlev')
      cycle
    end if
    if ( invarid(i) >= 0 ) then
      if ( .not. lchnameflag(i) ) then
        allocate(avar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'avar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, invarid(i), avar, &
                               istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, outvarid(i), avar, &
                               istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(avar)
      else
        allocate(txtvar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'txtvar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, invarid(i), txtvar, &
                               istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, outvarid(i), &
                               txtvar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(txtvar)
      end if
    end if
  end do

  ! Write time dependent variables

  do it = 1, nt
    istart(1) = it
    icount(1) = 1
    istatus = nf90_put_var(ncout, outvarid(itvarid), times(it:it), &
                           istart(1:1), icount(1:1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error writing time.')
    istart(1) = 1
    istart(2) = 1
    istart(3) = it
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncid, ipsvarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading ps.')
    istatus = nf90_put_var(ncout, outvarid(ipsvarid), ps, &
                           istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error writing ps.')
    if ( iodyn == 2 .and. ppvarid >= 0 ) then
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      icount(1) = jx
      icount(2) = iy
      icount(3) = kz
      icount(4) = 1
      istatus = nf90_get_var(ncid, ppvarid, pp, istart(1:4), icount(1:4))
      call checkncerr(istatus,__FILE__,__LINE__,'Error reading pp.')
      do k = 1, kz
        press(:,:,k) = ps0(:,:)*real(sigma(k)) + real(ptop) + pp(:,:,k)
      end do
    else if ( iodyn == 3 .and. paivarid >= 0) then
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      icount(1) = jx
      icount(2) = iy
      icount(3) = kz
      icount(4) = 1
      istatus = nf90_get_var(ncid, paivarid, pai, istart(1:4), icount(1:4))
      call checkncerr(istatus,__FILE__,__LINE__,'Error reading pai.')
      press = p00 * (pai**cpovr)
    end if
    if ( iodyn /= 3 ) then
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      icount(1) = jx
      icount(2) = iy
      icount(3) = kz
      icount(4) = 1
      istatus = nf90_get_var(ncid, tvarid, tazvar, istart(1:4), icount(1:4))
      call checkncerr(istatus,__FILE__,__LINE__,'Error reading temp.')
      if ( iodyn == 1 ) then
        call htsig_o(tazvar,hzvar,ps,topo,sigma,ptop,jx,iy,kz)
      else if ( iodyn == 2 ) then
        call nonhydrost(hzvar,tazvar,press,ps,topo,jx,iy,kz)
      end if
    end if
    do i = 1, nvars
      if (.not. ltvarflag(i)) cycle
      if (i == avarid) cycle
      if (i == bvarid) cycle
      if (i == itvarid) cycle
      if (i == ipsvarid) cycle
      if (i == ircm_map) cycle

      if (lkvarflag(i)) then

        ! Do interpolation

        iv = nvdims(i)
        if ( iv == 4 ) then
          istart(iv) = it
          icount(iv) = 1
          istart(1:iv-1) = 1
          icount(1:iv-1) = dimsize(1:iv-1,i)
          allocate(avar(varsize(i)),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'avar')
          istatus = nf90_get_var(ncid, invarid(i), avar, &
                                 istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error reading var to interpolate.')
          n3d = varsize(i) / i3d
          iz3d = p3d*n3d
          allocate(azvar(iz3d),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'azvar')
!$OMP PARALLEL DO
          do ii = 1, n3d
            xvar = reshape(avar((ii-1)*i3d+1:ii*i3d),[jx,iy,kz])
            call intlinz(zvar,xvar,hzvar,jx,iy,kz,zlevs,nz)
            azvar((ii-1)*iz3d+1:ii*iz3d) = reshape(zvar,[iz3d])
          end do
!$OMP END PARALLEL DO
          if ( i == qvarid .and. make_rh ) then
            qazvar = xvar
            if ( has_sph ) then
              call sph2mxr(qazvar,jx,iy,kz)
            end if
          end if
          icount(iv-1) = nz
          istatus = nf90_put_var(ncout, outvarid(i), azvar, &
                                 istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error writing interp variable.')
          deallocate(avar)
          deallocate(azvar)
        else if ( iv == 5 ) then
          istart(iv) = it
          icount(iv) = 1
          istart(1:iv-1) = 1
          icount(1:iv-1) = dimsize(1:iv-1,i)
          allocate(avar(varsize(i)),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'avar')
          istatus = nf90_get_var(ncid, invarid(i), avar, &
                                 istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error reading var to interpolate.')
          n3d = varsize(i) / dimsize(iv-1,i) / i3d
          iz3d = p3d*n3d
          do ich = 1, dimsize(iv-1,i)
            istart(iv) = it
            icount(iv) = 1
            istart(iv-1) = ich
            icount(iv-1) = 1
            istart(1:iv-2) = 1
            icount(1:iv-2) = dimsize(1:iv-2,i)
            icount(iv-2) = nz
            allocate(azvar(iz3d),stat=istatus)
            call checkalloc(istatus,__FILE__,__LINE__,'azvar')
!$OMP PARALLEL DO
            do ii = 1, n3d
              xvar = reshape(avar((ii-1)*i3d+(ich-1)*i3d+1:(ii+ich-1)*i3d), &
                             [jx,iy,kz])
              call intlinz(zvar,xvar,hzvar,jx,iy,kz,zlevs,nz)
              azvar((ii-1)*iz3d+1:ii*iz3d) = reshape(zvar,[iz3d])
            end do
!$OMP END PARALLEL DO
            istatus = nf90_put_var(ncout, outvarid(i), &
                                   azvar, istart(1:iv), icount(1:iv))
            call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error writing interp variable.')
            deallocate(azvar)
          end do
          deallocate(avar)
        end if

      else

        ! No interpolation

        iv = nvdims(i)
        istart(iv) = it
        icount(iv) = 1
        istart(1:iv-1) = 1
        icount(1:iv-1) = dimsize(1:iv-1,i)
        allocate(avar(varsize(i)),stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'avar')
        istatus = nf90_get_var(ncid, invarid(i), avar, &
                               istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error reading variable to pass')
        istatus = nf90_put_var(ncout, outvarid(i), avar, &
                               istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error writing variable to pass')
        deallocate(avar)

      end if

    end do

    if ( make_rh ) then
      if ( iodyn == 2 ) then
        call mxr2rh(tazvar,qazvar,press,jx,iy,kz)
      else
        call mxr2rh(tazvar,qazvar,ps,sigma,ptop,jx,iy,kz)
      end if
      call intlinz(zvar,qazvar,hzvar,jx,iy,kz,zlevs,nz)
      zvar = zvar * 100.0 ! Put in %
      iv = 4
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      icount(iv-1) = nz
      istatus = nf90_put_var(ncout, irhvar, zvar, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing rh variable.')
    end if
    if ( make_mslp ) then
      call mslp(tazvar,ps,topo,mslpr,jx,iy,kz)
      call gs_filter(mslpr,ps,jx,iy)
      iv = 3
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      istatus = nf90_put_var(ncout, imslzvar, mslpr, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing mslp variable.')
    end if
  end do

  deallocate(dimlen)
  deallocate(dimids)
  deallocate(lkvarflag)
  deallocate(ltvarflag)
  deallocate(varsize)
  deallocate(invarid)
  deallocate(outvarid)
  deallocate(nvdims)
  deallocate(dimsize)
  deallocate(sigma)
  deallocate(times)
  deallocate(ps)
  deallocate(xvar)
  deallocate(zvar)
  deallocate(tazvar)
  if ( make_rh ) then
    deallocate(qazvar)
  end if
  if ( make_mslp ) then
    deallocate(topo)
    deallocate(mslpr)
  end if

  if ( iodyn == 2 ) then
    deallocate(ps0)
    deallocate(pp)
    deallocate(press)
  else if ( iodyn == 3 ) then
    deallocate(press)
    deallocate(pai)
    deallocate(ak)
    deallocate(bk)
  end if
  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close input file '//trim(ncsfile))
  istatus = nf90_close(ncout)
  call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close output file '//trim(ncpfile))

  deallocate(zlevs)

  contains

  subroutine get_nz(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    integer(ik4) :: np ! Unused in sigma2z
    namelist /pp_param/ nz, np
    rewind(iu)
    read(iu, nml=pp_param, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading pp_param namelist in ',trim(ncsfile)
      write (stderr,*) 'Exiting...'
      stop
    end if
  end subroutine get_nz

  subroutine get_zlevs(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    namelist /height/ zlevs
    rewind(iu)
    read(iu, nml=height, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading height namelist in ',trim(ncsfile)
      write (stderr,*) 'Exiting...'
      stop
    end if
  end subroutine get_zlevs

#ifdef NETCDF4_HDF5
  subroutine get_ncfilter(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    namelist /ncfilters/ ncfilter, ncfilter_nparams, ncfilter_params
    rewind(iu)
    read(iu, nml=ncfilters, iostat=iresult)
  end subroutine get_ncfilter
#endif

end program sigma2z
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
