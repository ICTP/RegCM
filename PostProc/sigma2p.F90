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

program sigma2p

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam, only : iomode, dsmax, dsmin
  use mod_stdio, only : stderr
#ifdef NETCDF4_HDF5
  use mod_dynparam, only : ncfilter, ncfilter_nparams, ncfilter_params
#endif
  use mod_message
  use mod_vertint
  use mod_memutil
  use mod_sigma, only : init_sigma, half_sigma_coordinate
  use mod_zita
  use mod_hgt
  use mod_humid
  use mod_nchelper
  use netcdf

  implicit none

  character(256) :: prgname, ncsfile, ncpfile
  character(128) :: attname, dimname, varname, psunit
  integer(ik4) :: numarg, istatus, ncid, ncout

  integer(ik4), allocatable, dimension(:) :: dimids, dimlen
  real(rk4), allocatable, dimension(:) :: sigma, ak, bk
  real(rk4), allocatable, dimension(:,:,:) :: tmpvar, qvar, hzvar
  real(rk4), allocatable, dimension(:,:,:) :: pp, press, zeta, pai
  real(rk4), allocatable, dimension(:,:,:) :: fm, z0, tvar, qvvar
  real(rk4), allocatable, save, dimension(:,:,:) :: pvar, xvar
  real(rk4), allocatable, dimension(:,:) :: ps, topo, mslpr, ps0
  real(rk4), allocatable, dimension(:) :: avar
  character, allocatable, dimension(:) :: txtvar
  real(rk4), allocatable, dimension(:) :: apvar
  real(rkx), allocatable, dimension(:) :: times
  real(rkx), allocatable, dimension(:) :: sigfix
  logical, allocatable, dimension(:) :: lkvarflag, ltvarflag, lchnameflag
  integer(ik4), allocatable, dimension(:) :: invarid
  integer(ik4), allocatable, dimension(:) :: outvarid
  integer(ik4), allocatable, dimension(:) :: varsize
  integer(ik4), allocatable, dimension(:) :: intscheme
  integer(ik4), allocatable, dimension(:) :: nvdims
  integer(ik4), allocatable, dimension(:,:) :: dimsize
  integer(ik4), allocatable, dimension(:) :: istart, icount
  integer(ik4) :: ndims, nvars, natts, udimid, nvatts
  integer(ik4) :: ivarid, idimid, xtype
  integer(ik4) :: jxdimid, iydimid, kzdimid, itdimid, itvarid, ikvarid
  integer(ik4) :: avarid, bvarid, ipsvarid, ishvarid, ppvarid, ip0varid
  integer(ik4) :: paivarid
  integer(ik4) :: jx, iy, kz, nt
  real(rkx) :: ptop, dzita, ztop, zh, a0, htg
  integer(ik4), dimension(4) :: tdimids
  integer(ik4), dimension(3) :: psdimids
  integer(ik4) :: i, j, k, it, iv, iid1, iid2, ii, i3d, p3d, ich
  integer(ik4) :: tvarid, qvarid, irhvar, ihgvar, imslpvar, ircm_map
  logical :: has_t, has_q, has_rh, is_icbc
  logical :: make_rh, make_hgt, has_sph
  integer(ik4) :: n3d, ip3d, iodyn, np, ipunit, iresult
  real(rk4), allocatable, dimension(:) :: plevs
  real(rkx), allocatable, dimension(:) :: zita, zitah

  data has_t /.false./
  data has_q /.false./
  data has_rh /.false./
  data make_rh /.false./
  data make_hgt /.false./

!$OMP THREADPRIVATE(xvar,pvar)

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
    call get_np(ipunit)
    allocate(plevs(np))
    call get_plevs(ipunit)
    call get_command_argument(2,value=ncsfile)
    istatus = nf90_open(ncsfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
           'Error Opening Input file '//trim(ncsfile))
  else
    call checkncerr(istatus,__FILE__,__LINE__, &
           'Error Opening Input file '//trim(ncsfile))
    np = 12
    allocate(plevs(np))
    plevs = [1000.,925.,850.,700.,600.,500.,400.,300.,250.,200.,150.,100.]
  end if

  iid1 = scan(ncsfile, '/', .true.)
  iid2 = scan(ncsfile, '.', .true.)
  ncpfile = trim(ncsfile(iid1+1:iid2-1))//'_pressure.nc'

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
      istatus = nf90_def_dim(ncout, 'plev', np, idimid)
    else
      istatus = nf90_def_dim(ncout, dimname, dimlen(i), idimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Error define dimension '//trim(dimname))
  end do

  i3d = jx*iy*kz
  p3d = jx*iy*np

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
  allocate(pvar(jx,iy,np), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'pvar')
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
    allocate(ak(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'ak')
    allocate(bk(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'bk')
    allocate(zeta(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'zeta')
    allocate(pai(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'pai')
    allocate(press(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'press')
    allocate(pai(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'pai')
  end if

  ppvarid = -1
  ip0varid = -1
  paivarid = -1

  has_sph = .false.
  is_icbc = .false.
  ircm_map = -1
  invarid = -1
  outvarid = -1

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
    else if (varname == 'kz' .or. varname == 'sigma' .or. varname == 'lev') then
      varname = 'plev'
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
    else if ( varname == 'qas' .or. &
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
    if (varname == 'plev') then
      istatus = nf90_put_att(ncout, ivarid, 'standard_name', 'pressure')
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error adding pressure standard name')
      istatus = nf90_put_att(ncout, ivarid, 'long_name', &
                             'Interpolated pressure')
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error adding pressure long name')
      istatus = nf90_put_att(ncout, ivarid, 'units', 'hPa')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding pressure units')
      istatus = nf90_put_att(ncout, ivarid, 'axis', 'Z')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding pressure axis')
      istatus = nf90_put_att(ncout, ivarid, 'positive', 'down')
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error adding pressure positive')
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

  if ( is_icbc .and. iodyn == 3 ) then
    istatus = nf90_get_att(ncid, nf90_global, 'zita_height_top', ztop)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error reading attribute zita_height_top')
    istatus = nf90_get_att(ncid, nf90_global, 'zita_atmosphere_h', zh)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error reading attribute zita_atmosphere_h')
    istatus = nf90_get_att(ncid, nf90_global, 'zita_factor_a0', a0)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error reading attribute zita_factor_a0')
    allocate(fm(jx,iy,kz+1), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'fm')
    allocate(z0(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'z0')
    allocate(tvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'tvar')
    allocate(qvvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'qvvar')
    dzita = ztop/real(kz,rkx)
    allocate(zita(kz+1), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'zita')
    allocate(zitah(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'zitah')
    call model_zitaf(zita,ztop)
    call model_zitah(zitah,ztop)
  end if

  if ( has_t ) then
    make_hgt = .true.
    istart(1) = 1
    istart(2) = 1
    icount(1) = jx
    icount(2) = iy
    allocate(mslpr(jx,iy), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'mslpr')
    allocate(tmpvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'tmpvar')
    allocate(hzvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'hzvar')
    istatus = nf90_def_var(ncout, 'hgt', nf90_float, tdimids, ihgvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable hgt')
#ifdef NETCDF4_HDF5
#ifdef NCFILTERS_AVAIL
    istatus = nf90_def_var_filter(ncout, ihgvar, &
                  ncfilter,ncfilter_nparams,ncfilter_params)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set filter for hgt')
#endif
#endif
    istatus = nf90_put_att(ncout, ihgvar, 'standard_name', 'height')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, ihgvar,  &
            'long_name', 'Vertical distance above surface')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, ihgvar, 'units', 'm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
    istatus = nf90_put_att(ncout, ihgvar, '_FillValue', smissval)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding missval')
    istatus = nf90_put_att(ncout, ihgvar, 'coordinates', 'xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding coordinates')
    istatus = nf90_put_att(ncout, ihgvar, 'grid_mapping', 'crs')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding grid_mapping')
    istatus = nf90_def_var(ncout, 'mslp', nf90_float, psdimids, imslpvar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error define variable mslp')
#ifdef NETCDF4_HDF5
#ifdef NCFILTERS_AVAIL
    istatus = nf90_def_var_filter(ncout, imslpvar, &
                  ncfilter,ncfilter_nparams,ncfilter_params)
    call checkncerr(istatus,__FILE__,__LINE__,'Error set filter for mslp')
#endif
#endif
    istatus = nf90_put_att(ncout, imslpvar, 'standard_name', &
                     'air_pressure_at_sea_level')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding standard name')
    istatus = nf90_put_att(ncout, imslpvar, 'long_name', &
                     'Sea Level pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding long name')
    istatus = nf90_put_att(ncout, imslpvar, 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding units')
    istatus = nf90_put_att(ncout, imslpvar, 'coordinates', 'xlat xlon')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding coordinates')
    istatus = nf90_put_att(ncout, imslpvar, 'grid_mapping', 'crs')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding grid_mapping')
  end if
  if ( has_t .and. has_q .and. .not. has_rh ) then
    make_rh = .true.
    allocate(qvar(jx,iy,kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'qvar')
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
      ! height above the surface here
      zeta(:,:,k) = ak(k) + (bk(k) - d_one) * topo(:,:)
    end do
    if ( is_icbc ) then
      do k = 1, kz+1
        do i = 1, iy
          do j = 1, jx
            htg = topo(j,i)
            fm(j,i,k) = md_fmz_h(zita(k),htg,ztop,zh,a0)
          end do
        end do
      end do
      do k = 1, kz
        do i = 1, iy
          do j = 1, jx
            htg = topo(j,i)
            z0(j,i,k) = md_zeta_h(zitah(k),htg,ztop,zh,a0)
          end do
        end do
      end do
    end if
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
    if (i == ikvarid) then
      istatus = nf90_put_var(ncout, outvarid(i), plevs)
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable plev')
      cycle
    end if
    if ( invarid(i) >= 0 ) then
      if ( .not. lchnameflag(i) ) then
        allocate(avar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'avar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, invarid(i), &
                               avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, outvarid(i), &
                               avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(avar)
      else
        allocate(txtvar(varsize(i)), stat=istatus)
        call checkalloc(istatus,__FILE__,__LINE__,'txtvar')
        iv = nvdims(i)
        istart(:) = 1
        icount(1:iv) = dimsize(1:iv,i)
        istatus = nf90_get_var(ncid, invarid(i), &
                               txtvar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading variable.')
        istatus = nf90_put_var(ncout, outvarid(i), &
                               txtvar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__,'Error writing variable.')
        deallocate(txtvar)
      end if
    end if
  end do

  plevs = plevs * 100.0

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
    istatus = nf90_get_att(ncid, ipsvarid, 'units', psunit)
    call checkncerr(istatus,__FILE__,__LINE__,'Error reading ps.')
    if ( psunit == 'hPa' .or. psunit == 'mb' ) then
      ps = ps * d_100
    end if
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
    else if ( iodyn == 3 ) then
      if ( paivarid >= 0 ) then
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = it
        icount(1) = jx
        icount(2) = iy
        icount(3) = kz
        icount(4) = 1
        istatus = nf90_get_var(ncid, paivarid, pai, istart(1:4), icount(1:4))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading ta.')
      else
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = it
        icount(1) = jx
        icount(2) = iy
        icount(3) = kz
        icount(4) = 1
        istatus = nf90_get_var(ncid, tvarid, tvar, istart(1:4), icount(1:4))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading ta.')
        istatus = nf90_get_var(ncid, qvarid, qvvar, istart(1:4), icount(1:4))
        call checkncerr(istatus,__FILE__,__LINE__,'Error reading ta.')
        call paicompute(ps,z0,tvar,qvvar,pai,fm,dzita,jx,iy,kz)
      end if
      press = p00 * (pai**cpovr)
    end if
    do i = 1, nvars
      if (.not. ltvarflag(i)) cycle
      if (i == avarid) cycle
      if (i == bvarid) cycle
      if (i == paivarid) cycle
      if (i == ppvarid) cycle
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
          istatus = nf90_get_var(ncid, invarid(i), &
                                 avar, istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error reading var to interpolate.')
          n3d = varsize(i) / i3d
          ip3d = p3d*n3d
          allocate(apvar(ip3d),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'apvar')
!$OMP PARALLEL DO
          do ii = 1, n3d
            xvar = reshape(avar((ii-1)*i3d+1:ii*i3d),[jx,iy,kz])
            if ( iodyn == 2 ) then
              if (intscheme(i) == 1) then
                call intlin(pvar,xvar,press,jx,iy,kz,plevs,np)
              else if (intscheme(i) == 2) then
                call intlog(pvar,xvar,press,jx,iy,kz,plevs,np)
              end if
            else if ( iodyn == 3 ) then
              if (intscheme(i) == 1) then
                call intlin(pvar,xvar,press,jx,iy,kz,plevs,np)
              else if (intscheme(i) == 2) then
                call intlog(pvar,xvar,press,jx,iy,kz,plevs,np)
              end if
            else
              if (intscheme(i) == 1) then
                call intlin(pvar,xvar,ps,sigma,ptop,jx,iy,kz,plevs,np)
              else if (intscheme(i) == 2) then
                call intlog(pvar,xvar,ps,sigma,ptop,jx,iy,kz,plevs,np)
              end if
            end if
            apvar((ii-1)*ip3d+1:ii*ip3d) = reshape(pvar,[ip3d])
          end do
!$OMP END PARALLEL DO
          if ( i == tvarid ) then
            tmpvar = xvar
          else if ( i == qvarid .and. ( make_rh .or. iodyn == 3 ) ) then
            qvar = xvar
            if ( has_sph ) then
              call sph2mxr(qvar,jx,iy,kz)
            end if
          end if
          icount(iv-1) = np
          istatus = nf90_put_var(ncout, outvarid(i), apvar, &
                                 istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error writing interp variable.')
          deallocate(avar)
          deallocate(apvar)
        else if ( iv == 5 ) then
          istart(iv) = it
          icount(iv) = 1
          istart(1:iv-1) = 1
          icount(1:iv-1) = dimsize(1:iv-1,i)
          allocate(avar(varsize(i)),stat=istatus)
          call checkalloc(istatus,__FILE__,__LINE__,'avar')
          istatus = nf90_get_var(ncid, invarid(i), &
                                 avar, istart(1:iv), icount(1:iv))
          call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error reading var to interpolate.')
          n3d = varsize(i) / dimsize(iv-1,i) / i3d
          ip3d = p3d*n3d
          do ich = 1, dimsize(iv-1,i)
            istart(iv) = it
            icount(iv) = 1
            istart(iv-1) = ich
            icount(iv-1) = 1
            istart(1:iv-2) = 1
            icount(1:iv-2) = dimsize(1:iv-2,i)
            icount(iv-2) = np
            allocate(apvar(ip3d),stat=istatus)
            call checkalloc(istatus,__FILE__,__LINE__,'apvar')
!$OMP PARALLEL DO
            do ii = 1, n3d
              xvar = reshape(avar((ii-1)*i3d+(ich-1)*i3d+1:(ii+ich-1)*i3d), &
                             [jx,iy,kz])
              if ( iodyn == 2 ) then
                if (intscheme(i) == 1) then
                  call intlin(pvar,xvar,press,jx,iy,kz,plevs,np)
                else if (intscheme(i) == 2) then
                  call intlog(pvar,xvar,press,jx,iy,kz,plevs,np)
                end if
              else if ( iodyn == 3 ) then
                if (intscheme(i) == 1) then
                  call intlin(pvar,xvar,press,jx,iy,kz,plevs,np)
                else if (intscheme(i) == 2) then
                  call intlog(pvar,xvar,press,jx,iy,kz,plevs,np)
                end if
              else
                if (intscheme(i) == 1) then
                  call intlin(pvar,xvar,ps,sigma,ptop,jx,iy,kz,plevs,np)
                else if (intscheme(i) == 2) then
                  call intlog(pvar,xvar,ps,sigma,ptop,jx,iy,kz,plevs,np)
                end if
              end if
              apvar((ii-1)*ip3d+1:ii*ip3d) = reshape(pvar,[ip3d])
            end do
!$OMP END PARALLEL DO
            istatus = nf90_put_var(ncout, outvarid(i), apvar, &
                                   istart(1:iv), icount(1:iv))
            call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error writing interp variable.')
            deallocate(apvar)
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
        istatus = nf90_get_var(ncid, invarid(i), &
                               avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error reading variable to pass')
        istatus = nf90_put_var(ncout, outvarid(i), &
                               avar, istart(1:iv), icount(1:iv))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error writing variable to pass')
        deallocate(avar)
      end if
    end do

    if ( make_rh ) then
      if ( iodyn == 2 ) then
        call mxr2rh(tmpvar,qvar,press,jx,iy,kz)
        call intlin(pvar,qvar,press,jx,iy,kz,plevs,np)
      else if ( iodyn == 3 ) then
        call mxr2rh(tmpvar,qvar,press,jx,iy,kz)
        call intlin(pvar,qvar,press,jx,iy,kz,plevs,np)
      else
        call mxr2rh(tmpvar,qvar,ps,sigma,ptop,jx,iy,kz)
        call intlin(pvar,qvar,ps,sigma,ptop,jx,iy,kz,plevs,np)
      end if
      pvar = pvar * 100.0 ! Put in %
      iv = 4
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      icount(iv-1) = np
      istatus = nf90_put_var(ncout, irhvar, pvar, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing rh variable.')
    end if
    if ( make_hgt ) then
      if ( iodyn == 2 ) then
        call nonhydrost(hzvar,tmpvar,press,ps,topo,jx,iy,kz)
        call height_o(pvar,hzvar,tmpvar,ps,topo,press,jx,iy,kz,plevs,np)
      else if ( iodyn == 3 ) then
        do k = 1, kz
          zeta(:,:,k) = zeta(:,:,k) + topo
        end do
        call height_o(pvar,zeta,tmpvar,ps,topo,press,jx,iy,kz,plevs,np)
        do k = 1, kz
          zeta(:,:,k) = zeta(:,:,k) - topo
        end do
      else
        call htsig_o(tmpvar,hzvar,ps,topo,sigma,ptop,jx,iy,kz)
        call height_o(pvar,hzvar,tmpvar,ps,topo,sigma,ptop,jx,iy,kz,plevs,np)
      end if
      iv = 4
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      icount(iv-1) = np
      istatus = nf90_put_var(ncout, ihgvar, pvar, istart(1:iv), icount(1:iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error writing hgt variable.')
      call mslp(tmpvar,ps,topo,mslpr,jx,iy,kz)
      call gs_filter(mslpr,ps,jx,iy)
      iv = 3
      istart(iv) = it
      icount(iv) = 1
      istart(1:iv-1) = 1
      icount(1:iv-1) = dimsize(1:iv-1,tvarid)
      istatus = nf90_put_var(ncout, imslpvar, mslpr, istart(1:iv), icount(1:iv))
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
  deallocate(topo)
  deallocate(ps)
  deallocate(xvar)
  deallocate(pvar)
  if ( allocated(qvar) ) then
    deallocate(qvar)
  end if
  if ( has_t ) then
    deallocate(hzvar)
    deallocate(tmpvar)
    deallocate(mslpr)
  end if

  if ( iodyn == 2 ) then
    deallocate(ps0)
    deallocate(pp)
    deallocate(press)
  else if ( iodyn == 3 ) then
    deallocate(press)
    deallocate(zeta)
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

  deallocate(plevs)

  contains

  subroutine top2btm(x,nlon1,nlat1,nlev1)
    implicit none
    integer(ik4), intent(in) :: nlat1, nlev1, nlon1
    real(rk4), intent(inout), dimension(nlon1,nlat1,nlev1) :: x

    integer(ik4) :: i, j, k, kr
    real(rk4), dimension(nlev1) :: work

    do j = 1, nlat1
      do i = 1, nlon1
        do k = 1, nlev1
          work(k) = x(i,j,k)
        end do
        do k = 1, nlev1
          kr = nlev1 - k + 1
          x(i,j,k) = work(kr)
        end do
      end do
    end do
  end subroutine top2btm

  subroutine get_np(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    integer(ik4) :: nz ! Unused in sigma2p
    namelist /pp_param/ nz, np
    rewind(iu)
    read(iu, nml=pp_param, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading pp_param namelist in ',trim(ncsfile)
      write (stderr,*) 'Exiting...'
      stop
    end if
  end subroutine get_np

  subroutine get_plevs(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    namelist /pressure/ plevs
    rewind(iu)
    read(iu, nml=pressure, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading pressure namelist in ',trim(ncsfile)
      write (stderr,*) 'Exiting...'
      stop
    end if
  end subroutine get_plevs

#ifdef NETCDF4_HDF5
  subroutine get_ncfilter(iu)
    implicit none
    integer(ik4), intent(in) :: iu
    namelist /ncfilters/ ncfilter, ncfilter_nparams, ncfilter_params
    rewind(iu)
    read(iu, nml=ncfilters, iostat=iresult)
  end subroutine get_ncfilter
#endif

  subroutine paicompute(ps,z,t,q,pai,fm,dzita,nx,ny,nz)
    implicit none
    integer(ik4), intent(in) :: nx, ny, nz
    real(rkx), intent(in) :: dzita
    real(rk4), dimension(nx,ny), intent(in) :: ps
    real(rk4), dimension(nx,ny,nz), intent(in) :: z, t, q
    real(rk4), dimension(nx,ny,nz+1), intent(in) :: fm
    real(rk4), dimension(nx,ny,nz), intent(inout) :: pai
    real(rk4) :: tv, tv1, tv2, p, zb, zdelta, zz, lrt
    integer(ik4) :: i, j, k
    ! Hydrostatic initialization of pai
    do i = 1, ny
      do j = 1, nx
        zdelta = z(j,i,nz)*egrav
        tv1 = t(j,i,nz) * (d_one + ep1*q(j,i,nz))
        tv2 = t(j,i,nz-1) * (d_one + ep1*q(j,i,nz-1))
        lrt = (tv1-tv2)/(z(j,i,nz-1)-z(j,i,nz))
        tv = tv1 + 0.5_rk4*z(j,i,nz)*lrt
        zz = d_one/(rgas*tv)
        p = ps(j,i) * exp(-zdelta*zz)
        pai(j,i,nz) = (p/p00)**rovcp
      end do
    end do
    do k = nz-1, 1, -1
      do i = 1, ny
        do j = 1, nx
          tv1 = t(j,i,k) * (d_one + ep1*q(j,i,k))
          tv2 = t(j,i,k+1) * (d_one + ep1*q(j,i,k+1))
          zb = d_two*egrav*dzita/(fm(j,i,k+1)*cpd) + tv1 - tv2
          zdelta = sqrt(zb**2 + d_four * tv2 * tv1)
          pai(j,i,k) = -pai(j,i,k+1) / (d_two * tv2) * (zb - zdelta)
        end do
      end do
    end do
  end subroutine paicompute

end program sigma2p
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
