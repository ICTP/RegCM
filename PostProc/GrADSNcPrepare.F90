!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort

program ncprepare

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_projections
  use mod_stdio
  use mod_date
  use mod_message
  use mod_nchelper
  use mod_memutil
  use netcdf

  implicit none

  character(256) :: prgname , ncfile , tmpctl , tmpcoord , experiment , &
                    clmfile
  character(512) :: levels
  character(32) :: lvformat , varname
  character(64) :: vardesc , timeunit , timecal
  character(16) :: varunit
  character(16) :: dimdesc
  integer(ik4) :: numarg , istatus , ncid , ncid_clm
  type(rcm_time_and_date) :: idate1 , idate2
  type(rcm_time_interval) :: tdif
  integer(ik4) :: delta
  character(3) , dimension(12) :: cmon

  character(256) :: charatt
  character(6) :: iproj
  real(rkx) :: clat , clon , plat , plon , ds , centeri , centerj
  real(rkx) :: minlat , minlon , maxlat , maxlon , rlatinc , rloninc
  real(rkx) , dimension(2) :: icntr
  real(rkx) , dimension(2) :: trlat
  real(rkx) , allocatable , dimension(:,:) :: xlat , xlon
  real(rkx) , allocatable , dimension(:) :: level , tmplon
  real(rkx) , allocatable , dimension(:) :: times
  real(rkx) :: time1
  real(rk4) , allocatable , dimension(:,:) :: r4in , r4jn , r4uv
  real(rkx) , allocatable , dimension(:,:) :: rin , rjn , ruv
  logical , allocatable , dimension(:) :: lvarflag
  integer(ik4) , allocatable , dimension(:) :: dimids
  integer(ik4) :: ndims , nvars , natts , udimid , totvars
  integer(ik4) :: ivarid , idimid , xtype , ip1 , ip2
  integer(ik4) :: jxdimid , iydimid , kzdimid , itdimid , dptdimid
  integer(ik4) :: jx , iy , kz , nd , nt , nlat , nlon , ilat , ilon , isplit
  real(rk8) , dimension(:) , allocatable :: alon , alat
  real(rkx) :: flat , flon
  integer(ik4) :: i , j , iid
  integer(ik4) :: year , month , day , hour
  logical :: lvarsplit , existing , lsigma , ldepth , lu , lua , luas , lclm
  logical :: is_model_output = .false.
  logical :: uvrotate = .false.

  type(anyprojparams) :: pjpara
  type(regcm_projection) :: pj

  data cmon /'jan','feb','mar','apr','may','jun', &
             'jul','aug','sep','oct','nov','dec'/
  data lsigma /.true./
  data ldepth /.false./
  data lu /.false./
  data lua /.false./
  data luas /.false./
  data lclm /.false./
  data kzdimid  /-1/
  data itdimid  /-1/
  data dptdimid /-1/

  call get_command_argument(0,value=prgname)
  numarg = command_argument_count()
  if (numarg < 1) then
    write (stderr,*) 'Not enough arguments.'
    write (stderr,*) ' '
    write (stderr,*) 'Usage : ', trim(prgname), ' Rcmfile.nc'
    write (stderr,*) '        ', trim(prgname), ' clmoutput.nc DOMAIN.nc'
    write (stderr,*) '        ', trim(prgname), ' CHEMISS.nc DOMAIN.nc'
    write (stderr,*) ' '
    stop
  end if

  call get_command_argument(1,value=ncfile)

  call memory_init( )

  iid = scan(ncfile, '/', .true.)
  tmpctl = trim(ncfile)//'.ctl'

  if ( numarg == 2 ) then
    lclm = .true.
    clmfile = ncfile
    call get_command_argument(2,value=ncfile)
    iid = scan(clmfile, '/', .true.)
  end if

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error Open file '//trim(ncfile))

  if ( lclm ) then
    istatus = nf90_open(clmfile, nf90_nowrite, ncid_clm)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error Open file '//trim(clmfile))
    istatus = nf90_inquire(ncid_clm,ndims,nvars,natts,udimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire file '//trim(ncfile))
  else
    istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire file '//trim(ncfile))
  end if

  allocate(lvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__, &
                  'lvarflag')
  allocate(dimids(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__, &
                  'dimids')

  open(newunit=ip1, file=tmpctl, form='formatted', status='replace')
  if ( lclm ) then
    write(ip1, '(a)') 'dset ^'//trim(clmfile(iid+1:))
    write(ip1, '(a)') 'dtype netcdf'
    write(ip1, '(a)') 'undef 1e+20_FillValue'
  else
    write(ip1, '(a)') 'dset ^'//trim(ncfile(iid+1:))
    write(ip1, '(a)') 'dtype netcdf'
    write(ip1, '(a)') 'undef 1e+20_FillValue'
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'title', charatt)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error reading title attribute')
  write(ip1, '(a)') 'title '//trim(charatt)

  istatus = nf90_get_att(ncid, nf90_global, 'ipcc_scenario_code', charatt)
  if ( istatus == nf90_noerr ) then
    is_model_output = .true.
  end if

  istatus = nf90_get_att(ncid, nf90_global, &
                    'wind_rotated_eastward_northward', charatt)
  if ( istatus == nf90_noerr ) then
    uvrotate = .true.
  end if

  if ( lclm ) then
    istatus = nf90_inq_dimid(ncid_clm, "lon", jxdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid_clm, "x", jxdimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension x missing')
    istatus = nf90_inquire_dimension(ncid_clm, jxdimid, len=jx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dimension x')
    istatus = nf90_inq_dimid(ncid_clm, "lat", iydimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid_clm, "y", iydimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension y missing')
    istatus = nf90_inquire_dimension(ncid_clm, iydimid, len=iy)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dimension y')
  else
    istatus = nf90_inq_dimid(ncid, "jx", jxdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "x", jxdimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension x missing')
    istatus = nf90_inquire_dimension(ncid, jxdimid, len=jx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dimension x')
    istatus = nf90_inq_dimid(ncid, "iy", iydimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "y", iydimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension y missing')
    istatus = nf90_inquire_dimension(ncid, iydimid, len=iy)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dimension y')
  end if
  istatus = nf90_inq_dimid(ncid, "kz", kzdimid)
  if (istatus /= nf90_noerr) then
    lsigma = .false.
    istatus = nf90_inq_dimid(ncid, "lev", kzdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "plev", kzdimid)
      if (istatus /= nf90_noerr) then
        istatus = nf90_inq_dimid(ncid, "zlev", kzdimid)
      end if
    end if
  end if
  if (istatus == nf90_noerr) then
    istatus = nf90_inquire_dimension(ncid, kzdimid, len=kz)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dimension z')
  else
    kz = 0
  end if
  if ( lclm ) then
    istatus = nf90_inq_dimid(ncid_clm, "time", itdimid)
    if (istatus == nf90_noerr) then
      istatus = nf90_inquire_dimension(ncid_clm, itdimid, len=nt)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dimension time')
    else
      nt = 0
    end if
    istatus = nf90_inq_dimid(ncid_clm, "levsoi", dptdimid)
    if (istatus == nf90_noerr) then
      ldepth = .true.
      istatus = nf90_inquire_dimension(ncid_clm, dptdimid, len=nd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dimension depth')
    end if
  else
    istatus = nf90_inq_dimid(ncid, "time", itdimid)
    if (istatus == nf90_noerr) then
      istatus = nf90_inquire_dimension(ncid, itdimid, len=nt)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dimension time')
    else
      nt = 0
    end if
    istatus = nf90_inq_dimid(ncid, "depth", dptdimid)
    if (istatus == nf90_noerr) then
      ldepth = .true.
      istatus = nf90_inquire_dimension(ncid, dptdimid, len=nd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dimension depth')
    end if
  end if

#ifdef NETCDF4_HDF5
  if (ldepth) then
    if (iy*jx*64*4 > 524288) then
      write(ip1, '(a)') 'cachesize ', 524288 ! 1MB
    else
      write(ip1, '(a,i10)') 'cachesize ', iy*jx*nd*4
    end if
  else
    write(ip1, '(a,i10)') 'cachesize ', iy*jx*kz*4
  end if
#endif

  allocate(xlat(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__, &
                  'xlat')
  allocate(xlon(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__, &
                  'xlon')
  allocate(tmplon(iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__, &
                  'tmplon')

  istatus = nf90_inq_varid(ncid, "xlat", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lat", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find variable xlat')
  end if
  istatus = nf90_get_var(ncid, ivarid, xlat)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read variable xlat')
  istatus = nf90_inq_varid(ncid, "xlon", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lon", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find variable xlon')
  end if
  istatus = nf90_get_var(ncid, ivarid, xlon)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read variable xlon')

  istatus = nf90_get_att(ncid, nf90_global, 'projection', iproj)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute projection')

  istatus = nf90_get_att(ncid, nf90_global, 'experiment', experiment)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute experiment')
  istatus = nf90_get_att(ncid, nf90_global, 'index_of_projection_origin', &
                         icntr)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute index_of_projection_origin')
  istatus = nf90_get_att(ncid, nf90_global, 'latitude_of_projection_origin', &
                         clat)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute latitude_of_projection_origin')
  istatus = nf90_get_att(ncid, nf90_global, 'longitude_of_projection_origin', &
                         clon)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute longitude_of_projection_origin')
  istatus = nf90_get_att(ncid, nf90_global, 'grid_size_in_meters', ds)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read attribute grid_size_in_meters')

  minlat = rounder(minval(xlat),.false.)
  maxlat = rounder(maxval(xlat),.true.)
  if (abs(minlat+90.0_rkx)<0.001_rkx .or. abs(maxlat-90.0_rkx)<0.001_rkx) then
    minlon = -180.0_rkx
    maxlon = 180.0_rkx
  else if ( (xlon(jx,iy) - xlon(1,1)) < dlowval ) then
    minlon = -180.0_rkx
    maxlon = 180.0_rkx
  else
    tmplon(:) = xlon(jx,:)
    if ( (tmplon(1 ) > 0.0_rkx .and. tmplon(iy) < 0.0_rkx) .or. &
         (tmplon(iy) > 0.0_rkx .and. tmplon(1 ) < 0.0_rkx) ) then
      where ( tmplon < 0.0_rkx )
        tmplon = tmplon + 360.0_rkx
      endwhere
    end if
    maxlon = rounder(maxval(tmplon),.true.)
    tmplon(:) = xlon(1,:)
    if ( (tmplon(1 ) > 180.0_rkx .and. tmplon(iy) < 0.0_rkx) .or. &
         (tmplon(iy) > 180.0_rkx .and. tmplon(1 ) < 0.0_rkx) ) then
      where ( tmplon > 180.0_rkx )
        tmplon = tmplon - 360.0_rkx
      endwhere
    end if
    minlon = rounder(minval(tmplon),.false.)
  end if
  rlatinc = max(rounder(ds/111000.0_rkx,.false.),0.01_rkx)
  rloninc = max(rounder(ds/111000.0_rkx,.false.),0.01_rkx)
  nlat = nint(abs(maxlat-minlat)/rlatinc)
  if (minlon > 0.0_rkx .and. maxlon < 0.0_rkx) then
    nlon = nint(abs((maxlon+360.0_rkx)-minlon)/rloninc) + 1
  else if (minlon > 0.0_rkx .and. maxlon < 1e-30_rkx) then
    nlon = nint(360.0_rkx/rloninc) + 1
  else
    nlon = nint(abs(maxlon-minlon)/rloninc) + 1
  end if
  centeri = icntr(2)
  centerj = icntr(1)
  if ( is_model_output ) then
    centeri = centeri + 0.5_rkx
    centerj = centerj + 0.5_rkx
  end if
  deallocate(xlat)
  deallocate(xlon)

  if (iproj == 'LAMCON') then
    istatus = nf90_get_att(ncid, nf90_global, 'standard_parallel', trlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute standard_parallel')
  else if (iproj == 'ROTLLR' .or. iproj == 'ROTMER' ) then
    istatus = nf90_get_att(ncid, nf90_global, 'grid_north_pole_latitude', &
                           plat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_latitude')
    istatus = nf90_get_att(ncid, nf90_global, 'grid_north_pole_longitude', &
                           plon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_longitude')
  end if

  pjpara%pcode = iproj
  pjpara%ds = ds
  pjpara%clat = clat
  pjpara%clon = clon
  pjpara%plat = plat
  pjpara%plon = plon
  pjpara%trlat1 = trlat(1)
  pjpara%trlat2 = trlat(2)
  pjpara%nlon = jx
  pjpara%nlat = iy
  pjpara%staggerx = .false.
  pjpara%staggery = .false.
  pjpara%rotparam = .true.
  call pj%initialize(pjpara)

  allocate(alon(nlon),alat(nlat))
  if ( iproj == 'ROTLLR' ) then
    call pj%rl00(flat,flon)
    write(ip1, '(a,i8,i8,a,6f8.2)') 'pdef ', jx , iy ,   &
           ' rotll ',plon, plat, raddeg*ds/earthrad, &
           raddeg*ds/earthrad,flon,flat
  else
    if ( lclm ) then
      tmpcoord = trim(ncfile)//'.coord'
    else
      tmpcoord = trim(experiment)//'.coord'
    end if
    inquire (file=tmpcoord, exist=existing)
    if (.not. existing) then
      allocate(rin(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'rin')
      allocate(r4in(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'r4in')
      allocate(rjn(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'rjn')
      allocate(r4jn(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'r4jn')
      allocate(ruv(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'ruv')
      allocate(r4uv(nlon,nlat), stat=istatus)
      call checkalloc(istatus,__FILE__,__LINE__, &
                      'r4uv')
      do ilon = 1 , nlon
        alon(ilon) = minlon + (ilon-1) * rloninc
      end do
      do ilat = 1 , nlat
        alat(ilat) = minlat + (ilat-1) * rlatinc
      end do
      do ilon = 1 , nlon
        do ilat = 1 , nlat
          call pj%llij(alat(ilat),alon(ilon),rin(ilon,ilat),rjn(ilon,ilat))
        end do
      end do
      if ( iproj /= 'ROTLLR' ) then
        call pj%rotation_angle(alon,alat,ruv)
      end if
      r4in = real(rin)
      r4jn = real(rjn)
      r4uv = real(ruv)
      open(newunit=ip2, file=tmpcoord, form='unformatted', status='replace')
      write(ip2) r4in
      write(ip2) r4jn
      write(ip2) r4uv
      close(ip2)
      deallocate(rin,rjn,ruv)
      deallocate(r4in,r4jn,r4uv)
      write(ip1, '(a,i8,i8,a,a)') 'pdef ', jx , iy ,       &
             ' bilin sequential binary-big ', trim(tmpcoord)
    else
      write(stdout,*) 'Coordinate file exist, not recreating it'
      write(ip1, '(a,i8,i8,a,a)') 'pdef ', jx , iy ,       &
             ' bilin sequential binary-big ', trim(tmpcoord)
    end if
  end if

  write(ip1, '(a,i8,a,f7.2,f8.3)') 'xdef ', nlon , ' linear ',           &
         minlon, rloninc
  write(ip1, '(a,i8,a,f7.2,f8.3)') 'ydef ', nlat , ' linear ',           &
         minlat, rlatinc

  if (.not. ldepth .and. kz /= 0) then
    allocate(level(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__, &
                    'level')
    if (lsigma) then
      istatus = nf90_inq_varid(ncid, "kz", ivarid)
      if ( istatus /= nf90_noerr) then
        istatus = nf90_inq_varid(ncid, "sigma", ivarid)
      end if
    else
      istatus = nf90_inq_varid(ncid, "lev", ivarid)
      if ( istatus /= nf90_noerr) then
        istatus = nf90_inq_varid(ncid, "plev", ivarid)
        if (istatus /= nf90_noerr) then
          istatus = nf90_inq_varid(ncid, "zlev", ivarid)
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Vertical var not present')
    istatus = nf90_get_var(ncid, ivarid, level)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Read Z var')

    if (lsigma) level = level * 1000.0
    write (lvformat, '(a,i4,a)') '(a,i4,a,',kz,'f9.2)'
    write (levels, lvformat) 'zdef ', kz , ' levels ', level
    write (ip1, '(a)') trim(levels)
    deallocate(level)
  else if (ldepth) then
    write (ip1, '(a,i5,a,i5,i5)') 'zdef ', nd, ' linear ', 1, nd
  else
    write (ip1, '(a)') 'zdef 1 levels 1000.0'
  end if

  if (nt > 1) then
    allocate(times(nt), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__, &
                    'times')
    if ( lclm ) then
      istatus = nf90_inq_varid(ncid_clm, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time units not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time calendar not present')
      istatus = nf90_get_var(ncid_clm, ivarid, times)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Read time variable')
    else
      istatus = nf90_inq_varid(ncid, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable not present')
      istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time units not present')
      istatus = nf90_get_att(ncid, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time calendar not present')
      istatus = nf90_get_var(ncid, ivarid, times)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Read time variable')
    end if
    idate1 = timeval2date(times(1), timeunit, timecal)
    idate2 = timeval2date(times(2), timeunit, timecal)
    tdif = idate2-idate1
    delta = nint(tohours(tdif))
    deallocate(times)
    call split_idate(idate1,year,month,day,hour)
    if (delta == 24) then
      write (ip1, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month), year , ' 1dy'
    else if (delta == 168) then
      write (ip1, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month), year , ' 7dy'
    else if (delta >= 672 .and. delta <= 744) then
      write (ip1, '(a,i8,a,a,a1,a,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', '12', 'Z', '15', &
              cmon(month), year , ' 1mo'
    else if (delta > 8640) then
      write (ip1, '(a,i8,a,a,a1,a,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', '12' , 'Z', '15' , &
             '06' , year , ' 1yr'
    else
      write (ip1, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day,  &
             cmon(month), year , delta, 'hr'
    end if
  else if (nt == 1) then
    if ( lclm ) then
      istatus = nf90_inq_varid(ncid_clm, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time units not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time calendar not present')
      istatus = nf90_get_var(ncid_clm, ivarid, time1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable read')
    else
      istatus = nf90_inq_varid(ncid, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable not present')
      istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time units not present')
      istatus = nf90_get_att(ncid, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time calendar not present')
      istatus = nf90_get_var(ncid, ivarid, time1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Time variable read')
    end if
    idate1 = timeval2date(time1, timeunit, timecal)
    delta = 6
    call split_idate(idate1,year,month,day,hour)
    write (ip1, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
           'tdef ', nt, ' linear ', hour, 'Z', day,  &
           cmon(month), year , delta, 'hr'
  else
    write (ip1, '(a)') 'tdef 1 linear 00Z31dec1999 1yr'
  end if

  totvars = 0
  do i = 1 , nvars
    lvarflag(i) = .false.
    if ( lclm ) then
      istatus = nf90_inquire_variable(ncid_clm,i,name=varname,xtype=xtype, &
                                      ndims=idimid,dimids=dimids)
    else
      istatus = nf90_inquire_variable(ncid,i,name=varname,xtype=xtype, &
                                      ndims=idimid,dimids=dimids)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Inquire variable error')

    if ( varname == 'u' ) lu = .true.
    if ( varname == 'ua' ) lua = .true.
    if ( varname == 'uas' ) luas = .true.
    if (idimid > 1) then
      lvarflag(i) = .true.
      if (idimid == 2) then
        if (dimids(2) == iydimid) then
          totvars = totvars + 1
        else
          lvarflag(i) = .false.
        end if
      else if (idimid == 3) then
        if ((dimids(3) == kzdimid .or. dimids(3) == dptdimid) .and. &
             dimids(2) == iydimid) then
          totvars = totvars + 1
        else if (dimids(3) == itdimid .and. dimids(2) == iydimid) then
          totvars = totvars + 1
        else if (dimids(2) == iydimid .and. dimids(1) == jxdimid) then
          if ( lclm ) then
            istatus = nf90_inquire_dimension(ncid_clm, dimids(3), len=isplit)
          else
            istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
          end if
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Inquire split dim error')
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else if (idimid == 4) then
        if (dimids(4) == itdimid .and. &
            (dimids(3) == kzdimid .or. dimids(3) == dptdimid)) then
          totvars = totvars + 1
        else if (dimids(4) == itdimid .and. dimids(2) == iydimid) then
          if ( lclm ) then
            istatus = nf90_inquire_dimension(ncid_clm, dimids(3), len=isplit)
          else
            istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
          end if
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Inquire split dim error')
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else if (idimid == 5) then
        if (dimids(5) == itdimid .and. dimids(3) == kzdimid .and. &
            dimids(2) == iydimid) then
          if ( lclm ) then
            istatus = nf90_inquire_dimension(ncid_clm, dimids(4), len=isplit)
          else
            istatus = nf90_inquire_dimension(ncid, dimids(4), len=isplit)
          end if
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Inquire split dim error')
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else
        lvarflag(i) = .false.
      end if
    end if
  end do

  if ( .not. uvrotate ) then
    if ( lu ) then
      write (ip1, '(a)') 'vectorpairs u,v'
    else if ( lua ) then
      write (ip1, '(a)') 'vectorpairs ua,va'
    else if ( luas ) then
      write (ip1, '(a)') 'vectorpairs s01uas,s01vas'
    end if
  end if
  write (ip1, '(a,i4)') 'vars ', totvars

  do i = 1 , nvars
    if (lvarflag(i) .eqv. .false.) cycle
    lvarsplit = .false.
    if ( lclm ) then
      istatus = nf90_inquire_variable(ncid_clm,i,name=varname,ndims=idimid, &
                                      dimids=dimids)
    else
      istatus = nf90_inquire_variable(ncid,i,name=varname,ndims=idimid, &
                                      dimids=dimids)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Inquire variable error')
    if (idimid == 2) then
      if (dimids(2) == iydimid .and. dimids(1) == jxdimid) then
        dimdesc = ' 0 y,x'
      else
        cycle
      end if
    else if (idimid == 3) then
      if ((dimids(3) == kzdimid .or. dimids(3) == dptdimid) .and. &
           dimids(2) == iydimid) then
        if (ldepth) then
          write (dimdesc, '(a,i4,a)') ' ', nd, ' z,y,x'
        else
          write (dimdesc, '(a,i4,a)') ' ', kz, ' z,y,x'
        end if
      else if (dimids(3) == itdimid .and. dimids(2) == iydimid) then
        dimdesc = ' 0 t,y,x'
      else if (dimids(2) == iydimid) then
        lvarsplit = .true.
      else
        cycle
      end if
    else if (idimid == 4) then
      if (dimids(4) == itdimid .and. &
          (dimids(3) == kzdimid .or. dimids(3) == dptdimid)) then
        if (ldepth) then
          write (dimdesc, '(a,i4,a)') ' ', nd, ' t,z,y,x'
        else
          write (dimdesc, '(a,i4,a)') ' ', kz, ' t,z,y,x'
        end if
      else if (dimids(4) == itdimid .and. dimids(2) == iydimid) then
        lvarsplit = .true.
      else
        cycle
      end if
    else if (idimid == 5) then
      lvarsplit = .true.
    else
      cycle
    end if
    if ( lclm ) then
      istatus = nf90_get_att(ncid_clm,i,'long_name', vardesc)
      if ( istatus /= nf90_noerr ) then
        vardesc = 'unknown'
      end if
      istatus = nf90_get_att(ncid_clm,i,'units', varunit)
      if ( istatus /= nf90_noerr ) then
        varunit = '1'
      end if
    else
      istatus = nf90_get_att(ncid,i,'long_name', vardesc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Inquire variable '//trim(varname)//' long_name')
      istatus = nf90_get_att(ncid,i,'units', varunit)
      if ( istatus /= nf90_noerr ) then
        varunit = '1'
      end if
    end if

    if (.not. lvarsplit) then
      write (ip1, '(a,a,a,a,a,a,a,a,a)') trim(varname),'=>',trim(varname), &
                              trim(dimdesc), ' ', trim(vardesc) , ' (',   &
                              trim(varunit), ')'
    else if (idimid <= 4 .and. lvarsplit) then
      if ( lclm ) then
        istatus = nf90_inquire_dimension(ncid_clm, dimids(3), len=isplit)
      else
        istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Inquire split dimension')
      if (idimid == 3) then
        do j = 1 , isplit
          write (ip1, '(a,a,i0.2,a,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                                j, trim(varname), ' 0 ', j-1, ',y,x ',   &
                                trim(vardesc) , ' (', trim(varunit), ')'
        end do
      else if (idimid == 4) then
        do j = 1 , isplit
          write (ip1, '(a,a,i0.2,a,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                                j, trim(varname), ' 0 t,', j-1, ',y,x ',   &
                                trim(vardesc) , ' (', trim(varunit), ')'
        end do
      end if
    else if (idimid == 5 .and. lvarsplit) then
      if ( lclm ) then
        istatus = nf90_inquire_dimension(ncid_clm, dimids(4), len=isplit)
      else
        istatus = nf90_inquire_dimension(ncid, dimids(4), len=isplit)
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Inquire split dimension')
      do j = 1 , isplit
        if (ldepth) then
          write (ip1, '(a,a,i0.2,a,i2,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                    j, trim(varname)//' ', nd , ' t,', j-1, ',z,y,x ',   &
                    trim(vardesc) , ' (', trim(varunit), ')'
        else
          write (ip1, '(a,a,i0.2,a,i2,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                    j, trim(varname)//' ', kz, ' t,', j-1, ',z,y,x ',   &
                    trim(vardesc) , ' (', trim(varunit), ')'
        end if
      end do
    else
      cycle
    end if
  end do

  deallocate(alon,alat)
  deallocate(lvarflag)
  deallocate(dimids)

  write (ip1, '(a)') 'endvars'
  close(ip1)

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Close file error')
  if ( lclm ) then
    istatus = nf90_close(ncid_clm)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Close file error')
  end if

  call memory_destroy( )

end program ncprepare
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
