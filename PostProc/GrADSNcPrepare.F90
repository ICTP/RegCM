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
  call abort
end subroutine myabort

program ncprepare

  use mod_intkinds
  use mod_realkinds
  use mod_projections
  use mod_stdio
  use mod_date
  use mod_message
  use mod_nchelper
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
  real(rk8) :: clat , clon , plat , plon , ds , centeri , centerj
  real(rk8) :: minlat , minlon , maxlat , maxlon , rlatinc , rloninc
  real(rk8) , dimension(2) :: trlat
  real(rk8) , allocatable , dimension(:,:) :: xlat , xlon
  real(rk8) , allocatable , dimension(:) :: level , tmplon
  real(rk8) , allocatable , dimension(:) :: times
  real(rk8) :: time1
  real(rk4) , allocatable , dimension(:,:) :: r4in , r4jn , r4uv
  real(rk8) , allocatable , dimension(:,:) :: rin , rjn , ruv
  logical , allocatable , dimension(:) :: lvarflag
  integer(ik4) , allocatable , dimension(:) :: dimids
  integer(ik4) :: ndims , nvars , natts , udimid , totvars
  integer(ik4) :: ivarid , idimid , xtype
  integer(ik4) :: jxdimid , iydimid , kzdimid , itdimid , dptdimid
  integer(ik4) :: jx , iy , kz , nd , nt , nlat , nlon , ilat , ilon , isplit
  real(rk8) :: alat , alon , angle
  integer(ik4) :: i , j , iid
  integer(ik4) :: year , month , day , hour
  logical :: lvarsplit , existing , lsigma , ldepth , lu , lua , luas , lclm
  logical :: is_model_output = .false.

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
  iid = scan(ncfile, '/', .true.)
  tmpctl = trim(ncfile)//'.ctl'

  if ( numarg == 2 ) then
    lclm = .true.
    clmfile = ncfile
    call get_command_argument(2,value=ncfile)
    iid = scan(clmfile, '/', .true.)
  end if

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error Open file '//trim(ncfile))

  if ( lclm ) then
    istatus = nf90_open(clmfile, nf90_nowrite, ncid_clm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error Open file '//trim(clmfile))
    istatus = nf90_inquire(ncid_clm,ndims,nvars,natts,udimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire file '//trim(ncfile))
  else
    istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire file '//trim(ncfile))
  end if

  allocate(lvarflag(nvars), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'lvarflag')
  allocate(dimids(ndims), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'dimids')

  open(11, file=tmpctl, form='formatted', status='replace')
  if ( lclm ) then
    write(11, '(a)') 'dset ^'//trim(clmfile(iid+1:))
    write(11, '(a)') 'dtype netcdf'
    write(11, '(a)') 'undef 1e+36_FillValue'
  else
    write(11, '(a)') 'dset ^'//trim(ncfile(iid+1:))
    write(11, '(a)') 'dtype netcdf'
    write(11, '(a)') 'undef 1e+20_FillValue'
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'title', charatt)
  call checkncerr(istatus,__FILE__,__LINE__,'Error reading title attribute')
  write(11, '(a)') 'title '//trim(charatt)

  istatus = nf90_get_att(ncid, nf90_global, 'ipcc_scenario_code', charatt)
  if ( istatus == nf90_noerr ) then
    is_model_output = .true.
  end if

  if ( lclm ) then
    istatus = nf90_inq_dimid(ncid_clm, "lon", jxdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid_clm, "x", jxdimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension x missing')
    istatus = nf90_inquire_dimension(ncid_clm, jxdimid, len=jx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension x')
    istatus = nf90_inq_dimid(ncid_clm, "lat", iydimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid_clm, "y", iydimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension y missing')
    istatus = nf90_inquire_dimension(ncid_clm, iydimid, len=iy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension y')
  else
    istatus = nf90_inq_dimid(ncid, "jx", jxdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "x", jxdimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension x missing')
    istatus = nf90_inquire_dimension(ncid, jxdimid, len=jx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension x')
    istatus = nf90_inq_dimid(ncid, "iy", iydimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "y", iydimid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension y missing')
    istatus = nf90_inquire_dimension(ncid, iydimid, len=iy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension y')
  end if
  istatus = nf90_inq_dimid(ncid, "kz", kzdimid)
  if (istatus /= nf90_noerr) then
    lsigma = .false.
    istatus = nf90_inq_dimid(ncid, "lev", kzdimid)
    if (istatus /= nf90_noerr) then
      istatus = nf90_inq_dimid(ncid, "plev", kzdimid)
    end if
  end if
  if (istatus == nf90_noerr) then
    istatus = nf90_inquire_dimension(ncid, kzdimid, len=kz)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension z')
  else
    kz = 0
  end if
  if ( lclm ) then
    istatus = nf90_inq_dimid(ncid_clm, "time", itdimid)
    if (istatus == nf90_noerr) then
      istatus = nf90_inquire_dimension(ncid_clm, itdimid, len=nt)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension time')
    else
      nt = 0
    end if
    istatus = nf90_inq_dimid(ncid_clm, "levsoi", dptdimid)
    if (istatus == nf90_noerr) then
      ldepth = .true.
      istatus = nf90_inquire_dimension(ncid_clm, dptdimid, len=nd)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension depth')
    end if
  else
    istatus = nf90_inq_dimid(ncid, "time", itdimid)
    if (istatus == nf90_noerr) then
      istatus = nf90_inquire_dimension(ncid, itdimid, len=nt)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension time')
    else
      nt = 0
    end if
    istatus = nf90_inq_dimid(ncid, "depth", dptdimid)
    if (istatus == nf90_noerr) then
      ldepth = .true.
      istatus = nf90_inquire_dimension(ncid, dptdimid, len=nd)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension depth')
    end if
  end if

#ifdef NETCDF4_HDF5
  if (ldepth) then
    if (iy*jx*64*4 > 524288) then
      write(11, '(a)') 'cachesize ', 524288 ! 1MB
    else
      write(11, '(a,i10)') 'cachesize ', iy*jx*nd*4
    end if
  else
    write(11, '(a,i10)') 'cachesize ', iy*jx*kz*4
  end if
#endif

  allocate(xlat(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'xlat')
  allocate(xlon(jx,iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'xlon')
  allocate(tmplon(iy), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'tmplon')

  istatus = nf90_inq_varid(ncid, "xlat", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lat", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find variable xlat')
  end if
  istatus = nf90_get_var(ncid, ivarid, xlat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read variable xlat')
  istatus = nf90_inq_varid(ncid, "xlon", ivarid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_varid(ncid, "lon", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find variable xlon')
  end if
  istatus = nf90_get_var(ncid, ivarid, xlon)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read variable xlon')

  istatus = nf90_get_att(ncid, nf90_global, 'projection', iproj)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read attribute projection')

  istatus = nf90_get_att(ncid, nf90_global, 'experiment', experiment)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read attribute experiment')

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
  if (abs(minlat+90.0)<0.001 .or. abs(maxlat-90.0)<0.001) then
    minlon = -180.0
    maxlon = 180.0
  else if ( (xlon(jx,iy) - xlon(1,1)) < dlowval ) then
    minlon = -180.0
    maxlon = 180.0
  else
    tmplon(:) = xlon(jx,:)
    if ( (tmplon(1 ) > 0.0 .and. tmplon(iy) < 0.0) .or. &
         (tmplon(iy) > 0.0 .and. tmplon(1 ) < 0.0) ) then
      where ( tmplon < 0.0 )
        tmplon = tmplon + 360.0
      endwhere
    end if
    maxlon = rounder(maxval(tmplon),.true.)
    tmplon(:) = xlon(1,:)
    if ( (tmplon(1 ) > 180.0 .and. tmplon(iy) < 0.0) .or. &
         (tmplon(iy) > 180.0 .and. tmplon(1 ) < 0.0) ) then
      where ( tmplon > 180.0 )
        tmplon = tmplon - 360.0
      endwhere
    end if
    minlon = rounder(minval(tmplon),.false.)
  end if
  rlatinc = rounder(ds/111000.0/2.0,.false.)
  rloninc = rounder(ds/111000.0/2.0,.false.)
  nlat = nint(abs(maxlat-minlat)/rlatinc)
  if (minlon > 0.0 .and. maxlon < 0.0) then
    nlon = nint(abs((maxlon+360.0)-minlon)/rloninc) + 1
  else if (minlon > 0.0 .and. maxlon < 1e-30) then
    nlon = nint(360.0/rloninc) + 1
  else
    nlon = nint(abs(maxlon-minlon)/rloninc) + 1
  end if
  if ( is_model_output ) then
    centeri = dble(iy)/2.0D0+0.5
    centerj = dble(jx)/2.0D0+0.5
  else
    centeri = dble(iy)/2.0D0
    centerj = dble(jx)/2.0D0
  end if
  deallocate(xlat)
  deallocate(xlon)

  allocate(rin(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'rin')
  allocate(r4in(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'r4in')
  allocate(rjn(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'rjn')
  allocate(r4jn(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'r4jn')
  allocate(ruv(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'ruv')
  allocate(r4uv(nlon,nlat), stat=istatus)
  call checkalloc(istatus,__FILE__,__LINE__,'r4uv')

  if (iproj == 'LAMCON') then
    istatus = nf90_get_att(ncid, nf90_global, 'standard_parallel', trlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute standard_parallel')
    call setup_lcc(clat,clon,centerj,centeri,ds,clon,trlat(1),trlat(2))
    do ilon = 1 , nlon
      alon = minlon + (ilon-1) * rloninc
      call uvrot_lc(alon,angle)
      do ilat = 1 , nlat
        alat = minlat + (ilat-1) * rlatinc
        call llij_lc(alat,alon,rin(ilon,ilat),rjn(ilon,ilat))
        ruv(ilon,ilat) = angle
      end do
    end do
  else if (iproj == 'POLSTR') then
    call setup_plr(clat,clon,centerj,centeri,ds,clon)
    do ilon = 1 , nlon
      alon = minlon + (ilon-1) * rloninc
      call uvrot_ps(alon,angle)
      do ilat = 1 , nlat
        alat = minlat + (ilat-1) * rlatinc
        call llij_ps(alat,alon,rin(ilon,ilat),rjn(ilon,ilat))
        ruv(ilon,ilat) = angle
      end do
    end do
  else if (iproj == 'NORMER') then
    call setup_mrc(clat,clon,centerj,centeri,ds)
    do ilon = 1 , nlon
      alon = minlon + (ilon-1) * rloninc
      do ilat = 1 , nlat
        alat = minlat + (ilat-1) * rlatinc
        call llij_mc(alat,alon,rin(ilon,ilat),rjn(ilon,ilat))
        ruv(ilon,ilat) = 1.0
      end do
    end do
  else if (iproj == 'ROTMER') then
    istatus = nf90_get_att(ncid, nf90_global, 'grid_north_pole_latitude', &
                           plat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_latitude')
    istatus = nf90_get_att(ncid, nf90_global, 'grid_north_pole_longitude', &
                           plon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_longitude')
    call setup_rmc(clat,clon,centerj,centeri,ds,plon,plat)
    do ilon = 1 , nlon
      alon = minlon + (ilon-1) * rloninc
      do ilat = 1 , nlat
        alat = minlat + (ilat-1) * rlatinc
        call llij_rc(alat,alon,rin(ilon,ilat),rjn(ilon,ilat))
        call uvrot_rc(alat,alon,angle)
        ruv(ilon,ilat) = angle
      end do
    end do
  else
    write (stderr,*) 'Unknown Projection : ', iproj
    stop
  end if

  tmpcoord = ncfile(1:iid)//trim(experiment)//'.coord'
  inquire (file=tmpcoord, exist=existing)
  if (.not. existing) then
    r4in = real(rin)
    r4jn = real(rjn)
    r4uv = real(ruv)
    open(12, file=tmpcoord, form='unformatted', status='replace')
    write(12) r4in
    write(12) r4jn
    write(12) r4uv
    close(12)
  else
    write(stdout,*) 'Coordinate file exist, not recreating it'
  end if

  deallocate(rin,rjn,ruv)
  deallocate(r4in,r4jn,r4uv)

  if ( lclm ) then
    write(11, '(a,i8,i8,a,a)') 'pdef ', jx , iy ,                         &
           ' bilin sequential binary-big ', trim(tmpcoord)
  else
    write(11, '(a,i8,i8,a,a)') 'pdef ', jx , iy ,                         &
           ' bilin sequential binary-big ^', trim(experiment)//'.coord'
  end if
  write(11, '(a,i8,a,f7.2,f7.2)') 'xdef ', nlon , ' linear ',           &
         minlon, rloninc
  write(11, '(a,i8,a,f7.2,f7.2)') 'ydef ', nlat , ' linear ',           &
         minlat, rlatinc

  if (.not. ldepth .and. kz /= 0) then
    allocate(level(kz), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'level')
    if (lsigma) then
      istatus = nf90_inq_varid(ncid, "sigma", ivarid)
    else
      istatus = nf90_inq_varid(ncid, "lev", ivarid)
      if ( istatus /= nf90_noerr) then
        istatus = nf90_inq_varid(ncid, "plev", ivarid)
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, 'Vertical var not present')
    istatus = nf90_get_var(ncid, ivarid, level)
    call checkncerr(istatus,__FILE__,__LINE__, 'Read Z var')

    if (lsigma) level = level * 1000.0
    write (lvformat, '(a,i4,a)') '(a,i4,a,',kz,'f7.1)'
    write (levels, lvformat) 'zdef ', kz , ' levels ', level
    write (11, '(a)') trim(levels)
    deallocate(level)
  else if (ldepth) then
    write (11, '(a,i5,a,i5,i5)') 'zdef ', nd, ' linear ', 1, nd
  else
    write (11, '(a)') 'zdef 1 levels 1000.0'
  end if

  if (nt > 1) then
    allocate(times(nt), stat=istatus)
    call checkalloc(istatus,__FILE__,__LINE__,'times')
    if ( lclm ) then
      istatus = nf90_inq_varid(ncid_clm, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time units not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time calendar not present')
      istatus = nf90_get_var(ncid_clm, ivarid, times)
      call checkncerr(istatus,__FILE__,__LINE__, 'Read time variable')
    else
      istatus = nf90_inq_varid(ncid, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable not present')
      istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time units not present')
      istatus = nf90_get_att(ncid, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time calendar not present')
      istatus = nf90_get_var(ncid, ivarid, times)
      call checkncerr(istatus,__FILE__,__LINE__, 'Read time variable')
    end if
    idate1 = timeval2date(times(1), timeunit, timecal)
    idate2 = timeval2date(times(2), timeunit, timecal)
    tdif = idate2-idate1
    delta = idnint(tohours(tdif))
    deallocate(times)
    call split_idate(idate1,year,month,day,hour)
    if (delta == 24) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month), year , ' 1dy'
    else if (delta == 168) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month), year , ' 7dy'
    else if (delta >= 672 .and. delta <= 744) then
      write (11, '(a,i8,a,a,a1,a,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', '12', 'Z', '15', &
              cmon(month), year , ' 1mo'
    else if (delta > 8640) then
      write (11, '(a,i8,a,a,a1,a,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', '12' , 'Z', '15' , &
             '06' , year , ' 1yr'
    else
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day,  &
             cmon(month), year , delta, 'hr'
    end if
  else if (nt == 1) then
    if ( lclm ) then
      istatus = nf90_inq_varid(ncid_clm, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time units not present')
      istatus = nf90_get_att(ncid_clm, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time calendar not present')
      istatus = nf90_get_var(ncid_clm, ivarid, time1)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable read')
    else
      istatus = nf90_inq_varid(ncid, "time", ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable not present')
      istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time units not present')
      istatus = nf90_get_att(ncid, ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time calendar not present')
      istatus = nf90_get_var(ncid, ivarid, time1)
      call checkncerr(istatus,__FILE__,__LINE__, 'Time variable read')
    end if
    idate1 = timeval2date(time1, timeunit, timecal)
    delta = 6
    call split_idate(idate1,year,month,day,hour)
    write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
           'tdef ', nt, ' linear ', hour, 'Z', day,  &
           cmon(month), year , delta, 'hr'
  else
    write (11, '(a)') 'tdef 1 linear 00Z31dec1999 1yr'
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
    call checkncerr(istatus,__FILE__,__LINE__, 'Inquire variable error')

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
          call checkncerr(istatus,__FILE__,__LINE__, 'Inquire split dim error')
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
          call checkncerr(istatus,__FILE__,__LINE__, 'Inquire split dim error')
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
          call checkncerr(istatus,__FILE__,__LINE__, 'Inquire split dim error')
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else
        lvarflag(i) = .false.
      end if
    end if
  end do

  if ( lu ) then
    write (11, '(a)') 'vectorpairs u,v'
  else if ( lua ) then
    write (11, '(a)') 'vectorpairs ua,va'
  else if ( luas ) then
    write (11, '(a)') 'vectorpairs s01uas,s01vas'
  end if
  write (11, '(a,i4)') 'vars ', totvars

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
    call checkncerr(istatus,__FILE__,__LINE__, 'Inquire variable error')
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
      call checkncerr(istatus,__FILE__,__LINE__, 'Inquire variable long_name')
      istatus = nf90_get_att(ncid,i,'units', varunit)
      call checkncerr(istatus,__FILE__,__LINE__, 'Inquire variable units')
    end if

    if (.not. lvarsplit) then
      write (11, '(a,a,a,a,a,a,a,a,a)') trim(varname),'=>',trim(varname), &
                              trim(dimdesc), ' ', trim(vardesc) , ' (',   &
                              trim(varunit), ')'
    else if (idimid <= 4 .and. lvarsplit) then
      if ( lclm ) then
        istatus = nf90_inquire_dimension(ncid_clm, dimids(3), len=isplit)
      else
        istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
      end if
      call checkncerr(istatus,__FILE__,__LINE__, 'Inquire split dimension')
      if (idimid == 3) then
        do j = 1 , isplit
          write (11, '(a,a,i0.2,a,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                                j, trim(varname), ' 0 ', j-1, ',y,x ',   &
                                trim(vardesc) , ' (', trim(varunit), ')'
        end do
      else if (idimid == 4) then
        do j = 1 , isplit
          write (11, '(a,a,i0.2,a,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
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
      call checkncerr(istatus,__FILE__,__LINE__, 'Inquire split dimension')
      do j = 1 , isplit
        if (ldepth) then
          write (11, '(a,a,i0.2,a,i2,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                    j, trim(varname)//' ', nd , ' t,', j-1, ',z,y,x ',   &
                    trim(vardesc) , ' (', trim(varunit), ')'
        else
          write (11, '(a,a,i0.2,a,i2,a,i2,a,a,a,a,a)') trim(varname),'=>s', &
                    j, trim(varname)//' ', kz, ' t,', j-1, ',z,y,x ',   &
                    trim(vardesc) , ' (', trim(varunit), ')'
        end if
      end do
    else
      cycle
    end if
  end do

  deallocate(lvarflag)
  deallocate(dimids)

  write (11, '(a)') 'endvars'
  close(11)

  istatus = nf90_close(ncid)
  call checkncerr(istatus,__FILE__,__LINE__, 'Close file error')
  if ( lclm ) then
    istatus = nf90_close(ncid_clm)
    call checkncerr(istatus,__FILE__,__LINE__, 'Close file error')
  end if

end program ncprepare
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
