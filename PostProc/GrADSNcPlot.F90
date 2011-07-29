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

program ncplot
  use netcdf
  use mod_projections
  use mod_date
  implicit none

  character(256) :: prgname , ncfile , tmpctl , tmpcoord
  character(512) :: command , levels
  character(32) :: lvformat , varname
  character(64) :: vardesc , timeunit
  character(16) :: varunit
  character(16) :: dimdesc
  integer :: numarg , istatus , ncid
  integer :: year , month , day , hour , delta , idate
  character(3) , dimension(12) :: cmon
  character(12) :: cdum

  character(256) :: charatt
  character(6) :: iproj
  real(8) :: clat , clon , plat , plon , ds , centeri , centerj
  real(4) :: minlat , minlon , maxlat , maxlon , rlatinc , rloninc
  real(8) , dimension(2) :: trlat
  real(4) , allocatable , dimension(:,:) :: xlat , xlon
  real(4) , allocatable , dimension(:) :: level
  real(8) , allocatable , dimension(:) :: times
  real(8) :: time1
  real(4) , allocatable , dimension(:,:) :: rin , rjn , ruv
  logical , allocatable , dimension(:) :: lvarflag
  integer , allocatable , dimension(:) :: dimids
  integer :: ndims , nvars , natts , udimid , totvars
  integer :: ivarid , idimid , xtype
  integer :: jxdimid , iydimid , kzdimid , itdimid , dptdimid
  integer :: jx , iy , kz , nd , nt , nlat , nlon , ilat , ilon , isplit
  real(4) :: alat , alon , angle
  integer :: i , j
  logical :: lvarsplit , lsigma , ldepth
#ifdef __PGI
  integer , external :: iargc
#endif
#ifdef IBM
  integer , external :: iargc
#endif

  data cmon /'jan','feb','mar','apr','may','jun', &
             'jul','aug','sep','oct','nov','dec'/
  data lsigma /.true./
  data ldepth /.false./
  data kzdimid  /-1/
  data itdimid  /-1/
  data dptdimid /-1/

  call getarg(0, prgname)
  numarg = iargc( )
  if (numarg < 1) then
    write (6,*) 'Not enough arguments.'
    write (6,*) ' '
    write (6,*) 'Usage : ', trim(prgname), ' Rcmfile.nc'
    write (6,*) ' '
    stop
  end if

  call getarg(1, ncfile)
  tmpctl = trim(ncfile)//'.ctl'
  tmpcoord = trim(ncfile)//'.coord'

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error Opening NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error Reading NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  allocate(lvarflag(nvars))
  allocate(dimids(ndims))

  open(11, file=tmpctl, form='formatted', status='replace')
  open(12, file=tmpcoord, form='unformatted', status='replace')
  write(11, '(a)') 'dset '//trim(ncfile)
  write(11, '(a)') 'dtype netcdf'
  write(11, '(a)') 'undef -1e+34_FillValue'

  istatus = nf90_get_att(ncid, nf90_global, 'title', charatt)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading title attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  write(11, '(a)') 'title '//trim(charatt)

  istatus = nf90_inq_dimid(ncid, "jx", jxdimid)
  if (istatus /= nf90_noerr) then
    istatus = nf90_inq_dimid(ncid, "x", jxdimid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Dimension jx missing'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  end if
  istatus = nf90_inquire_dimension(ncid, jxdimid, len=jx)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error dimension jx'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inq_dimid(ncid, "iy", iydimid)
  if (istatus /= nf90_noerr) then
    istatus = nf90_inq_dimid(ncid, "y", iydimid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Dimension iy missing'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  end if
  istatus = nf90_inquire_dimension(ncid, iydimid, len=iy)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error dimension iy'
    write (6,*) nf90_strerror(istatus)
    stop
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
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error dimension kz'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  else
    kz = 0
  end if
  istatus = nf90_inq_dimid(ncid, "depth", dptdimid)
  if (istatus == nf90_noerr) then
    ldepth = .true.
    istatus = nf90_inquire_dimension(ncid, dptdimid, len=nd)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error dimension depth'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  end if
  istatus = nf90_inq_dimid(ncid, "time", itdimid)
  if (istatus == nf90_noerr) then
    istatus = nf90_inquire_dimension(ncid, itdimid, len=nt)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error dimension time'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  else
    nt = 0
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
  if (istatus /= 0) then
    write (6,*) 'Memory error allocating xlat'
    stop
  end if
  allocate(xlon(jx,iy), stat=istatus)
  if (istatus /= 0) then
    write (6,*) 'Memory error allocating xlon'
    stop
  end if

  istatus = nf90_inq_varid(ncid, "xlat", ivarid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error : xlat variable undefined'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, xlat)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error reading xlat variable'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inq_varid(ncid, "xlon", ivarid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error : xlon variable undefined'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, xlon)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error reading xlon variable'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'projection', iproj)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading projection attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'latitude_of_projection_origin', &
                         clat)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading latitude_of_projection_origin attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'longitude_of_projection_origin', &
                         clon)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading longitude_of_projection_origin attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_get_att(ncid, nf90_global, 'grid_size_in_meters', ds)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading grid_size_in_meters attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  minlat = rounder(minval(xlat),.false.)
  maxlat = rounder(maxval(xlat),.true.)
  if (abs(minlat+90.0)<0.001 .or. abs(maxlat-90.0)<0.001) then
    minlon = -180.0
    maxlon = 180.0
  else
    if (xlon(1,1) < 0 .and. xlon(1,iy) > 0.0) then
      minlon = rounder(xlon(1,iy),.false.)
    else if (xlon(1,1) > 0 .and. xlon(1,iy) < 0.0) then
      minlon = rounder(xlon(1,1),.false.)
    else
      minlon = rounder(minval(xlon(1,:)),.false.)
    end if
    if (xlon(jx,1) > 0 .and. xlon(jx,iy) < 0.0) then
      maxlon = rounder(xlon(jx,iy),.true.)
    else if (xlon(jx,1) < 0 .and. xlon(jx,iy) > 0.0) then
      maxlon = rounder(xlon(jx,1),.true.)
    else
      maxlon = rounder(maxval(xlon(jx,:)),.true.)
    end if
  end if
  rlatinc = rounder(real(ds/111000.0/2.0),.false.)
  rloninc = rounder(real(ds/111000.0/2.0),.false.)
  nlat = nint(abs(maxlat-minlat)/rlatinc)
  if (minlon > 0.0 .and. maxlon < 0.0) then
    nlon = nint(abs((maxlon+360.0)-minlon)/rloninc) + 1
  else if (minlon > 0.0 .and. maxlon < 1e-30) then
    nlon = nint(360.0/rloninc) + 1
  else
    nlon = nint(abs(maxlon-minlon)/rloninc) + 1
  end if
  centeri = dble(iy)/2.0D0
  centerj = dble(jx)/2.0D0
  deallocate(xlat)
  deallocate(xlon)

  allocate(rin(nlon,nlat))
  allocate(rjn(nlon,nlat))
  allocate(ruv(nlon,nlat))

  if (iproj == 'LAMCON') then
    istatus = nf90_get_att(ncid, nf90_global, 'standard_parallel', trlat)
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading standard_parallel attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
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
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading grid_north_pole_latitude attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_att(ncid, nf90_global, 'grid_north_pole_longitude', &
                           plon)
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading grid_north_pole_longitude attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
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
    write (6,*) 'Unknown Projection : ', iproj
    stop
  end if

  write(11, '(a,i4,i4,a,a)') 'pdef ', jx , iy ,                         &
         ' bilin sequential binary-big ', trim(tmpcoord)
  write(11, '(a,i5,a,f7.2,f7.2)') 'xdef ', nlon , ' linear ',           &
         minlon, rloninc 
  write(11, '(a,i5,a,f7.2,f7.2)') 'ydef ', nlat , ' linear ',           &
         minlat, rlatinc 

  write(12) rin
  write(12) rjn
  write(12) ruv
  close(12)

  if (.not. ldepth .and. kz /= 0) then
    allocate(level(kz), stat=istatus)
    if (istatus /= 0) then
      write (6,*) 'Memory error allocating level'
      stop
    end if
    if (lsigma) then
      istatus = nf90_inq_varid(ncid, "sigma", ivarid)
    else
      istatus = nf90_inq_varid(ncid, "lev", ivarid)
      if (istatus /= nf90_noerr) then
        istatus = nf90_inq_varid(ncid, "plev", ivarid)
      end if
    end if
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error : level variable undefined'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_var(ncid, ivarid, level)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading level variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if

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
    if (istatus /= 0) then
      write (6,*) 'Memory error allocating times'
      stop
    end if
    istatus = nf90_inq_varid(ncid, "time", ivarid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error : time variable undefined'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    read (timeunit,'(a12,i4,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour
    istatus = nf90_get_var(ncid, ivarid, times)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    delta = nint(times(2)-times(1))
    idate = year*1000000+month*10000+day*100+hour
    call addhours(idate,idnint(times(1)))
    call split_idate(idate,year,month,day,hour)
    deallocate(times)
    if (delta == 24) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
             year , ' 1dy'
    else if (delta == 168) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
             year , ' 7dy'
    else if (delta >= 672 .and. delta <= 744) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
             year , ' 1mo'
    else if (delta > 8640) then
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
             year , ' 1yr'
    else
      write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
             'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
             year , delta, 'hr'
    end if
  else if (nt == 1) then
    istatus = nf90_inq_varid(ncid, "time", ivarid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error : time variable undefined'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    read (timeunit,'(a12,i4,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour
    istatus = nf90_get_var(ncid, ivarid, time1)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    idate = year*1000000+month*10000+day*100+hour
    call addhours(idate,idnint(time1))
    call split_idate(idate,year,month,day,hour)
    delta = 6
    write (11, '(a,i8,a,i0.2,a1,i0.2,a3,i0.4,i5,a)') &
           'tdef ', nt, ' linear ', hour, 'Z', day, cmon(month),&
           year , delta, 'hr'
  else
    write (11, '(a)') 'tdef 1 linear 00Z31dec1999 1yr'
  end if

  totvars = 0
  do i = 1 , nvars
    lvarflag(i) = .false.
    istatus = nf90_inquire_variable(ncid,i,xtype=xtype,ndims=idimid, &
                                    dimids=dimids)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error inquire variable ', i
      write (6,*) nf90_strerror(istatus)
      stop
    end if
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
          istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension splitting'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else if (idimid == 4) then
        if (dimids(4) == itdimid .and. &
            (dimids(3) == kzdimid .or. dimids(3) == dptdimid)) then
          totvars = totvars + 1
        else if (dimids(4) == itdimid .and. dimids(2) == iydimid) then
          istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension splitting'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else if (idimid == 5) then
        if (dimids(5) == itdimid .and. dimids(3) == kzdimid .and. &
            dimids(2) == iydimid) then
          istatus = nf90_inquire_dimension(ncid, dimids(4), len=isplit)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error dimension splitting'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          totvars = totvars + isplit
        else
          lvarflag(i) = .false.
        end if
      else
        lvarflag(i) = .false.
      end if
    end if
  end do

  write (11, '(a)') 'vectorpairs u,v s01u10m,s01v10m'
  write (11, '(a,i4)') 'vars ', totvars

  do i = 1 , nvars
    if (lvarflag(i) .eqv. .false.) cycle
    lvarsplit = .false.
    istatus = nf90_inquire_variable(ncid,i,name=varname,ndims=idimid, &
                                    dimids=dimids)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error inquire variable ', i
      write (6,*) nf90_strerror(istatus)
      stop
    end if
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
    istatus = nf90_get_att(ncid,i,'long_name', vardesc)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error variable ', i, ' : missing long_name attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_att(ncid,i,'units', varunit)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error variable ', i, ' : missing units attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if

    if (.not. lvarsplit) then
      write (11, '(a,a,a,a,a,a,a,a,a)') trim(varname),'=>',trim(varname), &
                              trim(dimdesc), ' ', trim(vardesc) , ' (',   &
                              trim(varunit), ')'
    else if (idimid <= 4 .and. lvarsplit) then
      istatus = nf90_inquire_dimension(ncid, dimids(3), len=isplit)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error dimension splitting'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
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
      istatus = nf90_inquire_dimension(ncid, dimids(4), len=isplit)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error dimension splitting'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
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
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error closing NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  command = 'grads -lc '//char(39)//'open '//trim(tmpctl)//char(39)

  call system(command)
  call unlink(tmpctl)
  call unlink(tmpcoord)

end program ncplot
