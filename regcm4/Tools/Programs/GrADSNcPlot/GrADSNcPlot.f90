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
  implicit none

  character(256) :: prgname , ncfile , tmpctl , tmpcoord
  character(512) :: command , levels
  character(32) :: lvformat , varname
  character(64) :: vardesc
  character(16) :: varunit
  character(16) :: dimdesc
  integer :: numarg , istatus , ncid

  character(256) :: charatt
  character(6) :: iproj
  real(4) :: clat , clon , plat , plon , ds , centeri , centerj
  real(4) :: minlat , minlon , maxlat , maxlon , rlatinc , rloninc
  real(4) , dimension(2) :: trlat
  real(4) , allocatable , dimension(:,:) :: xlat , xlon
  real(4) , allocatable , dimension(:) :: sigma
  real(8) , allocatable , dimension(:) :: times
  real(4) , allocatable , dimension(:,:) :: rin , rjn , ruv
  logical , allocatable , dimension(:) :: lvarflag
  integer , allocatable , dimension(:) :: dimids
  integer :: ndims , nvars , natts , udimid , totvars
  integer :: ivarid , idimid , xtype
  integer :: jxdimid , iydimid , kzdimid , itdimid
  integer :: jx , iy , kz, nt , nlat , nlon , ilat , ilon
  real(4) :: alat , alon , angle
  integer :: i

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
    write (6,*) 'Dimension jx missing'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inquire_dimension(ncid, jxdimid, len=jx)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error dimension jx'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inq_dimid(ncid, "iy", iydimid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Dimension iy missing'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inquire_dimension(ncid, iydimid, len=iy)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error dimension iy'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inq_dimid(ncid, "kz", kzdimid)
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
  istatus = nf90_inq_dimid(ncid, "time", itdimid)
  if (istatus == nf90_noerr) then
    istatus = nf90_inquire_dimension(ncid, itdimid, len=nt)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error dimension kz'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  else
    nt = 0
  end if

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

  minlat = rounder(minval(xlat),.false.)
  maxlat = rounder(maxval(xlat),.true.)
  minlon = rounder(minval(xlon),.false.)
  maxlon = rounder(maxval(xlon),.true.)

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

  rlatinc = rounder(ds/111000.0/2.0,.false.)
  rloninc = rounder(ds/111000.0/2.0,.false.)
  nlat = nint(abs(maxlat-minlat)/rlatinc)
  nlon = nint(abs(maxlon-minlon)/rloninc)
  centeri = iy/2
  centerj = jx/2
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
    istatus = nf90_get_att(ncid, nf90_global, 'latitude_of_projection_pole', &
                           plat)
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading latitude_of_projection_pole attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_att(ncid, nf90_global, 'longitude_of_projection_pole', &
                           plon)
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading longitude_of_projection_pole attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    call setup_rmc(clat,clon,centerj,centeri,ds,plon,plat)
    do ilon = 1 , nlon
      alon = minlon + (ilon-1) * rloninc
      do ilat = 1 , nlat
        alat = minlat + (ilat-1) * rlatinc
        call llij_rc(alat,alon,rin(ilon,ilat),rjn(ilon,ilat))
        ruv(ilon,ilat) = 1.0
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

  if (kz /= 0) then
    allocate(sigma(kz), stat=istatus)
    if (istatus /= 0) then
      write (6,*) 'Memory error allocating sigma'
      stop
    end if
    istatus = nf90_inq_varid(ncid, "sigma", ivarid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error : sigma variable undefined'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_var(ncid, ivarid, sigma)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading sigma variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if

    sigma = sigma * 1000.0
    write (lvformat, '(a,i4,a)') '(a,i4,a,',kz,'f7.1)'
    write (levels, lvformat) 'zdef ', kz , ' levels ', sigma
    write (11, '(a)') trim(levels)
    deallocate(sigma)
  else
    write (11, '(a)') 'zdef 1 levels 1000.0'
  end if

  if (nt /= 0) then
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
    istatus = nf90_get_var(ncid, ivarid, times)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
  else
    write (11, '(a)') 'tdef 1 linear 00Z31dec1999 1yr'
  endif

  totvars = 0
  do i = 1 , nvars
    lvarflag(i) = .false.
    istatus = nf90_inquire_variable(ncid,i,xtype=xtype,ndims=idimid)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error inquire variable ', i
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    if (idimid > 1) then
      totvars = totvars + 1
      lvarflag(i) = .true.
    end if
  end do

  write (11, '(a,i4)') 'vars ', totvars

  do i = 1 , nvars
    if (lvarflag(i) .eqv. .false.) cycle
    istatus = nf90_inquire_variable(ncid,i,name=varname,dimids=dimids)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error inquire variable ', i
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    if (dimids(2) == iydimid .and. dimids(1) == jxdimid) then
      dimdesc = ' 0 y,x'
    else if (dimids(2) == kzdimid .and. dimids(1) == iydimid) then
      write (dimdesc, '(a,i4,a)') ' ', kz, ' z,y,x'
    else if (dimids(2) == itdimid .and. dimids(1) == iydimid) then
      dimdesc = ' 0 t,y,x'
    else if (dimids(2) == itdimid .and. dimids(1) == kzdimid) then
      write (dimdesc, '(a,i4,a)') ' ', kz, ' t,z,y,x'
    else
      cycle
    endif
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

    write (11, '(a,a,a,a,a,a,a,a,a)') trim(varname),'=>',trim(varname), &
                            trim(dimdesc), ' ', trim(vardesc) , ' (',   &
                            trim(varunit), ')'
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
