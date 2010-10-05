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

program readregcm
  use mod_date
  use netcdf
  implicit none

  character(256) :: prgname , ncfile
  character(32) :: varname
  character(64) :: vardesc , timeunit
  character(16) :: varunit
  integer :: numarg , istatus , ncid

  character(256) :: charatt
  character(6) :: iproj
  real(4) :: clat , clon , plat , plon , ds
  real(4) :: minlat , minlon , maxlat , maxlon , rlatinc , rloninc
  real(4) , dimension(2) :: trlat
  real(4) , allocatable , dimension(:,:) :: xlat , xlon , var
  real(4) , allocatable , dimension(:) :: sigma
  real(8) , allocatable , dimension(:) :: times
  integer , allocatable , dimension(:) :: dimids
  integer :: ndims , nvars , natts , udimid , totvars
  integer :: ivarid , idimid , xtype
  integer :: jxdimid , iydimid , kzdimid , itdimid
  integer :: jx , iy , kz, nt , nlat , nlon
  integer :: i , j
  integer , dimension(4) :: istart , icount
#ifdef IBM
  integer , external :: iargc
#endif

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

! Open the file

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error Opening NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

! Ask general numbers inside

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error Reading NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  allocate(dimids(ndims))

! Example on how to read a global text attribute

  istatus = nf90_get_att(ncid, nf90_global, 'title', charatt)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading title attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  print *, trim(charatt)

! Read dimensions: can now be used to allocate space

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

  print *, 'Horizontal dimensions (JXxIY) : ', jx, 'x', iy

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

  print *, 'Vertical dimension (KZ)       : ', kz

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

  print *, 'Total timesteps               : ', nt

! Example on how to read not time dependent variables.

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

! Schema is ask by name and get

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

  print *, 'Projection is                 : ', iproj

! Getting more local attributes. Just to printout.

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

  print *, 'Centered on (lat,lon)         : (', clat, ',', clon, ')' 

  istatus = nf90_get_att(ncid, nf90_global, 'grid_size_in_meters', ds)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error reading grid_size_in_meters attribute'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  print *, 'Grid size is                  : ', ds, ' m'

  minlat = rounder(minval(xlat),.false.)
  maxlat = rounder(maxval(xlat),.true.)
  if (abs(minlat+90.0)<0.001 .or. abs(maxlat-90.0)<0.001) then
    minlon = -180.0
    maxlon = 180.0
  else
    minlon = rounder(minval(xlon(1,:)),.false.)
    maxlon = rounder(maxval(xlon(jx,:)),.true.)
  end if
  rlatinc = rounder(ds/111000.0/2.0,.false.)
  rloninc = rounder(ds/111000.0/2.0,.false.)
  nlat = nint(abs(maxlat-minlat)/rlatinc)
  if (minlon > 0.0 .and. maxlon < 0.0) then
    nlon = nint(abs((maxlon+360.0)-minlon)/rloninc) + 1
  else if (minlon > 0.0 .and. maxlon < 1e-30) then
    nlon = nint(360.0/rloninc) + 1.0
  else
    nlon = nint(abs(maxlon-minlon)/rloninc) + 1
  end if
  deallocate(xlat)
  deallocate(xlon)

! This is the 'window' our model has done / will dosimulation

  print *, 'Minimum latitude is           : ', minlat
  print *, 'Maximum latitude is           : ', maxlat
  print *, 'Minimum longitude is          : ', minlon
  print *, 'Maximum longitude is          : ', maxLon

  if (iproj == 'LAMCON') then
    istatus = nf90_get_att(ncid, nf90_global, 'standard_parallel', trlat)
    if ( istatus /= nf90_noerr) then
      write (6,*) 'Error reading standard_parallel attribute'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    
    print *, 'Lambert parallels are         : ', trlat

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

    print *, 'Rotated mercator pole         : (', plat, ',', plon, ')'

  end if

! Read vertical levels

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

    print *, 'Sigma coordinate values       : ', sigma
    
    deallocate(sigma)

  end if

! Read times. Look at function timeval2idate in ICBC mod_date.f90

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
    istatus = nf90_get_att(ncid, ivarid, 'units', timeunit)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_get_var(ncid, ivarid, times)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error reading time variable'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    print *, 'Time units                    : ', trim(timeunit)
    do i = 1 , nt
      print *, '    Time                      : ', &
                timeval2idate(times(i),timeunit)
    end do

    deallocate(times)

  end if

  do i = 1 , nvars
    istatus = nf90_inquire_variable(ncid,i,name=varname,ndims=idimid, &
                                    dimids=dimids)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error inquire variable ', i
      write (6,*) nf90_strerror(istatus)
      stop
    end if

    istatus = nf90_get_att(ncid,i,'long_name', vardesc)
    if (istatus /= nf90_noerr) then
      vardesc = 'Not specified'
    end if
    istatus = nf90_get_att(ncid,i,'units', varunit)
    if (istatus /= nf90_noerr) then
      varunit = 'Not specified'
    end if

    print *, 'Variable                      : ', trim(varname)
    print *, '   Description                : ', trim(vardesc)
    print *, '   Units                      : ', trim(varunit)

!   EXAMPLE READING T2M to obtain statistics for this var at second timestep

    if (varname == 't2m') then
      allocate(var(jx,iy))
      if (nt >= 2) then
        istart(4) = 2 ! Second timestep : this is record number
      else
        istart(4) = 1 ! Ok, first if second not available
      end if
      istart(3) = 1  ! This var is 4d, as to be compliant with CF it
                     ! has this 2m level dimension. Just another element in
                     ! the strt/count arrays
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1  ! Time
      icount(3) = 1  ! 2m
      icount(2) = iy ! Num in SN direction
      icount(1) = jx ! Num in WE direction

! That's it. read the variable

      istatus = nf90_get_var(ncid,i,var,istart,icount)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error Reading T2M .'
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      print *, 'T2M at timestep 2'
      print *, 'Maximum value                 : ', maxval(var)
      print *, 'Minimum value                 : ', minval(var)
      print *, 'Mean value                    : ', sum(var)/(jx*iy)
      deallocate(var)
    end if

  end do

  deallocate(dimids)

  istatus = nf90_close(ncid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error closing NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  contains

  function rounder(value,ltop)
    implicit none
    real(4) , intent(in) :: value
    logical, intent(in) :: ltop
    real(4) :: rounder
    integer :: tmpval
    if (ltop) then
      tmpval = ceiling(value*100.0)
    else
      tmpval = floor(value*100.0)
    end if
    rounder = real(tmpval)/100.0
  end function rounder

end program readregcm
