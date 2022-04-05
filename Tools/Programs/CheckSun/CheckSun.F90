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
  call abort
endsubroutine myabort

program checksun
  use mod_date
  use mod_constants
  use netcdf
  implicit none

  character(256) :: chararg , ncfile
  integer :: numarg , istatus , ncid
  real(4) , allocatable , dimension(:,:) :: xlat , xlon , solin
  integer , dimension(3) :: idims , istart , icount
  integer , dimension(1) :: ixtime
  integer :: ivarid , itimid , ilonid , ilatid
  integer(8) :: idate0 , idate1 , idate2
  integer :: ifrq
  integer :: year , month , day , hour
  type(rcm_time_and_date) :: xidate0 , xidate1 , xidate2 , xidate
  type(rcm_time_interval) :: tdif
  integer :: jxdimid , iydimid
  integer :: jx , iy
  integer :: it , nt , ibase
  real(8) :: xtime , gmt
  character(32) :: csdate , calendar
  character(256) :: ofname
  real(8) , parameter :: dayspy = 365.2422D0
  real(8) , parameter :: solcon = 1367.0D0

  call get_command_argument(0,value=chararg)
  numarg = command_argument_count( )
  if (numarg < 5) then
    write (6,*) 'Not enough arguments.'
    write (6,*) ' '
    write (6,*) 'Usage : ', trim(chararg), &
              ' XXXX_DOMAIN000.nc CALENDAR IDATE0 IDATE1 IDATE2 IFRQ'
    write (6,*) ' '
    write (6,*) 'Where IDATEs are dates in the format YYYYMMDDHH,', &
                ' and IFRQ is delta in hours'
    write (6,*) 'to calculate SOLIN'
    write (6,*) ' '
    write (6,*) 'Example:'
    write (6,*) trim(chararg), &
             ' EUROPE_DOMAIN000.nc gregorian 1999010100 2001010100 2010010100 3'
    write (6,*) ' '
    stop
  end if

  call get_command_argument(1,value=ncfile)
  call get_command_argument(2,value=calendar)
  call get_command_argument(3,value=chararg)
  read(chararg, '(i10)',iostat=istatus) idate0
  if (istatus /= 0) then
    print *, 'Cannot parse idate0'
    stop
  end if
  xidate0 = i8wcal(idate0,calendar_int(calendar))
  call get_command_argument(4,value=chararg)
  read(chararg, '(i10)',iostat=istatus) idate1
  if (istatus /= 0) then
    print *, 'Cannot parse idate1'
    stop
  end if
  xidate1 = i8wcal(idate1,calendar_int(calendar))
  call get_command_argument(5,value=chararg)
  read(chararg, '(i10)',iostat=istatus) idate2
  if (istatus /= 0) then
    print *, 'Cannot parse idate2'
    stop
  end if
  xidate2 = i8wcal(idate2,calendar_int(calendar))
  call get_command_argument(6,value=chararg)
  read(chararg, '(i5)',iostat=istatus) ifrq
  if (istatus /= 0) then
    print *, 'Cannot parse ifrq'
    stop
  end if

  call split_idate(xidate0,year,month,day,hour)
  gmt = dble(hour)

! Open the file

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  if ( istatus /= nf90_noerr) then
    write (6,*) 'Error Opening NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

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
  allocate(solin(jx,iy), stat=istatus)
  if (istatus /= 0) then
    write (6,*) 'Memory error allocating solin'
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

  istatus = nf90_close(ncid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error closing NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  write (ofname,'(a,i10,a)') 'solin_',toint10(xidate0),'.nc'
  istatus = nf90_create(ofname, nf90_clobber, ncid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create NetCDF file '//ofname
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create dimension IY'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create dimension JX'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create dimension time'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_def_var(ncid, 'time', nf90_int, idims(3:3), itimid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create variable time'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_put_att(ncid, itimid, 'units', &
                         'hours since '//tochar(xidate0))
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error adding time units'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), ilatid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create variable xlat'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_put_att(ncid, ilatid, 'units', 'degrees_north')
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error adding xlat units'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), ilonid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create variable xlon'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_put_att(ncid, ilonid, 'units', 'degrees_east')
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error adding xlon units'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_def_var(ncid, 'solin', nf90_float, idims, ivarid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error create variable solin'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_put_att(ncid, ivarid, 'units', 'W/m^2')
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error adding solin units'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_enddef(ncid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error end of definitions'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  xtime = 0.0D0
  xidate = xidate1

  tdif = xidate2-xidate1
  nt = tohours(tdif)/ifrq+1
  tdif = xidate1-xidate0
  ibase =tohours(tdif)/ifrq

  istatus = nf90_put_var(ncid,ilatid,xlat)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error variable xlat write'
    write (6,*) nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_put_var(ncid,ilonid,xlon)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error variable xlon write'
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  istart(1) = 1
  istart(2) = 1
  icount(1) = jx
  icount(2) = iy
  icount(3) = 1
  do it = 1 , nt
    call calcsolin
    istart(3) = it
    print *, 'Doing ', tochar(xidate)
    ixtime(1) = (it-1)*ifrq
    istatus = nf90_put_var(ncid,itimid,ixtime,istart(3:3),icount(3:3))
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error variable time write'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    istatus = nf90_put_var(ncid,ivarid,solin,istart,icount)
    if (istatus /= nf90_noerr) then
      write (6,*) 'Error variable solin write'
      write (6,*) nf90_strerror(istatus)
      stop
    end if
    tdif = rcm_time_interval(ifrq,uhrs)
    xidate = xidate + tdif
  end do

  istatus = nf90_close(ncid)
  if (istatus /= nf90_noerr) then
    write (6,*) 'Error closing NetCDF file ', trim(ncfile)
    write (6,*) nf90_strerror(istatus)
    stop
  end if

  contains

  subroutine calcsolin
    implicit none
    integer :: i , j , idiff
    real(8) :: eccf , theta , calday , xt24 , tlocap , omga , xxlat , coszrs
    real(8) :: delta , decdeg, xday , lhour
    call split_idate(xidate,year,month,day,hour)
    lhour = dble(hour)
    tdif = xidate-xidate0
    idiff = (tohours(tdif)/ifrq)-ibase
    xday = dble(idiff)/24.0D0
    calday = yeardayfrac(xidate)
    theta = twopi*calday/dayspy
    delta = 0.006918D0 - 0.399912D0*dcos(theta) + &
            0.070257D0*dsin(theta) -              &
            0.006758D0*dcos(2.0D0*theta) +        &
            0.000907D0*dsin(2.0D0*theta) -        &
            0.002697D0*dcos(3.0D0*theta) +        &
            0.001480D0*dsin(3.0D0*theta)
    decdeg = delta/degrad
    eccf = 1.000110D0 + 0.034221D0*dcos(theta) +  &
           0.001280D0 * dsin(theta) +             &
           0.000719D0 * dcos(d_two*theta) +       &
           0.000077D0 * dsin(d_two*theta)
    xt24 = dmod(lhour*minph+xtime,minpd)
    do i = 1 , iy
      do j = 1 , jx
        tlocap = xt24/minph + xlon(j,i)/15.0D0
        tlocap = dmod(tlocap+houpd,houpd)
        omga = 15.0D0*(tlocap-12.0D0)*degrad
        xxlat = xlat(j,i)*degrad
        coszrs = dsin(delta)*dsin(xxlat) +           &
                 dcos(delta)*dcos(xxlat)*dcos(omga)
        coszrs = dmax1(0.0D0,coszrs)
        coszrs = dmin1(1.0D0,coszrs)
        solin(j,i) = real(solcon*eccf*coszrs)
      end do
    end do
    xtime = xtime + dble(ifrq)
    if (xtime >= 23.999D0) xtime = 0.0D0
  end subroutine calcsolin

end program checksun
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
