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
end subroutine myabort

program timeseries

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_projections
  use mod_nchelper
  use mod_stdio
  use mod_date

  implicit none

  character(len=8) , parameter :: tasname = 'tas     '
  character(len=8) , parameter :: prname  = 'pr      '
  character(len=8) , parameter :: prcname = 'prc     '
  !character(len=8) , parameter :: tasname = 't2m     '
  !character(len=8) , parameter :: prname  = 'pr      '
  !character(len=8) , parameter :: prcname = 'prc     '

  character(len=256) :: infile , prgname , outname , line
  integer , parameter :: outunit = 100
  character(len=256) :: charatt
  character(len=64) :: timeunit , timecal
  character(len=8) :: c1 , c2
  real(rk8) :: lat , lon , iiy , jjx , centeri , centerj
  real(rk8) , dimension(2) :: trlat
  real(rk8) , dimension(:) , allocatable :: xtimes
  integer(ik4) :: istat , imodel , ifpd , iext
  integer(ik4) :: ncid , ivarid
  integer(ik4) :: jxdimid , iydimid , itdimid , it , nt
  integer(ik4) , dimension(4) :: istart , icount
  integer(ik4) :: itempid , iprecid , icprid
  real(rk8) , dimension(2,2) :: xbuf
  real(rk8) :: temp , precip , cprc
  real(rk8) :: wt1 , wt2
  type(rcm_time_and_date) :: idate
  !
  ! Read input file
  !
  call get_command_argument(0,value=prgname)
  if ( command_argument_count( ) < 3 ) then
    call usage
  end if
  call get_command_argument(1,value=infile)
  istat = nf90_open(infile, nf90_nowrite, ncid)
  call checkncerr(istat,__FILE__,__LINE__,'Error Open file '//trim(infile))
  ifpd = scan(infile, '/', .true.)
  iext = scan(infile, '.', .true.)

  istat = nf90_get_att(ncid, nf90_global, 'ipcc_scenario_code', charatt)
  if ( istat == nf90_noerr ) then
    imodel = 1
  end if

  istat = nf90_inq_dimid(ncid, "jx", jxdimid)
  if (istat /= nf90_noerr) then
    istat = nf90_inq_dimid(ncid, "x", jxdimid)
  end if
  call checkncerr(istat,__FILE__,__LINE__,'Dimension x missing')
  istat = nf90_inquire_dimension(ncid, jxdimid, len=jx)
  call checkncerr(istat,__FILE__,__LINE__,'Error inquire dimension x')
  istat = nf90_inq_dimid(ncid, "iy", iydimid)
  if (istat /= nf90_noerr) then
    istat = nf90_inq_dimid(ncid, "y", iydimid)
  end if
  call checkncerr(istat,__FILE__,__LINE__,'Dimension y missing')
  istat = nf90_inquire_dimension(ncid, iydimid, len=iy)
  call checkncerr(istat,__FILE__,__LINE__,'Error inquire dimension y')

  istat = nf90_inq_dimid(ncid, "time", itdimid)
  if (istat == nf90_noerr) then
    istat = nf90_inquire_dimension(ncid, itdimid, len=nt)
    call checkncerr(istat,__FILE__,__LINE__,'Error inquire dimension time')
  else
    nt = 0
  end if

  if ( nt > 0 ) then
    allocate(xtimes(nt), stat=istat)
    if ( istat /= 0 ) then
      write(stderr,*) 'Memory error...'
    end if
    istat = nf90_inq_varid(ncid, "time", ivarid)
    call checkncerr(istat,__FILE__,__LINE__, 'Time variable not present')
    istat = nf90_get_att(ncid, ivarid, 'units', timeunit)
    call checkncerr(istat,__FILE__,__LINE__, 'Time units not present')
    istat = nf90_get_att(ncid, ivarid, 'calendar', timecal)
    call checkncerr(istat,__FILE__,__LINE__, 'Time calendar not present')
    istat = nf90_get_var(ncid, ivarid, xtimes)
    call checkncerr(istat,__FILE__,__LINE__, 'Read time variable')
  end if

  istat = nf90_inq_varid(ncid, tasname, itempid)
  call checkncerr(istat,__FILE__,__LINE__, tasname//' variable not present')
  istat = nf90_inq_varid(ncid, prname, iprecid)
  call checkncerr(istat,__FILE__,__LINE__, prname//' variable not present')
  istat = nf90_inq_varid(ncid, prcname, icprid)
  call checkncerr(istat,__FILE__,__LINE__, prcname//' variable not present')

  istat = nf90_get_att(ncid, nf90_global, 'projection', iproj)
  call checkncerr(istat,__FILE__,__LINE__,'Error read attribute projection')

  istat = nf90_get_att(ncid, nf90_global, 'latitude_of_projection_origin', &
                         clat)
  call checkncerr(istat,__FILE__,__LINE__, &
                  'Error read attribute latitude_of_projection_origin')
  istat = nf90_get_att(ncid, nf90_global, 'longitude_of_projection_origin', &
                         clon)
  call checkncerr(istat,__FILE__,__LINE__, &
                  'Error read attribute longitude_of_projection_origin')
  istat = nf90_get_att(ncid, nf90_global, 'grid_size_in_meters', ds)
  call checkncerr(istat,__FILE__,__LINE__, &
                  'Error read attribute grid_size_in_meters')

  lon = -500.0D0
  lat = -500.0D0
  call get_command_argument(2,value=line)
  read(line,*,iostat=istat) lon
  if ( istat /= 0 .or. lon < -360.0 .and. lon > 360.0 ) then
    write (stderr,*) 'LONGITUDE NOT IN RANGE -360,360'
    call usage
  end if
  call get_command_argument(3,value=line)
  read(line,*,iostat=istat) lat
  if ( istat /= 0 .or. lat < -90.0 .and. lat > 90.0 ) then
    write (stderr,*) 'LATITUDE NOT IN RANGE -90,90'
    call usage
  end if

  write(c1,'(f8.3)') lon
  write(c2,'(f8.3)') lat
  write(outname,'(a)') infile(ifpd+1:iext-1)//'_'//trim(adjustl(c1))// &
            '_'//trim(adjustl(c2))//'.dat'

  if ( imodel == 1 ) then
    centeri = dble(iy)/2.0D0+0.5
    centerj = dble(jx)/2.0D0+0.5
  else
    centeri = dble(iy)/2.0D0
    centerj = dble(jx)/2.0D0
  end if
  !
  ! calling the map projection subroutine
  !
  if ( iproj=='LAMCON' ) then
    istat = nf90_get_att(ncid, nf90_global, 'standard_parallel', trlat)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read attribute standard_parallel')
    truelatl = trlat(1)
    truelath = trlat(2)
    call setup_lcc(clat,clon,centerj,centeri,ds,clon,truelatl,truelath)
    call llij_lc(lat,lon,jjx,iiy)
  else if ( iproj=='POLSTR' ) then
    call setup_plr(clat,clon,centerj,centeri,ds,clon)
    call llij_ps(lat,lon,jjx,iiy)
  else if ( iproj=='NORMER' ) then
    call setup_mrc(clat,clon,centerj,centeri,ds)
    call llij_mc(lat,lon,jjx,iiy)
  else if ( iproj=='ROTMER' ) then
    istat = nf90_get_att(ncid, nf90_global, 'grid_north_pole_latitude', &
                           plat)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_latitude')
    istat = nf90_get_att(ncid, nf90_global, 'grid_north_pole_longitude', &
                           plon)
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error read attribute grid_north_pole_longitude')
    call setup_rmc(clat,clon,centerj,centeri,ds,plon,plat)
    call llij_rc(lat,lon,jjx,iiy)
  else
    write (stderr,*) 'iproj = ', iproj
    write (stderr,*) 'Unrecognized or unsupported projection'
    write (stderr,*) 'Set iproj in LAMCON, POLSTR, NORMER, ROTMER'
    stop
  end if

  write(stdout,*) 'LAT = ', lat, ', LON = ', lon
  write(stdout,*) 'JX  = ', jjx, ', IY  = ', iiy

  open(unit=outunit, file=outname, form='formatted', status='replace', &
       action='write', iostat=istat)
  if ( istat /= 0 ) then
    write (stderr,*) 'Cannot open output file !'
    stop
  end if

  istart(1) = int(floor(jjx))
  istart(2) = int(floor(iiy))

  wt1 = jjx - real(istart(1))
  wt2 = iiy - real(istart(2))

  icount(1) = 2
  icount(2) = 2
  do it = 1 , nt
    istart(3) = 1
    icount(3) = 1
    istart(4) = it
    icount(4) = 1
    istat = nf90_get_var(ncid, itempid, xbuf, istart, icount)
    call checkncerr(istat,__FILE__,__LINE__, 'Read '//tasname//' variable')
    temp = bilinear()
    istart(3) = it
    icount(3) = 1
    istat = nf90_get_var(ncid, iprecid, xbuf, istart(1:3), icount(1:3))
    call checkncerr(istat,__FILE__,__LINE__, 'Read '//prname//' variable')
    precip = bilinear() * 3600.0D0 ! Change in mm/h
    istat = nf90_get_var(ncid, icprid, xbuf, istart(1:3), icount(1:3))
    call checkncerr(istat,__FILE__,__LINE__, 'Read '//prcname//' variable')
    cprc = bilinear() * 3600.0D0 ! Change in mm/h
    idate = timeval2date(xtimes(it), timeunit, timecal)
    write(outunit,'(a,4f9.3)') toiso8601(idate) , temp , precip , cprc , &
      precip-cprc
  end do

  close(outunit)
  istat = nf90_close(ncid)
  call checkncerr(istat,__FILE__,__LINE__,'Cannot close input file!')

  deallocate(xtimes)

  contains

  real(rk8) function bilinear() result(res)
    implicit none
    res = (d_one-wt2)*((d_one-wt1)*xbuf(1,1) + wt1*xbuf(2,1)) + &
          wt2*((d_one-wt1)*xbuf(1,2) + wt1*xbuf(2,2))
  end function bilinear

  subroutine usage
    implicit none
    write(stdout,*) 'Extracts temp and precip timeseries from SRF files'
    write(stdout,*) ''
    write(stdout,*) 'Usage:'
    write(stdout,*) ''
    write(stdout,*) '       ',trim(prgname),' filename.nc lon lat'
    write(stdout,*) ''
    stop
  end subroutine usage

end program timeseries

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
