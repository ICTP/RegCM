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

program nearpoint

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_projections
  use mod_stdio

  implicit none
!
  character(len=256) :: namelistfile , prgname , outname , line
  real(rk8) :: lat , lon , i , j , centeri , centerj
  integer(ik4) :: istat , ierr
!
!     Read input global namelist
!
  call getarg(0, prgname)
  if ( iargc( ) < 3 ) then
    call usage
  end if
  call getarg(1, namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    call usage
  end if

  lon = -500.0D0
  lat = -500.0D0
  call getarg(2, line)
  read(line,fmt='(f)',iostat=istat) lon
  if ( istat /= 0 .or. lon < -360.0 .and. lon > 360.0 ) then
    write (stderr,*) 'LONGITUDE NOT IN RANGE -360,360'
    call usage
  end if
  call getarg(3, line)
  read(line,fmt='(f)',iostat=istat) lat
  if ( istat /= 0 .or. lat < -90.0 .and. lat > 90.0 ) then
    write (stderr,*) 'LATITUDE NOT IN RANGE -90,90'
    call usage
  end if
  centeri = dble(iy)/2.0D0
  centerj = dble(jx)/2.0D0
  ds = ds*1000.0

  write(stdout,*) 'LON = ',lon,', LAT = ',lat
!
!-----calling the map projection subroutine
!
  if ( iproj=='LAMCON' ) then
    call setup_lcc(clat,clon,centerj,centeri,ds,clon,truelatl,truelath)    
    call llij_lc(lat,lon,i,j)
  else if ( iproj=='POLSTR' ) then
    call setup_plr(clat,clon,centerj,centeri,ds,clon)
    call llij_ps(lat,lon,i,j)
  else if ( iproj=='NORMER' ) then
    call setup_mrc(clat,clon,centerj,centeri,ds)
    call llij_mc(lat,lon,i,j)
  else if ( iproj=='ROTMER' ) then
    call setup_rmc(clat,clon,centerj,centeri,ds,plon,plat)
    call llij_rc(lat,lon,i,j)
  else
    write (stderr,*) 'iproj = ', iproj
    write (stderr,*) 'Unrecognized or unsupported projection'
    write (stderr,*) 'Set iproj in LAMCON, POLSTR, NORMER, ROTMER'
    stop
  end if
  write(stdout,*) 'Geo mapping done: J = ', j, ', I = ', i
!
  contains

  subroutine usage
    implicit none

    write(stdout,*) 'Usage:'
    write(stdout,*) ''
    write(stdout,*) '       ',trim(prgname),' regcm.in lon lat'
    write(stdout,*) ''
    stop
  end subroutine usage

end program nearpoint
