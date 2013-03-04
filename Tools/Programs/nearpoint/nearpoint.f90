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
  real(rk8) :: lat , lon , iiy , jjx , jmax , imax , centeri , centerj
  integer(ik4) :: istat , ierr , imodel , iver
!
!     Read input global namelist
!
  call getarg(0, prgname)
  if ( iargc( ) < 5 ) then
    call usage
  end if
  call getarg(1, namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    call usage
  end if
  jmax = jx
  imax = iy

  call getarg(2, line)
  read(line,*,iostat=istat) iver
  if ( iver < 1 .or. iver > 2 ) then
    write (stderr,*) 'iver == 1  =>  lon/lat to jx,iy'
    write (stderr,*) 'iver == 2  =>  jx,iy to lon/lat'
    call usage
  end if

  call getarg(3, line)
  read(line,*,iostat=istat) imodel
  if ( imodel < 1 .or. imodel > 2 ) then
    write (stderr,*) 'model == 1  =>  model output'
    write (stderr,*) 'model == 2  =>  pre-proc output'
    call usage
  end if
  if ( imodel == 1 ) then
    imax = iy - 3
    jmax = jx - 3
  end if

  if ( iver == 1 ) then
    lon = -500.0D0
    lat = -500.0D0
    call getarg(4, line)
    read(line,*,iostat=istat) lon
    if ( istat /= 0 .or. lon < -360.0 .and. lon > 360.0 ) then
      write (stderr,*) 'LONGITUDE NOT IN RANGE -360,360'
      call usage
    end if
    call getarg(5, line)
    read(line,*,iostat=istat) lat
    if ( istat /= 0 .or. lat < -90.0 .and. lat > 90.0 ) then
      write (stderr,*) 'LATITUDE NOT IN RANGE -90,90'
      call usage
    end if
    write(stdout,*) 'LON = ',lon,', LAT = ',lat
  else
    jjx = -1
    iiy = -1
    call getarg(4, line)
    read(line,*,iostat=istat) jjx
    if ( jjx < 1 .or. jjx > jmax ) then
      write (stderr,*) 'JX NOT IN RANGE 1 to ',jmax
      call usage
    end if
    call getarg(5, line)
    read(line,*,iostat=istat) iiy
    if ( iiy < 1 .or. iiy > imax ) then
      write (stderr,*) 'IY NOT IN RANGE 1 to ',imax
      call usage
    end if
    write(stdout,*) 'JX  = ',jjx,', IY  = ',iiy
  end if
  if ( imodel == 1 ) then
    centeri = dble(iy-3)/2.0D0+0.5
    centerj = dble(jx-3)/2.0D0+0.5
  else
    centeri = dble(iy)/2.0D0
    centerj = dble(jx)/2.0D0
  end if
  ds = ds*1000.0
!
!-----calling the map projection subroutine
!
  if ( iproj=='LAMCON' ) then
    call setup_lcc(clat,clon,centerj,centeri,ds,clon,truelatl,truelath)    
    if ( iver == 1 ) then
      call llij_lc(lat,lon,jjx,iiy)
    else
      call ijll_lc(jjx,iiy,lat,lon)
    end if
  else if ( iproj=='POLSTR' ) then
    call setup_plr(clat,clon,centerj,centeri,ds,clon)
    if ( iver == 1 ) then
      call llij_ps(lat,lon,jjx,iiy)
    else
      call ijll_ps(jjx,iiy,lat,lon)
    end if
  else if ( iproj=='NORMER' ) then
    call setup_mrc(clat,clon,centerj,centeri,ds)
    if ( iver == 1 ) then
      call llij_mc(lat,lon,jjx,iiy)
    else
      call ijll_mc(jjx,iiy,lat,lon)
    end if
  else if ( iproj=='ROTMER' ) then
    call setup_rmc(clat,clon,centerj,centeri,ds,plon,plat)
    if ( iver == 1 ) then
      call llij_rc(lat,lon,jjx,iiy)
    else
      call ijll_rc(jjx,iiy,lat,lon)
    end if
  else
    write (stderr,*) 'iproj = ', iproj
    write (stderr,*) 'Unrecognized or unsupported projection'
    write (stderr,*) 'Set iproj in LAMCON, POLSTR, NORMER, ROTMER'
    stop
  end if
  if ( iver == 1 ) then
    write(stdout,*) 'JX  = ', jjx, ', IY  = ', iiy
  else
    write(stdout,*) 'LON = ', lon, ', LAT = ', lat
  end if
!
  contains

  subroutine usage
    implicit none

    write(stdout,*) 'Usage:'
    write(stdout,*) ''
    write(stdout,*) '       ',trim(prgname),' regcm.in iver model c1 c2'
    write(stdout,*) ''
    write(stdout,*) 'Where :'
    write(stdout,*) '        iver : 1 from lon,lat to jx,iy'
    write(stdout,*) '               2 from jx,iy to lon,lat'
    write(stdout,*) '       model : 1 model output'
    write(stdout,*) '               2 pre-processing output'
    write(stdout,*) '          c1 : coordinate W-E (lon or jx)'
    write(stdout,*) '          c2 : coordinate S-N (lat or iy)'
    write(stdout,*) ''
    stop
  end subroutine usage

end program nearpoint
