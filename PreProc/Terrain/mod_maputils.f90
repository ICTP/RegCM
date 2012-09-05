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

module mod_maputils

  use mod_constants

  contains

  subroutine lambrt(xlon,xlat,smap,coriol,iy,jx,clon,clat,ds,idot,  &
                    xn,truelatl,truelath)

  use mod_projections
  implicit none
!
  real(rk8) :: clat , clon , ds , truelath , truelatl , xn
  integer(ik4) :: idot , iy , jx
  real(rk4) , dimension(iy,jx) :: coriol , smap , xlat , xlon
  intent (in) clat , clon , ds , idot , iy , jx , truelath ,        &
              truelatl
  intent (out) coriol , smap , xlat , xlon , xn
!
  real(rk8) :: cntri , cntrj
  integer(ik4) :: i , j
!
  cntrj = dble(jx+idot)/d_two
  cntri = dble(iy+idot)/d_two
  call setup_lcc(clat,clon,cntrj,cntri,ds,clon,truelath,truelatl)
!
  do j = 1 , jx
    do i = 1 , iy
      call ijll_lc(real(j),real(i),xlat(i,j),xlon(i,j))
      call mapfac_lc(xlat(i,j), smap(i,j))
    end do
  end do
  if ( idot==1 ) then
    do i = 1 , iy
      do j = 1 , jx
        coriol(i,j) = real(eomeg2*dsin(dble(xlat(i,j))*degrad))
      end do
    end do
  end if
  xn = conefac
  end subroutine lambrt

  subroutine mappol(xlon,xlat,xmap,coriol,iy,jx,clon,clat,delx,idot)

  use mod_projections
  implicit none
!
  real(rk8) :: clat , clon , delx
  integer(ik4) :: idot , iy , jx
  real(rk4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
  intent (in) clat , clon , delx , idot , iy , jx
  intent (out) coriol , xlat , xlon , xmap
!
  real(rk8) :: cntrj , cntri
  integer(ik4) :: i , j
!
  cntrj = dble(jx+idot)/d_two
  cntri = dble(iy+idot)/d_two
  call setup_plr(clat,clon,cntrj,cntri,delx,clon)
!
  do i = 1 , iy
    do j = 1 , jx
      call ijll_ps(real(j),real(i),xlat(i,j),xlon(i,j))
      call mapfac_ps(xlat(i,j), xmap(i,j))
    end do
  end do
 
  if ( idot==1 ) then
    do i = 1 , iy
      do j = 1 , jx
        coriol(i,j) = real(eomeg2*dsin(dble(xlat(i,j))*degrad))
      end do
    end do
  end if
  end subroutine mappol

  subroutine normer(xlon,xlat,xmap,coriol,iy,jx,clon,clat,delx,idot)

  use mod_projections
  implicit none
!
  real(rk8) :: clat , clon , delx
  integer(ik4) :: idot , iy , jx
  real(rk4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
  intent (in) clat , clon , delx , idot , iy , jx
  intent (out) coriol , xlat , xlon , xmap
!
  real(rk8) :: cntri , cntrj
  integer(ik4) :: i , j
!
  cntrj = dble(jx+idot)/d_two
  cntri = dble(iy+idot)/d_two
  call setup_mrc(clat,clon,cntrj,cntri,delx)
!
  do i = 1 , iy
    do j = 1 , jx
      call ijll_mc(real(j),real(i),xlat(i,j),xlon(i,j))
      call mapfac_mc(xlat(i,j), xmap(i,j))
    end do
  end do
 
  if ( idot==1 ) then
    do i = 1 , iy
      do j = 1 , jx
        coriol(i,j) = real(eomeg2*dsin(dble(xlat(i,j))*degrad))
      end do
    end do
  end if
 
  end subroutine normer

  subroutine rotmer(xlon,xlat,xmap,coriol,iy,jx,clon,clat,pollon,   &
                    pollat,ds,idot)

  use mod_projections
  implicit none
!
  real(rk8) :: clat , clon , ds , pollat , pollon
  integer(ik4) :: idot , iy , jx
  real(rk4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
  intent (in) clat , clon , ds , idot , iy , jx
  intent (out) coriol , xlat , xlon , xmap
!
  real(rk8) :: cntri , cntrj
  integer(ik4) :: i , j
!
  cntrj = dble(jx+idot)/d_two
  cntri = dble(iy+idot)/d_two
  call setup_rmc(clat,clon,cntrj,cntri,ds,pollon,pollat)
!
  do i = 1 , iy
    do j = 1 , jx
      call ijll_rc(real(j),real(i),xlat(i,j),xlon(i,j))
      call mapfac_rc(real(i), xmap(i,j))
    end do
  end do
 
  if ( idot==1 ) then
    do i = 1 , iy
      do j = 1 , jx
        coriol(i,j) = real(eomeg2*dsin(dble(xlat(i,j))*degrad))
      end do
    end do
  end if
 
  end subroutine rotmer
!
end module mod_maputils
