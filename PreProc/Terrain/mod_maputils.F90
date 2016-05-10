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
  use mod_intkinds
  use mod_realkinds

  private

  public :: lambrt , mappol , normer , rotmer

  contains

  subroutine lambrt(xlon,xlat,smap,coriol,jx,iy,clon,clat,ds,idot,  &
                    xn,truelatl,truelath)
    use mod_projections
    implicit none
    integer(ik4) , intent(in) :: idot , iy , jx
    real(rkx) , intent(in) :: clat , clon , ds , truelath , truelatl
    real(rkx) , dimension(jx,iy) , intent(out) :: coriol , smap , xlat , xlon
    real(rkx) , intent(out) :: xn
    real(rkx) :: cntri , cntrj
    integer(ik4) :: i , j

    cntrj = real(jx+idot,rkx)/d_two
    cntri = real(iy+idot,rkx)/d_two
    call setup_lcc(clat,clon,cntrj,cntri,ds,clon,truelath,truelatl)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_lc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
        call mapfac_lc(xlat(j,i), smap(j,i))
      end do
    end do
    if ( idot==1 ) then
      do i = 1 , iy
        do j = 1 , jx
          coriol(j,i) = eomeg2*sin(xlat(j,i)*degrad)
        end do
      end do
    end if
    xn = conefac
  end subroutine lambrt

  subroutine mappol(xlon,xlat,xmap,coriol,jx,iy,clon,clat,delx,idot)
    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , delx
    integer(ik4) , intent(in) :: idot , iy , jx
    real(rkx) , intent(out) , dimension(jx,iy) :: coriol , xlat , xlon , xmap
    real(rkx) :: cntrj , cntri
    integer(ik4) :: i , j

    cntrj = real(jx+idot,rkx)/d_two
    cntri = real(iy+idot,rkx)/d_two
    call setup_plr(clat,clon,cntrj,cntri,delx,clon)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_ps(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
        call mapfac_ps(xlat(j,i), xmap(j,i))
      end do
    end do
    if ( idot==1 ) then
      do i = 1 , iy
        do j = 1 , jx
          coriol(j,i) = eomeg2*sin(xlat(j,i)*degrad)
        end do
      end do
    end if
  end subroutine mappol

  subroutine normer(xlon,xlat,xmap,coriol,jx,iy,clon,clat,delx,idot)
    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , delx
    integer(ik4) , intent(in) :: idot , iy , jx
    real(rkx) , dimension(jx,iy) , intent(out) :: coriol , xlat , xlon , xmap
    real(rkx) :: cntri , cntrj
    integer(ik4) :: i , j

    cntrj = real(jx+idot,rkx)/d_two
    cntri = real(iy+idot,rkx)/d_two
    call setup_mrc(clat,clon,cntrj,cntri,delx)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_mc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
        call mapfac_mc(xlat(j,i), xmap(j,i))
      end do
    end do
    if ( idot==1 ) then
      do i = 1 , iy
        do j = 1 , jx
          coriol(j,i) = eomeg2*sin(xlat(j,i)*degrad)
        end do
      end do
    end if
  end subroutine normer

  subroutine rotmer(xlon,xlat,xmap,coriol,jx,iy,clon,clat,pollon,   &
                    pollat,ds,idot)

    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , ds , pollat , pollon
    integer(ik4) , intent(in) :: idot , iy , jx
    real(rkx) , dimension(jx,iy) , intent(out) :: coriol , xlat , xlon , xmap
    real(rkx) :: cntri , cntrj
    integer(ik4) :: i , j

    cntrj = real(jx+idot,rkx)/d_two
    cntri = real(iy+idot,rkx)/d_two
    call setup_rmc(clat,clon,cntrj,cntri,ds,pollon,pollat)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_rc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
        call mapfac_rc(real(i,rkx), xmap(j,i))
      end do
    end do
    if ( idot == 1 ) then
      do i = 1 , iy
        do j = 1 , jx
          coriol(j,i) = eomeg2*sin(xlat(j,i)*degrad)
        end do
      end do
    end if
  end subroutine rotmer

end module mod_maputils

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
