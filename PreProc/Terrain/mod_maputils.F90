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

  public :: lambrt , mappol , normer , rotmer , corpar , mapfac

  contains

  subroutine corpar(dlat,coriol,jx,iy)
    implicit none
    integer(ik4) , intent(in) :: iy , jx
    real(rkx) , dimension(jx,iy) , intent(in) :: dlat
    real(rkx) , dimension(jx,iy) , intent(out) :: coriol
    integer(ik4) :: i , j
    do i = 1 , iy
      do j = 1 , jx
        coriol(j,i) = eomeg2*sin(dlat(j,i)*degrad)
      end do
    end do
  end subroutine corpar

  subroutine mapfac(xlat,mapf,jx,iy,iproj)
    use mod_projections
    implicit none
    character(len=6) , intent(in) :: iproj
    integer(ik4) , intent(in) :: iy , jx
    real(rkx) , dimension(jx,iy) , intent(in) :: xlat
    real(rkx) , dimension(jx,iy) , intent(out) :: mapf
    integer(ik4) :: i , j
    if ( iproj=='LAMCON' ) then
      do i = 1 , iy
        do j = 1 , jx
          call mapfac_lc(xlat(j,i), mapf(j,i))
        end do
      end do
    else if ( iproj=='POLSTR' ) then
      do i = 1 , iy
        do j = 1 , jx
          call mapfac_ps(xlat(j,i), mapf(j,i))
        end do
      end do
    else if ( iproj=='NORMER' ) then
      do i = 1 , iy
        do j = 1 , jx
          call mapfac_mc(xlat(j,i), mapf(j,i))
        end do
      end do
    else if ( iproj=='ROTMER' ) then
      do i = 1 , iy
        do j = 1 , jx
          call mapfac_rc(real(i,rkx), mapf(j,i))
        end do
      end do
    else
      ! Not a known projection, assume regular
      mapf(:,:) = d_one
    end if
  end subroutine mapfac

  subroutine lambrt(xlon,xlat,jx,iy,clon,clat,ds,idoti,idotj,  &
                    xn,truelatl,truelath,cntri,cntrj)
    use mod_projections
    implicit none
    integer(ik4) , intent(in) :: idoti , idotj , iy , jx
    real(rkx) , intent(in) :: clat , clon , ds , truelath , truelatl
    real(rkx) , intent(in) :: cntri , cntrj
    real(rkx) , dimension(jx,iy) , intent(out) :: xlat , xlon
    real(rkx) , intent(out) :: xn
    real(rkx) :: xcntri , xcntrj
    integer(ik4) :: i , j

    xcntrj = cntrj + real(idotj,rkx)/2.0_rkx
    xcntri = cntri + real(idoti,rkx)/2.0_rkx
    call setup_lcc(clat,clon,xcntrj,xcntri,ds,clon,truelath,truelatl)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_lc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
      end do
    end do
    xn = conefac
  end subroutine lambrt

  subroutine mappol(xlon,xlat,jx,iy,clon,clat,delx,idoti,idotj,cntri,cntrj)
    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , delx , cntri , cntrj
    integer(ik4) , intent(in) :: idoti , idotj , iy , jx
    real(rkx) , intent(out) , dimension(jx,iy) :: xlat , xlon
    real(rkx) :: xcntrj , xcntri
    integer(ik4) :: i , j

    xcntrj = cntrj + real(idotj,rkx)/2.0_rkx
    xcntri = cntri + real(idoti,rkx)/2.0_rkx
    call setup_plr(clat,clon,xcntrj,xcntri,delx,clon)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_ps(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
      end do
    end do
  end subroutine mappol

  subroutine normer(xlon,xlat,jx,iy,clon,clat,delx,idoti,idotj,cntri,cntrj)
    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , delx , cntri , cntrj
    integer(ik4) , intent(in) :: idoti , idotj , iy , jx
    real(rkx) , dimension(jx,iy) , intent(out) :: xlat , xlon
    real(rkx) :: xcntri , xcntrj
    integer(ik4) :: i , j

    xcntrj = cntrj + real(idotj,rkx)/2.0_rkx
    xcntri = cntri + real(idoti,rkx)/2.0_rkx
    call setup_mrc(clat,clon,xcntrj,xcntri,delx)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_mc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
      end do
    end do
  end subroutine normer

  subroutine rotmer(xlon,xlat,jx,iy,clon,clat,pollon,   &
                    pollat,ds,idoti,idotj,cntri,cntrj)

    use mod_projections
    implicit none
    real(rkx) , intent(in) :: clat , clon , ds , pollat , pollon , cntri , cntrj
    integer(ik4) , intent(in) :: idoti , idotj , iy , jx
    real(rkx) , dimension(jx,iy) , intent(out) :: xlat , xlon
    real(rkx) :: xcntri , xcntrj
    integer(ik4) :: i , j

    xcntrj = cntrj + real(idotj,rkx)/2.0_rkx
    xcntri = cntri + real(idoti,rkx)/2.0_rkx
    call setup_rmc(clat,clon,xcntrj,xcntri,ds,pollon,pollat)
    do i = 1 , iy
      do j = 1 , jx
        call ijll_rc(real(j,rkx),real(i,rkx),xlat(j,i),xlon(j,i))
      end do
    end do
  end subroutine rotmer

end module mod_maputils

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
