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

module mod_derived

  use mod_constants
  use mod_realkinds

  private

  real(sp) , parameter :: rt0 = real(rtzero)
  real(sp) , parameter :: t0 = real(tzero)
  real(sp) , parameter :: srgas = real(rgas)
  real(sp) , parameter :: slh0 = real(lh0)
  real(sp) , parameter :: slh1 = real(lh1)
  real(sp) , parameter :: slsvp1 = real(lsvp1)
  real(sp) , parameter :: slsvp2 = real(lsvp2)
  real(sp) , parameter :: sep2 = real(ep2)
  real(sp) , parameter :: segrav = real(egrav)
  real(sp) , parameter :: srovg = real(rovg)
  real(sp) , parameter :: slrate = real(lrate)

  public :: calc_rh , calc_hgt , calc_slpres

  contains

  subroutine calc_rh(t,q,preslv,im,jm,kp)
    implicit none
    integer , intent(in) :: im , jm , kp
    real(sp) , dimension(im,jm,kp) , intent(in) :: t
    real(sp) , dimension(kp) , intent(in) :: preslv
    real(sp) , dimension(im,jm,kp) , intent(inout) :: q
    integer :: i , j , k
    real(sp) :: hl , satvp , qs
    real(sp) , parameter :: qmin = 0.0 ! minimum value of specific humidity
    do k = 1 , kp
      do j = 1 , jm
        do i = 1 , im
          hl = slh0-slh1*(t(i,j,k)-t0)
          satvp = slsvp1*exp(slsvp2*hl*(rt0-1.0/t(i,j,k)))
          qs = sep2*satvp/(preslv(k)-satvp)  ! preslv (hpa)
          q(i,j,k) = amin1(amax1(q(i,j,k)/qs,qmin),1.0)*100.0
        end do
      end do
    end do
  end subroutine calc_rh
!
! Just use hydrostatic equation for pressure grater than surface p
! Extrapolation is performed below model surface.
!
  subroutine calc_hgt(hp,tp,ps,topo,plev,im,jm,kp)
    implicit none
    integer , intent(in) :: im , jm , kp
    real(sp) , dimension(im,jm,kp) , intent(out) :: hp
    real(sp) , dimension(im,jm,kp) , intent(in) :: tp
    real(sp) , dimension(im,jm) , intent(in) :: ps , topo
    real(sp) , dimension(kp) , intent(in) :: plev
    integer :: i , j , k
    do k = 1 , kp
      do j = 1 , jm
        do i = 1 , im
          if ( ps(i,j) < plev(k) ) then
            hp(i,j,k) = topo(i,j) - (tp(i,j,k)/slrate) * &
                  (1.0-exp(-srovg*slrate*log(plev(k)/ps(i,j))))
          else
            hp(i,j,k) = topo(i,j) + srovg*tp(i,j,k)*log(ps(i,j)/plev(k))
          end if
        end do
      end do
    end do
  end subroutine calc_hgt
!
! Extrapolate surface air temperature, extrapolate pressure at SL with
! hydrostatic equation.
!
  subroutine calc_slpres(h,t,ps,ht,slp,plev,im,jm,kp)
    implicit none
    integer , intent(in) :: im , jm , kp
    real(sp) , dimension(im,jm,kp) , intent(in) :: t , h
    real(sp) , dimension(im,jm) , intent(in) :: ht , ps
    real(sp) , dimension(im,jm) , intent(out) :: slp
    real(sp) , dimension(kp) , intent(in) :: plev
    integer :: i , j , kbc
    real(sp) :: tsfc
!
    do j = 1 , jm
      do i = 1 , im
        kbc = 1
        do while ( plev(kbc) >= ps(i,j) )
          kbc = kbc + 1
        end do
        tsfc = t(i,j,kbc)+slrate*(h(i,j,kbc)-ht(i,j))
        slp(i,j) = ps(i,j) * &
                   exp(-segrav/(srgas*slrate)*log(1.0-ht(i,j)*slrate/tsfc))
      end do
    end do
  end subroutine calc_slpres
!
end module mod_derived
