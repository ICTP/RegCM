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

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  private

  real(rk4) , parameter :: rt0 = real(rtzero)
  real(rk4) , parameter :: t0 = real(tzero)
  real(rk4) , parameter :: srgas = real(rgas)
  real(rk4) , parameter :: slh0 = real(lh0)
  real(rk4) , parameter :: slh1 = real(lh1)
  real(rk4) , parameter :: slsvp1 = real(lsvp1)
  real(rk4) , parameter :: slsvp2 = real(lsvp2)
  real(rk4) , parameter :: sep2 = real(ep2)
  real(rk4) , parameter :: segrav = real(egrav)
  real(rk4) , parameter :: srovg = real(rovg)
  real(rk4) , parameter :: slrate = real(lrate)

  public :: calc_rh , calc_hgt , calc_slpres

  contains

  subroutine calc_rh(t,q,preslv,im,jm,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , kp
    real(rk4) , dimension(im,jm,kp) , intent(in) :: t
    real(rk4) , dimension(kp) , intent(in) :: preslv
    real(rk4) , dimension(im,jm,kp) , intent(inout) :: q
    integer(ik4) :: i , j , k
    real(rk4) :: hl , satvp , qs
    real(rk4) , parameter :: qmin = 0.0 ! minimum value of specific humidity
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
    integer(ik4) , intent(in) :: im , jm , kp
    real(rk4) , dimension(im,jm,kp) , intent(out) :: hp
    real(rk4) , dimension(im,jm,kp) , intent(in) :: tp
    real(rk4) , dimension(im,jm) , intent(in) :: ps , topo
    real(rk4) , dimension(kp) , intent(in) :: plev
    integer(ik4) :: i , j , k
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
    integer(ik4) , intent(in) :: im , jm , kp
    real(rk4) , dimension(im,jm,kp) , intent(in) :: t , h
    real(rk4) , dimension(im,jm) , intent(in) :: ht , ps
    real(rk4) , dimension(im,jm) , intent(out) :: slp
    real(rk4) , dimension(kp) , intent(in) :: plev
    integer(ik4) :: i , j
    real(rk4) :: tstar , saval , sraval
!
    saval = real(lrate*rgas/egrav)
    sraval = real(egrav/(lrate*rgas))
    do j = 1 , jm
      do i = 1 , im
        tstar = t(i,j,1) + saval*t(i,j,1)*(ps(j,i)/plev(1)-1.0)
        slp(i,j) = ps(i,j) * ( 1.0 + (saval*ht(i,j))/(srgas*tstar) )**sraval
      end do
    end do
  end subroutine calc_slpres
!
end module mod_derived
