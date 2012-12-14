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

  public :: calc_rh , calc_hgt , calc_mslpres

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
    real(rk4) , dimension(im,jm) :: psmask
    real(rk4) , dimension(im,jm,kp) :: hp1
    integer(ik4) :: i , j , k , n
    do k = 1 , kp
      do j = 1 , jm
        do i = 1 , im
          if ( ps(i,j) < plev(k) ) then
            hp(i,j,k) = topo(i,j) - (tp(i,j,k)/slrate) * &
                  (1.0-exp(srovg*slrate*log(ps(i,j)/plev(k))))
          else
            hp(i,j,k) = topo(i,j) + srovg*tp(i,j,k)*log(ps(i,j)/plev(k))
          end if
        end do
      end do
    end do
    ! Now apply Gauss-Seidel filtering
    hp1(:,1,:) = hp(:,1,:)
    hp1(1,:,:) = hp(1,:,:)
    hp1(im,:,:) = hp(im,:,:)
    hp1(:,jm,:) = hp(:,jm,:)
    psmask(:,:) = 0.0
    do j = 2 , jm-1
      do i = 2 , im-1
        psmask(i,j) = (0.001*ps(i-1,j)+0.001*ps(i+1,j) + &
                       0.001*ps(i,j-1)+0.001*ps(i,j+1) - 0.004*ps(i,j))
      end do
    end do
    do n = 1 , 20
      do k = 1 , kp
        do j = 2 , jm-1
          do i = 2 , im-1
            hp1(i,j,k) = &
              ((hp1(i-1,j,k)+hp(i+1,j,k)+ &
                hp1(i,j-1,k)+hp(i,j+1,k))-psmask(i,j))/4.0
          end do
        end do
      end do
      hp(:,:,:) = hp1(:,:,:)
    end do
  end subroutine calc_hgt
!
! Extrapolate surface air temperature, extrapolate pressure at SL with
! hydrostatic equation.
!
  subroutine calc_mslpres(t,ps,ht,slp,plev,im,jm,kp)
    implicit none
    integer(ik4) , intent(in) :: im , jm , kp
    real(rk4) , dimension(im,jm,kp) , intent(in) :: t
    real(rk4) , dimension(im,jm) , intent(in) :: ht , ps
    real(rk4) , dimension(im,jm) , intent(out) :: slp
    real(rk4) , dimension(kp) , intent(in) :: plev
    integer(ik4) :: i , j , n
    real(rk4) :: tstar , hstar , alpha , sraval
    real(rk4) , dimension(im,jm) :: slp1
    real(rk4) , dimension(im,jm) :: psmask
!
    ! Follow Kallen 1996
    alpha = real(lrate*rgas/egrav)
    do j = 1 , jm
      do i = 1 , im
        tstar = t(i,j,1)*(1.0 + alpha*(ps(i,j)/plev(1)-1.0))
        if ( tstar < 255.0 ) then
          tstar = (tstar+255.0)*0.5
        else if ( tstar > 290.5 ) then
          tstar = 290.5 + (0.005*(tstar-290.5))**2
        end if
        hstar = ht(i,j)*segrav/(srgas*tstar)
        sraval = 0.5*alpha*hstar
        slp(i,j) = ps(i,j) * exp(hstar*(1.0 - sraval + (sraval*sraval)/3.0))
      end do
    end do
    ! Now apply Gauss-Seidel filtering
    slp1(:,1) = slp(:,1)
    slp1(1,:) = slp(1,:)
    slp1(im,:) = slp(im,:)
    slp1(:,jm) = slp(:,jm)
    psmask(:,:) = 0.0
    do j = 2 , jm-1
      do i = 2 , im-1
        psmask(i,j) = (0.001*ps(i-1,j)+0.001*ps(i+1,j) + &
                       0.001*ps(i,j-1)+0.001*ps(i,j+1) - 0.004*ps(i,j))
      end do
    end do
    do n = 1 , 20
      do j = 2 , jm-1
        do i = 2 , im-1
          if ( ht(i,j) > 10.0 ) then
            slp1(i,j) = &
              ((slp1(i-1,j)+slp(i+1,j)+slp1(i,j-1)+slp(i,j+1))-psmask(i,j))/4.0
          else
            slp1(i,j) = slp(i,j)
          end if
        end do
      end do
      slp(:,:) = slp1(:,:)
    end do
  end subroutine calc_mslpres
!
end module mod_derived
