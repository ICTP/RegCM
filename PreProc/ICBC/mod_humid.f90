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

module mod_humid

  use mod_intkinds
  use mod_realkinds
  use mod_constants
!
  contains

  subroutine humid1(t,q,ps,ptop,sigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) :: ps , ptop
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rk8) , intent(in) , dimension(nk) :: sigma
!
    real(rk8) :: lh , p , qs , satvp
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = (ptop+sigma(k)*ps)*d_10       ! PRESSURE AT LEVEL K
          lh = lh0 - lh1*(t(i,j,k)-tzero)
          satvp = lsvp1*dexp(lsvp2*lh*(d_one/tzero-d_one/t(i,j,k)))
          qs = ep2*satvp/(p-satvp)        ! SAT. MIXING RATIO
          q(i,j,k) = max(q(i,j,k)/qs,d_zero)
        end do
      end do
    end do
  end subroutine humid1
!
!-----------------------------------------------------------------------
!
  subroutine humid1_o(t,q,ps,sigma,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(im,jm) :: ps
    real(rk8) , intent(in) , dimension(im,jm,km) :: t
    real(rk8) , intent(in) , dimension(km) :: sigma
    real(rk8) , intent(inout) , dimension(im,jm,km) :: q
!
    real(rk8) :: lh , p , qs , satvp
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    ! DATA ON SIGMA LEVELS
    !
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          p = sigma(k)*(ps(i,j)-ptop) + ptop
          lh = lh0 - lh1*(t(i,j,k)-tzero)       ! LATENT HEAT OF EVAP.
          satvp = lsvp1*dexp(lsvp2*lh*(rtzero-d_one/t(i,j,k)))
                                                ! SATURATION VAP PRESS.
          qs = ep2*satvp/(p-satvp)              ! SAT. MIXING RATIO
          q(i,j,k) = max(q(i,j,k)/qs,d_zero)
        end do
      end do
    end do
  end subroutine humid1_o
!
!-----------------------------------------------------------------------
!
  subroutine humid1fv(t,q,p3d,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: p3d , t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q
!
    real(rk8) :: lh , qs , satvp
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          if ( p3d(i,j,k) > -9990.0D0 ) then
            lh = lh0 - lh1*(t(i,j,k)-tzero)  ! LATENT HEAT OF EVAP.
            satvp = lsvp1*dexp(lsvp2*lh*(rtzero-d_one/t(i,j,k)))
                                                 ! SATURATION VAP PRESS.
            qs = ep2*satvp/(p3d(i,j,k)-satvp)   ! SAT. MIXING RATIO
            q(i,j,k) = max(q(i,j,k)/qs,d_zero)    !ALREADY MIXING RATIO
          else
            q(i,j,k) = -9999.D0
          end if
        end do
      end do
    end do
  end subroutine humid1fv
!
!-----------------------------------------------------------------------
!
  subroutine humid2(t,q,ps,ptop,sigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(ni,nj) :: ps
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rk8) , intent(in) , dimension(nk) :: sigma
!
    real(rk8) :: lh , p , qs , satvp
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = (ptop+sigma(k)*ps(i,j))*d_10
          lh = lh0 - lh1*(t(i,j,k)-tzero)
          satvp = lsvp1*dexp(lsvp2*lh*(rtzero-d_one/t(i,j,k)))
          qs = ep2*satvp/(p-satvp)
          q(i,j,k) = max(q(i,j,k)*qs,d_zero)
        end do
      end do
    end do
  end subroutine humid2
!
end module mod_humid
