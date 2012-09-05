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
  real(rk4) , parameter :: rt0 = real(rtzero)
  real(rk4) , parameter :: t0 = real(tzero)
  real(rk4) , parameter :: slh0 = real(lh0)
  real(rk4) , parameter :: slh1 = real(lh1)
  real(rk4) , parameter :: slsvp1 = real(lsvp1)
  real(rk4) , parameter :: slsvp2 = real(lsvp2)
  real(rk4) , parameter :: sep2 = real(ep2)

  contains

  subroutine humid1(t,q,ps,ptop,sigma,ni,nj,nk)
  implicit none
!
  integer(ik4) :: ni , nj , nk
  real(rk4) :: ps
  real(rk8) :: ptop
  real(rk4) , dimension(ni,nj,nk) :: q , t
  real(rk4) , dimension(nk) :: sigma
  intent (in) ni , nj , nk , ps , ptop , sigma , t
  intent (inout) q
!
  real(rk4) :: lh , p , qs , pt , satvp
  integer(ik4) :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
  pt = real(ptop)
  do i = 1 , ni
    do j = 1 , nj
      do k = 1 , nk
        p = (pt+sigma(k)*ps)*10.        ! PRESSURE AT LEVEL K
        lh = slh0 - slh1*(t(i,j,k)-t0)
        satvp = slsvp1*exp(slsvp2*lh*(1./t0-1./t(i,j,k)))
        qs = sep2*satvp/(p-satvp)        ! SAT. MIXING RATIO
        q(i,j,k) = amax1(q(i,j,k)/qs,0.0)
      end do
    end do
  end do
  end subroutine humid1
!
!-----------------------------------------------------------------------
!
  subroutine humid1_o(t,q,ps,sigma,ptop,im,jm,km)
  implicit none
!
  integer(ik4) :: im , jm , km
  real(rk8) :: ptop
  real(rk4) , dimension(im,jm) :: ps
  real(rk4) , dimension(im,jm,km) :: q , t
  real(rk4) , dimension(km) :: sigma
  intent (in) im , jm , km , ps , ptop , sigma , t
  intent (inout) q
!
  real(rk4) :: lh , p , qs , satvp , pt
  integer(ik4) :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!     DATA ON SIGMA LEVELS
!
  pt = real(ptop)
  do k = 1 , km
    do j = 1 , jm
      do i = 1 , im
        p = sigma(k)*(ps(i,j)-pt) + pt
        lh = slh0 - slh1*(t(i,j,k)-t0)       ! LATENT HEAT OF EVAP.
        satvp = slsvp1*exp(slsvp2*lh*(rt0-1./t(i,j,k)))
                                                  ! SATURATION VAP PRESS.
        qs = sep2*satvp/(p-satvp)                 ! SAT. MIXING RATIO
        q(i,j,k) = amax1(q(i,j,k)/qs,0.0)
      end do
    end do
  end do
  end subroutine humid1_o
!
!-----------------------------------------------------------------------
!
  subroutine humid1fv(t,q,p3d,ni,nj,nk)
  implicit none
!
  integer(ik4) :: ni , nj , nk
  real(rk4) , dimension(ni,nj,nk) :: p3d , q , t
  intent (in) ni , nj , nk , p3d , t
  intent (inout) q
!
  real(rk4) :: lh , qs , satvp
  integer(ik4) :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
  do i = 1 , ni
    do j = 1 , nj
      do k = 1 , nk
        if ( p3d(i,j,k) > -9990. ) then
          lh = slh0 - slh1*(t(i,j,k)-t0)  ! LATENT HEAT OF EVAP.
          satvp = slsvp1*exp(slsvp2*lh*(rt0-1./t(i,j,k)))
                                               ! SATURATION VAP PRESS.
          qs = sep2*satvp/(p3d(i,j,k)-satvp)   ! SAT. MIXING RATIO
          q(i,j,k) = amax1(q(i,j,k)/qs,0.0)    !ALREADY MIXING RATIO
        else
          q(i,j,k) = -9999.
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
!
  integer(ik4) :: ni , nj , nk
  real(rk8) :: ptop
  real(rk4) , dimension(ni,nj) :: ps
  real(rk4) , dimension(ni,nj,nk) :: q , t
  real(rk4) , dimension(nk) :: sigma
  intent (in) ni , nj , nk , ps , ptop , sigma , t
  intent (inout) q
!
  real(rk4) :: lh , p , qs , pt , satvp
  integer(ik4) :: i , j , k
!
!     THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
!
  pt = real(ptop)
  do i = 1 , ni
    do j = 1 , nj
      do k = 1 , nk
        p = (pt+sigma(k)*ps(i,j))*10.
        lh = slh0 - slh1*(t(i,j,k)-t0)
        satvp = slsvp1*exp(slsvp2*lh*(rt0-1./t(i,j,k)))
        qs = sep2*satvp/(p-satvp)
        q(i,j,k) = amax1(q(i,j,k)*qs,0.0)
      end do
    end do
  end do
!
  end subroutine humid2
!
end module mod_humid
