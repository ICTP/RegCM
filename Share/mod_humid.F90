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

  implicit none

  private

  interface humid1_o
    module procedure humid1_o_double
    module procedure humid1_o_single
  end interface humid1_o

  public :: humid1 , humid1_o , humid1fv , humid2

  contains

  subroutine humid1(t,q,ps,ptop,sigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) :: ps , ptop
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rk8) , intent(in) , dimension(nk) :: sigma

    real(rk8) :: p , qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = sig2p(ps,sigma(k),ptop)
          qs = pfqsat(t(i,j,k),p)
          q(i,j,k) = max(q(i,j,k)/qs,d_zero)
        end do
      end do
    end do
  end subroutine humid1

  subroutine humid1_o_double(t,q,ps,sigma,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk8) , intent(in) , dimension(im,jm) :: ps
    real(rk8) , intent(in) , dimension(im,jm,km) :: t
    real(rk8) , intent(in) , dimension(km) :: sigma
    real(rk8) , intent(inout) , dimension(im,jm,km) :: q

    real(rk8) :: p , qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    ! DATA ON SIGMA LEVELS
    !
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! PS in output file is ps + ptop and bot are hPa !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          p = (sigma(k)*ps(i,j) + ptop) * d_100
          qs = pfqsat(t(i,j,k),p)
          q(i,j,k) = max(q(i,j,k)/qs,d_zero)
        end do
      end do
    end do
  end subroutine humid1_o_double

  subroutine humid1_o_single(t,q,ps,sigma,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rk8) , intent(in) :: ptop
    real(rk4) , intent(in) , dimension(im,jm) :: ps
    real(rk4) , intent(in) , dimension(im,jm,km) :: t
    real(rk4) , intent(in) , dimension(km) :: sigma
    real(rk4) , intent(inout) , dimension(im,jm,km) :: q

    real(rk4) :: qs
    real(rk8) :: p
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    ! DATA ON SIGMA LEVELS
    !
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! PS in output file is ps + ptop and bot are hPa !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          p = (dble(sigma(k))*dble(ps(i,j)) + ptop) * d_100
          qs = real(pfqsat(dble(t(i,j,k)),p))
          q(i,j,k) = max(q(i,j,k)/qs,0.0)
        end do
      end do
    end do
  end subroutine humid1_o_single
!
!-----------------------------------------------------------------------
!
  subroutine humid1fv(t,q,p3d,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: p3d , t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q

    real(rk8) :: qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          if ( p3d(i,j,k) > -9990.0D0 ) then
            qs = pfqsat(t(i,j,k),p3d(i,j,k)*d_100) ! P in mb -> Pa
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

    real(rk8) :: p , qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = (ptop + sigma(k)*ps(i,j))*d_1000
          qs = pfqsat(t(i,j,k),p)
          q(i,j,k) = max(q(i,j,k)*qs,d_zero)
        end do
      end do
    end do
  end subroutine humid2
!
end module mod_humid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
