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

  interface mxr2rh
    module procedure mxr2rh
    module procedure mxr2rhfv
    module procedure mxr2rh_o_double
    module procedure mxr2rh_o_single
    module procedure mxr2rh_o_double_nonhydro
    module procedure mxr2rh_o_single_nonhydro
  end interface mxr2rh

  interface sph2mxr
    module procedure sph2mxr_double
    module procedure sph2mxr_single
  end interface sph2mxr

  public :: mxr2rh , rh2mxr , clwfromt
  public :: sph2mxr , mxr2sph

  contains

#include <pfesat.inc>
#include <pfwsat.inc>
#include <sig2p.inc>

  subroutine sph2mxr_double(q,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q
    integer(ik4) :: i , j , k
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          q(i,j,k) = q(i,j,k) / (1.0D0 - q(i,j,k))
        end do
      end do
    end do
  end subroutine sph2mxr_double

  subroutine sph2mxr_single(q,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk4) , intent(inout) , dimension(ni,nj,nk) :: q
    integer(ik4) :: i , j , k
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          q(i,j,k) = q(i,j,k) / (1.0 - q(i,j,k))
        end do
      end do
    end do
  end subroutine sph2mxr_single

  subroutine mxr2sph(q,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rkx) , intent(inout) , dimension(ni,nj,nk) :: q
    integer(ik4) :: i , j , k
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          q(i,j,k) = q(i,j,k) / (d_one + q(i,j,k))
        end do
      end do
    end do
  end subroutine mxr2sph

  subroutine mxr2rh(t,q,ps,ptop,sigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rkx) , intent(in) :: ps , ptop
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: t
    real(rkx) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rkx) , intent(in) , dimension(nk) :: sigma
    real(rkx) :: p , qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = sig2p(ps,sigma(k),ptop)
          qs = pfwsat(t(i,j,k),p)
          q(i,j,k) = max(q(i,j,k)/qs,d_zero)
        end do
      end do
    end do
  end subroutine mxr2rh

  subroutine mxr2rh_o_double(t,q,ps,sigma,ptop,im,jm,km)
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
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    ! DATA ON SIGMA LEVELS
    !
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! PS in output file is ps + ptop
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          p = (sigma(k)*(ps(i,j)-ptop*100.0D0)) + ptop * 100.0D0
          qs = pfwsat(real(t(i,j,k),rkx),real(p,rkx))
          q(i,j,k) = max(q(i,j,k)/qs,0.0D0)
        end do
      end do
    end do
  end subroutine mxr2rh_o_double

  subroutine mxr2rh_o_single(t,q,ps,sigma,ptop,im,jm,km)
    implicit none
    integer(ik4) , intent(in) :: im , jm , km
    real(rkx) , intent(in) :: ptop
    real(rk4) , intent(in) , dimension(im,jm) :: ps
    real(rk4) , intent(in) , dimension(im,jm,km) :: t
    real(rk4) , intent(in) , dimension(km) :: sigma
    real(rk4) , intent(inout) , dimension(im,jm,km) :: q

    real(rk4) :: qs
    real(rkx) :: p
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    ! DATA ON SIGMA LEVELS
    !
    do k = 1 , km
      do j = 1 , jm
        do i = 1 , im
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! PS in output file is ps + ptop
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          p = real(sigma(k),rkx)*(real(ps(i,j),rkx) - ptop*d_100) + ptop * d_100
          qs = real(pfwsat(real(t(i,j,k),rkx),p))
          q(i,j,k) = max(q(i,j,k)/qs,0.0)
        end do
      end do
    end do
  end subroutine mxr2rh_o_single

  subroutine mxr2rhfv(t,q,p3d,ni,nj,nk,mval)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: p3d , t
    real(rkx) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rkx) , intent(in) :: mval

    real(rkx) :: qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          if ( p3d(i,j,k) > mval ) then
            qs = pfwsat(t(i,j,k),p3d(i,j,k)*d_100) ! P in mb -> Pa
            q(i,j,k) = max(q(i,j,k)/qs,d_zero)
          else
            q(i,j,k) = mval
          end if
        end do
      end do
    end do
  end subroutine mxr2rhfv

  subroutine mxr2rh_o_double_nonhydro(t,q,p3d,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk8) , intent(in) , dimension(ni,nj,nk) :: p3d , t
    real(rk8) , intent(inout) , dimension(ni,nj,nk) :: q

    real(rk8) :: qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          qs = pfwsat(real(t(i,j,k),rkx),real(p3d(i,j,k),rkx))
          q(i,j,k) = max(q(i,j,k)/qs,0.0D0)
        end do
      end do
    end do
  end subroutine mxr2rh_o_double_nonhydro

  subroutine mxr2rh_o_single_nonhydro(t,q,p3d,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rk4) , intent(in) , dimension(ni,nj,nk) :: p3d , t
    real(rk4) , intent(inout) , dimension(ni,nj,nk) :: q

    real(rk4) :: qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES MIXING RATIO BY RELATIVE HUMIDITY
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          qs = real(pfwsat(real(t(i,j,k),rkx),real(p3d(i,j,k),rkx)))
          q(i,j,k) = max(q(i,j,k)/qs,0.0)
        end do
      end do
    end do
  end subroutine mxr2rh_o_single_nonhydro

  subroutine rh2mxr(t,q,ps,ptop,sigma,ni,nj,nk)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    real(rkx) , intent(in) :: ptop
    real(rkx) , intent(in) , dimension(ni,nj) :: ps
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: t
    real(rkx) , intent(inout) , dimension(ni,nj,nk) :: q
    real(rkx) , intent(in) , dimension(nk) :: sigma

    real(rkx) :: p , qs
    integer(ik4) :: i , j , k
    !
    ! THIS ROUTINE REPLACES RELATIVE HUMIDITY BY MIXING RATIO
    !
    do k = 1 , nk
      do j = 1 , nj
        do i = 1 , ni
          p = (ptop + sigma(k)*ps(i,j))*d_1000
          qs = pfwsat(t(i,j,k),p)
          q(i,j,k) = max(q(i,j,k)*qs,d_zero)
        end do
      end do
    end do
  end subroutine rh2mxr

  real(rkx) function clwfromt(t) result(clw)
    implicit none
    real(rkx) , intent(in) :: t
    real(rkx) :: tcel
    ! Temperature dependency for cloud water content
    ! in g/m3 (Lemus et al., 1997)
    ! NOTE : THIS IS IN-CLOUD VARIABLE.
    tcel = t - tzero
    if ( tcel < -50.0_rkx ) then
      clw = 0.001_rkx
    else
      clw = 0.127_rkx + 6.78e-03_rkx * tcel +    &
                        1.29e-04_rkx * tcel**2 + &
                        8.36e-07_rkx * tcel**3
    end if
  end function clwfromt

end module mod_humid

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
