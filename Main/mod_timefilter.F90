!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_timefilter
  !
  ! Implementation of time filtering technique
  !
  ! RA  Robert-Asselin filter Eq 2.4.6 in the MM5 manual
  !
  ! RAW Robert-Asselin-Williams filter - MWR 2009
  !              http://dx.doi.org/10.1175/2009MWR2724.1
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams, only : iqv

  implicit none

  private

  real(rkx), public, parameter :: betaraw = 0.53_rkx

  public :: timefilter_apply

  interface timefilter_apply
    module procedure filter_raw_2d
    module procedure filter_raw_3d
    module procedure filter_raw_4d
    module procedure filter_raw_qv
    module procedure filter_raw_uv
    module procedure filter_ra_2d
    module procedure filter_ra_3d
    module procedure filter_ra_4d
    module procedure filter_ra_qv
    module procedure filter_ra_uv
  end interface timefilter_apply

  contains

  subroutine filter_ra_2d(phin,phinm1,phinp1,alpha)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha
    integer(ik4) :: i, j
    real(rkx) :: d
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      d = alpha * (phinp1(j,i) + phinm1(j,i) - d_two * phin(j,i))
      phinm1(j,i) = phin(j,i) + d
      phin(j,i) = phinp1(j,i)
    end do
  end subroutine filter_ra_2d

  subroutine filter_raw_2d(phin,phinm1,phinp1,alpha,beta)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha, beta
    integer(ik4) :: i, j
    real(rkx) :: d
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      d = alpha * (phinp1(j,i) + phinm1(j,i) - d_two * phin(j,i))
      phinm1(j,i) = phin(j,i) + beta * d
      phin(j,i) = phinp1(j,i) + (beta - d_one) * d
    end do
  end subroutine filter_raw_2d

  subroutine filter_ra_3d(phin,phinm1,phinp1,alpha)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
      d = alpha * (phinp1(j,i,k) + phinm1(j,i,k) - d_two * phin(j,i,k))
      phinm1(j,i,k) = phin(j,i,k) + d
      phin(j,i,k) = phinp1(j,i,k)
    end do
  end subroutine filter_ra_3d

  subroutine filter_raw_3d(phin,phinm1,phinp1,alpha,beta)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha, beta
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
      d = alpha * (phinp1(j,i,k) + phinm1(j,i,k) - d_two * phin(j,i,k))
      phinm1(j,i,k) = phin(j,i,k) + beta * d
      phin(j,i,k) = phinp1(j,i,k) + (beta - d_one) * d
    end do
  end subroutine filter_raw_3d

  subroutine filter_ra_uv(un,unm1,unp1,vn,vnm1,vnp1,alpha)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: unm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: un
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: unp1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vnm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vn
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: vnp1
    real(rkx), intent(in) :: alpha
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(un,3)
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:nk )
      d = alpha * (unp1(j,i,k) + unm1(j,i,k) - d_two * un(j,i,k))
      unm1(j,i,k) = un(j,i,k) + d
      un(j,i,k) = unp1(j,i,k)
      d = alpha * (vnp1(j,i,k) + vnm1(j,i,k) - d_two * vn(j,i,k))
      vnm1(j,i,k) = vn(j,i,k) + d
      vn(j,i,k) = vnp1(j,i,k)
    end do
  end subroutine filter_ra_uv

  subroutine filter_raw_uv(un,unm1,unp1,vn,vnm1,vnp1,alpha,beta)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: unm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: un
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: unp1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vnm1
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vn
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: vnp1
    real(rkx), intent(in) :: alpha, beta
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(un,3)
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:nk )
      d = alpha * (unp1(j,i,k) + unm1(j,i,k) - d_two * un(j,i,k))
      unm1(j,i,k) = un(j,i,k) + beta * d
      un(j,i,k) = unp1(j,i,k) + (beta - d_one) * d
      d = alpha * (vnp1(j,i,k) + vnm1(j,i,k) - d_two * vn(j,i,k))
      vnm1(j,i,k) = vn(j,i,k) + beta * d
      vn(j,i,k) = vnp1(j,i,k) + (beta - d_one) * d
    end do
  end subroutine filter_raw_uv

  subroutine filter_ra_4d(phin,phinm1,phinp1,alpha,n1,n2,low)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha, low
    integer(ik4), intent(in) :: n1, n2
    integer(ik4) :: i, j, k, n, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk, n = n1:n2 )
      d = alpha * (phinp1(j,i,k,n) + phinm1(j,i,k,n) - &
          d_two * phin(j,i,k,n))
      phinm1(j,i,k,n) = phin(j,i,k,n) + d
      if ( phinm1(j,i,k,n) < low ) phinm1(j,i,k,n) = d_zero
      phin(j,i,k,n) = phinp1(j,i,k,n)
    end do
  end subroutine filter_ra_4d

  subroutine filter_ra_qv(phin,phinm1,phinp1,alpha,ps)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(in) :: phinp1
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: ps
    real(rkx), intent(in) :: alpha
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
      d = alpha * (phinp1(j,i,k,iqv) + phinm1(j,i,k,iqv) - &
          d_two * phin(j,i,k,iqv))
      phinm1(j,i,k,iqv) = max(phin(j,i,k,iqv) + d, minqq*ps(j,i))
      phin(j,i,k,iqv) = phinp1(j,i,k,iqv)
    end do
  end subroutine filter_ra_qv

  subroutine filter_raw_4d(phin,phinm1,phinp1,alpha,beta,n1,n2,low)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(in) :: phinp1
    real(rkx), intent(in) :: alpha, beta, low
    integer(ik4), intent(in) :: n1, n2
    integer(ik4) :: i, j, k, n, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk, n = n1:n2 )
      d = alpha * (phinp1(j,i,k,n) + phinm1(j,i,k,n) - &
          d_two * phin(j,i,k,n))
      phinm1(j,i,k,n) = phin(j,i,k,n) + beta * d
      phin(j,i,k,n) = phinp1(j,i,k,n) + (beta - d_one) * d
      if ( phinm1(j,i,k,n) < low ) phinm1(j,i,k,n) = d_zero
      if ( phin(j,i,k,n) < low ) phin(j,i,k,n) = d_zero
    end do
  end subroutine filter_raw_4d

  subroutine filter_raw_qv(phin,phinm1,phinp1,alpha,beta,psa,psb)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phinm1
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: phin
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(in) :: phinp1
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: psa, psb
    real(rkx), intent(in) :: alpha, beta
    integer(ik4) :: i, j, k, nk
    real(rkx) :: d
    nk = size(phin,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
      d = alpha * (phinp1(j,i,k,iqv) + phinm1(j,i,k,iqv) - &
          d_two * phin(j,i,k,iqv))
      phinm1(j,i,k,iqv) = max(phin(j,i,k,iqv) + beta * d, minqq*psa(j,i))
      phin(j,i,k,iqv) = max(phinp1(j,i,k,iqv) + &
                      (beta - d_one) * d, minqq*psb(j,i))
    end do
  end subroutine filter_raw_qv

end module mod_timefilter
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
