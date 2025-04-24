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

module mod_cloud_xuran

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: xuran_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud fraction
  ! using the semi-empirical formula of Xu and Randall (1996, JAS)
  !
  subroutine xuran_cldfrac(p,qc,qs,rh,qcrit,fcc)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: p, rh
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: qc, qs
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: qcrit
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: fcc
    integer(ik4) :: i, j, k
    real(rkx), parameter :: parm_p = 0.25_rkx
    real(rkx), parameter :: parm_gamma = 0.49_rkx
    real(rkx), parameter :: parm_alpha0 = 100.0_rkx
    real(rkx) :: botm, rm, qcld, rhrng

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      if ( qc(j,i,k) > qcrit(j,i) ) then
        qcld = qc(j,i,k)
        rhrng = max(0.0_rkx,min(1.0_rkx,rh(j,i,k)))
        if ( rhrng > 0.99999 ) then
          fcc(j,i,k) = d_one
        else
          botm = rhrng ** parm_p
          rm = -min((parm_alpha0 * qcld) / &
            ((d_one-rhrng)*qs(j,i,k))**parm_gamma,25.0_rkx)
          fcc(j,i,k) = botm * (1.0_rkx - exp(rm))
        end if
      else
        fcc(j,i,k) = d_zero
      end if
    end do

  end subroutine xuran_cldfrac

end module mod_cloud_xuran
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
