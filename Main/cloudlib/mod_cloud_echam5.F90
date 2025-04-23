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

module mod_cloud_echam5

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: echam5_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! See Johannes Quaas
  ! Evaluating the “critical relative humidity” as a measure
  ! of subgrid-scale variability of humidity in general circulation
  ! model cloud cover parameterizations using satellite data
  !
  subroutine echam5_cldfrac(qc,rh,p,ps,qcrit,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qc , rh , p
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps , qcrit
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    real(rkx) , parameter :: ct = 0.35_rkx
    real(rkx) , parameter :: cs = 0.85_rkx
    real(rkx) , parameter :: nx = 4.0_rkx
    integer(ik4) :: i , j , k
    real(rkx) :: rhrng , rhcrit , sig

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      if ( qc(j,i,k) > qcrit(j,i) ) then
        ! Relative humidity
        rhrng = min(max(rh(j,i,k),0.001_rkx),0.999_rkx)
        sig = ps(j,i)/p(j,i,k)
        rhcrit = ct + (ct-cs)*exp(1.0_rkx-sig**nx)
        if ( rhrng < rhcrit ) then
          fcc(j,i,k) = d_zero
        else if ( rhrng > 0.99999_rkx ) then
          fcc(j,i,k) = d_one
        else
          ! Sundqvist formula
          fcc(j,i,k) = 1.0_rkx - sqrt((1.0_rkx-rhrng) / &
                                      (1.0_rkx-rhcrit))
        end if
      else
        fcc(j,i,k) = d_zero
      end if
    end do

  end subroutine echam5_cldfrac

end module mod_cloud_echam5
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
