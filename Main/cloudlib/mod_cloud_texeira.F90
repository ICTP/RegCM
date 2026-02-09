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

module mod_cloud_texeira

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none (type, external)

  private

  public :: texeira_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! See Cloud Fraction and Relative Humidity in a Prognostic Cloud
  ! Fraction Scheme, João Teixeira,
  ! https://doi.org/10.1175/1520-0493(2001)129<1750:CFARHI>2.0.CO;2
  subroutine texeira_cldfrac(qc,qs,rh,rh0,qcrit,fcc)
    implicit none (type, external)
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: qs, qc, rh
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: rh0, qcrit
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: fcc
    integer(ik4) :: i, j, k

    real(rkx), parameter :: kappa = 1.0e-6 ! sec-1
    real(rkx), parameter :: d = 4.0e-6 ! sec-1
    real(rkx) :: rhrng, liq, spq

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      ! Adjusted relative humidity threshold
      rhrng = min(max(rh(j,i,k),0.001_rkx),0.999_rkx)
      liq = qc(j,i,k)
      spq = qs(j,i,k) / (d_one + qs(j,i,k))
      if ( liq > qcrit(j,i) .and. rhrng > rh0(j,i) ) then
        fcc(j,i,k) = d*liq / (d_two*spq*(d_one-rhrng)*kappa) * &
          ( -d_one + sqrt(d_one + &
          (d_four*spq*((d_one-rhrng)*kappa))/(d*liq)) )
      else
        fcc(j,i,k) = d_zero
      end if
    end do

  end subroutine texeira_cldfrac

end module mod_cloud_texeira
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
