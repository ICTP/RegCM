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

module mod_cloud_tompkins

  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: tompkins_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! See A cloud scheme for data assimilation purposes
  ! ECMWF Technical Memorandum 410
  subroutine tompkins_cldfrac(qc,qs,rh,p,ps,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qs , qc , rh , p
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    real(rkx) :: rhrng , kappa , rhcrit
    real(rkx) :: spq , spqs , sig
    integer(ik4) :: i , j , k

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Adjusted relative humidity threshold
          rhrng = min(max(rh(j,i,k),0.001_rkx),0.999_rkx)
          spq = qc(j,i,k) / (d_one + qc(j,i,k))
          spqs = qs(j,i,k) / (d_one + qs(j,i,k))
          if ( spq > 1.0e-8_rkx ) then
            sig = p(j,i,k)/ps(j,i)
            kappa = max(0.0_rkx,0.9_rkx*(sig-0.2_rkx)**0.2_rkx)
            rhcrit = 0.7_rkx * sig * (1.0_rkx-sig) * &
                   (1.85_rkx + 0.95_rkx*(sig-0.5_rkx))
            fcc(j,i,k) = 1.0_rkx - sqrt((1.0_rkx-rhrng) / &
                   (1.0_rkx - rhcrit - kappa*(rhrng-rhcrit)))
          else
            fcc(j,i,k) = d_zero
          end if
        end do
      end do
    end do

  end subroutine tompkins_cldfrac

end module mod_cloud_tompkins
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
