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

module mod_cloud_texeira

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

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
  subroutine texeira_cldfrac(qc,qs,rh,rh0,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qs , qc , rh
    real(rkx) , pointer , dimension(:,:) , intent(in) :: rh0
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    real(rkx) :: rhrng
    integer(ik4) :: i , j , k

    real(rkx) , parameter :: kappa = 1.0e-6 ! sec-1
    real(rkx) , parameter :: d = 4.0e-6 ! sec-1
    real(rkx) :: liq , spq

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Adjusted relative humidity threshold
          rhrng = min(max(rh(j,i,k),0.001_rkx),0.999_rkx)
          liq = qc(j,i,k)
          spq = qs(j,i,k) / (d_one + qs(j,i,k))
          if ( liq > 1.0e-12_rkx .and. rhrng > rh0(j,i) ) then
            fcc(j,i,k) = d*liq / (d_two*spq*(d_one-rhrng)*kappa) * &
              ( -d_one + sqrt(d_one + &
                 (d_four*spq*((d_one-rhrng)*kappa))/(d*liq)) )
          else
            fcc(j,i,k) = d_zero
          end if
        end do
      end do
    end do

  end subroutine texeira_cldfrac

end module mod_cloud_texeira
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
