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

module mod_cloud_guli2007

  use mod_realkinds
  use mod_constants
  use mod_dynparam

  implicit none

  private

  public :: gulisa_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! Details on the formulation is in:
  !  Gultepe and Isaac,
  !   Cloud fraction parametrization as a function of mean cloud water
  !   content and its variance using in-situ observations
  !   GRL vol 34, L07801, 2007
  !
  subroutine gulisa_cldfrac(qv,qs,qt,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qt
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qv
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qs
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    integer(ik4) :: i , j , k
    real(rkx) :: qgkg

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          qgkg = (qt(j,i,k)+max(qv(j,i,k)-qs(j,i,k),d_zero))*1.0e3_rkx
          fcc(j,i,k) = 1.0_rkx/90.0_rkx * &
                      ((100_rkx-ds)  * 5.57_rkx * qgkg**0.78_rkx + &
                       (ds-10.0_rkx) * 4.82_rkx * qgkg**0.94_rkx)
        end do
      end do
    end do

  end subroutine gulisa_cldfrac

end module mod_cloud_guli2007
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
