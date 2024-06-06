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

  use mod_intkinds
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
  subroutine gulisa_cldfrac(qt,z,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qt
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    integer(ik4) :: i , j , k
    real(rkx) :: qgkg , stddev

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

#ifdef STDPAR
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz ) &
      local(qgkg,stddev)
#else
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          stddev = 1.0_rkx + (0.04_rkx/(z(j,i,k)/80000.0_rkx*0.625_rkx)) * &
            exp(-log(0.0005_rkx*z(j,i,k))**2/0.625_rkx)
          qgkg = qt(j,i,k)*1.0e3_rkx/stddev
          if ( qgkg < 0.18_rkx ) then
            fcc(j,i,k) = d_zero
          else if ( qgkg > 2.0_rkx ) then
            fcc(j,i,k) = d_one
          else
            fcc(j,i,k) = -0.1754_rkx + 0.9811_rkx*qgkg - &
                                       0.2223_rkx*qgkg*qgkg - &
                                       0.0104_rkx*qgkg*qgkg*qgkg
          end if
#ifndef STDPAR
        end do
      end do
#endif
    end do

  end subroutine gulisa_cldfrac

end module mod_cloud_guli2007
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
