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
  subroutine gulisa_cldfrac(qt,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qt
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: fcc
    integer(ik4) :: i , j , k
    real(rkx) :: qgkg , s10 , s100 , r10 , r100
    real(rkx) :: cs10 , cs100 , w1 , w2

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          qgkg = qt(j,i,k)*d_100
          s10 = 0.28_rkx + qgkg ** 0.49_rkx
          r10 = qgkg/s10
          s100 = 0.18_rkx + qgkg ** 0.48_rkx
          r100 = qgkg/s100
          if ( r10 < 0.18_rkx ) then
            cs10 = d_zero
          else if ( r10 > 2.0_rkx ) then
            cs10 = d_one
          else
            cs10 = -0.1754_rkx + 0.9811_rkx * r10 -      &
                                 0.2223_rkx * r10*r10 +  &
                                 0.0104_rkx * r10*r10*r10
          end if
          if ( r100 < 0.12_rkx ) then
            cs100 = d_zero
          else if ( r10 > 1.85_rkx ) then
            cs100 = d_one
          else
            cs100 = -0.0913_rkx + 0.7213_rkx * r10 +     &
                                  0.1060_rkx * r10*r10 - &
                                  0.0946_rkx * r10*r10*r10
          end if
          if ( ds <= 10.0_rkx ) then
            fcc(j,i,k) = cs10
          else if ( ds >= 100.0_rkx ) then
            fcc(j,i,k) = cs100
          else
            w2 = (ds-10.0_rkx)/90.0_rkx
            w1 = d_one - w2
            fcc(j,i,k) = cs10 * w1 + cs100 * w2
          end if
        end do
      end do
    end do

  end subroutine gulisa_cldfrac

end module mod_cloud_guli2007
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
