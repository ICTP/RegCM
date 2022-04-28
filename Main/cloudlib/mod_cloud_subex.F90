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

module mod_cloud_subex

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: subex_cldfrac

  contains
  !
  ! This subroutine computes the fractional cloud coverage
  ! which is used in radiation.
  !
  ! The fractional coverage of large scale clouds is a function of
  ! relative humidity, using the relationship of sundqvist et
  ! al., 1989.  The relative humidity at which clouds begin to
  ! form is lower over land than ocean, due to the greater number
  ! of cloud condensation nucleii.
  !
  ! See Pal et al (2000) for more info.
  !
  subroutine subex_cldfrac(t,p,qv,qc,rh,tc0,rh0,fcc)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: t , p , qv , qc , rh
    real(rkx) , pointer , dimension(:,:) , intent(in) :: rh0
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: fcc
    real(rkx) , intent(in) :: tc0
    real(rkx) :: rh0adj , rhrng
    integer(ik4) :: i , j , k

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( qc(j,i,k) > dlowval ) then
            ! Adjusted relative humidity threshold
            rhrng = min(max(rh(j,i,k),rhmin),1.0_rkx)
            if ( t(j,i,k) > tc0 ) then
              rh0adj = rh0(j,i)
            else ! high cloud (less subgrid variability)
              rh0adj = 0.99999_rkx - &
                  (d_one-rh0(j,i))/(d_one+0.15_rkx*(tc0-t(j,i,k)))
            end if
            if ( rhrng <= rh0adj ) then
              fcc(j,i,k) = d_zero
            else if ( rhrng > 0.99999_rkx ) then
              fcc(j,i,k) = d_one
            else
              ! Use Sundqvist (1989) formula
              fcc(j,i,k) = d_one-sqrt((d_one-rhrng)/(d_one-rh0adj))
            end if
          else
            fcc(j,i,k) = d_zero
          end if
        end do
      end do
    end do
    !
    ! Correction:
    !   Ivan Guettler, 14.10.2010.
    ! Based on: Vavrus, S. and Waliser D., 2008,
    ! An Improved Parameterization for Simulating Arctic Cloud Amount
    ! in the CCSM3 Climate Model, J. Climate
    !
    if ( larcticcorr ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! clouds below 750hPa, extremely cold conditions,
            !  when no cld microphy
            if ( p(j,i,k) >= 75000.0_rkx .and. qv(j,i,k) <= 0.003_rkx ) then
              fcc(j,i,k) = fcc(j,i,k) * &
                    max(0.15_rkx,min(d_one,qv(j,i,k)/0.003_rkx))
            end if
          end do
        end do
      end do
    end if

  end subroutine subex_cldfrac

end module mod_cloud_subex
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
