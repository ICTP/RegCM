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

module mod_cloud_subex

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: subex_cldfrac

  logical, parameter :: use_dynamic_rhcrit = .true.

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
  subroutine subex_cldfrac(t,p,qv,rh,tc0,rh0,fcc)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: t, p
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: qv, rh
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: rh0
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: fcc
    real(rkx), intent(in) :: tc0
    integer(ik4) :: i, j, k
    real(rkx) :: rh0adj, rhrng

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      ! Use Pal et al. formula
      ! rhrng = rh(j,i,k)
      ! Adjusted relative humidity threshold
      rhrng = min(max(rh(j,i,k),rhmin),rhmax)
      if ( use_dynamic_rhcrit ) then
        rh0adj = dynamic_rhcrit(rh0(j,i),t(j,i,k))
      else
        if ( t(j,i,k) > tc0 ) then
          rh0adj = rh0(j,i)
        else ! high cloud (less subgrid variability)
          ! Use Pal et al. formula
          rh0adj = rhmax - (rhmax-rh0(j,i))/(d_one+0.15_rkx*(tc0-t(j,i,k)))
          ! Adjusted for Sundqvist
          !rh0adj = d_one - &
          !  (d_one-rh0(j,i))/(d_one+0.15_rkx*(tc0-t(j,i,k)))
        end if
      end if
      if ( rhrng <= rh0adj ) then
        fcc(j,i,k) = d_zero
      else if ( rhrng > 0.99999_rkx ) then
        fcc(j,i,k) = d_one
      else
        ! Use Pal et al. (2000) formula
        ! fcc(j,i,k) = sqrt((rhrng-rh0adj)/(rhmax-rh0adj))
        ! Use Sundqvist (1989) formula
        fcc(j,i,k) = rhmax - sqrt((rhmax-rhrng)/(rhmax-rh0adj))
      end if
    end do
    !
    ! Correction:
    !   Ivan Guettler, 14.10.2010.
    ! Based on: Vavrus, S. and Waliser D., 2008,
    ! An Improved Parameterization for Simulating Arctic Cloud Amount
    ! in the CCSM3 Climate Model, J. Climate
    !
    if ( larcticcorr ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ! clouds below 750hPa, extremely cold conditions,
        !  when no cld microphy
        if ( p(j,i,k) >= 75000.0_rkx .and. qv(j,i,k) <= 0.003_rkx ) then
          fcc(j,i,k) = fcc(j,i,k) * &
                    max(0.15_rkx,min(d_one,qv(j,i,k)/0.003_rkx))
        end if
      end do
    end if
  end subroutine subex_cldfrac
  !
  ! A dynamic critical relative humidity based on temperature in
  ! cloud parameterization to improve low cloud in an AGCM
  ! https://doi.org/10.1016/j.atmosres.2025.108183
  !
  pure real(rkx) function dynamic_rhcrit(rh0,tk) result(rhc)
    !$acc routine seq
    implicit none
    real(rkx), intent(in) :: rh0, tk
    real(rkx) :: tku
    tku = min(301.0_rkx,max(216.0_rkx,tk))
    rhc = 0.01_rkx * min(60.0_rkx, max(100.0_rkx, rh0*100.0_rkx + &
      -21760.24_rkx-0.000004969437_rkx*tku**4+0.005182194_rkx*tku**3 - &
      2.009108_rkx*tku**2+343.0437_rkx*tku))
  end function dynamic_rhcrit

end module mod_cloud_subex
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
