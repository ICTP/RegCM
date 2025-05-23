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

module mod_che_ccn

  use mod_intkinds
  use mod_constants
  use mod_realkinds
  use mod_che_common
  use mod_che_carbonaer

  private

  ! Parameters for calculating ccn concentrations form aerosol mass
  ! densities in kg/m3
  ! other aerosol densities are passed by module already
  real(rkx), parameter :: rhos = 1800.0_rkx
  real(rkx), parameter :: rhosslt = 1000.0_rkx
  ! coef_ccn, abulk are parameter adjustable from the namelist
  ! now particule to ccn convesrion following Lee et al.,
  ! 2016 update from Jones et al., 1998
  real(rkx), parameter :: c1 = 4.34e8_rkx
  real(rkx), parameter :: c2 = -1.8e-9_rkx

  ! minimal number concentration of ccn = 30 cm-3
  real(rkx), parameter :: ccnmin = 30.e6_rkx

  public :: ccn, calc_ccn

  contains

  subroutine ccn
    implicit none
    integer(ik4) :: i, j, k, n
    cccn(:,:,:) = d_zero
    if ( nbchl > 0 ) then
      do n = 1, nbchl
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              cccn(j,i,k) = cccn(j,i,k) + &
                 calc_ccn(rhobchl(n),chib3d(j,i,k,ibchl(n))*crhob3d(j,i,k))
            end do
          end do
        end do
      end do
    end if
    if ( nochl > 0 ) then
      do n = 1, nochl
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              cccn(j,i,k) = cccn(j,i,k) + &
                 calc_ccn(rhoochl(n),chib3d(j,i,k,iochl(n))*crhob3d(j,i,k))
            end do
          end do
        end do
      end do
    end if
    if ( iso4 > 0 ) then
      do k = 1, kz
        do i = ici1, ici2
          do j = jci1, jci2
            cccn(j,i,k) = cccn(j,i,k) + &
                 calc_ccn(rhos,chib3d(j,i,k,iso4)*crhob3d(j,i,k))
          end do
        end do
      end do
    end if
    if ( isslt(1) > 0 ) then
      do k = 1, kz
        do i = ici1, ici2
          do j = jci1, jci2
            cccn(j,i,k) = cccn(j,i,k) + &
                 calc_ccn(rhosslt,chib3d(j,i,k,isslt(1))*crhob3d(j,i,k))
          end do
        end do
      end do
    end if
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          ! now calculate ccn number from particle number
          ! cccn = cccn * abulk
          cccn(j,i,k) = c1 * (d_one - exp(c2 * cccn(j,i,k)))
          !
          ! finally consider a minimal concentration of CCN
          cccn(j,i,k) = max(ccnmin, cccn(j,i,k))
        end do
      end do
    end do
  end subroutine ccn

  pure elemental real(rkx) function calc_ccn(denx,mrat) result(res)
    implicit none
    ! Calculate the particle number from the mass distriubution of
    ! hydrophilic particle
    ! coef_ccn is a coefficient detremined by assuming a lognormal mass
    ! distributionand calculated as (1/Dm**3)/exp(-9(logsigma)**2/2)
    ! passed in nameliste
    real(rkx), intent(in) :: denx, mrat
    !
    res = 6.0_rkx / (mathpi * denx ) * mrat * coef_ccn
  end function calc_ccn

end module mod_che_ccn

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
