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

module mod_che_carbonaer

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  implicit none

  private

  ! Parameter usefull for wet and dry deposition of carbon aerosol
  ! densities in kg/m3

  real(rkx) , public , parameter :: rhobc   = 2000.0_rkx
  real(rkx) , public , parameter :: rhooc   = 1200.0_rkx
  real(rkx) , public , parameter :: rhobchl = 1600.0_rkx
  real(rkx) , public , parameter :: rhoochl = 1200.0_rkx
  real(rkx) , public , parameter :: rhoso4 = 1830.0_rkx  !(95-96% H2SO4)
  real(rkx) , public , parameter :: rhosm1 = 1200.0_rkx
  real(rkx) , public , parameter :: rhosm2 = 1200.0_rkx

  ! effctive diameters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? )
  real(rkx) , public , parameter :: reffbc   = 0.05_rkx
  real(rkx) , public , parameter :: reffbchl = 0.3_rkx
  real(rkx) , public , parameter :: reffoc   = 0.2_rkx
  real(rkx) , public , parameter :: reffochl = 0.3_rkx
  ! (75% H2SO4,25% H2O; Kiehl and Briegleb, 1993)
  real(rkx) , public , parameter :: reffso4 = 0.1_rkx
  real(rkx) , public , parameter :: reffsm1 = 0.3_rkx
  real(rkx) , public , parameter :: reffsm2 = 0.3_rkx

  ! aging efolding time (s), from hydrophobic to hydrophilic
  ! Cooke et al.
  real(rkx) :: chagct
  real(rkx) , parameter :: const_chagct = 1.15_rkx * 86400.0_rkx
  real(rkx) , parameter :: chsmct = 1.15_rkx * 86400.0_rkx
  real(rkx) , parameter :: kcond = 0.01_rkx ! nm-1
  real(rkx) , parameter :: kcoag = 6e-6_rkx ! cm^3/h

  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rkx) , parameter :: solbc = 0.05_rkx
  real(rkx) , parameter :: solbchl = 0.8_rkx
  real(rkx) , parameter :: soloc = 0.05_rkx
  real(rkx) , parameter :: solochl = 0.8_rkx
  real(rkx) , parameter :: solsm1 = 0.05_rkx
  real(rkx) , parameter :: solsm2 = 0.8_rkx

  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(rkx) , public , dimension(9) :: carbed

  real(rkx) , pointer , dimension(:,:,:) :: ncon
  ! variable for surface area of aerosol
  real(rkx) , pointer , dimension(:,:,:) :: surf
  real(rkx) , pointer , dimension(:,:,:) :: so4chagct

  public :: solbc , solbchl , soloc , solochl , solsm1 , solsm2
  ! unit of so4chagct = kg kg-1 sec-1
  public :: so4chagct
  public :: carb_init , carb_prepare , aging_carb

  contains

    subroutine carb_init( )
      implicit none
      if ( carb_aging_control ) then
        call getmem3d(ncon,jci1,jci2,ici1,ici2,1,kz,'carbonaer:ncon')
        call getmem3d(surf,jci1,jci2,ici1,ici2,1,kz,'carbonaer:surf')
        call getmem3d(so4chagct,jci1,jci2, &
                      ici1,ici2,1,kz,'che_common:so4chagct')
      end if
    end subroutine carb_init

    subroutine carb_prepare( )
      implicit none
      integer(ik4) :: i , j , k
      real(rkx) :: kav , pm
      ncon(:,:,:) = d_zero
      surf(:,:,:) = d_zero
      if ( ibchb > 0 ) then
        ! particle mass in kg
        pm = rhobc * mathpi / d_six * reffbc * reffbc * reffbc*1e-18
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! mixing ratio in kg/kg
              kav = max(chib(j,i,k,ibchb)-mintr,d_zero)
              ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
              ! unit of surf m2/kg
              surf(j,i,k) = surf(j,i,k)+mathpi*reffbc*reffbc*kav/pm*1e-12
            end do
          end do
        end do
      end if
      if ( ibchl > 0 ) then
        ! particle mass in kg
        pm = rhobchl * mathpi / d_six * reffbchl*reffbchl*reffbchl*1e-18
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! mixing ratio in kg/kg
              kav = max(chib(j,i,k,ibchl)-mintr,d_zero)
              ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
              ! unit of surf m2/kg
              surf(j,i,k) = surf(j,i,k)+mathpi*reffbchl*reffbchl*kav/pm*1e-12
            end do
          end do
        end do
      end if
      if ( iochb > 0 ) then
        ! particle mass in kg
        pm = rhooc * mathpi / d_six * reffoc*reffoc*reffoc*1e-18
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! mixing ratio in kg/kg
              kav = max(chib(j,i,k,iochb)-mintr,d_zero)
              ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
              ! unit of surf m2/kg
              surf(j,i,k) = surf(j,i,k)+mathpi*reffoc*reffoc*kav/pm*1e-12
            end do
          end do
        end do
      end if
      if ( iochl > 0 ) then
        ! particle mass in kg
        pm = rhoochl * mathpi / d_six * reffochl*reffochl*reffochl*1e-18
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! mixing ratio in kg/kg
              kav = max(chib(j,i,k,iochl)-mintr,d_zero)
              ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
              ! unit of surf m2/kg
              surf(j,i,k) = surf(j,i,k)+mathpi*reffochl*reffochl*kav/pm*1e-12
            end do
          end do
        end do
      end if
      if ( iso4 > 0 ) then
        ! particle mass in kg
        pm = rhoso4 * mathpi / d_six * reffso4 * reffso4 * reffso4*1e-18
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! mixing ratio in kg/kg
              kav = max(chib(j,i,k,iso4)-mintr,d_zero)
              ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
              ! unit of surf m2/kg
              surf(j,i,k) = surf(j,i,k)+mathpi*reffso4*reffso4*kav/pm*1e-12
            end do
          end do
        end do
      end if
    end subroutine carb_prepare

    subroutine aging_carb(i)
      implicit none
      integer, intent(in) :: i
      integer(ik4) :: j , k
      real(rkx) :: agingtend1 , agingtend2 , ksp , icon , arg
      !
      ! aging o carbon species : Conversion from hydrophobic to
      ! hydrophilic: Carbonaceopus species time constant
      ! ( 1.15 day Cooke et al.,1999 )
      !
      if ( ibchb > 0 .and. ibchl > 0 ) then
        if ( carb_aging_control ) then
          do k = 1 , kz
            do j = jci1 , jci2
              !in nm/hr from m/sec
              if ( surf(j,i,k) > d_zero ) then
                arg = so4chagct(j,i,k)/rhoso4/surf(j,i,k)
                icon = arg*3.6e+12_rkx !?????????????
              else
                icon = d_zero
              end if
              arg = min(max(1.0e-3, kcond*icon + kcoag*ncon(j,i,k)),d_one)
              chagct = 3600.0_rkx * d_one/arg
              ksp = max(chib(j,i,k,ibchb)-mintr,d_zero)
              arg = max(min(dt/chagct,25.0_rkx),d_zero)
              agingtend1 = -ksp*(d_one-exp(-arg))/dt
              agingtend2 = -agingtend1
              chiten(j,i,k,ibchb) = chiten(j,i,k,ibchb) + agingtend1
              chiten(j,i,k,ibchl) = chiten(j,i,k,ibchl) + agingtend2
              if ( ichdiag > 0 ) then
                chemdiag(j,i,k,ibchb) = chemdiag(j,i,k,ibchb)+agingtend1*cfdout
                chemdiag(j,i,k,ibchl) = chemdiag(j,i,k,ibchl)+agingtend2*cfdout
              end if
            end do
          end do
        else
          do k = 1 , kz
            do j = jci1 , jci2
              ksp = max(chib(j,i,k,ibchb)-mintr,d_zero)
              arg = max(min(dt/const_chagct,25.0_rkx),d_zero)
              agingtend1 = -ksp*(d_one-exp(-arg))/dt
              agingtend2 = -agingtend1
              chiten(j,i,k,ibchb) = chiten(j,i,k,ibchb) + agingtend1
              chiten(j,i,k,ibchl) = chiten(j,i,k,ibchl) + agingtend2
              if ( ichdiag > 0 ) then
                chemdiag(j,i,k,ibchb) = chemdiag(j,i,k,ibchb)+agingtend1*cfdout
                chemdiag(j,i,k,ibchl) = chemdiag(j,i,k,ibchl)+agingtend2*cfdout
              end if
            end do
          end do
        end if
      end if
      if ( iochb > 0  .and. iochl > 0 ) then
        if ( carb_aging_control ) then
          do k = 1 , kz
            do j = jci1 , jci2
              ! in nm/hr from m/sec
              if ( surf(j,i,k) > d_zero ) then
                arg = so4chagct(j,i,k)/rhoso4/surf(j,i,k)
                icon = arg*3.6e+12_rkx !????????????
              else
                icon = d_zero
              end if
              arg = min(max(1.0e-3, kcond*icon + kcoag*ncon(j,i,k)),d_one)
              chagct = 3600.0_rkx * d_one/arg
              ksp = max(chib(j,i,k,iochb)-mintr,d_zero)
              arg = max(min(dt/chagct,25.0_rkx),d_zero)
              agingtend1 = -ksp*(d_one-exp(-arg))/dt
              agingtend2 = -agingtend1
              chiten(j,i,k,iochb) = chiten(j,i,k,iochb) + agingtend1
              chiten(j,i,k,iochl) = chiten(j,i,k,iochl) + agingtend2
              if ( ichdiag > 0 ) then
                chemdiag(j,i,k,iochb) = chemdiag(j,i,k,iochb)+agingtend1*cfdout
                chemdiag(j,i,k,iochl) = chemdiag(j,i,k,iochl)+agingtend2*cfdout
              end if
            end do
          end do
        else
          do k = 1 , kz
            do j = jci1 , jci2
              ksp = max(chib(j,i,k,iochb)-mintr,d_zero)
              arg = max(min(dt/const_chagct,25.0_rkx),d_zero)
              agingtend1 = -ksp*(d_one-exp(-arg))/dt
              agingtend2 = -agingtend1
              chiten(j,i,k,iochb) = chiten(j,i,k,iochb) + agingtend1
              chiten(j,i,k,iochl) = chiten(j,i,k,iochl) + agingtend2
              if ( ichdiag > 0 ) then
                chemdiag(j,i,k,iochb) = chemdiag(j,i,k,iochb)+agingtend1*cfdout
                chemdiag(j,i,k,iochl) = chemdiag(j,i,k,iochl)+agingtend2*cfdout
              end if
            end do
          end do
        end if
      end if
      if ( ism1 > 0  .and. ism2 > 0 ) then
        do k = 1 , kz
          do j = jci1 , jci2
            ksp = max(chib(j,i,k,ism1)-mintr,d_zero)
            agingtend1 = -ksp*(d_one-exp(-dt/chsmct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,ism1) = chiten(j,i,k,ism1) + agingtend1
            chiten(j,i,k,ism2) = chiten(j,i,k,ism2) + agingtend2
            if ( ichdiag > 0 ) then
              chemdiag(j,i,k,ism1) = chemdiag(j,i,k,ism1) + agingtend1*cfdout
              chemdiag(j,i,k,ism2) = chemdiag(j,i,k,ism2) + agingtend2*cfdout
            end if
          end do
        end do
      end if

    end subroutine aging_carb

end module mod_che_carbonaer

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
