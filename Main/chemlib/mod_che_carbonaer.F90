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
  use mod_runparams , only : chechgact

  implicit none

  private

  ! Parameter usefull for wet and dry deposition of carbon aerosol
  ! densities in kg/m3

  real(rkx) , public , parameter :: rhobc   = 2000.0_rkx
  real(rkx) , public , parameter :: rhooc   = 1200.0_rkx
  ! higher modes
  real(rkx) , public , parameter , dimension(nchlmax) :: rhobchl = &
                      [ 1600.0_rkx , 1600.0_rkx , 1600.0_rkx ]
  real(rkx) , public , parameter , dimension(nchlmax) :: rhoochl = &
                      [ 1200.0_rkx , 1200.0_rkx , 1200.0_rkx ]
  real(rkx) , public , parameter :: rhoso4 = 1830.0_rkx  !(95-96% H2SO4)
  real(rkx) , public , parameter :: rhosm1 = 1200.0_rkx
  real(rkx) , public , parameter :: rhosm2 = 1200.0_rkx

  ! effctive diameters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? )
  real(rkx) , public , parameter :: reffbc   = 0.05_rkx
  real(rkx) , public , parameter :: reffoc   = 0.2_rkx
  ! main mode
  real(rkx) , public , parameter , dimension(nchlmax) :: reffbchl = &
                         [ 0.3_rkx , 0.3_rkx , 0.3_rkx ]
  real(rkx) , public , parameter , dimension(nchlmax) :: reffochl = &
                         [ 0.3_rkx , 0.3_rkx , 0.3_rkx ]
  ! higher modes
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
  real(rkx) , parameter :: soloc = 0.05_rkx
  ! Higher modes
  real(rkx) , parameter , dimension(nchlmax) :: solbchl = &
                      [ 0.8_rkx , 0.8_rkx , 0.8_rkx ]
  real(rkx) , parameter , dimension(nchlmax) :: solochl = &
                      [ 0.8_rkx , 0.8_rkx , 0.8_rkx ]
  real(rkx) , parameter :: solsm1 = 0.05_rkx
  real(rkx) , parameter :: solsm2 = 0.8_rkx

  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(rkx) , public , dimension(12) :: carbed

  real(rkx) , pointer , dimension(:,:,:) :: ncon
  real(rkx) , pointer , dimension(:,:,:,:) :: save_chagct
  real(rkx) , pointer , dimension(:,:,:,:) :: save_ncon
  ! variable for surface area of aerosol
  real(rkx) , pointer , dimension(:,:,:) :: surf
  real(rkx) , pointer , dimension(:,:,:) :: so4chagct

  public :: solbc , solbchl , soloc , solochl
  public :: solsm1 , solsm2
  public :: ncon , surf
  ! unit of so4chagct = kg kg-1 sec-1
  public :: so4chagct , save_chagct , save_ncon
  public :: carb_init , carb_prepare , aging_carb
  real(rkx) , pointer , dimension(:,:,:,:) :: pp

  contains

    subroutine carb_init( )
      implicit none
      if ( carb_aging_control ) then
        call getmem3d(ncon,jci1,jci2,ici1,ici2,1,kz,'carbonaer:ncon')
        call getmem3d(surf,jci1,jci2,ici1,ici2,1,kz,'carbonaer:surf')
        call getmem3d(so4chagct,jci1,jci2, &
                      ici1,ici2,1,kz,'che_common:so4chagct')
        if ( chechgact ) then
          call getmem4d(save_chagct,jci1,jci2, &
                        ici1,ici2,1,kz,1,ntr,'che_common:save_chagct')
          call getmem4d(save_ncon,jci1,jci2, &
                        ici1,ici2,1,kz,1,ntr,'che_common:save_ncon')
        end if
      end if
      if ( idynamic == 3 ) then
        call assignpnt(chemt,pp)
      else
        call assignpnt(chib,pp)
      end if
    end subroutine carb_init

    subroutine carb_prepare( )
      implicit none
      integer(ik4) :: n
      ncon(:,:,:) = d_zero
      surf(:,:,:) = d_zero
      if ( ibchb > 0 ) then
        call addto_ncon_surf(reffbc,rhobc,ibchb)
      end if
      if ( nbchl > 0 ) then
        ! particle mass in kg
        do n = 1 , nbchl
          call addto_ncon_surf(reffbchl(n),rhobchl(n),ibchl(n))
          end do
      end if
      if ( iochb > 0 ) then
        call addto_ncon_surf(reffoc,rhooc,iochb)
      end if
      if ( nochl > 0 ) then
        ! particle mass in kg
        do n = 1 , nochl
          call addto_ncon_surf(reffochl(n),rhoochl(n),iochl(n))
        end do
      end if
      if ( iso4 > 0 ) then
        ! particle mass in kg
        call addto_ncon_surf(reffso4,rhoso4,iso4)
      end if
    end subroutine carb_prepare

    subroutine aging_carb
      implicit none
      integer(ik4) :: n
      integer(ik4) , dimension(nchlmax+1) :: ids
      !
      ! aging o carbon species : Conversion from hydrophobic to
      ! hydrophilic: Carbonaceopus species time constant
      ! ( 1.15 day Cooke et al.,1999 )
      !
      if ( ibchb > 0 .and. nbchl > 0 ) then
        ids(1) = ibchb
        do n = 1 , nbchl
          ids(n+1) = ibchl(n)
        end do
        if ( carb_aging_control ) then
          call doagingdyn(ids(1:nbchl+1))
        else
          call doaging(ids(1:nbchl+1),const_chagct)
        end if
      end if
      if ( iochb > 0 .and. nochl > 0 ) then
        ids(1) = iochb
        do n = 1 , nochl
          ids(n+1) = iochl(n)
        end do
        if ( carb_aging_control ) then
          call doagingdyn(ids(1:nochl+1))
        else
          call doaging(ids(1:nochl+1),const_chagct)
        end if
      end if
      if ( ism1 > 0 .and. ism2 > 0 ) then
        call doaging([ism1,ism2],chsmct)
      end if

    end subroutine aging_carb

    subroutine addto_ncon_surf(refx,rhox,id)
      implicit none
      real(rkx) , intent(in) :: refx , rhox
      integer(ik4) , intent(in) :: id
      integer(ik4) :: i , j , k
      real(rkx) :: kav , pm
        ! particle mass in kg
        pm = rhox * mathpi / d_six * refx * refx * refx*1e-18
        if ( idynamic == 3 ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! mixing ratio in kg/kg
                kav = max(pp(j,i,k,id)-mintr,d_zero)
                ncon(j,i,k) = ncon(j,i,k)+kav*crhob3d(j,i,k)/pm*1e-6
                if ( chechgact ) then
                  save_ncon(j,i,k,id) = kav*crhob3d(j,i,k)/pm*1e-6
                end if
                ! unit of surf m2/kg
                surf(j,i,k) = surf(j,i,k)+mathpi*refx*refx*kav/pm*1e-12
              end do
            end do
          end do
        else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! mixing ratio in kg/kg
                kav = max(pp(j,i,k,id)-mintr,d_zero)
                ncon(j,i,k) = ncon(j,i,k)+kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
                if ( chechgact ) then
                  save_ncon(j,i,k,id) = kav/cpsb(j,i)*crhob3d(j,i,k)/pm*1e-6
                end if
                ! unit of surf m2/kg
                surf(j,i,k) = surf(j,i,k)+mathpi*reffso4*reffso4*kav/pm*1e-12
              end do
            end do
          end do
        end if
      end subroutine addto_ncon_surf

      subroutine doagingdyn(ids)
        implicit none
        integer(ik4) , dimension(:) , intent(in) :: ids
        integer(ik4) :: i , j , k , n , b1 , b2
        real(rkx) , dimension(jci1:jci2,ici1:ici2) :: agingtend
        real(rkx) :: kav , arg , icon

        do n = 1 , size(ids) - 1
          b1 = ids(n)
          b2 = ids(n+1)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! in nm/hr from m/sec
                if ( surf(j,i,k) > d_zero ) then
                  arg = so4chagct(j,i,k)/rhoso4/surf(j,i,k)
                  icon = arg*3.6e+12_rkx
                else
                  icon = d_zero
                end if
                arg = min(max(1.0e-3_rkx, &
                          kcond*icon + kcoag*ncon(j,i,k)),d_one)
                chagct = 3600.0_rkx * d_one/arg
                if ( chechgact ) then
                  save_chagct(j,i,k,b1) = chagct
                end if
                kav = max(pp(j,i,k,b1)-mintr,d_zero)
                arg = max(min(dt/chagct,25.0_rkx),d_zero)
                agingtend(j,i) = -kav*(d_one-exp(-arg))/dt
                chiten(j,i,k,b1) = chiten(j,i,k,b1) + agingtend(j,i)
                chiten(j,i,k,b2) = chiten(j,i,k,b2) - agingtend(j,i)
              end do
            end do
            if ( ichdiag > 0 ) then
              do i = ici1 , ici2
                do j = jci1 , jci2
                  chemdiag(j,i,k,b1) = chemdiag(j,i,k,b1) + &
                                agingtend(j,i)*cfdout
                  chemdiag(j,i,k,b2) = chemdiag(j,i,k,b2) - &
                                agingtend(j,i)*cfdout
                end do
              end do
            end if
          end do
        end do
      end subroutine doagingdyn

      subroutine doaging(ids,sm)
        implicit none
        integer(ik4) , dimension(:) , intent(in) :: ids
        real(rkx) , intent(in) :: sm
        integer(ik4) :: i , j , k , n , b1 , b2
        real(rkx) , dimension(jci1:jci2,ici1:ici2) :: agingtend
        real(rkx) :: kav

        do n = 1 , size(ids) - 1
          b1 = ids(n)
          b2 = ids(n+1)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                kav = max(pp(j,i,k,b1)-mintr,d_zero)
                agingtend(j,i) = -kav*(d_one-exp(-dt/sm))/dt
                chiten(j,i,k,b1) = chiten(j,i,k,b1) + agingtend(j,i)
                chiten(j,i,k,b2) = chiten(j,i,k,b2) - agingtend(j,i)
              end do
            end do
            if ( ichdiag > 0 ) then
              do i = ici1 , ici2
                do j = jci1 , jci2
                  chemdiag(j,i,k,b1) = chemdiag(j,i,k,b1) + &
                                agingtend(j,i)*cfdout
                  chemdiag(j,i,k,b2) = chemdiag(j,i,k,b2) - &
                                agingtend(j,i)*cfdout
                end do
              end do
            end if
          end do
        end do
      end subroutine doaging

end module mod_che_carbonaer

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
