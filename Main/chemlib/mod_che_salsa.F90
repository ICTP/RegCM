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

module mod_che_salsa

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams , only : dtsec 
  !use mod_runparams , only : rcmtimer
  use mod_che_common
  use mod_che_indices
  use mod_che_species
! FAB problem with ncld parameter, rename 
!  USE mo_submctl, ONLY : nbins,ncld,nprc,pi6,          &
!                               rhowa, rhosu, rhobc, rhooc,   &
!                               rhono, rhonh, rhoss, rhodu, fn2a
  USE mo_submctl
  USE mo_salsa, ONLY : salsa
  USE mo_salsa_init
  USE mo_salsa_properties, ONLY : equilibration
  USE class_componentIndex, ONLY : ComponentIndex, GetIndex, GetNcomp, IsUsed

  implicit none

  private
    ! grid points for SALSA
  integer(ik4)   :: kproma 
  INTEGER(ik4)  :: kbdim 
  INTEGER(ik4)  :: klev 
  INTEGER(ik4)  :: krow

!  REAL(rkx), PARAMETER :: init_rh(kbdim,klev) = 0.3_dp

  ! -- Local hydrometeor properties (set up in aero initialize)
  TYPE(t_section), ALLOCATABLE, SAVE :: cloud(:,:,:) ! cloud properties
  TYPE(t_section), ALLOCATABLE, SAVE :: aero(:,:,:)  ! Aerosol properties
  TYPE(t_section), ALLOCATABLE, SAVE :: precp(:,:,:) ! Precipitation properties

  TYPE(ComponentIndex) :: prtcl ! Contains "getIndex" which gives the index for a


  ! -- Local gas compound tracers [# m-3]
  REAL(rkx), allocatable :: zgso4(:,:),   &
              zghno3(:,:),  &
              zgnh4(:,:),   &
              zgocnv(:,:),  &
              zgocsv(:,:),  &
              in_t(:,:),    & 
              in_p(:,:),    &
              in_rv(:,:),   &
              in_rs(:,:),   &
              in_w(:,:) ,   &
              actd(:,:)  


  public :: run_salsa, init_salsa 

  contains 

  subroutine init_salsa()
    implicit none

    kproma = (jci2-jci1+1)*(ici2-ici1+1)
    kbdim = kproma
    klev = kz 
    allocate(in_t(1:kbdim,1:klev))    
    allocate(in_p(1:kbdim,1:klev))    
    allocate(in_rv(1:kbdim,1:klev))   
    allocate(in_rs(1:kbdim,1:klev))   
    allocate(in_w(1:kbdim,1:klev)) 
    
    call salsa_initialize()
  end subroutine init_salsa

  subroutine run_salsa()
    implicit none
    integer(ik4) :: i,j,k,n,nc,vc
! set input
   do k = 1 , kz
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        in_t(n,k) = ctb3d(j,i,k)
        in_p(n,k) = cpb3d(j,i,k)
        in_rv(n,k)= cqxb3d(j,i,k,iqv)
        in_rs(n,k) =  cqsb3d(j,i,k) 
        in_w(n,k) = cwpx3d(j,i,k) 

        IF (IsUsed(prtcl,'SO4')) THEN
                nc = GetIndex(prtcl,'SO4')
                vc = 1
             !   str = (nc-1)*nbins+1
             !   end = nc*nbins
                aero(n,k,1:nbins)%volc(vc) = 0.
        END IF

        n = n + 1
      end do
    end do
   end do

   print* , 'FAB', in_t (50, kz), in_p(50,kz), in_rv(50,kz),in_rs(50,kz), in_w(50,kz) 
 

!    actd(1:kproma,:,:)%numc = 0._dp
!    aero(1:kproma,:,:)%numc = 0._dp
!    cloud(1:kproma,:,:)%numc = 0._dp
!    precp(1:kproma,:,:)%numc = 0._dp
!    DO ss = 1,8 !GetNcomp(prtcl)+1  !!!! FIXED, should be 1,8

!       actd(1:kproma,:,:)%volc(ss) = 0._dp
!       aero(1:kproma,:,:)%volc(ss) = 0._dp
!       cloud(1:kproma,:,:)%volc(ss) = 0._dp
!       precp(1:kproma,:,:)%volc(ss) = 0._dp

!    END DO

    ! Set the SALSA runtime config (saisiko hoidettua tehokkaammin?)
!    CALL set_salsa_runtime(prunmode)
!
!            CALL salsa(kproma, kbdim,  klev,   krow,          &
!                        in_p,   in_rv,  in_rs,  in_t, tstep,   &
!                        zgso4,  zgocnv, zgocsv, zghno3,        &
!                        zgnh3,  aero,   cloud,  precp,         &
!                        actd,   in_w,   dbg3,   prtcl          )



  end subroutine run_salsa

end module mod_che_salsa
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
