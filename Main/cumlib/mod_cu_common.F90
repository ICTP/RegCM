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

module mod_cu_common

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams , only : clfrcv , icumcloud , icup , ichem , dt , nqx
  use mod_mppparam , only : ma
  use mod_humid , only : clwfromt
  use mod_memutil
  use mod_regcm_types

  implicit none

  public

  real(rkx) , pointer , dimension(:,:,:) :: cu_tten
  real(rkx) , pointer , dimension(:,:,:) :: avg_tten
  real(rkx) , pointer , dimension(:,:,:) :: cu_uten
  real(rkx) , pointer , dimension(:,:,:) :: cu_vten
  real(rkx) , pointer , dimension(:,:,:) :: cu_utenx
  real(rkx) , pointer , dimension(:,:,:) :: cu_vtenx
  real(rkx) , pointer , dimension(:,:,:) :: avg_uten
  real(rkx) , pointer , dimension(:,:,:) :: avg_vten
  real(rkx) , pointer , dimension(:,:,:,:) :: cu_qten
  real(rkx) , pointer , dimension(:,:,:,:) :: avg_qten
  real(rkx) , pointer , dimension(:,:,:,:) :: cu_chiten
  real(rkx) , pointer , dimension(:,:,:,:) :: avg_chiten
  real(rkx) , pointer , dimension(:,:) :: cu_prate
  real(rkx) , pointer , dimension(:,:,:) :: cu_qdetr
  real(rkx) , pointer , dimension(:,:,:) :: cu_raincc
  real(rkx) , pointer , dimension(:,:,:) :: cu_convpr
  real(rkx) , pointer , dimension(:,:,:) :: cu_cldfrc
  real(rkx) , pointer , dimension(:,:,:) :: cu_cldlwc
  integer(ik4) , pointer , dimension(:,:) :: cu_ktop
  integer(ik4) , pointer , dimension(:,:) :: cu_kbot

  real(rkx) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
  real(rkx) :: rdt

  integer(ik4) , pointer , dimension(:,:) :: cuscheme ! which scheme to use
  integer(ik4) :: total_precip_points

  real(rkx) , dimension(10) :: cld_profile
  real(rkx) , dimension(10) :: fixed_cld_profile
  real(rkx) , dimension(10) :: rnum

  real(rkx) , parameter :: maxcloud_dp =  60000.0_rkx ! In Pa
  logical , parameter :: addnoise = .false.

  contains

  subroutine init_mod_cumulus
    implicit none
    integer(ik4) , dimension(:) , allocatable:: iseed
    integer :: k , nseed
    real(rk4) :: cputime

    rdt = d_one / dt

    call getmem3d(cu_tten,jci1,jci2,ici1,ici2,1,kz,'cumulus:tten')
    call getmem3d(avg_tten,jci1,jci2,ici1,ici2,1,kz,'cumulus:avg_tten')
    call getmem3d(cu_uten,jci1,jci2,ici1,ici2,1,kz,'cumulus:uten')
    call getmem3d(cu_vten,jci1,jci2,ici1,ici2,1,kz,'cumulus:vten')
    if ( any(icup == 5) ) then
      call getmem3d(avg_uten,jci1,jci2,ici1,ici2,1,kz,'cumulus:avg_uten')
      call getmem3d(avg_vten,jci1,jci2,ici1,ici2,1,kz,'cumulus:avg_vten')
    end if
    call getmem3d(cu_utenx,jci1ga,jci2ga,ici1ga,ici2ga,1,kz,'cumulus:utenx')
    call getmem3d(cu_vtenx,jci1ga,jci2ga,ici1ga,ici2ga,1,kz,'cumulus:vtenx')
    call getmem4d(cu_qten,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cumulus:qten')
    call getmem4d(avg_qten,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cumulus:avg_qten')
    call getmem3d(cu_cldfrc,jci1,jci2,ici1,ici2,1,kz,'cumulus:cldfrc')
    call getmem3d(cu_cldlwc,jci1,jci2,ici1,ici2,1,kz,'cumulus:cldlwc')
    call getmem2d(cu_prate,jci1,jci2,ici1,ici2,'cumulus:prate')
    call getmem2d(cu_ktop,jci1,jci2,ici1,ici2,'cumulus:ktop')
    call getmem2d(cu_kbot,jci1,jci2,ici1,ici2,'cumulus:kbot')
    if ( ichem == 1 ) then
      call getmem4d(cu_chiten,jci1,jci2,ici1,ici2,1,kz,1,ntr,'cumulus:chiten')
      call getmem4d(avg_chiten,jci1,jci2, &
                               ici1,ici2,1,kz,1,ntr,'cumulus:avgchiten')
      call getmem3d(cu_convpr,jci1,jci2,ici1,ici2,1,kz,'cumulus:convpr')
    end if
    if ( any(icup == 5) ) then
      call getmem3d(cu_qdetr,jdi1,jdi2,idi1,idi2,1,kz,'cumulus:qdetr')
      call getmem3d(cu_raincc,jdi1,jdi2,idi1,idi2,1,kz,'cumulus:raincc')
    end if

    if ( icumcloud == 2 ) then
      !
      ! Free hand draw of a generic ten layer cumulus cloud shape.
      !
      fixed_cld_profile(1)  = 0.130_rkx
      fixed_cld_profile(2)  = 0.125_rkx
      fixed_cld_profile(3)  = 0.120_rkx
      fixed_cld_profile(4)  = 0.080_rkx
      fixed_cld_profile(5)  = 0.080_rkx
      fixed_cld_profile(6)  = 0.080_rkx
      fixed_cld_profile(7)  = 0.085_rkx
      fixed_cld_profile(8)  = 0.085_rkx
      fixed_cld_profile(9)  = 0.105_rkx
      fixed_cld_profile(10) = 0.110_rkx
      if ( addnoise ) then
        call random_seed(size=nseed)
        call cpu_time(cputime)
        allocate(iseed(nseed))
        iseed = int(cputime) + 37*(/(k-1,k=1,nseed)/)
        call random_seed(put=iseed)
        deallocate(iseed)
      else
        cld_profile = fixed_cld_profile
      end if
    end if
  end subroutine init_mod_cumulus

  subroutine model_cumulus_cloud(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    real(rkx) :: akclth , scalep , scalef
    integer(ik4):: i , j , k , ktop , kbot , kclth , ikh
    scalef = (d_one-clfrcv)
    if ( icumcloud <= 1 ) then
      iloop1: &
      do i = ici1 , ici2
        jloop1: &
        do j = jci1 , jci2
          if ( cuscheme(j,i) /= 1 .and. cuscheme(j,i) /= 3 ) cycle jloop1
          ! The regcm model is top to bottom
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop1
          akclth = d_one/real(kclth,rkx)
          do k = ktop , kbot
            cu_cldfrc(j,i,k) = d_one - scalef**akclth
          end do
        end do jloop1
      end do iloop1
    else if ( icumcloud == 2 ) then
      if ( addnoise ) then
        ! Put 25% noise level. Update cld_profile each time.
        call random_number(rnum)
        cld_profile = (0.75_rkx+(rnum/2.0_rkx))*fixed_cld_profile
      end if
      iloop3: &
      do i = ici1 , ici2
        jloop3: &
        do j = jci1 , jci2
          if ( cuscheme(j,i) /= 1 .and. cuscheme(j,i) /= 3 ) cycle jloop3
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop3
          scalep = min((m2c%pas(j,i,kbot)-m2c%pas(j,i,ktop)) / &
                  maxcloud_dp,d_one)
          do k = ktop , kbot
            ikh = max(1,min(10,int((real(k-ktop+1,rkx)/real(kclth,rkx))*d_10)))
            cu_cldfrc(j,i,k) = cld_profile(ikh)*clfrcv*scalep
          end do
        end do jloop3
      end do iloop3
    end if
    if ( icumcloud == 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          if ( ktop > 1 .and. kbot > 1 ) then
            do k = ktop , kbot
              if ( cuscheme(j,i) == 4 ) then
                cu_cldlwc(j,i,k) = 0.5e-4_rkx
              else
                cu_cldlwc(j,i,k) = 0.3e-3_rkx
              end if
            end do
          end if
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( cu_cldfrc(j,i,k) > 0.001_rkx ) then
              cu_cldlwc(j,i,k) = clwfromt(m2c%tas(j,i,k))
            else
              cu_cldlwc(j,i,k) = d_zero
            end if
          end do
        end do
      end do
    end if
  end subroutine model_cumulus_cloud

end module mod_cu_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
