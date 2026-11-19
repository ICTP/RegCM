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

module mod_cu_common

  use mod_intkinds, only : ik4
  use mod_realkinds, only : rkx
  use mod_constants, only : d_10, d_1000, d_half, d_one, d_two, d_zero, egrav
  use mod_dynparam, only : ice1ga, ice1gb, ice2ga, ice2gb, ici1, ici2, idi1, &
    idi2, idynamic, jce1ga, jce1gb, jce2ga, jce2gb, jci1, jci2, jdi1, jdi2, &
    kz, ntr
  use mod_runparams, only : clfrcv, nqx
  use mod_runparams, only : icup, ipptls, icumcloud, ichem
  use mod_memutil, only : getmem
  use mod_regcm_types, only : mod_2_cum

  implicit none

  private

  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_tten
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_uten
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_vten
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: cu_qten
  real(rkx), pointer, contiguous, public, dimension(:,:,:,:) :: cu_chiten
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: avg_ww
  real(rkx), pointer, contiguous, public, dimension(:,:) :: cu_prate
  real(rkx), pointer, contiguous, public, dimension(:,:) :: cu_srate
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_qdetr
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_raincc
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_convpr
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_cldfrc
  ! In cloud convective condensate for radiation [g/m3], zero where the
  ! scheme does not provide it
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cu_cldlwc
  integer(ik4), pointer, contiguous, public, dimension(:,:) :: cu_ktop
  integer(ik4), pointer, contiguous, public, dimension(:,:) :: cu_kbot

  ! which scheme to use
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: cuscheme

  integer(ik4), public :: total_precip_points

  real(rkx), dimension(10) :: cld_profile
  real(rkx), dimension(10) :: fixed_cld_profile
  real(rkx), dimension(10) :: rnum

  real(rkx), parameter :: maxcloud_dp =  60000.0_rkx ! In Pa
  logical, parameter :: addnoise = .false.

  public :: init_mod_cumulus, model_cumulus_cloud
  public :: cu_cloud_fraction, cu_cloud_water

  contains

  subroutine init_mod_cumulus
    implicit none

    if ( any(icup == 4) .or. any(icup == 5) ) then
      if ( idynamic == 3 ) then
        call getmem(cu_uten,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'cumulus:uten')
        call getmem(cu_vten,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'cumulus:vten')
      else
        call getmem(cu_uten,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'cumulus:uten')
        call getmem(cu_vten,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'cumulus:vten')
      end if
    end if
    call getmem(cu_tten,jci1,jci2,ici1,ici2,1,kz,'cumulus:tten')
    call getmem(cu_qten,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cumulus:qten')
    call getmem(cu_cldfrc,jci1,jci2,ici1,ici2,1,kz,'cumulus:cldfrc')
    call getmem(cu_cldlwc,jci1,jci2,ici1,ici2,1,kz,'cumulus:cldlwc')
    call getmem(cu_prate,jci1,jci2,ici1,ici2,'cumulus:prate')
    if ( ipptls > 1 .and. any(icup == 5) ) then
      call getmem(cu_srate,jci1,jci2,ici1,ici2,'cumulus:srate')
    end if
    call getmem(cu_ktop,jci1,jci2,ici1,ici2,'cumulus:ktop')
    call getmem(cu_kbot,jci1,jci2,ici1,ici2,'cumulus:kbot')
    if ( ichem == 1 ) then
      call getmem(cu_chiten,jci1,jci2,ici1,ici2,1,kz,1,ntr,'cumulus:chiten')
      call getmem(cu_convpr,jci1,jci2,ici1,ici2,1,kz,'cumulus:convpr')
    end if
    if ( any(icup == 5) ) then
      call getmem(cu_qdetr,jdi1,jdi2,idi1,idi2,1,kz,'cumulus:qdetr')
      call getmem(cu_raincc,jdi1,jdi2,idi1,idi2,1,kz,'cumulus:raincc')
    end if
    if ( any(icup == 5) .or. any(icup == 6) ) then
      call getmem(avg_ww,jci1,jci2,ici1,ici2,1,kz,'cumulus:avg_ww')
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
      if ( .not. addnoise ) then
        cld_profile = fixed_cld_profile
      end if
    end if
  end subroutine init_mod_cumulus

  subroutine model_cumulus_cloud(m2c)
    implicit none
    type(mod_2_cum), intent(in) :: m2c
    real(rkx) :: akclth, scalep, scalef
    integer(ik4):: i, j, k, ktop, kbot, kclth, ikh
    scalef = (d_one-clfrcv)
    if ( icumcloud <= 1 ) then
      do i = ici1, ici2
        do j = jci1, jci2
          ! The regcm model is top to bottom
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle
          akclth = d_one/real(kclth,rkx)
          do k = ktop, kbot
            cu_cldfrc(j,i,k) = d_one - scalef**akclth
          end do
        end do
      end do
    else if ( icumcloud == 2 ) then
      if ( addnoise ) then
        ! Put 25% noise level. Update cld_profile each time.
        call random_number(rnum)
        cld_profile = (0.75_rkx+(rnum/2.0_rkx))*fixed_cld_profile
      end if
      do i = ici1, ici2
        do j = jci1, jci2
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle
          scalep = min((m2c%pas(j,i,kbot)-m2c%pas(j,i,ktop)) / &
                  maxcloud_dp,d_one)
          do k = ktop, kbot
            ikh = max(1,min(10,int((real(k-ktop+1,rkx)/real(kclth,rkx))*d_10)))
            cu_cldfrc(j,i,k) = cld_profile(ikh)*clfrcv*scalep
          end do
        end do
      end do
    end if
  end subroutine model_cumulus_cloud
  !
  ! Convective cloud fraction seen by radiation, the same for all the mass
  ! flux cumulus schemes, selected by iconvcf (all arguments, no module
  ! state, for the device):
  !   icf    : 0 kfac*ln(1+k2*mfu), 1 Xu and Krueger (1991) fit, 2 as 1
  !            converted to the area fraction of a layer of depth dz
  !   ctype  : 1 deep, 2 shallow, 3 mid level convection
  !   mfu    : grid mean updraft mass flux (kg m-2 s-1)
  !   zag    : height of the level above the ground (m)
  !   dz     : depth of the layer (m)
  !   kdeep, kshal, k2 : kfac_deep, kfac_shal, k2_const (icf = 0)
  ! The fraction is at most 0.6 for deep, 0.2 for shallow and 0.4 for mid
  ! level convection; the grid spacing dependent maximum (clfrcv) is
  ! applied by the cumulus interface.
  !
  pure real(rkx) function cu_cloud_fraction(icf,ctype,mfu,zag,dz, &
                                            kdeep,kshal,k2) result(cf)
!$acc routine seq
    implicit none
    integer(ik4), intent(in) :: icf, ctype
    real(rkx), intent(in) :: mfu, zag, dz, kdeep, kshal, k2
    ! Xu and Krueger (1991), Table 5: convective cloud fraction of thin
    ! (30 mb) layers at 0.8, 3.9 and 7.9 km, fit to the log of the grid
    ! mean cumulus mass flux in mb/h (valid above 0.01 mb/h).
    real(rkx), parameter :: xkz1 = 800.0_rkx
    real(rkx), parameter :: xkz2 = 3900.0_rkx
    real(rkx), parameter :: xkz3 = 7900.0_rkx
    real(rkx), parameter :: xkc01 = 0.0116_rkx
    real(rkx), parameter :: xkc11 = 0.0103_rkx
    real(rkx), parameter :: xkc21 = 0.0025_rkx
    real(rkx), parameter :: xkc02 = 0.0363_rkx
    real(rkx), parameter :: xkc12 = 0.0226_rkx
    real(rkx), parameter :: xkc22 = 0.0054_rkx
    real(rkx), parameter :: xkc03 = 0.0131_rkx
    real(rkx), parameter :: xkc13 = 0.0123_rkx
    real(rkx), parameter :: xkc23 = 0.0030_rkx
    ! kg m-2 s-1 to mb/h
    real(rkx), parameter :: xkmbh = egrav*3600.0_rkx/100.0_rkx
    ! Neggers et al. (2011), JGR 116, D21202: the cover projected by a
    ! cumulus layer of depth dz is (1 + b dz) times its volume fraction
    ! (BOMEX value of b, m-1).
    real(rkx), parameter :: cuovlb = 0.0064_rkx
    real(rkx) :: zmc, zw, zc0, zc1, zc2, zx, zcov
    cf = d_zero
    if ( icf == 1 .or. icf == 2 ) then
      ! Fit coefficients linear in height between the three levels,
      ! constant outside. Below the minimum of the quadratic (about
      ! 0.13 mb/h) the minimum value is used.
      zmc = max(mfu,d_zero)*xkmbh
      if ( zmc > 0.01_rkx ) then
        if ( zag < xkz2 ) then
          zw = max(zag-xkz1,d_zero)/(xkz2-xkz1)
          zc0 = xkc01 + zw*(xkc02-xkc01)
          zc1 = xkc11 + zw*(xkc12-xkc11)
          zc2 = xkc21 + zw*(xkc22-xkc21)
        else
          zw = min((zag-xkz2)/(xkz3-xkz2),d_one)
          zc0 = xkc02 + zw*(xkc03-xkc02)
          zc1 = xkc12 + zw*(xkc13-xkc12)
          zc2 = xkc22 + zw*(xkc23-xkc22)
        end if
        zx = max(log(zmc),-zc1/(d_two*zc2))
        cf = zc0 + zc1*zx + zc2*zx*zx
        if ( icf == 2 ) then
          ! Volume fraction to the area fraction seen by the radiation,
          ! in the symmetric form of Brooks et al. (2005), JAS 62, 2248,
          ! which tends to 1 for full cover.
          zcov = d_one + cuovlb*dz
          cf = d_one/(d_one+(d_one/max(cf,1.0e-6_rkx)-d_one)/zcov)
        end if
      end if
    else
      if ( ctype == 1 ) then
        cf = kdeep*log(d_one+k2*max(mfu,d_zero))
      else if ( ctype == 2 ) then
        cf = kshal*log(d_one+k2*max(mfu,d_zero))
      else
        cf = d_half*(kdeep+kshal)*log(d_one+k2*max(mfu,d_zero))
      end if
    end if
    if ( ctype == 1 ) then
      cf = min(max(d_zero,cf),0.6_rkx)
    else if ( ctype == 2 ) then
      cf = min(max(d_zero,cf),0.2_rkx)
    else
      cf = min(max(d_zero,cf),0.4_rkx)
    end if
  end function cu_cloud_fraction
  !
  ! In cloud convective condensate seen by radiation (g/m3), for the mass
  ! flux cumulus schemes: the updraft condensate per unit volume of the
  ! grid box, rho*a*lu with the updraft area a = mfu/(rho*wu), spread over
  ! the convective cloud fraction cf.
  !   mfu : grid mean updraft mass flux (kg m-2 s-1)
  !   lu  : updraft condensate, liquid and ice (kg/kg)
  !   wu  : updraft vertical velocity (m/s), at least 1 m/s
  !   cf  : convective cloud fraction of the layer
  !
  pure real(rkx) function cu_cloud_water(mfu,lu,wu,cf) result(cwc)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: mfu, lu, wu, cf
    if ( cf > d_zero ) then
      cwc = d_1000 * max(mfu,d_zero) * max(lu,d_zero) / (max(wu,d_one) * cf)
    else
      cwc = d_zero
    end if
  end function cu_cloud_water

end module mod_cu_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
