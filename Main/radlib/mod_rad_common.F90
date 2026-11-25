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

module mod_rad_common

  use mod_intkinds, only : ik4
  use mod_realkinds, only : is_nan, rk4, rk8, rkx
  use mod_dynparam, only : nspi

  implicit none

  public

  ! absnxt  - Nearest layer absorptivities
  ! abstot  - Non-adjacent layer absorptivites
  ! emstot  - Total emissivity

  ! Those need to be saved in output file
  real(rk8), pointer, contiguous, dimension(:,:,:) :: o3prof => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:,:)  :: gasabsnxt => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:,:)  :: gasabstot => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: gasemstot => null( )
  ! Column inputs of the last absorptivity/emissivity update, saved in
  ! the restart file instead of the tables, which are rebuilt from them
  ! (5*kz+8 values per column, against (kz+1)**2+5*kz+1 for the tables).
  ! Layout of the last index: surface temperature, tropopause pressure,
  ! co2, ch4, n2o, cfc11, cfc12 before the stratospheric scaling, then
  ! temperature, water vapour, mid level pressure and ozone (kz values
  ! each) and interface pressure (kz+1 values).
  real(rk8), pointer, contiguous, dimension(:,:,:) :: absemin => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: taucldsp => null( )
  logical, save :: doabsems, dolw, dosw
  integer(ik4), save :: ichso4, ichbc, ichoc
  integer(ik4), save :: kclimh, kth, ktf, ksf, kclimf
  real(rk8) :: chfrovrradfr ! chfrq/rafrq

  real(rk8), dimension(nspi), parameter :: wavmin = &
    [ 0.200_rk8, 0.245_rk8, 0.265_rk8, 0.275_rk8, 0.285_rk8 ,&
      0.295_rk8, 0.305_rk8, 0.350_rk8, 0.640_rk8, 0.700_rk8 ,&
      0.701_rk8, 0.701_rk8, 0.701_rk8, 0.701_rk8, 0.702_rk8 ,&
      0.702_rk8, 2.630_rk8, 4.160_rk8, 4.160_rk8 ]

  real(rk8), dimension(nspi), parameter :: wavmax = &
    [ 0.245_rk8, 0.265_rk8, 0.275_rk8, 0.285_rk8, 0.295_rk8 ,&
      0.305_rk8, 0.350_rk8, 0.640_rk8, 0.700_rk8, 5.000_rk8 ,&
      5.000_rk8, 5.000_rk8, 5.000_rk8, 5.000_rk8, 5.000_rk8 ,&
      5.000_rk8, 2.860_rk8, 4.550_rk8, 4.550_rk8 ]

  real(rk8), dimension(14), parameter :: wavnm1 = &
    [ 2600.0_rk8, 3250.0_rk8, 4000.0_rk8, 4650.0_rk8, 5150.0_rk8, &
      6150.0_rk8, 7700.0_rk8, 8050.0_rk8,12850.0_rk8,16000.0_rk8, &
     22650.0_rk8,29000.0_rk8,38000.0_rk8,  820.0_rk8 ]

  real(rk8), dimension(14), parameter :: wavnm2 = &
    [  3250.0_rk8, 4000.0_rk8, 4650.0_rk8, 5150.0_rk8, 6150.0_rk8, &
       7700.0_rk8, 8050.0_rk8,12850.0_rk8,16000.0_rk8,22650.0_rk8, &
      29000.0_rk8,38000.0_rk8,50000.0_rk8, 2600.0_rk8 ]

  ! Minimum in cloud water path (g/m2) of a layer for a cloud to be seen
  ! by the radiation, the same for all the schemes
  real(rk8), parameter :: cldwpmin = 0.01_rk8

  ! Ice to liquid conversion of the generalized effective size De of
  ! Fu (1996) and Sun (2001) to the ice effective radius: re = 3*sqrt(3)/8 De
  real(rk8), parameter :: dge2rei = 0.64952_rk8

  contains

  !
  ! Cloud optical properties shared by the radiation schemes (CCM and
  ! RRTMG), so that the short and long wave see the same cloud.
  !
  ! Ice fraction of the cloud water seen by radiation. With the
  ! microphysics condensate (lgrid) the ratio of grid ice to total; where
  ! the grid has no condensate (convective only cloud) the temperature
  ! ramp of the mixed phase saturation (satliqfrac in saturation.inc).
  !
  pure real(rk8) function cloud_fice(lgrid,ql,qi,t) result(fice)
!$acc routine seq
    use mod_constants, only : rtice, rtwat, rtwat_rtice_r
    implicit none
    logical, intent(in) :: lgrid
    real(rk8), intent(in) :: ql, qi, t
    real(rk8), parameter :: qmin = 1.0e-11_rk8
    if ( lgrid .and. ql+qi > qmin ) then
      fice = max(qi,0.0_rk8)/(max(ql,0.0_rk8)+max(qi,0.0_rk8))
    else
      fice = 1.0_rk8 - ((max(real(rtice,rk8),min(real(rtwat,rk8),t)) - &
                     real(rtice,rk8))*real(rtwat_rtice_r,rk8))**2
    end if
  end function cloud_fice
  !
  ! Liquid effective radius (micron) from the in cloud liquid water
  ! content lwc (kg/m3), Martin et al. (1994): droplet number 300 cm-3
  ! and k = 0.67 over land, 100 cm-3 and k = 0.80 over sea, or the
  ! given droplet number cdnc (cm-3, > 0).
  !
  pure real(rk8) function cloud_rel(lwc,island,cdnc) result(rel)
!$acc routine seq
    use mod_constants, only : rhoh2o, mathpi
    implicit none
    real(rk8), intent(in) :: lwc, cdnc
    logical, intent(in) :: island
    real(rk8) :: kp, nd
    if ( island ) then
      kp = 0.67_rk8
      nd = 300.0e6_rk8
    else
      kp = 0.80_rk8
      nd = 100.0e6_rk8
    end if
    if ( cdnc > 0.0_rk8 ) nd = cdnc*1.0e6_rk8
    rel = ((3.0_rk8*max(lwc,1.0e-10_rk8)) / &
           (4.0_rk8*real(mathpi,rk8)*real(rhoh2o,rk8)*kp*nd))**(1.0_rk8/3.0_rk8)
    rel = max(4.0_rk8,min(20.0_rk8,rel*1.0e6_rk8))
  end function cloud_rel
  !
  ! Liquid effective radius (micron) in the form of the NCAR CAM3
  ! (Collins et al., 2004): 14 micron over sea and sea ice, over land
  ! ramping to the sea value between 0 and -20 C. The land value is
  ! 11 micron (CAM3 8 micron), closer to the MODIS land radius.
  !
  pure real(rk8) function cloud_rel_cam(t,island) result(rel)
!$acc routine seq
    implicit none
    real(rk8), intent(in) :: t
    logical, intent(in) :: island
    real(rk8), parameter :: rliqland = 11.0_rk8
    real(rk8), parameter :: rliqocean = 14.0_rk8
    if ( island ) then
      rel = rliqland + (rliqocean-rliqland) * &
            min(1.0_rk8,max(0.0_rk8,(273.16_rk8-t)*0.05_rk8))
    else
      rel = rliqocean
    end if
  end function cloud_rel_cam
  !
  ! Ice effective radius (micron) from the in cloud ice water content iwc
  ! (g/m3), temperature (K) and latitude (degrees): Sun (2001) generalized
  ! effective size with the minimum size 20 + 40 cos(lat) micron, as in
  ! the ECMWF IFS, converted to the radius with dge2rei.
  !
  pure real(rk8) function cloud_rei(iwc,t,lat) result(rei)
!$acc routine seq
    use mod_constants, only : degrad
    implicit none
    real(rk8), intent(in) :: iwc, t, lat
    real(rk8) :: tk, x, de
    tk = min(273.15_rk8,max(193.15_rk8,t))
    x = max(iwc,1.0e-7_rk8)
    de = (1.2351_rk8+0.0105_rk8*(tk-273.15_rk8)) * &
         (45.8966_rk8*x**0.2214_rk8 + 0.7957_rk8*x**0.2535_rk8*(tk-83.15_rk8))
    de = max(20.0_rk8+40.0_rk8*cos(lat*real(degrad,rk8)),min(155.0_rk8,de))
    rei = dge2rei*de
  end function cloud_rei
  !
  ! The cloud seen by the radiation schemes, CCM and RRTMG, so that they
  ! differ only in the optics. For one layer, from the cloud fraction cfr
  ! and the in cloud condensate cwc (g/m3) of the large scale and
  ! convective clouds (cldfrac), the layer depth dz (m), temperature t
  ! (K), grid mean cloud liquid ql and ice qi (kg/kg), land flag island,
  ! latitude lat (degrees) and droplet number cdnc (cm-3, 0 for the
  ! climatological radius):
  !   cf   : cloud fraction, zero if the in cloud water path is below
  !          cldwpmin
  !   cwp  : in cloud water path, liquid and ice (g/m2)
  !   fice : ice fraction of cwp
  !   rel  : liquid effective radius (micron), CAM3, or Martin et al.
  !          (1994) with the droplet number cdnc
  !   rei  : ice effective radius (micron), Sun (2001)
  ! The radii have filler values where the cloud fraction is at most
  ! lowcld.
  !
  pure subroutine cloud_layer(cfr,cwc,dz,t,ql,qi,island,lat,cdnc, &
                              cf,cwp,fice,rel,rei)
!$acc routine seq
    use mod_constants, only : lowcld
    implicit none
    real(rk8), intent(in) :: cfr, cwc, dz, t, ql, qi, lat, cdnc
    logical, intent(in) :: island
    real(rk8), intent(out) :: cf, cwp, fice, rel, rei
    real(rk8) :: lwc
    fice = cloud_fice(.true.,ql,qi,t)
    cwp = cwc*dz
    if ( cwp > cldwpmin ) then
      cf = cfr
    else
      cwp = 0.0_rk8
      cf = 0.0_rk8
    end if
    if ( cf > real(lowcld,rk8) ) then
      rel = cloud_rel_cam(t,island)
      rei = cloud_rei(cwp*fice/dz,t,lat)
      if ( cdnc > 0.0_rk8 ) then
        ! in cloud liquid water content (kg/m3)
        lwc = cwp*(1.0_rk8-fice)/dz*0.001_rk8
        if ( lwc >= 1.0e-6_rk8 ) rel = cloud_rel(lwc,island,cdnc)
      end if
    else
      ! filler for no clouds
      rel = 8.5_rk8
      rei = 20.0_rk8
    end if
  end subroutine cloud_layer
  !
  ! Cloud droplet number (cm-3) from the sulfate mass mixing ratio so4
  ! (kg/kg) and the air density rho (kg/m3), Qian et al. (1999), the
  ! simple sulfate indirect effect; zero (climatological radius) below
  ! 0.1 microgram/m3.
  !
  pure real(rk8) function cloud_cdnc_so4(rho,so4) result(cdnc)
!$acc routine seq
    implicit none
    real(rk8), intent(in) :: rho, so4
    real(rk8) :: aerc
    aerc = rho*so4*1.0e9_rk8
    if ( aerc > 0.1_rk8 ) then
      cdnc = 90.7_rk8*aerc**0.45_rk8 + 23.0_rk8
    else
      cdnc = 0.0_rk8
    end if
  end function cloud_cdnc_so4

end module mod_rad_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
