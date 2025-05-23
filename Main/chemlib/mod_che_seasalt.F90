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

module mod_che_seasalt

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common

  implicit none

  private

  ! sea-salt density
  real(rkx), parameter :: rhosslt = 1020.0_rkx
  real(rkx), dimension(sbin,2) :: ssltbsiz

  data ssltbsiz /0.05_rkx, 1.0_rkx, 1.0_rkx, 10.0_rkx/

  real(rkx), dimension(sbin) :: ssltbed

  ! solubility of od dust aer for param of giorgi and chameides
  real(rkx), dimension(sbin) :: solsslt

  data ssltbed /0.6_rkx, 6.0_rkx/
  data solsslt /0.8_rkx, 0.8_rkx/

  public :: sea_salt, rhosslt, ssltbsiz, solsslt, ssltbed

  contains
!
!  FAB : NB THIS ROUTINE IS NOT OPTIMIZED AT ALL
!
  !************************************************************
  !*                                                        ***
  !* Purpose:                                               ***
  !* Compute ocean to atmosphere emissions of seasalt mass, ***
  !* number, and surface area for multiple modes            ***
  !*                                                        ***
  !*                                                        ***
  !* Method:                                                ***
  !* Integrates Monahan seasalt emission parameterization   ***
  !* over size range appropriate for each mode.             ***
  !* Applies empirical correction (based on ODowd           ***
  !* measurements) at Dp < 0.2 micrometer.                  ***
  !************************************************************

  subroutine sea_salt(wind10,ivegcov)

    implicit none

    integer(ik4), dimension(jci1:jci2,ici1:ici2), intent(in) :: ivegcov
    real(rkx), dimension(jci1:jci2,ici1:ici2), intent(in) :: wind10

    real(rkx), dimension(sbin) :: seasalt_flx
    real(rkx) :: dumu10
    real(rkx) :: dplo_acc, dphi_acc, dplo_cor, dphi_cor
    real(rkx) :: seasalt_emfac_mascor, seasalt_emfac_numcor
    real(rkx) :: seasalt_emfac_masacc, seasalt_emfac_numacc
    integer(ik4) :: i, j

    real(rkx) :: specmw_seasalt_amode
    real(rkx) :: dum_mw

    specmw_seasalt_amode = 1.0_rkx

    dplo_acc = ssltbsiz(1,1)*1.0e-4_rkx
    dphi_acc = ssltbsiz(1,2)*1.0e-4_rkx
    dplo_cor = ssltbsiz(2,1)*1.0e-4_rkx
    dphi_cor = ssltbsiz(2,2)*1.0e-4_rkx

    !************************************************************
    !* initialization -- compute seasalt emissions factors    ***
    !* for accumulation and coarse modes                      ***
    !************************************************************

    seasalt_emfac_numacc = d_zero
    seasalt_emfac_masacc = d_zero
    seasalt_emfac_numcor = d_zero
    seasalt_emfac_mascor = d_zero

    call seasalt_emit(1,dplo_acc,dphi_acc,seasalt_emfac_numacc, &
                       seasalt_emfac_masacc)

    call seasalt_emit(1,dplo_cor,dphi_cor,seasalt_emfac_numcor, &
                      seasalt_emfac_mascor)

    !************************************************************
    !* convert number factors from (#/m2/s) to ((k#)/m2/s)    ***
    !************************************************************

    seasalt_emfac_numacc = seasalt_emfac_numacc * 1.0e-3_rkx
    seasalt_emfac_numcor = seasalt_emfac_numcor * 1.0e-3_rkx

    !************************************************************
    !* convert mass factors from (g/m2/s) to (kmol/m2/s)     ****
    !************************************************************

    dum_mw = specmw_seasalt_amode
    seasalt_emfac_masacc = seasalt_emfac_masacc * (1.0e-3_rkx/dum_mw)
    seasalt_emfac_mascor = seasalt_emfac_mascor * (1.0e-3_rkx/dum_mw)

    ! print*, seasalt_emfac_masacc, seasalt_emfac_mascor

    do i = ici1, ici2
      do j = jci1, jci2
        if ( ivegcov(j,i) > 0 ) cycle

        seasalt_flx(1) = seasalt_emfac_masacc
        seasalt_flx(2) = seasalt_emfac_mascor

        if ( wind10(j,i) > d_zero ) then
          dumu10 = max( d_zero, min( 100.0_rkx, wind10(j,i) ))
          dumu10 = dumu10**3.41_rkx
          seasalt_flx(1) = seasalt_flx(1) * dumu10
          seasalt_flx(2) = seasalt_flx(2) * dumu10
        else
          seasalt_flx(1) = d_zero
          seasalt_flx(2) = d_zero
        end if

        if ( idynamic == 3 ) then
          ! calculate the source tendancy
          chiten(j,i,kz,isslt(1)) = chiten(j,i,kz,isslt(1)) + &
                  seasalt_flx(1)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
          chiten(j,i,kz,isslt(2)) = chiten(j,i,kz,isslt(2)) + &
                  seasalt_flx(2)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
        else
          ! calculate the source tendancy
          chiten(j,i,kz,isslt(1)) = chiten(j,i,kz,isslt(1)) + &
                  seasalt_flx(1)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
          chiten(j,i,kz,isslt(2)) = chiten(j,i,kz,isslt(2)) + &
                  seasalt_flx(2)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
        end if

        ! diagnostic source
        cemtrac(j,i,isslt(1)) = cemtrac(j,i,isslt(1)) + &
                 seasalt_flx(1)*cfdout
        cemtrac(j,i,isslt(2)) = cemtrac(j,i,isslt(2)) + &
                 seasalt_flx(2)*cfdout
      end do
    end do

  end subroutine sea_salt
!
!===================================================================
!
  subroutine seasalt_emit(ireduce_smallr_emit,dpdrylo_cm, &
                          dpdryhi_cm,emitfact_numb,emitfact_mass)
    implicit none

    integer(ik4), intent(in) :: ireduce_smallr_emit
    real(rkx), intent(in) :: dpdrylo_cm, dpdryhi_cm
    real(rkx), intent(out) :: emitfact_numb, emitfact_mass

    integer(ik4) :: nsub_bin, isub_bin
    real(rkx) :: drydens, drydens_f
    real(rkx) :: relhum
    real(rkx) :: rdry_star, sigmag_star
    real(rkx) :: dumsum_na, dumsum_ma
    real(rkx) :: rdrylowermost, rdryuppermost
    real(rkx) :: rdrylo, rdryhi
    real(rkx) :: alnrdrylo, dlnrdry
    real(rkx) :: rdrybb, rwetbb
    real(rkx) :: rdryaa, rwetaa
    real(rkx) :: rdry_cm, rwet_cm
    real(rkx) :: rdry, rwet, drwet
    real(rkx) :: dum,xmdry, dumb,dumexpb
    real(rkx) :: dumadjust, df0dlnrdry
    real(rkx) :: df0drwet
    real(rkx) :: df0dlnrdry_star

    !************************************************************
    !*  c1-c4 are constants for sea-salt hygroscopi!growth  ****
    !*  parametrization in Eqn3 and table 2 of Gong et al.   ****
    !* (1997)                                                ****
    !************************************************************

    real(rkx), parameter :: c1 = 0.7674_rkx
    real(rkx), parameter :: c2 = 3.079_rkx
    real(rkx), parameter :: c3 = 2.573e-11_rkx
    real(rkx), parameter :: c4 = -1.424_rkx

    !************************************************************
    !*  dry particle density (g/cm3)                          ***
    !************************************************************

    drydens = 2.165_rkx

    !************************************************************
    !* factor for radius (micrometers) to mass (g)            ***
    !************************************************************

    drydens_f = drydens * fourt * mathpi * 1.0e-12_rkx

    !************************************************************
    !*   bubble emissions formula assume RH=80%              ****
    !************************************************************

    relhum = 0.80_rkx

    !************************************************************
    !* rdry star = dry radius (micrometers) below which the *****
    !* dF0/dr emission formula is adjusted downwards        *****
    !************************************************************
    rdry_star = 0.1_rkx

    if ( ireduce_smallr_emit <= 0 ) rdry_star = -1.0e20_rkx

    !************************************************************
    !* sigmag_star  = geometri!standard deviation used for  ****
    !*   rdry < rdry_star                                    ****
    !************************************************************

    sigmag_star = 1.9_rkx

    !************************************************************
    !*                  i n i t i a l i z e    s u m s       ****
    !************************************************************

    dumsum_na = d_zero
    dumsum_ma = d_zero

    !************************************************************
    !* rdrylowermost, rdryuppermost = lower and upper dry   *****
    !* radii (micrometers) for overall integration          *****
    !************************************************************

    rdrylowermost = dpdrylo_cm*0.5e4_rkx
    rdryuppermost = dpdryhi_cm*0.5e4_rkx

    !************************************************************
    !*                    "S E C T I O N  1"               ******
    !************************************************************

    !************************************************************
    !* integrate over rdry > rdry_star, where the dF0/dr    *****
    !* emissions formula is applicable                      *****
    !* (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20, *****
    !* and the entire integration is done here)             *****
    !************************************************************

    if ( rdryuppermost > rdry_star ) then

      !************************************************************
      !*  rdrylo ,  rdryhi  =  lower and upper dry radii     *****
      !*  (micrometers) for this part of the integration      *****
      !************************************************************

      rdrylo = max( rdrylowermost, rdry_star )
      rdryhi = rdryuppermost

      nsub_bin = 1000

      alnrdrylo = log( rdrylo )
      dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

      !************************************************************
      !*  compute rdry, rwet (micrometers) at lowest size       ***
      !************************************************************

      rdrybb = exp( alnrdrylo )
      rdry_cm = rdrybb*1.0e-4_rkx
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onet
      rwetbb = rwet_cm*1.0e4_rkx

      do isub_bin = 1, nsub_bin

        !************************************************************
        !*  rdry, rwet at sub_bin lower boundary are those       ****
        !*  at upper boundary of previous sub_bin                ****
        !************************************************************

        rdryaa = rdrybb
        rwetaa = rwetbb

        !************************************************************
        !* compute rdry, rwet (micrometers) at sub_bin upper    *****
        !* boundary                                             *****
        !************************************************************

        dum = alnrdrylo + isub_bin*dlnrdry
        rdrybb = exp( dum )

        rdry_cm = rdrybb*1.0e-4_rkx
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                 ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onet
        rwetbb = rwet_cm*1.0e4_rkx

        !************************************************************
        !* geometric mean rdry, rwet (micrometers) for sub_bin   ****
        !************************************************************

        rdry = sqrt(rdryaa * rdrybb)
        rwet = sqrt(rwetaa * rwetbb)
        drwet = rwetbb - rwetaa

        !************************************************************
        !*   xmdry is dry mass in g                             *****
        !************************************************************

        xmdry = drydens_f * (rdry**3)

        !************************************************************
        !* dumb is "B" in Gong's Eqn 5a                          ****
        !* df0drwet is "dF0/dr" in Gong's Eqn 5a                 ****
        !************************************************************

        dumb = ( 0.380_rkx - log10(rwet) ) / 0.650_rkx
        dumexpb = exp( -dumb*dumb)
        df0drwet = 1.373_rkx * (rwet**(-d_three)) *       &
                  (1.0_rkx + 0.057_rkx*(rwet**1.05_rkx)) *    &
                  (10.0_rkx**(1.19_rkx*dumexpb))

        dumsum_na = dumsum_na + drwet*df0drwet
        dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
      end do
    end if

    !************************************************************
    !*                "S E C T I O N  2"                   ******
    !************************************************************

    !************************************************************
    !*  integrate over rdry < rdry_star, where the dF0/dr     ***
    !*  emissions formula is just an extrapolation and        ***
    !*  predicts too many emissions                           ***
    !*                                                        ***
    !* 1. compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry) ***
    !* at rdry_star                                           ***
    !* 2.  for rdry < rdry_star, assume dF0/dln(rdry) is      ***
    !* lognormal,with the same lognormal parameters observed  ***
    !*  in O Dowd et al. (1997)                               ***
    !************************************************************

    if ( rdrylowermost < rdry_star ) then

      !************************************************************
      !* compute dF0/dln(rdry) at rdry_star                    ****
      !************************************************************

      rdryaa = 0.99_rkx*rdry_star
      rdry_cm = rdryaa*1.0e-4_rkx
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onet
      rwetaa = rwet_cm*1.0e4_rkx
      rdrybb = 1.01_rkx*rdry_star
      rdry_cm = rdrybb*1.0e-4_rkx
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onet
      rwetbb = rwet_cm*1.0e4_rkx
      rwet = d_half*(rwetaa + rwetbb)
      dumb = ( 0.380_rkx - log10(rwet) ) / 0.650_rkx
      dumexpb = exp( -dumb*dumb)
      df0drwet = 1.373_rkx * (rwet**(-d_three)) *     &
                 (1.0_rkx + 0.057_rkx*(rwet**1.05_rkx)) * &
                 (10.0_rkx**(1.19_rkx*dumexpb))
      drwet = rwetbb - rwetaa
      dlnrdry = log( rdrybb/rdryaa )
      df0dlnrdry_star = df0drwet * (drwet/dlnrdry)

      !************************************************************
      !* rdrylo, rdryhi = lower and upper dry radii (micrometers)**
      !* for this part of the integration                        **
      !************************************************************

      rdrylo = rdrylowermost
      rdryhi = min( rdryuppermost, rdry_star )
      nsub_bin = 1000
      alnrdrylo = log( rdrylo )
      dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

      do isub_bin = 1, nsub_bin

        !************************************************************
        !*    geometric mean rdry (micrometers) for sub_bin      ****
        !************************************************************

        dum = alnrdrylo + (isub_bin-d_half)*dlnrdry
        rdry = exp( dum )

        !************************************************************
        !*  xmdry is dry mass in g                               ****
        !************************************************************

        xmdry = drydens_f * (rdry**3)

        !************************************************************
        !* dumadjust is adjustment factor to reduce dF0/dr      *****
        !************************************************************

        dum = log( rdry/rdry_star ) / log( sigmag_star )
        dumadjust = exp( -d_half*dum*dum )
        df0dlnrdry = df0dlnrdry_star * dumadjust
        dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
        dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
      end do
    end if

    !************************************************************
    !*                     A L L     D O N E                *****
    !************************************************************

    emitfact_numb = dumsum_na
    emitfact_mass = dumsum_ma
    !  print *, emitfact_numb, emitfact_mass

  end subroutine seasalt_emit

end module mod_che_seasalt
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
