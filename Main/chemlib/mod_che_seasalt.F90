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

module mod_che_seasalt

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common

  implicit none

  private

  ! sea-salt density
  real(rk8) , parameter :: rhosslt = 1000.0D0
  real(rk8) , dimension(sbin,2) :: ssltbsiz

  data ssltbsiz /0.05D0, 1.0D0, 1.0D0, 10.0D0/
!
  real(rk8) , dimension(sbin) :: ssltbed


 ! solubility of od dust aer for param of giorgi and chameides
  real(rk8) , dimension(sbin) :: solsslt

  data ssltbed /0.6D0, 6.0D0/
  data solsslt /0.8D0, 0.8D0/

  public :: sea_salt , rhosslt , ssltbsiz , solsslt , ssltbed

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

  subroutine sea_salt(j,wind10,ivegcov,seasalt_flx)

    implicit none

    integer(ik4) , intent(in) :: j
    integer(ik4) , dimension(ici1:ici2) :: ivegcov
    real(rk8) , dimension(ici1:ici2) :: wind10
    real(rk8) , dimension(ici1:ici2,sbin) , intent(out) :: seasalt_flx

    real(rk8) :: dumu10
    real(rk8) :: dplo_acc , dphi_acc , dplo_cor , dphi_cor
    real(rk8) , dimension(ici1:ici2,2) :: qflxm , qflxn
    real(rk8) :: seasalt_emfac_mascor , seasalt_emfac_numcor
    real(rk8) :: seasalt_emfac_masacc , seasalt_emfac_numacc
    integer(ik4) :: modeptr_coarseas
    integer(ik4) :: modeptr_accum
    integer(ik4) :: i , iflg , ib

    real(rk8) :: specmw_seasalt_amode
    real(rk8) :: dum_mw

    specmw_seasalt_amode = 1.0D0
    modeptr_coarseas = 1
    modeptr_accum = 1

    dplo_acc = ssltbsiz(1,1)*1.0D-04
    dphi_acc = ssltbsiz(1,2)*1.0D-04
    dplo_cor = ssltbsiz(2,1)*1.0D-04
    dphi_cor = ssltbsiz(2,2)*1.0D-04

    !************************************************************
    !* initialization -- compute seasalt emissions factors    ***
    !* for accumulation and coarse modes                      ***
    !************************************************************

    seasalt_flx(:,1) = d_zero
    seasalt_flx(:,2) = d_zero

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

    seasalt_emfac_numacc = seasalt_emfac_numacc * 1.0D-3
    seasalt_emfac_numcor = seasalt_emfac_numcor * 1.0D-3

    !************************************************************
    !* convert mass factors from (g/m2/s) to (kmol/m2/s)     ****
    !************************************************************

    dum_mw = specmw_seasalt_amode
    seasalt_emfac_masacc = seasalt_emfac_masacc * (1.0D-3/dum_mw)
    seasalt_emfac_mascor = seasalt_emfac_mascor * (1.0D-3/dum_mw)

    ! print*, seasalt_emfac_masacc, seasalt_emfac_mascor

    do i = ici1 , ici2
      if ( ivegcov(i) == 0 ) then
        iflg = 1
      else
        iflg = 0
      end if

      seasalt_flx(i,1) = seasalt_emfac_masacc * iflg
      seasalt_flx(i,2) = seasalt_emfac_mascor * iflg

      if ( wind10(i) > d_zero ) then
        dumu10 = dmax1( d_zero, dmin1( 100.0D0, wind10(i) ))
        dumu10 = dumu10**3.41D0
        qflxm(i,1) = seasalt_emfac_mascor * dumu10
        qflxn(i,1) = seasalt_emfac_numcor * dumu10
        qflxm(i,2) = seasalt_emfac_masacc * dumu10
        qflxn(i,2) = seasalt_emfac_numacc * dumu10

        seasalt_flx(i,1) = seasalt_flx(i,1) * dumu10
        seasalt_flx(i,2) = seasalt_flx(i,2) * dumu10
      else
        seasalt_flx(i,1) = d_zero
        seasalt_flx(i,2) = d_zero
      end if
    end do

    ! Update the emission tendency
    do ib = 1 , sbin
      ! calculate the source tendancy
      do i = ici1 , ici2
        ! chemsrc(i,j,lmonth,isslt(ib)) = seasalt_flx(i,ib)

        chiten(j,i,kz,isslt(ib)) = chiten(j,i,kz,isslt(ib)) + &
                seasalt_flx(i,ib)*egrav/(dsigma(kz)*1.D3)
        ! diagnostic source
        cemtrac(j,i,isslt(ib)) = cemtrac(j,i,isslt(ib)) + &
                 seasalt_flx(i,ib)*cfdout

      end do
    end do

  end subroutine sea_salt
!
!===================================================================
!
  subroutine seasalt_emit(ireduce_smallr_emit,dpdrylo_cm, &
                          dpdryhi_cm,emitfact_numb,emitfact_mass)
    implicit none

    integer(ik4) , intent(in) :: ireduce_smallr_emit
    real(rk8) , intent(in) :: dpdrylo_cm , dpdryhi_cm
    real(rk8) , intent(out) :: emitfact_numb , emitfact_mass

    integer(ik4) :: nsub_bin , isub_bin
    real(rk8) :: drydens , drydens_f
    real(rk8) :: relhum
    real(rk8) :: rdry_star , sigmag_star
    real(rk8) :: dumsum_na , dumsum_ma
    real(rk8) :: rdrylowermost , rdryuppermost
    real(rk8) :: rdrylo , rdryhi
    real(rk8) :: alnrdrylo , dlnrdry
    real(rk8) :: rdrybb , rwetbb
    real(rk8) :: rdryaa , rwetaa
    real(rk8) :: rdry_cm , rwet_cm
    real(rk8) :: rdry , rwet , drwet
    real(rk8) :: dum,xmdry , dumb,dumexpb
    real(rk8) :: dumadjust , df0dlnrdry
    real(rk8) :: df0drwet
    real(rk8) :: df0dlnrdry_star

    !************************************************************
    !*  c1-c4 are constants for sea-salt hygroscopi!growth  ****
    !*  parametrization in Eqn3 and table 2 of Gong et al.   ****
    !* (1997)                                                ****
    !************************************************************

    real(rk8) , parameter :: c1 = 0.7674D0
    real(rk8) , parameter :: c2 = 3.079D0
    real(rk8) , parameter :: c3 = 2.573D-11
    real(rk8) , parameter :: c4 = -1.424D0

    !************************************************************
    !*  dry particle density (g/cm3)                          ***
    !************************************************************

    drydens = 2.165D0

    !************************************************************
    !* factor for radius (micrometers) to mass (g)            ***
    !************************************************************

    drydens_f = drydens * fourt * mathpi * 1.0D-12

    !************************************************************
    !*   bubble emissions formula assume RH=80%              ****
    !************************************************************

    relhum = 0.80D0

    !************************************************************
    !* rdry star = dry radius (micrometers) below which the *****
    !* dF0/dr emission formula is adjusted downwards        *****
    !************************************************************
    rdry_star = 0.1D0

    if ( ireduce_smallr_emit <= 0 ) rdry_star = -1.0D20

    !************************************************************
    !* sigmag_star  = geometri!standard deviation used for  ****
    !*   rdry < rdry_star                                    ****
    !************************************************************

    sigmag_star = 1.9D0

    !************************************************************
    !*                  i n i t i a l i z e    s u m s       ****
    !************************************************************

    dumsum_na = d_zero
    dumsum_ma = d_zero

    !************************************************************
    !* rdrylowermost, rdryuppermost = lower and upper dry   *****
    !* radii (micrometers) for overall integration          *****
    !************************************************************

    rdrylowermost = dpdrylo_cm*0.5D4
    rdryuppermost = dpdryhi_cm*0.5D4

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
      !*  rdrylo  ,  rdryhi  =  lower and upper dry radii     *****
      !*  (micrometers) for this part of the integration      *****
      !************************************************************

      rdrylo = max( rdrylowermost, rdry_star )
      rdryhi = rdryuppermost

      nsub_bin = 1000

      alnrdrylo = dlog( rdrylo )
      dlnrdry = (dlog( rdryhi ) - alnrdrylo)/nsub_bin

      !************************************************************
      !*  compute rdry, rwet (micrometers) at lowest size       ***
      !************************************************************

      rdrybb = dexp( alnrdrylo )
      rdry_cm = rdrybb*1.0D-4
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - dlog10(relhum) ) )**onet
      rwetbb = rwet_cm*1.0D4

      do isub_bin = 1 , nsub_bin

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
        rdrybb = dexp( dum )

        rdry_cm = rdrybb*1.0D-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                 ( (c3*(rdry_cm**c4)) - dlog10(relhum) ) )**onet
        rwetbb = rwet_cm*1.0D4

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

        dumb = ( 0.380D0 - dlog10(rwet) ) / 0.650D0
        dumexpb = dexp( -dumb*dumb)
        df0drwet = 1.373D0 * (rwet**(-d_three)) *       &
                  (1.0D0 + 0.057D0*(rwet**1.05D0)) *    &
                  (10.0D0**(1.19D0*dumexpb))

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

      rdryaa = 0.99D0*rdry_star
      rdry_cm = rdryaa*1.0D-4
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - dlog10(relhum) ) )**onet
      rwetaa = rwet_cm*1.0D4
      rdrybb = 1.01D0*rdry_star
      rdry_cm = rdrybb*1.0D-4
      rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2)) / &
                ( (c3*(rdry_cm**c4)) - dlog10(relhum) ) )**onet
      rwetbb = rwet_cm*1.0D4
      rwet = d_half*(rwetaa + rwetbb)
      dumb = ( 0.380D0 - dlog10(rwet) ) / 0.650D0
      dumexpb = dexp( -dumb*dumb)
      df0drwet = 1.373D0 * (rwet**(-d_three)) *     &
                 (1.0D0 + 0.057D0*(rwet**1.05D0)) * &
                 (10.0D0**(1.19D0*dumexpb))
      drwet = rwetbb - rwetaa
      dlnrdry = dlog( rdrybb/rdryaa )
      df0dlnrdry_star = df0drwet * (drwet/dlnrdry)

      !************************************************************
      !* rdrylo, rdryhi = lower and upper dry radii (micrometers)**
      !* for this part of the integration                        **
      !************************************************************

      rdrylo = rdrylowermost
      rdryhi = min( rdryuppermost, rdry_star )
      nsub_bin = 1000
      alnrdrylo = dlog( rdrylo )
      dlnrdry = (dlog( rdryhi ) - alnrdrylo)/nsub_bin

      do isub_bin = 1 , nsub_bin

        !************************************************************
        !*    geometric mean rdry (micrometers) for sub_bin      ****
        !************************************************************

        dum = alnrdrylo + (isub_bin-d_half)*dlnrdry
        rdry = dexp( dum )

        !************************************************************
        !*  xmdry is dry mass in g                               ****
        !************************************************************

        xmdry = drydens_f * (rdry**3)

        !************************************************************
        !* dumadjust is adjustment factor to reduce dF0/dr      *****
        !************************************************************

        dum = dlog( rdry/rdry_star ) / dlog( sigmag_star )
        dumadjust = dexp( -d_half*dum*dum )
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
    !  print *, emitfact_numb , emitfact_mass

  end subroutine seasalt_emit

end module mod_che_seasalt
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
