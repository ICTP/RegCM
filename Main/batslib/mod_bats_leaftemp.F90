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

module mod_bats_leaftemp
  !
  ! Calculate leaf temperature, leaf fluxes, and net transpiration.
  ! documented in NCAR Tech Note, Dickinson et al., 1986.
  ! modifications by Klaus Blumel, 1988.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_bats_internal
  use mod_bats_param
  use mod_service

  implicit none

  private

  public :: lftemp , satur , fseas

  contains

!=======================================================================
!     based on: bats version 1e          copyright 18 august 1989
!=======================================================================
!
!     Calculate leaf temperature, leaf fluxes, and net transpiration.
!     documented in NCAR Tech Note, Dickinson et al., 1986.
!     modifications by Klaus Blumel, 1988.
!
!        f l o w   d i a g r a m   f o r   l e f t e m
!
!              lftemp ===> stomat
!                          frawat
!                            root
!                           satur
!                          lfdrag
!                          condch
!                          condcq
!                           deriv
!
!     cf = heat transfer coeff. from leaves; assumes same for moisture
!          (see e.g. d. gates' book) for laminar flow past leaf;
!          from ewing paper, cf has dimensions t**1/2 l**-1;
!          so to make it dimensionless, premultiply by 0.01 (si).
!     cgrnd = deriv. of soil energy flux with respect to soil temp.
!                             (used in tgrund)
!     delt = ts-taf
!     ef = transpiration rate
!     efpot = potential evaporation rate (kg/m**2/s)
!     fevpg = evaporative heat flux from ground
!     flnet = (temp gradient)*(d/dt(sig t**4)); hence, t**3 term
!     fseng = sensible heat flux from ground
!     iter = leaf temp iteration counter; runs from 1 to max of 100
!     ra = leaf aerodynamic resistance factor
!     taf = air temperature within foliage canopy
!     tbef = leaf temp at time step before current one
!     ts = air temperature of lowest model layer
!     uaf = mean wind within canopy
!
!     taf(i) = temperature of air in canopy
!     delt(i)= difference between temperature of overlying air
!                          and that in canopy
!     delq(i)= difference between humidity of overlying air
!                          and that in canopy
!
!     convergence of leaf temperature calculation is declared if
!     enough iterations (itmin) and change of temp small enough and
!     change of latent heat fluxes small enough, or if
!     maximum iteration reached (itmax).
!
!=======================================================================

  subroutine lftemp
    implicit none
    real(rkx) :: dcn , delmax , efeb , eg1 , epss , fbare , qbare , &
               qcan , qsatdg , rppdry , sf1 , sf2 , sgtg3 , vakb ,  &
               xxkb , efpot , tbef , dels
    integer(ik4) :: iter , itfull , itmax , i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lftemp'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !=======================================================================
    ! 1.   setup information
    !=======================================================================
    !
    ! 1.1  get stress-free stomatal resistance
    !      (1st guess at vapor pressure deficit)
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        vpdc(i) = d_10
        sgtg3 = emiss(i)*(sigm*tgrd(i)**3)
        flneto(i) = d_four*sgtg3*(tlef(i)-tgrd(i))
      else
        flneto(i) = d_zero
      end if
    end do

    call stomat
    !
    ! 1.3  determine fraction of total and green canopy surface
    !      covered by water
    call frawat
    !
    ! 1.4  establish root function in terms of etrc = maximum
    !      sustainable transpiration rate
    ! (routine also returns efpr, used in subr. water to define upper soil
    !  layer transpiration)
    !
    call root
    ! 1.5  saturation specific humidity of leaf
    call satur(qsatl,tlef,sfcp)

    !==========================================================
    !l    2.   begin iteration for leaf temperature calculation
    !==========================================================
    iter = 0
    efeb = d_zero
    delmax = d_one
    itmax = 10
    itfull = itmax
    ! itmax = 40
    ! itfull = 20

    do iter = 0 , itmax
      !
      ! 2.1  recalc stability dependent canopy & leaf drag coeffs
      !
      if ( iter == 0 ) call condch

      call lfdrag
      call condch

      do i = ilndbeg , ilndend
        if ( sigf(i) > minsigf ) then
          lftra(i) = d_one/(cf(i)*uaf(i))
          cn1(i) = wtlh(i)*rhs(i)
          df(i) = cn1(i)*cpd
          !
          ! 2.2  decrease foliage conductance for stomatal
          !      resistance
          rppdry = lftra(i)*fdry(i)/(lftrs(i)+lftra(i))
          rpp(i) = rppdry + fwet(i)
          !
          ! 2.3  recalculate saturation vapor pressure
          !
          eg1 = eg(i)
          eg(i) = pfesat(tlef(i))
          ! call bats_psat(tlef(i),eg(i))
          qsatl(i) = qsatl(i)*eg(i)/eg1
        end if
      end do

      ! 2.4  canopy evapotranspiration
      if ( iter == 0 ) call condcq

      epss = 1.0e-10_rkx
      do i = ilndbeg , ilndend
        if ( sigf(i) > minsigf ) then
          efpot = cn1(i)*(wtgaq(i)*qsatl(i) - wtgq0(i)*qgrd(i) - wtaq0(i)*qs(i))
          if ( efpot > d_zero ) then
            etr(i) = efpot*lftra(i)*fdry(i)/(lftrs(i)+lftra(i))
            rpp(i) = min(rpp(i),(etr(i)+ldew(i)/dtbat)/efpot-epss)
          else
            etr(i) = d_zero
            rpp(i) = d_one
          end if
          if ( ( efpot >= d_zero ) .and. ( etr(i) >= etrc(i) ) )  then
            ! transpiration demand exceeds supply, stomat adjust
            ! demand
            rppdry = lftra(i)*fdry(i)/(lftrs(i)+lftra(i))
            rppdry = rppdry/(etr(i)/etrc(i))
            etr(i) = etrc(i)
            ! recalculate stomatl resistance and rpp
            lftrs(i) = lftra(i)*(fdry(i)/rppdry-d_one)
            rpp(i) = rppdry + fwet(i)
            rpp(i) = min(rpp(i),(etr(i)+ldew(i)/dtbat)/efpot-epss)
          end if
          rppq(i) = htvp(i)*rpp(i)
          efe(i) = rppq(i)*efpot
          if ( efe(i)*efeb < d_zero ) efe(i) = 0.1_rkx*efe(i)
        else
          etr(i) = etrc(i)
        end if
      end do
      !=====================================
      !      3.   solve for leaf temperature
      !=====================================
      !  3.1  update conductances for changes in rpp and cdr
      call condcq
      !
      !  3.2  derivatives of energy fluxes with respect to leaf
      !       temperature for newton-raphson calculation of leaf temperature.
      !  subr.  ii: rs,ra,cdrd,rppq,efe.
      !  subr. output: qsatld,dcd.
      if ( iter <= itfull ) call deriv
      !
      !  3.3  compute dcn from dcd, output from subr. deriv
      do i = ilndbeg , ilndend
        if ( sigf(i) > minsigf ) then
          dcn = dcd(i)*tlef(i)
          ! 1.2  radiative forcing for leaf temperature calculation
          sgtg3 = emiss(i)*(sigm*tgrd(i)**3)
          sf1 = sigf(i)*(abswveg(i)-lwflx(i)-(d_one-sigf(i))* &
                flneto(i)+d_four*sgtg3*tgrd(i))
          sf2 = d_four*sigf(i)*sgtg3 + df(i)*wtga(i) + dcd(i)
          ! 3.4  iterative leaf temperature calculation
          tbef = tlef(i)
          tlef(i) = (sf1+df(i)*(wta0(i)*sts(i)+wtg0(i)*tgrd(i))-efe(i)+dcn)/sf2
          !
          ! 3.5  chk magnitude of change; limit to max allowed value
          dels = tlef(i) - tbef
          if ( abs(dels) > delmax ) then
            tlef(i) = tbef + delmax*dels/abs(dels)
          end if
          ! 3.6  update dependence of stomatal resistance
          !      on vapor pressure deficit
          qcan = wtlq0(i)*qsatl(i) + qgrd(i)*wtgq0(i) + qs(i)*wtaq0(i)
          vpdc(i) = (d_one-rpp(i))*(qsatl(i)-qcan)*d_1000/ep2
        end if
      end do
      call stomat
      ! 3.8  end iteration
    end do

    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        !=========================================
        ! 4.   update dew accumulation (kg/m**2/s)
        !=========================================
        ldew(i) = ldew(i) + (etr(i) - efe(i)/htvp(i))*dtbat
        !===========================================
        ! 5.   collect parameters needed to evaluate
        !      sensible and latent fluxes
        !===========================================
        !  5.1  canopy properties
        taf(i) = wtg0(i)*tgrd(i) + wta0(i)*sts(i) + wtl0(i)*tlef(i)
        delt(i) = wtgl(i)*sts(i) - (wtl0(i)*tlef(i) + wtg0(i)*tgrd(i))
        delq(i) = wtglq(i)*qs(i) - (wtlq0(i)*qsatl(i) + wtgq0(i)*qgrd(i))
        sgtg3 = emiss(i)*(sigm*tgrd(i)**3)
        flnet(i) = sgtg3*(tlef(i)-tgrd(i))*d_four
        xxkb = min(rough(lveg(i)),d_one)
        vakb = (d_one-sigf(i))*vspda(i) + sigf(i) * &
               (xxkb*uaf(i)+(d_one-xxkb)*vspda(i))
        wtg2(i) = (d_one-sigf(i))*cdr(i)*vakb
        fbare = wtg2(i)*(tgrd(i)-sts(i))
        qbare = wtg2(i)*(qgrd(i)-qs(i))
        !  5.2  fluxes from soil
        fseng(i) = cpd*rhs(i)*(wtg(i)*((wta0(i) +    &
                   wtl0(i))*tgrd(i)-wta0(i)*sts(i) - &
                   wtl0(i)*tlef(i))+fbare)
        fevpg(i) = rhs(i)*rgr(i) * (wtg(i)*((wtaq0(i) +  &
                   wtlq0(i))*qgrd(i)-wtaq0(i)*qs(i) -    &
                   wtlq0(i)*qsatl(i))+qbare)
        !  5.3  deriv of soil energy flux with respect to soil temp
        qsatdg = pfqsdt(tgrd(i),sfcp(i))
        qsatdg = qsatdg * rgr(i)
        ! call bats_qsdt(tgrd(i),qgrd(i),qsatdg)
        ! qsatdg = qsatdg * rgr(i)
        cgrnds(i) = rhs(i)*cpd*(wtg(i)*(wta0(i)+wtl0(i))+wtg2(i))
        cgrndl(i) = rhs(i)*qsatdg*((wta(i)+wtlq(i)) * &
                    wtg(i)*wtsqi(i)+wtg2(i))
        cgrnd(i) = cgrnds(i) + cgrndl(i)*htvp(i)
        !  5.4  reinitialize cdrx
        !     shuttleworth mods #3 removed here !!!!!!
        cdrx(i) = cdr(i)
        !
        !  5.5  fluxes from canopy and soil to overlying air
        fbare = wtg2(i)*(tgrd(i)-sts(i))
        qbare = wtg2(i)*(qgrd(i)-qs(i))
        sent(i) = cpd*rhs(i)*(-wta(i)*delt(i)+fbare)
        evpr(i) = rhs(i)*(-wta(i)*delq(i) + rgr(i)*qbare)
        if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
        if ( abs(evpr(i)) < dlowval ) evpr(i) = d_zero
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>
#include <pfdesatdt.inc>
#include <pqderiv.inc>

  end subroutine lftemp

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     gives leaf stomatal resistance from environmental parameters
!             under conditions of no moisture stress
!
!     standard lai from xla=max & xlai0=min lai
!     fseas is a seasonal factor for reduced winter lai
!           and root water uptake
!        fc = light sensitivity for crops and grasses and has inverse
!             radiation units (m**2/watt)
!      rlai = sum of leaf and stem area indices
!     rmax0 = 20000. s/m (maximum resistance, in mod_constants)
!     radu & radl = visible light intensity in upper & lower canopy
!     ft & fb = the fractional intercepted photo-active radiation
!             per unit (leaf & stem) area in the top (upper) and
!             bottom (lower) canopies, respectively
!     radfi = average of upper and lower canopy light factors
!        rs = stomatal resistance = min.res. * rad.factor * leaf factor
!      trup = transmission of the upper canopy, assumed to be the same
!             for the lower canopy,i.e., trup=exp(-0.5*g*rlai/czenith),
!             where g = attenuation factor
!
!     documented in NCAR Tech Note, Dickinson et al., 1986
!     improved stomatal shading, Dickinson, nov 88.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine stomat
    implicit none
    real(rkx) :: difzen , g , radf , radfi , seas , vpdf , rilmax , &
                 rmini , fsol0 , fsold , maxarg
    real(rkx) :: trup , trupd
    integer(ik4) :: il , ilmax , i
    real(rkx) , dimension(10) :: rad , radd
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'stomat'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! seasonal temperature factor
    ! g is average leaf crosssection per unit lai
    ! difzen is ave of inverse of cos of angle of diffuse vis light
    ! ilmax is number of canopy layers
    ! czenith is cosine solar zenith angle for incident light
    ! (to spec from input data need a good treatment of diffuse rad)
    ! trup is transmission of direct beam light in one canopy layer
    ! trupd is transmission of diffuse light in one canopy layer
    !
    g = d_half
    difzen = d_two
    ilmax = 4
    rilmax = d_four
    maxarg = -log(epsilon(d_one))
    call fseas(tlef,bseas)
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        if ( czenith(i) > 1.e-3_rkx ) then
          fsold = fracd(i)*swsi(i)*fc(lveg(i))
          fsol0 = (d_one-fracd(i))*swsi(i)*fc(lveg(i))
          if ( g*rlai(i)/(rilmax*czenith(i)) > maxarg ) then
            rad(1) = fsol0*rilmax/rlai(i)
            do il = 2 , ilmax
              rad(il) = dlowval
            end do
          else
            trup = exp(-g*rlai(i)/(rilmax*czenith(i)))
            rad(1)  = (d_one-trup) *fsol0*rilmax/rlai(i)
            do il = 2 , ilmax
              rad(il)  = trup* rad(il-1)
            end do
          end if
          if ( difzen*g*rlai(i)/rilmax > maxarg ) then
            radd(1) = fsold*rilmax/rlai(i)
            do il = 2 , ilmax
              radd(il) = dlowval
            end do
          else
            trupd = exp(-difzen*g*rlai(i)/rilmax)
            radd(1) = (d_one-trupd)*fsold*rilmax/rlai(i)
            do il = 2 , ilmax
              radd(il) = trupd*radd(il-1)
            end do
          end if
          rmini = rsmin(lveg(i))/rmax0
          radfi = d_zero
          do il = 1 , ilmax
            radfi = radfi+(rmini+rad(il)+radd(il))/(d_one+rad(il)+radd(il))
          end do
          radf = rilmax/radfi
          vpdf = d_one/max(0.3_rkx,d_one-vpdc(i)*0.025_rkx)
          seas = d_one/(rmini+bseas(i))
          lftrs(i) = rsmin(lveg(i))*radf*seas*vpdf
          lftrs(i) = min(lftrs(i),rmax0)
        else
          lftrs(i) = rsmin(lveg(i))
        end if
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine stomat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  **** determines fraction of foliage covered by water fwet, and
!  **** the fraction of foliage that is dry transpiring leaf fdry.
!  note: their defns differ - fwet is the fraction of all veg surfaces
!  which are wet because mod_stems can evaporate, fdry is the fraction
!  of lai which is dry because only leaves can transpire
!
!  ldew(i) is in kg/m**2/s
!  fwet   = ratio of dew to max value to 2/3 power
!           ( 2/3 power comes from deardorff (1978) )
!              ** keep fwet le 1.0 **
!  dewmaxi = inverse of max allowed dew depth on leaf in mm
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine frawat
    implicit none
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'frawat'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        fwet(i) = d_zero
        if ( ldew(i) > d_zero ) then
          fwet(i) = ((dewmaxi/vegt(i))*ldew(i))**twot
          fwet(i) = min(fwet(i),d_one)
        end if
        fdry(i) = (d_one-fwet(i))*xlai(i)/xlsai(i)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine frawat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     this subroutine provides root function in terms of maximum
!     transpiration rate plants can sustain depending on soil moisture.
!
!     trsmx0 is a prescribed constant (kg/m**2/s).
!     trsmx is the maximum transpiration rate,
!        including a low temperature correction (=seasb)
!        and a correction for fractional vegetation (=sigf).
!
!     rotf is ratio of moisture extracton from top to total when
!           fully saturated
!     rootf is ratio of roots in upper soil layer
!                    to roots in root soil layer
!     bsw is the b param in clapp and hornberger
!
!     "wlt  " are ratios factors controlling the saturation
!                 cf wilting (see ewing paper)
!     wlttb (total) & wltub (upper) become 1 at the wilting point
!     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
!
!     etrc= max poss transpiration given the soil moisture distributions
!     efpr = the relative contribution of upper soil layer to
!     evapotranspiration - need soil moist. budget (subrout water)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine root
    implicit none
    real(rkx) :: bneg , rotf , trsmx , wlttb , wltub , wmli
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'root'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        ! trsmx = trsmx0*sigf(i)*seasb(i)
        trsmx = trsmx0*sigf(i)
        rotf = rootf(lveg(i))
        bneg = -bsw(i)
        wmli = d_one/(wiltr(i)**bneg-d_one)
        wlttb = (watr(i)**bneg-d_one)*wmli
        wltub = (watu(i)**bneg-d_one)*wmli
        wlttb = min(wlttb,d_one)
        wltub = min(wltub,d_one)
        etrc(i) = trsmx*(d_one-(d_one-rotf)*wlttb-rotf*wltub)
        efpr(i) = trsmx*rotf*(d_one-wltub)
        if ( etrc(i) < 1.0e-12_rkx ) then
          etrc(i) = 1.0e-12_rkx
          efpr(i) = d_one
        else
          efpr(i) = efpr(i)/etrc(i)
        end if
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine root

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****  calculates saturation vapor pressure (eg)
!           and qsat = saturated specific humidity (dimensionless)
!
!           uses tetens formula (1930) (ref. riegel,1974,jam,p606
!                                                 equation 1)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine satur(qsat,t,p)
    implicit none
    real(rkx) , pointer , dimension(:) , intent(in) :: p , t
    real(rkx) , pointer , dimension(:) , intent(inout) :: qsat
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'satur'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ilndbeg , ilndend
      eg(i) = pfesat(t(i))
      qsat(i) = pfqsat(t(i),p(i),eg(i))
      ! call bats_satur(t(i),p(i),eg(i),qsat(i))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>

  end subroutine satur

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     recalculate stability dependent drag coefficient for vegetation,
!     given the neutral drag coefficient.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine lfdrag
    implicit none
    real(rkx) :: dthdz , ribi , sqrtf , tkb , u1 , u2 , zatild , cdrmin
    real(rkx) :: dlstaf , rib1
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lfdrag'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        tkb = wta0(i)*sts(i) + wtl0(i)*tlef(i) + wtg0(i)*tgrd(i)
        dlstaf = sts(i) - sigf(i)*tkb - (d_one-sigf(i))*tgrd(i)
        if ( dlstaf <= d_zero ) then
          dthdz = (d_one-sigf(i))*tgrd(i) + sigf(i)*tkb-sts(i)
          if ( dthdz > dlowval ) then
            u1 = wtur + d_two*sqrt(dthdz)
          else
            u1 = wtur
          end if
          ribd(i) = usw(i)**2 + vsw(i)**2 + u1**2
        else
          u2 = wtur
          ribd(i) = usw(i)**2 + vsw(i)**2 + u2**2
        end if
        vspda(i) = sqrt(ribd(i))
        if ( vspda(i) < d_one ) then
          vspda(i) = d_one
          ribd(i) = d_one
        end if
        zatild = (zh(i)-displa(lveg(i)))*sigf(i) + zh(i)*(d_one-sigf(i))
        rib1 = egrav*zatild/(ribd(i)*sts(i))
        rib(i) = rib1*dlstaf
        if ( rib(i) < d_zero ) then
          cdr(i) = cdrn(i)*(d_one+24.5_rkx * sqrt(-cdrn(i)*rib(i)))
          sqrtf = min(sqrt(-cdrn(i)/rib(i)),11.5_rkx/12.25_rkx)
          cdrd(i) = cdrn(i)*12.25_rkx*wtl0(i)*rib1*sigf(i)*sqrtf
        else
          ribi = d_one/(d_one+11.5_rkx*rib(i))
          cdr(i) = cdrn(i)*ribi
          cdrd(i) = cdr(i)*ribi*11.5_rkx*rib1*wtl0(i)*sigf(i)
          cdrmin = max(cdrn(i)*d_rfour,6.0e-4_rkx)
          if ( (cdr(i) < cdrmin) ) then
            cdr(i) = cdrmin
            cdrd(i) = d_zero
          end if
        end if
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine lfdrag

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional sensible heat conductances
!               for canopy and soil flux calculations
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine condch
    implicit none
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'condch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !     csoilc = constant drag coefficient for soil under canopy
    !     symbols used for weights are:   wt : weight
    !     a : air
    !     l : leaf
    !     i : inverse
    !     s : sum
    !     h : sensible heat
    !     q : water vapor
    !     0 : normalized (sums to one)
    !     g : ground
    !
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        uaf(i) = vspda(i)*sqrt(cdr(i))
        cf(i) = 0.01_rkx*sqrtdi(lveg(i))/sqrt(uaf(i))
        wta(i) = sigf(i)*cdr(i)*vspda(i)
        wtlh(i) = cf(i)*uaf(i)*vegt(i)
        wtg(i) = csoilc*uaf(i)*sigf(i)
        wtshi(i) = d_one/(wta(i)+wtlh(i)+wtg(i))
        wtl0(i) = wtlh(i)*wtshi(i)
        wtg0(i) = wtg(i)*wtshi(i)
        wtgl(i) = wtl0(i) + wtg0(i)
        wta0(i) = d_one - wtgl(i)
        wtga(i) = wta0(i) + wtg0(i)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine condch

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional latent heat conductances
!               for canopy and soil flux calculations
!
!     latent fluxes differ from sensible due to stomatal resistance
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine condcq
    implicit none
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'condcq'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !     symbols used for weights are:   wt : weight
    !     a : air
    !     l : leaf
    !     i : inverse
    !     s : sum
    !     h : sensible heat
    !     q : water vapor
    !     0 : normalized (sums to one)
    !     g : ground
    !
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        rgr(i) = gwet(i)
        wtlq(i) = wtlh(i)*rpp(i)
        wtgq(i) = wtg(i)*rgr(i)
        wtsqi(i) = d_one/(wta(i)+wtlq(i)+wtgq(i))
        wtgq0(i) = wtgq(i)*wtsqi(i)
        wtlq0(i) = wtlq(i)*wtsqi(i)
        wtglq(i) = wtgq0(i) + wtlq0(i)
        wtaq0(i) = d_one - wtglq(i)
        wtgaq(i) = wtaq0(i) + wtgq0(i)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine condcq

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     derivatives of energy fluxes with respect to leaf temperature for
!     newton-raphson calculation of leaf temperature.
!     input: rs,ra,cdrd,rppq,efe.    output: dcd.
!
!     approximate by derivatives of cdr and ef.  many weaker
!     dependences on leaf temperature are omitted, as convergence
!     rate is not affected.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine deriv
    implicit none
    real(rkx) :: hfl , xkb , qsatld
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'deriv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        qsatld = pfqsdt(tlef(i),sfcp(i))
        ! call bats_qsdt(tlef(i),qsatl(i),qsatld)
        xkb = cdrd(i)/cdr(i)
        hfl = df(i)*(wtga(i)*tlef(i) - wtg0(i)*tgrd(i) - wta0(i)*sts(i))
        dcd(i) = cn1(i)*rppq(i)*wtgaq(i) * qsatld + (d_one-wtgaq(i)) *   &
                 efe(i) * xkb + (d_one-wtga(i)) * hfl * xkb
        dcd(i) = max(dcd(i),d_zero)
        dcd(i) = min(dcd(i),500.0_rkx)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>
#include <pfdesatdt.inc>
#include <pqderiv.inc>

  end subroutine deriv

  subroutine fseas(temp,ffsea)
    ! The seasonal function is a number between 0 and 1
    ! If the temperature is greater than 298.0, it is 1
    ! If the temperature is less than 298.0, it is less than 1
    ! If the temperature is less than 273.0, it is zero
    implicit none
    real(rkx) , pointer , dimension(:) , intent(in) :: temp
    real(rkx) , pointer , dimension(:) , intent(inout) :: ffsea
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'fseas'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( lcrop_cutoff ) then
      do i = ilndbeg , ilndend
        if ( lveg(i) == 1 ) then
          ffsea(i) = max(d_zero,d_one - 0.0016_rkx * &
                     (max(298.0_rkx-temp(i),d_zero)**4))
        else
          ffsea(i) = max(d_zero,d_one - 0.0016_rkx * &
                     (max(298.0_rkx-temp(i),d_zero)**2))
        end if
      end do
    else
      do i = ilndbeg , ilndend
        ffsea(i) = max(d_zero,d_one - 0.0016_rkx * &
                   (max(298.0_rkx-temp(i),d_zero)**2))
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine fseas

end module mod_bats_leaftemp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
