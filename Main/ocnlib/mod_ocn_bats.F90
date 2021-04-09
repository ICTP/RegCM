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

module mod_ocn_bats
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_ocn_internal
  use mod_runparams , only : icetriggert
  use mod_constants

  implicit none

  private

  public :: ocnbats , seaice

  contains

  subroutine ocnbats
    implicit none
    real(rkx) :: ribd , cdrn , qgrd
    real(rkx) :: qs , delq , delt , fact , factuv
    real(rkx) :: cdrmin , cdrx , ribn , vspda
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'ocnbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = iocnbeg , iocnend
      if ( mask(i) /= 1 ) cycle

      ! Update surface temperature from the input SST
      tgrd(i) = tgb(i)
      tgbrd(i) = tgb(i)

      ! Compute delt and delq
      qs = qv(i)/(d_one+qv(i))
      qgrd = pfqsat(tgrd(i),sfps(i))
      delt = sfta(i) - tgrd(i)
      ! Specific humidities
      delq = qs - qgrd
      ! Comnpute drag coefficient over ocean
      ribd = usw(i)**2 + vsw(i)**2 + wtur**2
      vspda = sqrt(ribd)
      cdrn = (vonkar/log(ht(i)/zoce))**2
      ribn = ht(i)*egrav*(d_one - tgrd(i)/tatm(i))
      br(i) = ribn/ribd
      if ( br(i) < d_zero ) then
        cdrx = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*br(i)))
      else
        cdrx = cdrn/(d_one+11.5_rkx*br(i))
      end if
      cdrmin = max(0.25_rkx*cdrn,6.0e-4_rkx)
      if ( cdrx < cdrmin ) cdrx = cdrmin
      ram1(i) = d_one/cdrx
      rah1(i) = d_one/cdrx
      drag(i) = cdrx*sqrt(ribd)*rhox(i)
      ustr(i) = sqrt((vspda*drag(i))/rhox(i))
      zoo(i) = 0.01_rkx*regrav*ustr(i)*ustr(i)

      ! Update output variables
      evpr(i) = max(-drag(i)*delq,d_zero)
      sent(i) = -drag(i)*cpd*delt
      if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( evpr(i) < dlowval ) evpr(i) = d_zero
      fact = log(ht(i)*d_half)/log(ht(i)/zoce)
      factuv = log(ht(i)*d_r10)/log(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      t2m(i) = tatm(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>

  end subroutine ocnbats

  subroutine seaice
    implicit none
    real(rkx) :: age , u1 , ribd
    real(rkx) :: cdrn , cdr , qgrd
    real(rkx) :: ps , qs , delq , delt , rhosw , ribl
    real(rkx) :: bb , fact , fss , hrl , hs , hsl
    real(rkx) :: rhosw3 , rsd1 , smc4 , smt , tg , tgrnd , qice
    real(rkx) :: ksnow , rsi , uv995 , sficemm , tau , xdens
    real(rkx) :: arg , arg2 , age1 , age2
    real(rkx) :: cdrmin , cdrx , clead , dela , dela0 , dels
    real(rkx) :: factuv , fevpg , fseng , qgrnd , ribn , sold , vspda
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'seaice'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = iocnbeg , iocnend
      if ( mask(i) /= 2 ) cycle

      ! Update surface temperature from the input SST
      if ( tgb(i) >= icetriggert ) then
        tgrd(i) = tgb(i)
      else
        if ( tatm(i) > icetriggert ) then
          tgrd(i) = icetriggert
        else
          tgrd(i) = tatm(i) - 0.01_rkx
        end if
      end if
      tgbrd(i) = icetriggert

      uv995 = max(sqrt(usw(i)**2+vsw(i)**2),wtur)

      ! Update Snow Cover
      delt = tatm(i) - tgrd(i)
      sold = sncv(i)
      if ( tatm(i) < tzero ) then
        ps = prcp(i)
      else
        ps = d_zero
      end if
      sncv(i) = sncv(i) + dtocn*ps

      ! Update Snow age
      if ( sncv(i) < dlowval ) then
        sncv(i) = d_zero
        snag(i) = d_zero
      else
        arg = 5.0e3_rkx*(d_one/tzero-d_one/tgrd(i))
        age1 = exp(arg)
        arg2 = min(d_zero,d_10*arg)
        age2 = exp(arg2)
        age = age1 + age2 + age3
        dela0 = 1.0e-6_rkx*dtocn
        dela = dela0*age
        dels = d_r10*max(d_zero,sncv(i)-sold)
        snag(i) = (snag(i)+dela)*(d_one-dels)
        if ( snag(i) < dlowval ) snag(i) = d_zero
      end if
      if ( sncv(i) > 800.0_rkx ) snag(i) = d_zero

      age = (d_one-d_one/(d_one+snag(i)))
      ! Compute drag velocity over seaice
      cdrn = (vonkar/log(ht(i)/zlnd))**2
      if ( delt < d_zero ) then
        u1 = wtur + d_two*sqrt(-delt)
      else
        u1 = wtur
      end if
      ribd = usw(i)**2 + vsw(i)**2 + u1**2
      vspda = sqrt(ribd)
      ribn = ht(i)*egrav*(delt/tatm(i))
      br(i) = ribn/ribd
      if ( br(i) < d_zero ) then
        cdr = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*br(i)))
      else
        cdr = cdrn/(d_one+11.5_rkx*br(i))
      end if
      cdrmin = max(0.25_rkx*cdrn,6.0e-4_rkx)
      if ( cdr < cdrmin ) cdr = cdrmin
      ! Over snow
      ! rhosw = density of snow relative to water
      rhosw = 0.10_rkx*(d_one+d_three*age)
      rhosw3 = rhosw**3
      cdrn = (vonkar/log(ht(i)/zsno))**2
      ribl = (d_one-icetriggert/tatm(i))*ht(i)*egrav/ribd
      if ( ribl < d_zero ) then
        clead = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*ribl))
      else
        clead = cdrn/(d_one+11.5_rkx*br(i))
      end if
      cdrx = (d_one-aarea)*cdr + aarea*clead
      ram1(i) = d_one/cdrx
      rah1(i) = d_one/cdrx
      drag(i) = cdrx*vspda*rhox(i)
      ustr(i) = sqrt((vspda*drag(i))/rhox(i))
      zoo(i) = 0.01_rkx*regrav*ustr(i)*ustr(i)

      ! Update now the other variables
      ! The ground temperature and heat fluxes for lake are computed
      ! in the lake model
      qs = qv(i)/(d_one+qv(i))
      ! shice = specific heat of sea-ice per unit volume
      sficemm = sfice(i)*d_1000
      rsd1 = shice*sficemm*d_r1000
      if ( sncv(i) > d_zero ) then
        ! include snow heat capacity
        rsd1 = rsd1 + csnw*sncv(i)*d_r1000
        ! subsurface heat flux through ice
        ! Following Maykut and Untersteiner (1971) and Semtner (1976)
        rsi = 1.4_rkx*rhosw3*sficemm/sncv(i)
        ksnow = 7.0e-4_rkx*rhosw3/sncv(i)
        fss = ksnow * (tgbrd(i)-tgrd(i)) / (d_one + rsi)
      else
        ! Slack, 1980
        fss = 2.14_rkx*(tgbrd(i)-tgrd(i))/sficemm
      end if
      if ( icpl(i) == 0 ) then
        sfice(i) = (sficemm + 1.087_rkx*(fss/wlhf)*dtocn) * d_r1000
      end if
      ! set sea ice parameter for melting if seaice less than 2 cm
      if ( sfice(i) <= iceminh ) then
        qgrd = pfqsat(271.36_rkx,sfps(i))
        sncv(i) = d_zero
        snag(i) = d_zero
        if ( icpl(i) == 0 ) then
          tgrd(i) = tgb(i)
          tgbrd(i) = tgb(i)
          sfice(i) = d_zero
          mask(i) = 1
        end if
        delq = qs - qgrd
        evpr(i) = max(-drag(i)*delq,d_zero)
        sent(i) = -drag(i)*cpd*delt
      else
        ! assume lead ocean temp is icetriggert
        ! flux of heat and moisture through leads
        fact = -drag(i)
        qgrd = pfqsat(fact*tatm(i),sfps(i))
        qice = pfqsat(icetriggert,sfps(i))
        !
        qgrnd = ((d_one-aarea)*cdr*qgrd + aarea*clead*qice)/cdrx
        tgrnd = ((d_one-aarea)*cdr*tgrd(i) + aarea*clead*icetriggert)/cdrx
        delt = tatm(i) - tgrnd
        delq = qs - qgrnd
        ! output fluxes, averaged over leads and ice
        evpr(i) = max(fact*delq,d_zero)
        sncv(i) = sncv(i) - dtocn*evpr(i)
        sent(i) = fact*cpd*delt
        hrl = rhox(i)*vspda*clead * (qice-qs)
        hsl = rhox(i)*vspda*clead * (icetriggert-tatm(i))*cpd
        ! get fluxes over ice for sublimation (subrout snow)
        ! and melt (below) calculation
        fseng = (sent(i)-aarea*hsl)/(d_one-aarea)
        fevpg = (evpr(i)-aarea*hrl)/(d_one-aarea)
        hs = rswf(i) - rlwf(i) - fseng - wlhs*fevpg
        bb = dtocn*(hs+fss)/rsd1
        ! snow melt
        if ( tgrd(i) >= tzero ) sm(i) = sm(i) + (hs+fss)/wlhf
        if ( sm(i) <= d_zero ) sm(i) = d_zero
        if ( sm(i) * dtocn > sncv(i) ) sm(i) = sncv(i)/dtocn
        smc4 = sm(i)*dtocn
        sncv(i) = sncv(i) - sm(i) * dtocn
        ! all snow removed, melt ice
        if ( sncv(i) < smc4 ) then
          smt = (sncv(i)/dtocn)
          ! rho(h2o)/rho(ice) = 1.087
          if ( icpl(i) == 0 ) then
            sfice(i) = sfice(i) + dtocn*(smt-sm(i))*1.087_rkx*d_r1000
          end if
          sm(i) = sm(i)+smt
          ! set sea ice parameter for melting if less than 2 cm
          if ( sfice(i) <= iceminh ) then
            sncv(i) = d_zero
            snag(i) = d_zero
            if ( icpl(i) == 0 ) then
              tgrd(i) = tgb(i)
              tgbrd(i) = tgb(i)
              sfice(i) = d_zero
              mask(i) = 1
            end if
          end if
        else
          ! snow or ice with no snow melting
          tg = tgrd(i) + bb
          tg = min(tg,icetriggert)
          tgrd(i) = tg
          tgb(i) = tg
        end if
      end if
      if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
      fact = log(ht(i)*d_half)/log(ht(i)/zoce)
      factuv = log(ht(i)*d_r10)/log(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      rhoa(i) = rhox(i)
      xdens = sfps(i)/(rgas*tgrd(i)*(d_one+ep1*qv(i)))
      tau = xdens*ustr(i)*ustr(i)
      taux(i) = tau*(usw(i)/uv995)
      tauy(i) = tau*(vsw(i)/uv995)
      t2m(i) = tatm(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>

  end subroutine seaice

end module mod_ocn_bats
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
