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
  use mod_constants

  implicit none

  private

  public :: ocnbats , seaice

  contains

  subroutine ocnbats
    implicit none
    real(rk8) :: ribd , cdrn , rib , qgrd
    real(rk8) :: qs , delq , delt , fact , factuv
    real(rk8) :: cdrmin , cdrx , ribn
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
      qs = qv(i)
      qgrd = pfqsat(tgrd(i),sfps(i))
      delt = tatm(i) - tgrd(i)
      ! Specific humidities
      delq = (qs/(d_one+qs) - qgrd/(d_one+qgrd))
      ! Comnpute drag coefficient over ocean
      ribd = usw(i)**2 + vsw(i)**2 + wtur**2
      cdrn = (vonkar/dlog(ht(i)/zoce))**2
      ribn = ht(i)*egrav*(d_one - tgrd(i)/sts(i))
      rib = ribn/ribd
      if ( rib < d_zero ) then
        cdrx = cdrn*(d_one+24.5D0*dsqrt(-cdrn*rib))
      else
        cdrx = cdrn/(d_one+11.5D0*rib)
      end if
      cdrmin = dmax1(0.25D0*cdrn,6.0D-4)
      if ( cdrx < cdrmin ) cdrx = cdrmin
      drag(i) = cdrx*dsqrt(ribd)*rhox(i)

      ! Update output variables
      evpr(i) = -drag(i)*delq
      sent(i) = -drag(i)*cpd*delt
      if ( dabs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( dabs(evpr(i)) < dlowval ) evpr(i) = d_zero
      fact = dlog(ht(i)*d_half)/dlog(ht(i)/zoce)
      factuv = dlog(ht(i)*d_r10)/dlog(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      t2m(i) = sts(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do
  end subroutine ocnbats

  subroutine seaice
    implicit none
    real(rk8) :: age , scrat , u1 , ribd
    real(rk8) :: cdrn , rib , cdr , qgrd
    real(rk8) :: ps , qs , delq , delt , rhosw , ribl
    real(rk8) :: bb , fact , fss , hrl , hs , hsl
    real(rk8) :: rhosw3 , rsd1 , smc4 , smt , tg , tgrnd , qice
    real(rk8) :: ksnow , rsi , uv995 , sficemm
    real(rk8) :: arg , arg2 , age1 , age2 , tage
    real(rk8) :: cdrmin , cdrx , clead , dela , dela0 , dels
    real(rk8) :: factuv , fevpg , fseng , qgrnd , ribn , sold , vspda
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'ocnbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = iocnbeg , iocnend
      if ( mask(i) /= 2 ) cycle

      ! Update surface temperature from the input SST
      if ( tgb(i) >= icetemp ) then
        tgrd(i) = tgb(i)
      end if

      uv995 = dsqrt(usw(i)**2+vsw(i)**2)

      ! Update Snow Cover
      delt = sts(i) - tgrd(i)
      sold = sncv(i)
      if ( sts(i) < tzero ) then
        ps = prcp(i)
      else
        ps = d_zero
        ! All is melting
        sm(i) = prcp(i)
      end if
      sncv(i) = sncv(i) + dtocn*(ps-evpr(i))

      ! Update Snow age
      if ( sncv(i) < dlowval ) then
        sncv(i) = d_zero
        snag(i) = d_zero
      else
        arg = 5.0D3*(d_one/tzero-d_one/tgrd(i))
        age1 = dexp(arg)
        arg2 = dmin1(d_zero,d_10*arg)
        age2 = dexp(arg2)
        tage = age1 + age2 + age3
        dela0 = 1.0D-6*dtocn
        dela = dela0*tage
        dels = d_r10*dmax1(d_zero,sncv(i)-sold)
        snag(i) = (snag(i)+dela)*(d_one-dels)
        if ( snag(i) < dlowval ) snag(i) = d_zero
      end if
      if ( sncv(i) > 800.0D0 ) snag(i) = d_zero

      ! Scvk is used to compute albedo over seaice
      ! Is the fraction of seaice covered by snow
      age = (d_one-d_one/(d_one+snag(i)))
      scrat = sncv(i)*0.01D0/(d_one+d_three*age)
      scvk(i) = scrat/(0.1D0+scrat)

      ! Compute drag velocity over seaice
      cdrn = (vonkar/dlog(ht(i)/zlnd))**2
      if ( delt < d_zero ) then
        u1 = wtur + d_two*dsqrt(-delt)
      else
        u1 = wtur
      end if
      ribd = usw(i)**2 + vsw(i)**2 + u1**2
      vspda = dsqrt(ribd)
      ribn = ht(i)*egrav*(delt/sts(i))
      rib = ribn/ribd
      if ( rib < d_zero ) then
        cdr = cdrn*(d_one+24.5D0*dsqrt(-cdrn*rib))
      else
        cdr = cdrn/(d_one+11.5D0*rib)
      end if
      cdrmin = dmax1(0.25D0*cdrn,6.0D-4)
      if ( cdr < cdrmin ) cdr = cdrmin
      ! Over snow
      ! rhosw = density of snow relative to water
      rhosw = 0.10D0*(d_one+d_three*age)
      rhosw3 = rhosw**3
      cdrn = (vonkar/dlog(ht(i)/zsno))**2
      ribl = (d_one-271.5D0/sts(i))*ht(i)*egrav/ribd
      if ( ribl < d_zero ) then
        clead = cdrn*(d_one+24.5D0*dsqrt(-cdrn*ribl))
      else
        clead = cdrn/(d_one+11.5D0*rib)
      end if
      cdrx = (d_one-aarea)*cdr + aarea*clead
      drag(i) = cdrx*vspda*rhox(i)

      ! Update now the other variables
      ! The ground temperature and heat fluxes for lake are computed
      ! in the lake model
      qs = qv(i)
      qgrd = pfqsat(tgrd(i),sfps(i))
      ! Move to specific humidity
      qs = qs/(d_one+qs)
      qgrd = qgrd/(d_one+qgrd)
      tgbrd(i) = -d_two + tzero
      ! shice = specific heat of sea-ice per unit volume
      sficemm = sfice(i)*d_1000
      rsd1 = shice*sficemm*d_r1000
      if ( sncv(i) > d_zero ) then
        ! include snow heat capacity
        rsd1 = rsd1 + csnw*sncv(i)*d_r1000
        ! subsurface heat flux through ice
        ! Following Maykut and Untersteiner (1971) and Semtner (1976)
        rsi = 1.4D0*rhosw3*sficemm/sncv(i)
        ksnow = 7.0D-4*rhosw3/sncv(i)
        fss = ksnow * (tgbrd(i)-tgrd(i)) / (d_one + rsi)
      else
        ! Slack, 1980
        fss = 2.14D0*(tgbrd(i)-tgrd(i))/sficemm
      end if
      if ( icpl(i) == 0 ) then
        sfice(i) = (sficemm + 1.087D0*(fss/wlhf)*dtocn) * d_r1000
      end if
      ! set sea ice parameter for melting if seaice less than 2 cm
      if ( sfice(i) <= iceminh ) then
        sncv(i) = d_zero
        scvk(i) = d_zero
        snag(i) = d_zero
        delq = qs - qgrd
        evpr(i) = -drag(i)*delq
        sent(i) = -drag(i)*cpd*delt
        if ( icpl(i) == 0 ) then
          tgrd(i) = tzero
          sfice(i) = d_zero
          mask(i) = 1
        end if
      else
        ! assume lead ocean temp is -1.8c
        ! flux of heat and moisture through leads
        ! sat. mixing ratio at t=-1.8c is 3.3e-3
        qice = 3.3D-3 * stdp/sfps(i)
        !
        qgrnd = ((d_one-aarea)*cdr*qgrd + aarea*clead*qice)/cdrx
        tgrnd = ((d_one-aarea)*cdr*tgrd(i) + aarea*clead*(tzero-1.8D0))/cdrx
        fact = -drag(i)
        delt = sts(i) - tgrnd
        delq = qs - qgrnd
        ! output fluxes, averaged over leads and ice
        evpr(i) = fact*delq
        sent(i) = fact*cpd*delt
        hrl = rhox(i)*vspda*clead * (qice-qs)
        hsl = rhox(i)*vspda*clead * (tzero-1.8D0-sts(i))*cpd
        ! get fluxes over ice for sublimation (subrout snow)
        ! and melt (below) calculation
        fseng = (sent(i)-aarea*hsl)/(d_one-aarea)
        fevpg = (evpr(i)-aarea*hrl)/(d_one-aarea)
        hs = rswf(i) - rlwf(i) - fseng - wlhs*fevpg
        bb = dtocn*(hs+fss)/rsd1
        ! snow melt
        if ( tgrd(i) >= tzero ) sm(i) = sm(i) + (hs+fss)/wlhf
        if ( sm(i) <= d_zero ) sm(i) = d_zero
        smc4 = sm(i)*dtocn
        sncv(i) = sncv(i) - sm(i) * dtocn
        ! all snow removed, melt ice
        if ( sncv(i) < smc4 ) then
          smt = (sncv(i)/dtocn)
          ! rho(h2o)/rho(ice) = 1.087
          if ( icpl(i) == 0 ) then
            sfice(i) = sfice(i) + dtocn*(smt-sm(i))*1.087D0*d_r1000
          end if
          sm(i) = sm(i)+smt
          ! set sea ice parameter for melting if less than 2 cm
          if ( sfice(i) <= iceminh ) then
            sncv(i) = d_zero
            scvk(i) = d_zero
            snag(i) = d_zero
            if ( icpl(i) == 0 ) then
              tgrd(i) = tzero
              sfice(i) = d_zero
              mask(i) = 1
            end if
          end if
        else
          ! snow or ice with no snow melting
          tg = tgrd(i) + bb
          if ( tg >= tzero ) tgrd(i) = tzero
          if ( tg < tzero ) tgrd(i) = tg
          tgb(i) = tgrd(i)
        end if
      end if
      if ( dabs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( dabs(evpr(i)) < dlowval ) evpr(i) = d_zero
      fact = dlog(ht(i)*d_half)/dlog(ht(i)/zoce)
      factuv = dlog(ht(i)*d_r10)/dlog(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      taux(i) = dmissval
      tauy(i) = dmissval
      t2m(i) = sts(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do
  end subroutine seaice

end module mod_ocn_bats
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
