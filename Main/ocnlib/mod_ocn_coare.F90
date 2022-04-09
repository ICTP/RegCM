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
!
module mod_ocn_coare
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_service
  use mod_ocn_internal
  use mod_runparams , only : iocnrough

  implicit none

  private

  public :: coare3_drv

  contains

    !
    !-----------------------------------------------------------------------
    ! This routine computes the bulk parameterization of surface
    ! atmospheric variables, wind stress and surface net heat fluxes
    !
    ! Adapted from COARE version 3.0 code for snow and ice which is
    ! written originally by C. Fairall (04/2005)
    !-----------------------------------------------------------------------
    !
    subroutine coare3_drv()
      implicit none
      real(rkx) :: ts , qs , us , uv995 , t995 , q995 , z995 , ta
      real(rkx) :: zu , zt , zq , zi , du , dt , dq , ut , dter
      real(rkx) :: ug , zogs , u10 , cdhg , zo10 , zot10
      real(rkx) :: usr , qsr , tsr , zetu , l10 , wetc , zet
      real(rkx) :: cd10 , ch10 , ct10 , cc , cd , ct , ribcu
      real(rkx) :: rr , rt , rq , zo , zot , zoq , dels , bigc , Al
      real(rkx) :: l , Bf , tkt , qout , qcol , alq , xlamx , dqer
      real(rkx) :: le , visa , rhoa , cpv , Rns , Rnl
      real(rkx) :: hsb , hlb , tau , uv10 , facttq
      integer(ik4) :: i , k , niter
      logical :: iflag

      real(rkx) , parameter :: beta = 1.25_rkx   ! gustiness coeff.
      real(rkx) , parameter :: fdg  = 1.0_rkx    ! ratio of thermal to wind VonKarman
      real(rkx) , parameter :: visw = 1.0e-6_rkx ! water kinematic viscosity
      real(rkx) , parameter :: tcw  = 0.6_rkx    ! water thermal diffusivity
      real(rkx) , parameter :: rhow = 1022.0     ! water density
      real(rkx) , parameter :: be   = 0.026_rkx  ! sal. expans. coef. of water
      real(rkx) , parameter :: cpw  = 4.0e3_rkx  ! spec. heat of water

      do i = iocnbeg , iocnend
        if ( mask(i) /= 1 ) cycle

        tgrd(i) = tgb(i)
        tgbrd(i) = tgb(i)

        !
        !-----------------------------------------
        ! Input bulk parameterization fields
        !-----------------------------------------
        !
        iflag = .false.
        if ( lseaice .or. llake ) then
          if (sfice(i) > 0.0_rkx) iflag = .true.
        end if
        ts = tgrd(i) - tzero
        ! comes from coupled model
        us = 0.0_rkx ! current speed
        uv995 = sqrt(usw(i)**2+vsw(i)**2)
        t995 = tatm(i)-tzero
        q995 = qv(i)
        z995 = ht(i)
        ta = sfta(i)
        ! height of the atmospheric data
        zu = z995
        zt = z995
        zq = z995

        ! height (m) of atmospheric boundary layer
        zi = hpbl(i)
        !
        !---------------------------------------
        ! Air constants
        !---------------------------------------
        !
        ! specific heat of moist air
        cpv = cpmf(q995)

        ! latent heat of vaporization (J/kg) at sea surface
        le = wlh(tgrd(i))

        ! moist air density (kg/m3)
        rhoa = sfps(i)/(rgas*ta*(d_one+ep1*q995))

        ! kinematic viscosity of dry air (m2/s), Andreas (1989)
        visa = 1.326e-5_rkx*(d_one+6.542e-3_rkx*t995+8.301e-6_rkx*t995*t995- &
               4.84e-9_rkx*t995*t995*t995)

        bigc = 16.0*egrav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)

        ! water thermal expansion coefft
        if (ts > -2.0_rkx) then
          Al = 2.1e-5_rkx*(ts+3.2_rkx)**0.79_rkx
        else
          Al = 2.4253e-05_rkx
        end if
        !
        !-----------------------------------------------------
        ! Compute net longwave and shortwave radiation (W/m2)
        !-----------------------------------------------------
        !
        Rns = rswf(i)
        Rnl = rlwf(i)
        !
        !--------------------------------
        ! Begin bulk loop - first guess
        !--------------------------------
        !
        qs = pfqsat(tgrd(i),sfps(i))*0.98_rkx
        wetc = pfqsdt(tgrd(i),sfps(i))

        ! Move to specific humidities
        q995 = q995/(d_one+q995)

        ! Deltas
        dt = ta - t995 - tzero
        dq = qs - q995
        du = uv995 - us

        ! assume that wind is measured relative to sea surface
        ! and include gustiness
        ug = d_half
        ut = sqrt(du*du+ug*ug)

        if (iflag) then
          dter = d_zero
        else
          dter = 0.3_rkx
        end if
        !
        !-----------------------
        ! Neutral coefficients
        !-----------------------
        !
        zogs = 1e-4_rkx
        if ( iflag ) then
          zogs = 4.5e-4_rkx
        end if

        u10 = ut*log(10.0_rkx/zogs)/log(zu/zogs)
        cdhg = vonkar/log(10.0_rkx/zogs)
        usr = cdhg*u10

        ! initial guess for the friction velocity
        if ( iflag ) then
          zo10 = zogs
        else
          zo10 = 0.011_rkx*usr*usr*regrav+0.11_rkx*visa/usr
        end if

        ! initial guess for drag coefficents
        cd10 = (vonkar/log(10.0_rkx/zo10))**2
        if ( iflag ) then
          ch10 = 0.0015_rkx
        else
          ch10 = 0.00115_rkx
        end if
        ct10 = ch10/sqrt(cd10)
        zot10 = 10.0_rkx/exp(vonkar/ct10)
        cd = (vonkar/log(zu/zo10))**2
        !
        !----------------------------
        ! Compute Richardson number
        !----------------------------
        !
        ct = vonkar/log(zt/zot10)
        cc = vonkar*ct/cd
        ribcu = -zu/zi/0.004_rkx/beta**3
        br(i) = -egrav*zu/ta*((dt-dter)+ep1*ta*dq)/ut**2

        niter = 3
        if ( br(i) < 0.0_rkx ) then
          zetu = cc*br(i)/(1.0_rkx+br(i)/ribcu) ! unstable
        else
          zetu = cc*br(i)*(1.0_rkx+3.0_rkx*br(i)/cc) ! stable
        end if
        l10 = zu/zetu
        if (zetu > 50.0) niter = 1
        !
        !---------------------------------------------------
        ! First guesses for Monon-Obukhov similarity scales
        !---------------------------------------------------
        !
        usr = ut*vonkar/(log(zu/zo10)-psiuo(zu/l10))
        tsr = -(dt-dter)*vonkar*fdg/(log(zt/zot10)-psit(zt/l10))
        qsr = -(dq-wetc*dter)*vonkar*fdg/(log(zq/zot10)-psit(zq/l10))
        tkt = 0.001_rkx ! cool skin thickness
        !
        !------------------------------------------
        ! Main loop - limited with three iteration
        !------------------------------------------
        !
        do k = 1, niter
          zet = vonkar*egrav*zu/ta*(tsr+ep1*ta*qsr)/(usr*usr)
          if ( iflag ) then
            zo = zogs
          else
            zo = 0.011_rkx*usr*usr*regrav+0.11_rkx*visa/usr
          end if

          rr = zo*usr/visa

          ! for snow/ice, Andreas, 1987
          if ( iflag ) then
            if ( rr <= 0.135_rkx ) then
              rt = rr*exp(1.250_rkx)
              rq = rr*exp(1.610_rkx)
            else if ( rr <= 2.5_rkx ) then
              rt = rr*exp(0.149_rkx-0.550_rkx*log(rr))
              rq = rr*exp(0.351_rkx-0.628_rkx*log(rr))
            else if ( rr <= 1000_rkx ) then
              rt = rr*exp(0.317_rkx-0.565_rkx*log(rr)-0.183_rkx*log(rr)*log(rr))
              rq = rr*exp(0.396_rkx-0.512_rkx*log(rr)-0.180_rkx*log(rr)*log(rr))
            else
              rt = 1e-10_rkx
              rq = 1e-10_rkx
            end if
            ! for ocean, Lui et al., 1979
          else
            if ( rr <= 0.11_rkx ) then
              rt = 0.177_rkx
              rq = 0.292_rkx
            else if ( rr <= 0.8_rkx ) then
              rt = 1.376_rkx*rr**0.929_rkx
              rq = 1.808_rkx*rr**0.826_rkx
            else if ( rr <= 3.0_rkx ) then
              rt = 1.026_rkx*rr**(-0.599_rkx)
              rq = 1.393_rkx*rr**(-0.528_rkx)
            else if ( rr <= 10.0_rkx ) then
              rt = 1.625_rkx*rr**(-1.018_rkx)
              rq = 1.956_rkx*rr**(-0.870_rkx)
            else if ( rr <= 30.0_rkx ) then
              rt = 4.661_rkx*rr**(-1.475_rkx)
              rq = 4.994_rkx*rr**(-1.297_rkx)
            else if (rr <= 100.0_rkx) then
              rt = 34.904_rkx*rr**(-2.067_rkx)
              rq = 30.709_rkx*rr**(-1.845_rkx)
            else if ( rr <= 300.0_rkx ) then
              rt = 1667.19_rkx*rr**(-2.907_rkx)
              rq = 1448.68_rkx*rr**(-2.682_rkx)
            else if ( rr <= 1000.0_rkx ) then
              rt = 5.88e5_rkx*rr**(-3.935_rkx)
              rq = 2.98e5_rkx*rr**(-3.616_rkx)
            else
              rt = 1e-10_rkx
              rq = 1e-10_rkx
            end if
          end if

          ! compute Monin-Obukhov stability parameter, Z/L
          L = zu/zet
          zot = rt*visa/usr
          zoq = rq*visa/usr

          ! update scaling params
          ram1(i) = (log(zu/zo)-psiuo(zu/L))
          rah1(i) = (log(zt/zot)-psit(zt/L))
          usr = ut*vonkar/ram1(i)
          tsr = -(dt-dter)*vonkar*fdg/rah1(i)
          qsr = -(dq-wetc*dter)*vonkar*fdg/(log(zq/zoq)-psit(zq/L))

          ! compute gustiness in wind speed
          Bf = -egrav/ta*usr*(tsr+ep1*ta*qsr)
          if ( Bf > 0.0_rkx ) then
            ug = beta*(Bf*zi)**0.333_rkx
          else
            ug = 0.2_rkx
          end if
          ut = sqrt(du*du+ug*ug)

          ! background sensible and latent heat flux
          hsb = -rhoa*cpd*usr*tsr
          hlb = -rhoa*Le*usr*qsr

          ! total cooling at the interface
          qout = Rnl+hsb+hlb
          dels = Rns*(0.137_rkx+11.0_rkx*tkt-6.6e-5_rkx / &
                   tkt*(1.0_rkx-exp(-tkt/8.0e-4_rkx)))
          qcol = qout-dels

          if ( qcol > 0 ) then
            ! buoy flux water, Eq. 7
            alq = Al*qcol+be*hlb*cpw/Le
            ! Saunders, Eq. 13
            if ( alq < 0.0_rkx ) then
              dter = d_zero
            else
              xlamx = 6.0_rkx/(1.0_rkx+(bigc*alq/usr**4)**0.75_rkx)**0.333_rkx
              ! sub. thk, Eq. 11
              tkt = xlamx*visw/(sqrt(rhoa/rhow)*usr)
              ! cool skin, Eq. 12
              dter = qcol*tkt/tcw
            end if
          else
            dter = d_zero
          end if

          if ( iflag ) then
            dter = d_zero
          end if

          dqer = wetc*dter
        end do
        !
        if ( zetu < d_zero ) then
          uv10 = uv995+(usr/vonkar)*(log(10.0_rkx/zu)- &
                       (psiuo(zetu)-psiuo(zu/l10)))
        else
          uv10 = uv995+(usr/vonkar)*(log(10.0_rkx/zu)+ &
                       d_five*zetu-d_five*zu/l10)
        end if
        !
        !-------------------------------------
        ! Compute atmosphere/ocean/ice fluxes
        !-------------------------------------
        !
        ! heat fluxes
        sent(i) = hsb
        evpr(i) = max(hlb/Le,d_zero)
        if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
        if ( evpr(i) < dlowval ) evpr(i) = d_zero

        ! drag coefficents
        facttq = log(z995*d_half)/log(z995/zo)
        ustr(i) = usr
        zoo(i) = zo
        drag(i) = usr**2*rhox(i)/uv995

        ! wind stress components
        tau = rhoa*usr*usr*du/ut
        taux(i) = tau*(usw(i)/uv995)
        tauy(i) = tau*(vsw(i)/uv995)

        ! wind components
        u10m(i) = usw(i)*uv10/uv995
        v10m(i) = vsw(i)*uv10/uv995

        ! surface atmospheric variables
        t2m(i) = t995+tzero-dt*facttq
        q2m(i) = q995-dq*facttq
      end do

      contains

#include <pfesat.inc>
#include <pfqsat.inc>
#include <pfdesatdt.inc>
#include <pqderiv.inc>
#include <wlh.inc>
#include <cpmf.inc>

      pure real(rkx) function psiuo(zet)
        implicit none
        real(rkx) , intent (in) :: zet
        real(rkx) :: x, psik, f, psic, c
        if (zet < 0.0_rkx) then
          x = (1.0_rkx-15.0_rkx*zet)**0.25_rkx
          psik = 2.0_rkx*log((1.0_rkx+x)/2.0_rkx)+log((1.0_rkx+x*x)/2.0_rkx)- &
                 2.0_rkx*atan(x)+2.0_rkx*atan(1.0_rkx)
          x = (1.0_rkx-10.15_rkx*zet)**0.3333_rkx
          psic = 1.5_rkx*log((1.0_rkx+x+x*x)/3.0_rkx)-sqrt(3.0_rkx)*        &
                 atan((1.0_rkx+2.0_rkx*x)/sqrt(3.0_rkx))+4.0_rkx*           &
                 atan(1.0_rkx)/sqrt(3.0_rkx)
          f = zet*zet/(1.0_rkx+zet*zet)
          psiuo = (1.0_rkx-f)*psik+f*psic
        else
          c = min(50.0_rkx,0.35_rkx*zet)
          psiuo = -((1.0_rkx+1.0_rkx*zet)**1.0_rkx+0.667_rkx*                 &
                  (zet-14.28_rkx)/exp(c)+8.525_rkx)
        end if
      end function psiuo

      pure real(rkx) function psit(zet)
        implicit none
        real(rkx) , intent (in) :: zet
        real(rkx) :: x, psik, f, psic, c
        if (zet < 0.0_rkx) then
          x = (1.0_rkx-15.0_rkx*zet)**0.5_rkx
          psik = 2.0_rkx*log((1.0_rkx+x)/2.0_rkx)
          x = (1.0_rkx-34.15_rkx*zet)**0.3333_rkx
          psic = 1.5_rkx*log((1.0_rkx+x+x*x)/3.0_rkx)-sqrt(3.0_rkx)* &
                 atan((1.0_rkx+2.0_rkx*x)/sqrt(3.0_rkx))+4.0_rkx*    &
                 atan(1.0_rkx)/sqrt(3.0_rkx)
          f = zet*zet/(1.0_rkx+zet*zet)
          psit = (1.0_rkx-f)*psik+f*psic
        else
          c = min(50.0_rkx,0.35_rkx*zet)
          psit = -((1.0_rkx+2.0_rkx/3.0_rkx*zet)**1.5_rkx+0.667_rkx*     &
                  (zet-14.28_rkx)/exp(c)+8.525_rkx)
        end if
      end function psit

    end subroutine coare3_drv

end module mod_ocn_coare
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
