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
      real(rk8) :: ts , qs , us , uv995 , t995 , q995 , z995 , ta
      real(rk8) :: zu , zt , zq , zi , du , dt , dq , ut , dter
      real(rk8) :: ug , zogs , u10 , cdhg , zo10 , zot10
      real(rk8) :: usr , qsr , tsr , zetu , l10 , wetc , zet
      real(rk8) :: cd10 , ch10 , ct10 , cc , cd , ct , ribcu , ribu
      real(rk8) :: rr , rt , rq , zo , zot , zoq , dels , bigc , Al
      real(rk8) :: L , Bf , tkt , qout , qcol , alq , xlamx , dqer
      real(rk8) :: Le , visa , rhoa , cpv , Rns , Rnl
      real(rk8) :: hsb , hlb , tau , uv10 , facttq
      integer(ik4) :: i , k , niter
      logical :: iflag

      real(rk8) , parameter :: beta = 1.25D0 ! gustiness coeff.
      real(rk8) , parameter :: fdg  = 1.0D0 ! ratio of thermal to wind VonKarman
      real(rk8) , parameter :: visw = 1.0d-6 ! water kinematic viscosity
      real(rk8) , parameter :: tcw  = 0.6D0     ! water thermal diffusivity
      real(rk8) , parameter :: rhow = 1022.0    ! water density
      real(rk8) , parameter :: be   = 0.026D0   ! sal. expans. coef. of water
      real(rk8) , parameter :: cpw  = 4.0d3     ! spec. heat of water
      real(rk8) , parameter :: cpa  = 1004.67D0 ! spec. heat of moist. air


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
          if (sfice(i) > 0.0D0) iflag = .true.
        end if
        ts = tgrd(i) - tzero
        ! comes from coupled model
        us = 0.0D0 ! current speed
        uv995 = dsqrt(usw(i)**2+vsw(i)**2)
        t995 = sts(i)-tzero
        q995 = qv(i)
        z995 = ht(i)
        ta = t995 + tzero
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
        cpv = cpa*(1.0D0+0.84D0*q995)

        ! latent heat of vaporization (J/kg) at sea surface
        Le = wlh(tgrd(i))

        ! moist air density (kg/m3)
        rhoa = sfps(i)/(rgas*ta*(d_one+0.61D0*q995))

        ! kinematic viscosity of dry air (m2/s), Andreas (1989)
        visa = 1.326d-5*(d_one+6.542d-3*t995+8.301d-6*t995*t995- &
               4.84d-9*t995*t995*t995)

        bigc = 16.0*egrav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)

        ! water thermal expansion coefft
        if (ts > -2.0D0) then
          Al = 2.1d-5*(ts+3.2D0)**0.79D0
        else
          Al = 2.4253d-05
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
        qs = pfqsat(tgrd(i),sfps(i))*0.98D0
        wetc = pfqsdt(tgrd(i),sfps(i))

        ! Move all to specific humidities
        qs = qs/(d_one+qs)
        q995 = q995/(d_one+q995)
        wetc = wetc/(d_one+wetc)

        ! Deltas
        dt = ts - t995 - 0.0098D0*zt
        dq = qs - q995
        du = uv995 - us

        ! assume that wind is measured relative to sea surface
        ! and include gustiness
        ug = d_half
        ut = dsqrt(du*du+ug*ug)

        if (iflag) then
          dter = d_zero
        else
          dter = 0.3D0
        end if
        !
        !-----------------------
        ! Neutral coefficients
        !-----------------------
        !
        zogs = 1D-4
        if ( iflag ) then
          zogs = 4.5D-4
        end if

        u10 = ut*dlog(10.0D0/zogs)/dlog(zu/zogs)
        cdhg = vonkar/dlog(10.0D0/zogs)
        usr = cdhg*u10

        ! initial guess for the friction velocity
        if ( iflag ) then
          zo10 = zogs
        else
          zo10 = 0.011D0*usr*usr*regrav+0.11D0*visa/usr
        end if

        ! initial guess for drag coefficents
        cd10 = (vonkar/dlog(10.0D0/zo10))**2
        if ( iflag ) then
          ch10 = 0.0015D0
        else
          ch10 = 0.00115D0
        end if
        ct10 = ch10/dsqrt(cd10)
        zot10 = 10.0D0/dexp(vonkar/ct10)
        cd = (vonkar/dlog(zu/zo10))**2
        !
        !----------------------------
        ! Compute Richardson number
        !----------------------------
        !
        ct = vonkar/dlog(zt/zot10)
        cc = vonkar*ct/cd
        ribcu = -zu/zi/0.004D0/beta**3
        ribu = -egrav*zu/ta*((dt-dter)+0.61D0*ta*dq)/ut**2

        niter = 3
        if (ribu < 0.0D0) then
          zetu = cc*ribu/(1.0D0+ribu/ribcu) ! unstable
        else
          zetu = cc*ribu*(1.0D0+3.0D0*ribu/cc) ! stable
        end if
        l10 = zu/zetu
        if (zetu > 50.0) niter = 1
        !
        !---------------------------------------------------
        ! First guesses for Monon-Obukhov similarity scales
        !---------------------------------------------------
        !
        usr = ut*vonkar/(dlog(zu/zo10)-psiuo(zu/l10))
        tsr = -(dt-dter)*vonkar*fdg/(dlog(zt/zot10)-psit(zt/l10))
        qsr = -(dq-wetc*dter)*vonkar*fdg/(dlog(zq/zot10)-psit(zq/l10))
        tkt = 0.001D0 ! cool skin thickness
        !
        !------------------------------------------
        ! Main loop - limited with three iteration
        !------------------------------------------
        !
        do k = 1, niter
          zet = vonkar*egrav*zu/ta*(tsr+0.61D0*ta*qsr)/(usr*usr)
          if ( iflag ) then
            zo = zogs
          else
            zo = 0.011D0*usr*usr*regrav+0.11D0*visa/usr
          end if

          rr = zo*usr/visa

          ! for snow/ice, Andreas, 1987
          rt = d_one
          rq = d_one
          if ( iflag ) then
            if ( rr <= 0.135D0 ) then
              rt = rr*dexp(1.250D0)
              rq = rr*dexp(1.610D0)
            else if ( rr <= 2.5D0 ) then
              rt = rr*dexp(0.149D0-0.550D0*dlog(rr))
              rq = rr*dexp(0.351D0-0.628D0*dlog(rr))
            else if ( rr <= 1000D0 ) then
              rt = rr*dexp(0.317D0-0.565D0*dlog(rr)-0.183D0*dlog(rr)*dlog(rr))
              rq = rr*dexp(0.396D0-0.512D0*dlog(rr)-0.180D0*dlog(rr)*dlog(rr))
            end if
            ! for ocean, Lui et al., 1979
          else
            if ( rr <= 0.11D0 ) then
              rt = 0.177D0
              rq = 0.292D0
            else if ( rr <= 0.8D0 ) then
              rt = 1.376D0*rr**0.929D0
              rq = 1.808D0*rr**0.826D0
            else if ( rr <= 3.0D0 ) then
              rt = 1.026D0*rr**(-0.599D0)
              rq = 1.393D0*rr**(-0.528D0)
            else if ( rr <= 10.0D0 ) then
              rt = 1.625D0*rr**(-1.018D0)
              rq = 1.956D0*rr**(-0.870D0)
            else if ( rr <= 30.0D0 ) then
              rt = 4.661D0*rr**(-1.475D0)
              rq = 4.994D0*rr**(-1.297D0)
            else if (rr <= 100.0D0) then
              rt = 34.904D0*rr**(-2.067D0)
              rq = 30.709D0*rr**(-1.845D0)
            else if ( rr <= 300.0D0 ) then
              rt = 1667.19D0*rr**(-2.907D0)
              rq = 1448.68D0*rr**(-2.682D0)
            else if ( rr <= 1000.0D0 ) then
              rt = 5.88d5*rr**(-3.935D0)
              rq = 2.98d5*rr**(-3.616D0)
            end if
          end if

          ! compute Monin-Obukhov stability parameter, Z/L
          L = zu/zet
          zot = rt*visa/usr
          zoq = rq*visa/usr

          ! update scaling params
          usr = ut*vonkar/(dlog(zu/zo)-psiuo(zu/L))
          tsr = -(dt-dter)*vonkar*fdg/(dlog(zt/zot)-psit(zt/L))
          qsr = -(dq-wetc*dter)*vonkar*fdg/(dlog(zq/zoq)-psit(zq/L))

          ! compute gustiness in wind speed
          Bf = -egrav/ta*usr*(tsr+0.61D0*ta*qsr)
          if ( Bf > 0.0D0 ) then
            ug = beta*(Bf*zi)**0.333D0
          else
            ug = 0.2D0
          end if
          ut = dsqrt(du*du+ug*ug)

          ! background sensible and latent heat flux
          hsb = -rhoa*cpa*usr*tsr
          hlb = -rhoa*Le*usr*qsr

          ! total cooling at the interface
          qout = Rnl+hsb+hlb
          dels = Rns*(0.137D0+11.0D0*tkt-6.6d-5/tkt*(1.0D0-dexp(-tkt/8.0d-4)))
          qcol = qout-dels

          if ( qcol > 0 ) then
            ! buoy flux water, Eq. 7
            alq = Al*qcol+be*hlb*cpw/Le
            ! Saunders, Eq. 13
            if ( alq < 0.0D0 ) then
              dter = d_zero
            else
              xlamx = 6.0D0/(1.0D0+(bigc*alq/usr**4)**0.75D0)**0.333D0
              ! sub. thk, Eq. 11
              tkt = xlamx*visw/(dsqrt(rhoa/rhow)*usr)
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
          uv10 = uv995+(usr/vonkar)*(dlog(10.0D0/zu)- &
                       (psiuo(zetu)-psiuo(zu/l10)))
        else
          uv10 = uv995+(usr/vonkar)*(dlog(10.0D0/zu)+ &
                       d_five*zetu-d_five*zu/l10)
        end if
        !
        !-------------------------------------
        ! Compute atmosphere/ocean/ice fluxes
        !-------------------------------------
        !
        ! heat fluxes
        sent(i) = hsb
        evpr(i) = hlb/Le

        ! drag coefficents
        facttq = dlog(z995*d_half)/dlog(z995/zo)
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
    end subroutine coare3_drv

    pure real(rk8) function psiuo(zet)
      implicit none
      real(rk8) , intent (in) :: zet
      real(rk8) :: x, psik, f, psic, c
      if (zet < 0.0D0) then
        x = (1.0D0-15.0D0*zet)**0.25D0
        psik = 2.0D0*dlog((1.0D0+x)/2.0D0)+dlog((1.0D0+x*x)/2.0D0)- &
               2.0D0*datan(x)+2.0D0*datan(1.0D0)
        x = (1.0D0-10.15D0*zet)**0.3333D0
        psic = 1.5D0*dlog((1.0D0+x+x*x)/3.0D0)-dsqrt(3.0D0)*        &
               datan((1.0D0+2.0D0*x)/dsqrt(3.0D0))+4.0D0*           &
               datan(1.0D0)/dsqrt(3.0D0)
        f = zet*zet/(1.0D0+zet*zet)
        psiuo = (1.0D0-f)*psik+f*psic
      else
        c = dmin1(50.0D0,0.35D0*zet)
        psiuo = -((1.0D0+1.0D0*zet)**1.0D0+0.667D0*                 &
                (zet-14.28D0)/dexp(c)+8.525D0)
      end if
    end function psiuo

    pure real(rk8) function psit(zet)
      implicit none
      real(rk8) , intent (in) :: zet
      real(rk8) :: x, psik, f, psic, c
      if (zet < 0.0D0) then
        x = (1.0D0-15.0D0*zet)**0.5D0
        psik = 2.0D0*dlog((1.0D0+x)/2.0D0)
        x = (1.0D0-34.15D0*zet)**0.3333D0
        psic = 1.5D0*dlog((1.0D0+x+x*x)/3.0D0)-dsqrt(3.0D0)* &
               datan((1.0D0+2.0D0*x)/dsqrt(3.0D0))+4.0D0*    &
               datan(1.0D0)/dsqrt(3.0D0)
        f = zet*zet/(1.0D0+zet*zet)
        psit = (1.0D0-f)*psik+f*psic
      else
        c = dmin1(50.0D0,0.35D0*zet)
        psit = -((1.0D0+2.0D0/3.0D0*zet)**1.5D0+0.667D0*     &
                (zet-14.28D0)/dexp(c)+8.525D0)
      end if
    end function psit

end module mod_ocn_coare
