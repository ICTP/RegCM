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
!     LAKE MODEL
!
      module mod_lake

      use mod_constants
!
      private
!
      public :: initlake , lake
!
      integer , parameter :: ndpmax = 400
!
      real(8) , dimension(ndpmax) :: tprof
      integer :: idep
      real(8) :: eta
      real(8) :: hi
      real(8) :: aveice
      real(8) :: hsnow
      real(8) :: tgl
!
      real(8) , dimension(ndpmax) :: de , dnsty , tt
!
      public :: tprof , idep , aveice , hsnow , tgl
!
!     surface thickness
      real(8) , parameter :: surf = 1.0D0
!     vertical grid spacing in m
      real(8) , parameter :: dz = surf
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine initlake(depth)
      implicit none
! 
      real(8) , intent(in) :: depth
      integer :: i, j, n

      idep   = int(max(2.D0,min(depth,dble(ndpmax))))
      hi     = 0.01D0
      aveice = 0.0D0
      hsnow  = 0.0D0
      eta    = 0.5D0
      tprof(1:idep) = 6.0D0
      tprof(idep+1:ndpmax) = -1D+34

      end subroutine initlake
!
!-----------------------------------------------------------------------
!
      subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen, &
                      xl,prec,evl)
! 
      implicit none
!
      real(8) :: dtlake , evl , hsen , flw , &
               & prec , ql , fsw , tl , vl , zl , xl
      intent (in) hsen , ql , tl , vl , zl
      intent (inout) evl
!
      real(8) :: ai , ea , ev , hs , ld , lu , qe , qh , tac , tk , u2
!
!***  dtlake:  time step in seconds
!***  zo:      surface roughness length
!
      real(8) , parameter :: zo = 0.001D0
      real(8) , parameter :: z2 = 2.0D0
      real(8) , parameter :: tcutoff = -0.001D0
      real(8) , parameter :: twatui = 1.78D0
      logical , parameter :: lfreeze = .false.
      integer , parameter :: kmin = 1
!
!     interpolate winds at z1 m to 2m via log wind profile
      u2 = vl*log(z2/zo)/log(zl/zo)
      if ( u2 < d_half ) u2 = d_half
 
!     ****** Check if conditions not exist for lake ice
      if ( (aveice < 1.0D-8) .and. (tprof(1) > tcutoff) ) then
 
        qe = -d_one*evl*wlhv
        qh = hsen

!       ******    Calculate eddy diffusivities
        call eddy(idep,dtlake,u2,xl,tprof)
 
!       ******    Lake temperature calc using sensible and latent heats
        call temp(idep,dtlake,fsw,flw,qe,qh,eta,tprof)
 
!       ******    Convective mixer
        call mixer(kmin,idep,tprof)

        hi     = 0.01D0
        aveice = d_zero
        hsnow  = d_zero

!     ****** Lake ice
      else
 
!       convert mixing ratio to air vapor pressure
        ea  = ql*88.0D0/(ep2+0.378D0*ql)
        tac = tl - tzero
        tk  = tzero + tprof(1)
        lu  = -emsw*sigm*tk**d_four
        ld  = flw - lu
        ev  = evl*secph         ! convert to mm/hr
        ai  = aveice * d_r1000  ! convert to m
        hs  = hsnow * d_r100    ! convert to m

        call ice(dtlake,fsw,ld,tac,u2,ea,hs,hi,ai,ev,prec,tprof)
        if ( .not. lfreeze ) tprof(1) = twatui

        evl    = ev/secph       ! convert evl  from mm/hr to mm/sec
        aveice = ai*d_1000      ! convert ice  from m to mm
        hsnow  = hs*d_100       ! convert snow from m depth to mm h20
        if (aveice < dlowval) aveice = d_zero
        if (hsnow < dlowval) hsnow = d_zero
 
      end if
 
      tgl = tprof(1) + tzero
 
      end subroutine lake
!
!-----------------------------------------------------------------------
!
      subroutine eddy(ndpt,dtlake,u2,xl,tprof)
 
! Computes density and eddy diffusivity
 
      implicit none
!
      integer , intent (in) :: ndpt
      real(8) , intent (in) :: dtlake , u2 , xl
      real(8) , dimension(ndpmax) , intent (in) :: tprof
!
      real(8) :: demax , demin , dpdz , ks , n2 , po
      real(8) :: zmax , rad , ri , ws , z
      integer :: k
!
!     demin molecular diffusion of heat in water
      demin = hdmw
!
!     Added to keep numerical stability of code
      demax = .50D0*dz**d_two/dtlake
      demax = .99D0*demax
!
      do k = 1 , ndpt
        dnsty(k) = d_1000*(d_one-1.9549D-05 * &
                      (dabs((tprof(k)+tzero)-277.0D0))**1.68D0)
      end do
! 
! Compute eddy diffusion profile
!
! Reference:
!
! B. Henderson-Sellers
!  New formulation of eddy diffusion thermocline models.
!  Appl. Math. Modelling, 1985, Vol. 9 December, pp. 441-446
!
 
!     Decay constant of shear velocity - Ekman profile parameter
      ks = 6.6D0*dsqrt(dsin(xl*degrad))*u2**(-1.84D0)

!     Ekman layer depth where eddy diffusion happens
      zmax = dble(ceiling(surf+40.0D0/(vonkar*ks)))

!     Surface shear velocity
      ws = 0.0012D0*u2

!     Inverse of turbulent Prandtl number
      po = d_one
 
      do k = 1 , ndpt - 1

!       Actual depth from surface
        z = surf + dble(k-1)*dz
        if (z >= zmax) then
          de(k) = demin
          cycle
        end if

        if ( k == 1 ) then
          dpdz = (dnsty(k+1)-dnsty(k))/surf
        else
          dpdz = (dnsty(k+1)-dnsty(k))/dz
        end if

!       Brunt Vaisala frequency squared : we do not mind stability,
!       we just look for energy here.
!        n2 = dabs((dpdz/dnsty(k))*egrav)
        n2 = (dpdz/dnsty(k))*egrav
        if (dabs(n2) < dlowval) then
          de(k) = demin
          cycle
        end if

!       Richardson number estimate
        rad = d_one+40.0D0*n2*((vonkar*z)/(ws*dexp(-ks*z)))**d_two
        if (rad < d_zero) rad = d_zero
        ri = (-d_one+dsqrt(rad))/20.0D0

!       Total diffusion coefficient for heat: molecular + eddy (Eqn 42)
        de(k) = demin + vonkar*ws*z*po*dexp(-ks*z) / &
                        (d_one+37.0D0*ri**d_two)
        if ( de(k) < demin ) de(k) = demin
        if ( de(k) > demax ) de(k) = demax

      end do
      de(ndpt) = demin
 
      end subroutine eddy
!
!-----------------------------------------------------------------------
!
      subroutine temp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
!
!*****************BEGIN SUBROUTINE TEMP********************
!             COMPUTES TEMPERATURE PROFILE                *
!**********************************************************
!
      implicit none
!
      integer , intent(in) :: ndpt
      real(8) , intent(in) :: dtlake , eta , flw , qe , qh , fsw
      real(8) , dimension(ndpmax) , intent(inout) :: tprof
!
      real(8) :: bot , dt1 , dt2 , top
      integer :: k
 
!******    solve differential equations of heat transfer

      tt(1:ndpt) = tprof(1:ndpt)
 
      dt1 = (fsw*(d_one-dexp(-eta*surf))+(flw+qe+qh)) / &
              (surf*dnsty(1)*cpw)
      dt2 = -de(1)*(tprof(1)-tprof(2))/surf
      tt(1) = tt(1) + (dt1+dt2)*dtlake
 
      do k = 2 , ndpt - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        dt1 = fsw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cpw)
        dt2 = (de(k-1)*(tprof(k-1)-tprof(k))    -    &
               de(k)  *(tprof(k)  -tprof(k+1))) / dz
        tt(k) = tt(k) + (dt1+dt2)*dtlake
      end do
 
      top = (surf+(ndpt-2)*dz)
      dt1 = fsw*dexp(-eta*top)/(dz*dnsty(ndpt)*cpw)
      dt2 = de(ndpt-1)*(tprof(ndpt-1)-tprof(ndpt))/dz
      tt(ndpt) = tt(ndpt) + (dt1+dt2)*dtlake
 
      do k = 1 , ndpt
        tprof(k) = tt(k)
        dnsty(k) = d_1000*(d_one-1.9549D-05 * &
                   (dabs((tprof(k)+tzero)-277.0D0))**1.68D0)
      end do

      end subroutine temp
!
!-----------------------------------------------------------------------
!
      subroutine mixer(kmin,ndpt,tprof)
!
! Simulates convective mixing
!
      implicit none
!
      integer , intent(in) :: ndpt , kmin
      real(8) , intent(inout) , dimension(ndpmax) :: tprof
!
      real(8) :: avet , avev , tav , vol
      integer :: k , k2
! 
      tt(kmin:ndpt) = tprof(kmin:ndpt)
 
      do k = kmin , ndpt - 1
        avet = d_zero
        avev = d_zero
 
        if ( dnsty(k) > dnsty(k+1) ) then
 
          do k2 = kmin , k + 1
            if ( k2 == 1 ) then
              vol = surf
            else
              vol = dz
            end if
            avet = avet + tt(k2)*vol
            avev = avev + vol
          end do
 
          tav = avet/avev

          do k2 = kmin , k + 1
            tt(k2) = tav
            dnsty(k2) = d_1000*(d_one-1.9549D-05 * &
                        (dabs((tav+tzero)-277.0D0))**1.68D0)
          end do
        end if
 
      end do ! K loop
 
      tprof(kmin:ndpt) = tt(kmin:ndpt)
 
      end subroutine mixer
!
!-----------------------------------------------------------------------
!
      subroutine ice(dtx,fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,tprof)

      implicit none
      real(8) :: ea , evl , hi , aveice , hs , fsw , &
                 ld , prec , tac , u2 , dtx
      real(8) , dimension(ndpmax) :: tprof
      intent (in) dtx , ea , ld , prec , tac , u2
      intent (out) evl
      intent (inout) hi , aveice , hs , fsw , tprof
!
      real(8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               & t2 , tf , theta , rho , xlexpc
      real(8) :: xea , xeb , xec
      integer :: nits
!
      real(8) , parameter :: isurf = 0.6D0
      ! attenuation coeff for ice in visible band (m-1)
      real(8) , parameter :: lami1 = 1.5D0
      ! attenuation coeff for ice in infrared band (m-1)
      real(8) , parameter :: lami2 = 20.0D0
      ! attenuation coeff for snow in visible band (m-1)
      real(8) , parameter :: lams1 = 6.0D0
      ! attenuation coeff for snow in infrared band (m-1)
      real(8) , parameter :: lams2 = 20.0D0
      ! thermal conductivity of ice (W/m/C)
      real(8) , parameter :: ki = 2.3D0
      ! thermal conductivity of snow (W/m/C)
      real(8) , parameter :: ks = 0.31D0
      ! standard atmospheric pressure (hPa) ????
      real(8) , parameter :: atm = 950.0D0
      ! heat flux from water to ice (w/m2) ???
      real(8) , parameter :: qw = 1.389D0
      ! latent heat of fusion (J/kg)
      real(8) , parameter :: li = 334.0D03
      ! drag coefficient for the turbulent momentum flux.
      real(8) , parameter :: cd = 0.001D0
      ! Maximum exponent
      real(8) , parameter :: minexp = -25.0D0
!
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           
!**********************************************************************
 
      if ( (tac <= d_zero) .and. (aveice > d_zero) ) &
        hs = hs + prec*d_r100  ! convert prec(mm) to depth(m)
      if ( hs < dlowval ) hs = d_zero
 
      ! temperature of ice/snow surface
      t0 = tprof(1)
      ! freezing temp of water
      tf = d_zero
      ! approximate density of air (1 kg/m3)
      rho = rhoh2o*d_r1000
 
      khat = (ki*hs+ks*hi)/(ki*ks)
      theta = cpd*rho*cd*u2
      psi = wlhv*rho*cd*u2*ep2/atm
      evl = d_100*psi*(eomb(t0)-ea)/(wlhv*rho)
      ! amount of radiation that penetrates through the ice (W/m2)
      xea = -lams1*hs
      xeb = -lami1*hi
      xec = -lami2*hi
      if ( xea > minexp ) then
        xea = dexp(xea)
      else
        xea = d_zero
      end if
      if ( xeb > minexp ) then
        xeb = dexp(xeb)
      else
        xeb = d_zero
      end if
      if ( xec > minexp ) then
        xec = dexp(xec)
      else
        xec = d_zero
      end if

      qpen = fsw*0.7D0*((d_one-xea)/(ks*lams1) +            &
                        (xea*(d_one-xeb)/(ki*lami1))) +     &
             fsw*0.3D0*((d_one-dexp(-lams2))/(ks*lams2)+    &
                        (-lams2*hs)*(d_one-xec)/(ki*lami2))
      ! radiation absorbed at the ice surface
      fsw = fsw - qpen
 
      ! test qpen sensitivity
      !qpen = qpen * 0.5

      nits = 0
      t1 = -50.0D0
      f0 = f(t0)
      f1 = f(t1)
      do
        nits = nits + 1
        t2 = t1 - (t1-t0)*f1/(f1-f0)
        if ( dabs((t2-t1)/t1) >= 0.001D0 ) then
          t0 = t1
          t1 = t2
          f0 = f1
          f1 = f(t1)
          cycle
        end if
 
        t0 = t2
        if ( t0 >= tf ) then
 
          if ( hs > d_zero ) then
            ds = dtx*                                            &
               & ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea)+     &
               &  theta*(tf-tac)-fsw)-d_one/khat*(tf-t0+qpen)) / &
               & (rhosnow*li)
            if ( ds > d_zero ) ds = d_zero
            hs = hs + ds
            if ( hs < d_zero ) then
              hs = d_zero
              tprof(1) = (aveice*t0+(isurf-aveice)*tprof(2))/isurf
            end if
          end if
          if ( (dabs(hs) < dlowval) .and. (aveice > d_zero) ) then
            di = dtx*                                        &
              & ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea) + &
                 theta*(tf-tac)-fsw)-d_one/khat*(tf-t0+qpen))/ &
                 (rhoice*li)
            if ( di > d_zero ) di = d_zero
            hi = hi + di
          end if
 
        else if ( t0 < tf ) then
 
          q0 = -ld + 0.97D0*sigm*t4(t0) + psi*(eomb(t0)-ea) + &
               theta*(t0-tac) - fsw
          xlexpc = -(lams1*hs+lami1*hi)
          ! Graziano : limit exponential
          if (xlexpc > minexp ) then
            qpen = fsw*0.7D0*(d_one-dexp(-(lams1*hs+lami1*hi))) + &
                   fsw*0.3D0*(d_one-dexp(-(lams2*hs+lami2*hi)))
          else
            qpen = fsw
          end if
          di = dtx*(q0-qw-qpen)/(rhoice*li)
 
          hi = hi + di
        end if
 
        if ( hi <= 0.01D0 ) then
          hi = 0.01D0
          aveice = d_zero
          hs = d_zero
          tprof(1) = (hi*t0+(isurf-hi)*tprof(2))/isurf
        else
          aveice = hi
          tprof(1) = t0
        end if
        exit
      end do
 
      contains

      function t4(x)
        implicit none
        real(8) :: t4
        real(8) , intent(in) :: x
        t4 = (x+tzero)**d_four
      end function t4
      ! Computes air vapor pressure as a function of temp (in K)
      function tr1(x)
        implicit none
        real(8) :: tr1
        real(8) , intent(in) :: x
        tr1 = d_one - (tboil/(x+tzero))
      end function tr1
      function eomb(x)
        implicit none
        real(8) :: eomb
        real(8) , intent(in) :: x
        eomb = stdpmb*dexp(13.3185D0*tr1(x)-1.976D0*tr1(x)**d_two   &
           &   -0.6445D0*tr1(x)**d_three- 0.1299D0*tr1(x)**d_four)
       end function eomb
      function f(x)
        implicit none
        real(8) :: f
        real(8) , intent(in) :: x
        f = (-ld+0.97D0*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-tac)-fsw)  &
            - d_one/khat*(qpen+tf-x)
      end function f
 
      end subroutine ice
!
      end module mod_lake
