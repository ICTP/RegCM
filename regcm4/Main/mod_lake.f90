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

      use mod_dynparam , only : lkpts

      implicit none

      private

      public :: lakedrv , initlk

      integer :: iin , iout , lcount

      contains

      subroutine lakedrv(jslc)
      use mod_dynparam
      use mod_param1
      use mod_iunits
      use mod_bats
      use mod_date
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      intent (in) jslc
!
! Local variables
!
      real(8) :: aveice , dayl , evl , flw , fsw , hlat , hsen , prec , &
               & ql , hsnow , tgl , tl , vl , zl
      integer :: ilake , n
!
      do ilake = 1 , iym1
        do n = 1 , nnsg
          if ( veg2d1(n,ilake,jslc).eq.14. ) then
            dayl = (nnnnnn-nstrt0)/4. + (xtime+dtmin)/1440.
            tl = ts1d(n,ilake)
            vl = dsqrt(us1d(ilake)**2+vs1d(ilake)**2)
            zl = z1d(n,ilake)
            ql = qs1d(n,ilake)
            fsw = fsw1d(ilake)
            flw = -1.*flw1d(ilake)
            prec = prca2d(ilake,jslc)      !  units of prec = mm
            hsen = -1.*sent1d(n,ilake)
            hlat = -1.*evpr1d(n,ilake)

            call lake(iutlak,dayl,dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat, &
                    & tgl,evl,prec,aveice,hsnow)
 
            tg1d(n,ilake) = tgl
            tgb1d(n,ilake) = tgl
            if ( aveice.le.10. ) then
              ldoc1d(n,ilake) = 0.
              sice1d(n,ilake) = 0.
              scv1d(n,ilake) = 0.
              sag1d(n,ilake) = 0.
            else
              ldoc1d(n,ilake) = 2.
              sice1d(n,ilake) = aveice  !  units of ice = mm
              scv1d(n,ilake) = hsnow   !  units of snow = mm h2o
            end if
 
          end if
        end do
      end do
 
      end subroutine lakedrv

!
!-----------------------------------------------------------------------
!
      subroutine initlk(veg2d,ix1,jx1)
 
      implicit none
!
! Dummy arguments
!
      integer :: ix1 , jx1
      real(8) , dimension(ix1,jx1) :: veg2d
      intent (in) ix1 , jx1
      intent (out) veg2d
!
! Local variables
!
      integer :: depth , freeze , idata , ilake , j , jlake , lakeset , &
               & n
      real(8) :: eta , hi , hice , hsnow
      real(8) , dimension(400) :: t
!
 
!     ******  lkpts = number of lake points desired

      data hi , hice , hsnow , eta/.01 , 0. , 0. , .5/
 
!     ******  unit number containing lake pt locations, depths
      idata = 40
 
      lcount = 0
      iin = 41
      iout = 42
 
!     read in vegetation type desired for lake points
!     (18 = mixed woodland (no lake model);  14 = inland water (lake
!     model)
      read (idata,99001) lakeset
      print * , '*** lake points set to bats surface type ' , lakeset
 
!     initialize data for lake model
      do n = 1 , lkpts
        read (idata,99002) ilake , jlake , depth
        freeze = 1
 
!       ******     Initially set eta to mean values
        if ( depth.lt.50 ) then
          eta = .7
        else if ( depth.gt.100 ) then
          eta = .3
        else
          eta = .5
        end if
 
!       ******     Initially set lake points isothermal at 6.0 C (June)
        do j = 1 , depth
          t(j) = 6.0
        end do
 
        veg2d(ilake,jlake) = dble(lakeset)
        write (iin) ilake , jlake , depth , freeze , hi , hice , hsnow ,&
                  & eta , (t(j),j=1,depth)
        print * , ilake , jlake , depth , freeze
 
      end do
 
      rewind (iin)
99001 format (i2)
99002 format (3I4)
 
      end subroutine initlk
!
!-----------------------------------------------------------------------
!
      subroutine lake(iutlak,day,dt,ta,ua,za,q,sw,lnet,hsen,hlat,ts,    &
                    & evap,prec,hice,hsnow)
 
      use mod_constants , only : ep2 , tzero , sigm , wlhv , vonkar
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: depmax = 400
!
! Dummy arguments
!
      real(8) :: day , dt , evap , hice , hlat , hsen , hsnow , lnet ,  &
               & prec , q , sw , ta , ts , ua , za
      integer :: iutlak
      intent (in) day , hlat , hsen , iutlak , q , ta , ua , za
      intent (out) ts
      intent (inout) evap , hice , hsnow
!
! Local variables
!
      real(8) , dimension(depmax) :: de , dnsty
      real(8) , dimension(depmax,2) :: t
      integer :: depth , freeze , ilake , j , jlake , k , kmin
      real(8) :: ea , eta , hi , hs , ld , lu , qe , qh , tac , surf ,  &
              &  tcutoff , tk , u2
!
!***  dt:  time step in seconds
!***  surf:surface thickness
!***  dz:  vertical grid spacing in m
!***  zo:  surface roughness length
      real(8) , parameter :: dz = 1.0
      real(8) , parameter :: zo = 0.001
      real(8) , parameter :: z2 = 2.0
!
      surf = 1.0
!
!     interpolate winds at z1 m to 2m via log wind profile
      u2 = ua*dlog(z2/zo)/dlog(za/zo)
 
!******    depth: 1-m slices of lake depth
      read (iin) ilake , jlake , depth , freeze , hi , hice , hsnow ,   &
               & eta , (t(j,1),j=1,depth)
      do k = 1 , depth
        t(k,2) = t(k,1)
      end do
 
      tac = ta - tzero
      tk = tzero + t(1,1)
      lu = -0.97*sigm*tk**4
      ld = lnet - lu
      qe = hlat*wlhv
      qh = hsen
 
!     convert mixing ratio to air vapor pressure
      ea = q*88.0/(ep2+0.378*q)
 
!     ******    Check if conditions exist for lake ice
      tcutoff = -0.001
      if ( (hice.eq.0.0) .and. (t(1,1).gt.tcutoff) ) then
 
!       ******    Calculate eddy diffusivities
        call eddy(dt,surf,dz,vonkar,u2,t,dnsty,de,depth)
 
!       ******    Lake temperature calc using BATS sensible and latent
!       heats
        call temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,depmax,depth)
 
!       ******    Convective mixer
        kmin = 1
        call mixer(kmin,surf,dz,t,dnsty,depmax,depth)
 
      else
 
        call ice(sw,ld,tac,u2,ea,hs,hi,hice,evap,t,depth,prec)
        if ( freeze.eq.0 ) t(1,1) = t(1,2)
 
      end if
 
      write (iout) ilake , jlake , depth , freeze , hi , hice , hsnow , &
                 & eta , (t(j,1),j=1,depth)
 
      write (iutlak) day , ilake , jlake , depth , evap , hi , hice ,   &
                   & hsnow , (t(j,1),j=1,depth)
 
      ts = t(1,1) + tzero
      evap = evap/3600.          !  convert evap from mm/hr to mm/sec
      hice = hice*1000.          !  convert ice  from m to mm
      hsnow = hsnow*100.         !  convert snow from m depth to mm h20
 
      lcount = lcount + 1
      if ( lcount.eq.lkpts ) then
        lcount = 0
        iin = 83 - iin
        iout = 83 - iout
        rewind (iin)
        rewind (iout)
      end if
 
      end subroutine lake
!
!-----------------------------------------------------------------------
!
      subroutine eddy(dt,surf,dz,kv,u2,t,dnsty,de,depth)
 
! Computes density, eddy diffusivity and variable time step
 
      use mod_constants , only : gti , tzero
      implicit none
!
! PARAMETER definitions
!
!cc   dm=5.148e-04       ! value used in prev simulations
      real(8) , parameter :: dm = 1.38889E-07
!
! Dummy arguments
!
      integer :: depth
      real(8) :: dt , dz , kv , surf , u2
      real(8) , dimension(depth) :: de , dnsty
      real(8) , dimension(depth,2) :: t
      intent (in) depth , dt , dz , kv , surf , t
      intent (inout) de , dnsty , u2
!
! Local variables
!
      real(8) :: demax , dpdz , ks , n2 , po , rad , ri , rimax , ws , z
      integer :: k
!
      demax = .5*dz**2/dt
      demax = .99*demax
      rimax = 0.0
      do k = 1 , depth
        dnsty(k) = 1000.0*(1.0-1.9549E-05*(dabs((t(k,1)+tzero)-277.0))  &
                 & **1.68)
      end do
 
      if ( u2.lt.0.5 ) u2 = 0.5
 
!******     compute eddy diffusion profile
!     N2 Brunt-Vaisala number
!     Ri gradient Richardson number
!     dm molecular diffusion of water
!******     compute eddy diffusion profile
 
      ks = 0.745*u2**(-1.84)
      ws = 0.0012*u2
      po = 1.0
 
      do k = 1 , depth - 1
        dpdz = (dnsty(k+1)-dnsty(k))/dz   ! gtb removed /2.0
        n2 = dpdz/dnsty(k)*gti            ! gtb removed minus
        z = surf + dble(k-1)              ! gtb: k was k-1
        rad = 1. + 40.*n2*(kv*z*dexp(ks*z)/ws)**2
        if ( rad.lt.0 ) rad = 0.0
        ri = (-1.0+dsqrt(rad))/20.0
        de(k) = dm + kv*ws*z*po*dexp(-ks*z)/(1.0+37.0*ri**2)
        if ( de(k).gt.demax ) de(k) = demax
        if ( dabs(ri).gt.rimax ) rimax = dabs(ri)
      end do
      de(depth) = 0.0
 
      end subroutine eddy
!
!-----------------------------------------------------------------------
!
      subroutine temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,depmax,   &
                    & depth)
!*****************BEGIN SUBROUTINE TEMP********************
!             COMPUTES TEMPERATURE PROFILE                *
!**********************************************************
      use mod_constants , only : tzero , cpw
      implicit none
!
! Dummy arguments
!
      integer :: depmax , depth
      real(8) :: dt , dz , eta , lnet , qe , qh , surf , sw
      real(8) , dimension(depmax) :: de , dnsty
      real(8) , dimension(depmax,2) :: t
      intent (in) de , depmax , depth , dt , dz , eta , lnet , qe , qh ,&
                & sw
      intent (inout) dnsty , surf , t
!
! Local variables
!
      real(8) :: bot , t1 , t2 , tdiff , top
      integer :: k
 
      surf = 1.0
 
!******    solve differential equations of heat transfer
      do k = 1 , depth
        t(k,2) = t(k,1)
      end do
 
      k = 1
      t1 = sw*(1.-dexp(-eta*surf))/(surf*dnsty(k)*cpw) + (lnet+qe+qh)   &
         & /(surf*dnsty(k)*cpw)
      t2 = -de(k)*(t(k,1)-t(k+1,1))/surf
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      do k = 2 , depth - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        t1 = sw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cpw)
        t2 = (de(k-1)*(t(k-1,1)-t(k,1))-de(k)*(t(k,1)-t(k+1,1)))/dz
        t(k,2) = t(k,2) + (t1+t2)*dt
      end do
 
      k = depth
      top = (surf+(k-2)*dz)
      t1 = sw*dexp(-eta*top)/(dz*dnsty(k)*cpw)
      t2 = de(k-1)*(t(depth-1,1)-t(depth,1))/dz
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      tdiff = 0.
      do k = 1 , depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if ( k.eq.1 ) tdiff = tdiff*surf
        t(k,1) = t(k,2)
        dnsty(k) = 1000.0*(1.0-1.9549E-05*(dabs((t(k,2)+tzero)-277.0))  &
                 & **1.68)
      end do

!sb   print *, 'TEMP: Total temp change = ', tdiff
 
      end subroutine temp
!
!-----------------------------------------------------------------------
!
      subroutine mixer(kmin,surf,dz,t,dnsty,depmax,depth)
!
! Simulates convective mixing
!
      use mod_constants , only : tzero
      implicit none
!
! Dummy arguments
!
      integer :: depmax , depth , kmin
      real(8) :: dz , surf
      real(8) , dimension(depmax) :: dnsty
      real(8) , dimension(depmax,2) :: t
      intent (in) depmax , depth , dz , kmin , surf
      intent (inout) dnsty , t
!
! Local variables
!
      real(8) :: avet , avev , tav , tdiff , vol
      integer :: k , k2 , m
! 
      do k = kmin , depth - 1
        avet = 0.0
        avev = 0.0
 
        if ( dnsty(k).gt.dnsty(k+1) ) then
 
          do m = kmin , k + 1
            if ( m.eq.1 ) then
              vol = surf
            else
              vol = dz
            end if
            avet = avet + t(m,2)*vol
            avev = avev + vol
          end do
 
          tav = avet/avev
          do k2 = kmin , k + 1
            t(k2,2) = tav
            dnsty(k2) = 1000.0*(1.0-1.9549E-05*(dabs((t(k2,2)+tzero)-   &
                      & 277.0))**1.68)
          end do
        end if
 
      end do ! K loop
 
      tdiff = 0.0
      do k = kmin , depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if ( k.eq.1 ) tdiff = tdiff*surf
        t(k,1) = t(k,2)
      end do
 
      end subroutine mixer
!
!-----------------------------------------------------------------------
!
      subroutine ice(kd,ld,ta,u2,ea,hs,hi,hii,evap,t,depth,precip)

      use mod_constants , only : ep2 , tzero , sigm , wlhv
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: surf = 0.6 , lami1 = 1.5 , lami2 = 20 ,    &
                           & lams1 = 6.0 , lams2 = 20.0 , ki = 2.3 ,    &
                           & ks = 0.31 , rhoi = 917.0 , atm = 950 ,     &
                           & qw = 1.389 , rhos = 330.0 , li = 334.0E03 ,&
                           & cd = 0.001 , cp = 1000.0 , rho = 1.0 ,     &
                           & sec = 3600
!
! Dummy arguments
!
      integer :: depth
      real(8) :: ea , evap , hi , hii , hs , kd , ld , precip , ta , u2
      real(8) , dimension(depth,2) :: t
      intent (in) depth , ea , ld , precip , ta , u2
      intent (out) evap
      intent (inout) hi , hii , hs , kd , t
!
! Local variables
!
      real(8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               & t2 , tf , theta , x
      real(8) :: f , t4
      integer :: nits
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           *
!***********************************************************************
 
      t4(x) = (x+tzero)**4
 
!     ****** g. bates changed air to ta, qpen1 to qpen (4/92)
      f(x) = (-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-ta)-kd)    &
           & - 1./khat*(qpen+tf-x)
 
!CC   f(x)=(-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta
!CC   +      *(x-air)-kd)-1/khat*(qpen1+tf-x)
 
      if ( (ta.le.0.0) .and. (hii.gt.0.0) ) hs = hs + precip*10./1000.
                                     ! convert precip(mm) to depth(m)
 
      t0 = t(1,1)
      tf = 0.0
 
      khat = (ki*hs+ks*hi)/(ki*ks)
      theta = cp*rho*cd*u2
      psi = wlhv*rho*cd*u2*ep2/atm
      evap = 100.*psi*(eomb(t0)-ea)/(wlhv*rho)
      qpen = kd*0.7*((1.0-dexp(-lams1*hs))/(ks*lams1)+(dexp(-lams1*hs)) &
           & *(1.0-dexp(-lami1*hi))/(ki*lami1))                         &
           & + kd*0.3*((1.0-dexp(-lams2))/(ks*lams2)+(-lams2*hs)        &
           & *(1.0-dexp(-lami2*hi))/(ki*lami2))
      kd = kd - qpen
 
      nits = 0
      t1 = -50
      f0 = f(t0)
      f1 = f(t1)
      do
        nits = nits + 1
        t2 = t1 - (t1-t0)*f1/(f1-f0)
        if ( dabs((t2-t1)/t1).ge.0.001 ) then
          t0 = t1
          t1 = t2
          f0 = f1
          f1 = f(t1)
          cycle
        end if
 
        t0 = t2
        if ( t0.ge.tf ) then
 
          if ( hs.gt.0. ) then
            ds = sec*                                                   &
               & ((-ld+0.97*sigm*t4(tf)+psi*(eomb(tf)-ea)+theta*(tf-ta)&
               & -kd)-1/khat*(t0-tf+qpen))/(rhos*li)
            if ( ds.gt.0.0 ) ds = 0.0
            hs = hs + ds
            if ( hs.lt.0 ) then
              hs = 0.0
              t(1,1) = (hii*t0+(surf-hii)*t(2,1))/surf
            end if
          end if
          if ( (hs.eq.0.) .and. (hii.gt.0.0) ) then
            di = sec*                                                   &
               & ((-ld+0.97*sigm*t4(tf)+psi*(eomb(tf)-ea)+theta*(tf-ta)&
               & -kd)-1/khat*(t0-tf+qpen))/(rhoi*li)
            if ( di.gt.0 ) di = 0.0
            hi = hi + di
          end if
 
        else if ( t0.lt.tf ) then
 
          q0 = -ld + 0.97*sigm*t4(t0) + psi*(eomb(t0)-ea)              &
             & + theta*(t0-ta) - kd
          qpen = kd*0.7*(1.0-dexp(-(lams1*hs+lami1*hi)))                &
               & + kd*0.3*(1.0-dexp(-(lams2*hs+lami2*hi)))
          di = sec*(q0-qw-qpen)/(rhoi*li)
 
          hi = hi + di
 
        else
        end if
 
        if ( hi.le.0.01 ) then
          hi = 0.01
          hii = 0.0
          hs = 0.0
          t(1,1) = (hi*t0+(surf-hi)*t(2,1))/surf
        else
          hii = hi
          t(1,1) = t0
        end if
        exit
      end do
 
      end subroutine ice
!
! Computes air vapor pressure as a function of temp (in K)
!
      function eomb(x)

      use mod_constants , only : stdpmb , tboil , tzero
      implicit none
!
! Dummy arguments
!
      real(8) :: x
      intent (in) x
      real(8) :: eomb
!
! Local variables
!
      real(8) :: tr1
!
      tr1 = 1.0 - (tboil/(x+tzero))
      eomb = stdpmb*dexp(13.3185*tr1-1.976*tr1**2-0.6445*tr1**3-       &
           & 0.1299*tr1**4)

      end function eomb
!
      end module mod_lake
