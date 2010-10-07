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

      use mod_dynparam

      implicit none
      real(8) ,allocatable, dimension(:,:,:) :: dhlake1
      integer ,allocatable, dimension(:,:,:) :: depth2d
      real(8) ,allocatable, dimension(:,:,:) :: eta2d
      real(8) ,allocatable, dimension(:,:,:) :: hi2d
      real(8) ,allocatable, dimension(:,:,:) :: aveice2d
      real(8) ,allocatable, dimension(:,:,:) :: hsnow2d
      real(8) ,allocatable, dimension(:,:,:) :: evl2d
      real(8) ,allocatable, dimension(:,:,:,:,:) :: tlak3d

      contains
!
!-----------------------------------------------------------------------
!
      subroutine allocate_lake
      implicit none
#ifdef MPP1
      allocate(dhlake1(nnsg,iy,jxp))
      allocate(depth2d(nnsg,iym1,jxp))
      allocate(eta2d(nnsg,iym1,jxp))
      allocate(hi2d(nnsg,iym1,jxp))
      allocate(aveice2d(nnsg,iym1,jxp))
      allocate(hsnow2d(nnsg,iym1,jxp))
      allocate(evl2d(nnsg,iym1,jxp))
      allocate(tlak3d(400,2,nnsg,iym1,jxp))
#else
      allocate(dhlake1(nnsg,iy,jx))
      allocate(depth2d(nnsg,iym1,jx))
      allocate(eta2d(nnsg,iym1,jx))
      allocate(hi2d(nnsg,iym1,jx))
      allocate(aveice2d(nnsg,iym1,jx))
      allocate(hsnow2d(nnsg,iym1,jx))
      allocate(evl2d(nnsg,iym1,jx))
      allocate(tlak3d(400,2,nnsg,iym1,jx))
#endif
      end subroutine allocate_lake

      subroutine initlake
      use mod_dynparam
      use mod_bats, only : satbrt1
#ifdef MPP1
      use mod_mppio, only : depth2d_io
#ifndef IBM
      use mpi
#endif
#endif
      implicit none
#ifdef MPP1
#ifdef IBM
      include 'mpif.h'
#endif
#endif
! 
! Local variables
!
      integer :: i, j, n
#ifdef MPP1
      integer :: ierr
#endif

      hi2d = 0.01
      aveice2d = 0.0
      hsnow2d = 0.0
      eta2d = 0.5
      tlak3d  = 6.
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        do i = 2 , iym1
          do n = 1 , nnsg

!     ******  initialize hostetler lake model
            if( (satbrt1(n,i,j).gt.13.9.and.satbrt1(n,i,j)  &
               & .lt.15.1) .and. dhlake1(n,i,j).gt.-9990.) then
              depth2d(n,i,j) = dmax1(2.d0,dmin1(dhlake1(n,i,j),400.d0))
              if(depth2d(n,i,j).lt.50) then
                eta2d(n,i,j) = .7
              else if(depth2d(n,i,j).gt.100) then
                eta2d(n,i,j) = .3
              else
                eta2d(n,i,j) = .5
              endif
            else
              depth2d(n,i,j) = 0
              eta2d(n,i,j) = 0.5
            endif
          enddo
        enddo
      enddo
#ifdef MPP1
      call mpi_gather(depth2d,nnsg*iym1*jxp,mpi_integer, &
                    & depth2d_io,nnsg*iym1*jxp,mpi_integer, &
                    & 0, mpi_comm_world,ierr)
#endif
      end subroutine initlake
!
      subroutine lakedrv(jslc)
      use mod_dynparam
      use mod_runparams
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
      real(8) :: evl , flw , fsw , hlat , hsen , prec , &
               & ql , tgl , tl , vl , zl
      integer :: i , n
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( depth2d(n,i,jslc).gt.0 ) then
            tl = ts1d(n,i)
            vl = dsqrt(us1d(i)**2+vs1d(i)**2)
            zl = z1d(n,i)
            ql = qs1d(n,i)
            fsw = fsw1d(i)
            flw = -1.*flw1d(i)
            prec = prca2d(i,jslc)      !  units of prec = mm
            hsen = -1.*sent1d(n,i)
            hlat = -1.*evpr1d(n,i)

            call lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat, &
                    & tgl,prec,depth2d(n,i,jslc),eta2d(n,i,jslc),   &
                    & hi2d(n,i,jslc),aveice2d(n,i,jslc),  &
                    & hsnow2d(n,i,jslc),evl2d(n,i,jslc),tlak3d(1,1,n,i,jslc) )
 
            tg1d(n,i) = tgl
            tgb1d(n,i) = tgl
            if ( aveice2d(n,i,jslc).le.10. ) then
              ldoc1d(n,i) = 0.
              sice1d(n,i) = 0.
              scv1d(n,i) = 0.
              sag1d(n,i) = 0.
            else
              ldoc1d(n,i) = 2.
              sice1d(n,i) = aveice2d(n,i,jslc)  !  units of ice = mm
              scv1d(n,i) = hsnow2d(n,i,jslc)    !  units of snow = mm h2o
            end if
 
          end if
        end do
      end do
 
      end subroutine lakedrv
!
!-----------------------------------------------------------------------
!
      subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat,tgl,  &
                    & prec,depth,eta,hi,aveice,hsnow,evl,t)
 
      use mod_constants , only : ep2 , tzero , sigm , wlhv
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: depmax = 400
!
! Dummy arguments
!
      real(8) :: dtlake , evl , aveice , hlat , hsen , hsnow , flw , &
               & prec , ql , fsw , tl , tgl , vl , zl , eta,hi
      integer :: depth
      intent (in) hlat , hsen , ql , tl , vl , zl
      intent (in) depth,eta
      intent (out) tgl
      intent (inout) evl , aveice , hsnow
      real(8) , dimension(depmax,2) :: t
      intent (inout) t 
!
! Local variables
!
      real(8) , dimension(depmax) :: de , dnsty
      integer :: freeze , j , k , kmin
      real(8) :: ea , hs , ld , lu , qe , qh , tac , surf ,  &
              &  tcutoff , tk , u2
!
!***  dtlake:  time step in seconds
!***  surf:surface thickness
!***  dz:  vertical grid spacing in m
!***  zo:  surface roughness length
      real(8) , parameter :: dz = 1.0
      real(8) , parameter :: zo = 0.001
      real(8) , parameter :: z2 = 2.0
!
      freeze = 1
      surf = 1.0
!
!     interpolate winds at z1 m to 2m via log wind profile
      u2 = vl*dlog(z2/zo)/dlog(zl/zo)
 
!******    depth: 1-m slices of lake depth
      do k = 1 , depth
        t(k,2) = t(k,1)
      end do
 
      tac = tl - tzero
      tk = tzero + t(1,1)
      lu = -0.97*sigm*tk**4
      ld = flw - lu
      qe = hlat*wlhv
      qh = hsen
 
!     convert mixing ratio to air vapor pressure
      ea = ql*88.0/(ep2+0.378*ql)
 
!     ******    Check if conditions exist for lake ice
      tcutoff = -0.001
      if ( (aveice.eq.0.0) .and. (t(1,1).gt.tcutoff) ) then
 
!       ******    Calculate eddy diffusivities
        call eddy(dtlake,surf,dz,u2,t,dnsty,de,depth)
 
!       ******    Lake temperature calc using BATS sensible and latent
!       heats
        call temp(dtlake,surf,dz,t,fsw,flw,qe,qh,dnsty,de,eta,depmax,depth)
 
!       ******    Convective mixer
        kmin = 1
        call mixer(kmin,surf,dz,t,dnsty,depmax,depth)
 
      else
 
        call ice(fsw,ld,tac,u2,ea,hs,hi,aveice,evl,t,depth,prec)
        if ( freeze.eq.0 ) t(1,1) = t(1,2)
 
      end if
 
      tgl = t(1,1) + tzero
      evl = evl/3600.          !  convert evl from mm/hr to mm/sec
      aveice = aveice*1000.          !  convert ice  from m to mm
      hsnow = hsnow*100.         !  convert snow from m depth to mm h20
 
      end subroutine lake
!
!-----------------------------------------------------------------------
!
      subroutine eddy(dtlake,surf,dz,u2,t,dnsty,de,depth)
 
! Computes density, eddy diffusivity and variable time step
 
      use mod_constants , only : gti , tzero , vonkar
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
      real(8) :: dtlake , dz , surf , u2
      real(8) , dimension(depth) :: de , dnsty
      real(8) , dimension(depth,2) :: t
      intent (in) depth , dtlake , dz , surf , t
      intent (inout) de , dnsty , u2
!
! Local variables
!
      real(8) :: demax , dpdz , ks , n2 , po , rad , ri , rimax , ws , z
      integer :: k
!
      demax = .5*dz**2/dtlake
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
        rad = 1. + 40.*n2*(vonkar*z*dexp(ks*z)/ws)**2
        if ( rad.lt.0 ) rad = 0.0
        ri = (-1.0+dsqrt(rad))/20.0
        de(k) = dm + vonkar*ws*z*po*dexp(-ks*z)/(1.0+37.0*ri**2)
        if ( de(k).gt.demax ) de(k) = demax
        if ( dabs(ri).gt.rimax ) rimax = dabs(ri)
      end do
      de(depth) = 0.0
 
      end subroutine eddy
!
!-----------------------------------------------------------------------
!
      subroutine temp(dtlake,surf,dz,t,fsw,flw,qe,qh,dnsty,de,eta,depmax,   &
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
      real(8) :: dtlake , dz , eta , flw , qe , qh , surf , fsw
      real(8) , dimension(depmax) :: de , dnsty
      real(8) , dimension(depmax,2) :: t
      intent (in) de , depmax , depth , dtlake , dz , eta , flw , qe , qh ,&
                & fsw
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
      t1 = fsw*(1.-dexp(-eta*surf))/(surf*dnsty(k)*cpw) + (flw+qe+qh)   &
         & /(surf*dnsty(k)*cpw)
      t2 = -de(k)*(t(k,1)-t(k+1,1))/surf
      t(k,2) = t(k,2) + (t1+t2)*dtlake
 
      do k = 2 , depth - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        t1 = fsw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cpw)
        t2 = (de(k-1)*(t(k-1,1)-t(k,1))-de(k)*(t(k,1)-t(k+1,1)))/dz
        t(k,2) = t(k,2) + (t1+t2)*dtlake
      end do
 
      k = depth
      top = (surf+(k-2)*dz)
      t1 = fsw*dexp(-eta*top)/(dz*dnsty(k)*cpw)
      t2 = de(k-1)*(t(depth-1,1)-t(depth,1))/dz
      t(k,2) = t(k,2) + (t1+t2)*dtlake
 
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
      subroutine ice(fsw,ld,tac,u2,ea,hs,hi,aveice,evl,t,depth,prec)

      use mod_constants , only : ep2 , tzero , sigm , wlhv,stdpmb,tboil
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
      real(8) :: ea , evl , hi , aveice , hs , fsw , ld , prec , tac , u2
      real(8) , dimension(depth,2) :: t
      intent (in) depth , ea , ld , prec , tac , u2
      intent (out) evl
      intent (inout) hi , aveice , hs , fsw , t
!
! Local variables
!
      real(8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               & t2 , tf , theta , x
      real(8) :: f , t4 , tr1 , eomb
      integer :: nits
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           
!**********************************************************************
 
      t4(x) = (x+tzero)**4
 
! Computes air vapor pressure as a function of temp (in K)
      tr1(x) = 1.0 - (tboil/(x+tzero))
      eomb(x) = stdpmb*dexp(13.3185*tr1(x)-1.976*tr1(x)**2   &
           &   -0.6445*tr1(x)**3- 0.1299*tr1(x)**4)

!     ****** g. bates changed air to tac, qpen1 to qpen (4/92)
      f(x) = (-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-tac)-fsw)    &
           & - 1./khat*(qpen+tf-x)
 
!CC   f(x)=(-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta
!CC   +      *(x-air)-fsw)-1/khat*(qpen1+tf-x)
 
      if ( (tac.le.0.0) .and. (aveice.gt.0.0) ) hs = hs + prec*10./1000.
                                     ! convert prec(mm) to depth(m)
 
      t0 = t(1,1)
      tf = 0.0
 
      khat = (ki*hs+ks*hi)/(ki*ks)
      theta = cp*rho*cd*u2
      psi = wlhv*rho*cd*u2*ep2/atm
      evl = 100.*psi*(eomb(t0)-ea)/(wlhv*rho)
      qpen = fsw*0.7*((1.0-dexp(-lams1*hs))/(ks*lams1)+(dexp(-lams1*hs)) &
           & *(1.0-dexp(-lami1*hi))/(ki*lami1))                         &
           & + fsw*0.3*((1.0-dexp(-lams2))/(ks*lams2)+(-lams2*hs)        &
           & *(1.0-dexp(-lami2*hi))/(ki*lami2))
      fsw = fsw - qpen
 
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
               & ((-ld+0.97*sigm*t4(tf)+psi*(eomb(tf)-ea)+theta*(tf-tac)&
               & -fsw)-1/khat*(t0-tf+qpen))/(rhos*li)
            if ( ds.gt.0.0 ) ds = 0.0
            hs = hs + ds
            if ( hs.lt.0 ) then
              hs = 0.0
              t(1,1) = (aveice*t0+(surf-aveice)*t(2,1))/surf
            end if
          end if
          if ( (hs.eq.0.) .and. (aveice.gt.0.0) ) then
            di = sec*                                                   &
              & ((-ld+0.97*sigm*t4(tf)+psi*(eomb(tf)-ea)+theta*(tf-tac) &
              & -fsw)-1/khat*(t0-tf+qpen))/(rhoi*li)
            if ( di.gt.0 ) di = 0.0
            hi = hi + di
          end if
 
        else if ( t0.lt.tf ) then
 
          q0 = -ld + 0.97*sigm*t4(t0) + psi*(eomb(t0)-ea)              &
             & + theta*(t0-tac) - fsw
          qpen = fsw*0.7*(1.0-dexp(-(lams1*hs+lami1*hi)))               &
               & + fsw*0.3*(1.0-dexp(-(lams2*hs+lami2*hi)))
          di = sec*(q0-qw-qpen)/(rhoi*li)
 
          hi = hi + di
 
        else
        end if
 
        if ( hi.le.0.01 ) then
          hi = 0.01
          aveice = 0.0
          hs = 0.0
          t(1,1) = (hi*t0+(surf-hi)*t(2,1))/surf
        else
          aveice = hi
          t(1,1) = t0
        end if
        exit
      end do
 
      end subroutine ice
!
!-----------------------------------------------------------------------
!
      subroutine outlake

      use mod_dynparam
      use mod_runparams
      use mod_date
#ifdef MPP1
      use mod_mppio, only : depth2d_io, eta2d_io, hi2d_io, &
                   &        aveice2d_io, hsnow2d_io, evl2d_io, &
                   &        tlak3d_io 
#ifndef IBM
      use mpi
#endif
#endif
      implicit none
#ifdef MPP1
#ifdef IBM
      include 'mpif.h'
#endif
#endif
!
! Local variables
!
      real(8) :: dayl
      integer :: i, j, k, n
#ifdef MPP1
      integer :: ierr
#endif
!
      dayl = (nnnnnn-nstrt0)/4. + (xtime+dtmin)/1440.

#ifdef MPP1
      call mpi_gather(eta2d,nnsg*iym1*jxp,mpi_real8, &
                    & eta2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(hi2d,nnsg*iym1*jxp,mpi_real8, &
                    & hi2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(aveice2d,nnsg*iym1*jxp,mpi_real8, &
                    & aveice2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(hsnow2d,nnsg*iym1*jxp,mpi_real8, &
                    & hsnow2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(evl2d,nnsg*iym1*jxp,mpi_real8, &
                    & evl2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(tlak3d,400*2*nnsg*iym1*jxp,mpi_real8, &
                    & tlak3d_io,400*2*nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d_io(n,i,j).gt.0 ) then
              write(58) dayl, i, j, depth2d_io(n,i,j), evl2d_io(n,i,j), &
                   & hi2d_io(n,i,j), aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                   & (tlak3d_io(k,1,n,i,j),k=1,depth2d_io(n,i,j))  
            endif
          enddo
        enddo
      enddo
      endif

#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1 
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d(n,i,j).gt.0 ) then
              write(58) dayl, i, j, depth2d(n,i,j), evl2d(n,i,j), &
                   & hi2d(n,i,j), aveice2d(n,i,j), hsnow2d(n,i,j), &
                   & (tlak3d(k,1,n,i,j),k=1,depth2d(n,i,j))  
            endif
          enddo
        enddo
      enddo
#endif

      end subroutine outlake
!
!-----------------------------------------------------------------------
!
      subroutine lakesav0_o

      use mod_dynparam
#ifdef MPP1
      use mod_mppio, only : eta2d_io, hi2d_io, aveice2d_io, &
                       &    hsnow2d_io, evl2d_io, tlak3d_io 
#ifndef IBM
      use mpi
#endif
#endif
      implicit none
#ifdef MPP1
#ifdef IBM
      include 'mpif.h'
#endif
#endif
!
! Local variables
!
#ifdef MPP1
      integer :: ierr
#endif
!
#ifdef MPP1
      call mpi_gather(eta2d,nnsg*iym1*jxp,mpi_real8, &
                    & eta2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(hi2d,nnsg*iym1*jxp,mpi_real8, &
                    & hi2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(aveice2d,nnsg*iym1*jxp,mpi_real8, &
                    & aveice2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(hsnow2d,nnsg*iym1*jxp,mpi_real8, &
                    & hsnow2d_io,nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      call mpi_gather(tlak3d,400*2*nnsg*iym1*jxp,mpi_real8, &
                    & tlak3d_io,400*2*nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
#endif

      end subroutine lakesav0_o
!
!-----------------------------------------------------------------------
!
      subroutine lakesav_o(iutl)

      use mod_dynparam
#ifdef MPP1
      use mod_mppio, only : depth2d_io, eta2d_io, hi2d_io, &
                   &        aveice2d_io, hsnow2d_io, evl2d_io, &
                   &        tlak3d_io 
#endif
      implicit none
      integer :: iutl
      intent (in) iutl
!
! Local variables
!
      integer :: i, j, k, n, numpts
!
#ifdef MPP1
#ifdef BAND
      write (iutl) (((depth2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      write (iutl) (((depth2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      numpts = 0
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d_io(n,i,j).gt.0 ) then
              numpts = numpts+1
            endif
          enddo
        enddo
      enddo
      write (iutl) numpts
      print * , 'writing lake model restart file. numpts = ' , numpts
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d_io(n,i,j).gt.0 ) then
              write(iutl) n, i, j, depth2d_io(n,i,j), eta2d_io(n,i,j), &
                   & hi2d_io(n,i,j), &
                   & aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                   & (tlak3d_io(k,1,n,i,j),k=1,depth2d_io(n,i,j))  
            endif
          enddo
        enddo
      enddo

#else
#ifdef BAND
      write (iutl) (((depth2d(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      write (iutl) (((depth2d(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      numpts = 0
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d(n,i,j).gt.0 ) then
              numpts = numpts+1
            endif
          enddo
        enddo
      enddo
      write (iutl) numpts
      print * , 'writing lake model restart file. numpts = ' , numpts
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1 
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( depth2d(n,i,j).gt.0 ) then
              write(iutl) n, i, j, depth2d(n,i,j), eta2d(n,i,j), &
                   & hi2d(n,i,j), &
                   & aveice2d(n,i,j), hsnow2d(n,i,j), &
                   & (tlak3d(k,1,n,i,j),k=1,depth2d(n,i,j))  
            endif
          enddo
        enddo
      enddo
#endif

      end subroutine lakesav_o
!
!-----------------------------------------------------------------------
!
      subroutine lakesav_i(iutl)

      use mod_dynparam
#ifdef MPP1
      use mod_mppio, only : depth2d_io, eta2d_io, hi2d_io, &
                   &        aveice2d_io, hsnow2d_io, tlak3d_io 
#endif
      implicit none
      integer :: iutl
      intent (in) iutl
!
! Local variables
!
      integer :: i, j, k, l, n, numpts
!

#ifdef MPP1
      if ( myid.eq.0 ) then
#ifdef BAND
        read (iutl) (((depth2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
        read (iutl) (((depth2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif

        hi2d_io = 0.01
        aveice2d_io = 0.0
        hsnow2d_io = 0.0
        eta2d_io = 0.5
        tlak3d_io  = 6.

        read (iutl) numpts
        print * , 'reading lake model restart file. numpts = ' , numpts
        do l = 1, numpts
          read(iutl) n, i, j, depth2d_io(n,i,j), eta2d_io(n,i,j), &
                   & hi2d_io(n,i,j), &
                   & aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                   & (tlak3d_io(k,1,n,i,j),k=1,depth2d_io(n,i,j))  
        enddo
      endif
#else
#ifdef BAND
      read (iutl) (((depth2d(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      read (iutl) (((depth2d(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      read (iutl) numpts
      print * , 'reading lake model restart file. numpts = ' , numpts
      do l = 1, numpts
        read(iutl) n, i, j, depth2d(n,i,j), eta2d(n,i,j), &
                 & hi2d(n,i,j), &
                 & aveice2d(n,i,j), hsnow2d(n,i,j), &
                 & (tlak3d(k,1,n,i,j),k=1,depth2d(n,i,j))  
      enddo
#endif

      end subroutine lakesav_i
!
!-----------------------------------------------------------------------
!
      subroutine lakesav0_i

      use mod_dynparam
#ifdef MPP1
      use mod_mppio, only : depth2d_io, eta2d_io, hi2d_io, &
                   &        aveice2d_io, hsnow2d_io, tlak3d_io 
#ifndef IBM
      use mpi
#endif
#endif
      implicit none
#ifdef MPP1
#ifdef IBM
      include 'mpif.h'
#endif
#endif
!
! Local variables
!
#ifdef MPP1
      integer :: ierr
#endif
!

#ifdef MPP1
      call mpi_scatter(depth2d_io,nnsg*iym1*jxp,mpi_integer, &
                     & depth2d,nnsg*iym1*jxp,mpi_integer, &
                     & 0, mpi_comm_world,ierr)
      call mpi_scatter(eta2d_io,nnsg*iym1*jxp,mpi_real8, &
                     & eta2d,nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
      call mpi_scatter(hi2d_io,nnsg*iym1*jxp,mpi_real8, &
                     & hi2d,nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
      call mpi_scatter(aveice2d_io,nnsg*iym1*jxp,mpi_real8, &
                     & aveice2d,nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
      call mpi_scatter(hsnow2d_io,nnsg*iym1*jxp,mpi_real8, &
                     & hsnow2d,nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
      call mpi_scatter(tlak3d_io,400*2*nnsg*iym1*jxp,mpi_real8, &
                     & tlak3d,400*2*nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
#endif

      end subroutine lakesav0_i
!
      end module mod_lake
