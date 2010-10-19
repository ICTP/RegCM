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
      use mod_dynparam
      use mod_runparams
      use mod_bats
      use mod_date
#ifdef MPP1
      use mod_mppio
#endif
!
      private
!
      public :: allocate_lake , lakesav_i, lakesav_o , lakesav0_i
      public :: initlake , outlake , lakedrv
      public :: dhlake1
!
      real(8) , allocatable , dimension(:,:,:) :: dhlake1
      integer , allocatable , dimension(:,:,:) :: idep2d
      real(8) , allocatable , dimension(:,:,:) :: eta2d
      real(8) , allocatable , dimension(:,:,:) :: hi2d
      real(8) , allocatable , dimension(:,:,:) :: aveice2d
      real(8) , allocatable , dimension(:,:,:) :: hsnow2d
      real(8) , allocatable , dimension(:,:,:) :: evl2d
      real(8) , allocatable , dimension(:,:,:,:) :: tlak3d
!
      real(8) , dimension(ndpmax) :: de , dnsty , tt
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
      subroutine allocate_lake
      implicit none
#ifdef MPP1
      allocate(dhlake1(nnsg,iy,jxp))
      allocate(idep2d(nnsg,iym1,jxp))
      allocate(eta2d(nnsg,iym1,jxp))
      allocate(hi2d(nnsg,iym1,jxp))
      allocate(aveice2d(nnsg,iym1,jxp))
      allocate(hsnow2d(nnsg,iym1,jxp))
      allocate(evl2d(nnsg,iym1,jxp))
      allocate(tlak3d(ndpmax,nnsg,iym1,jxp))
#else
      allocate(dhlake1(nnsg,iy,jx))
      allocate(idep2d(nnsg,iym1,jx))
      allocate(eta2d(nnsg,iym1,jx))
      allocate(hi2d(nnsg,iym1,jx))
      allocate(aveice2d(nnsg,iym1,jx))
      allocate(hsnow2d(nnsg,iym1,jx))
      allocate(evl2d(nnsg,iym1,jx))
      allocate(tlak3d(ndpmax,nnsg,iym1,jx))
#endif
      dhlake1 = 0.0D0
      idep2d = 0
      eta2d = 0.0D0
      hi2d = 0.0D0
      aveice2d = 0.0D0
      hsnow2d = 0.0D0
      evl2d = 0.0D0
      tlak3d = 0.0D0
      end subroutine allocate_lake

      subroutine initlake
#ifdef MPP1
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

      hi2d = 0.01D0
      aveice2d = 0.0D0
      hsnow2d = 0.0D0
      eta2d = 0.5D0
      tlak3d  = 6.0D0
      idep2d = 0

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
            if ( (satbrt1(n,i,j).gt.13.9.and.satbrt1(n,i,j)  &
               & .lt.14.1) .and. dhlake1(n,i,j).gt.1.0) then
              idep2d(n,i,j) = int(max(2.D0,min(dhlake1(n,i,j), &
                                   dble(ndpmax)))/dz)
              if (idep2d(n,i,j).lt.50) then
                eta2d(n,i,j) = .7
              else if (idep2d(n,i,j).gt.100) then
                eta2d(n,i,j) = .3
              else
                eta2d(n,i,j) = .5
              end if
            else
              idep2d(n,i,j) = 0
              eta2d(n,i,j) = 0.5
            end if
          end do
        end do
      end do
#ifdef MPP1
      call mpi_gather(idep2d,   nnsg*iym1*jxp,mpi_integer, &
                    & idep2d_io,nnsg*iym1*jxp,mpi_integer, &
                    & 0, mpi_comm_world,ierr)
#endif
      end subroutine initlake
!
      subroutine lakedrv(jslc)
      implicit none
!
      integer , intent(in) :: jslc
!
      real(8) :: flw , fsw , hlat , hsen , prec , &
               & ql , tgl , tl , vl , zl
      integer :: i , n
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( idep2d(n,i,jslc).gt.1 ) then
            tl = ts1d(n,i)
            vl = sqrt(us1d(i)**2.0D0+vs1d(i)**2.0D0)
            zl = z1d(n,i)
            ql = qs1d(n,i)
            fsw = fsw1d(i)
            flw = -1.*flw1d(i)
            prec = prca2d(i,jslc)      !  units of prec = mm
            hsen = -1.0D0*sent1d(n,i)
            hlat = -1.0D0*evpr1d(n,i)

            call lake( dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat, &
                    &  tgl,prec,idep2d(n,i,jslc),eta2d(n,i,jslc),   &
                    &  hi2d(n,i,jslc),aveice2d(n,i,jslc),  &
                    &  hsnow2d(n,i,jslc),evl2d(n,i,jslc),  &
                    &  tlak3d(:,n,i,jslc) )

!           Feed back ground temperature
            tg1d(n,i) = tgl
            tgb1d(n,i) = tgl

            if ( aveice2d(n,i,jslc).le.10.0D0 ) then
              ldoc1d(n,i) = 0.0D0
              sice1d(n,i) = 0.0D0
              scv1d(n,i) = 0.0D0
              sag1d(n,i) = 0.0D0
            else
              ldoc1d(n,i) = 2.0D0
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
                    & prec,ndpt,eta,hi,aveice,hsnow,evl,tprof)
 
      implicit none
!
      real(8) :: dtlake , evl , aveice , hlat , hsen , hsnow , flw , &
               & prec , ql , fsw , tl , tgl , vl , zl , eta , hi
      real(8) , dimension(ndpmax) :: tprof
      integer :: ndpt
      intent (in) hlat , hsen , ql , tl , vl , zl
      intent (in) ndpt , eta
      intent (out) tgl
      intent (inout) evl , aveice , hsnow
      intent (inout) tprof
!
! Local variables
!
      integer :: k
      real(8) :: ea , hs , ld , lu , qe , qh , tac , tk , u2
!
!***  dtlake:  time step in seconds
!***  zo:      surface roughness length
!
      real(8) , parameter :: zo = 0.001D0
      real(8) , parameter :: z2 = 2.0D0
      real(8) , parameter :: tcutoff = -0.001D0
      logical , parameter :: lfreeze = .true.
      integer , parameter :: kmin = 1
!
!     interpolate winds at z1 m to 2m via log wind profile
      u2 = vl*log(z2/zo)/log(zl/zo)
      if ( u2.lt.0.5D0 ) u2 = 0.5D0
 
!******    depth: 1-m slices of lake depth

      tac = tl - tzero
      tk = tzero + tprof(1)
      lu = -0.97D0*sigm*tk**4.0D0
      ld = flw - lu
      qe = hlat*wlhv
      qh = hsen
 
!     ******    Check if conditions exist for lake ice
      if ( (aveice.eq.0.0D0) .and. (tprof(1).gt.tcutoff) ) then
 
!       ******    Calculate eddy diffusivities
        call eddy(ndpt,dtlake,u2,tprof)
        do k = 1 , ndpmax
          if ((tprof(k)/=tprof(k)) .or. ((tprof(k)>0.0).eqv.(tprof(k)<=0.0))) then
            print *, 'EDDY: At k = ',k
            print *, 'NAN'
            call fatal(__FILE__,__LINE__,'NaN in lake')
          end if
        end do
 
!       ******    Lake temperature calc using sensible and latent heats
        call temp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
        do k = 1 , ndpmax
          if ((tprof(k)/=tprof(k)) .or. ((tprof(k)>0.0).eqv.(tprof(k)<=0.0))) then
            print *, 'TEMP: At k = ',k
            print *, 'NAN'
            call fatal(__FILE__,__LINE__,'NaN in lake')
          end if
        end do
 
!       ******    Convective mixer
        call mixer(kmin,ndpt,tprof)
        do k = 1 , ndpmax
          if ((tprof(k)/=tprof(k)) .or. ((tprof(k)>0.0).eqv.(tprof(k)<=0.0))) then
            print *, 'MIX: At k = ',k
            print *, 'NAN'
            call fatal(__FILE__,__LINE__,'NaN in lake')
          end if
        end do
 
      else
 
!       convert mixing ratio to air vapor pressure
        ea = ql*88.0D0/(ep2+0.378D0*ql)
 
        call ice(fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,tprof)
        if ( lfreeze ) tprof(1) = tt(1)
 
      end if
 
      tgl = tprof(1) + tzero
      evl = evl/3600.0D0          !  convert evl from mm/hr to mm/sec
      aveice = aveice*1000.0D0    !  convert ice  from m to mm
      hsnow = hsnow*100.0D0       !  convert snow from m depth to mm h20
 
      end subroutine lake
!
!-----------------------------------------------------------------------
!
      subroutine eddy(ndpt,dtlake,u2,tprof)
 
! Computes density and eddy diffusivity
 
      implicit none
!
      integer , intent (in) :: ndpt
      real(8) , intent (in) :: dtlake , u2
      real(8) , dimension(ndpmax) , intent (in) :: tprof
!
      real(8) :: demax , demin , dpdz , ks , n2 , po , rad , ri , ws , z
      integer :: k
!
      demin = hdmw
      demax = .50D0*dz**2.0D0/dtlake
      demax = .99D0*demax

      do k = 1 , ndpt
        dnsty(k) = 1000.0D0*(1.0D0-1.9549D-05 * &
                      (abs((tprof(k)+tzero)-277.0D0))**1.68D0)
      end do
 
! compute eddy diffusion profile
!
!     n2    Brunt-Vaisala frequency squared
!     ri    gradient Richardson number
!     demin molecular diffusion of water
 
! Decay constant of shear velocity: Shouldn't it be function of latitude?
      ks = 0.745D0*u2**(-1.84D0)
! Surface shear velocity
      ws = 0.0012D0*u2
! Prandtl number
      po = 1.0D0
 
      do k = 1 , ndpt - 1
        dpdz = (dnsty(k+1)-dnsty(k))/dz
        n2 = (dpdz/dnsty(k))*gti
        z = surf + dble(k-1)*dz
!
! This line has problems:
!
!       rad = 1.0D0+40.0D0*n2*((vonkar*z)/(ws*exp(-ks*z)))**2.0D0
!
! I am proposing to change the previous to the below (change sign in
! the exponential).
! But I have no reference for this. Where this formula is taken from?
! Ciao!
        rad = 1.0D0+40.0D0*n2*((vonkar*z)/(ws*exp(ks*z)))**2.0D0
        ri = (-1.0D0+sqrt(rad))/20.0D0
        de(k) = demin + vonkar*ws*z*po*exp(-ks*z)/(1.0D0+37.0D0*ri**2.0D0)
        if ( de(k).lt.demin ) de(k) = demin
        if ( de(k).gt.demax ) de(k) = demax
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
 
      dt1 = (fsw*(1.0D0-exp(-eta*surf))+(flw+qe+qh)) / &
              (surf*dnsty(1)*cpw)
      dt2 = -de(1)*(tprof(1)-tprof(2))/surf
      tt(1) = tt(1) + (dt1+dt2)*dtlake
 
      do k = 2 , ndpt - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        dt1 = fsw*(exp(-eta*top)-exp(-eta*bot))/(dz*dnsty(k)*cpw)
        dt2 = (de(k-1)*(tprof(k-1)-tprof(k))    -    &
               de(k)  *(tprof(k)  -tprof(k+1))) / dz
        tt(k) = tt(k) + (dt1+dt2)*dtlake
      end do
 
      top = (surf+(ndpt-2)*dz)
      dt1 = fsw*exp(-eta*top)/(dz*dnsty(ndpt)*cpw)
      dt2 = de(ndpt-1)*(tprof(ndpt-1)-tprof(ndpt))/dz
      tt(ndpt) = tt(ndpt) + (dt1+dt2)*dtlake
 
      do k = 1 , ndpt
        tprof(k) = tt(k)
        dnsty(k) = 1000.0D0*(1.0D0-1.9549D-05 * &
                   (abs((tprof(k)+tzero)-277.0D0))**1.68D0)
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
        avet = 0.0D0
        avev = 0.0D0
 
        if ( dnsty(k).gt.dnsty(k+1) ) then
 
          do k2 = kmin , k + 1
            if ( k2.eq.1 ) then
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
            dnsty(k2) = 1000.0D0*(1.0D0-1.9549D-05 * &
                        (abs((tav+tzero)-277.0D0))**1.68D0)
          end do
        end if
 
      end do ! K loop
 
      tprof(kmin:ndpt) = tt(kmin:ndpt)
 
      end subroutine mixer
!
!-----------------------------------------------------------------------
!
      subroutine ice(fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,tprof)

      implicit none
      real(8) :: ea , evl , hi , aveice , hs , fsw , &
                 ld , prec , tac , u2
      real(8) , dimension(ndpmax) :: tprof
      intent (in) ea , ld , prec , tac , u2
      intent (out) evl
      intent (inout) hi , aveice , hs , fsw , tprof
!
      real(8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               & t2 , tf , theta , rho
      integer :: nits
!
      real(8) , parameter :: isurf = 0.6D0
      real(8) , parameter :: lami1 = 1.5D0
      real(8) , parameter :: lami2 = 20.0D0
      real(8) , parameter :: lams1 = 6.0D0
      real(8) , parameter :: lams2 = 20.0D0
      real(8) , parameter :: ki = 2.3D0
      real(8) , parameter :: ks = 0.31D0
      real(8) , parameter :: atm = 950.0D0
      real(8) , parameter :: qw = 1.389D0
      real(8) , parameter :: li = 334.0D03
      real(8) , parameter :: cd = 0.001D0
      real(8) , parameter :: sec = 3600.0D0
!
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           
!**********************************************************************
 
      if ( (tac.le.0.0D0) .and. (aveice.gt.0.0D0) ) &
        hs = hs + prec*10.0D0/1000.0D0  ! convert prec(mm) to depth(m)
 
      t0 = tprof(1)
      tf = 0.0D0
      rho = rhoh2o/1000.0D0
 
      khat = (ki*hs+ks*hi)/(ki*ks)
      theta = cpd*rho*cd*u2
      psi = wlhv*rho*cd*u2*ep2/atm
      evl = 100.0D0*psi*(eomb(t0)-ea)/(wlhv*rho)
      qpen = fsw*0.7D0*((1.0D0-exp(-lams1*hs))/(ks*lams1) +            &
                        (exp(-lams1*hs))*(1.0D0-exp(-lami1*hi)) /      &
                        (ki*lami1))+fsw*0.3D0*((1.0D0-exp(-lams2)) /   &
                        (ks*lams2)+(-lams2*hs)*(1.0D0-exp(-lami2*hi))/ &
                        (ki*lami2))
      fsw = fsw - qpen
 
      nits = 0
      t1 = -50.0D0
      f0 = f(t0)
      f1 = f(t1)
      do
        nits = nits + 1
        t2 = t1 - (t1-t0)*f1/(f1-f0)
        if ( abs((t2-t1)/t1).ge.0.001D0 ) then
          t0 = t1
          t1 = t2
          f0 = f1
          f1 = f(t1)
          cycle
        end if
 
        t0 = t2
        if ( t0.ge.tf ) then
 
          if ( hs.gt.0.0D0 ) then
            ds = sec*                                        &
               & ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea)+ &
               &  theta*(tf-tac)-fsw)-1.0D0/khat*(t0-tf+qpen))/(rhos*li)
            if ( ds.gt.0.0D0 ) ds = 0.0D0
            hs = hs + ds
            if ( hs.lt.0.0D0 ) then
              hs = 0.0D0
              tprof(1) = (aveice*t0+(isurf-aveice)*tprof(2))/isurf
            end if
          end if
          if ( (hs.eq.0.0D0) .and. (aveice.gt.0.0D0) ) then
            di = sec*                                        &
              & ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea) + &
                 theta*(tf-tac)-fsw)-1.0D0/khat*(t0-tf+qpen))/(rhoi*li)
            if ( di.gt.0.0D0 ) di = 0.0D0
            hi = hi + di
          end if
 
        else if ( t0.lt.tf ) then
 
          q0 = -ld + 0.97D0*sigm*t4(t0) + psi*(eomb(t0)-ea)             &
             & + theta*(t0-tac) - fsw
          qpen = fsw*0.7D0*(1.0D0-exp(-(lams1*hs+lami1*hi))) +          &
               & fsw*0.3D0*(1.0D0-exp(-(lams2*hs+lami2*hi)))
          di = sec*(q0-qw-qpen)/(rhoi*li)
 
          hi = hi + di
        end if
 
        if ( hi.le.0.01D0 ) then
          hi = 0.01D0
          aveice = 0.0D0
          hs = 0.0D0
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
        t4 = (x+tzero)**4.0D0
      end function t4
      ! Computes air vapor pressure as a function of temp (in K)
      function tr1(x)
        implicit none
        real(8) :: tr1
        real(8) , intent(in) :: x
        tr1 = 1.0D0 - (tboil/(x+tzero))
      end function tr1
      function eomb(x)
        implicit none
        real(8) :: eomb
        real(8) , intent(in) :: x
        eomb = stdpmb*exp(13.3185D0*tr1(x)-1.976D0*tr1(x)**2.D0   &
           &   -0.6445D0*tr1(x)**3.D0- 0.1299D0*tr1(x)**4.D0)
       end function eomb
      function f(x)
        implicit none
        real(8) :: f
        real(8) , intent(in) :: x
!     ****** g. bates changed air to tac, qpen1 to qpen (4/92)
!       f = (-ld+0.97D0*sigm*t4(x)+psi*(eomb(x)-ea)+theta
!              *(x-air)-fsw)-1.0D0/khat*(qpen1+tf-x)
        f = (-ld+0.97D0*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-tac)-fsw)  &
            - 1.0D0/khat*(qpen+tf-x)
      end function f
 
      end subroutine ice
!
!-----------------------------------------------------------------------
!
      subroutine outlake

#ifdef MPP1
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
      call mpi_gather(tlak3d,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                    & tlak3d_io,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
          do i = 2 , iym1
            do n = 1 , nnsg
              if ( idep2d_io(n,i,j).gt.1 ) then
                write(58) dayl, i, j, idep2d_io(n,i,j),  &
                        evl2d_io(n,i,j),hi2d_io(n,i,j), &
                        aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                        (tlak3d_io(k,n,i,j),k=1,idep2d_io(n,i,j))  
              end if
            end do
          end do
        end do
      end if

#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1 
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( idep2d(n,i,j).gt.1 ) then
              write(58) dayl, i, j, idep2d(n,i,j), evl2d(n,i,j), &
                   & hi2d(n,i,j), aveice2d(n,i,j), hsnow2d(n,i,j), &
                   & (tlak3d(k,n,i,j),k=1,idep2d(n,i,j))  
            end if
          end do
        end do
      end do
#endif

      end subroutine outlake
!
!-----------------------------------------------------------------------
!
      subroutine lakesav0_o

#ifdef MPP1
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
      call mpi_gather(tlak3d,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                    & tlak3d_io,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                    & 0, mpi_comm_world,ierr)
#endif

      end subroutine lakesav0_o
!
!-----------------------------------------------------------------------
!
      subroutine lakesav_o(iutl)

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
      write (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      write (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      numpts = 0
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( idep2d_io(n,i,j).gt.1 ) then
              numpts = numpts+1
            end if
          end do
        end do
      end do
      write (iutl) numpts
      print * , 'writing lake model restart file. numpts = ' , numpts
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( idep2d_io(n,i,j).gt.1 ) then
              write(iutl) n, i, j, idep2d_io(n,i,j), eta2d_io(n,i,j), &
                   & hi2d_io(n,i,j), &
                   & aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                   & (tlak3d_io(k,n,i,j),k=1,idep2d_io(n,i,j))  
            end if
          end do
        end do
      end do

#else
#ifdef BAND
      write (iutl) (((idep2d(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      write (iutl) (((idep2d(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      numpts = 0
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( idep2d(n,i,j).gt.1 ) then
              numpts = numpts+1
            end if
          end do
        end do
      end do
      write (iutl) numpts
      print * , 'writing lake model restart file. numpts = ' , numpts
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1 
#endif
        do i = 2 , iym1
          do n = 1 , nnsg
            if ( idep2d(n,i,j).gt.1 ) then
              write(iutl) n, i, j, idep2d(n,i,j), eta2d(n,i,j), &
                   & hi2d(n,i,j), &
                   & aveice2d(n,i,j), hsnow2d(n,i,j), &
                   & (tlak3d(k,n,i,j),k=1,idep2d(n,i,j))  
            end if
          end do
        end do
      end do
#endif

      end subroutine lakesav_o
!
!-----------------------------------------------------------------------
!
      subroutine lakesav_i(iutl)

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
        read (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
        read (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif

        hi2d_io = 0.01
        aveice2d_io = 0.0
        hsnow2d_io = 0.0
        eta2d_io = 0.5
        tlak3d_io  = 6.

        read (iutl) numpts
        print * , 'reading lake model restart file. numpts = ' , numpts
        do l = 1, numpts
          read(iutl) n, i, j, idep2d_io(n,i,j), eta2d_io(n,i,j), &
                   & hi2d_io(n,i,j), &
                   & aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
                   & (tlak3d_io(k,n,i,j),k=1,idep2d_io(n,i,j))  
        end do
      end if
#else
#ifdef BAND
      read (iutl) (((idep2d(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
      read (iutl) (((idep2d(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
      read (iutl) numpts
      print * , 'reading lake model restart file. numpts = ' , numpts
      do l = 1, numpts
        read(iutl) n, i, j, idep2d(n,i,j), eta2d(n,i,j), &
                 & hi2d(n,i,j), &
                 & aveice2d(n,i,j), hsnow2d(n,i,j), &
                 & (tlak3d(k,n,i,j),k=1,idep2d(n,i,j))  
      end do
#endif

      end subroutine lakesav_i
!
!-----------------------------------------------------------------------
!
      subroutine lakesav0_i

#ifdef MPP1
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
      call mpi_scatter(idep2d_io,nnsg*iym1*jxp,mpi_integer, &
                     & idep2d,nnsg*iym1*jxp,mpi_integer, &
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
      call mpi_scatter(tlak3d_io,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                     & tlak3d,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                     & 0, mpi_comm_world,ierr)
#endif

      end subroutine lakesav0_i
!
      end module mod_lake
