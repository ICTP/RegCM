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
      use mod_bats
      use mod_runparams
      use mod_message
#ifdef MPP1
      use mod_mppio
#endif

      private

      integer :: lkpts

      integer , dimension(:,:,:) , allocatable :: lakemask

      integer , dimension(:) , allocatable :: ndpts
      integer , dimension(:) , allocatable :: ifreeze
      real(8) , dimension(:) , allocatable :: xhi
      real(8) , dimension(:) , allocatable :: xhii
      real(8) , dimension(:) , allocatable :: xhs
      real(8) , dimension(:) , allocatable :: xeta

      integer , parameter :: maxdep = 400
      real(8) , dimension(:,:) , allocatable :: xt

      public :: lakedrv , initlk
      public :: lakesavread , lakesavwrite , allocate_lake
#ifdef MPP1
      public :: mpilake
#endif
!
      contains
!
!-----------------------------------------------------------------------
!
#ifdef MPP1
      subroutine mpilake
        use mpi
        implicit none
        integer :: ierr
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_bcast(lkpts,1,mpi_integer,0,mpi_comm_world,ierr)
        if (myid /= 0) then
#ifdef BAND
          allocate(lakemask(nnsg,iym1,jx))
#else
          allocate(lakemask(nnsg,iym1,jxm1))
#endif
          call allocate_lake
        endif
        call mpi_bcast(ndpts,lkpts,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ifreeze,lkpts,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(xhi,lkpts,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(xhii,lkpts,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(xhs,lkpts,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(xeta,lkpts,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(xt,lkpts*maxdep,mpi_real8,0,mpi_comm_world,ierr)
#ifdef BAND
        call mpi_bcast(lakemask,nnsg*iym1*jx, &
                       mpi_integer,0,mpi_comm_world,ierr)
#else
        call mpi_bcast(lakemask,nnsg*iym1*jxm1, &
                       mpi_integer,0,mpi_comm_world,ierr)
#endif

      end subroutine mpilake
#endif

      subroutine allocate_lake
        implicit none
        allocate(ndpts(lkpts))
        allocate(ifreeze(lkpts))
        allocate(xeta(lkpts))
        allocate(xhi(lkpts))
        allocate(xhii(lkpts))
        allocate(xhs(lkpts))
        allocate(xt(lkpts,maxdep))
        ndpts = 0
        ifreeze = 0
        xeta = 0.0D0
        xhi = 0.0D0
        xhii = 0.0D0
        xhs = 0.0D0
        xt = 0.0D0
      end subroutine

      subroutine initlk
        implicit none
        integer :: ilak , i , j , n , ie , je
        real(8) :: dpt
!
#ifdef BAND
        allocate(lakemask(nnsg,iym1,jx))
        ie = iym1
        je = jx
#else
        allocate(lakemask(nnsg,iym1,jxm1))
        ie = iym1
        je = jxm1
#endif
#ifdef MPP1
        where (satbrt1_io(:,1:ie,1:je) > 13.9 .and. &
               satbrt1_io(:,1:ie,1:je) < 14.1)
#else
        where (satbrt1(:,1:ie,1:je) > 13.9 .and. &
               satbrt1(:,1:ie,1:je) < 14.1)
#endif
          lakemask = 1
        elsewhere
          lakemask = 0
        end where

        lkpts = 0
        do n = 1 , nnsg
          do i = 1 , ie
            do j = 1 , je
              if (lakemask(n,i,j) > 0) then
                lkpts = lkpts + 1
                lakemask(n,i,j) = lkpts
              end if
            end do
          end do
        end do

        if (lkpts == 0) then
          write (6,*) 'lakemod : No points for lake model.'
        else
          write (6,*) 'lakemod : ', lkpts, ' points for lake model.'
        end if

        call allocate_lake

        ilak = 1
        do n = 1 , nnsg
          do i = 1 , ie
            do j = 1 , je
              if (lakemask(n,i,j) > 0) then
#ifdef MPP1
                dpt = lkdpth_io(n,i,j)
#else
                dpt = lkdpth(n,i,j)
#endif
                if (dpt < 50.0) then
                  xeta(ilak) = 7.0
                else if (dpt > 100.0) then
                  xeta(ilak) = 3.0
                else
                  xeta(ilak) = 5.0
                end if
                ndpts(ilak) = min(int(dpt), maxdep)
                ifreeze(ilak) = 1
                xt(ilak,:) = 6.0
                xhi(ilak) = 0.01
                xhii(ilak) = 0.0
                xhs(ilak) = 0.0
                ilak = ilak + 1
              end if
            end do
          end do
        end do
 
      end subroutine initlk
!
      subroutine lakesavwrite(iut)
        implicit none
        integer, intent(in) :: iut
        write(iut) lkpts
        if (lkpts > 0) then
          write(iut) ndpts , ifreeze , xhi , xhii , xhs , xeta , xt
        end if
      end subroutine lakesavwrite
!
      subroutine lakesavread(iut)
        implicit none
        integer, intent(in) :: iut
        integer :: ilak , i , j , n , ie , je
        read(iut) ilak
!
#ifdef BAND
        allocate(lakemask(nnsg,iym1,jx))
        ie = iym1
        je = jx
#else
        allocate(lakemask(nnsg,iym1,jxm1))
        ie = iym1
        je = jxm1
#endif
#ifdef MPP1
        where (satbrt1_io(:,1:ie,1:je) > 13.9 .and. &
               satbrt1_io(:,1:ie,1:je) < 14.1)
#else
        where (satbrt1(:,1:ie,1:je) > 13.9 .and. &
               satbrt1(:,1:ie,1:je) < 14.1)
#endif
          lakemask = 1
        elsewhere
          lakemask = 0
        end where
        lkpts = 0
        do n = 1 , nnsg
          do i = 1 , ie
            do j = 1 , je
              if (lakemask(n,i,j) > 0) then
                lkpts = lkpts + 1
                lakemask(n,i,j) = lkpts
              end if
            end do
          end do
        end do

        if (ilak /= lkpts) then
          write (6,*) 'Error Reading SAV for lake model.'
          write (6,*) 'Number of points saved does not match calculated'
          write (6,*) 'Saved      : ilak'
          write (6,*) 'Calculated : lkpts'
          call fatal(__FILE__,__LINE__,'SAV LAKE ERROR')
        end if

        if (lkpts == 0) then
          write (6,*) 'No points for lake model.'
          return
        end if

        call allocate_lake
        read(iut) ndpts , ifreeze , xhi , xhii , xhs , xeta , xt

      end subroutine lakesavread
!
      subroutine lakedrv(jslc)
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      intent (in) jslc
!
! Local variables
!
      real(8) :: aveice , evl , flw , fsw , hlat , hsen , prec , &
               & ql , hsnow , tgl , tl , vl , zl
      integer :: ilake , i , n , jj
!

      if (lkpts == 0) return

      do i = 1 , iym1
        do n = 1 , nnsg
          jj = (jxp*myid) + jslc
          if ( lakemask(n,i,jj) > 0 ) then
            ilake = lakemask(n,i,jj)
            tl = ts1d(n,i)
            vl = dsqrt(us1d(i)**2+vs1d(i)**2)
            zl = z1d(n,i)
            ql = qs1d(n,i)
            fsw = fsw1d(i)
            flw = -1.*flw1d(i)
            prec = prca2d(i,jslc)      !  units of prec = mm
            hsen = -1.*sent1d(n,i)
            hlat = -1.*evpr1d(n,i)

            call lake(ilake,dtlake,tl,vl,zl,ql,fsw,flw,hsen,hlat, &
                    & tgl,evl,prec,aveice,hsnow)
 
            tg1d(n,i) = tgl
            tgb1d(n,i) = tgl
            if ( aveice.le.10. ) then
              ldoc1d(n,i) = 0.
              sice1d(n,i) = 0.
              scv1d(n,i) = 0.
              sag1d(n,i) = 0.
            else
              ldoc1d(n,i) = 2.
              sice1d(n,i) = aveice  !  units of ice = mm
              scv1d(n,i) = hsnow   !  units of snow = mm h2o
            end if
 
          end if
        end do
      end do
 
      end subroutine lakedrv
!
!-----------------------------------------------------------------------
!
      subroutine lake(ilak,dt,ta,ua,za,q,sw,lnet,hsen,hlat,ts,    &
                    & evap,prec,hice,hsnow)
 
      implicit none
!
! Dummy arguments
!
      real(8) :: dt , evap , hice , hlat , hsen , hsnow , lnet ,  &
               & prec , q , sw , ta , ts , ua , za
      integer :: ilak
      intent (in) hlat , hsen , ilak , q , ta , ua , za
      intent (out) ts , evap , hice , hsnow
!
! Local variables
!
      real(8) , dimension(maxdep) :: de , dnsty
      real(8) , dimension(maxdep,2) :: t
      integer :: ndep , nfrz , kmin
      real(8) :: ea , eta , hi , hs , ld , lu , qe , qh , tac , surf ,  &
              &  tcutoff , tk , u2
!
!***  dt:  time step in seconds
!***  dz:  vertical grid spacing in m
!***  zo:  surface roughness length
      real(8) , parameter :: dz = 1.0
      real(8) , parameter :: zo = 0.001
      real(8) , parameter :: z2 = 2.0
!
!     surf:surface thickness
!
      surf = 1.0
!
!     interpolate winds at z1 m to 2m via log wind profile
!
      u2 = ua*dlog(z2/zo)/dlog(za/zo)
 
!     Set current values for lake point

      ndep = ndpts(ilak)
      nfrz = ifreeze(ilak)

      hi = xhi(ilak)
      hice = xhii(ilak)
      hsnow = xhs(ilak)
      eta = xeta(ilak)

      t(:,1) = xt(ilak,:)
      t(:,2) = t(:,1)
 
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
        call eddy(dt,surf,dz,vonkar,u2,t,dnsty,de,ndep)
!       ******    Lake temperature calc using BATS sensible and latent
!       heats
        call temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,maxdep,ndep)
!       ******    Convective mixer
        kmin = 1
        call mixer(kmin,surf,dz,t,dnsty,maxdep,ndep)
      else
        call ice(sw,ld,tac,u2,ea,hs,hi,hice,evap,t,ndep,prec)
        if ( nfrz.eq.0 ) t(1,1) = t(1,2)
      end if
 
      ndpts(ilak) = ndep
      ifreeze(ilak) = nfrz
      xhi(ilak) = hi
      xhii(ilak) = hice
      xhs(ilak) = hsnow
      xeta(ilak) = eta
      xt(ilak,:) = t(:,1)

      ts = t(1,1) + tzero
      evap = evap/3600.          !  convert evap from mm/hr to mm/sec
      hice = hice*1000.          !  convert ice  from m to mm
      hsnow = hsnow*100.         !  convert snow from m depth to mm h2o
 
      end subroutine lake
!
!-----------------------------------------------------------------------
!
      subroutine eddy(dt,surf,dz,kv,u2,t,dnsty,de,ndep)
 
! Computes density, eddy diffusivity and variable time step
 
      implicit none
!
! PARAMETER definitions
!
!cc   dm=5.148e-04       ! value used in prev simulations
      real(8) , parameter :: dm = 1.38889E-07
!
! Dummy arguments
!
      integer :: ndep
      real(8) :: dt , dz , kv , surf , u2
      real(8) , dimension(ndep) :: de , dnsty
      real(8) , dimension(ndep,2) :: t
      intent (in) ndep , dt , dz , kv , surf , t
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
      do k = 1 , ndep
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
 
      do k = 1 , ndep - 1
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
      de(ndep) = 0.0
 
      end subroutine eddy
!
!-----------------------------------------------------------------------
!
      subroutine temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,maxdep,   &
                    & ndep)
!*****************BEGIN SUBROUTINE TEMP********************
!             COMPUTES TEMPERATURE PROFILE                *
!**********************************************************
      implicit none
!
! Dummy arguments
!
      integer :: maxdep , ndep
      real(8) :: dt , dz , eta , lnet , qe , qh , surf , sw
      real(8) , dimension(maxdep) :: de , dnsty
      real(8) , dimension(maxdep,2) :: t
      intent (in) de , maxdep , ndep , dt , dz , eta , lnet , qe , qh ,&
                & sw
      intent (inout) dnsty , surf , t
!
! Local variables
!
      real(8) :: bot , t1 , t2 , tdiff , top
      integer :: k
 
      surf = 1.0
 
!******    solve differential equations of heat transfer
      do k = 1 , ndep
        t(k,2) = t(k,1)
      end do
 
      k = 1
      t1 = sw*(1.-dexp(-eta*surf))/(surf*dnsty(k)*cpw) + (lnet+qe+qh)   &
         & /(surf*dnsty(k)*cpw)
      t2 = -de(k)*(t(k,1)-t(k+1,1))/surf
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      do k = 2 , ndep - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        t1 = sw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cpw)
        t2 = (de(k-1)*(t(k-1,1)-t(k,1))-de(k)*(t(k,1)-t(k+1,1)))/dz
        t(k,2) = t(k,2) + (t1+t2)*dt
      end do
 
      k = ndep
      top = (surf+(k-2)*dz)
      t1 = sw*dexp(-eta*top)/(dz*dnsty(k)*cpw)
      t2 = de(k-1)*(t(ndep-1,1)-t(ndep,1))/dz
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      tdiff = 0.
      do k = 1 , ndep
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
      subroutine mixer(kmin,surf,dz,t,dnsty,maxdep,ndep)
!
! Simulates convective mixing
!
      implicit none
!
! Dummy arguments
!
      integer :: maxdep , ndep , kmin
      real(8) :: dz , surf
      real(8) , dimension(maxdep) :: dnsty
      real(8) , dimension(maxdep,2) :: t
      intent (in) maxdep , ndep , dz , kmin , surf
      intent (inout) dnsty , t
!
! Local variables
!
      real(8) :: avet , avev , tav , tdiff , vol
      integer :: k , k2 , m
! 
      do k = kmin , ndep - 1
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
      do k = kmin , ndep
        tdiff = tdiff + t(k,2) - t(k,1)
        if ( k.eq.1 ) tdiff = tdiff*surf
        t(k,1) = t(k,2)
      end do
 
      end subroutine mixer
!
! Computes air vapor pressure as a function of temp (in K)
!
      function eomb(x) result(e)

      implicit none
!
! Dummy arguments
!
      real(8) :: x
      intent (in) x
      real(8) :: e
!
! Local variables
!
      real(8) :: tr1
!
      tr1 = 1.0 - (tboil/(x+tzero))
      e = stdpmb*dexp(13.3185*tr1-1.976*tr1**2-0.6445*tr1**3-       &
           & 0.1299*tr1**4)
      end function eomb

!	
!
!       function t4(x) result(t4) 
!       implicit none 
!       real(8), intent(in) :: x 
!       real(8) :: t4 

!        t4 = (x+tzero)**4
!       end function t4

!       function f(x) result(f)
!       real(8), intent(in) :: x 
!       real(8) :: f

 
!-----------------------------------------------------------------------
!
      subroutine ice(kd,ld,ta,u2,ea,hs,hi,hii,evap,t,ndep,precip)

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
      integer :: ndep
      real(8) :: ea , evap , hi , hii , hs , kd , ld , precip , ta , u2
      real(8) , dimension(ndep,2) :: t
      intent (in) ndep , ea , ld , precip , ta , u2
      intent (out) evap
      intent (inout) hi , hii , hs , kd , t
!
! Local variables
!
      real(8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               & t2 , tf , theta , x
      real(8) :: f , t4 ,tr1,eomb
      integer :: nits
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           *
!***********************************************************************
 
      t4(x) = (x+tzero)**4
 
      tr1(x) = 1.0 - (tboil/(x+tzero))
      eomb(x) = stdpmb*dexp(13.3185*tr1(x)-1.976*tr1(x)**2-0.6445* &
                tr1(x)**3-0.1299*tr1(x)**4)

!     ****** g. bates changed air to ta, qpen1 to qpen (4/92)
      f(x) = (-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-ta)-kd)    &
           & - 1./khat*(qpen+tf-x)
 
!CC   f(x)=(-ld+0.97*sigm*t4(x)+psi*(eomb(x)-ea)+theta
!CC   +      *(x-air)-kd)-1/khat*(qpen1+tf-x)
 
      if ( (ta.le.0.0) .and. (hii.gt.0.0) ) hs = hs + precip*10./1000.
                                     ! convert precip(mm) to ndep(m)
 
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
      end module mod_lake
