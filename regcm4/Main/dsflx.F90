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
 
      subroutine sfflux(ilg,il1,il2,jloop,luc,ivegcov,vegfrac,isoiltex, &
                      & ustarnd, z0,soilw,surfwd,roarow,trsize,rsfrow)
 
!  **********************************************************
!  *  dust emission scheme                             ******
!  *                                                   ******
!  * this scheme based on marticorena and bergametti,  ******
!  * 1995; gong et al.,(2003); alfaro et al.,(1997)    ******
!  *                                                   ******
!  * the modification coded by:                        ******
!  * ashraf s. zakey                                   ******
!  **********************************************************
 
      use mod_dynparam
      use mod_dust
      implicit none
!
! Dummy arguments
!
      integer :: il1 , il2 , ilg , jloop , luc
      integer , dimension(ilg) :: isoiltex , ivegcov
      real(8) , dimension(ilg) :: roarow , soilw , surfwd , vegfrac , &
                       &        z0,ustarnd
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(nbin,2) :: trsize
      intent (in) il1 , il2 , isoiltex , ivegcov , jloop , roarow ,     &
                & soilw , surfwd , vegfrac , z0, ustarnd
      intent (out) rsfrow
!
! Local variables
!
      integer :: i , ieff , ieffmax , n , ns
      real(8) , dimension(ilg) :: xclayrow , xroarow , xsoilw ,         &
                                & xsurfwd , xvegfrac , xz0, xustarnd
      real(8) , dimension(ilg,20) :: xfland
      integer , dimension(ilg) :: xisoiltex
      real(8) , dimension(ilg,nbin) :: xrsfrow
      real(8) , dimension(ilg,nats) :: xsand2row
      real(8) , dimension(ilg,nsoil,nats) :: xsrel2d
! 
      rsfrow = 0.0
!     effective emitter cell ( depending on ivegcov)
      xvegfrac = 0.
      xisoiltex = 0.
      xsoilw = 0.
      xsurfwd = 0.
      xz0 = 0.
      xclayrow = 0.
      xroarow = 0.
      xsrel2d = 0.
      xsand2row = 0.
      xustarnd=0.
      xfland = 0.0
      xrsfrow = 0.
 
      ieff = 0
      ieffmax = 0
      do i = il1 , il2
        if ( (ivegcov(i).eq.8 .or. ivegcov(i).eq.11) .and.              &
           & (isoiltex(i).gt.0 .and. isoiltex(i).le.12) ) then
          ieff = ieff + 1
          xvegfrac(ieff) = vegfrac(i)
          xisoiltex(ieff) = isoiltex(i)
          xsoilw(ieff) = soilw(i)
          xsurfwd(ieff) = surfwd(i)
          xz0(ieff) = z0(i)
          xroarow(ieff) = roarow(i)
          xustarnd(ieff) = ustarnd(i) 
!         soil parameters
!         ok if one texture per grid cell
          do ns = 1 , nsoil
            xsrel2d(ieff,ns,xisoiltex(ieff)) = srel2d(i,jloop,ns)
          end do
 
          xclayrow(ieff) = clayrow2(i,jloop)
          do n = 1 , nats
            xsand2row(ieff,n) = sand2row2(i,n,jloop)
          end do
          do n = 1 , 20
            xfland(ieff,n) = 1
          end do
        end if
      end do
 
      ieffmax = ieff
!      if (ieffmax>0. ) print*,&
!         & maxval(xustarnd)

      if ( ieffmax.gt.0 ) call dust_module(1,ieffmax,ilg,trsize,xsoilw, &
         & xvegfrac,xsurfwd,xfland,xclayrow,xsand2row,xroarow,xz0,      &
         & xsrel2d,xustarnd,xrsfrow,luc)
        
!        if (ieffmax>0. ) print*,'FLUX',maxval(xrsfrow)
!     put back the dust flux on the right grid
 
      ieff = 0
 
      do i = il1 , il2
        if ( (ivegcov(i).eq.8 .or. ivegcov(i).eq.11) .and.              &
           & (isoiltex(i).gt.0 .and. isoiltex(i).le.12) ) then
          ieff = ieff + 1
          do n = 1 , nbin
            rsfrow(i,n) = xrsfrow(ieff,n)
          end do
        end if
      end do
 
      end subroutine sfflux
! 
      subroutine dust_module(il1,il2,ilg,trsize,soilw,vegfrac,surfwd,   &
                           & fland,clayrow,sand2row,roarow,z0,srel,ustarnd,     &
                           & rsfrow,luc)
 
      use mod_dynparam
      use mod_dust
      use mod_message
      use mod_constants , only : vonkar
      use mod_aero_param , only : rhop
      implicit none
!
! Dummy arguments
!
      integer :: il1 , il2 , ilg , luc
      real(8) , dimension(ilg) :: clayrow , roarow , soilw , surfwd ,   &
                                & vegfrac , z0, ustarnd
      real(8) , dimension(ilg,20) :: fland
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(ilg,nats) :: sand2row
      real(8) , dimension(ilg,nsoil,nats) :: srel
      real(8) , dimension(nbin,2) :: trsize
      intent (in) clayrow , soilw , surfwd , z0,ustarnd
!
! Local variables
!
      real(8) , dimension(ilg) :: alamda , hc , rc , srl , wprim
      real(8) :: arc1 , arc2 , br , cly1 , cly2 , sigr , tempd ,        &
               & umin , ustarns , uth , utmin , x , xz , ym , z0s
      integer :: i
      real(8) , dimension(ilg,nats) :: ustar
      real(8) , dimension(ilg,nsoil) :: utheff
!
      data umin/15./
      data xz/0.25/ , br/202.0/ , ym/0.16/ , sigr/1.45/
      data z0s/3.E-3/ , x/10./

      do i = il1 , il2
 
        srl(i) = z0(i)*100.0
        rc(i) = 1.0
 
        if ( jfs.eq.0 ) then
 
! *****************************************************************
!         * raupach et al. (1993)                                     
!         ****
! *****************************************************************
          if ( vegfrac(i).lt.1.0 ) then
            alamda(i) = xz*(log(1.0-vegfrac(i)))*(-1.0)
            arc1 = sigr*ym*alamda(i)
            arc2 = br*ym*alamda(i)
            if ( arc1.le.1.0 .and. arc2.le.1.0 ) rc(i)                  &
               & = (sqrt(1.0-arc1)*sqrt(1.0+arc2))
          end if
 
        else if ( jfs.eq.1 ) then
! Marticorena et al., 1997: correction factor for non erodible elements
!  
          rc(i) = 1 - (dlog(0.5E-2/z0s)/(dlog(0.35*(x/z0s)**0.8)))
 
        end if
 
!       *************************************************************
!       threshold velocity correction for soil humidity hc
!***************************************************************
 
        if ( jsoilm.eq.0 ) then
 
          if ( soilw(i).lt.0.0 ) then
            write (aline,*) 'hc, rc =' , soilw(i) , ' less than zero'
            call say
            call fatal(__FILE__,__LINE__,'NEGATIVE SOILW')
          else if ( soilw(i).lt.0.03 ) then
            hc(i) = exp(22.7*soilw(i))
          else if ( soilw(i).ge.0.03 ) then
            hc(i) = exp(95.3*soilw(i)-2.029)
          else
            hc(i) = 1.0
          end if
 
        else if ( jsoilm.eq.1 ) then
 
          cly1 = clayrow(i)
          cly2 = cly1*cly1
          wprim(i) = 0.0014*cly2 + 0.17*cly1
          if ( soilw(i).lt.wprim(i) ) then
            hc(i) = sqrt(1.0+1.21*tempd**0.68)
          else
            hc(i) = 1.0
          end if
 
! no soil humidity correction facor if jsoilm > 1
        else
          hc(i)=1.0
        end if
 
! *****************************************************************
!       * total correction factor for both hc and rc                
!       ****
! *****************************************************************
        rc(i) = rc(i)/hc(i)
 
! *******************************************************************
!       *     computation of the wind friction velocity              
!       ***** *     accounting for the increase of the roughness length
!       ***** *     due to the saltation layer (gillette etal. jgr 103,
!       ***** *     no. d6, p6203-6209, 1998                           
!       *****
! *******************************************************************
!        ustarns = (vonkar*100.*surfwd(i))/(log(1000./srl(i)))

        ustarns = ustarnd(i)*100 !cm.s-1
        utmin = (umin/(100.*vonkar*rc(i)))*log(1000./srl(i))
 
! *******************************************************************
!       *     vonkar=karman constant, and 1000 cm =10 m:            
!       ***** *     the height of wind defined level. umin: 21 cm/s:   
!       ***** *    the minimal threshold wind friction velocity.       
!       *****
 
! *******************************************************************
        if ( surfwd(i).ge.utmin ) then
          ustar(i,:) = ustarns + 0.3*(surfwd(i)-utmin)*(surfwd(i)-utmin)
        else
          ustar(i,:) = ustarns
        end if
 
      end do       ! end i loop
 
      call uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop,dp)
 
      call emission(ilg,il1,il2,luc,rhop,nsoil,nbin,nats,fland,uth,     &
                  & roarow,dp,rc,utheff,ustar,srel,rsfrow,trsize,       &
                  & sand2row,vegfrac)
 
      end subroutine dust_module
! 
      subroutine uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop,dp)
      implicit none
!
! Dummy arguments
!
      integer :: il1 , il2 , ilg , nsoil , ust
      real(8) :: rhop
      real(8) , dimension(nsoil) :: dp
      real(8) , dimension(ilg) :: roarow
      real(8) , dimension(ilg,nsoil) :: utheff
      intent (in) il1 , il2 , ilg , nsoil , ust
      intent (out) utheff
!
! Local variables
!
      integer :: i , j
      real(8) , external :: ustart0, ustart01
!
      do i = 1 , nsoil
        do j = il1 , il2
          if ( ust.eq.0 ) utheff(j,i) = ustart0(rhop,dp(i),roarow(j))
          if ( ust.eq.1 ) utheff(j,i) = ustart01(rhop,dp(i),roarow(j))
        end do
      end do
 
      end subroutine uthefft
! 
      function ustart01(rhop,dum,rhair)
      use mod_constants , only : gti
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: a2 = 0.129 , c1 = 0.006 , c2 = 1.928 ,     &
                           & c3 = 0.0858 , c4 = -0.0617 , c5 = 2.5 ,    &
                           & y1 = 1331.647 , y2 = 1.561228 ,            &
                           & y3 = 0.38194
!
! Dummy arguments
!
      real(8) :: dum , rhair , rhop
      real(8) :: ustart01
      intent (in) dum , rhair , rhop
!
! Local variables
!
      real(8) :: dm , rep , term , term1 , term2
      real(8) , external :: cvmgt
!
!     *****************************************************************
!     * calculate of ustar01(d) using iversen and white (1982)     ****
!     * for smoth surface:                                         ****
!     * coded by :                                                 ****
!     * ashraf s. zakey, 2003                                      ****
!     * dum    : particle diameter [um]                            ****
!     * ustar0 : threshold frication velocity [m/s]                ****
!     *****************************************************************
 
      dm = dum  !* 1.0e-4      ! cm
      rep = y1*(dm**y2) + y3
      term1 = sqrt(1.0+(c1/(rhop*gti*0.1*(dm**c5))))
      term2 = sqrt(rhop*gti*100.0*dm/rhair)
      term = term1*term2
      ustart01 = cvmgt(a2*term*(1.0-c3*exp(c4*(rep-10.0))),             &
               & a2*term/sqrt(c2*(rep**0.092)-1.0),rep.gt.10.0)
 
      end function ustart01
! 
      function ustart0(rhop,dum,rhoa)
      use mod_constants , only : gti
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: agamma = 3.0E-4 , f = 0.0123
!
! Dummy arguments
!
      real(8) :: dum , rhoa , rhop
      real(8) :: ustart0
      intent (in) dum , rhoa , rhop
!
! Local variables
!
      real(8) :: dm , sigma
!
!     *****************************************************************
!     *                                                            ****
!     * modified by a.s.zakey, nov.2003                            ****
!     * y. shao, 13 june 2000                                      ****
!     * calculate ustar0(d) using shao and lu (2000) for uncovered ****
!     * dry surface                                                ****
!     * dum:    particle diameter                   [um]           ****
!     * ustar0: threshold friction velocity       [cm/s]           ****
!     *****************************************************************
                                   ! a constant
 
      sigma = rhop/rhoa        ! gravity parameter    [m s^-2]
      dm = dum*1.0E-2    !* 1.0e-6
      ustart0 = f*(sigma*gti*dm+agamma/(rhoa*dm))
      ustart0 = sqrt(ustart0)
      ustart0 = ustart0*100.0
      end function ustart0
! 
      function cvmgt(val1,val2,cond)
      implicit none
!
! Dummy arguments
!
      logical :: cond
      real(8) :: val1 , val2
      real(8) :: cvmgt
      intent (in) cond , val1 , val2
!
      if ( cond ) then
        cvmgt = val1
      else
        cvmgt = val2
      end if
!
      end function cvmgt
!
      subroutine emission(ilg,il1,il2,luc,rhop,nsoil,nbin,nats,fland,   &
                        & uth,roarow,dp,rc,utheff,ustar,srel,rsfrow,    &
                        & trsize,sand2row,vegfrac)
 
      use mod_constants , only : rgti , mathpi
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: isize = 12
!
! Dummy arguments
!
      integer :: il1 , il2 , ilg , luc , nats , nbin , nsoil
      real(8) :: rhop , uth
      real(8) , dimension(nsoil) :: dp
      real(8) , dimension(ilg,luc) :: fland
      real(8) , dimension(ilg) :: rc , roarow , vegfrac
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(ilg,nats) :: sand2row , ustar
      real(8) , dimension(ilg,nsoil,nats) :: srel
      real(8) , dimension(nbin,2) :: trsize
      real(8) , dimension(ilg,nsoil) :: utheff
      intent (in) dp , fland , il1 , il2 , ilg , luc , nats , nbin ,    &
                & nsoil , rc , rhop , roarow , sand2row , srel ,        &
                & trsize , ustar , utheff , vegfrac
      intent (inout) rsfrow , uth
!
! Local variables
!
      real(8) :: aeffect , alogdi , amean1 , amean2 , amean3 , asigma1 ,&
               & asigma2 , asigma3 , beffect , beta , const , d1 , d2 , &
               & d3 , dec , e1 , e2 , e3 , ec , f , fdp1 , fdp2 , p1 ,  &
               & p2 , p3 , rwi , sigma1 , sigma2 , sigma3 , totv1 ,     &
               & totv2 , totv3
      real(8) , dimension(2,isize) :: aerosize
      real(8) , dimension(isize) :: frac1 , frac2 , frac3
      real(8) , dimension(ilg) :: fsoil , fsoil1 , fsoil2 , fsoil3
      integer :: i , j , k , n
      real(8) , dimension(ilg,isize) :: rsfrowsub
!
      data const/2.3/ , beta/16300./
!     alfaro 's values
      data e1/3.61/ , e2/3.52/ , e3/3.46/
      data d1/1.5/ , d2/6.7/ , d3/14.2/
      data sigma1/1.7/ , sigma2/1.6/ , sigma3/1.5/
 
!     emissions  bins (sub-bins)
      data aerosize/1.0E-08 , 2.0E-08 , 2.0E-08 , 4.0E-08 , 4.0E-08 ,   &
         & 8.0E-08 , 8.0E-08 , 1.6E-07 , 1.6E-07 , 3.2E-07 , 3.2E-07 ,  &
         & 6.4E-07 , 6.4E-07 , 1.28E-06 , 1.28E-06 , 2.56E-06 ,         &
         & 2.56E-06 , 5.12E-06 , 5.12E-06 , 10.4E-06 , 10.24E-06 ,      &
         & 20.48E-06 , 20.48E-06 , 40.6E-06/
 
      p1 = 0.0
      p2 = 0.0
      p3 = 0.0
      fsoil(:) = 0.
      fsoil1(:) = 0.
      fsoil2(:) = 0.
      fsoil3(:) = 0.
 
      do i = 1 , nats
        do j = 1 , nsoil
          do k = il1 , il2
 
            if ( rc(k).gt.0.0 .and. ustar(k,i).ne.0. ) then
              uth = utheff(k,j)/(rc(k)*ustar(k,i))
 
              if ( uth.le.1.0 ) then
 
                fdp1 = ustar(k,i)**3*(1.0-uth*uth)
                fdp2 = (1.0+uth)*const*(1.E-5)*roarow(k)*rgti
 
                if ( fdp2.le.0.0 ) fdp2 = 0.
 
                f = 0.0D0
                aeffect = (1-f)*(1-vegfrac(k))
                beffect = 0.01*fland(k,i)*sand2row(k,i)
 
!                fsoil(k) = srel(k,j,i)*fdp1*fdp2*aeffect*beffect
! FAB 
                fsoil(k) = srel(k,j,i)*fdp1*fdp2 
 
!               size-distributed kinetic energy flux
                dec = fsoil(k)*beta
!               individual kinetic energy for an aggregate of size dp (
!               g cm2 s-2) cf alfaro (dp) is in cm
!               ec=(pi/3.)*1.e-1*rhop*(dp(j)**3.0)*(ustar(k,i)**2.0)
                ec = (mathpi/12)*rhop*1E-3*(dp(j)**3.0)*                &
                    & (20*ustar(k,i))**2.0
 
                if ( ec.gt.e1 ) then
                  p1 = (ec-e1)/(ec-e3)
                  p2 = (1-p1)*(ec-e2)/(ec-e3)
                  p3 = 1 - p1 - p2
                else if ( ec.gt.e2 .and. ec.le.e1 ) then
                  p1 = 0.
                  p2 = (ec-e2)/(ec-e3)
                  p3 = 1 - p2
                else if ( ec.gt.e3 .and. ec.le.e2 ) then
                  p1 = 0.
                  p2 = 0.
                  p3 = 1.
                else if ( ec.le.e3 ) then
                  p1 = 0.
                  p2 = 0.
                  p3 = 0.
                else
                end if
 
                fsoil1(k) = fsoil1(k) + 1.E-2*p1*(dec/e1)*(mathpi/6.)   &
                          & *rhop*((d1*1.E-04)**3.)
                fsoil2(k) = fsoil2(k) + 1.E-2*p2*(dec/e2)*(mathpi/6.)   &
                          & *rhop*((d2*1.E-04)**3.)
                fsoil3(k) = fsoil3(k) + 1.E-2*p3*(dec/e3)*(mathpi/6.)   &
                          & *rhop*((d3*1.E-04)**3.)
              end if
            end if
          end do
        end do
      end do

      totv1 = 0.0
      totv2 = 0.0
      totv3 = 0.0
 
      do n = 1 , isize
        rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.E6
        alogdi = log10(rwi)
        amean1 = log10(d1)
        amean2 = log10(d2)
        amean3 = log10(d3)
 
        asigma1 = log10(sigma1)
        asigma2 = log10(sigma2)
        asigma3 = log10(sigma3)
 
        frac1(n) = exp(-(alogdi-amean1)**2./(2*asigma1**2))
        frac2(n) = exp(-(alogdi-amean2)**2./(2*asigma2**2))
        frac3(n) = exp(-(alogdi-amean3)**2./(2*asigma3**2))
 
        totv1 = totv1 + frac1(n)
        totv2 = totv2 + frac2(n)
        totv3 = totv3 + frac3(n)
      end do


 
       do n = 1 , isize
        frac1(n) = frac1(n)/totv1
        frac2(n) = frac2(n)/totv2
        frac3(n) = frac3(n)/totv3
        if ( frac1(n).lt.1.E-9 ) frac1(n) = 0.0
        if ( frac2(n).lt.1.E-9 ) frac2(n) = 0.0
        if ( frac3(n).lt.1.E-9 ) frac3(n) = 0.0
      end do

      do n = 1 , isize
        do i = il1 , il2
 
!         discretisation of the modal emission in isize emission sub bin
          rsfrowsub(i,n) = fsoil1(i)*frac1(n) + fsoil2(i)*frac2(n)      &
                         & + fsoil3(i)*frac3(n)
 
!         and in tranport bins (nbin)
          rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.E6

          
          do k = 1 , nbin
            if ( rwi.ge.trsize(k,1) .and. rwi.lt.trsize(k,2) )          &
               & rsfrow(i,k) = rsfrow(i,k) + rsfrowsub(i,n)
          end do
        end do
      end do

 
      end subroutine emission
