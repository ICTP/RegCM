C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine sfflux(ilg,il1,il2,jloop,luc,ivegcov,
     &                  vegfrac, isoiltex, ustarnd, z0,
     &                  soilw, surfwd, roarow, trsize, rsfrow)

c  **********************************************************
c  *  dust emission scheme                             ******
c  *                                                   ******
c  * this scheme based on marticorena and bergametti,  ******
c  * 1995; gong et al.,(2003); alfaro et al.,(1997)    ******
c  *                                                   ******
c  * the modification coded by:                        ******
c  * ashraf s. zakey                                   ******
c  **********************************************************

      implicit none

#     include "parameter.inc"
#     include "include/dust.cb"

      integer il1,il2,ilg,jloop
      integer n, ns
      integer luc
      integer i,ieff, ieffmax
c  
      real(kind=8)  rsfrow(ilg,nbin)
      real(kind=8)  soilw(ilg),ustarnd(ilg)
      real(kind=8)  surfwd(ilg),z0(ilg) 
      integer ivegcov(ilg),isoiltex(ilg) 
      real(kind=8)  vegfrac(ilg)
      real(kind=8)  roarow(ilg),sflt(ilg)       
      real(kind=8)  trsize(nbin,2)

c
c        effective work array  ( maybe think about allocatbale array)

      real(kind=8)  xrsfrow(ilg,nbin)
      real(kind=8)  xsoilw(ilg)
      real(kind=8)  xsurfwd(ilg) , xz0(ilg)
      real(kind=8)  xvegfrac(ilg), xfland(ilg,20)
      real(kind=8)  xroarow(ilg),xsflt(ilg)       
      real(kind=8)  xsrel2d(ilg,nsoil,nats)
      real(kind=8)  xustarnd(ilg) 
      real(kind=8) xsandrow(ilg),  xclayrow(ilg)
      real(kind=8) xsand2row(ilg,nats), xclay2row(ilg,nats)

      integer xivegcov(ilg), xisoiltex(ilg), xirb(ilg)

      rsfrow = 0.0

c  effective emitter cell ( depending on ivegcov) 
      xivegcov = 0.
      xvegfrac = 0. 
      xisoiltex= 0. 
      xsoilw   = 0.
      xsurfwd  = 0.
      xz0      = 0.
      xclayrow = 0.
      xsandrow = 0.
      xirb = 0
      xustarnd = 0. 
      xsflt = 0.
      xroarow =0.
      xsrel2d =0.
      xclay2row =0.
      xsand2row =0.
      xfland =0.0
      xrsfrow =0.

      ieff = 0  
      ieffmax=0
      do i=il1, il2
        if ((ivegcov(i) .eq.8 .or.  ivegcov(i) .eq.11) .and.
     &     (isoiltex(i).gt.0 .and. isoiltex(i).le.12) )then
          ieff = ieff + 1 
          xivegcov(ieff) =  ivegcov(i)
          xvegfrac(ieff) =  vegfrac(i) 
          xisoiltex(ieff) = isoiltex(i) 
          xsoilw(ieff)  =  soilw(i)
          xsurfwd(ieff) = surfwd(i)
          xz0(ieff)    = z0(i)
          xroarow(ieff) = roarow(i)
          xustarnd(ieff)   =  ustarnd(i)

c soil parameters
c ok if one texture per grid cell
c (remember srel2d,4 =1 because of memory problem for big domains)
          do ns=1,nsoil
            xsrel2d(ieff,ns,xisoiltex(ieff) )=srel2d(i,jloop,ns)
          end do

          xsandrow(ieff)=sandrow2(i,jloop) 
          xclayrow(ieff)=clayrow2(i,jloop) 
          do n=1,nats
            xsand2row(ieff,n)=sand2row2(i,n,jloop)
            xclay2row(ieff,n)=clay2row2(i,n,jloop)
          end do
c  trivial !
          xirb(ieff) = 1
          do n=1,20
            xfland(ieff,n)  = 1
          end do
        end if
      end do 

      ieffmax = ieff 

      if (ieffmax .gt. 0 ) then 
        call dust_module(1,ieffmax,ilg, trsize,
     $                   xsoilw, xvegfrac,xsurfwd, xfland,
     $                   xclayrow,xsand2row,
     $                   xroarow, xz0, xsrel2d,
     $                   xrsfrow,luc)
      end if

c  put back the dust flux on the right grid

      ieff = 0         
             
      do i=il1, il2
         if ((ivegcov(i) .eq.8 .or.  ivegcov(i) .eq.11) .and.
     &      (isoiltex(i).gt.0 .and. isoiltex(i).le.12) )then
          ieff = ieff + 1
          sflt(i)  = xsflt(ieff)  
          do n=1,nbin
            rsfrow(i,n) =  xrsfrow(ieff,n)  
          end do
        end if
      end do 

      return
      end subroutine sfflux

      subroutine dust_module(il1,il2,ilg,trsize,
     $                       soilw, vegfrac, surfwd, fland,
     &                       clayrow,sand2row,
     $                       roarow,z0,srel,
     $                       rsfrow,luc)

      implicit none

#     include "parameter.inc"
#     include "include/dust.cb"
#     include "include/message.cb"

      integer il1, il2, ilg
      integer i, luc
      real(kind=8)  fland(ilg,20), soilw(ilg)
      real(kind=8)  clayrow(ilg)
      real(kind=8)  sand2row(ilg,nats)
      real(kind=8)  roarow(ilg), surfwd(ilg)
      real(kind=8)  vegfrac(ilg),z0(ilg)

      real(kind=8)  ustar(ilg,nats),utheff(ilg,nsoil)
      real(kind=8)  uth
      real(kind=8)  srl(ilg), rc(ilg), alamda(ilg)
      real(kind=8)  srel(ilg,nsoil,nats)
      real(kind=8)  umin, rhop, hc(ilg)
      real(kind=8)  ustarns, utmin, vk
      real(kind=8)  xz, br,ym
      real(kind=8)  sigr,xc,z0s,x
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  rsfrow(ilg,nbin)
      real(kind=8)  wprim(ilg)

      real(kind=8)  tempd
      real(kind=8)  cly1,cly2
      real(kind=8)  arc1,arc2

      data umin/15./,rhop/2500.000/, vk/0.4/

      data xz/0.25/, br/202.0/, ym/0.16/, sigr/1.45/
      data z0s/3.e-3/, x/10./

      do i = il1,il2

       srl(i) = z0(i)*100.0
       rc(i) = 1.0
       
       if (jfs .eq. 0) then

c   *****************************************************************
c   * raupach et al. (1993)                                      ****
c   *****************************************************************               

       xc=1.0/(sigr*ym)     

       if (vegfrac(i) .lt. 1.0) then
         alamda(i)= xz*(log(1.0 - vegfrac(i)))*(-1.0)
         arc1= sigr*ym*alamda(i)
         arc2= br*ym*alamda(i)

       if (arc1 .le. 1.0 .and. arc2 .le. 1.0 ) then

         rc(i)= (sqrt(1.0 - arc1)*sqrt(1.0 + arc2))

       end if
       end if
     
       else if (jfs .eq. 1) then

       rc(i) = 1 -
     &       (dlog (0.5e-2/z0s)/
     &       (dlog( 0.35*(x/z0s)**0.8) )  )

      end if

c *************************************************************
c threshold velocity correction for soil humidity hc
c***************************************************************

       if (jsoilm .eq. 0) then

           if (soilw(i) .lt. 0.0) then
            write(aline, *) 'hc, rc =', soilw(i), ' less than zero'
            call say
            call fatal(__FILE__,__LINE__,'NEGATIVE SOILW')
          else if (soilw(i) .lt. 0.03 ) then
            hc(i) = exp(22.7*soilw(i))
          else if (soilw(i) .ge. 0.03 ) then
            hc(i) = exp(95.3 * soilw(i) - 2.029)
          else
             hc(i) = 1.0   
          end if

       else if (jsoilm .eq. 1) then   

          cly1 = clayrow(i)
          cly2 = cly1 * cly1


          wprim(i) = 0.0014* cly2 + 0.17 * cly1

c          print*,'-----------',wprim(i),clayrow(i)
c          if (isoiltex(i).eq.1) then
c            wprim(i) = 1.
c          end if

c          if ((isoiltex(i).eq.2).or.(isoiltex(i).eq.3)) then
c            wprim(i) = 1.84
c          end if

c          if ((isoiltex(i).eq.4).or.(isoiltex(i).eq.5)
c     &                          .or.(isoiltex(i).eq.6)) then
c            wprim(i) = 2.5
c          end if

c          if (isoiltex(i).eq.7) then
c            wprim(i) = 3.5
c          end if

c          if ((isoiltex(i).eq.8).or.(isoiltex(i).eq.9)) then
c            wprim(i) = 6.8
c          end if

c          if (isoiltex(i).eq.10) then
c            wprim(i) = 10.
c          end if

c          if (isoiltex(i).eq.11) then
c            wprim(i) = 10.5
c          end if

c          if (isoiltex(i).eq.12) then
c            wprim(i) = 11.5
c          end if

          tempd    = dmax1(0.00001d0, soilw(i)*100.0 - wprim(i))

          if (soilw(i) .lt. wprim(i)) then
            hc(i) = sqrt(1.0+ 1.21* tempd**0.68)
          else
            hc(i) = 1.0
          end if

c          print*,'hc= ', i,hc(i),hc1(i)
       end if

c   *****************************************************************
c   * total correction factor for both hc and rc                 **** 
c   *****************************************************************

         rc(i) = rc(i)/ hc(i)
    
c *******************************************************************
c *     computation of the wind friction velocity               *****
c *     accounting for the increase of the roughness length     *****
c *     due to the saltation layer (gillette etal. jgr 103,     *****
c *     no. d6, p6203-6209, 1998                                *****
c *******************************************************************

      ustarns=(vk*100.*surfwd(i))/(log(1000./srl(i)))
c     ustarns = ustarnd(i)*100 ! cm.s-1
      utmin=(umin/(100.*vk*rc(i)))*log(1000./srl(i))

c *******************************************************************
c *     vk=0.4: karman constant, and 1000 cm =10 m:             *****
c *     the height of wind defined level. umin: 21 cm/s:        *****
c *    the minimal threshold wind friction velocity.            *****
c *******************************************************************

          if (surfwd(i) .ge. utmin) then
            ustar(i,:)=ustarns+0.3*(surfwd(i)-utmin)
     &                *(surfwd(i)-utmin)
          else
            ustar(i,:)=ustarns
          end if

          end do   ! end i loop

      call uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop,dp)

      call emission(ilg, il1, il2, luc, rhop, nsoil,nbin,
     %              nats, fland, uth, roarow,dp, rc,utheff,
     %              ustar, srel,rsfrow,trsize,
     %              sand2row, vegfrac)

      return
      end subroutine dust_module

      subroutine uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,
     %                    rhop,dp)

      integer il1, il2, ilg
      integer i, j
      integer nsoil, ust
      real(kind=8)  roarow(ilg), dp(nsoil)
      real(kind=8)  utheff(ilg,nsoil)
      real(kind=8)  rhop,  ustart0, ustart01

      do i= 1, nsoil
      do j= il1,il2
        if (ust .eq. 0) then
          utheff(j,i)= ustart0(rhop,dp(i),roarow(j))
        end if
        if (ust .eq. 1) then
          utheff(j,i)=ustart01(rhop,dp(i),roarow(j))
        end if
      end do
      end do

      return
      end subroutine uthefft

      real(kind=8) function ustart01(rhop,dum,rhair)

c   *****************************************************************
c   * calculate of ustar01(d) using iversen and white (1982)     ****
c   * for smoth surface:                                         ****
c   * coded by :                                                 ****
c   * ashraf s. zakey, 2003                                      **** 
c   * dum    : particle diameter [um]                            ****
c   * ustar0 : threshold frication velocity [m/s]                ****
c   *****************************************************************

      real(kind=8) dm, dum
      real(kind=8) a1,a2,a3,c1,c2,c3,c4,c5
      parameter(a1=0.2, a2=0.129, a3=0.12)
      parameter(c1=0.006, c2=1.928, c3=0.0858)
      parameter(c4=-0.0617,c5=2.5)
      real(kind=8)  y1,y2,y3,rhop,rhair,g
      parameter(y1=1331.647, y2=1.561228, y3=0.38194)
      parameter(g=9.806160)
      real(kind=8)  rep, term1,term2
      real(kind=8)  term, cvmgt

      dm = dum  !* 1.0e-4      ! cm
      rep = y1 * (dm**y2) + y3
      term1 = sqrt(1.0 + (c1/(rhop*g*0.1*(dm**c5))))
      term2 = sqrt(rhop*g*100.0*dm/rhair)        
      term= term1*term2
      ustart01= cvmgt(a2*term*(1.0-c3*exp(c4*(rep-10.0))),
     %                a2*term/sqrt(c2 * (rep**0.092)-1.0),
     %                rep .gt. 10.0)

      return
      end function ustart01


      real(kind=8) function ustart0(rhop,dum,rhoa)

c   *****************************************************************
c   *                                                            **** 
c   * modified by a.s.zakey, nov.2003                            ****
c   * y. shao, 13 june 2000                                      ****
c   * calculate ustar0(d) using shao and lu (2000) for uncovered ****
c   * dry surface                                                ****
c   * dum:    particle diameter                   [um]           ****
c   * ustar0: threshold friction velocity       [cm/s]           **** 
c   *****************************************************************
      real(kind=8)  dm, dum,gamma 
      parameter(gamma=3.0e-4)      ! a constant
      real(kind=8) pi,g,f,rhoa,rhop,sigma
      parameter(pi = 3.1415927)
      parameter(g  = 9.81)     ! gravity parameter    [m s^-2]
      parameter(f  = 0.0123)

      sigma = rhop/rhoa
      dm = dum*1.0e-2    !* 1.0e-6 
      ustart0 = f* (sigma * g * dm + gamma/ (rhoa*dm)) 
      ustart0 = sqrt(ustart0 )
      ustart0 = ustart0 * 100.0  
      return 
      end function ustart0

      real(kind=8) function cvmgt(val1, val2, cond)
c
      logical cond
      real(kind=8) val1, val2
c
      if (cond) then
        cvmgt=val1
      else
        cvmgt=val2
      end if
c
      return
      end function cvmgt

c ======================================================================
      subroutine emission(ilg, il1, il2, luc, rhop, nsoil,nbin,
     %                    nats, fland, uth, roarow,dp, rc,utheff,
     %                    ustar, srel,rsfrow,trsize,
     %                    sand2row, vegfrac) 

      implicit none

      integer isize,luc
      integer il1, il2, ilg
      integer nsoil, nats, nbin
      integer i, j, k
      integer n
      parameter(isize=12)
      real(kind=8)  srel(ilg,nsoil,nats)
      real(kind=8)  fsoil(ilg),fsoil1(ilg),fsoil2(ilg),fsoil3(ilg)
      real(kind=8)  ustar(ilg,nats),sand2row(ilg,nats)
      real(kind=8)  fland(ilg,luc),utheff(ilg,nsoil)
      real(kind=8)  roarow(ilg),rc(ilg)
      real(kind=8)  frac1(isize),frac2(isize),frac3(isize)
      real(kind=8)  aerosize (2,isize)
      real(kind=8)  dp(nsoil)
      real(kind=8)  rsfrowsub(ilg, isize)
      real(kind=8)  rsfrow(ilg, nbin)
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  const, beta
      real(kind=8)  e1, e2, e3
      real(kind=8)  fdp1,fdp2,dec,ec
      real(kind=8)  uth,p1,p2,p3
      real(kind=8)  d1, d2, d3
      real(kind=8)  pi, g, rwi
      real(kind=8)  alogdi
      real(kind=8)  amean1, amean2, amean3
      real(kind=8)  rhop,totv1,totv2,totv3
      real(kind=8)  sigma1, sigma2, sigma3
      real(kind=8)  asigma1, asigma2, asigma3
      real(kind=8)  vegfrac(ilg)
      real(kind=8)  aeffect,beffect
      real(kind=8)  f 

      data  const/2.3/, beta/16300./
      data  pi/3.1415926535897/, g/9.806160/
c alfaro 's values 
      data  e1/3.61/, e2/3.52/, e3/3.46/
      data  d1/1.5/, d2/6.7/, d3/14.2/
      data sigma1/1.7/,sigma2/1.6/,sigma3/1.5/

c  emissions  bins (sub-bins)       
      data aerosize/ 1.0e-08, 2.0e-08,
     &         2.0e-08, 4.0e-08,
     &         4.0e-08, 8.0e-08,
     &         8.0e-08, 1.6e-07,
     &         1.6e-07, 3.2e-07,
     &         3.2e-07, 6.4e-07,
     &         6.4e-07, 1.28e-06,
     &         1.28e-06, 2.56e-06,
     &         2.56e-06,5.12e-06,
     &         5.12e-06,10.4e-06,
     &         10.24e-06,20.48e-06,
     &         20.48e-06,40.6e-06/ 

c      call putzero(fsoil1,ilg)
c      call putzero(fsoil2,ilg)
c      call putzero(fsoil3,ilg)

       fsoil (:) =0.
       fsoil1(:) =0.
       fsoil2(:) =0.
       fsoil3(:) =0.



c15    format(3(f12.4,1x))
      do i=1,nats
      do j=1, nsoil
        do k=il1,il2

          if (rc(k) .gt. 0.0 .and. ustar(k,i) .ne. 0.) then
            uth=utheff(k,j)/(rc(k)*ustar(k,i))
                  
            if (uth .le. 1.0) then

              fdp1=ustar(k,i)**3*(1.0-uth*uth)
              fdp2=(1.0+uth)*const*(1.e-5)*roarow(k)/g

              if (fdp2 .le. 0.0) fdp2=0.

              f = 0.0d0
              aeffect = (1-f) * (1-vegfrac(k))
              beffect = 0.01 * fland(k,i) * sand2row(k,i)

              fsoil(k)=srel(k,j,i)*fdp1*fdp2*aeffect*beffect

c              fsoil(k)=srel(k,j,i)*fdp1*fdp2*(1-f)*(1-vegfrac(k))

!     size-distributed kinetic energy flux
              dec=fsoil(k)*beta
!     individual kinetic energy for an aggregate of size dp ( g cm2 s-2)
! cf alfaro (dp) is in cm              
c              ec=(pi/3.)*1.e-1*rhop*(dp(j)**3.0)*(ustar(k,i)**2.0)
              ec= (pi/12)*rhop*1e-3*(dp(j)**3.0)*(20*ustar(k,i))**2.0
                        
              if (ec .gt.e1) then 
                p1=(ec-e1)/(ec-e3)
                p2=(1-p1)*(ec-e2)/(ec-e3)
                p3=1-p1-p2
              else if ( ec .gt. e2 .and. ec .le. e1) then
                p1 =0.
                p2=(ec-e2)/(ec-e3)
                p3 = 1-p2
              else if (  ec .gt. e3 .and. ec .le. e2)then 
                p1=0.
                p2=0.
                p3=1.
              else if ( ec .le. e3) then 
                p1 =0.
                p2 =0.
                p3 =0.
              end if   

              fsoil1(k)=fsoil1(k)+1.e-2*p1*(dec/e1)*(pi/6.)*rhop*
     &                  ((d1*1.e-04)**3.)
              fsoil2(k)=fsoil2(k)+1.e-2*p2*(dec/e2)*(pi/6.)*rhop*
     &                  ((d2*1.e-04)**3.)
              fsoil3(k)=fsoil3(k)+1.e-2*p3*(dec/e3)*(pi/6.)*rhop*
     &                  ((d3*1.e-04)**3.)
            end if
          end if
        end do
      end do
      end do

      totv1=0.0
      totv2=0.0
      totv3=0.0

      do n=1,isize
        rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.e6
        alogdi = log10(rwi)
        amean1=log10(d1)
        amean2=log10(d2)
        amean3=log10(d3)

        asigma1 = log10(sigma1)
        asigma2 = log10(sigma2) 
        asigma3 = log10(sigma3)

        frac1(n) = exp (-(alogdi - amean1)**2./ (2*asigma1**2))
        frac2(n) = exp (-(alogdi - amean2)**2./ (2*asigma2**2))
        frac3(n) = exp (-(alogdi - amean3)**2./ (2*asigma3**2))

        totv1=totv1+frac1(n)
        totv2=totv2+frac2(n)
        totv3=totv3+frac3(n)
      end do
   
      do n=1,isize
        frac1(n)=frac1(n)/totv1
        frac2(n)=frac2(n)/totv2
        frac3(n)=frac3(n)/totv3
        if (frac1(n) .lt. 1.e-9)  frac1(n)=0.0
        if (frac2(n) .lt. 1.e-9)  frac2(n)=0.0
        if (frac3(n) .lt. 1.e-9)  frac3(n)=0.0
      end do

      do 842 n = 1,isize
      do 842 i=il1,il2

c discretisation of the modal emission in isize emission sub bin         
        rsfrowsub(i,n) = fsoil1(i) *frac1(n) + fsoil2(i)*frac2(n)
     &                 + fsoil3(i) *frac3(n)
        


c and in tranport bins (nbin)
        rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.e6
        do k = 1,nbin
          if ( rwi .ge.trsize(k,1) .and. rwi .lt. trsize(k,2)) then
            rsfrow(i,k) = rsfrow (i,k) +  rsfrowsub(i,n)
          end if 
        end do
842   continue

      return
      end subroutine emission
