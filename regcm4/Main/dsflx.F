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

c  **********************************************************
c  *  Dust Emission Scheme                             ******
c  *                                                   ******
c  * This scheme based on Marticorena and Bergametti,  ******
c  * 1995; Gong et al.,(2003); Alfaro et al.,(1997)    ******
c  *                                                   ******
c  * The Modification coded by:                        ******
c  * Ashraf S. Zakey                                   ******
C  **********************************************************

      subroutine SFFLUX(ilg,il1,il2,jloop,luc,ivegcov,
     &                  vegfrac, isoiltex, ustarnd, Z0,
     &                  soilw, surfwd, roarow, trsize, RSFROW)

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/dust.cb"

      integer IL1,IL2,ILG,jloop
      integer N, NS
      integer LUC
      integer I,ieff, ieffmax
c  
      real(kind=8)  RSFROW(ILG,nbin)
      real(kind=8)  SOILW(ILG),ustarnd(ilg)
      real(kind=8)  SURFWD(ILG),Z0(ILG) 
      integer ivegcov(ilg),isoiltex(ilg) 
      real(kind=8)  vegfrac(ilg)
      real(kind=8)  ROAROW(ilg),SFLT(ilg)       
      real(kind=8)  trsize(nbin,2)

c
c        effective work array  ( maybe think about allocatbale array)

      real(kind=8)  XRSFROW(ILG,nbin)
      real(kind=8)  XSOILW(ILG)
      real(kind=8)  XSURFWD(ILG) , XZ0(ILG)
      real(kind=8)  Xvegfrac(ilg), XFLAND(ilg,20)
      real(kind=8)  XROAROW(ilg),XSFLT(ilg)       
      real(kind=8)  XSREL2D(ilg,nsoil,nats)
      real(kind=8)  Xustarnd(ilg) 
      real(kind=8) XSANDROW(ilg),  XCLAYROW(ilg)
      real(kind=8) XSAND2ROW(ilg,nats), XCLAY2ROW(ilg,nats)

      integer Xivegcov(ilg), Xisoiltex(ilg), XIRB(ilg)

      RSFROW = 0.0

c  effective emitter cell ( depending on ivegcov) 
      Xivegcov = 0.
      Xvegfrac = 0. 
      Xisoiltex= 0. 
      XSOILW   = 0.
      XSURFWD  = 0.
      XZ0      = 0.
      XCLAYROW = 0.
      XSANDROW = 0.
      XIRB = 0
      Xustarnd = 0. 
      XSFLT = 0.
      XROAROW =0.
      XSREL2D =0.
      XCLAY2ROW =0.
      XSAND2ROW =0.
      XFLAND =0.0
      XRSFROW =0.

      ieff = 0  
      ieffmax=0
      do i=il1, il2
        if((ivegcov(i) .eq.8 .or.  ivegcov(i) .eq.11) .and.
     &     (isoiltex(i).gt.0 .and. isoiltex(i).le.12) )then
          ieff = ieff + 1 
          Xivegcov(ieff) =  ivegcov(i)
          Xvegfrac(ieff) =  vegfrac(i) 
          Xisoiltex(ieff) = isoiltex(i) 
          XSOILW(ieff)  =  SOILW(i)
          XSURFWD(ieff) = SURFWD(i)
          XZ0(ieff)    = Z0(i)
          XROAROW(ieff) = ROAROW(i)
          Xustarnd(ieff)   =  ustarnd(i)

c soil parameters
c ok if one texture per grid cell
c (remember SREL2D,4 =1 because of memory problem for big domains)
          do ns=1,nsoil
            XSREL2D(ieff,ns,Xisoiltex(ieff) )=SREL2D(i,jloop,ns)
          enddo

          XSANDROW(ieff)=SANDROW2(i,jloop) 
          XCLAYROW(ieff)=CLAYROW2(i,jloop) 
          do n=1,nats
            XSAND2ROW(ieff,n)=SAND2ROW2(i,n,jloop)
            XCLAY2ROW(ieff,n)=CLAY2ROW2(i,n,jloop)
          enddo
c  trivial !
          XIRB(ieff) = 1
          do n=1,20
            XFLAND(ieff,n)  = 1
          enddo
        end if
      end do 

      ieffmax = ieff 

      if(ieffmax .gt. 0 ) then 
        call dust_module(1,ieffmax,ilg, trsize,
     $                   XSOILW, XVEGFRAC,XSURFWD, XFLAND,
     $                   XCLAYROW,XSAND2ROW,
     $                   XROAROW, XZ0, XSREL2D,
     $                   XRSFROW,LUC)
      end if

c  put back the dust flux on the right grid

      ieff = 0         
             
      do i=il1, il2
         if((ivegcov(i) .eq.8 .or.  ivegcov(i) .eq.11) .and.
     &      (isoiltex(i).gt.0 .and. isoiltex(i).le.12) )then
          ieff = ieff + 1
          SFLT(i)  = XSFLT(ieff)  
          do n=1,nbin
            RSFROW(i,n) =  XRSFROW(ieff,n)  
          enddo
        end if
      end do 

      return
      end

      subroutine dust_module(IL1,IL2,ILG,trsize,
     $                       SOILW, VEGFRAC, SURFWD, FLAND,
     &                       CLAYROW,SAND2ROW,
     $                       ROAROW,Z0,SREL,
     $                       RSFROW,LUC)

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/dust.cb"

      integer IL1, IL2, ILG
      integer I, LUC
      real(kind=8)  fland(ilg,20), soilw(ilg)
      real(kind=8)  clayrow(ilg)
      real(kind=8)  SAND2ROW(ilg,nats)
      real(kind=8)  roarow(ilg), surfwd(ilg)
      real(kind=8)  vegfrac(ilg),Z0(ilg)

      real(kind=8)  ustar(ilg,nats),utheff(ilg,nsoil)
      real(kind=8)  uth
      real(kind=8)  srl(ilg), rc(ilg), alamda(ilg)
      real(kind=8)  srel(ilg,nsoil,nats)
      real(kind=8)  UMIN, RHOP, HC(ilg)
      real(kind=8)  USTARNS, UTMIN, vk
      real(kind=8)  xz, br,ym
      real(kind=8)  sigr,xc,z0s,X
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  RSFROW(ILG,nbin)
      real(kind=8)  wprim(ilg)

      real(kind=8)  TEMPD
      real(kind=8)  cly1,cly2
      real(kind=8)  arc1,arc2

      data UMIN/15./,rhop/2500.000/, vk/0.4/

      data xz/0.25/, br/202.0/, ym/0.16/, sigr/1.45/
      data z0s/3.E-3/, X/10./

      do i = il1,il2

       srl(i) = z0(i)*100.0
       rc(i) = 1.0
       
       if (jfs .eq. 0) then

c   *****************************************************************
c   * Raupach et al. (1993)                                      ****
c   *****************************************************************               

       xc=1.0/(sigr*ym)     

       if(vegfrac(i) .lt. 1.0) then
         alamda(i)= xz*(log(1.0 - vegfrac(i)))*(-1.0)
         arc1= sigr*ym*alamda(i)
         arc2= br*ym*alamda(i)

       if(arc1 .le. 1.0 .and. arc2 .le. 1.0 ) then

         rc(i)= (sqrt(1.0 - arc1)*sqrt(1.0 + arc2))

       end if
       end if
     
       else if(jfs .eq. 1) then

       rc(i) = 1 -
     &       (dlog (0.5E-2/z0s)/
     &       (dlog( 0.35*(X/z0s)**0.8) )  )

      end if

c *************************************************************
c threshold velocity correction for soil humidity hc
c***************************************************************

       if (jSOILM .EQ. 0) THEN

           IF(SOILW(I) .lt. 0.0) THEN
            WRITE(6,*)'HC, RC =', SOILW(I), ' LESS THAN ZERO'
            STOP
          ELSEIF(SOILW(I) .lt. 0.03 ) THEN
            hc(i) = EXP(22.7*SOILW(I))
          ELSEIF(SOILW(I) .ge. 0.03 ) THEN
            hc(i) = EXP(95.3 * SOILW(I) - 2.029)
          ELSE
             hc(i) = 1.0   
          END IF

       elseif(jSOILM .EQ. 1) THEN   

          cly1 = clayrow(i)
          cly2 = cly1 * cly1


          wprim(i) = 0.0014* cly2 + 0.17 * cly1

c          print*,'-----------',WPRIM(I),CLAYROW(I)
c          if (isoiltex(i).eq.1) then
c            WPRIM(i) = 1.
c          endif

c          if ((isoiltex(i).eq.2).or.(isoiltex(i).eq.3)) then
c            WPRIM(i) = 1.84
c          endif

c          if ((isoiltex(i).eq.4).or.(isoiltex(i).eq.5)
c     &                          .or.(isoiltex(i).eq.6)) then
c            WPRIM(i) = 2.5
c          endif

c          if (isoiltex(i).eq.7) then
c            WPRIM(i) = 3.5
c          endif

c          if ((isoiltex(i).eq.8).or.(isoiltex(i).eq.9)) then
c            WPRIM(i) = 6.8
c          endif

c          if (isoiltex(i).eq.10) then
c            WPRIM(i) = 10.
c          endif

c          if (isoiltex(i).eq.11) then
c            WPRIM(i) = 10.5
c          endif

c          if (isoiltex(i).eq.12) then
c            WPRIM(i) = 11.5
c          endif

          tempd    = dmax1(0.00001d0, soilw(i)*100.0 - wprim(i))

          IF (SOILW(i) .lt. WPRIM(i)) THEN
            hc(i) = SQRT(1.0+ 1.21* tempd**0.68)
          ELSE
            hc(i) = 1.0
          END IF

c          print*,'HC= ', i,hc(i),hc1(i)
       END IF

c   *****************************************************************
c   * Total correction factor for both hc and rc                 **** 
c   *****************************************************************

         rc(i) = rc(i)/ hc(i)
    
c *******************************************************************
c *     COMPUTATION OF THE WIND FRICTION VELOCITY               *****
c *     ACCOUNTING FOR THE INCREASE OF THE ROUGHNESS LENGTH     *****
c *     DUE TO THE SALTATION LAYER (GILLETTE ETAL. JGR 103,     *****
c *     NO. D6, P6203-6209, 1998                                *****
c *******************************************************************

      ustarns=(vk*100.*surfwd(i))/(log(1000./srl(i)))
c     ustarns = ustarnd(i)*100 ! cm.s-1
      utmin=(umin/(100.*vk*rc(i)))*log(1000./srl(i))

c *******************************************************************
c *     VK=0.4: KARMAN CONSTANT, AND 1000 CM =10 M:             *****
c *     THE HEIGHT OF WIND DEFINED LEVEL. UMIN: 21 CM/S:        *****
c *    THE MINIMAL THRESHOLD WIND FRICTION VELOCITY.            *****
c *******************************************************************

          if (surfwd(i) .ge. utmin) then
            ustar(i,:)=ustarns+0.3*(surfwd(i)-utmin)
     &                *(surfwd(i)-utmin)
          else
            ustar(i,:)=ustarns
          end if

          end do   ! end i loop

      call uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop,dp)

      call emission(ilg, il1, il2, luc, RHOP, nsoil,nbin,
     %              nats, fland, uth, roarow,dp, rc,utheff,
     %              ustar, srel,RSFROW,trsize,
     %              sand2row, vegfrac)

      return
      end

      subroutine uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,
     %                    rhop,dp)

      integer IL1, IL2, ILG
      integer I, J
      integer NSOIL, UST
      real(kind=8)  roarow(ilg), dp(nsoil)
      real(kind=8)  utheff(ilg,nsoil)
      real(kind=8)  RHOP,  ustart0, USTART01

      do i= 1, nsoil
      do j= il1,il2
        if(ust .eq. 0) then
          utheff(j,i)= ustart0(rhop,dp(i),roarow(j))
        end if
        if(ust .eq. 1) then
          utheff(j,i)=USTART01(rhop,dp(i),ROAROW(j))
        end if
      end do
      end do

      return
      end     

      real(kind=8) FUNCTION USTART01(rhop,DUM,rhair)

c   *****************************************************************
c   * CALCULATE OF USTAR01(D) USING IVERSEN AND WHITE (1982)     ****
c   * FOR SMOTH SURFACE:                                         ****
c   * CODED BY :                                                 ****
c   * ASHRAF S. ZAKEY, 2003                                      **** 
c   * DUM    : PARTICLE DIAMETER [UM]                            ****
c   * USTAR0 : THRESHOLD FRICATION VELOCITY [M/S]                ****
c   *****************************************************************

      real(kind=8) DM, DUM
      real(kind=8) a1,a2,a3,c1,c2,c3,c4,c5
      PARAMETER(A1=0.2, A2=0.129, A3=0.12)
      PARAMETER(C1=0.006, C2=1.928, C3=0.0858)
      PARAMETER(C4=-0.0617,C5=2.5)
      real(kind=8)  y1,y2,y3,rhop,rhair,g
      PARAMETER(Y1=1331.647, Y2=1.561228, Y3=0.38194)
      parameter(g=9.806160)
      real(kind=8)  REP, TERM1,TERM2
      real(kind=8)  TERM, cvmgt

      DM = DUM  !* 1.0E-4      ! cm
      REP = Y1 * (DM**Y2) + Y3
      TERM1 = SQRT(1.0 + (C1/(RHOP*G*0.1*(DM**C5))))
      term2 = sqrt(rhop*g*100.0*DM/rhair)        
      term= term1*term2
      USTART01= cvmgt(a2*term*(1.0-c3*exp(c4*(rep-10.0))),
     %                a2*term/sqrt(C2 * (REP**0.092)-1.0),
     %                rep .gt. 10.0)

      return
      END     


      real(kind=8) FUNCTION USTART0(rhop,DUM,rhoa)

c   *****************************************************************
c   *                                                            **** 
c   * MODIFIED BY A.S.ZAKEY, NOV.2003                            ****
c   * Y. Shao, 13 June 2000                                      ****
c   * Calculate ustar0(d) using Shao and Lu (2000) for uncovered ****
c   * dry surface                                                ****
c   * dum:    particle diameter                   [um]           ****
c   * ustar0: threshold friction velocity       [cm/s]           **** 
c   *****************************************************************
      real(kind=8)  DM, DUM,gamma 
      PARAMETER(GAMMA=3.0E-4)      ! a constant
      real(kind=8) pi,g,f,rhoa,rhop,sigma
      PARAMETER(pi = 3.1415927)
      PARAMETER(g  = 9.81)     ! gravity parameter    [m s^-2]
      PARAMETER(F  = 0.0123)

      sigma = rhop/rhoa
      DM = DUM*1.0E-2    !* 1.0E-6 
      USTART0 = F* (SIGMA * G * DM + GAMMA/ (RHOA*DM)) 
      USTART0 = SQRT(USTART0 )
      USTART0 = USTART0 * 100.0  
      return 
      END 

      real(kind=8) FUNCTION CVMGT(VAL1, VAL2, COND)
C
      LOGICAL COND
      real(kind=8) VAL1, VAL2
C
      IF(COND) THEN
        CVMGT=VAL1
      ELSE
        CVMGT=VAL2
      ENDIF
C
      RETURN
      END

c ======================================================================
      subroutine emission(ilg, il1, il2, luc, RHOP, nsoil,nbin,
     %                    nats, fland, uth, roarow,dp, rc,utheff,
     %                    ustar, srel,RSFROW,trsize,
     %                    sand2row, vegfrac) 

      implicit none

      integer isize,luc
      integer il1, il2, ilg
      integer nsoil, nats, nbin
      integer i, j, k
      integer n
      parameter(ISIZE=12)
      real(kind=8)  srel(ilg,nsoil,nats)
      real(kind=8)  fsoil(ilg),fsoil1(ilg),fsoil2(ilg),fsoil3(ilg)
      real(kind=8)  ustar(ilg,nats),sand2row(ilg,nats)
      real(kind=8)  fland(ilg,luc),utheff(ilg,nsoil)
      real(kind=8)  roarow(ilg),rc(ilg)
      real(kind=8)  frac1(ISIZE),frac2(ISIZE),frac3(ISIZE)
      real(kind=8)  AEROSIZE (2,ISIZE)
      real(kind=8)  dp(nsoil)
      real(kind=8)  RSFROWSUB(ilg, ISIZE)
      real(kind=8)  RSFROW(ilg, nbin)
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  CONST, BETA
      real(kind=8)  e1, e2, e3
      real(kind=8)  FDP1,FDP2,DEC,EC
      real(kind=8)  uth,P1,P2,P3
      real(kind=8)  d1, d2, d3
      real(kind=8)  pi, g, RWI
      real(kind=8)  ALOGDI
      real(kind=8)  AMEAN1, AMEAN2, AMEAN3
      real(kind=8)  RHOP,TOTV1,TOTV2,TOTV3
      real(kind=8)  SIGMA1, SIGMA2, SIGMA3
      real(kind=8)  ASIGMA1, ASIGMA2, ASIGMA3
      real(kind=8)  vegfrac(ilg)
      real(kind=8)  Aeffect,Beffect
      real(kind=8)  f 

      data  CONST/2.3/, BETA/16300./
      data  pi/3.1415926535897/, g/9.806160/
c Alfaro 's values 
      data  e1/3.61/, e2/3.52/, e3/3.46/
      data  d1/1.5/, d2/6.7/, d3/14.2/
      data SIGMA1/1.7/,SIGMA2/1.6/,SIGMA3/1.5/

c  EMISSIONS  Bins (sub-bins)       
      data AEROSIZE/ 1.0E-08, 2.0E-08,
     &         2.0E-08, 4.0E-08,
     &         4.0E-08, 8.0E-08,
     &         8.0E-08, 1.6E-07,
     &         1.6E-07, 3.2E-07,
     &         3.2E-07, 6.4E-07,
     &         6.4E-07, 1.28E-06,
     &         1.28E-06, 2.56E-06,
     &         2.56E-06,5.12E-06,
     &         5.12E-06,10.4E-06,
     &         10.24E-06,20.48E-06,
     &         20.48E-06,40.6E-06/ 

c      CALL PUTZERO(FSOIL1,ILG)
c      CALL PUTZERO(FSOIL2,ILG)
c      CALL PUTZERO(FSOIL3,ILG)

       FSOIL (:) =0.
       FSOIL1(:) =0.
       FSOIL2(:) =0.
       FSOIL3(:) =0.



c15    format(3(F12.4,1x))
      do i=1,nats
      do j=1, NSOIL
        do k=il1,il2

          if (rc(k) .gt. 0.0 .and. ustar(k,i) .ne. 0.) then
            uth=utheff(k,j)/(rc(k)*ustar(k,i))
                  
            if (uth .le. 1.0) then

              FDP1=USTAR(k,i)**3*(1.0-uth*uth)
              FDP2=(1.0+uth)*CONST*(1.E-5)*roarow(k)/g

              if(FDP2 .le. 0.0) FDP2=0.

              f = 0.0d0
              Aeffect = (1-f) * (1-vegfrac(k))
              Beffect = 0.01 * fland(k,i) * sand2row(k,i)

              FSOIL(K)=srel(k,j,i)*FDP1*FDP2*aeffect*beffect

c              FSOIL(K)=srel(k,j,i)*FDP1*FDP2*(1-f)*(1-vegfrac(k))

!     SIZE-DISTRIBUTED KINETIC ENERGY FLUX
              DEC=FSOIL(k)*BETA
!     INDIVIDUAL KINETIC ENERGY FOR AN AGGREGATE OF SIZE DP ( g cm2 s-2)
! cf Alfaro (DP) is in cm              
c              EC=(PI/3.)*1.E-1*RHOP*(DP(j)**3.0)*(USTAR(k,i)**2.0)
              EC= (PI/12)*RHOP*1E-3*(DP(j)**3.0)*(20*USTAR(k,i))**2.0
                        
              if (EC .gt.E1) then 
                P1=(EC-E1)/(EC-E3)
                P2=(1-P1)*(EC-E2)/(EC-E3)
                P3=1-P1-P2
              else if ( EC .gt. E2 .and. EC .le. E1) then
                P1 =0.
                P2=(EC-E2)/(EC-E3)
                P3 = 1-P2
              else if (  EC .gt. E3 .and. EC .le. E2)then 
                P1=0.
                P2=0.
                P3=1.
              else if( EC .le. E3) then 
                P1 =0.
                P2 =0.
                P3 =0.
              end if   

              FSOIL1(k)=FSOIL1(k)+1.E-2*P1*(DEC/E1)*(PI/6.)*RHOP*
     &                  ((D1*1.E-04)**3.)
              FSOIL2(k)=FSOIL2(k)+1.E-2*P2*(DEC/E2)*(PI/6.)*RHOP*
     &                  ((D2*1.E-04)**3.)
              FSOIL3(k)=FSOIL3(k)+1.E-2*P3*(DEC/E3)*(PI/6.)*RHOP*
     &                  ((D3*1.E-04)**3.)
            END IF
          end if
        end do
      end do
      end do

      TOTV1=0.0
      TOTV2=0.0
      TOTV3=0.0

      DO N=1,ISIZE
        RWI = (AEROSIZE(1,N)+AEROSIZE(2,N))/2.0*1.E6
        ALOGDI = LOG10(RWI)
        AMEAN1=LOG10(D1)
        AMEAN2=LOG10(D2)
        AMEAN3=LOG10(D3)

        ASIGMA1 = LOG10(SIGMA1)
        ASIGMA2 = LOG10(SIGMA2) 
        ASIGMA3 = LOG10(SIGMA3)

        FRAC1(N) = EXP (-(ALOGDI - AMEAN1)**2./ (2*ASIGMA1**2))
        FRAC2(N) = EXP (-(ALOGDI - AMEAN2)**2./ (2*ASIGMA2**2))
        FRAC3(N) = EXP (-(ALOGDI - AMEAN3)**2./ (2*ASIGMA3**2))

        TOTV1=TOTV1+FRAC1(N)
        TOTV2=TOTV2+FRAC2(N)
        TOTV3=TOTV3+FRAC3(N)
      END DO
   
      DO N=1,ISIZE
        FRAC1(N)=FRAC1(N)/TOTV1
        FRAC2(N)=FRAC2(N)/TOTV2
        FRAC3(N)=FRAC3(N)/TOTV3
        IF (FRAC1(N) .LT. 1.E-9)  FRAC1(N)=0.0
        IF (FRAC2(N) .LT. 1.E-9)  FRAC2(N)=0.0
        IF (FRAC3(N) .LT. 1.E-9)  FRAC3(N)=0.0
      END DO

      DO 842 N = 1,ISIZE
      DO 842 I=IL1,IL2

c Discretisation of the modal emission in ISIZE emission sub bin         
        RSFROWSUB(I,N) = FSOIL1(I) *FRAC1(N) + FSOIL2(I)*FRAC2(N)
     &                 + FSOIL3(I) *FRAC3(N)
        


c and in tranport bins (nbin)
        RWI = (AEROSIZE(1,N)+AEROSIZE(2,N))/2.0*1.E6
        do k = 1,nbin
          if ( RWI .ge.trsize(k,1) .and. RWI .lt. trsize(k,2)) then
            RSFROW(I,k) = RSFROW (I,k) +  RSFROWSUB(I,N)
          end if 
        end do
842   CONTINUE

      RETURN
      END
