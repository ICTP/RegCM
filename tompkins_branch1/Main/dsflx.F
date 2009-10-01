c  **********************************************************
c  *  Dust Emission Scheme                             ******
c  *                                                   ******
c  * This scheme based on Gong et al., (2003); Alfaro  ******
c  * et al.,(1997); Marticorena and Bergametti, 1995   ******
c  *                                                   ******
c  * The Modification coded by:                        ******
c  * Ashraf S. Zakey                                   ******
C  **********************************************************
      subroutine SFFLUX(ilg,il1,il2,jloop,ilev,luc,ivegcov,
     &                  vegfrac, isoiltex, sm,Z0,
     &                  soilw, surfwd, roarow, trsize,RSFROW)

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'Common2/dust.cb'
#else
      include 'Commons/dust.cb'
#endif

      INTEGER IL1,IL2,ILG,jloop
      INTEGER ILEV,LEV, N, NS
      INTEGER LUC
      INTEGER J,L,I,ieff, ieffmax
c  
      real*8  RSFROW(ILG,nbin)
      real*8  SOILW(ILG),sm(ilg)
      real*8  SURFWD(ILG),Z0(ILG) 
      integer ivegcov(ilg),isoiltex(ilg) 
      real*8  vegfrac(ilg), FLAND(ilg,20)
      real*8  ROAROW(ilg),SFLT(ilg)       
      real*8  trsize(nbin,2)

      INTEGER  IRB(ilg),ism(ilg)
c
c        effective work array  ( maybe think about allocatbale array)

      real*8  XRSFROW(ILG,nbin)
      real*8  XSOILW(ILG)
      real*8  XSURFWD(ILG) , XZ0(ILG)
      integer Xivegcov(ilg),Xisoiltex(ilg) 
      real*8  Xvegfrac(ilg), XFLAND(ilg,20)
      real*8  XROAROW(ilg),XSFLT(ilg)       
      real*8  XSREL2D(ilg,nsoil,nats)

      real*8 XSANDROW(ilg),  XCLAYROW(ilg)
      real*8 XSAND2ROW(ilg,nats), XCLAY2ROW(ilg,nats) 

      INTEGER  XIRB(ilg),Xism(ilg)      

      DO N = 1, 12
      DO I = IL1, IL2
         RSFROW(I,N) = 0.0
      END DO
      END DO

c  effective emitter cell ( depending on ivegcov) 
      do i=1,ilg
        Xivegcov(i) = 0.
        Xvegfrac(i) = 0. 
        Xisoiltex(i)= 0. 
        XSOILW(i)   = 0.
        XSURFWD(i)  = 0.
        XZ0(i)      = 0.
        XCLAYROW(i) = 0.
        XSANDROW(i) = 0.
        XIRB(i) = 0
        Xism(i) = 0  
        XSFLT(i) = 0.
        XROAROW(i) =0. 
      enddo

      do n=1,nats
      do ns=1,nsoil
      do i=1,ilg
        XSREL2D(i,ns,n)=0.
      enddo
      enddo
      enddo
      do n=1,nats
      do i=1,ilg
        XCLAY2ROW(i,n) =0.
        XSAND2ROW(i,n) =0.
      enddo
      enddo

      do n=1,20
      do i=1,ilg
        XFLAND(i,n)  =0. 
      enddo
      enddo
         
      do n=1,nbin
      do i=1,ilg
        XRSFROW(i,n) =0.
      enddo
      enddo

      ieff = 0  
      ieffmax=0
      do i=il1, il2
        if((ivegcov(i) .eq.8 .or.  ivegcov(i) .eq.11) .and. 
     &    (isoiltex(i).gt.0 .and. isoiltex(i).le.12)) then
          ieff = ieff + 1 
          Xivegcov(ieff) =  ivegcov(i)
          Xvegfrac(ieff) =  vegfrac(i) 
          Xisoiltex(ieff) = isoiltex(i) 
          XSOILW(ieff)  =  SOILW(i)
          XSURFWD(ieff) = SURFWD(i)
          XZ0(ieff)    = Z0(i)
          XROAROW(ieff) = ROAROW(i)

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
          Xism(ieff) = 1  
          do n=1,20
            XFLAND(ieff,n)  = 1
          enddo
        end if
      end do 

      ieffmax = ieff 

      if(ieffmax .gt. 0 ) then 
        call dust_module(1,ieffmax,ilg, ILEV,trsize,
     $                   XSOILW, XIVEGCOV, XVEGFRAC,XSURFWD, XFLAND,
     $                   XSANDROW,XCLAYROW,XSAND2ROW,XCLAY2ROW,
     $                   XROAROW, XZ0, XISOILTEX,XSREL2D,
     $                   XSFLT,XRSFROW,XiRB,Xism,LUC)
      end if


c  put back the dust flux on the right grid

      ieff = 0         
             
      do i=il1, il2
        if(ivegcov(i) .eq. 8 .or. ivegcov(i) .eq. 11) then
          ieff = ieff + 1
          SFLT(i)  = XSFLT(ieff)  
          do n=1,nbin
            RSFROW(i,n) =  XRSFROW(ieff,n)  
          enddo
        end if
      end do 

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine dust_module(IL1,IL2,ILG,ILEV,trsize,
     $                       SOILW, IVEGCOV, VEGFRAC, SURFWD, FLAND,
     &                       SANDROW,CLAYROW,SAND2ROW,CLAY2ROW,
     $                       ROAROW,Z0,ISOILTEX,SREL,SFLT, 
     $                       RSFROW,iRB,ism,LUC)

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'Common2/dust.cb'
#else
      include 'Commons/dust.cb'
#endif

      INTEGER IL1, IL2, ILG
      INTEGER I,L, LUC, ilev
      integer j,k 
      real*8  fland(ilg,20), soilw(ilg),silt2row(Ilg,nats)
      real*8  sandrow(ilg), clayrow(ilg),clay2row(ilg,nats)
      real*8  SAND2ROW(ilg,nats)
      real*8  roarow(ilg), surfwd(ilg)
      real*8  vegfrac(ilg),Z0(ilg)
      integer ivegcov(ilg),isoiltex(ilg)
      real*8  ustar(ilg,nats),utheff(ilg,nsoil)
      real*8  HT(ilg),uth
      real*8  srl(ilg), rc(ilg), alamda(ilg)
      real*8  srel(ilg,nsoil,nats),sflt(ilg)
      real*8  BCLY(12),BSND(12),BSLT(12),hh(20)
      real*8  BETA(ilg,nats,3),azz(Ilg),zz(ilg)
      REAL*8  UMIN, RHOP, HC
      REAL*8  USTARNS, UTMIN, vk
      real*8  xz, br,ym
      real*8  sigr,xc
      real*8  trsize(nbin,2)
      real*8  RSFROW(ILG,nbin)

      integer irb(ilg),ism(ilg)
  
      data UMIN/15./,rhop/2650.000/, vk/0.4/

      data xz/0.2/, br/10.0/, ym/0.16/, sigr/0.90/

      call smfv(il1,il2,ilg,soilw,clayrow,hc,jsoilm)

c **********************************************************
c * Raupach et al. (1993)                               ****
c **********************************************************
      xc=1.0/(sigr*ym)

      do i = il1,il2
        if(VEGFRAC(i).lt.1.0) then
          alamda(i)= xz*log(1.0 - VEGFRAC(i))*(-1.0)
             
          if(sigr*ym*alamda(i).LE.1..and. ym*br*alamda(i).le.1.) then
c
            rc(i)= sqrt(1.0-sigr*ym*alamda(i))*sqrt(1.0-ym*br*alamda(i))
          end if
        end if
      end do

CCCCCCC  Please if the value of roughness length calculated inside the model

      do i= il1, il2
        SRL(I) = Z0(I)*100.
      end do

c update the impact of soil moisture
             
      do i=1, nats
      do j=il1, il2
        rc(j) = rc(j)/hc
        if( rc(j) .le. 0.) then
          ustar(j,i) = 0
        else

c *******************************************************************
c *     COMPUTATION OF THE WIND FRICTION VELOCITY               ***** 
c *     ACCOUNTING FOR THE INCREASE OF THE ROUGHNESS LENGTH     *****
c *     DUE TO THE SALTATION LAYER (GILLETTE ETAL. JGR 103,     *****
c *     NO. D6, P6203-6209, 1998                                *****
c *******************************************************************

          ustarns=(vk*100.*surfwd(j))/(log(1000./srl(j)))
          utmin=(umin/(100.*vk*rc(j)))*log(1000./srl(j))

c *******************************************************************
c *     VK=0.4: KARMAN CONSTANT, AND 1000 CM =10 M:             *****
c *     THE HEIGHT OF WIND DEFINED LEVEL. UMIN: 21 CM/S:        *****
c *    THE MINIMAL THRESHOLD WIND FRICTION VELOCITY.            *****
c *******************************************************************

          if (surfwd(j) .ge. utmin) then
            ustar(j,i)=ustarns+0.3*(surfwd(j)-utmin)
     &                *(surfwd(j)-utmin)
          else
            ustar(j,i)=ustarns
          end if
        end if
      end do
      end do

      call uthefft(il1,il2,ilg,ilev,ust,nsoil,roarow,utheff,rhop,dp)

      call emission(ilg, il1, il2, ilev, luc, RHOP, nsoil,nbin,
     %              nats, fland, uth, roarow,dp, rc,utheff,
     %              ustar, srel,RSFROW,ivegcov,trsize,
     %              sand2row, irb,ism,sflt)

      return
      end

c=================================================================
c=================================================================

      subroutine smfv(il1,il2,ilg,soilw,clayrow,hc,jsoilm) 
      INTEGER IL1, IL2, ILG
      INTEGER K,JSOILM

      real*8  soilw(ilg), clayrow(ilg)
      real*8  wprim(ilg) 
      REAL*8  HC, TEMPD

      IF (jSOILM .EQ. 0) THEN
        DO K =IL1,IL2
          IF(SOILW(K) .lt. 0.0) THEN
            WRITE(6,*)'HC, RC =', SOILW(K), ' LESS THAN ZERO'
            STOP
          ELSEIF(SOILW(K) .lt. 0.03 ) THEN
            HC = EXP(22.7*SOILW(K))
          ELSEIF(SOILW(K) .ge. 0.03 ) THEN
            HC = EXP(95.3 * SOILW(K) - 2.029)
          END IF
          hc=hc*0.01
          if(hc .gt. 1.0) hc=1.0
        END DO
      end if

      IF(jSOILM .EQ. 1) THEN
        DO K =IL1,IL2
          WPRIM(K) = (0.0014*CLAYROW(K) + 0.17) * CLAYROW(K)
          TEMPD    = DMAX1(0.00001d0,SOILW(K)*100.0 -WPRIM(K))
          IF (SOILW(K) .gt.WPRIM(K)) THEN
            HC = SQRT(1.0+1.21*TEMPD**0.68)
          ELSE
            HC = 1.0
          END IF
        END DO
      END IF
c  
      return
      end 

c===========================================================

      subroutine uthefft(il1,il2,ilg,ilev,ust,nsoil,roarow,utheff,
     %                    rhop,dp)

      INTEGER IL1, IL2, ILG
      INTEGER I, J, K
      INTEGER ILEV, NSOIL, UST
      real*8  roarow(ilg), dp(nsoil)
      real*8  utheff(ilg,nsoil)
      REAL*8  RHOP,  ustart0, USTART01
      real*8  deldp

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

c=========================================================================

      REAL*8 FUNCTION USTART01(rhop,DUM,rhair)
c !=========================================================
c ! CALCULATE OF USTAR01(D) USING IVERSEN AND WHITE (1982)
c ! FOR SMOTH SURFACE:
c ! CODED BY :
c ! ASHRAF S. ZAKEY, 2003
c ! DUM    : PARTICLE DIAMETER [UM]
c ! USTAR0 : THRESHOLD FRICATION VELOCITY [M/S]
c !==========================================================

      REAL*8 DM, DUM
      REAL*8 a1,a2,a3,c1,c2,c3,c4,c5
      PARAMETER(A1=0.2, A2=0.129, A3=0.12)
      PARAMETER(C1=0.006, C2=1.928, C3=0.0858)
      PARAMETER(C4=-0.0617,C5=2.5)
      REAL*8  y1,y2,y3,rhop,rhair,g
      PARAMETER(Y1=1331.647, Y2=1.561228, Y3=0.38194)
      parameter(g=9.806160)
      REAL*8  REP, TERM1,TERM2
      REAL*8  TERM, cvmgt

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

c===============================================================

      REAL*8 FUNCTION USTART0(rhop,DUM,rhoa)
!==========================================================
! 
! MODIFIED BY A.S.ZAKEY, NOV.2003
! Y. Shao, 13 June 2000
!
! Calculate ustar0(d) using Shao and Lu (2000) for uncovered
! dry surface
!
! dum:    particle diameter			[um]
! ustar0: threshold friction velocity   	[cm/s]
!==========================================================
      REAL*8  DM, DUM,gamma 
      PARAMETER(GAMMA=3.0E-4)      ! a constant
      REAL*8 pi,g,f,rhoa,rhop,sigma
      PARAMETER(pi = 3.1415927)
      PARAMETER(g  = 9.81)     	! gravity parameter    [m s^-2]
      PARAMETER(F  = 0.0123)

      sigma = rhop/rhoa
      DM = DUM*1.0E-2    !* 1.0E-6 
      USTART0 = F* (SIGMA * G * DM + GAMMA/ (RHOA*DM)) 
      USTART0 = SQRT(USTART0 )
      USTART0 = USTART0 * 100.0  
      return 
      END 
c============================================================ 
      SUBROUTINE PUTZERO(F,IFINAL)

C     * ANDRE PLANTE MARCH 1993
C     * PUT ZERO IN NT VECTORS F
      INTEGER IFINAL ,I
      REAL*8 F(IFINAL)

      DO 100 I=1,IFINAL
        F(I)=0.
 100  CONTINUE

      RETURN
      END

      REAL*8 FUNCTION CVMGT(VAL1, VAL2, COND)
C
      LOGICAL COND
      REAL*8 VAL1, VAL2
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

c ======================================================================
      subroutine emission(ilg, il1, il2, ilev, luc, RHOP, nsoil,nbin,
     %                    nats, fland, uth, roarow,dp, rc,utheff,
     %                    ustar, srel,RSFROW,ivegcov,trsize,
     %                    sand2row, irb, ism,sflt) 

      implicit none

      integer isize,luc,np
      integer il1, il2, ilg
      integer nsoil, nats, nbin
      integer ilev, i, j, k
      integer n, irb(ilg),ism(ilg),ivegcov(ilg)
      parameter(ISIZE=12)
      real*8  srel(ilg,nsoil,nats)
      real*8  fsoil(ilg),fsoil1(ilg),fsoil2(ilg),fsoil3(ilg)
      real*8  ustar(ilg,nats),sand2row(ilg,nats)
      real*8  fland(ilg,luc),utheff(ilg,nsoil)
      real*8  roarow(ilg),rc(ilg),ASFLT(ilg)
      real*8  frac1(ISIZE),frac2(ISIZE),frac3(ISIZE)
      real*8  AEROSIZE (2,ISIZE)
      real*8  TRAC(ILG,ILEV,ISIZE), DSHJ(ILG,ILEV)
      real*8  PRESSG(ILG),dp(nsoil),sflt(ilg)
      real*8  RSFROWSUB(ilg, ISIZE)
      real*8  RSFROW(ilg, nbin)
      real*8  trsize(nbin,2)
      real*8  CONST, BETA
      real*8  e1, e2, e3
      real*8  FDP1,FDP2,DEC,EC
      real*8  uth,P1,P2,P3
      real*8  d1, d2, d3
      real*8  pi, g, RWI
      real*8  ALOGDI,DELT
      real*8  AMEAN1, AMEAN2, AMEAN3
      real*8  RHOP,TOTV1,TOTV2,TOTV3
      real*8  SUMF3,SUMF1,avsum
      real*8  SIGMA1, SIGMA2, SIGMA3
      real*8  ASIGMA1, ASIGMA2, ASIGMA3
      integer  IAE1, icont

      data  CONST/2.61/, BETA/16300./, IAE1/5/
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

      CALL PUTZERO(FSOIL1,ILG)
      CALL PUTZERO(FSOIL2,ILG)
      CALL PUTZERO(FSOIL3,ILG)

15    format(3(F12.4,1x))
      do i=1,nats
      do j=1, NSOIL
        do k=il1,il2

          if (rc(k) .gt. 0.0 .and. ustar(k,i) .ne. 0.) then
            uth=utheff(k,j)/(rc(k)*ustar(k,i))
                  
            if (uth .le. 1.0) then
              FDP1=USTAR(k,i)**3*(1.0+uth*uth)
              FDP2=(1.0-uth)*CONST*(1.E-5)*roarow(k)/g
              if(FDP2 .le. 0.0) FDP2=0.

              FSOIL(K)=srel(k,j,i)*FDP1*FDP2*fland(k,i)
     %                    *sand2row(k,i)*0.01 

!     SIZE-DISTRIBUTED KINETIC ENERGY FLUX
              DEC=FSOIL(k)*BETA
!     INDIVIDUAL KINETIC ENERGY FOR AN AGGREGATE OF SIZE DP ( g cm2 s-2)
! cf Alfaro (DP) is in cm              
              EC=(PI/3.)*1.E-1*RHOP*(DP(j)**3.0)*(USTAR(k,i)**2.0)
                        
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
