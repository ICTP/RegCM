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

      SUBROUTINE CHDRYDEP(ILG,IL1,IL2,ILEV,LUC,NBIN, IVEGCOV,
     &         THROW, ROAROW, SHJ, PRESSG,TEMP2, SUTEMP,
     &         SRAD, RH10,WIND10,ZEFF,TRSIZE,
     &         PDEPV )   
C     ***************************************************************
C     *  DRY DEPOSTION SCHEME FOR DUST PARTICLES                *****
C     *  THIS SCHEME BASED ON :                                 *****
C     * - ZHANG ET AL,(2001) : A SIZE-SEGREGATED PARTICLE       *****
C     *   DRY DEPOSITION SCHEME FOR AN ATMOSPHERIC AEROSOL      *****
C     *   MODULE, ATMOS. ENV. 35, 549-560                       *****
C     *                                                         *****
C     * - GIORGI, F. (1986): A PARTICLE DRY DEPOSITION          *****
C     *      PARAMETERIZATION SCHEME FOR USE IN TRACER          *****
C     *      TRANSPORT MODELS. JGR,91, 9794-9806                *****
C     *                                                         *****
C     *  INPUT:                                                 *****
C     *  =====                                                  *****
C     * - THROW : Temperature in K                              *****
C     * - ROAROW : Air density                                  *****
C     * - SHJ    : Local mid-layer sigma value                  *****
C     * - PRESSG : Grid row of surface pressure [Pa]            *****
C     * - TEMP2  : TEMPERATURE AT 10M. (DEG K)                  *****
C     * - SUTEMP : SURFACE TEMPERATURE (DEG K)                  *****
C     * - SRAD   : SOLAR IRRADIANCE AT THE GROUND(W/M**2)       *****
C     * - RH10   : RELATIVE HUMIDITY OF AIR AT 10M.             *****
C     * -VEGCOVER: VEGETATION COVER                             *****
C     * - WIND10 : WIND at 10 M                                 *****
C     * - XROW   : Dust concentration      [kg/kg]              *****
C     * - USTAR  : U*                                           *****
C     * - ILEV   : Number of model level                        *****
C     *                                                         *****
C     * OUTPUT                                                  *****
C     * ======                                                  *****
C     * - RTDRY  : dry deposition tendency                      *****
C     *                                                         *****
C     ***************************************************************

      implicit none
      INTEGER N,L,I,K,J,JC
      INTEGER ISIZE,ILEV,IL1,IL2,ILG
      INTEGER LUC,LEV,nbin
      parameter(ISIZE = 12)    
c subbbin
      real(kind=8)  SHJ(ILEV),PRESSG(ILG),THROW(ILG,ILEV)
      real(kind=8)  ROAROW(ILG,ILEV)
      real(kind=8)  PDEPVSUB(ILG,ILEV,ISIZE)
      real(kind=8)  TAUREL(ILG,ILEV,ISIZE)
      real(kind=8)  AEROSIZE (2,ISIZE)
      real(kind=8)  USTAR(ILG,LUC)
      real(kind=8)  CFAC(ILG,ILEV,ISIZE)
      real(kind=8)  AMU(ILG,ILEV)
      real(kind=8)  SCHM(ILG)
      real(kind=8)  RS(ILG, LUC,ISIZE) 
      real(kind=8)  TEMP2(ILG), SUTEMP(ILG)
      real(kind=8)  WIND10(ILG), RH10(ILG)
      real(kind=8)  SRAD(ILG)
      real(kind=8)  RA(ILG,LUC),VEGCOVER(ILG,LUC)
      real(kind=8)  DELT, ADT2
!      real(kind=8)  arw
      real(kind=8)  AVESIZE(ISIZE)  
      INTEGER IVEGCOV(ILG) ,tot
      real(kind=8)  ZEFF(ILG) 
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  PDEPV(ILG,ILEV,nbin)

      data    AEROSIZE/ 1.0E-08,   2.0E-08,
     &                  2.0E-08,   4.0E-08,
     &                  4.0E-08,   8.0E-08,
     &                  8.0E-08,   1.6E-07,
     &                  1.6E-07,   3.2E-07,
     &                  3.2E-07,   6.4E-07,
     &                  6.4E-07,   1.28E-06,
     &                  1.28E-06,  2.56E-06,
     &                  2.56E-06,  5.12E-06,
     &                  5.12E-06, 10.4E-06,
     &                  10.24E-06,20.48E-06,
     &                  20.48E-06,40.6E-06/ 

      real(kind=8)  PRE,PRII,PRIIV,AMFP
      real(kind=8)  A1,A2,A3,C1,C2,C3,C4,AA1,AA2,AA3,G,BOLTZK,RHOP,PI
      DATA    A1/145.8/, A2/1.5/, A3/110.4/
      DATA    C1/6.54E-8/,C2/1.818E-5/,C3/1.013E5/,C4/293.15/
      DATA    AA1/1.257/, AA2/0.4/, AA3/1.1/
      DATA    G/9.806160/,BOLTZK/1.3806044E-23/
      DATA    RHOP/2650.0/,PI/3.14159265/

      real(kind=8)  ARYE(20),AEST(20)
      DATA    AEST/0.80,  0.80, 0.8, 0.8, 1.2,
     1             1.20,  2.00, 1.5, 1.5, 2.0,
     2            15.0 , 15.0 , 1.5, 1.5, 1.5,
     5            15.0 ,  1.20, 1.2, 1.2, 1.2 /
      DATA    ARYE/0.5   ,5.0   ,0.5 , 5.0 , 1.0,
     1             1.0   ,0.0001,5.0 ,10.0 ,10.0,
     2             0.0001,0.0001,0.56, 0.56, 0.56,
     3             0.56  ,0.56  ,0.56, 0.56, 0.56/
      INTEGER KCOV
      real(kind=8)  FRX1,AMOB,EB,ST, EIM, EIN, R1
      real(kind=8)  ANU(ILG),RHSIZE(ILG,ILEV,ISIZE)
      real(kind=8)  PDIFF(ILG,ILEV,ISIZE)

      real(kind=8)  WW, MOL,PTEMP2,TBAR,RATIOZ,LOGRATIO,ASQ,DTEMP,VP,ES
      real(kind=8)  WVPM,USTARSQ,UTSTAR,THSTAR,KUI,AA,CM,CH,FM,FH,DELZ
      real(kind=8)  RIB,VPTEMP,TSW,QS,TSV,Z0WATER,DTHV,CUN,ZDL
      real(kind=8)  PSIU,X,Z,ZL,Y,PSIT,ZZ0(ILG)

      real(kind=8)  KRM,ZR,Z10,PMB
      DATA    KRM/0.4/,ZR/1.0/,Z10/10.0/,PMB /1000.0/

      DO N=1,ISIZE
        AVESIZE(N) = (AEROSIZE(1,N)+AEROSIZE(2,N))/2.0
      end do

      DELT=1200.0
      ADT2=2.0*DELT

c======================================================================
C  ********************************************************
C  *   AEROSIZE - DRY RADIUS                           ****
C  *   RHOP  - DENSITY FOR EACH AEROSOL TYPE           ****
C  ********************************************************
      DO N=1,ISIZE
         DO L=1,ILEV
           DO I=IL1,IL2  

C  ********************************************************
C  *  AEROSOL GRAVITATIONAL SETTLING VELOCITY          ****
C  *  AND DIFFUSION COEFFICIENT                        ****
C  *                                                   ****
C  * AIR'S DYNAMIC VISCOSITY                           ****
C  ********************************************************

             AMU(I,L)=A1*1.E-8*THROW(I,L)**A2/(THROW(I,L)+A3)

C . . . . MID LAYER PRESSURE IN [PASCAL].
             PRE = PRESSG(I)*SHJ(L)
C  ********************************************************
C  * MEAN MOLECULAR FREE PATH.                         ****
C  *     K.V. BEARD [1976], J ATM. SCI., 33            ****
C  ********************************************************

             AMFP=C1*(AMU(I,L)/C2)*(C3/PRE)*(THROW(I,L)/C4)**(1./2.)
             PRII=2./9.*G/AMU(I,L)
             PRIIV=PRII*(RHOP-ROAROW(I,L))

C  ********************************************************
C  * CUNNINGHAM SLIP CORRECTION FACTOR AND             ****
C  * RELAXATION TIME = vg/Grav.                        ****
C  ********************************************************

             CFAC(I,L,N) = 1.+AMFP/AVESIZE(N)*(AA1+AA2*
     &            EXP(-AA3*AVESIZE(N)/AMFP))
             TAUREL(I,L,N)=DMAX1(PRIIV*AVESIZE(N)**2*CFAC(I,L,N)/G,0.D0)

C  ********************************************************
C  * STOKES FRICTION                                  *****  
c    PDEPVSUB(I,L,N) ' sellting dep. velocity = '
C  ********************************************************

             PDEPVSUB(I,L,N) =TAUREL(I,L,N)*G
           end do
         end do
      end do
c======================================================================

c======================================================================
C     ****************************************************
C     * RA : IS THE AERODYNAMIC RESISTANCE ABOVE THE  ****
C     *      CANOPY AND IT IS FUNCTION IN U* AND      ****
C     *      Z0: ROUGHNESS LENGTH AND THE STABILITY   ****
C     *      FUNCTION                                 ****
C     * MOL  - MONIN OBUKHOV LENGTH (M) - CALCULATED  ****
C     *           FOR EACH LAND USE CATEGORY          ****
C     * PTEMP2 -POTENTIAL TEMPERATURE AT Z2  (DEG. K) ****
C     * TEMP2 - TEMPERATURE AT 10M. (DEG K)           ****
C     * Z10   - 10 M.                                 ****
C     * SUTEMP -SURFACE TEMPERATURE (DEG K)           ****
C     * SRAD   -SOLAR IRRADIANCE AT THE GROUND(W/M**2)****
C     * RH10  - RELATIVE HUMIDITY OF AIR AT 10M.      ****
C     *           (0.0-1.0)                           ****
C     * PMB   - SEA LEVEL PRESSURE (MB)               ****
C     ****************************************************
      DELZ=Z10-ZR
      DO J=1,LUC
        DO I=IL1,IL2
          WW = DMAX1(WIND10(I),1.0D0)
          ZZ0(I) = ZEFF(I)              

C     ***************************************************************
C     *  POTENTIAL TEMPERATURE AT Z2  (DEG. K)                  ***** 
C     ***************************************************************

          PTEMP2 = TEMP2(I) + Z10 * 0.0098

C     ***************************************************************
C     * FOR CALCULATIONS OVER WATER COMPUTE VALUES OF CRITICAL  *****
C     * PROFILE VARIABLES: L AND USTAR                          *****
C     *           ******BEGIN FOR WATER***                      *****
C     ***************************************************************

          IF ( IVEGCOV(I) .EQ. 0) THEN 
   
C     **************************************************************
C     *  VP  - VAPOUR PRESSURE AT Z2                           *****
C     *  WVPM- WATER VAPOUR MIXING RATIO AT  Z2                *****
C     *  VPTEMP- VIRTUAL POTENTIAL TEMPERATURE AT Z2 (DEG. K)  *****
C     **************************************************************
            ES = 6.108*EXP(17.27*(TEMP2(I)-273.16)/(TEMP2(I)-35.86))
            VP = RH10(I) * ES
            WVPM = 0.622 * VP/(PMB-VP)           
            VPTEMP = PTEMP2 * (1.0 + 0.61 * WVPM)

C     **************************************************************
C     *  ASSUME RH10 AT WATER SURFACE IS 100%                  *****
C     *   VP = ES(TSW-273.15) !SAT. VAP PRESS AT SURFACE       *****
C     *   SATURATED VAPOUR PRESSURE AT SURFACE                 *****
C     *   SATURATED MIXING RATIO AT SURFACE                    *****
C     *   TSV - VIRTUAL POTENTIAL TEMPERATURE AT SURFACE       *****
C     *           (DEG. K)                                     *****
C     **************************************************************
      
            TSW= SUTEMP(I)
            VP = 6.108*EXP(17.27*(TSW - 273.16)/(TSW - 35.86))
            QS = 0.622 * VP / (PMB - VP)
            TSV = TSW * (1. + 0.61 * QS)
            Z0WATER = 1.0E-4
C     **************************************************************
C     * SCALET  :  NOT REQUIRED IF  Z2 = 10M                   *****
C     **************************************************************
            DTHV = (VPTEMP - TSV)
C     **************************************************************
C     * CALCULATE DRAG COEFFICIENT CUN WITH NEUTRAL CONDITION  *****
C     * ASSUMPTION  GARRATT (1977)                             ***** 
C     **************************************************************
                
            CUN=7.5E-4+6.7E-5*WW
            MOL=9999.0

      IF(ABS(DTHV).GT.1.0E-6)MOL=VPTEMP*CUN**1.5*WW**2/(5.096E-3*DTHV)
            IF(MOL.GT.0..AND.MOL.LT.5.0) MOL = 5.0
            IF(MOL.GT.-5.0.AND.MOL.LT.0) MOL =-5.0
            ZDL = Z10/MOL
C
            IF(ZDL.LT.0.0) THEN

C     **************************************************************
C     *                        WIND SPEED                      *****
C     ************************************************************** 

              X=(1.0 - 15.0*ZDL)**0.25
              PSIU=2.*DLOG(0.5*(1.0+X))+DLOG(0.5*(1.0+X*X))-
     &                     2.0*ATAN(X) + 0.5*PI

C     **************************************************************
C     *                       POT TEMP                         *****
C     **************************************************************
              Y=SQRT(1.-9.*ZDL)                   
              PSIT=2.*0.74*DLOG((1+Y)/2.0)
            ELSE
              PSIU = -4.7*ZDL
              PSIT = PSIU
            end if
            Z0WATER = 0.000002 *WW**2.5
C
            USTAR(I,J)= KRM*WW/(DLOG(Z10/Z0WATER) - PSIU)
            THSTAR = KRM*(PTEMP2-SUTEMP(I))/(0.74*
     &                   DLOG(Z10/Z0WATER)-PSIT)
            ZZ0(I)=Z0WATER
C
          ELSE

C     **************************************************************
C     * COMPUTE USTAR AND L FOR LAND USE CATEGORIES OTHER THAN  ****
C     * WATER USE LOUIS METHOD. !PKK 7/16/85, FIND BULK         **** 
C     * RICHARDSON NUMBER.                                      ****
C     **************************************************************

            RIB = G*Z10*(PTEMP2 - SUTEMP(I))/(SUTEMP(I)*WW**2)

C     ***************************************************************
C     * ENSURE THAT CONDITIONS OVER LAND ARE NEVER STABLE WHEN  ***** 
C     * THERE IS INCOMING SOLAR RADIATIOM                       *****
C     *************************************************************** 
            IF (SRAD(I).GT.0.0.AND.RIB.GT.0.0) RIB = 1.E-15
C
            DTEMP = PTEMP2 - SUTEMP(I)
            IF(DABS(DTEMP).LT.1.E-10) DTEMP=DSIGN(1.D-10,DTEMP)
            TBAR = 0.5*(PTEMP2 + SUTEMP(I))
C
            RATIOZ = Z10/ZZ0(I)
            LOGRATIO = DLOG(RATIOZ)
            ASQ = 0.16/(LOGRATIO**2)
C
            IF(RIB.LE.0.0)THEN
              AA = ASQ*9.4*SQRT(RATIOZ)
              CM = 7.4*AA
              CH = 5.3*AA
              FM = 1. - (9.4*RIB/(1. + CM*SQRT(ABS(RIB))))
              FH = 1. - (9.4*RIB/(1. + CH*SQRT(ABS(RIB))))
            ELSE
              FM = 1./((1. + 4.7*RIB)**2)
              FH = FM
            end if
C
            USTARSQ = ASQ*WW**2*FM
            UTSTAR = ASQ*WW*DTEMP*FH/0.74
            USTAR(I,J) = SQRT(USTARSQ)
            THSTAR = UTSTAR/USTAR(I,J)
C
            MOL = TBAR*USTARSQ/(0.4*9.8*THSTAR)
          end if

          KUI=1.0/(KRM*USTAR(I,J))

C     **************************************************************  
C     * COMPUTE THE VALUES OF  RA                            *******
C     **************************************************************
C         Z=0.36788*(DELZ)+ZR
 
          Z=Z10
          ZL=Z/MOL

          IF(ZL.GE.0.) THEN
            RA(I,J)=KUI*(0.74*DLOG(Z/ZZ0(I))+4.7*ZL)
          ELSE
            RA(I,J)=KUI*0.74*(DLOG(Z/ZZ0(I))-
     &                         2.0*DLOG((1+SQRT(1-9.*ZL))*0.5))
          end if
          RA(I,J)=DMAX1(RA(I,J),0.99D0)
          RA(I,J)=DMIN1(RA(I,J),999.9D0)
        end do
      end do 
c======================================================================

c======================================================================
C find the right table index for the cell cover ( ocean and lake are 0
C      in the ivegcov and 14-15 in the table )

      if (ivegcov(I) .eq. 0) then 
        KCOV = 14
      else 
        KCOV = IVEGCOV(I)
      end if
 
C     *****************************************************
C     * THE SCHMIDT NUMBER IS THE RATIO OF THE         ****
C     * KINEMATIC VISCOSITY OF AIR TO THE PARTICLE     ****
C     * BROWNIAN DIFFUSIVITY ===> SC=V/D               ****
C     *****************************************************
         
      DO 555 N=1,ISIZE
        DO 666 L=1,ILEV
          DO 777 I=IL1,IL2

C     *****************************************************
C     * FOR NOW WE WILL NOT CONSIDER THE HUMIDITY      ****
C     * IMPACT SO WE WILL SET THE VARIABLE FRX1=1.0    ****
C     * I.E. ONLY DRY PARTICLES                        ****
C     *****************************************************

            FRX1 = 1.0
            RHSIZE(I,L,N)=AVESIZE(N)*FRX1
            ANU(I)=AMU(I,L)/ROAROW(I,L)
            AMOB=6.*PI*AMU(I,L)*RHSIZE(I,L,N)/CFAC(I,L,N)
            PDIFF(I,L,N)=BOLTZK*THROW(I,L)/AMOB
            SCHM(I)=ANU(I)/PDIFF(I,L,N)

C    ******************************************************
C    * FOR BROWNIAN DIFFUSION, THERE IS EVIDENCE THAT  ****
C    * ITS FROMULA DEPEND ON SCHMIDT NUMBER AS :       ****
c    * EB= SCHM x C^GAMA                               ****
C    * WHERE GAMA IS EFFICIENCY FACTOR AND ITS VALUE   ****
C    * BETWEEN 1/2 AND 2/3 WITH LARGER VALUES FOR      ****
C    * ROUGHER SURFACES                                ****
C    * ****************************************************
              
 777      CONTINUE
          IF (L .EQ. ILEV) THEN
            DO 400 K=1,LUC ! LUC  = 1 for the moment
              DO 200 I = IL1, IL2
         
C    ******************************************************
C    * The parameter governing impaction processes is *****
C    * the Stokes number,st, which has the form of    *****
C    * 1) ST = VG x U* /G A for vegetated surefaces   *****
C    *    (Slinn, 1982)                               *****
C    * 2) ST = VG x U*2/ANU for smothed surfaces or   *****
C    *    surfaces with bluff roughness elements      *****
C    *    (Giorgi,1988)                               *****
C    ******************************************************

                ST = TAUREL(I,L,N)* USTAR(I,K)*USTAR(I,K)/ANU(I)

                EB = SCHM(I)**(-0.666667)
c               EIM=(ST/(ST+AEST(K)))**2 
                EIM=(ST/(ST+AEST(KCOV)))**2 

              
                EIM=DMIN1(EIM,0.6D0)
                EIN = 0.0
c               IF (ARYE(K) .GT. 0.0001) THEN
c                  EIN = (1000.0*2.0*AVESIZE(N)/ARYE(K))**1.5
c               end if

                IF (ARYE(KCOV) .GT. 0.0001) THEN
                  EIN = (1000.0*2.0*AVESIZE(N)/ARYE(KCOV))**1.5
                end if

                EIN = DMIN1(EIN,0.5D0)

C     *****************************************************
C     * PARTICKES LARGER THAN 5 MICRO MAY REBOUNDED   *****
C     * AFTER HITTING A SURFACE,THIS PROCESS MAY BE   *****
C     * INCLUDED BY MODIFYING THE TOTAL COLLECTION    *****
C     * BY THE FACTOR OF R1, WHICH REPRESENTS THE     *****
C     * FRACTION OF PARTICLES STICKING TO THE SURFACE *****
C     * SLINN (1982) SUGGESTED THE FOLLOWING:         *****
C     * R = EXP (- ST^0.2)                            *****
C     ***************************************************** 
              
C               R1 = EXP (-ST**0.5)
                R1= DMAX1(0.5D0,  EXP(-ST**0.5))
c               IF(K .GE. 11 .AND. R1 .LT. 0.5 ) R1=0.5
                IF( KCOV .GE. 11 .AND. R1 .LT. 0.5 ) R1=0.5
                IF(R1 .LT. 0.4) R1=0.4
                 
C      ***************************************************
C      * CALCULATION OF RS: THE SURFACE RESISTANCE   *****
C      * WHICH DEPENDS ON THE COLLECTION EFFICIENCY  *****
C      * OF THE SURFACE AND IS DETERMINED BY THE     *****
C      * VARIOUS DEPOSITION PROCESSES                *****
C      ***************************************************

C               RS= 1.0/USTAR(I,K)/(EB+EIM+EIN)/R1
                RS(I,K,N) = 1.0/3.0/USTAR(I,K)/(EB+EIM+EIN)/R1
 200          CONTINUE                                    
 400        CONTINUE
          end if
 666    CONTINUE
 555  CONTINUE
c======================================================================

c NOTE for the moment VEGCOVER = 1 ( one type of cover IVEGcov PER GRIDE CELL) 

      DO N = 1,ISIZE 
        DO 12 JC = 1,LUC
          DO 13  I= IL1, IL2

c modify the surface layer 

            VEGCOVER(I,JC) =1.
            PDEPVSUB(I,ILEV,N) =  PDEPVSUB(I,ILEV,N) +
     &        VEGCOVER(I,JC) *  1.0/(RA(I,JC)+RS(I,JC,N) )
 
13        CONTINUE
12      CONTINUE      

      end do  ! end NSIZE loop

c average deposition velocities on bin 
      do k = 1,nbin
        tot = 0
        do lev=1,ilev
        do i=1,ilg
          PDEPV(i,lev,k) =0.
        enddo
        enddo
        do N=1,ISIZE                      
          if (  AVESIZE(N)*1.E6 .ge.trsize(k,1) .and. 
     &          AVESIZE(N)*1.E6 .lt.trsize(k,2)) then
            do lev=1,ilev
            do i=1,ilg
              PDEPV(i,lev,k)=PDEPV(i,lev,k)+PDEPVSUB(i,lev,N) 
            enddo
            enddo
            tot = tot + 1                 
          end if 
        end do
        if(tot.gt.0) then
          do lev=1,ilev
          do i=1,ilg
            PDEPV(i,lev,k) = PDEPV(i,lev,k) / tot
          enddo
          enddo
        end if
      end do

      return
      end


