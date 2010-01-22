compiler on linux: pgf77 -byteswapio -o SIGMAtoP SIGMAtoP.f
compiler on sun  : f77 -o SIGMAtoP SIGMAtoP.f
C
C Note: this SIGMAtoP.f should be put under RegCM/Input/tools/
      implicit none
      include '../../PreProc/ICBC/icbc.param'
      integer npp
      parameter(npp=11)
      character*14 fn,fz
      data    fn/'ICBC1993010100'/
      data    fz/'IC_P1993010100'/
      integer i,j,k,nrec,mrec
      real       fin(jx,iy)
      real      fout(jx-2,iy-2)
      real      xlat(jx,iy),xlon(jx,iy)
      real         u(jx,iy,kz)
      real         v(jx,iy,kz)
      real         t(jx,iy,kz)
      real        qv(jx,iy,kz)
      real         w(jx,iy,kz+1),pp(jx,iy,kz)
      real        ps(jx,iy)
      real       tgb(jx,iy) ! temperature of lower soil layer
      real        ht(jx,iy),h(jx,iy,kz)
      real    satbrt(jx,iy)
      real      slp1(jx,iy),slp2(jx,iy)
      real        up(jx,iy,npp)
      real        vp(jx,iy,npp)
      real        tp(jx,iy,npp)
      real        hp(jx,iy,npp)
      real        qp(jx,iy,npp)
      real         a(jx,iy)
      real         b(jx,iy,86)
      real   sigf(kz+1),sig(kz),plev(npp)
      COMMON /ARRAY/ u,v,t,qv,w,pp,ps,tgb
     &             ,ht,h
     &             ,satbrt
     &             ,slp1,slp2
     &             ,up,vp,tp,hp,qp
     &             ,SIG,PLEV
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE

      data plev/1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100./
      INTEGER iyy,jxx,kzz,IDATE,NI,NJ,NK
      REAL    delx,clat,clon,plat,plon,grdfac
      CHARACTER*6 cgtype
      INTEGER igrads,ibigend
      INTEGER nday,nslice,numb
c
      open(61,file='../DOMAIN.INFO',form='unformatted'
     &         ,access='direct',recl=iy*jx*4)
      READ(61,rec=1) IYY,JXX,KZZ,DELX,CLAT,CLON,PLAT,PLON,GRDFAC
     &              ,CGTYPE,(SIGF(K),K=1,KZ+1),PTOP,igrads,ibigend
      IF(IYY.NE.IY.OR.JXX.NE.JX.OR.KZZ.NE.KZ) THEN
         WRITE(*,*)'There is inconsistence among parameters:'
         WRITE(*,*)'IY,JX,KZ,IYY,JXX,KZZ',IY,JX,KZ,IYY,JXX,KZZ
         STOP
      ENDIF
      do k=1,kz
         sig(k) = 0.5*(sigf(k)+sigf(k+1))
      enddo
      read(61,rec=2) ht
      read(61,rec=5) xlat
      read(61,rec=6) xlon
      close(61)
      open(10,file=fz,form='unformatted'
     &       ,recl=(jx-2)*(iy-2)*4,access='direct')
      nrec=0
      open(64,file='../'//fn,form='unformatted'
     &       ,access='direct',recl=iy*jx*4)
      MREC=0
      do nday=1,31
         NSLICE=4
         if(nday.eq.1) NSLICE=5
         do numb=1,NSLICE
            MREC=MREC+1
            READ(64,rec=MREC) IDATE,NI,NJ,NK
            if(NI.ne.JX.or.NJ.ne.IY.or.NK.ne.KZ) then
               write(*,*)'IDATE,NI,NJ,NK = ',IDATE,NI,NJ,NK
               STOP
            else
               write(*,*)'IDATE = ',IDATE
            endif
            DO K=NK,1,-1
               MREC=MREC+1
               READ(64,rec=MREC)((U(I,J,K),I=1,NI),J=1,NJ)
            ENDDO
            DO K=NK,1,-1
               MREC=MREC+1
               READ(64,rec=MREC)((V(I,J,K),I=1,NI),J=1,NJ)
            ENDDO
            DO K=NK,1,-1
               MREC=MREC+1
               READ(64,rec=MREC)((T(I,J,K),I=1,NI),J=1,NJ)
            ENDDO
            DO K=NK,1,-1
               MREC=MREC+1
               READ(64,rec=MREC)((QV(I,J,K),I=1,NI),J=1,NJ)
            ENDDO
            MREC=MREC+1
            READ(64,rec=MREC) PS
            MREC=MREC+1
            READ(64,rec=MREC) TGB

            do j=1,iy
            do i=1,jx
               ps(i,j)=ps(i,j)*10.
               tgb(i,j)=t(i,j,kz)
            enddo
            enddo
C
C           to calculate Heights on sigma surfaces.
            CALL HTSIG(T,H,PS,HT,SIG,jx,iy,kz)
C
C           to calculate Sea-Level Pressure using
C            1. ERRICO's solution described in height
C            2. a simple formulae
C            3. MM5 method
            CALL SLPRES(H,T,PS,HT,TGB,SLP1,SLP2,SIG,jx,iy,kz)

C           to interpolate H,U,V,T,Q and QC
C              1. For Heights
            CALL HEIGHT(HP,H,T,PS,HT,SIG,jx,iy,kz,PLEV,NPP)
C              2. For Zonal and Meridional Winds
            CALL INTLIN(UP,U,PS,SIG,jx,iy,kz,PLEV,NPP)
            CALL INTLIN(VP,V,PS,SIG,jx,iy,kz,PLEV,NPP)
C              3. For Temperatures
            CALL INTLOG(TP,T,PS,SIG,jx,iy,kz,PLEV,NPP)
C              4. For Moisture qva & qca
            CALL HUMID1(T,QV,PS,SIG,jx,iy,kz)
            CALL INTLIN(QP,QV,PS,SIG,jx,iy,kz,PLEV,NPP)
            CALL HUMID2(TP,QP,PLEV,jx,iy,NPP)
            do k=1,npp
               do j=1,iy-2
               do i=1,jx-2
                  fout(i,j) = hp(i+1,j+1,k)
               enddo
               enddo
               nrec=nrec+1
               write(10,rec=nrec) fout
            enddo
            do k=1,npp
               do j=1,iy-2
               do i=1,jx-2
                  fout(i,j) = tp(i+1,j+1,k)
               enddo
               enddo
               nrec=nrec+1
               write(10,rec=nrec) fout
            enddo
            do k=1,npp
               do j=1,iy-2
               do i=1,jx-2
                  fout(i,j) = up(i+1,j+1,k)
               enddo
               enddo
               nrec=nrec+1
               write(10,rec=nrec) fout
            enddo
            do k=1,npp
               do j=1,iy-2
               do i=1,jx-2
                  fout(i,j) = vp(i+1,j+1,k)
               enddo
               enddo
               nrec=nrec+1
               write(10,rec=nrec) fout
            enddo
            do k=1,npp
               do j=1,iy-2
               do i=1,jx-2
                  fout(i,j) = qp(i+1,j+1,k)
               enddo
               enddo
               nrec=nrec+1
               write(10,rec=nrec) fout
            enddo
            do j=1,iy-2
            do i=1,jx-2
               fout(i,j) = ps(i+1,j+1)
            enddo
            enddo
            nrec=nrec+1
            write(10,rec=nrec) fout
            do j=1,iy-2
            do i=1,jx-2
               fout(i,j) = slp1(i+1,j+1)
            enddo
            enddo
            nrec=nrec+1
            write(10,rec=nrec) fout
c           do j=1,iy-2
c           do i=1,jx-2
c              fout(i,j) = slp2(i+1,j+1)
c           enddo
c           enddo
c           nrec=nrec+1
c           write(10,rec=nrec) fout
         enddo
      enddo
c
      stop
      end
      SUBROUTINE HUMID2(T,Q,PRESLV,IM,JM,KP)
      implicit none
      INTEGER IM,JM,KP
      REAL    TR,QMIN
      PARAMETER (TR=1./273.16)
      PARAMETER (QMIN=0.0)   ! MINIMUM VALUE OF SPECIFIC HUMIDITY
      REAL    T(IM,JM,KP),Q(IM,JM,KP)
      REAL    PRESLV(KP)
      INTEGER I,J,K
      REAL    HL,SATVP,QS
C
C  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
C  DATA ON SIGMA LEVELS
C
      DO K=1,KP
        DO J=1,JM
          DO I=1,IM
            HL=597.3-.566*(T(I,J,K)-273.16)
            SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))
            QS=.622*SATVP/(PRESLV(K)-SATVP)           ! preslv (hPa)
            IF (Q(I,J,K).LT.QMIN) Q(I,J,K)=QMIN       ! SPECIFIED MINIMUM
            Q(I,J,K)=Q(I,J,K)*QS
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE HUMID1(T,Q,PS,SIGMA,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    TR,QMIN
      PARAMETER (TR=1./273.16)
      PARAMETER (QMIN=0.0)   ! MINIMUM VALUE OF SPECIFIC HUMIDITY
      REAL    T(IM,JM,KM),Q(IM,JM,KM)
      REAL    PS(IM,JM)
      REAL    SIGMA(KM)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K
      REAL    HL,SATVP,QS,P
C
C  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
C  DATA ON SIGMA LEVELS
C
      DO K=1,KM
        DO J=1,JM
          DO I=1,IM
            P=SIGMA(K)*(PS(I,J)-PTOP)+ptop
            HL=597.3-.566*(T(I,J,K)-273.16)           ! LATENT HEAT OF EVAP.
            SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K))) ! SATURATION VAP PRESS.
            QS=.622*SATVP/(P-SATVP)                   ! SAT. MIXING RATIO
            Q(I,J,K)=Q(I,J,K)/QS
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE INTLOG(FP,F,PSTAR,SIG,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P,KBC
      REAL    SIGP,WP,W1
C
C  INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
C        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
C        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
C        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
C        CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
C        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
C        TWO EXTREME TEMPERATURES IN THE LAYER.

C
C** FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      DO K=1,KM
        IF(SIG(K).LT.BLTOP) KBC = K
      ENDDO
      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = LOG(SIGP/SIG(K1)) / LOG(SIG(K1P)/SIG(K1))
              W1  = 1. - WP
              FP(I,J,N)= W1*F(I,J,K1) + WP*F(I,J,K1P)
            ELSE IF((SIGP.GE.SIG(KM)).AND.(SIGP.LE.1.))THEN
              FP(I,J,N)= F(I,J,KM)
            ELSE IF(SIGP.GT.1.) THEN
              FP(I,J,N) = F(I,J,KBC) 
     &        * EXP(-RGAS*TLAPSE*LOG(SIGP/SIG(KBC))/GRAV)
C                  ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
      RETURN
      END
C
      SUBROUTINE INTLIN(FP,F,PSTAR,SIG,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P
      REAL    SIGP,WP,W1
C
C  INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
C        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
C        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.

      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = (SIGP-SIG(K1))/(SIG(K1P)-SIG(K1))
              W1  = 1.-WP
              FP(I,J,N)  = W1*F(I,J,K1)+WP*F(I,J,K1P)
            ELSE IF(SIGP.GE.SIG(KM)) THEN
              FP(I,J,N)  = F(I,J,KM)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE HEIGHT(HP,H,T,PSTAR,HT,SIG,IM,JM,KM,P,KP)

C  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
C     ON INPUT:
C        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
C        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
C        SIG = SIGMA LEVELS.
C        P = PRESSURE LEVELS DESIRED.
C     ON OUTPUT:
C        ALL FIELDS EXCEPT H ARE UNCHANGED.
C        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
C
C  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
C  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
C     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
C     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
C     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
C     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
C      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
C
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    T(IM,JM,KM),H(IM,JM,KM),HP(IM,JM,KP)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      REAL    PSIG(100)
      INTEGER I,J,K,KBC,N,KT,KB
      REAL    PSFC,TEMP,WT,WB
C
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
      PRINT *,'FIRST SIGMA LEVEL ABOVE BNDY LAYER:', SIG(KBC)
C
      DO J=1,JM
        DO I=1,IM
          DO K=1,KM
            PSIG(K) = SIG(K) * (PSTAR(I,J)-PTOP) + PTOP
          ENDDO
          PSFC = PSTAR(I,J)
          DO N = 1,KP
            KT = 1
            DO K=1,KM
              IF (PSIG(K).LT.P(N)) KT=K
            ENDDO
            KB = KT + 1
            IF(P(N).LE.PSIG(1)) THEN
              TEMP = T(I,J,1)
              HP(I,J,N) =H(I,J,1)+RGAS*TEMP*LOG(PSIG(1)/P(N))/GRAV
            ELSE IF((P(N).GT.PSIG(1)).AND.(P(N).LT.PSIG(KM))) THEN
              WT = LOG(PSIG(KB)/P(N)) / LOG(PSIG(KB)/PSIG(KT))
              WB = LOG(P(N)/PSIG(KT)) / LOG(PSIG(KB)/PSIG(KT))
              TEMP = WT * T(I,J,KT) + WB * T(I,J,KB)
              TEMP = ( TEMP + T(I,J,KB) ) / 2.
              HP(I,J,N) =H(I,J,KB)+RGAS*TEMP*LOG(PSIG(KB)/P(N))/GRAV
            ELSE IF((P(N).GE.PSIG(KM)).AND.(P(N).LE.PSFC)) THEN
              TEMP = T(I,J,KM)
              HP(I,J,N) =HT(I,J)+RGAS*TEMP*LOG(PSFC/P(N))/GRAV
            ELSE IF(P(N).GT.PSFC) THEN
              TEMP = T(I,J,KBC) - TLAPSE * (H(I,J,KBC)-HT(I,J))
              HP(I,J,N) =HT(I,J)-(TEMP/TLAPSE)
     &              * ( 1.-EXP(-RGAS*TLAPSE*LOG(P(N)/PSFC)/GRAV))
C
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE SLPRES(H,T,PSTAR,HT,TG,SLP1,SLP2,SIG,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM)
      REAL    PSTAR(IM,JM),HT(IM,JM),TG(IM,JM)
      REAL    SLP1(IM,JM),SLP2(IM,JM)
      REAL    SIG(KM)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER KBC,I,J,K
      REAL    TSFC
C
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
      DO J=1,JM
         DO I=1,IM
            TSFC = T(I,J,KBC)-TLAPSE*(H(I,J,KBC)-HT(I,J))
            SLP1(I,J) = PSTAR(I,J)
     &      * EXP( -GRAV/(RGAS*TLAPSE)*LOG(1.-HT(I,J)*TLAPSE/TSFC))
         ENDDO
      ENDDO

      DO J=1,JM
         DO I=1,IM
            SLP2(I,J) = PSTAR(I,J)
     &      * EXP( GRAV*HT(I,J)/(RGAS*0.5*(TG(I,J)+288.15)))
         ENDDO
      ENDDO

      RETURN
      END
C
      SUBROUTINE HTSIG(T,H,PSTAR,HT,SIG,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K
      REAL    TBAR
C
      DO J=1,JM
      DO I=1,IM
         H(I,J,KM) = HT(I,J) + RGAS/GRAV*T(I,J,KM)
     &             * LOG(PSTAR(I,J)/((PSTAR(I,J)-PTOP)*SIG(KM)+PTOP))
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM
      DO I=1,IM
         TBAR = 0.5*( T(I,J,K)+T(I,J,K+1) )
         H(I,J,K) = H(I,J,K+1) +RGAS/GRAV*TBAR
     &            * LOG(((PSTAR(I,J)-PTOP)*SIG(K+1)+PTOP)
     &                 /((PSTAR(I,J)-PTOP)*SIG(K)+PTOP))
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
C
      BLOCK DATA
      implicit none
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      DATA    RGAS/287.04/
      DATA    GRAV/9.80616/
      DATA    BLTOP/.96/
      DATA    TLAPSE/-6.5E-3/
      END
