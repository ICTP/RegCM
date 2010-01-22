compiler on linux: pgf77 -byteswapio -o SIGMAtoP SIGMAtoP.f
compiler on sun  : f77 -o SIGMAtoP SIGMAtoP.f
C
C Note: this SIGMAtoP.f should be put under RegCM/???Run/output/
      implicit none
      integer iy,jx,kx,np
      parameter(iy=32,jx=49,kx=18)
      parameter(np=11)
      real    plev(np)
      data    plev/1000.,925.,850.,700.,500.,400.,300.,250.,200.
     &             ,150.,100./
      integer i,j,k,mrec,lrec
c
      integer nfile
      parameter (nfile=2)
      character*14  inout(nfile)
      data inout/'ATM.1994070100','ATM.1994070600'/
      integer   number(nfile)
      data      number/20,40/

      character*8  allout
      data allout /'PLEV_VAR'/

      integer   nf,mslice,idate

      real       fin(jx,iy)
      real      fout(jx,iy)
      real     xlong(jx,iy),difflon

      real         u(jx,iy,kx)
      real         v(jx,iy,kx)
      real         w(jx,iy,kx)
      real         t(jx,iy,kx)
      real         q(jx,iy,kx)
      real        qc(jx,iy,kx)
      real        ps(jx,iy)
      real      rain(jx,iy)
      real       tgb(jx,iy) ! temperature of lower soil layer
      real       swt(jx,iy) ! total soil water in mm H2O
      real       rno(jx,iy) ! accumulated infiltration
      real        ht(jx,iy),h(jx,iy,kx)
      real    satbrt(jx,iy)
      real      slp1(jx,iy),slp2(jx,iy)
      real        up(jx,iy,np)
      real        vp(jx,iy,np)
      real        wp(jx,iy,np)
      real        tp(jx,iy,np)
      real        hp(jx,iy,np)
      real        qp(jx,iy,np)
      real        cp(jx,iy,np)
      real      sigf(kx+1),sig(kx)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE

      integer mdate0,ibltyp,icup,ipptls,iboudy,il,jl,kl
      real    sp1d(kx+1),dxsp,ptsp,clat,clon,plat,plon
      character*6 proj
      real    dto,dtb,dtr,dtc
      integer iotyp

      open(49,file='OUT_HEAD',form='unformatted'
     &         ,access='direct',recl=iy*jx*4)
      read(49,rec=1) mdate0,ibltyp,icup,ipptls,iboudy
     &    , il,jl,kl,sp1d,dxsp,ptsp,clat,clon,plat,plon,proj
     &    , dto,dtb,dtr,dtc,iotyp
      if(il-2.ne.iy.or.jl-2.ne.jx.or.kl.ne.kx) then
         write(*,*) 'domain parameters are not consistent'
         write(*,*) 'il,jl,kl,iy,jx,kx',il,jl,kl,iy,jx,kx
         stop
      endif
      do k=1,kx+1
         sigf(k)=sp1d(kx+2-k)
      enddo
      ptop=ptsp*10.
      read(49,rec=2) ht
      read(49,rec=5) satbrt
      close(49)

      do k=1,kx
         sig(k) = 0.5*(sigf(k)+sigf(k+1))
      enddo

      open(50,file=allout,form='unformatted'
     &       ,access='direct',recl=iy*jx*4)

      mrec=0

      do nf=1,nfile
         if(iotyp.eq.1) then
            open(55,file=inout(nf),form='unformatted'
     &             ,recl=iy*jx*4  ,access='direct')
            lrec=0
         else if(iotyp.eq.2) then
            open(55,file=inout(nf),form='unformatted')
         endif
      do mslice=1,number(nf)
         write(*,*) ' mrec = ', mrec

         if(iotyp.eq.1) then
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((u(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((v(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((w(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((t(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((q(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               lrec=lrec+1
               read(55,rec=lrec) ((qc(i,j,k),i=1,jx),j=1,iy)
            enddo
            lrec=lrec+1
            read(55,rec=lrec) ((ps(i,j),i=1,jx),j=1,iy)
            lrec=lrec+1
            read(55,rec=lrec) ((rain(i,j),i=1,jx),j=1,iy)
            lrec=lrec+1
            read(55,rec=lrec) ((tgb(i,j),i=1,jx),j=1,iy)
            lrec=lrec+1
            read(55,rec=lrec) ((swt(i,j),i=1,jx),j=1,iy)
            lrec=lrec+1
            read(55,rec=lrec) ((rno(i,j),i=1,jx),j=1,iy)
         else if(iotyp.eq.2) then
            read(55) idate
            print*,' IDATE = ',idate
            do k=kx,1,-1
               read(55) ((u(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               read(55) ((v(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               read(55) ((w(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               read(55) ((t(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               read(55) ((q(i,j,k),i=1,jx),j=1,iy)
            enddo
            do k=kx,1,-1
               read(55) ((qc(i,j,k),i=1,jx),j=1,iy)
            enddo
            read(55) ((ps(i,j),i=1,jx),j=1,iy)
            read(55) ((rain(i,j),i=1,jx),j=1,iy)
            read(55) ((tgb(i,j),i=1,jx),j=1,iy)
            read(55) ((swt(i,j),i=1,jx),j=1,iy)
            read(55) ((rno(i,j),i=1,jx),j=1,iy)
         endif
C
C     to calculate Heights on sigma surfaces.
      CALL HTSIG(T,H,PS,HT,SIG,jx,iy,kx)
C
C     to calculate Sea-Level Pressure using
C         1. ERRICO's solution described in height
C         2. a simple formulae
C         3. MM5 method
      CALL SLPRES(H,T,PS,HT,TGB,SLP1,SLP2,SIG,jx,iy,kx)

C     to interpolate H,U,V,T,Q and QC
C        1. For Heights
      CALL HEIGHT(HP,H,T,PS,HT,SIG,jx,iy,kx,PLEV,NP)
C        2. For Zonal and Meridional Winds
      CALL INTLIN(UP,U,PS,SIG,jx,iy,kx,PLEV,NP)
      CALL INTLIN(VP,V,PS,SIG,jx,iy,kx,PLEV,NP)
      CALL INTLIN(WP,W,PS,SIG,jx,iy,kx,PLEV,NP)
C        3. For Temperatures
      CALL INTLOG(TP,T,PS,SIG,jx,iy,kx,PLEV,NP)
C        4. For Moisture qva & qca
      CALL HUMID1(T,Q,PS,SIG,jx,iy,kx)
      CALL INTLIN(QP,Q,PS,SIG,jx,iy,kx,PLEV,NP)
      CALL HUMID2(TP,QP,PLEV,jx,iy,NP)
      CALL INTLIN(CP,QC,PS,SIG,jx,iy,kx,PLEV,NP)

      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = hp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = tp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = up(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = vp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = wp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = qp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do k=1,np
        do j=1,iy
          do i=1,jx
            fout(i,j) = cp(i,j,k)
          enddo
        enddo
        mrec=mrec+1
        write(50,rec=mrec) fout
      enddo
      do j=1,iy
        do i=1,jx
          fout(i,j) = ps(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout
      do j=1,iy
        do i=1,jx
          fout(i,j) = slp1(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout
      do j=1,iy
        do i=1,jx
          fout(i,j) = rain(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout
      do j=1,iy
        do i=1,jx
          fout(i,j) = tgb(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout
      do j=1,iy
        do i=1,jx
          fout(i,j) = swt(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout
      do j=1,iy
        do i=1,jx
          fout(i,j) = rno(i,j)
        enddo
      enddo
      mrec=mrec+1
      write(50,rec=mrec) fout

      enddo
      close(55)
      enddo
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
