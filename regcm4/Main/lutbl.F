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

      SUBROUTINE LUTBL(PTOP)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
      real(kind=8)  PTOP
      INTEGER ITB,JTB
      real(kind=8)  A2,A3,A4
      PARAMETER(ITB=100,JTB=150,A2=17.2693882E0,A3=273.16E0,A4=35.86E0)
      real(kind=8)  PQ0,EPS
      PARAMETER(PQ0=379.90516E0,EPS=2.E-12)
      real(kind=8)  ELIWV,CP,R,ROVCP
      PARAMETER(ELIWV=2.72E6,CP=1004.E0,R=287.04,ROVCP=R/CP)
C
      real(kind=8)  PL,THL,RDQ,RDTH,RDP,RDTHE,QS0,SQS,THE0,STHE,PTBL
      real(kind=8)  TTBL
      COMMON /TBLS/ PL,THL,RDQ,RDTH,RDP,RDTHE,QS0(JTB),SQS(JTB),
     1             THE0(JTB),STHE(JTB),PTBL(ITB,JTB),TTBL(JTB,ITB)
      real(kind=8)  QSOLD(JTB),POLD(JTB),QSNEW(JTB),Y2P(JTB),PNEW(JTB)
      real(kind=8)  TOLD(JTB),THEOLD(JTB),THENEW(JTB),Y2T(JTB),TNEW(JTB)
      real(kind=8)  APP(JTB),AQP(JTB),APT(JTB),AQT(JTB)
      INTEGER KTHM,KPM,KTHM1,KPM1,KTH,KP
      real(kind=8)  THH,PH,DTH,DP,TH,P,APE,QS0K,SQSK,DQS,QS,THE0K,STHEK
      real(kind=8)  DTHE,PT
C
C--------------COARSE LOOK-UP TABLE FOR SATURATION POINT----------------
C
      PT=PTOP*1000.
! PTOP IN PASCAL
      KTHM=JTB
      KPM=ITB
      KTHM1=KTHM-1
      KPM1=KPM-1
C
      THL=210.
      THH=385.
      PL=PT
      PH=105000.
C
      DTH=(THH-THL)/dble(KTHM-1)
      DP=(PH-PL)/dble(KPM-1)
C
      RDTH=1./DTH
      RDP=1./DP
      RDQ=KPM-1
C
      TH=THL-DTH
C-----------------------------------------------------------------------
      DO 500 KTH=1,KTHM
        TH=TH+DTH
        P=PL-DP
        DO 510 KP=1,KPM
          P=P+DP
          APE=(100000./P)**(ROVCP)
          QSOLD(KP)=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
          POLD(KP)=P
  510   CONTINUE
C
        QS0K=QSOLD(1)
        SQSK=QSOLD(KPM)-QSOLD(1)
        QSOLD(1)=0.
        QSOLD(KPM)=1.
C
        DO 520 KP=2,KPM1
          QSOLD(KP)=(QSOLD(KP)-QS0K)/SQSK
CWWWWWWWWWWWWWW FIX DUE TO CYBER HALF PREC. LIMITATION WWWWWWWWWWWWWWWWW
          IF((QSOLD(KP)-QSOLD(KP-1)).LT.EPS)QSOLD(KP)=QSOLD(KP-1)+EPS
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  520   CONTINUE
C
        QS0(KTH)=QS0K
        SQS(KTH)=SQSK
C-----------------------------------------------------------------------
        QSNEW(1)=0.
        QSNEW(KPM)=1.
        DQS=1./dble(KPM-1)
C
        DO 530 KP=2,KPM1
  530   QSNEW(KP)=QSNEW(KP-1)+DQS
C
        Y2P(1)=0.
        Y2P(KPM)=0.
C
        CALL SPLINE(KPM,QSOLD,POLD,Y2P,KPM,QSNEW,PNEW,APP,AQP)
C
        DO 540 KP=1,KPM
  540   PTBL(KP,KTH)=PNEW(KP)
C-----------------------------------------------------------------------
  500 CONTINUE
C--------------COARSE LOOK-UP TABLE FOR T(P) FROM CONSTANT THE----------
      P=PL-DP
      DO 550 KP=1,KPM
        P=P+DP
        TH=THL-DTH
        DO 560 KTH=1,KTHM
          TH=TH+DTH
          APE=(100000./P)**(ROVCP)
          QS=PQ0/P*EXP(A2*(TH-A3*APE)/(TH-A4*APE))
          TOLD(KTH)=TH/APE
          THEOLD(KTH)=TH*EXP(ELIWV*QS/(CP*TOLD(KTH)))
  560   CONTINUE
C
        THE0K=THEOLD(1)
        STHEK=THEOLD(KTHM)-THEOLD(1)
        THEOLD(1)=0.
        THEOLD(KTHM)=1.
C
        DO 570 KTH=2,KTHM1
          THEOLD(KTH)=(THEOLD(KTH)-THE0K)/STHEK
CWWWWWWWWWWWWWW FIX DUE TO CYBER HALF PREC. LIMITATION WWWWWWWWWWWWWWWWW
          IF((THEOLD(KTH)-THEOLD(KTH-1)).LT.EPS)THEOLD(KTH)=
     +      THEOLD(KTH-1)+EPS
CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  570   CONTINUE
C
        THE0(KP)=THE0K
        STHE(KP)=STHEK
C-----------------------------------------------------------------------
        THENEW(1)=0.
        THENEW(KTHM)=1.
        DTHE=1./dble(KTHM-1)
        RDTHE=1./DTHE
C
        DO 580 KTH=2,KTHM1
  580   THENEW(KTH)=THENEW(KTH-1)+DTHE
C
        Y2T(1)=0.
        Y2T(KTHM)=0.
C
        CALL SPLINE(KTHM,THEOLD,TOLD,Y2T,KTHM,THENEW,TNEW,APT,AQT)
C
        DO 590 KTH=1,KTHM
          TTBL(KTH,KP)=TNEW(KTH)
  590   CONTINUE
C-----------------------------------------------------------------------
  550 CONTINUE
C
      RETURN
      END
