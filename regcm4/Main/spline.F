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

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SPLINE(NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
      IMPLICIT NONE
      INTEGER NOLD,NNEW
C     ******************************************************************
C     *                                                                *
C     *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
C     *  PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
C     *                                                                *
C     *  PROGRAMER\ Z. JANJIC, YUGOSLAV FED. HYDROMET. INST., BEOGRAD  *
C     *                                                                *
C     *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
C     *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
C     *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
C     *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
C     *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
C     *         SPECIFIED.                                             *
C     *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
C     *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
C     *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
C     *         AND LE XOLD(NOLD).                                     *
C     *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
C     *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
C     *                                                                *
C     ******************************************************************
      INTEGER JTB
      PARAMETER (JTB=150)
      real(kind=8)  XOLD(JTB),YOLD(JTB),Y2(JTB),P(JTB),Q(JTB),XNEW(JTB)
     &       ,YNEW(JTB)
      INTEGER NOLDM1,K,K1,K2,KOLD
      real(kind=8)  DXL,DXR,DYDXL,DYDXR,RTDXC
      real(kind=8)  DXC,DEN,XK,Y2K,Y2KP1,DX,RDX,AK,BK,CK,X,XSQ
C-----------------------------------------------------------------------
      NOLDM1=NOLD-1
C
      DXL=XOLD(2)-XOLD(1)
      DXR=XOLD(3)-XOLD(2)
      DYDXL=(YOLD(2)-YOLD(1))/DXL
      DYDXR=(YOLD(3)-YOLD(2))/DXR
      RTDXC=.5/(DXL+DXR)
C
      P(1)=RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
      Q(1)=-RTDXC*DXR
C
      IF(NOLD.EQ.3)GOTO 700
C-----------------------------------------------------------------------
      K=3
C
  100 DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
C
      P(K-1)=DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
C
      K=K+1
      IF(K.LT.NOLD)GOTO 100
C-----------------------------------------------------------------------
  700 K=NOLDM1
C
  200 Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
C
      K=K-1
      IF(K.GT.1)GOTO 200
C-----------------------------------------------------------------------
      K1=1
C
  300 XK=XNEW(K1)
C
      DO 400 K2=2,NOLD
        IF(XOLD(K2).LE.XK)GOTO 400
        KOLD=K2-1
        GOTO 450
  400 CONTINUE
      YNEW(K1)=YOLD(NOLD)
      GOTO 600
C
  450 IF(K1.EQ.1)GOTO 500
      IF(K.EQ.KOLD)GOTO 550
C
  500 K=KOLD
C
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
C
CVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C     WRITE(6,5000) K,Y2K,Y2KP1,DX,RDX,YOLD(K),YOLD(K+1)
C5000 FORMAT(' K=',I4,' Y2K=',E12.4,' Y2KP1=',E12.4,' DX=',E12.4,' RDX='
C    2,E12.4,' YOK=',E12.4,' YOP1=',E12.4)
CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      AK=.1666667*RDX*(Y2KP1-Y2K)
      BK=.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)
C
  550 X=XK-XOLD(K)
      XSQ=X*X
C
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
C
  600 K1=K1+1
      IF(K1.LE.NNEW)GOTO 300
C-----------------------------------------------------------------------
      RETURN
      END
