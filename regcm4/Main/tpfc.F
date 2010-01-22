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

      FUNCTION TPFC(PRESS,THETAE,TGS,D273,RL,QS,PI)
      IMPLICIT NONE
      real(kind=8)  PRESS,THETAE,TGS,D273,RL,QS,PI,TPFC
      real(kind=8)  RL461,RL1004,RP,ES,FO,T1,TGUESS,F1,DT
C
C...ITERATIVELY EXTRACT TEMPERATURE FROM EQUIVALENT POTENTIAL
C...TEMPERATURE.
C
      RL461=RL/461.5
      RL1004=RL/1004.
      RP=THETAE/PI
      ES=611.*DEXP(RL461*(D273-1./TGS))
      QS=.622*ES/(PRESS-ES)
      FO=TGS*DEXP(RL1004*QS/TGS)-RP
      T1=TGS-.5*FO
      TGUESS=TGS
   10 ES=611.*DEXP(RL461*(D273-1./T1))
      QS=.622*ES/(PRESS-ES)
      F1=T1*EXP(RL1004*QS/T1)-RP
      IF(ABS(F1).LT..1)GOTO 20
      DT=F1*(T1-TGUESS)/(F1-FO)
      TGUESS=T1
      FO=F1
      T1=T1-DT
      GOTO 10
C
   20 TPFC=T1
      RETURN
      END
