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

      function tpfc(press,thetae,tgs,d273,rl,qs,pi)

      implicit none

      real(kind=8)  press,thetae,tgs,d273,rl,qs,pi,tpfc
      real(kind=8)  rl461,rl1004,rp,es,fo,t1,tguess,f1,dt
c
c...iteratively extract temperature from equivalent potential
c...temperature.
c
      rl461  = rl/461.5
      rl1004 = rl/1004.
      rp     = thetae/pi
      es     = 611.*dexp(rl461*(d273-1./tgs))
      qs     = 0.622*es/(press-es)
      fo     = tgs*dexp(rl1004*qs/tgs)-rp
      t1     = tgs-0.5*fo
      tguess = tgs
   10 es     = 611.*dexp(rl461*(d273-1./t1))
      qs     = 0.622*es/(press-es)
      f1     = t1*exp(rl1004*qs/t1)-rp
      if (abs(f1).lt..1) goto 20
      dt     = f1*(t1-tguess)/(f1-fo)
      tguess = t1
      fo     = f1
      t1     = t1-dt
      goto 10
c
   20 tpfc = t1
      return
      end function tpfc
