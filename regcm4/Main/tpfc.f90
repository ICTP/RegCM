!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      function tpfc(press,thetae,tgs,d273,rl,qs,pi)
 
      use mod_constants , only : rcpd , ep2 , rwat
      implicit none
!
! Dummy arguments
!
      real(8) :: d273 , pi , press , qs , rl , tgs , thetae
      real(8) :: tpfc
      intent (in) d273 , pi , press , rl , tgs , thetae
      intent (inout) qs
!
! Local variables
!
      real(8) :: dt , es , f1 , fo , rlocpd , rlorw , rp , t1 , tguess
!
!...iteratively extract temperature from equivalent potential
!...temperature.
!
      rlorw = rl/rwat
      rlocpd = rl*rcpd
      rp = thetae/pi
      es = 611.*dexp(rlorw*(d273-1./tgs))
      qs = ep2*es/(press-es)
      fo = tgs*dexp(rlocpd*qs/tgs) - rp
      t1 = tgs - 0.5*fo
      tguess = tgs
 100  es = 611.*dexp(rlorw*(d273-1./t1))
      qs = ep2*es/(press-es)
      f1 = t1*exp(rlocpd*qs/t1) - rp
      if ( abs(f1).lt..1 ) then
!
        tpfc = t1
      else
        dt = f1*(t1-tguess)/(f1-fo)
        tguess = t1
        fo = f1
        t1 = t1 - dt
        go to 100
      end if
      end function tpfc
