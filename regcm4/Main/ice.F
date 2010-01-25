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

*****************************SUBROUINE ICE*****************************
*                         SIMULATES LAKE ICE                           *
************************************************************************
      subroutine ice(kd,ld,ta,u2,ea,hs,hi,hii,evap,t,depth, precip)
      implicit none

      real(kind=8)  ta,u2,ea,hs,hi,hii,evap,precip
      integer depth
      real(kind=8)  le,kd,ld,lami1,lami2,lams1,lams2,Ki,Ks,khat,Li
      real(kind=8)  surf
      parameter(surf=0.6)
      real(kind=8)  rhoi,atm,qw,rhos,cd,cp,rho,sec,delta
      parameter(lami1=1.5, lami2=20, lams1=6.0, lams2=20.0,
     +          ki=2.3, ks=0.31, rhoi=917.0, atm=950, qw=1.389,
     +          rhos=330.0,  li=334.0E03,  cd=0.001,  cp=1000.0,
     +          rho=1.0, sec=3600, delta=5.67E-08, le=2.47E06)
      real(kind=8)  t(depth,2)
      real(kind=8)  psi,theta,qpen,tf,t0,t1,f0,f1,t2,ds,di,q0
      integer nits

      real(kind=8)  eomb
      external eomb

      real(kind=8)  x,t4,f
      t4(x)=(x+273.15)**4

c ****** g. bates changed air to ta, qpen1 to qpen (4/92)
      f(x)=(-ld+0.97*delta*t4(x)+psi*(eomb(x)-ea)+theta
     +      *(x-ta)-kd)-1./khat*(qpen+tf-x)

CCC   f(x)=(-ld+0.97*delta*t4(x)+psi*(eomb(x)-ea)+theta
CCC  +      *(x-air)-kd)-1/khat*(qpen1+tf-x)

      if ((ta .le. 0.0) .and. (hii .gt. 0.0)) then
        hs = hs+precip*10. / 1000.   ! convert precip(mm) to depth(m)
      end if

      T0=T(1,1)
      tf =  0.0

      khat=(ki*hs+ks*hi)/(ki*ks)
      theta=cp*rho*cd*u2
      psi=le*rho*cd*u2*0.622/atm
      evap= 100.*psi*(eomb(T0)-ea) / (le*rho)
      qpen=kd*0.7*((1.0-dexp(-lams1*hs))/(ks*lams1)+(dexp(-lams1*hs))*
     +     (1.0-dexp(-lami1*hi))/(ki*lami1))+
     +     kd*0.3*((1.0-dexp(-lams2))/(ks*lams2)+(-lams2*hs)*
     +     (1.0-dexp(-lami2*hi))/(ki*lami2))
      kd=kd-qpen

      nits=0
      t1=-50
      f0=f(t0)
      f1=f(t1)
  100 nits=nits+1
      t2=t1-(t1-t0)*f1/(f1-f0)
      if (dabs((t2-t1)/t1) .ge. 0.001) then
        t0=t1
        t1=t2
        f0=f1
        f1=f(t1)
        goto 100
      end if

      t0=t2
      if (t0 .ge. tf) then

        if (hs .gt.0.) then
          ds= sec*((-ld+0.97*delta*t4(tf)+psi*(eomb(tf)-ea)+
     +       theta*(tf-ta)-kd)-1/khat*(t0-tf+qpen)) / (rhos*li)
          if (ds .gt. 0.0) then
            ds=0.0
          end if
          hs=hs+ds
          if (hs .lt. 0) then
            hs=0.0
            t(1,1) = ( hii*t0 + (surf-hii)*t(2,1) ) / surf
          end if
        end if
        if ((hs .eq. 0.) .and. (hii .gt. 0.0)) then
          di= sec*((-ld+0.97*delta*t4(tf)+psi*(eomb(tf)-ea)+
     +       theta*(tf-ta)-kd)-1/khat*(t0-tf+qpen)) / (rhoi*li)
          if (di .gt. 0) then
            di=0.0
          end if
          hi=hi+di
        end if

      else if (t0 .lt. tf) then

        q0=-ld+0.97*delta*t4(t0)+psi*(eomb(t0)-ea)+theta*(t0-ta)-kd
        qpen=kd*0.7*(1.0-dexp(-(lams1*hs+lami1*hi)))+
     +       kd*0.3*(1.0-dexp(-(lams2*hs+lami2*hi)))
        di= sec*(q0-qw-qpen)/(rhoi*li)

        hi=hi+di

      end if

      if (hi .le. 0.01) then
        hi=0.01
        hii=0.0
        hs=0.0
        t(1,1) = ( hi*t0 + (surf-hi)*t(2,1) ) / surf
      else
        hii = hi
        t(1,1) = t0
      end if

      return
      end
