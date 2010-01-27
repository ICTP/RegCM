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

      subroutine spline(nold,xold,yold,y2,nnew,xnew,ynew,p,q)

      implicit none
C     ******************************************************************
C     *                                                                *
C     *  this is a one-dimensional cubic spline fitting routine        *
C     *  programed for a small scalar machine.                         *
C     *                                                                *
C     *  programer\ z. janjic, yugoslav fed. hydromet. inst., beograd  *
C     *                                                                *
C     *  nold - number of given values of the function.  must be ge 3. *
C     *  xold - locations of the points at which the values of the     *
C     *         function are given.  must be in ascending order.       *
C     *  yold - the given values of the function at the points xold.   *
C     *  y2   - the second derivatives at the points xold.  if natural *
C     *         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
C     *         specified.                                             *
C     *  nnew - number of values of the function to be calculated.     *
C     *  xnew - locations of the points at which the values of the     *
C     *         function are calculated.  xnew(k) must be ge xold(1)   *
C     *         and le xold(nold).                                     *
C     *  ynew - the values of the function to be calculated.           *
C     *  p, q - auxiliary vectors of the length nold-2.                *
C     *                                                                *
C     ******************************************************************
      integer nold,nnew
      integer jtb
      parameter (jtb=150)
      real(kind=8)  xold(jtb),yold(jtb),y2(jtb),p(jtb),q(jtb),xnew(jtb)
     &       ,ynew(jtb)
      integer noldm1,k,k1,k2,kold
      real(kind=8)  dxl,dxr,dydxl,dydxr,rtdxc
      real(kind=8)  dxc,den,xk,y2k,y2kp1,dx,rdx,ak,bk,ck,x,xsq
C-----------------------------------------------------------------------
      noldm1=nold-1
C
      dxl=xold(2)-xold(1)
      dxr=xold(3)-xold(2)
      dydxl=(yold(2)-yold(1))/dxl
      dydxr=(yold(3)-yold(2))/dxr
      rtdxc=.5/(dxl+dxr)
C
      p(1)=rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1)=-rtdxc*dxr
C
      if(nold.eq.3)goto 700
C-----------------------------------------------------------------------
      k=3
C
  100 dxl=dxr
      dydxl=dydxr
      dxr=xold(k+1)-xold(k)
      dydxr=(yold(k+1)-yold(k))/dxr
      dxc=dxl+dxr
      den=1./(dxl*q(k-2)+dxc+dxc)
C
      p(k-1)=den*(6.*(dydxr-dydxl)-dxl*p(k-2))
      q(k-1)=-den*dxr
C
      k=k+1
      if(k.lt.nold)goto 100
C-----------------------------------------------------------------------
  700 k=noldm1
C
  200 y2(k)=p(k-1)+q(k-1)*y2(k+1)
C
      k=k-1
      if(k.gt.1)goto 200
C-----------------------------------------------------------------------
      k1=1
C
  300 xk=xnew(k1)
C
      do 400 k2=2,nold
        if(xold(k2).le.xk)goto 400
        kold=k2-1
        goto 450
  400 continue
      ynew(k1)=yold(nold)
      goto 600
C
  450 if(k1.eq.1)goto 500
      if(k.eq.kold)goto 550
C
  500 k=kold
C
      y2k=y2(k)
      y2kp1=y2(k+1)
      dx=xold(k+1)-xold(k)
      rdx=1./dx
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C     write(6,5000) k,y2k,y2kp1,dx,rdx,yold(k),yold(k+1)
C5000 format(' k=',i4,' y2k=',e12.4,' y2kp1=',e12.4,' dx=',e12.4,' rdx='
C    2,e12.4,' yok=',e12.4,' yop1=',e12.4)
Caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
      ak=.1666667*rdx*(y2kp1-y2k)
      bk=.5*y2k
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)
C
  550 x=xk-xold(k)
      xsq=x*x
C
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)
C
  600 k1=k1+1
      if(k1.le.nnew)goto 300
C-----------------------------------------------------------------------
      return
      end subroutine spline
