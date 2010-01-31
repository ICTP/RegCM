!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine spline(nold,xold,yold,y2,nnew,xnew,ynew)
 
      implicit none
!
!*****************************************************************
!                                                                *
!  this is a one-dimensional cubic spline fitting routine        *
!  programed for a small scalar machine.                         *
!                                                                *
!  programer\ z. janjic, yugoslav fed. hydromet. inst., beograd  *
!                                                                *
!  nold - number of given values of the function.  must be ge 3. *
!  xold - locations of the points at which the values of the     *
!         function are given.  must be in ascending order.       *
!  yold - the given values of the function at the points xold.   *
!  y2   - the second derivatives at the points xold.  if natural *
!         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
!         specified.                                             *
!  nnew - number of values of the function to be calculated.     *
!  xnew - locations of the points at which the values of the     *
!         function are calculated.  xnew(k) must be ge xold(1)   *
!         and le xold(nold).                                     *
!  ynew - the values of the function to be calculated.           *
!  p, q - auxiliary vectors of the length nold-2.                *
!                                                                *
!*****************************************************************
!
! Dummy arguments
!
      integer :: nnew , nold
      real(8) , dimension(nold) :: xold , yold , y2
      real(8) , dimension(nnew) :: xnew, ynew
      intent (in) nnew , nold , xnew , xold , yold
      intent (out) ynew
      intent (inout) y2
!
! Local variables
!
      real(8) , dimension(nold-2) :: p , q
      real(8) :: ak , bk , ck , den , dx , dxc , dxl , dxr , dydxl ,    &
               & dydxr , rdx , rtdxc , x , xk , xsq , y2k , y2kp1
      integer :: k , k1 , k2 , kold , noldm1
!-----------------------------------------------------------------------
      noldm1 = nold - 1
!
      dxl = xold(2) - xold(1)
      dxr = xold(3) - xold(2)
      dydxl = (yold(2)-yold(1))/dxl
      dydxr = (yold(3)-yold(2))/dxr
      rtdxc = .5/(dxl+dxr)
!
      p(1) = rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))
      q(1) = -rtdxc*dxr
!
      if ( nold.eq.3 ) then
!-----------------------------------------------------------------------
        k = noldm1
      else
!-----------------------------------------------------------------------
        k = 3
        do
!
          dxl = dxr
          dydxl = dydxr
          dxr = xold(k+1) - xold(k)
          dydxr = (yold(k+1)-yold(k))/dxr
          dxc = dxl + dxr
          den = 1./(dxl*q(k-2)+dxc+dxc)
!
          p(k-1) = den*(6.*(dydxr-dydxl)-dxl*p(k-2))
          q(k-1) = -den*dxr
!
          k = k + 1
          if ( k.ge.nold ) then
            k = noldm1
            exit
          end if
        end do
      end if
      do
!
        y2(k) = p(k-1) + q(k-1)*y2(k+1)
!
        k = k - 1
        if ( k.le.1 ) then
!-----------------------------------------------------------------------
          k1 = 1
          exit
        end if
      end do
!
 100  continue
      xk = xnew(k1)
!
      do k2 = 2 , nold
        if ( xold(k2).gt.xk ) then
          kold = k2 - 1
!
          if ( k1.eq.1 ) go to 200
          if ( k.ne.kold ) go to 200
          go to 300
        end if
      end do
      ynew(k1) = yold(nold)
      go to 400
!
 200  continue
      k = kold
!
      y2k = y2(k)
      y2kp1 = y2(k+1)
      dx = xold(k+1) - xold(k)
      rdx = 1./dx
      ak = .1666667*rdx*(y2kp1-y2k)
      bk = .5*y2k
      ck = rdx*(yold(k+1)-yold(k)) - .1666667*dx*(y2kp1+y2k+y2k)
!
 300  continue
      x = xk - xold(k)
      xsq = x*x
!
      ynew(k1) = ak*xsq*x + bk*xsq + ck*x + yold(k)
!
 400  continue
      k1 = k1 + 1
      if ( k1.le.nnew ) go to 100
!-----------------------------------------------------------------------
      end subroutine spline
