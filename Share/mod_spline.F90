!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! This file is part of ICTP RegCM.
!
! ICTP RegCM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ICTP RegCM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! how these work: consider 1D interpolation, for which you use the
! SPLINE and SPLINT.  The interpolation is performed in a two step
! process.  Step 1 is that you compute the 2nd derivatives of the Y
! array data at each X array data point.  This is accomplished using
! a single call to subroutine SPLINE.  Step 2 is the interpolation
! itself, which is accomplished by calling SPLINT.  Multiple
! successive calls to SPLINT can be made to interpolate the Y array
! data to obtain an interpolated y value at various x values.

module mod_spline

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  implicit none

  private

  public :: splini , spline , splint
  public :: splie2 , splin2

  contains
    !
    ! Calculate 2nd derivatives of cubic spline interp function
    ! Adapted from numerical recipes by press et al
    !
    subroutine spline(x,y,yp1,ypn,y2)
      implicit none
      ! Arrays of tabulated function in ascending order by x with y = f(x)
      real(rkx) , dimension(:) , intent(in)  :: x , y
      ! Specified derivatives at x(1) and x(n)
      ! Values > 1E30 signals second derivative zero
      real(rkx) , intent(in) :: yp1 , ypn
      ! Output array of second derivatives
      real(rkx) , dimension(:) , intent(out) :: y2

      ! SPLINE use: given an 1D array of X data and an array of Y data,
      ! both of length N, this routine computes the 2nd derivatives, Y2 at
      ! each X data point.  The user needs to specify the values of YP1
      ! and YP2, which flags how the Y2 are computed at the edges.  For
      ! natural spline fitting (recommended), set YP1 and YPN to numbers
      ! greater than 1.0E+30.

      ! this routine called once, prior to using routine SPLINT, as a set
      ! up for using routine SPLINT, which performs the actual
      ! interpolation

      ! IMPORTANT NOTE: the X data values in array X must be in ascending
      ! order or the interpolation will fail

      integer(ik4) :: n , i , k
      real(rkx) :: p , qn , sig , un
      real(rkx) , dimension(:) , allocatable :: u

      ! if YP1>1.0E+30 use natural spline, otherwise estimate Y2 at the
      ! first point

      n = size(x,1)
      allocate(u(n-1))

      if ( yp1 > 0.99e30_rkx ) then
        y2(1) = d_zero
        u(1) = d_zero
      else
        y2(1) = -d_half
        u(1) = (d_three/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if

      ! store intermediate values of terms in the expansion series

      do i = 2 , n - 1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1) + d_two
        y2(i) = (sig-d_one)/p
        u(i) = (d_six*((y(i+1)-y(i))/(x(i+1)-x(i))- &
                     (y(i)-y(i-1))/(x(i)-x(i-1))) / &
                     (x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      ! if YPN>1.0E+30 use natural spline, otherwise estimate Y2 at the
      ! last point point

      if ( ypn > 0.99e30_rkx ) then
        qn = d_zero
        un = d_zero
      else
        qn = d_half
        un = (d_three/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+d_one)

      ! compute the Y2 from the 2nd order expansion series

      do k = n - 1 , 1 , -1
        y2(k) = y2(k)*y2(k+1) + u(k)
      end do

      deallocate(u)
    end subroutine spline
    !
    ! Integrate cubic spline function from xa(1) to x
    !
    subroutine splini(xa,ya,y2a,x,yi)
      implicit none
      ! Arrays of tabulated function in ascending order by xa with ya = f(xa)
      real(rkx) , dimension(:) , intent(in) :: xa , ya
      ! Array of second derivatives
      real(rkx) , dimension(:) , intent(in) :: y2a
      ! Ascissa endpoint of integration
      real(rkx) , intent(in) :: x
      ! Output value
      real(rkx) , intent(out) :: yi

      real(rkx) :: a , a2 , b , b2 , h , xx
      integer(ik4) :: n , khi , klo

      n = size(xa,1)
      yi = d_zero
      klo = 1
      khi = 2
      do while ( x > xa(klo) .and. khi <= n )
        xx = x
        if ( khi<n ) xx = min(x,xa(khi))
        h = xa(khi) - xa(klo)
        a = (xa(khi)-xx)/h
        b = (xx-xa(klo))/h
        a2 = a*a
        b2 = b*b
        yi = yi + ((d_one-a2)*ya(klo)/d_two+b2*ya(khi)/d_two+ &
                   ((-(d_one+a2*a2)/d_four+a2/d_two)*y2a(klo)+ &
                   (b2*b2/d_four-b2/d_two)*y2a(khi))*h*h/d_six)*h
        klo = klo + 1
        khi = khi + 1
      end do
    end subroutine splini
    !
    ! Calculate cubic spline interp value
    ! Adapted from numerical recipes by press et al.
    !
    subroutine splint(xa,ya,y2a,x,y)
      implicit none

      ! Arrays of tabulated function values in ascending xa order
      real(rkx) , dimension(:) , intent(in) :: xa , ya
      ! Arrays of second derivatives
      real(rkx) , dimension(:) , intent(in) :: y2a
      ! Abscissa of interpolation
      real(rkx) , intent(in) :: x
      ! Output value
      real(rkx) , intent(out) :: y

      real(rkx) :: a , b , h
      integer(ik4) :: n , k , khi , klo

      ! SPLINT use: given an 1D array of XA data, an array of YA data, and
      ! an array of the 2nd derivatives Y2A, all of length N, this routine
      ! performs cubic spline interpolation, returning the interpolated
      ! value Y at the user input value X.  The Y2A are computed in
      ! routine SPLINE, which is called once before calling SPLINT.

      ! IMPORTANT NOTE: the X data values in array X must be in ascending
      ! order or the interpolation will fail

      n = size(xa,1)
      klo = 1
      khi = n
      ! determine the indices of array XA that bracket the input X value
      do while ( khi-klo > 1 )
        k = (khi+klo)/2
        if ( xa(k)>x ) then
          khi = k
        else
          klo = k
        end if
      end do

      ! determine the finite difference along the X dimension
      h = xa(khi) - xa(klo)

      if ( h==0 ) then
        y = 1.0e-30_rkx
        return
      end if

      ! interpolate

      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a*a*a-a)*y2a(klo) + &
                                   (b*b*b-b)*y2a(khi))*h*h/d_six
    end subroutine splint

    subroutine splie2(x1a,x2a,ya,m,n,y2a)
      implicit none
      !     SPLIE2 use: given an array of X1A data of length M, and an array
      !     of X2A data of length N, this routine computes the 2nd
      !     derivatives, Y2A, at each X1A,X2A data point.  Thus Y2A has
      !     dimensions Y2A(M,N).  Natural spline fitting is assumed.
      !     this routine called once, prior to using routine SPLIN2, as a set
      !     up for using routine SPLIN2, which performs the actual
      !     interpolation.
      !     Uses routines: SPLINE
      !     IMPORTANT NOTE: the X1A and X2A data values must both be in
      !     ascending order or the interpolation will fail
      real(rkx) , dimension(:) , intent(in) :: x1a , x2a
      real(rkx) , dimension(:,:) , intent(in) :: ya
      real(rkx) , dimension(:,:) , intent(out) :: y2a
      integer(ik4) :: m , n , j , k

      real(rkx) , allocatable , dimension(:) :: ytmp , y2tmp

      m = size(x1a,1)
      n = size(x2a,1)

      allocate(ytmp(n), y2tmp(n))
      do j = 1 , m
        do k = 1 , n
          ytmp(k) = ya(j,k)
        end do
        call spline(x2a,ytmp,1.e30_rkx,1.e30_rkx,y2tmp)
        do k = 1 , n
          y2a(j,k) = y2tmp(k)
        end do
      end do
      deallocate(ytmp,y2tmp)
    end subroutine splie2

    subroutine splin2(x1a,x2a,ya,y2a,x1,x2,y)
      implicit none
      ! SPLIN2 use: given an array of X1A data of length M, an array of
      ! X2A data of length N, and an array of 2nd derivatives, Y2A at each
      ! X1A,X2A data point, dimensioned Y2A(M,N), this routine performs 2D
      ! interpolation, returning the interpolated value Y at user input
      ! values X1 and X2.  Natural spline fitting is assumed.
      ! Uses routines: SPLINT, SPLINE
      ! IMPORTANT NOTE: the X1A and X2A data values must both be in
      ! ascending order or the interpolation will fail
      real(rkx) , dimension(:) , intent(in) :: x1a
      real(rkx) , dimension(:) , intent(in) :: x2a
      real(rkx) , dimension(:,:) , intent(in) :: ya
      real(rkx) , dimension(:,:) , intent(in) :: y2a
      real(rkx) , intent(in) :: x1
      real(rkx) , intent(in) :: x2
      real(rkx) , intent(out) :: y
      integer(ik4) :: j , k , m , n

      real(rkx) , allocatable , dimension(:) :: ytmp , y2tmp , yytmp

      m = size(x1a,1)
      n = size(x2a,1)
      allocate(ytmp(n), y2tmp(n))

      do j = 1 , m
        do k = 1 , n
          ytmp(k) = ya(j,k)
          y2tmp(k) = y2a(j,k)
        end do
        call splint(x2a,ytmp,y2tmp,x2,yytmp(j))
      end do
      call spline(x1a,yytmp,1.e30_rkx,1.e30_rkx,y2tmp)
      call splint(x1a,yytmp,y2tmp,x1,y)
    end subroutine splin2

end module mod_spline
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
