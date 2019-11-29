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

module mod_spbarcoord

  !
  ! Implementation of Spherical Barycentric Coordinates as described in:
  !
  !   Langer, Belyaev and Seidel
  !   Eurographics Symposium on Geometry Processing (2006)
  !   Polthier - Sheffer Editors
  !
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use iso_c_binding

  implicit none

  private

  public :: spherical_barycentric

  type , bind(C) :: vpoint
    integer(c_int) :: idx
    real(c_double) , dimension(3) :: v
  end type vpoint

  real(rk8) , dimension(3) :: centroid

  interface
    subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
      import
      type(c_ptr) , value       :: array
      integer(c_size_t) , value :: elem_count
      integer(c_size_t) , value :: elem_size
      type(c_funptr) , value    :: compare
    end subroutine qsort !standard C library qsort
  end interface

  contains

  subroutine spherical_barycentric(np,p,v,lambda)
    implicit none
    integer(ik4) , intent(in) :: np
    real(rk8) , intent(in) , dimension(3) :: p
    real(rk8) , intent(in) , dimension(3,np) :: v
    real(rk8) , intent(out) , dimension(np) :: lambda
    type(vpoint) , target , dimension(np) :: voc
    real(rk8) , dimension(np) :: alpha
    real(rk8) , dimension(np) :: theta
    real(rk8) , dimension(np) :: tansum
    real(rk8) :: norm
    integer(ik4) :: i

    ! Compute ordering centroid

    centroid = d_zero
    do i = 1 , np
      centroid(1) = centroid(1) + v(1,i)
      centroid(2) = centroid(2) + v(2,i)
      centroid(3) = centroid(3) + v(3,i)
    end do
    centroid = centroid / real(np,rk8)

    call set_clockwise_order
    call compute_angles
    tansum(1) = tan(alpha(np)*0.5_rk8) + tan(alpha(1)*0.5_rk8)
    do i = 2 , np
      tansum(i) = tan(alpha(i-1)*0.5_rk8) + tan(alpha(i)*0.5_rk8)
    end do
    norm = 0.0_rk8
    do i = 1 , np
      norm = norm + tansum(i) / tan(theta(i))
    end do
    do i = 1 , np
      lambda(voc(i)%idx) = ( tansum(i) / sin(theta(i)) ) / norm
    end do

    ! Double check sum weights is one
    norm = sum(lambda(1:np))-1.0_rk8
    if ( abs(norm) > epsilon(1.0_rk8) ) then
      lambda(:) = lambda(:) - norm * (1.0_rk8-lambda(:))
    end if

    contains

      subroutine set_clockwise_order
        implicit none
        integer(c_size_t) :: lnp , lsize
        integer(ik4) :: n
        do n = 1 , np
          voc(n)%idx = n
          voc(n)%v(:) = v(:,n)
        end do
        lnp = np
        lsize = sizeof(voc(1))
        call qsort(c_loc(voc(1)),lnp,lsize,c_funloc(compare))
      end subroutine set_clockwise_order

      pure real(rk8) function norma2(x) result(a)
        implicit none
        real(rk8) , dimension(3) , intent(in) :: x
        a = (x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      end function norma2

      pure real(rk8) function norma(x) result(a)
        implicit none
        real(rk8) , dimension(3) , intent(in) :: x
        a = sqrt(norma2(x))
      end function norma

      pure real(rk8) function dotprod(x,y) result(a)
        implicit none
        real(rk8) , dimension(3) , intent(in) :: x , y
        a = (x(1)*y(1) + x(2)*y(2) + x(3)*y(3))
      end function dotprod

      pure real(rk8) function angle_between(x,y) result(a)
        implicit none
        real(rk8) , dimension(3) , intent(in) :: x , y
        a = max(-1.0_rk8,min(1.0_rk8,dotprod(x,y) / (norma(x)*norma(y))))
        a = acos(a)
      end function angle_between

      subroutine vecprod(x,y,z)
        implicit none
        real(rk8) , dimension(3) , intent(in) :: x , y
        real(rk8) , dimension(3) , intent(out) :: z
        z(1) = x(2)*y(3) - x(3)*y(2)
        z(2) = x(3)*y(1) - x(1)*y(3)
        z(3) = x(1)*y(2) - x(2)*y(1)
      end subroutine vecprod

      subroutine compute_angles
        implicit none
        real(rk8) , dimension(3) :: vp1 , vp2
        do i = 1 , np
          theta(i) = angle_between(p,voc(i)%v)
        end do
        call vecprod(p,voc(np)%v,vp1)
        call vecprod(p,voc(1)%v,vp2)
        alpha(np) = angle_between(vp1,vp2)
        do i = 1 , np-1
          call vecprod(p,voc(i+1)%v,vp1)
          call vecprod(p,voc(i)%v,vp2)
          alpha(i) = angle_between(vp1,vp2)
        end do
      end subroutine compute_angles

  end subroutine spherical_barycentric

  integer(c_int) function compare(x1,x2) result(res) bind(C)
    implicit none
    type(vpoint) , intent(in) :: x1 , x2
    real(rk8) , dimension(3) :: p1 , p2 , n
    real(rk8) :: xres

    call vecdiff(centroid,x1%v,p1)
    call vecdiff(centroid,x2%v,p2)

    n(1) = p1(2)*p2(3) - p1(3)*p2(2)
    n(2) = p1(3)*p2(1) - p1(1)*p2(3)
    n(3) = p1(1)*p2(2) - p1(2)*p2(1)

    xres = ((n(1)*centroid(1) + n(2)*centroid(2) + n(3)*centroid(3)))

    if ( xres < 0 ) then
      res = -1
    else if ( xres > 0 ) then
      res = 1
    else
      res = 0
    end if

    contains

    subroutine vecdiff(a,b,c)
      implicit none
      real(rk8) , dimension(3) , intent(in) :: a , b
      real(rk8) , dimension(3) , intent(out) :: c
      c(1) = a(1) - b(1)
      c(2) = a(2) - b(2)
      c(3) = a(3) - b(3)
    end subroutine vecdiff

  end function compare

end module mod_spbarcoord

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
