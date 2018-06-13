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

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use iso_c_binding

  implicit none

  private

  public :: spherical_barycentric

  type , bind(C) :: vpoint
    integer(ik4) :: idx
    real(rkx) , dimension(3) :: v
  end type vpoint

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
    real(rkx) , intent(in) , dimension(3) :: p
    real(rkx) , intent(in) , dimension(3,np) :: v
    real(rkx) , intent(out) , dimension(np) :: lambda
    type(vpoint) , target , dimension(np) :: voc
    real(rkx) , dimension(np) :: alpha
    real(rkx) , dimension(np) :: theta
    real(rkx) , dimension(np) :: tansum
    real(rkx) :: norm
    integer(ik4) :: i , im1

    call set_clockwise_order
    call compute_angles
    do i = 1 , np
      im1 = i - 1
      if ( im1 == 0 ) im1 = np
      tansum(i) = tan(alpha(im1)*0.5_rkx) + tan(alpha(i)*0.5_rkx)
    end do
    norm = 0.0_rkx
    do i = 1 , np
      norm = norm + tansum(i) / tan(theta(i))
    end do
    do i = 1 , np
      lambda(voc(i)%idx) = ( tansum(i) / sin(theta(i)) ) / norm
    end do

    ! Double check sum weights is one
    norm = sum(lambda(1:np-1))
    lambda(np) = d_one - norm

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

      integer(c_int) function compare(x1,x2) result(res) bind(C)
        implicit none
        type(vpoint) , intent(in) :: x1 , x2
        if ( x1%v(3) > x2%v(3) ) then
          res = -1
          return
        end if
        if ( x1%v(3) < x2%v(3) ) then
          res = 1
          return
        end if
        if ( x1%v(1) > x2%v(1) ) then
          res = -1
        end if
        res = 1
      end function compare

      pure real(rkx) function norma(x) result(a)
        implicit none
        real(rkx) , dimension(3) , intent(in) :: x
        a = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      end function norma

      pure real(rkx) function dotprod(x,y) result(a)
        implicit none
        real(rkx) , dimension(3) , intent(in) :: x , y
        a = (x(1)*y(1) + x(2)*y(2) + x(3)*y(3))
      end function dotprod

      pure real(rkx) function angle_between(x,y) result(a)
        implicit none
        real(rkx) , dimension(3) , intent(in) :: x , y
        a = dotprod(x,y) / (norma(x)*norma(y))
        a = acos(a)
      end function angle_between

      subroutine compute_angles
        implicit none
        real(rkx) :: pn , tmp1 , di , pn2 , dip , vin , vipn , vvipn2
        integer(ik4) :: ip
        do i = 1 , np
          theta(i) = angle_between(p,voc(i)%v)
        end do
        do i = 1 , np
          ip = i+1
          if ( ip > np ) ip = 1
          pn = norma(p)
          pn2 = pn *pn
          di = pn * tan(theta(i))
          dip = pn * tan(theta(ip))
          tmp1 = angle_between(voc(i)%v,voc(ip)%v)
          vin = sqrt(di*di+pn2)
          vipn = sqrt(dip*dip+pn2)
          vvipn2 = (vipn*sin(tmp1))**2+(vin-vipn*cos(tmp1))**2
          alpha(i) = acos(0.5_rkx*(dip*dip+di*di-vvipn2)/(dip*di))
        end do
      end subroutine compute_angles

  end subroutine spherical_barycentric

end module mod_spbarcoord

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
