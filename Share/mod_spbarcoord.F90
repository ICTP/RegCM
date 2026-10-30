!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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

  implicit none

  private

  public :: spherical_barycentric

  type vpoint
    integer(ik4) :: idx
    real(rk8), dimension(3) :: v
  end type vpoint

  real(rk8), parameter :: notzero = epsilon(1.0_rk8)

  contains

  subroutine spherical_barycentric(np,p,v,lambda)
    implicit none
    integer(ik4), intent(in) :: np
    real(rk8), intent(in), dimension(3) :: p
    real(rk8), intent(in), dimension(3,np) :: v
    real(rk8), intent(out), dimension(np) :: lambda
    type(vpoint), dimension(np) :: voc
    real(rk8), dimension(np) :: alpha
    real(rk8), dimension(np) :: theta
    real(rk8), dimension(np) :: tansum
    real(rk8) :: norm, xs
    integer(ik4) :: i

    if ( np == 0 ) return

    call set_cyclic_order
    call compute_angles
    tansum(1) = tan(alpha(np)*0.5_rk8) + tan(alpha(1)*0.5_rk8)
    do i = 2, np
      tansum(i) = tan(alpha(i-1)*0.5_rk8) + tan(alpha(i)*0.5_rk8)
    end do
    norm = 0.0_rk8
    do i = 1, np
      xs = sin(theta(i))
      norm = norm + tansum(i) / sign(max(xs,notzero),xs)
    end do
    norm = max(norm,notzero)
    do i = 1, np
      xs = sin(theta(i))
      lambda(voc(i)%idx) = (tansum(i)/sign(max(xs,notzero),xs))/norm
    end do

    ! Double check sum weights is one
    norm = sum(lambda(1:np))
    if ( norm > notzero ) lambda(1:np) = lambda(1:np) / norm

    contains

      !
      ! Order the vertices by their angle around the query point p in the
      ! plane tangent to the sphere at p, so that they form a polygon
      ! around p (star-shaped as seen from p).
      !
      subroutine set_cyclic_order
        implicit none
        real(rk8), dimension(3) :: e1, e2, ax
        real(rk8), dimension(np) :: phi
        real(rk8) :: phin
        type(vpoint) :: vn
        integer(ik4) :: n, m
        ! Orthonormal basis (e1,e2) of the tangent plane at p, built from
        ! the coordinate axis least aligned with p
        ax = d_zero
        ax(minloc(abs(p),1)) = 1.0_rk8
        call vecprod(ax,p,e1)
        e1 = e1 / max(norma(e1),notzero)
        call vecprod(p,e1,e2)
        do n = 1, np
          voc(n)%idx = n
          voc(n)%v(:) = v(:,n)
          phi(n) = atan2(dotprod(v(:,n),e2),dotprod(v(:,n),e1))
        end do
        ! Insertion sort by angle: np is small
        do n = 2, np
          vn = voc(n)
          phin = phi(n)
          m = n - 1
          do while ( m >= 1 )
            if ( phi(m) <= phin ) exit
            voc(m+1) = voc(m)
            phi(m+1) = phi(m)
            m = m - 1
          end do
          voc(m+1) = vn
          phi(m+1) = phin
        end do
      end subroutine set_cyclic_order

      pure real(rk8) function norma2(x) result(a)
        implicit none
        real(rk8), dimension(3), intent(in) :: x
        a = (x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
      end function norma2

      pure real(rk8) function norma(x) result(a)
        implicit none
        real(rk8), dimension(3), intent(in) :: x
        a = sqrt(norma2(x))
      end function norma

      pure real(rk8) function dotprod(x,y) result(a)
        implicit none
        real(rk8), dimension(3), intent(in) :: x, y
        a = (x(1)*y(1) + x(2)*y(2) + x(3)*y(3))
      end function dotprod

      pure real(rk8) function angle_between(x,y) result(a)
        implicit none
        real(rk8), dimension(3), intent(in) :: x, y
        a = max(-1.0_rk8,min(1.0_rk8,dotprod(x,y) / &
                  max((norma(x)*norma(y)),notzero)))
        a = acos(a)
      end function angle_between

      subroutine vecprod(x,y,z)
        implicit none
        real(rk8), dimension(3), intent(in) :: x, y
        real(rk8), dimension(3), intent(out) :: z
        z(1) = x(2)*y(3) - x(3)*y(2)
        z(2) = x(3)*y(1) - x(1)*y(3)
        z(3) = x(1)*y(2) - x(2)*y(1)
      end subroutine vecprod

      subroutine compute_angles
        implicit none
        real(rk8), dimension(3) :: vp1, vp2
        integer(ik4) :: i
        do i = 1, np
          theta(i) = angle_between(p,voc(i)%v)
        end do
        call vecprod(p,voc(np)%v,vp1)
        call vecprod(p,voc(1)%v,vp2)
        alpha(np) = angle_between(vp1,vp2)
        do i = 1, np-1
          call vecprod(p,voc(i+1)%v,vp1)
          call vecprod(p,voc(i)%v,vp2)
          alpha(i) = angle_between(vp1,vp2)
        end do
      end subroutine compute_angles

  end subroutine spherical_barycentric

end module mod_spbarcoord

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
