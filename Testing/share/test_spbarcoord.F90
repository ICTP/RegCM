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

! Unit test of Share/mod_spbarcoord.F90: for a query point inside a convex
! spherical polygon of 3 to 8 vertices given in random order, the weights
! must be positive, sum to one, reproduce the query point and not depend
! on the order of the vertices.

program test_spbarcoord
  use mod_intkinds
  use mod_realkinds
  use mod_spbarcoord
  use share_test_utils
  implicit none
  integer(ik4), parameter :: ntest = 20000
  real(rk8), parameter :: d2r = acos(-1.0_rk8)/180.0_rk8
  real(rk8), dimension(3) :: p, q
  real(rk8), dimension(3,8) :: v, w
  real(rk8), dimension(8) :: lambda, mu
  integer(ik4), dimension(8) :: perm
  real(rk8) :: u, lat0, lon0, ang, r, err, errmax, permax, summax
  integer(ik4) :: it, np, i, nneg

  call fixed_seed
  errmax = 0.0_rk8
  permax = 0.0_rk8
  summax = 0.0_rk8
  nneg = 0
  do it = 1, ntest
    call random_number(u)
    np = 3 + int(u*6.0_rk8)
    call random_number(u)
    lat0 = (u-0.5_rk8)*170.0_rk8
    call random_number(u)
    lon0 = u*360.0_rk8
    ! Vertices in angular sectors around the centre, at nearly the same
    ! distance: a convex polygon around the centre.
    do i = 1, np
      call random_number(u)
      ang = (real(i-1,rk8)+0.3_rk8*u)*360.0_rk8/real(np,rk8)
      call random_number(u)
      r = 0.8_rk8 + 0.2_rk8*u
      call ll(lat0+r*sin(ang*d2r),lon0+r*cos(ang*d2r)/cos(lat0*d2r),v(:,i))
    end do
    call shuffle(v,np)
    call random_number(u)
    call ll(lat0+(u-0.5_rk8)*0.2_rk8,lon0,p)
    call spherical_barycentric(np,p,v(:,1:np),lambda(1:np))
    if ( any(lambda(1:np) < 0.0_rk8) ) nneg = nneg + 1
    summax = max(summax,abs(sum(lambda(1:np))-1.0_rk8))
    q = matmul(v(:,1:np),lambda(1:np))
    q = q/sqrt(sum(q*q))
    ! Distance from the cross product: acos loses precision near 1.
    err = sqrt((q(2)*p(3)-q(3)*p(2))**2 + (q(3)*p(1)-q(1)*p(3))**2 + &
               (q(1)*p(2)-q(2)*p(1))**2)*6371.0_rk8
    errmax = max(errmax,err)
    ! Same points in another order: the same weight for each point.
    do i = 1, np
      perm(i) = i
    end do
    call shuffle_index(perm,np)
    do i = 1, np
      w(:,i) = v(:,perm(i))
    end do
    call spherical_barycentric(np,p,w(:,1:np),mu(1:np))
    do i = 1, np
      permax = max(permax,abs(mu(i)-lambda(perm(i))))
    end do
  end do
  call check(nneg == 0,'weights are positive inside a convex polygon')
  call check(summax < 1.0e-12_rk8,'weights sum to one')
  call check(errmax < 1.0e-6_rk8,'weights reproduce the query point (1 mm)')
  call check(permax < 1.0e-12_rk8,'weights do not depend on the order')

  ! A query point on the bisector of two of four symmetric vertices.
  call ll(1.0_rk8,0.0_rk8,v(:,1))
  call ll(0.0_rk8,1.0_rk8,v(:,2))
  call ll(-1.0_rk8,0.0_rk8,v(:,3))
  call ll(0.0_rk8,-1.0_rk8,v(:,4))
  call ll(0.0_rk8,0.0_rk8,p)
  call spherical_barycentric(4,p,v(:,1:4),lambda(1:4))
  call check(all(abs(lambda(1:4)-0.25_rk8) < 1.0e-12_rk8), &
    'centre of a symmetric square has equal weights')

  call test_end('test_spbarcoord')

  contains

  subroutine ll(lat,lon,x)
    implicit none
    real(rk8), intent(in) :: lat, lon
    real(rk8), dimension(3), intent(out) :: x
    x = [cos(lat*d2r)*sin(lon*d2r), sin(lat*d2r), cos(lat*d2r)*cos(lon*d2r)]
  end subroutine ll

  subroutine shuffle(a,n)
    implicit none
    real(rk8), dimension(:,:), intent(inout) :: a
    integer(ik4), intent(in) :: n
    integer(ik4) :: k, j
    real(rk8), dimension(3) :: t
    real(rk8) :: uu
    do k = n, 2, -1
      call random_number(uu)
      j = 1 + int(uu*real(k,rk8))
      t = a(:,k)
      a(:,k) = a(:,j)
      a(:,j) = t
    end do
  end subroutine shuffle

  subroutine shuffle_index(a,n)
    implicit none
    integer(ik4), dimension(:), intent(inout) :: a
    integer(ik4), intent(in) :: n
    integer(ik4) :: k, j, t
    real(rk8) :: uu
    do k = n, 2, -1
      call random_number(uu)
      j = 1 + int(uu*real(k,rk8))
      t = a(k)
      a(k) = a(j)
      a(j) = t
    end do
  end subroutine shuffle_index

end program test_spbarcoord
