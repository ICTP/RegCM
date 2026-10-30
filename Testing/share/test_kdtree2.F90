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

! Unit test of Share/mod_kdtree2.F90 and Share/mod_kdtree2_priority.F90:
! the searches on random points of the unit sphere are compared to the
! brute force ones, and a tree with many coincident points is built.

program test_kdtree2
  use mod_kdtree2
  use share_test_utils
  implicit none
  integer, parameter :: n = 20000
  integer, parameter :: nq = 2000
  real(kdkind), allocatable, target, dimension(:,:) :: x
  real(kdkind), dimension(3) :: q
  real(kdkind), dimension(2) :: u
  real(kdkind) :: r2
  type(kdtree2), pointer :: tp
  type(kdtree2_result), allocatable, dimension(:) :: rt, rb
  integer :: i, k, m, nn, nf, nfb, nc, nbad, nbadr, nbadc, nunsort
  integer, dimension(3) :: nnv = [4, 8, 32]
  character(len=80) :: msg

  call fixed_seed
  allocate(x(3,n), rt(n), rb(n))
  do i = 1, n
    call random_number(u)
    call sph(u,x(:,i))
  end do
  tp => kdtree2_create(x,sort=.true.,rearrange=.true.)
  do m = 1, size(nnv)
    nbad = 0
    nunsort = 0
    do k = 1, nq
      call random_number(u)
      call sph(u,q)
      nn = nnv(m)
      call kdtree2_n_nearest(tp,q,nn,rt)
      call kdtree2_n_nearest_brute_force(tp,q,nnv(m),rb)
      if ( nn /= nnv(m) ) then
        nbad = nbad + 1
      else if ( maxval(abs(rt(1:nn)%dis-rb(1:nn)%dis)) > 1.0e-14_kdkind ) then
        nbad = nbad + 1
      end if
      if ( any(rt(2:nn)%dis < rt(1:nn-1)%dis) ) nunsort = nunsort + 1
    end do
    write(msg,'(a,i0,a)') 'n_nearest with ', nnv(m), &
      ' points equals brute force'
    call check(nbad == 0,trim(msg))
    write(msg,'(a,i0,a)') 'n_nearest with ', nnv(m), ' points is sorted'
    call check(nunsort == 0,trim(msg))
  end do

  nbadr = 0
  nbadc = 0
  r2 = 0.02_kdkind**2
  do k = 1, nq
    call random_number(u)
    call sph(u,q)
    nc = kdtree2_r_count(tp,q,r2)
    call kdtree2_r_nearest(tp,q,r2,nf,n,rt)
    call kdtree2_r_nearest_brute_force(tp,q,r2,nfb,rb)
    if ( nc /= nfb ) nbadc = nbadc + 1
    if ( nf /= nfb ) then
      nbadr = nbadr + 1
    else if ( nf > 0 ) then
      if ( maxval(abs(rt(1:nf)%dis-rb(1:nf)%dis)) > 1.0e-14_kdkind ) then
        nbadr = nbadr + 1
      end if
    end if
  end do
  call check(nbadr == 0,'r_nearest equals brute force')
  call check(nbadc == 0,'r_count equals brute force')
  call kdtree2_destroy(tp)

  ! Balanced tree.
  tp => kdtree2_create(x,balanced=.true.,sort=.true.,rearrange=.true.)
  nbad = 0
  do k = 1, nq
    call random_number(u)
    call sph(u,q)
    nn = 8
    call kdtree2_n_nearest(tp,q,nn,rt)
    call kdtree2_n_nearest_brute_force(tp,q,8,rb)
    if ( maxval(abs(rt(1:8)%dis-rb(1:8)%dis)) > 1.0e-14_kdkind ) then
      nbad = nbad + 1
    end if
  end do
  call check(nbad == 0,'balanced tree n_nearest equals brute force')
  ! More points in the radius than storage: warns and keeps the nearest.
  call kdtree2_r_nearest(tp,x(:,1),0.05_kdkind**2,nf,5,rt)
  call check(nf >= 5 .and. rt(1)%dis < 1.0e-14_kdkind, &
    'r_nearest with too little storage')
  call kdtree2_destroy(tp)
  deallocate(x)

  ! More coincident points than a bucket: used to recurse forever.
  allocate(x(3,200))
  do i = 1, 199
    x(:,i) = [1.0_kdkind, 0.0_kdkind, 0.0_kdkind]
  end do
  x(:,200) = [0.0_kdkind, 1.0_kdkind, 0.0_kdkind]
  tp => kdtree2_create(x,sort=.true.,rearrange=.true.)
  nn = 4
  call kdtree2_n_nearest(tp,[0.9_kdkind,0.1_kdkind,0.0_kdkind],nn,rt)
  call check(nn == 4 .and. abs(rt(1)%dis-0.02_kdkind) < 1.0e-12_kdkind, &
    'tree with 199 coincident points')
  call kdtree2_destroy(tp)

  call test_end('test_kdtree2')

  contains

  subroutine sph(u,v)
    implicit none
    real(kdkind), dimension(2), intent(in) :: u
    real(kdkind), dimension(3), intent(out) :: v
    real(kdkind) :: lat, lon
    lat = asin(2.0_kdkind*u(1)-1.0_kdkind)
    lon = 2.0_kdkind*acos(-1.0_kdkind)*u(2)
    v = [cos(lat)*sin(lon), sin(lat), cos(lat)*cos(lon)]
  end subroutine sph

end program test_kdtree2
