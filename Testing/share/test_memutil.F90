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

! Unit test of Share/mod_memutil.F90 and Share/mod_space.F90: bounds and
! initial values of getmem, release of any element of the pool list,
! reallocation of an associated pointer, zero size arrays, assignpnt and
! remappnt4, and a valgrind clean memory_destroy.

program test_memutil
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_date, only : rcm_time_and_date
  use share_test_utils
  implicit none
  real(rk8), pointer, contiguous, dimension(:,:) :: a => null()
  real(rk8), pointer, contiguous, dimension(:,:) :: b => null()
  real(rk8), pointer, contiguous, dimension(:,:) :: c => null()
  real(rk8), pointer, contiguous, dimension(:,:) :: z => null()
  real(rk8), pointer, contiguous, dimension(:,:) :: s2 => null()
  real(rk8), pointer, contiguous, dimension(:,:,:) :: t3 => null()
  real(rk4), pointer, contiguous, dimension(:,:,:,:) :: r4 => null()
  real(rk4), pointer, contiguous, dimension(:,:,:,:) :: m4 => null()
  integer(ik4), pointer, contiguous, dimension(:) :: iv => null()
  integer(ik2), pointer, contiguous, dimension(:,:,:) :: sv => null()
  logical, pointer, contiguous, dimension(:,:,:,:,:) :: l5 => null()
  type(rcm_time_and_date), pointer, contiguous, dimension(:) :: tv => null()
  integer(ik4) :: k

  call memory_init

  call getmem(a,1,3,0,4,'a')
  call check(all(lbound(a) == [1,0]) .and. all(ubound(a) == [3,4]), &
    'getmem bounds')
  call check(all(a == 0.0_rk8),'real arrays start at zero')
  call getmem(iv,1,5,'iv')
  call check(all(iv == -1),'integer arrays start at -1')
  call getmem(sv,1,2,1,2,1,2,'sv')
  call check(all(sv == -1_ik2),'short integer arrays start at -1')
  call getmem(l5,1,2,1,2,1,2,1,2,1,2,'l5')
  call check(.not. any(l5),'logical arrays start false')
  call getmem(tv,1,3,'tv')
  call check(size(tv) == 3,'date arrays')

  call getmem(b,1,3,0,4,'b')
  call getmem(c,-2,2,1,2,'c')
  call relmem(b)
  call check(.not. associated(b),'release in the middle of the list')
  call relmem(a)
  call check(.not. associated(a),'release at the head of the list')
  call relmem(c)
  call check(.not. associated(c),'release at the tail of the list')

  call getmem(a,1,2,1,2,'a')
  a = 5.0_rk8
  call getmem(a,1,4,1,4,'a')
  call check(all(shape(a) == [4,4]) .and. all(a == 0.0_rk8), &
    'getmem on an associated pointer reallocates it')
  call relmem(a)

  call getmem(z,1,0,1,3,'z')
  call check(associated(z) .and. size(z) == 0,'zero size array')
  call relmem(z)
  call check(.not. associated(z),'zero size array released')
  call getmem(t3,1,3,1,4,0,2,'t3')
  do k = 0, 2
    t3(:,:,k) = real(k,rk8)
  end do
  call assignpnt(t3,s2,2)
  call check(all(lbound(s2) == [1,1]) .and. all(ubound(s2) == [3,4]) &
    .and. all(s2 == 2.0_rk8),'assignpnt of a slice')
  s2 = 7.0_rk8
  call check(all(t3(:,:,2) == 7.0_rk8),'assignpnt shares the storage')
  call relmem(t3)
  call assignpnt(t3,s2,1)
  call check(.not. associated(s2),'assignpnt of a null pointer')

  call getmem(r4,0,3,1,2,0,4,1,2,'r4')
  call remappnt4(r4,m4)
  call check(all(lbound(m4) == [0,1,1,1]) .and. &
    all(ubound(m4) == [3,2,5,2]),'remappnt4 starts the third index at 1')

  call memory_destroy
  call test_end('test_memutil')
end program test_memutil
