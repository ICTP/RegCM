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

module mod_memutil

  use mod_space
  use mod_message
  use mod_realkinds
  use mod_constants
  use mod_date , only : rcm_time_and_date

  private

  public :: memory_init , memory_destroy
  public :: getmem1d , relmem1d
  public :: getmem2d , relmem2d
  public :: getmem3d , relmem3d
  public :: getmem4d , relmem4d
  public :: getmem5d , relmem5d
  public :: assignpnt

  interface assignpnt
    module procedure assignp1d_l
    module procedure assignp1d_s
    module procedure assignp1d_i
    module procedure assignp1d_r
    module procedure assignp1d_d
    module procedure assignp1d_t
    module procedure assignp2d_l
    module procedure assignp2d_s
    module procedure assignp2d_i
    module procedure assignp2d_r
    module procedure assignp2d_d
    module procedure assignp3d_l
    module procedure assignp3d_s
    module procedure assignp3d_i
    module procedure assignp3d_r
    module procedure assignp3d_d
    module procedure assignp4d_l
    module procedure assignp4d_s
    module procedure assignp4d_i
    module procedure assignp4d_r
    module procedure assignp4d_d
    module procedure assignp4d_3d_d
    module procedure assignp5d_l
    module procedure assignp5d_s
    module procedure assignp5d_i
    module procedure assignp5d_r
    module procedure assignp5d_d
  end interface assignpnt

  interface getmem1d
    module procedure getmem1d_l
    module procedure getmem1d_s
    module procedure getmem1d_i
    module procedure getmem1d_r
    module procedure getmem1d_d
    module procedure getmem1d_t
  end interface getmem1d

  interface relmem1d
    module procedure relmem1d_l
    module procedure relmem1d_s
    module procedure relmem1d_i
    module procedure relmem1d_r
    module procedure relmem1d_d
    module procedure relmem1d_t
  end interface relmem1d

  interface getmem2d
    module procedure getmem2d_l
    module procedure getmem2d_s
    module procedure getmem2d_i
    module procedure getmem2d_r
    module procedure getmem2d_d
  end interface getmem2d

  interface relmem2d
    module procedure relmem2d_l
    module procedure relmem2d_s
    module procedure relmem2d_i
    module procedure relmem2d_r
    module procedure relmem2d_d
  end interface relmem2d

  interface getmem3d
    module procedure getmem3d_l
    module procedure getmem3d_s
    module procedure getmem3d_i
    module procedure getmem3d_r
    module procedure getmem3d_d
  end interface getmem3d

  interface relmem3d
    module procedure relmem3d_l
    module procedure relmem3d_s
    module procedure relmem3d_i
    module procedure relmem3d_r
    module procedure relmem3d_d
  end interface relmem3d

  interface getmem4d
    module procedure getmem4d_l
    module procedure getmem4d_s
    module procedure getmem4d_i
    module procedure getmem4d_r
    module procedure getmem4d_d
  end interface getmem4d

  interface relmem4d
    module procedure relmem4d_l
    module procedure relmem4d_s
    module procedure relmem4d_i
    module procedure relmem4d_r
    module procedure relmem4d_d
  end interface relmem4d

  interface getmem5d
    module procedure getmem5d_l
    module procedure getmem5d_s
    module procedure getmem5d_i
    module procedure getmem5d_r
    module procedure getmem5d_d
  end interface getmem5d

  interface relmem5d
    module procedure relmem5d_l
    module procedure relmem5d_s
    module procedure relmem5d_i
    module procedure relmem5d_r
    module procedure relmem5d_d
  end interface relmem5d

  type pool1d_i
    type(pool1d_i) , pointer :: next => null()
    type(iarr1d) :: a
  end type pool1d_i

  type pool1d_s
    type(pool1d_s) , pointer :: next => null()
    type(sarr1d) :: a
  end type pool1d_s

  type pool1d_l
    type(pool1d_l) , pointer :: next => null()
    type(larr1d) :: a
  end type pool1d_l

  type pool1d_r
    type(pool1d_r) , pointer :: next => null()
    type(r4arr1d) :: a
  end type pool1d_r

  type pool1d_d
    type(pool1d_d) , pointer :: next => null()
    type(r8arr1d) :: a
  end type pool1d_d

  type pool1d_t
    type(pool1d_t) , pointer :: next => null()
    type(rtarr1d) :: a
  end type pool1d_t

  type pool2d_i
    type(pool2d_i) , pointer :: next => null()
    type(iarr2d) :: a
  end type pool2d_i

  type pool2d_s
    type(pool2d_s) , pointer :: next => null()
    type(sarr2d) :: a
  end type pool2d_s

  type pool2d_l
    type(pool2d_l) , pointer :: next => null()
    type(larr2d) :: a
  end type pool2d_l

  type pool2d_r
    type(pool2d_r) , pointer :: next => null()
    type(r4arr2d) :: a
  end type pool2d_r

  type pool2d_d
    type(pool2d_d) , pointer :: next => null()
    type(r8arr2d) :: a
  end type pool2d_d

  type pool3d_s
    type(pool3d_s) , pointer :: next => null()
    type(sarr3d) :: a
  end type pool3d_s

  type pool3d_i
    type(pool3d_i) , pointer :: next => null()
    type(iarr3d) :: a
  end type pool3d_i

  type pool3d_l
    type(pool3d_l) , pointer :: next => null()
    type(larr3d) :: a
  end type pool3d_l

  type pool3d_r
    type(pool3d_r) , pointer :: next => null()
    type(r4arr3d) :: a
  end type pool3d_r

  type pool3d_d
    type(pool3d_d) , pointer :: next => null()
    type(r8arr3d) :: a
  end type pool3d_d

  type pool4d_s
    type(pool4d_s) , pointer :: next => null()
    type(sarr4d) :: a
  end type pool4d_s

  type pool4d_i
    type(pool4d_i) , pointer :: next => null()
    type(iarr4d) :: a
  end type pool4d_i

  type pool4d_l
    type(pool4d_l) , pointer :: next => null()
    type(larr4d) :: a
  end type pool4d_l

  type pool4d_r
    type(pool4d_r) , pointer :: next => null()
    type(r4arr4d) :: a
  end type pool4d_r

  type pool4d_d
    type(pool4d_d) , pointer :: next => null()
    type(r8arr4d) :: a
  end type pool4d_d

  type pool5d_s
    type(pool5d_s) , pointer :: next => null()
    type(sarr5d) :: a
  end type pool5d_s

  type pool5d_i
    type(pool5d_i) , pointer :: next => null()
    type(iarr5d) :: a
  end type pool5d_i

  type pool5d_l
    type(pool5d_l) , pointer :: next => null()
    type(larr5d) :: a
  end type pool5d_l

  type pool5d_r
    type(pool5d_r) , pointer :: next => null()
    type(r4arr5d) :: a
  end type pool5d_r

  type pool5d_d
    type(pool5d_d) , pointer :: next => null()
    type(r8arr5d) :: a
  end type pool5d_d

  type (pool1d_i) , pointer :: r1di , l1di , c1di , n1di , p1di
  type (pool1d_s) , pointer :: r1ds , l1ds , c1ds , n1ds , p1ds
  type (pool1d_l) , pointer :: r1dl , l1dl , c1dl , n1dl , p1dl
  type (pool1d_r) , pointer :: r1dr , l1dr , c1dr , n1dr , p1dr
  type (pool1d_d) , pointer :: r1dd , l1dd , c1dd , n1dd , p1dd
  type (pool1d_t) , pointer :: r1dt , l1dt , c1dt , n1dt , p1dt

  type (pool2d_i) , pointer :: r2di , l2di , c2di , n2di , p2di
  type (pool2d_s) , pointer :: r2ds , l2ds , c2ds , n2ds , p2ds
  type (pool2d_l) , pointer :: r2dl , l2dl , c2dl , n2dl , p2dl
  type (pool2d_r) , pointer :: r2dr , l2dr , c2dr , n2dr , p2dr
  type (pool2d_d) , pointer :: r2dd , l2dd , c2dd , n2dd , p2dd

  type (pool3d_i) , pointer :: r3di , l3di , c3di , n3di , p3di
  type (pool3d_s) , pointer :: r3ds , l3ds , c3ds , n3ds , p3ds
  type (pool3d_l) , pointer :: r3dl , l3dl , c3dl , n3dl , p3dl
  type (pool3d_r) , pointer :: r3dr , l3dr , c3dr , n3dr , p3dr
  type (pool3d_d) , pointer :: r3dd , l3dd , c3dd , n3dd , p3dd

  type (pool4d_i) , pointer :: r4di , l4di , c4di , n4di , p4di
  type (pool4d_s) , pointer :: r4ds , l4ds , c4ds , n4ds , p4ds
  type (pool4d_l) , pointer :: r4dl , l4dl , c4dl , n4dl , p4dl
  type (pool4d_r) , pointer :: r4dr , l4dr , c4dr , n4dr , p4dr
  type (pool4d_d) , pointer :: r4dd , l4dd , c4dd , n4dd , p4dd

  type (pool5d_i) , pointer :: r5di , l5di , c5di , n5di , p5di
  type (pool5d_s) , pointer :: r5ds , l5ds , c5ds , n5ds , p5ds
  type (pool5d_l) , pointer :: r5dl , l5dl , c5dl , n5dl , p5dl
  type (pool5d_r) , pointer :: r5dr , l5dr , c5dr , n5dr , p5dr
  type (pool5d_d) , pointer :: r5dd , l5dd , c5dd , n5dd , p5dd

  integer :: ista

  contains

  subroutine memory_init
    implicit none
    allocate(r1di, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1d1')
    allocate(r1ds, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1ds')
    allocate(r1dl, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1dl')
    allocate(r1dr, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1dr')
    allocate(r1dd, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1dd')
    allocate(r1dt, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r1dt')
    l1di => r1di
    l1ds => r1ds
    l1dl => r1dl
    l1dr => r1dr
    l1dd => r1dd
    l1dt => r1dt
    allocate(r2di, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r2d1')
    allocate(r2ds, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r2ds')
    allocate(r2dl, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r2dl')
    allocate(r2dr, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r2dr')
    allocate(r2dd, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r2dd')
    l2di => r2di
    l2ds => r2ds
    l2dl => r2dl
    l2dr => r2dr
    l2dd => r2dd
    allocate(r3di, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r3d1')
    allocate(r3ds, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r3ds')
    allocate(r3dl, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r3dl')
    allocate(r3dr, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r3dr')
    allocate(r3dd, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r3dd')
    l3di => r3di
    l3ds => r3ds
    l3dl => r3dl
    l3dr => r3dr
    l3dd => r3dd
    allocate(r4di, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r4d1')
    allocate(r4ds, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r4ds')
    allocate(r4dl, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r4dl')
    allocate(r4dr, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r4dr')
    allocate(r4dd, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r4dd')
    l4di => r4di
    l4ds => r4ds
    l4dl => r4dl
    l4dr => r4dr
    l4dd => r4dd
    allocate(r5di, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r5d1')
    allocate(r5ds, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r5ds')
    allocate(r5dl, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r5dl')
    allocate(r5dr, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r5dr')
    allocate(r5dd, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'r5dd')
    l5di => r5di
    l5ds => r5ds
    l5dl => r5dl
    l5dr => r5dr
    l5dd => r5dd
  end subroutine memory_init

  subroutine getmem1d_l(a,l,h,vn)
    implicit none
    logical , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1dl => l1dl
    call getspc1d(c1dl%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1dl%a%space
    a(:) = .false.
    allocate(c1dl%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1dl%next')
    l1dl => c1dl%next
  end subroutine getmem1d_l

  subroutine relmem1d_l(a)
    implicit none
    logical , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1dl => null()
    c1dl => r1dl
    do while ( associated(c1dl) )
      n1dl => c1dl%next
      if ( associated(a,c1dl%a%space) ) then
        deallocate(c1dl%a%space)
        a => null()
        if ( associated(p1dl) ) then
          if ( associated(n1dl) ) then
            p1dl%next => n1dl
          else
            l1dl => p1dl
            p1dl%next => null()
          end if
        else
          r1dl => n1dl
        end if
        deallocate(c1dl)
        exit
      end if
      p1dl => c1dl
      c1dl => c1dl%next
    end do
  end subroutine relmem1d_l

  subroutine getmem1d_t(a,l,h,vn)
    implicit none
    type(rcm_time_and_date) , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1dt => l1dt
    call getspc1d(c1dt%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1dt%a%space
    allocate(c1dt%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1dt%next')
    l1dt => c1dt%next
  end subroutine getmem1d_t

  subroutine relmem1d_t(a)
    implicit none
    type(rcm_time_and_date) , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1dt => null()
    c1dt => r1dt
    do while ( associated(c1dt) )
      n1dt => c1dt%next
      if ( associated(a,c1dt%a%space) ) then
        deallocate(c1dt%a%space)
        a => null()
        if ( associated(p1dt) ) then
          if ( associated(n1dt) ) then
            p1dt%next => n1dt
          else
            l1dt => p1dt
            p1dt%next => null()
          end if
        else
          r1dt => n1dt
        end if
        deallocate(c1dt)
        exit
      end if
      p1dt => c1dt
      c1dt => c1dt%next
    end do
  end subroutine relmem1d_t

  subroutine getmem1d_s(a,l,h,vn)
    implicit none
    integer(2) , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1ds => l1ds
    call getspc1d(c1ds%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1ds%a%space
    a(:) = -1_2
    allocate(c1ds%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1ds%next')
    l1ds => c1ds%next
  end subroutine getmem1d_s

  subroutine relmem1d_s(a)
    implicit none
    integer(2) , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1ds => null()
    c1ds => r1ds
    do while ( associated(c1ds) )
      n1ds => c1ds%next
      if ( associated(a,c1ds%a%space) ) then
        deallocate(c1ds%a%space)
        a => null()
        if ( associated(p1ds) ) then
          if ( associated(n1ds) ) then
            p1ds%next => n1ds
          else
            l1ds => p1ds
            p1ds%next => null()
          end if
        else
          r1ds => n1ds
        end if
        deallocate(c1ds)
        exit
      end if
      p1ds => c1ds
      c1ds => c1ds%next
    end do
  end subroutine relmem1d_s

  subroutine getmem1d_i(a,l,h,vn)
    implicit none
    integer , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1di => l1di
    call getspc1d(c1di%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1di%a%space
    a(:) = -1
    allocate(c1di%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1di%next')
    l1di => c1di%next
  end subroutine getmem1d_i

  subroutine relmem1d_i(a)
    implicit none
    integer , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1di => null()
    c1di => r1di
    do while ( associated(c1di) )
      n1di => c1di%next
      if ( associated(a,c1di%a%space) ) then
        deallocate(c1di%a%space)
        a => null()
        if ( associated(p1di) ) then
          if ( associated(n1di) ) then
            p1di%next => n1di
          else
            l1di => p1di
            p1di%next => null()
          end if
        else
          r1di => n1di
        end if
        deallocate(c1di)
        exit
      end if
      p1di => c1di
      c1di => c1di%next
    end do
  end subroutine relmem1d_i

  subroutine getmem1d_r(a,l,h,vn)
    implicit none
    real(sp) , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1dr => l1dr
    call getspc1d(c1dr%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1dr%a%space
    a(:) = 0.0
    allocate(c1dr%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1dr%next')
    l1dr => c1dr%next
  end subroutine getmem1d_r

  subroutine relmem1d_r(a)
    implicit none
    real(sp) , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1dr => null()
    c1dr => r1dr
    do while ( associated(c1dr) )
      n1dr => c1dr%next
      if ( associated(a,c1dr%a%space) ) then
        deallocate(c1dr%a%space)
        a => null()
        if ( associated(p1dr) ) then
          if ( associated(n1dr) ) then
            p1dr%next => n1dr
          else
            l1dr => p1dr
            p1dr%next => null()
          end if
        else
          r1dr => n1dr
        end if
        deallocate(c1dr)
        exit
      end if
      p1dr => c1dr
      c1dr => c1dr%next
    end do
  end subroutine relmem1d_r

  subroutine getmem1d_d(a,l,h,vn)
    implicit none
    real(dp) , pointer , dimension(:) , intent(out) :: a
    integer , intent(in) :: l , h
    character (len=*) , intent(in) :: vn
    type (bounds) :: b
    if ( associated(a) ) call relmem1d(a)
    b = bounds(l,h)
    c1dd => l1dd
    call getspc1d(c1dd%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c1dd%a%space
    a(:) = d_zero
    allocate(c1dd%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c1dd%next')
    l1dd => c1dd%next
  end subroutine getmem1d_d

  subroutine relmem1d_d(a)
    implicit none
    real(dp) , pointer , dimension(:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p1dd => null()
    c1dd => r1dd
    do while ( associated(c1dd) )
      n1dd => c1dd%next
      if ( associated(a,c1dd%a%space) ) then
        deallocate(c1dd%a%space)
        a => null()
        if ( associated(p1dd) ) then
          if ( associated(n1dd) ) then
            p1dd%next => n1dd
          else
            l1dd => p1dd
            p1dd%next => null()
          end if
        else
          r1dd => n1dd
        end if
        deallocate(c1dd)
        exit
      end if
      p1dd => c1dd
      c1dd => c1dd%next
    end do
  end subroutine relmem1d_d

  subroutine finalize_pool1d_i(n)
    implicit none
    type(pool1d_i) , intent(inout) , pointer :: n
    type(pool1d_i) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_i

  subroutine finalize_pool1d_s(n)
    implicit none
    type(pool1d_s) , intent(inout) , pointer :: n
    type(pool1d_s) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_s

  subroutine finalize_pool1d_l(n)
    implicit none
    type(pool1d_l) , intent(inout) , pointer :: n
    type(pool1d_l) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_l

  subroutine finalize_pool1d_r(n)
    implicit none
    type(pool1d_r) , intent(inout) , pointer :: n
    type(pool1d_r) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_r

  subroutine finalize_pool1d_d(n)
    implicit none
    type(pool1d_d) , intent(inout) , pointer :: n
    type(pool1d_d) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_d

  subroutine finalize_pool1d_t(n)
    implicit none
    type(pool1d_t) , intent(inout) , pointer :: n
    type(pool1d_t) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool1d_t

  subroutine getmem2d_l(a,l1,h1,l2,h2,vn)
    implicit none
    logical , pointer , dimension(:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(2) :: b
    if ( associated(a) ) call relmem2d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    c2dl => l2dl
    call getspc2d(c2dl%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c2dl%a%space
    a(:,:) = .false.
    allocate(c2dl%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c2dl%next')
    l2dl => c2dl%next
  end subroutine getmem2d_l

  subroutine relmem2d_l(a)
    implicit none
    logical , pointer , dimension(:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p2dl => null()
    c2dl => r2dl
    do while ( associated(c2dl) )
      n2dl => c2dl%next
      if ( associated(a,c2dl%a%space) ) then
        deallocate(c2dl%a%space)
        a => null()
        if ( associated(p2dl) ) then
          if ( associated(n2dl) ) then
            p2dl%next => n2dl
          else
            l2dl => p2dl
            p2dl%next => null()
          end if
        else
          r2dl => n2dl
        end if
        deallocate(c2dl)
        exit
      end if
      p2dl => c2dl
      c2dl => c2dl%next
    end do
  end subroutine relmem2d_l

  subroutine getmem2d_s(a,l1,h1,l2,h2,vn)
    implicit none
    integer(2) , pointer , dimension(:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(2) :: b
    if ( associated(a) ) call relmem2d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    c2ds => l2ds
    call getspc2d(c2ds%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c2ds%a%space
    a(:,:) = -1_2
    allocate(c2ds%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c2ds%next')
    l2ds => c2ds%next
  end subroutine getmem2d_s

  subroutine relmem2d_s(a)
    implicit none
    integer(2) , pointer , dimension(:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p2ds => null()
    c2ds => r2ds
    do while ( associated(c2ds) )
      n2ds => c2ds%next
      if ( associated(a,c2ds%a%space) ) then
        deallocate(c2ds%a%space)
        a => null()
        if ( associated(p2ds) ) then
          if ( associated(n2ds) ) then
            p2ds%next => n2ds
          else
            l2ds => p2ds
            p2ds%next => null()
          end if
        else
          r2ds => n2ds
        end if
        deallocate(c2ds)
        exit
      end if
      p2ds => c2ds
      c2ds => c2ds%next
    end do
  end subroutine relmem2d_s

  subroutine getmem2d_i(a,l1,h1,l2,h2,vn)
    implicit none
    integer , pointer , dimension(:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(2) :: b
    if ( associated(a) ) call relmem2d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    c2di => l2di
    call getspc2d(c2di%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c2di%a%space
    a(:,:) = -1
    allocate(c2di%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c2di%next')
    l2di => c2di%next
  end subroutine getmem2d_i

  subroutine relmem2d_i(a)
    implicit none
    integer , pointer , dimension(:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p2di => null()
    c2di => r2di
    do while ( associated(c2di) )
      n2di => c2di%next
      if ( associated(a,c2di%a%space) ) then
        deallocate(c2di%a%space)
        a => null()
        if ( associated(p2di) ) then
          if ( associated(n2di) ) then
            p2di%next => n2di
          else
            l2di => p2di
            p2di%next => null()
          end if
        else
          r2di => n2di
        end if
        deallocate(c2di)
        exit
      end if
      p2di => c2di
      c2di => c2di%next
    end do
  end subroutine relmem2d_i

  subroutine getmem2d_r(a,l1,h1,l2,h2,vn)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(2) :: b
    if ( associated(a) ) call relmem2d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    c2dr => l2dr
    call getspc2d(c2dr%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c2dr%a%space
    a(:,:) = 0.0
    allocate(c2dr%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c2dr%next')
    l2dr => c2dr%next
  end subroutine getmem2d_r

  subroutine relmem2d_r(a)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p2dr => null()
    c2dr => r2dr
    do while ( associated(c2dr) )
      n2dr => c2dr%next
      if ( associated(a,c2dr%a%space) ) then
        deallocate(c2dr%a%space)
        a => null()
        if ( associated(p2dr) ) then
          if ( associated(n2dr) ) then
            p2dr%next => n2dr
          else
            l2dr => p2dr
            p2dr%next => null()
          end if
        else
          r2dr => n2dr
        end if
        deallocate(c2dr)
        exit
      end if
      p2dr => c2dr
      c2dr => c2dr%next
    end do
  end subroutine relmem2d_r

  subroutine getmem2d_d(a,l1,h1,l2,h2,vn)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(2) :: b
    if ( associated(a) ) call relmem2d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    c2dd => l2dd
    call getspc2d(c2dd%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c2dd%a%space
    a(:,:) = d_zero
    allocate(c2dd%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c2dd%next')
    l2dd => c2dd%next
  end subroutine getmem2d_d

  subroutine relmem2d_d(a)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p2dd => null()
    c2dd => r2dd
    do while ( associated(c2dd) )
      n2dd => c2dd%next
      if ( associated(a,c2dd%a%space) ) then
        deallocate(c2dd%a%space)
        a => null()
        if ( associated(p2dd) ) then
          if ( associated(n2dd) ) then
            p2dd%next => n2dd
          else
            l2dd => p2dd
            p2dd%next => null()
          end if
        else
          r2dd => n2dd
        end if
        deallocate(c2dd)
        exit
      end if
      p2dd => c2dd
      c2dd => c2dd%next
    end do
  end subroutine relmem2d_d

  subroutine finalize_pool2d_i(n)
    implicit none
    type(pool2d_i) , intent(inout) , pointer :: n
    type(pool2d_i) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool2d_i

  subroutine finalize_pool2d_s(n)
    implicit none
    type(pool2d_s) , intent(inout) , pointer :: n
    type(pool2d_s) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool2d_s

  subroutine finalize_pool2d_l(n)
    implicit none
    type(pool2d_l) , intent(inout) , pointer :: n
    type(pool2d_l) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool2d_l

  subroutine finalize_pool2d_r(n)
    implicit none
    type(pool2d_r) , intent(inout) , pointer :: n
    type(pool2d_r) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool2d_r

  subroutine finalize_pool2d_d(n)
    implicit none
    type(pool2d_d) , intent(inout) , pointer :: n
    type(pool2d_d) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool2d_d

  subroutine getmem3d_l(a,l1,h1,l2,h2,l3,h3,vn)
    implicit none
    logical , pointer , dimension(:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(3) :: b
    if ( associated(a) ) call relmem3d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    c3dl => l3dl
    call getspc3d(c3dl%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c3dl%a%space
    a(:,:,:) = .false.
    allocate(c3dl%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c3dl%next')
    l3dl => c3dl%next
  end subroutine getmem3d_l

  subroutine relmem3d_l(a)
    implicit none
    logical , pointer , dimension(:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p3dl => null()
    c3dl => r3dl
    do while ( associated(c3dl) )
      n3dl => c3dl%next
      if ( associated(a,c3dl%a%space) ) then
        deallocate(c3dl%a%space)
        a => null()
        if ( associated(p3dl) ) then
          if ( associated(n3dl) ) then
            p3dl%next => n3dl
          else
            l3dl => p3dl
            p3dl%next => null()
          end if
        else
          r3dl => n3dl
        end if
        deallocate(c3dl)
        exit
      end if
      p3dl => c3dl
      c3dl => c3dl%next
    end do
  end subroutine relmem3d_l

  subroutine getmem3d_s(a,l1,h1,l2,h2,l3,h3,vn)
    implicit none
    integer(2) , pointer , dimension(:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(3) :: b
    if ( associated(a) ) call relmem3d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    c3ds => l3ds
    call getspc3d(c3ds%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c3ds%a%space
    a(:,:,:) = -1_2
    allocate(c3ds%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c3ds%next')
    l3ds => c3ds%next
  end subroutine getmem3d_s

  subroutine relmem3d_s(a)
    implicit none
    integer(2) , pointer , dimension(:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p3ds => null()
    c3ds => r3ds
    do while ( associated(c3ds) )
      n3ds => c3ds%next
      if ( associated(a,c3ds%a%space) ) then
        deallocate(c3ds%a%space)
        a => null()
        if ( associated(p3ds) ) then
          if ( associated(n3ds) ) then
            p3ds%next => n3ds
          else
            l3ds => p3ds
            p3ds%next => null()
          end if
        else
          r3ds => n3ds
        end if
        deallocate(c3ds)
        exit
      end if
      p3ds => c3ds
      c3ds => c3ds%next
    end do
  end subroutine relmem3d_s

  subroutine getmem3d_i(a,l1,h1,l2,h2,l3,h3,vn)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(3) :: b
    if ( associated(a) ) call relmem3d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    c3di => l3di
    call getspc3d(c3di%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c3di%a%space
    a(:,:,:) = -1
    allocate(c3di%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c3di%next')
    l3di => c3di%next
  end subroutine getmem3d_i

  subroutine relmem3d_i(a)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p3di => null()
    c3di => r3di
    do while ( associated(c3di) )
      n3di => c3di%next
      if ( associated(a,c3di%a%space) ) then
        deallocate(c3di%a%space)
        a => null()
        if ( associated(p3di) ) then
          if ( associated(n3di) ) then
            p3di%next => n3di
          else
            l3di => p3di
            p3di%next => null()
          end if
        else
          r3di => n3di
        end if
        deallocate(c3di)
        exit
      end if
      p3di => c3di
      c3di => c3di%next
    end do
  end subroutine relmem3d_i

  subroutine getmem3d_r(a,l1,h1,l2,h2,l3,h3,vn)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(3) :: b
    if ( associated(a) ) call relmem3d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    c3dr => l3dr
    call getspc3d(c3dr%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c3dr%a%space
    a(:,:,:) = 0.0
    allocate(c3dr%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c3dr%next')
    l3dr => c3dr%next
  end subroutine getmem3d_r

  subroutine relmem3d_r(a)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p3dr => null()
    c3dr => r3dr
    do while ( associated(c3dr) )
      n3dr => c3dr%next
      if ( associated(a,c3dr%a%space) ) then
        deallocate(c3dr%a%space)
        a => null()
        if ( associated(p3dr) ) then
          if ( associated(n3dr) ) then
            p3dr%next => n3dr
          else
            l3dr => p3dr
            p3dr%next => null()
          end if
        else
          r3dr => n3dr
        end if
        deallocate(c3dr)
        exit
      end if
      p3dr => c3dr
      c3dr => c3dr%next
    end do
  end subroutine relmem3d_r

  subroutine getmem3d_d(a,l1,h1,l2,h2,l3,h3,vn)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(3) :: b
    if ( associated(a) ) call relmem3d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    c3dd => l3dd
    call getspc3d(c3dd%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c3dd%a%space
    a(:,:,:) = d_zero
    allocate(c3dd%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c3dd%next')
    l3dd => c3dd%next
  end subroutine getmem3d_d

  subroutine relmem3d_d(a)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p3dd => null()
    c3dd => r3dd
    do while ( associated(c3dd) )
      n3dd => c3dd%next
      if ( associated(a,c3dd%a%space) ) then
        deallocate(c3dd%a%space)
        a => null()
        if ( associated(p3dd) ) then
          if ( associated(n3dd) ) then
            p3dd%next => n3dd
          else
            l3dd => p3dd
            p3dd%next => null()
          end if
        else
          r3dd => n3dd
        end if
        deallocate(c3dd)
        exit
      end if
      p3dd => c3dd
      c3dd => c3dd%next
    end do
  end subroutine relmem3d_d

  subroutine finalize_pool3d_s(n)
    implicit none
    type(pool3d_s) , intent(inout) , pointer :: n
    type(pool3d_s) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool3d_s

  subroutine finalize_pool3d_i(n)
    implicit none
    type(pool3d_i) , intent(inout) , pointer :: n
    type(pool3d_i) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool3d_i

  subroutine finalize_pool3d_l(n)
    implicit none
    type(pool3d_l) , intent(inout) , pointer :: n
    type(pool3d_l) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool3d_l

  subroutine finalize_pool3d_r(n)
    implicit none
    type(pool3d_r) , intent(inout) , pointer :: n
    type(pool3d_r) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool3d_r

  subroutine finalize_pool3d_d(n)
    implicit none
    type(pool3d_d) , intent(inout) , pointer :: n
    type(pool3d_d) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool3d_d

  subroutine getmem4d_l(a,l1,h1,l2,h2,l3,h3,l4,h4,vn)
    implicit none
    logical , pointer , dimension(:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(4) :: b
    if ( associated(a) ) call relmem4d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    c4dl => l4dl
    call getspc4d(c4dl%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c4dl%a%space
    a(:,:,:,:) = .false.
    allocate(c4dl%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c4dl%next')
    l4dl => c4dl%next
  end subroutine getmem4d_l

  subroutine relmem4d_l(a)
    implicit none
    logical , pointer , dimension(:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p4dl => null()
    c4dl => r4dl
    do while ( associated(c4dl) )
      n4dl => c4dl%next
      if ( associated(a,c4dl%a%space) ) then
        deallocate(c4dl%a%space)
        a => null()
        if ( associated(p4dl) ) then
          if ( associated(n4dl) ) then
            p4dl%next => n4dl
          else
            l4dl => p4dl
            p4dl%next => null()
          end if
        else
          r4dl => n4dl
        end if
        deallocate(c4dl)
        exit
      end if
      p4dl => c4dl
      c4dl => c4dl%next
    end do
  end subroutine relmem4d_l

  subroutine getmem4d_s(a,l1,h1,l2,h2,l3,h3,l4,h4,vn)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(4) :: b
    if ( associated(a) ) call relmem4d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    c4ds => l4ds
    call getspc4d(c4ds%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c4ds%a%space
    a(:,:,:,:) = -1_2
    allocate(c4ds%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c4ds%next')
    l4ds => c4ds%next
  end subroutine getmem4d_s

  subroutine relmem4d_s(a)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p4ds => null()
    c4ds => r4ds
    do while ( associated(c4ds) )
      n4ds => c4ds%next
      if ( associated(a,c4ds%a%space) ) then
        deallocate(c4ds%a%space)
        a => null()
        if ( associated(p4ds) ) then
          if ( associated(n4ds) ) then
            p4ds%next => n4ds
          else
            l4ds => p4ds
            p4ds%next => null()
          end if
        else
          r4ds => n4ds
        end if
        deallocate(c4ds)
        exit
      end if
      p4ds => c4ds
      c4ds => c4ds%next
    end do
  end subroutine relmem4d_s

  subroutine getmem4d_i(a,l1,h1,l2,h2,l3,h3,l4,h4,vn)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(4) :: b
    if ( associated(a) ) call relmem4d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    c4di => l4di
    call getspc4d(c4di%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c4di%a%space
    a(:,:,:,:) = -1
    allocate(c4di%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c4di%next')
    l4di => c4di%next
  end subroutine getmem4d_i

  subroutine relmem4d_i(a)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p4di => null()
    c4di => r4di
    do while ( associated(c4di) )
      n4di => c4di%next
      if ( associated(a,c4di%a%space) ) then
        deallocate(c4di%a%space)
        a => null()
        if ( associated(p4di) ) then
          if ( associated(n4di) ) then
            p4di%next => n4di
          else
            l4di => p4di
            p4di%next => null()
          end if
        else
          r4di => n4di
        end if
        deallocate(c4di)
        exit
      end if
      p4di => c4di
      c4di => c4di%next
    end do
  end subroutine relmem4d_i

  subroutine getmem4d_r(a,l1,h1,l2,h2,l3,h3,l4,h4,vn)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(4) :: b
    if ( associated(a) ) call relmem4d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    c4dr => l4dr
    call getspc4d(c4dr%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c4dr%a%space
    a(:,:,:,:) = 0.0
    allocate(c4dr%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c4dr%next')
    l4dr => c4dr%next
  end subroutine getmem4d_r

  subroutine relmem4d_r(a)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p4dr => null()
    c4dr => r4dr
    do while ( associated(c4dr) )
      n4dr => c4dr%next
      if ( associated(a,c4dr%a%space) ) then
        deallocate(c4dr%a%space)
        a => null()
        if ( associated(p4dr) ) then
          if ( associated(n4dr) ) then
            p4dr%next => n4dr
          else
            l4dr => p4dr
            p4dr%next => null()
          end if
        else
          r4dr => n4dr
        end if
        deallocate(c4dr)
        exit
      end if
      p4dr => c4dr
      c4dr => c4dr%next
    end do
  end subroutine relmem4d_r

  subroutine getmem4d_d(a,l1,h1,l2,h2,l3,h3,l4,h4,vn)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(4) :: b
    if ( associated(a) ) call relmem4d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    c4dd => l4dd
    call getspc4d(c4dd%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c4dd%a%space
    a(:,:,:,:) = d_zero
    allocate(c4dd%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c4dd%next')
    l4dd => c4dd%next
  end subroutine getmem4d_d

  subroutine relmem4d_d(a)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p4dd => null()
    c4dd => r4dd
    do while ( associated(c4dd) )
      n4dd => c4dd%next
      if ( associated(a,c4dd%a%space) ) then
        deallocate(c4dd%a%space)
        a => null()
        if ( associated(p4dd) ) then
          if ( associated(n4dd) ) then
            p4dd%next => n4dd
          else
            l4dd => p4dd
            p4dd%next => null()
          end if
        else
          r4dd => n4dd
        end if
        deallocate(c4dd)
        exit
      end if
      p4dd => c4dd
      c4dd => c4dd%next
    end do
  end subroutine relmem4d_d

  subroutine finalize_pool4d_s(n)
    implicit none
    type(pool4d_s) , intent(inout) , pointer :: n
    type(pool4d_s) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool4d_s

  subroutine finalize_pool4d_i(n)
    implicit none
    type(pool4d_i) , intent(inout) , pointer :: n
    type(pool4d_i) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool4d_i

  subroutine finalize_pool4d_l(n)
    implicit none
    type(pool4d_l) , intent(inout) , pointer :: n
    type(pool4d_l) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool4d_l

  subroutine finalize_pool4d_r(n)
    implicit none
    type(pool4d_r) , intent(inout) , pointer :: n
    type(pool4d_r) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool4d_r

  subroutine finalize_pool4d_d(n)
    implicit none
    type(pool4d_d) , intent(inout) , pointer :: n
    type(pool4d_d) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool4d_d

  subroutine getmem5d_l(a,l1,h1,l2,h2,l3,h3,l4,h4,l5,h5,vn)
    implicit none
    logical , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4 , l5 , h5
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(5) :: b
    if ( associated(a) ) call relmem5d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    b(5) = bounds(l5,h5)
    c5dl => l5dl
    call getspc5d(c5dl%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c5dl%a%space
    a(:,:,:,:,:) = .false.
    allocate(c5dl%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c5dl%next')
    l5dl => c5dl%next
  end subroutine getmem5d_l

  subroutine relmem5d_l(a)
    implicit none
    logical , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p5dl => null()
    c5dl => r5dl
    do while ( associated(c5dl) )
      n5dl => c5dl%next
      if ( associated(a,c5dl%a%space) ) then
        deallocate(c5dl%a%space)
        a => null()
        if ( associated(p5dl) ) then
          if ( associated(n5dl) ) then
            p5dl%next => n5dl
          else
            l5dl => p5dl
            p5dl%next => null()
          end if
        else
          r5dl => n5dl
        end if
        deallocate(c5dl)
        exit
      end if
      p5dl => c5dl
      c5dl => c5dl%next
    end do
  end subroutine relmem5d_l

  subroutine getmem5d_s(a,l1,h1,l2,h2,l3,h3,l4,h4,l5,h5,vn)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4 , l5 , h5
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(5) :: b
    if ( associated(a) ) call relmem5d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    b(5) = bounds(l5,h5)
    c5ds => l5ds
    call getspc5d(c5ds%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c5ds%a%space
    a(:,:,:,:,:) = -1_2
    allocate(c5ds%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c5ds%next')
    l5ds => c5ds%next
  end subroutine getmem5d_s

  subroutine relmem5d_s(a)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p5ds => null()
    c5ds => r5ds
    do while ( associated(c5ds) )
      n5ds => c5ds%next
      if ( associated(a,c5ds%a%space) ) then
        deallocate(c5ds%a%space)
        a => null()
        if ( associated(p5ds) ) then
          if ( associated(n5ds) ) then
            p5ds%next => n5ds
          else
            l5ds => p5ds
            p5ds%next => null()
          end if
        else
          r5ds => n5ds
        end if
        deallocate(c5ds)
        exit
      end if
      p5ds => c5ds
      c5ds => c5ds%next
    end do
  end subroutine relmem5d_s

  subroutine getmem5d_i(a,l1,h1,l2,h2,l3,h3,l4,h4,l5,h5,vn)
    implicit none
    integer , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4 , l5 , h5
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(5) :: b
    if ( associated(a) ) call relmem5d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    b(5) = bounds(l5,h5)
    c5di => l5di
    call getspc5d(c5di%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c5di%a%space
    a(:,:,:,:,:) = -1
    allocate(c5di%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c5di%next')
    l5di => c5di%next
  end subroutine getmem5d_i

  subroutine relmem5d_i(a)
    implicit none
    integer , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p5di => null()
    c5di => r5di
    do while ( associated(c5di) )
      n5di => c5di%next
      if ( associated(a,c5di%a%space) ) then
        deallocate(c5di%a%space)
        a => null()
        if ( associated(p5di) ) then
          if ( associated(n5di) ) then
            p5di%next => n5di
          else
            l5di => p5di
            p5di%next => null()
          end if
        else
          r5di => n5di
        end if
        deallocate(c5di)
        exit
      end if
      p5di => c5di
      c5di => c5di%next
    end do
  end subroutine relmem5d_i

  subroutine getmem5d_r(a,l1,h1,l2,h2,l3,h3,l4,h4,l5,h5,vn)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4 , l5 , h5
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(5) :: b
    if ( associated(a) ) call relmem5d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    b(5) = bounds(l5,h5)
    c5dr => l5dr
    call getspc5d(c5dr%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c5dr%a%space
    a(:,:,:,:,:) = 0.0
    allocate(c5dr%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c5dr%next')
    l5dr => c5dr%next
  end subroutine getmem5d_r

  subroutine relmem5d_r(a)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p5dr => null()
    c5dr => r5dr
    do while ( associated(c5dr) )
      n5dr => c5dr%next
      if ( associated(a,c5dr%a%space) ) then
        deallocate(c5dr%a%space)
        a => null()
        if ( associated(p5dr) ) then
          if ( associated(n5dr) ) then
            p5dr%next => n5dr
          else
            l5dr => p5dr
            p5dr%next => null()
          end if
        else
          r5dr => n5dr
        end if
        deallocate(c5dr)
        exit
      end if
      p5dr => c5dr
      c5dr => c5dr%next
    end do
  end subroutine relmem5d_r

  subroutine getmem5d_d(a,l1,h1,l2,h2,l3,h3,l4,h4,l5,h5,vn)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    integer , intent(in) :: l1 , h1 , l2 , h2 , l3 , h3 , l4 , h4 , l5 , h5
    character (len=*) , intent(in) :: vn
    type (bounds) , dimension(5) :: b
    if ( associated(a) ) call relmem5d(a)
    b(1) = bounds(l1,h1)
    b(2) = bounds(l2,h2)
    b(3) = bounds(l3,h3)
    b(4) = bounds(l4,h4)
    b(5) = bounds(l5,h5)
    c5dd => l5dd
    call getspc5d(c5dd%a,b,ista)
    call checkalloc(ista,__FILE__,__LINE__,vn)
    a => c5dd%a%space
    a(:,:,:,:,:) = d_zero
    allocate(c5dd%next, stat=ista)
    call checkalloc(ista,__FILE__,__LINE__,'c5dd%next')
    l5dd => c5dd%next
  end subroutine getmem5d_d

  subroutine relmem5d_d(a)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:,:) , intent(out) :: a
    if ( .not. associated(a) ) return
    p5dd => null()
    c5dd => r5dd
    do while ( associated(c5dd) )
      n5dd => c5dd%next
      if ( associated(a,c5dd%a%space) ) then
        deallocate(c5dd%a%space)
        a => null()
        if ( associated(p5dd) ) then
          if ( associated(n5dd) ) then
            p5dd%next => n5dd
          else
            l5dd => p5dd
            p5dd%next => null()
          end if
        else
          r5dd => n5dd
        end if
        deallocate(c5dd)
        exit
      end if
      p5dd => c5dd
      c5dd => c5dd%next
    end do
  end subroutine relmem5d_d

  subroutine finalize_pool5d_s(n)
    implicit none
    type(pool5d_s) , intent(inout) , pointer :: n
    type(pool5d_s) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool5d_s

  subroutine finalize_pool5d_i(n)
    implicit none
    type(pool5d_i) , intent(inout) , pointer :: n
    type(pool5d_i) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool5d_i

  subroutine finalize_pool5d_l(n)
    implicit none
    type(pool5d_l) , intent(inout) , pointer :: n
    type(pool5d_l) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool5d_l

  subroutine finalize_pool5d_r(n)
    implicit none
    type(pool5d_r) , intent(inout) , pointer :: n
    type(pool5d_r) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool5d_r

  subroutine finalize_pool5d_d(n)
    implicit none
    type(pool5d_d) , intent(inout) , pointer :: n
    type(pool5d_d) , pointer :: p
    do while ( associated(n) )
      p => n%next
      if ( associated(n%a%space) ) then
        deallocate(n%a%space)
      end if
      deallocate(n)
      nullify(n)
      n => p
    end do
  end subroutine finalize_pool5d_d

  subroutine memory_destroy
    implicit none
    call finalize_pool1d_i(r1di)
    call finalize_pool1d_s(r1ds)
    call finalize_pool1d_l(r1dl)
    call finalize_pool1d_r(r1dr)
    call finalize_pool1d_d(r1dd)
    call finalize_pool1d_t(r1dt)
    call finalize_pool2d_i(r2di)
    call finalize_pool2d_s(r2ds)
    call finalize_pool2d_l(r2dl)
    call finalize_pool2d_r(r2dr)
    call finalize_pool2d_d(r2dd)
    call finalize_pool3d_i(r3di)
    call finalize_pool3d_s(r3ds)
    call finalize_pool3d_l(r3dl)
    call finalize_pool3d_r(r3dr)
    call finalize_pool3d_d(r3dd)
    call finalize_pool4d_i(r4di)
    call finalize_pool4d_s(r4ds)
    call finalize_pool4d_l(r4dl)
    call finalize_pool4d_r(r4dr)
    call finalize_pool4d_d(r4dd)
    call finalize_pool5d_i(r5di)
    call finalize_pool5d_s(r5ds)
    call finalize_pool5d_l(r5dl)
    call finalize_pool5d_r(r5dr)
    call finalize_pool5d_d(r5dd)
  end subroutine memory_destroy

  subroutine assignp1d_l(a,b)
    implicit none
    logical , pointer , dimension(:) , intent(in) :: a
    logical , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_l

  subroutine assignp1d_s(a,b)
    implicit none
    integer(2) , pointer , dimension(:) , intent(in) :: a
    integer(2) , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_s

  subroutine assignp1d_i(a,b)
    implicit none
    integer , pointer , dimension(:) , intent(in) :: a
    integer , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_i

  subroutine assignp1d_r(a,b)
    implicit none
    real(sp) , pointer , dimension(:) , intent(in) :: a
    real(sp) , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_r

  subroutine assignp1d_d(a,b)
    implicit none
    real(dp) , pointer , dimension(:) , intent(in) :: a
    real(dp) , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_d

  subroutine assignp1d_t(a,b)
    implicit none
    type(rcm_time_and_date) , pointer , dimension(:) , intent(in) :: a
    type(rcm_time_and_date) , pointer , dimension(:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):) => a
#endif
  end subroutine assignp1d_t

  subroutine assignp2d_l(a,b)
    implicit none
    logical , pointer , dimension(:,:) , intent(in) :: a
    logical , pointer , dimension(:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):) => a
#endif
  end subroutine assignp2d_l

  subroutine assignp2d_s(a,b)
    implicit none
    integer(2) , pointer , dimension(:,:) , intent(in) :: a
    integer(2) , pointer , dimension(:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):) => a
#endif
  end subroutine assignp2d_s

  subroutine assignp2d_i(a,b)
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: a
    integer , pointer , dimension(:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):) => a
#endif
  end subroutine assignp2d_i

  subroutine assignp2d_r(a,b)
    implicit none
    real(sp) , pointer , dimension(:,:) , intent(in) :: a
    real(sp) , pointer , dimension(:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):) => a
#endif
  end subroutine assignp2d_r

  subroutine assignp2d_d(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):) => a
#endif
  end subroutine assignp2d_d

  subroutine assignp3d_l(a,b)
    implicit none
    logical , pointer , dimension(:,:,:) , intent(in) :: a
    logical , pointer , dimension(:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):) => a
#endif
  end subroutine assignp3d_l

  subroutine assignp3d_s(a,b)
    implicit none
    integer(2) , pointer , dimension(:,:,:) , intent(in) :: a
    integer(2) , pointer , dimension(:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):) => a
#endif
  end subroutine assignp3d_s

  subroutine assignp3d_i(a,b)
    implicit none
    integer , pointer , dimension(:,:,:) , intent(in) :: a
    integer , pointer , dimension(:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):) => a
#endif
  end subroutine assignp3d_i

  subroutine assignp3d_r(a,b)
    implicit none
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: a
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):) => a
#endif
  end subroutine assignp3d_r

  subroutine assignp3d_d(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):) => a
#endif
  end subroutine assignp3d_d

  subroutine assignp4d_l(a,b)
    implicit none
    logical , pointer , dimension(:,:,:,:) , intent(in) :: a
    logical , pointer , dimension(:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a
#endif
  end subroutine assignp4d_l

  subroutine assignp4d_s(a,b)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:) , intent(in) :: a
    integer(2) , pointer , dimension(:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a
#endif
  end subroutine assignp4d_s

  subroutine assignp4d_i(a,b)
    implicit none
    integer , pointer , dimension(:,:,:,:) , intent(in) :: a
    integer , pointer , dimension(:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a
#endif
  end subroutine assignp4d_i

  subroutine assignp4d_r(a,b)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: a
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a
#endif
  end subroutine assignp4d_r

  subroutine assignp4d_d(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a
#endif
  end subroutine assignp4d_d

  subroutine assignp4d_3d_d(a,b,i)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: b
    integer , intent(in) :: i
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a(:,:,:,i)
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):) => a(:,:,:,i)
#endif
  end subroutine assignp4d_3d_d

  subroutine assignp5d_l(a,b)
    implicit none
    logical , pointer , dimension(:,:,:,:,:) , intent(in) :: a
    logical , pointer , dimension(:,:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):,lbound(a,5):) => a
#endif
  end subroutine assignp5d_l

  subroutine assignp5d_s(a,b)
    implicit none
    integer(2) , pointer , dimension(:,:,:,:,:) , intent(in) :: a
    integer(2) , pointer , dimension(:,:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):,lbound(a,5):) => a
#endif
  end subroutine assignp5d_s

  subroutine assignp5d_i(a,b)
    implicit none
    integer , pointer , dimension(:,:,:,:,:) , intent(in) :: a
    integer , pointer , dimension(:,:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):,lbound(a,5):) => a
#endif
  end subroutine assignp5d_i

  subroutine assignp5d_r(a,b)
    implicit none
    real(sp) , pointer , dimension(:,:,:,:,:) , intent(in) :: a
    real(sp) , pointer , dimension(:,:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):,lbound(a,5):) => a
#endif
  end subroutine assignp5d_r

  subroutine assignp5d_d(a,b)
    implicit none
    real(dp) , pointer , dimension(:,:,:,:,:) , intent(in) :: a
    real(dp) , pointer , dimension(:,:,:,:,:) , intent(out) :: b
    if ( .not. associated(a) ) then
      nullify(b)
      return
    end if
#ifndef __GFORTRAN__
    b => a
#else
    b(lbound(a,1):,lbound(a,2):,lbound(a,3):,lbound(a,4):,lbound(a,5):) => a
#endif
  end subroutine assignp5d_d

end module mod_memutil
