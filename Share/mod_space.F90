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

module mod_space

  use mod_intkinds
  use mod_realkinds
  use mod_date, only : rcm_time_and_date

  implicit none

  private

  public :: bounds
  public :: sarr1d, iarr1d, larr1d, r4arr1d, r8arr1d, rtarr1d
  public :: sarr2d, iarr2d, larr2d, r4arr2d, r8arr2d
  public :: sarr3d, iarr3d, larr3d, r4arr3d, r8arr3d
  public :: sarr4d, iarr4d, larr4d, r4arr4d, r8arr4d
  public :: sarr5d, iarr5d, larr5d, r4arr5d, r8arr5d

  public :: getspc

  type bounds
    integer(ik4) :: low
    integer(ik4) :: high
  end type bounds

  type sarr1d
    integer(ik2), pointer, contiguous, dimension(:) :: space => null()
  end type sarr1d

  type iarr1d
    integer(ik4), pointer, contiguous, dimension(:) :: space => null()
  end type iarr1d

  type larr1d
    logical, pointer, contiguous, dimension(:) :: space => null()
  end type larr1d

  type r4arr1d
    real(rk4), pointer, contiguous, dimension(:) :: space => null()
  end type r4arr1d

  type r8arr1d
    real(rk8), pointer, contiguous, dimension(:) :: space => null()
  end type r8arr1d

  type rtarr1d
    type(rcm_time_and_date), pointer, contiguous, dimension(:) :: space => null()
  end type rtarr1d

  type sarr2d
    integer(ik2), pointer, contiguous, dimension(:,:) :: space => null()
  end type sarr2d

  type iarr2d
    integer(ik4), pointer, contiguous, dimension(:,:) :: space => null()
  end type iarr2d

  type larr2d
    logical, pointer, contiguous, dimension(:,:) :: space => null()
  end type larr2d

  type r4arr2d
    real(rk4), pointer, contiguous, dimension(:,:) :: space => null()
  end type r4arr2d

  type r8arr2d
    real(rk8), pointer, contiguous, dimension(:,:) :: space => null()
  end type r8arr2d

  type sarr3d
    integer(ik2), pointer, contiguous, dimension(:,:,:) :: space => null()
  end type sarr3d

  type iarr3d
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: space => null()
  end type iarr3d

  type larr3d
    logical, pointer, contiguous, dimension(:,:,:) :: space => null()
  end type larr3d

  type r4arr3d
    real(rk4), pointer, contiguous, dimension(:,:,:) :: space => null()
  end type r4arr3d

  type r8arr3d
    real(rk8), pointer, contiguous, dimension(:,:,:) :: space => null()
  end type r8arr3d

  type sarr4d
    integer(ik2), pointer, contiguous, dimension(:,:,:,:) :: space => null()
  end type sarr4d

  type iarr4d
    integer(ik4), pointer, contiguous, dimension(:,:,:,:) :: space => null()
  end type iarr4d

  type larr4d
    logical, pointer, contiguous, dimension(:,:,:,:) :: space => null()
  end type larr4d

  type r4arr4d
    real(rk4), pointer, contiguous, dimension(:,:,:,:) :: space => null()
  end type r4arr4d

  type r8arr4d
    real(rk8), pointer, contiguous, dimension(:,:,:,:) :: space => null()
  end type r8arr4d

  type sarr5d
    integer(ik2), pointer, contiguous, dimension(:,:,:,:,:) :: space => null()
  end type sarr5d

  type iarr5d
    integer(ik4), pointer, contiguous, dimension(:,:,:,:,:) :: space => null()
  end type iarr5d

  type larr5d
    logical, pointer, contiguous, dimension(:,:,:,:,:) :: space => null()
  end type larr5d

  type r4arr5d
    real(rk4), pointer, contiguous, dimension(:,:,:,:,:) :: space => null()
  end type r4arr5d

  type r8arr5d
    real(rk8), pointer, contiguous, dimension(:,:,:,:,:) :: space => null()
  end type r8arr5d

  interface getspc
    module procedure getspc1d_larr
    module procedure getspc1d_sarr
    module procedure getspc1d_iarr
    module procedure getspc1d_r4arr
    module procedure getspc1d_r8arr
    module procedure getspc1d_rtarr
    module procedure getspc2d_larr
    module procedure getspc2d_sarr
    module procedure getspc2d_iarr
    module procedure getspc2d_r4arr
    module procedure getspc2d_r8arr
    module procedure getspc3d_larr
    module procedure getspc3d_sarr
    module procedure getspc3d_iarr
    module procedure getspc3d_r4arr
    module procedure getspc3d_r8arr
    module procedure getspc4d_larr
    module procedure getspc4d_sarr
    module procedure getspc4d_iarr
    module procedure getspc4d_r4arr
    module procedure getspc4d_r8arr
    module procedure getspc5d_larr
    module procedure getspc5d_sarr
    module procedure getspc5d_iarr
    module procedure getspc5d_r4arr
    module procedure getspc5d_r8arr
  end interface getspc

  contains

  subroutine getspc1d_larr(a,b,istat)
    type (larr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat, source=.false.)
  end subroutine getspc1d_larr

  subroutine getspc1d_sarr(a,b,istat)
    type (sarr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat, source=-1_ik2)
  end subroutine getspc1d_sarr

  subroutine getspc1d_iarr(a,b,istat)
    type (iarr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat, source=-1)
  end subroutine getspc1d_iarr

  subroutine getspc1d_r4arr(a,b,istat)
    type (r4arr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat, source=0.0)
  end subroutine getspc1d_r4arr

  subroutine getspc1d_r8arr(a,b,istat)
    type (r8arr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat, source=0.0d0)
  end subroutine getspc1d_r8arr

  subroutine getspc1d_rtarr(a,b,istat)
    type (rtarr1d), intent(inout) :: a
    type (bounds), intent(in) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat)
  end subroutine getspc1d_rtarr

  subroutine getspc2d_larr(a,b,istat)
    type (larr2d), intent(inout) :: a
    type (bounds), intent(in), dimension(2) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), &
             stat=istat, source=.false.)
  end subroutine getspc2d_larr

  subroutine getspc2d_sarr(a,b,istat)
    type (sarr2d), intent(inout) :: a
    type (bounds), intent(in), dimension(2) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), &
             stat=istat, source=-1_ik2)
  end subroutine getspc2d_sarr

  subroutine getspc2d_iarr(a,b,istat)
    type (iarr2d), intent(inout) :: a
    type (bounds), intent(in), dimension(2) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), &
             stat=istat, source=-1)
  end subroutine getspc2d_iarr

  subroutine getspc2d_r4arr(a,b,istat)
    type (r4arr2d), intent(inout) :: a
    type (bounds), intent(in), dimension(2) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), &
             stat=istat, source=0.0)
  end subroutine getspc2d_r4arr

  subroutine getspc2d_r8arr(a,b,istat)
    type (r8arr2d), intent(inout) :: a
    type (bounds), intent(in), dimension(2) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), &
             stat=istat, source=0.0d0)
  end subroutine getspc2d_r8arr

  subroutine getspc3d_larr(a,b,istat)
    type (larr3d), intent(inout) :: a
    type (bounds), intent(in), dimension(3) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat, source=.false.)
  end subroutine getspc3d_larr

  subroutine getspc3d_sarr(a,b,istat)
    type (sarr3d), intent(inout) :: a
    type (bounds), intent(in), dimension(3) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat, source=-1_ik2)
  end subroutine getspc3d_sarr

  subroutine getspc3d_iarr(a,b,istat)
    type (iarr3d), intent(inout) :: a
    type (bounds), intent(in), dimension(3) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat, source=-1)
  end subroutine getspc3d_iarr

  subroutine getspc3d_r4arr(a,b,istat)
    type (r4arr3d), intent(inout) :: a
    type (bounds), intent(in), dimension(3) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat, source=0.0)
  end subroutine getspc3d_r4arr

  subroutine getspc3d_r8arr(a,b,istat)
    type (r8arr3d), intent(inout) :: a
    type (bounds), intent(in), dimension(3) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat, source=0.0d0)
  end subroutine getspc3d_r8arr

  subroutine getspc4d_larr(a,b,istat)
    type (larr4d), intent(inout) :: a
    type (bounds), intent(in), dimension(4) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), &
                     stat=istat, source=.false.)
  end subroutine getspc4d_larr

  subroutine getspc4d_sarr(a,b,istat)
    type (sarr4d), intent(inout) :: a
    type (bounds), intent(in), dimension(4) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), &
                     stat=istat, source=-1_ik2)
  end subroutine getspc4d_sarr

  subroutine getspc4d_iarr(a,b,istat)
    type (iarr4d), intent(inout) :: a
    type (bounds), intent(in), dimension(4) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), &
                     stat=istat, source=-1)
  end subroutine getspc4d_iarr

  subroutine getspc4d_r4arr(a,b,istat)
    type (r4arr4d), intent(inout) :: a
    type (bounds), intent(in), dimension(4) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), &
                     stat=istat, source=0.0)
  end subroutine getspc4d_r4arr

  subroutine getspc4d_r8arr(a,b,istat)
    type (r8arr4d), intent(inout) :: a
    type (bounds), intent(in), dimension(4) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), &
                     stat=istat, source=0.0d0)
  end subroutine getspc4d_r8arr

  subroutine getspc5d_larr(a,b,istat)
    type (larr5d), intent(inout) :: a
    type (bounds), intent(in), dimension(5) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat, source=.false.)
  end subroutine getspc5d_larr

  subroutine getspc5d_sarr(a,b,istat)
    type (sarr5d), intent(inout) :: a
    type (bounds), intent(in), dimension(5) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat, source=-1_ik2)
  end subroutine getspc5d_sarr

  subroutine getspc5d_iarr(a,b,istat)
    type (iarr5d), intent(inout) :: a
    type (bounds), intent(in), dimension(5) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat, source=-1)
  end subroutine getspc5d_iarr

  subroutine getspc5d_r4arr(a,b,istat)
    type (r4arr5d), intent(inout) :: a
    type (bounds), intent(in), dimension(5) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat, source=0.0)
  end subroutine getspc5d_r4arr

  subroutine getspc5d_r8arr(a,b,istat)
    type (r8arr5d), intent(inout) :: a
    type (bounds), intent(in), dimension(5) ::  b
    integer(ik4), intent(out) :: istat
    if ( associated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat, source=0.0d0)
  end subroutine getspc5d_r8arr

end module mod_space
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
