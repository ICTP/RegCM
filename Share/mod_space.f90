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

module mod_space

  use m_realkinds

  private

  public :: bounds
  public :: iarr1d , larr1d , r4arr1d , r8arr1d
  public :: iarr2d , larr2d , r4arr2d , r8arr2d
  public :: iarr3d , larr3d , r4arr3d , r8arr3d
  public :: iarr4d , larr4d , r4arr4d , r8arr4d
  public :: iarr5d , larr5d , r4arr5d , r8arr5d

  public :: getspc1d , getspc2d , getspc3d , getspc4d , getspc5d

  type bounds
    integer :: low
    integer :: high
  end type bounds

  type iarr1d
    integer , allocatable , dimension(:) :: space
  end type iarr1d
 
  type larr1d
    logical , allocatable , dimension(:) :: space
  end type larr1d

  type r4arr1d
    real(sp) , allocatable , dimension(:) :: space
  end type r4arr1d

  type r8arr1d
    real(dp) , allocatable , dimension(:) :: space
  end type r8arr1d

  type iarr2d
    integer , allocatable , dimension(:,:) :: space
  end type iarr2d
 
  type larr2d
    logical , allocatable , dimension(:,:) :: space
  end type larr2d

  type r4arr2d
    real(sp) , allocatable , dimension(:,:) :: space
  end type r4arr2d

  type r8arr2d
    real(dp) , allocatable , dimension(:,:) :: space
  end type r8arr2d

  type iarr3d
    integer , allocatable , dimension(:,:,:) :: space
  end type iarr3d
 
  type larr3d
    logical , allocatable , dimension(:,:,:) :: space
  end type larr3d

  type r4arr3d
    real(sp) , allocatable , dimension(:,:,:) :: space
  end type r4arr3d

  type r8arr3d
    real(dp) , allocatable , dimension(:,:,:) :: space
  end type r8arr3d

  type iarr4d
    integer , allocatable , dimension(:,:,:,:) :: space
  end type iarr4d
 
  type larr4d
    logical , allocatable , dimension(:,:,:,:) :: space
  end type larr4d

  type r4arr4d
    real(sp) , allocatable , dimension(:,:,:,:) :: space
  end type r4arr4d

  type r8arr4d
    real(dp) , allocatable , dimension(:,:,:,:) :: space
  end type r8arr4d

  type iarr5d
    integer , allocatable , dimension(:,:,:,:,:) :: space
  end type iarr5d
 
  type larr5d
    logical , allocatable , dimension(:,:,:,:,:) :: space
  end type larr5d

  type r4arr5d
    real(sp) , allocatable , dimension(:,:,:,:,:) :: space
  end type r4arr5d

  type r8arr5d
    real(dp) , allocatable , dimension(:,:,:,:,:) :: space
  end type r8arr5d

  interface getspc1d
    module procedure getspc1d_larr
    module procedure getspc1d_iarr
    module procedure getspc1d_r4arr
    module procedure getspc1d_r8arr
  end interface getspc1d

  interface getspc2d
    module procedure getspc2d_larr
    module procedure getspc2d_iarr
    module procedure getspc2d_r4arr
    module procedure getspc2d_r8arr
  end interface getspc2d

  interface getspc3d
    module procedure getspc3d_larr
    module procedure getspc3d_iarr
    module procedure getspc3d_r4arr
    module procedure getspc3d_r8arr
  end interface getspc3d

  interface getspc4d
    module procedure getspc4d_larr
    module procedure getspc4d_iarr
    module procedure getspc4d_r4arr
    module procedure getspc4d_r8arr
  end interface getspc4d

  interface getspc5d
    module procedure getspc5d_larr
    module procedure getspc5d_iarr
    module procedure getspc5d_r4arr
    module procedure getspc5d_r8arr
  end interface getspc5d

  contains

  subroutine getspc1d_larr(a,b,istat)
    type (larr1d) , intent(inout) :: a
    type (bounds) , intent(in) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat)
  end subroutine getspc1d_larr

  subroutine getspc1d_iarr(a,b,istat)
    type (iarr1d) , intent(inout) :: a
    type (bounds) , intent(in) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat)
  end subroutine getspc1d_iarr

  subroutine getspc1d_r4arr(a,b,istat)
    type (r4arr1d) , intent(inout) :: a
    type (bounds) , intent(in) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat)
  end subroutine getspc1d_r4arr

  subroutine getspc1d_r8arr(a,b,istat)
    type (r8arr1d) , intent(inout) :: a
    type (bounds) , intent(in) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b%low:b%high), stat=istat)
  end subroutine getspc1d_r8arr

  subroutine getspc2d_larr(a,b,istat)
    type (larr2d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(2) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), stat=istat)
  end subroutine getspc2d_larr

  subroutine getspc2d_iarr(a,b,istat)
    type (iarr2d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(2) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), stat=istat)
  end subroutine getspc2d_iarr

  subroutine getspc2d_r4arr(a,b,istat)
    type (r4arr2d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(2) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), stat=istat)
  end subroutine getspc2d_r4arr

  subroutine getspc2d_r8arr(a,b,istat)
    type (r8arr2d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(2) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high), stat=istat)
  end subroutine getspc2d_r8arr

  subroutine getspc3d_larr(a,b,istat)
    type (larr3d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(3) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat)
  end subroutine getspc3d_larr

  subroutine getspc3d_iarr(a,b,istat)
    type (iarr3d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(3) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat)
  end subroutine getspc3d_iarr

  subroutine getspc3d_r4arr(a,b,istat)
    type (r4arr3d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(3) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat)
  end subroutine getspc3d_r4arr

  subroutine getspc3d_r8arr(a,b,istat)
    type (r8arr3d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(3) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high), stat=istat)
  end subroutine getspc3d_r8arr

  subroutine getspc4d_larr(a,b,istat)
    type (larr4d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(4) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), stat=istat)
  end subroutine getspc4d_larr

  subroutine getspc4d_iarr(a,b,istat)
    type (iarr4d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(4) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), stat=istat)
  end subroutine getspc4d_iarr

  subroutine getspc4d_r4arr(a,b,istat)
    type (r4arr4d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(4) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), stat=istat)
  end subroutine getspc4d_r4arr

  subroutine getspc4d_r8arr(a,b,istat)
    type (r8arr4d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(4) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high), stat=istat)
  end subroutine getspc4d_r8arr

  subroutine getspc5d_larr(a,b,istat)
    type (larr5d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(5) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat)
  end subroutine getspc5d_larr

  subroutine getspc5d_iarr(a,b,istat)
    type (iarr5d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(5) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat)
  end subroutine getspc5d_iarr

  subroutine getspc5d_r4arr(a,b,istat)
    type (r4arr5d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(5) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat)
  end subroutine getspc5d_r4arr

  subroutine getspc5d_r8arr(a,b,istat)
    type (r8arr5d) , intent(inout) :: a
    type (bounds) , intent(in) , dimension(5) ::  b
    integer , intent(out) :: istat
    if ( allocated(a%space) ) deallocate(a%space)
    allocate(a%space(b(1)%low:b(1)%high,b(2)%low:b(2)%high, &
                     b(3)%low:b(3)%high,b(4)%low:b(4)%high, &
                     b(5)%low:b(5)%high), stat=istat)
  end subroutine getspc5d_r8arr

end module mod_space
