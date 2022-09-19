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

module mod_sort
  use mod_realkinds
  use mod_intkinds

  implicit none

  private

  interface sort
    module procedure sort_int32
    module procedure sort_int64
    module procedure sort_real32
    module procedure sort_real64
  end interface sort

  interface quicksort
    module procedure quicksort_int32
    module procedure quicksort_int64
    module procedure quicksort_real32
    module procedure quicksort_real64
  end interface quicksort

  interface msi_index
    module procedure msi_r8
    module procedure msi_r4
    module procedure msi_i4
  end interface msi_index

  interface argsort
    module procedure argsort_r8
    module procedure argsort_r4
    module procedure argsort_i4
  end interface argsort

  public :: sort , argsort , msi_index

  contains

  subroutine sort_int32(x)
    implicit none
    integer(kind=ik4) , dimension(:) , intent(inout) :: x
    integer :: istart, istop
    istart = 1
    istop = size(x)
    if ( istop < 2 ) return
    call quicksort_int32(x,istart,istop)
  end subroutine sort_int32

  subroutine sort_int64(x)
    implicit none
    integer(kind=ik8) , dimension(:) , intent(inout) :: x
    integer :: istart, istop
    istart = 1
    istop = size(x)
    if ( istop < 2 ) return
    call quicksort_int64(x,istart,istop)
  end subroutine sort_int64

  subroutine sort_real32(x)
    implicit none
    real(kind=rk4) , dimension(:) , intent(inout) :: x
    integer :: istart, istop
    istart = 1
    istop = size(x)
    if ( istop < 2 ) return
    call quicksort_real32(x,istart,istop)
  end subroutine sort_real32

  subroutine sort_real64(x)
    implicit none
    real(kind=rk8) , dimension(:) , intent(inout) :: x
    integer :: istart, istop
    istart = 1
    istop = size(x)
    if ( istop < 2 ) return
    call quicksort_real64(x,istart,istop)
  end subroutine sort_real64

   recursive subroutine quicksort_int32(x,first,last)
    implicit none
    integer(kind=ik4) , dimension(:) , intent(inout) :: x
    integer , intent(in) :: first, last
    integer(kind=ik4) :: pivot, temp
    integer :: left, right

    if ( first >= last ) return
    pivot = x((first+last)/2)
    left = first
    right = last

    do while ( left <= right )
      do while ( x(left) < pivot )
        left = left + 1
      end do
      do while ( x(right) > pivot )
        right = right - 1
      end do
      if ( left <= right ) then
        temp = x(left)
        x(left) = x(right)
        x(right) = temp
        left = left + 1
        right = right - 1
      end if
      call quicksort_int32(x,first,right)
      call quicksort_int32(x,left,last)
    end do
  end subroutine quicksort_int32

  recursive subroutine quicksort_int64(x,first,last)
    implicit none
    integer(kind=ik8) , dimension(:) , intent(inout) :: x
    integer , intent(in) :: first, last
    integer(kind=ik8) :: pivot, temp
    integer :: left, right

    if ( first >= last ) return
    pivot = x((first+last)/2)
    left = first
    right = last

    do while ( left <= right )
      do while ( x(left) < pivot )
        left = left + 1
      end do
      do while ( x(right) > pivot )
        right = right - 1
      end do
      if ( left <= right ) then
        temp = x(left)
        x(left) = x(right)
        x(right) = temp
        left = left + 1
        right = right - 1
      end if
      call quicksort_int64(x,first,right)
      call quicksort_int64(x,left,last)
    end do
  end subroutine quicksort_int64

  recursive subroutine quicksort_real32(x,first,last)
    implicit none
    real(kind=rk4) , dimension(:) , intent(inout) :: x
    integer , intent(in) :: first, last
    real(kind=rk4) :: pivot, temp
    integer :: left, right

    if ( first >= last ) return
    pivot = x((first+last)/2)
    left = first
    right = last

    do while ( left <= right )
      do while ( x(left) < pivot )
        left = left + 1
      end do
      do while ( x(right) > pivot )
        right = right - 1
      end do
      if ( left <= right ) then
        temp = x(left)
        x(left) = x(right)
        x(right) = temp
        left = left + 1
        right = right - 1
      end if
      call quicksort_real32(x,first,right)
      call quicksort_real32(x,left,last)
    end do
  end subroutine quicksort_real32

  recursive subroutine quicksort_real64(x,first,last)
    implicit none
    real(kind=rk8) , dimension(:) , intent(inout) :: x
    integer , intent(in) :: first, last
    real(kind=rk8) :: pivot, temp
    integer :: left, right

    if ( first >= last ) return
    pivot = x((first+last)/2)
    left = first
    right = last

    do while ( left <= right )
      do while ( x(left) < pivot )
        left = left + 1
      end do
      do while ( x(right) > pivot )
        right = right - 1
      end do
      if ( left <= right ) then
        temp = x(left)
        x(left) = x(right)
        x(right) = temp
        left = left + 1
        right = right - 1
      end if
      call quicksort_real64(x,first,right)
      call quicksort_real64(x,left,last)
    end do
  end subroutine quicksort_real64

  subroutine msi_r8(a,idx,jdx)
    implicit none
    real(rk8) , dimension(:,:) , intent(in) :: a
    integer(ik4) , dimension(:,:) , intent(out) :: idx , jdx
    real(rk8) , dimension(:) , allocatable :: b
    integer(ik4), dimension(product(shape(a))) :: iord
    integer(ik4) :: nn , n , i , j , ii , jj , nx , ny
    integer(ik4) :: is , js
    nn = product(shape(a))
    b = reshape(a,(/nn/))
    iord = argsort(b)
    nx = size(a,1)
    ny = size(a,2)
    is = lbound(a,1)
    js = lbound(a,2)
    do n = 1 , nn
      j = js + (iord(n)-1)/nx
      i = is + (iord(n)-1)-(j-js)*nx
      jj = 1 + (n-1)/nx
      ii = 1 + (n-1)-(jj-1)*nx
      idx(ii,jj) = i
      jdx(ii,jj) = j
    end do
  end subroutine msi_r8

  subroutine msi_r4(a,idx,jdx)
    implicit none
    real(rk4) , dimension(:,:) , intent(in) :: a
    integer(ik4) , dimension(:,:) , intent(out) :: idx , jdx
    real(rk4) , dimension(:) , allocatable :: b
    integer(ik4), dimension(product(shape(a))) :: iord
    integer(ik4) :: nn , n , i , j , ii , jj , nx , ny
    integer(ik4) :: is , js
    nn = product(shape(a))
    b = reshape(a,(/nn/))
    iord = argsort(b)
    nx = size(a,1)
    ny = size(a,2)
    is = lbound(a,1)
    js = lbound(a,2)
    do n = 1 , nn
      j = js + (iord(n)-1)/nx
      i = is + (iord(n)-1)-(j-js)*nx
      jj = 1 + (n-1)/nx
      ii = 1 + (n-1)-(jj-1)*nx
      idx(ii,jj) = i
      jdx(ii,jj) = j
    end do
  end subroutine msi_r4

  subroutine msi_i4(a,idx,jdx)
    implicit none
    integer(ik4) , dimension(:,:) , intent(in) :: a
    integer(ik4) , dimension(:,:) , intent(out) :: idx , jdx
    integer(ik4) , dimension(:) , allocatable :: b
    integer(ik4), dimension(product(shape(a))) :: iord
    integer(ik4) :: nn , n , i , j , ii , jj , nx , ny
    integer(ik4) :: is , js
    nn = product(shape(a))
    b = reshape(a,(/nn/))
    iord = argsort(b)
    nx = size(a,1)
    ny = size(a,2)
    is = lbound(a,1)
    js = lbound(a,2)
    do n = 1 , nn
      j = js + (iord(n)-1)/nx
      i = is + (iord(n)-1)-(j-js)*nx
      jj = 1 + (n-1)/nx
      ii = 1 + (n-1)-(jj-1)*nx
      idx(ii,jj) = i
      jdx(ii,jj) = j
    end do
  end subroutine msi_i4

  function argsort_r8(a) result(b)
    implicit none
    real(rk8) , intent(in) :: a(:)
    integer(ik4) , dimension(size(a)) :: b
    integer :: n , i , imin , temp1
    real(rk8) :: temp2
    real(rk8) , dimension(size(a)) :: a2
    a2 = a
    n = size(a)
    do i = 1 , n
      b(i) = i
    end do
    if ( n == 1 ) return
    do i = 1 , n-1
      imin = minloc(a2(i:),1) + i - 1
      if ( imin /= i ) then
        temp2 = a2(i)
        a2(i) = a2(imin)
        a2(imin) = temp2
        temp1 = b(i)
        b(i) = b(imin)
        b(imin) = temp1
      end if
    end do
  end function argsort_r8

  function argsort_r4(a) result(b)
    implicit none
    real(rk4) , intent(in) :: a(:)
    integer(ik4) , dimension(size(a)) :: b
    integer :: n , i , imin , temp1
    real(rk4) :: temp2
    real(rk4) , dimension(size(a)) :: a2
    a2 = a
    n = size(a)
    do i = 1 , n
      b(i) = i
    end do
    if ( n == 1 ) return
    do i = 1 , n-1
      imin = minloc(a2(i:),1) + i - 1
      if ( imin /= i ) then
        temp2 = a2(i)
        a2(i) = a2(imin)
        a2(imin) = temp2
        temp1 = b(i)
        b(i) = b(imin)
        b(imin) = temp1
      end if
    end do
  end function argsort_r4

  function argsort_i4(a) result(b)
    implicit none
    integer(ik4) , intent(in) :: a(:)
    integer(ik4) , dimension(size(a)) :: b
    integer :: n , i , imin , temp1
    integer(ik4) :: temp2
    integer(ik4) , dimension(size(a)) :: a2
    a2 = a
    n = size(a)
    do i = 1 , n
      b(i) = i
    end do
    if ( n == 1 ) return
    do i = 1 , n-1
      imin = minloc(a2(i:),1) + i - 1
      if ( imin /= i ) then
        temp2 = a2(i)
        a2(i) = a2(imin)
        a2(imin) = temp2
        temp1 = b(i)
        b(i) = b(imin)
        b(imin) = temp1
      end if
    end do
  end function argsort_i4

end module mod_sort

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
