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

module mod_hash

  use mod_intkinds
  implicit none

  private

  !integer(ik4), parameter :: magic_numb = z'5d7a9f43'
  integer(ik4), parameter :: magic_numb = 1568317251

  public :: hash, bsearch

  contains

  !
  ! Search in ordered list
  !
  recursive integer(ik4) function bsearch(l,v,n) result(k)
    implicit none
    integer(ik4), dimension(:), intent(in) :: l
    integer(ik4), intent(in) :: v
    integer(ik4), intent(in) :: n
    integer(ik4) :: i, j, p
    i = 1
    j = n
    k = -1
    if ( n > 0 ) then
      do while ( i <= j )
        p = (i+j)/2
        if ( l(p) == v ) then
          k = p
          exit
        else if ( l(p) < v ) then
          i = p + 1
        else
          j = p - 1
        end if
      end do
    end if
  end function bsearch

  integer(ik4) function hash(text) result(hashed)
    implicit none
    character(len=*), intent(in) :: text
    integer(ik4) :: i, j
    hashed = 0
    do i = 1, len_trim(text)
      j = mod(i-1, 4) * 8
      hashed = ieor( hashed, ishft( ichar( text(i:i) ), j ) )
    end do
    hashed = abs( ieor( hashed, magic_numb ) )
  end function hash

end module mod_hash

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
