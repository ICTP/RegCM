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

module mod_posix

  use mod_stdio
  use, intrinsic :: iso_c_binding

  implicit none

  private

  type, bind(C) :: dirent
    integer(c_long) :: d_ino
    integer(c_long) :: d_off ! __off_t, check size
    integer(c_short) :: d_reclen
    character(len=1,kind=c_char) :: d_type
    character(len=1,kind=c_char) :: d_name(256)
  end type dirent

  interface
    function opendir(a) bind(C,name='opendir')
      import
      implicit none
      type(c_ptr) :: opendir
      ! A C string: the path terminated by C_NULL_CHAR
      character(kind=c_char), dimension(*), intent(in) :: a
    end function opendir
    function readdir(dir) bind(C,name='readdir')
      import
      implicit none
      type(c_ptr), value :: dir
      type(c_ptr) :: readdir
    end function readdir
    function closedir(dir) bind(C,name='closedir')
      import
      implicit none
      type(c_ptr), value :: dir
      integer(c_int) :: closedir
    end function closedir
    subroutine rewinddir(dir) bind(C,name='rewinddir')
      import
      implicit none
      type(c_ptr), value :: dir
    end subroutine rewinddir
  end interface

  type :: direntry
    character(len=256) :: ename
  end type direntry

  interface
    integer(kind=c_int) function gethostname(a,i) bind(C,name='gethostname')
      import
      implicit none
      type(c_ptr), value :: a
      integer(kind=c_size_t), intent(in), value :: i
    end function gethostname
  end interface

  public :: dirlist, direntry
  public :: hostname

  contains

  subroutine hostname(hname)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), intent(inout) :: hname
    integer(kind=c_size_t) :: ilen
    integer(kind=c_int) :: iret
    character(kind=c_char), target :: cname(256)
    type(c_ptr) :: a
    integer :: i
    a = c_loc(cname)
    ilen = 256
    cname(1:256) = C_NULL_CHAR
    iret = gethostname(a,ilen)
    hname = ''
    if ( iret /= 0 ) then
      hname = 'unknown'
      return
    end if
    do i = 1, min(256,len(hname))
      if ( cname(i) == C_NULL_CHAR ) exit
      hname(i:i) = cname(i)
    end do
  end subroutine hostname

  subroutine dirlist(path,dire)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), intent(in) :: path
    type(direntry), dimension(:), pointer, intent(inout) :: dire
    type(c_ptr) :: dir, dc
    integer(c_int) :: ires
    type(dirent), pointer :: d
    integer :: i, ic, fnum

    if ( associated(dire) ) then
      deallocate(dire)
      nullify(dire)
    end if
    dir = opendir(trim(path)//C_NULL_CHAR)
    if ( .not.c_associated(dir) ) then
      write(stderr,*) 'Cannot open directory ',trim(path)
      return
    end if
    fnum = 0
    do
      dc = readdir(dir)
      if ( .not. c_associated(dc) ) exit
      fnum = fnum + 1
    end do
    if ( fnum == 0 ) then
      write(stderr,*) 'Cannot read directory ',trim(path)
      ires = closedir(dir)
      return
    end if
    call rewinddir(dir)
    allocate(dire(fnum))
    do i = 1, fnum
      dire(i)%ename = ' '
    end do
    do i = 1, fnum
      dc = readdir(dir)
      if ( .not. c_associated(dc) ) exit
      call c_f_pointer(dc,d)
      charloop: &
      do ic = 1, 256
        if ( d%d_name(ic) == C_NULL_CHAR ) then
          dire(i)%ename(ic:) = ' '
          exit charloop
        end if
        dire(i)%ename(ic:ic) = d%d_name(ic)
      end do charloop
    end do
    ires = closedir(dir)
  end subroutine dirlist

end module mod_posix

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
