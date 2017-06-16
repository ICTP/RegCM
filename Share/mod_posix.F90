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

module mod_posix

  use mod_stdio
  use iso_c_binding

  implicit none

  private

  type , bind(C) :: dirent
    integer(c_long) :: d_ino
    integer(c_long) :: d_off ! __off_t, check size
    integer(c_short) :: d_reclen
    character(len=1,kind=c_char) :: d_type
    character(len=1,kind=c_char) :: d_name(256)
  end type

  interface
    function opendir(a) bind(C,name='opendir')
      import
      type(c_ptr) :: opendir
      character(len=1, kind=c_char) , intent(in) :: a(*)
    end function opendir
    function readdir(dir) bind(C,name='readdir')
      import
      type(c_ptr) , value :: dir
      type(c_ptr) :: readdir
    end function readdir
    function closedir(dir) bind(C,name='closedir')
      import
      type(c_ptr) , value :: dir
      integer(c_int) :: closedir
    end function closedir
    subroutine seekdir(dir,pos) bind(C,name='seekdir')
      import
      type(c_ptr) , value :: dir
      integer(c_long) :: pos
    end subroutine seekdir
    function telldir(dir) bind(C,name='telldir')
      import
      type(c_ptr) , value :: dir
      integer(c_long) :: telldir
    end function telldir
    subroutine rewinddir(dir) bind(C,name='rewinddir')
      import
      type(c_ptr) , value :: dir
    end subroutine rewinddir
  end interface

  type :: direntry
    character(len=256) :: ename
  end type direntry

  public :: dirlist , direntry

  contains

  subroutine dirlist(path,dire)
    use iso_c_binding
    implicit none
    character(len=*) , intent(in) :: path
    type(direntry) , dimension(:) , pointer :: dire
    type(c_ptr) :: dir , dc
    integer(c_int) :: ires
    type(dirent) , pointer :: d
    integer :: i , ic , fnum

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
    do i = 1 , fnum
      dc = readdir(dir)
      if ( .not. c_associated(dc) ) exit
      call c_f_pointer(dc,d)
      charloop: &
      do ic = 1 , 256
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
