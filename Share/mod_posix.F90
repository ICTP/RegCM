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
      character(len=*, kind=c_char), intent(in) :: a
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
    subroutine seekdir(dir,pos) bind(C,name='seekdir')
      import
      implicit none
      type(c_ptr), value :: dir
      integer(c_long), intent(in) :: pos
    end subroutine seekdir
    function telldir(dir) bind(C,name='telldir')
      import
      implicit none
      type(c_ptr), value :: dir
      integer(c_long) :: telldir
    end function telldir
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

  public :: seekdir, telldir
  public :: dirlist, direntry
  public :: replacestr, lower, upper, splitstr, basename
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
    cname(1:256) = ' '
    iret = gethostname(a,ilen)
    hname = adjustl("")
    do i = 1, min(256,len(hname))
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

  pure recursive function replacestr(string,search,sub) result(mstring)
    implicit none
    character(len=*), intent(in) :: string, search, sub
    character(len=:), allocatable :: mstring
    integer :: i, stringlen, searchlen
    stringlen = len(string)
    searchlen = len(search)
    if ( stringlen == 0 .or. searchlen == 0 ) then
      mstring = ""
      return
    else if ( stringlen < searchlen ) then
      mstring = string
      return
    end if
    i = 1
    do
      if ( string(i:i+searchlen-1) == search ) then
        mstring = string(1:i-1)//sub// &
          replacestr(string(i+searchlen:stringlen),search,sub)
        exit
      end if
      if ( i+searchlen > stringlen ) then
        mstring = string
        exit
      end if
      i = i + 1
    end do
  end function replacestr

  elemental pure function lower(str,istart,istop) result (string)
    implicit none
    character(*), intent(in) :: str
    character(len(str)) :: string
    integer, intent(in), optional :: istart, istop
    integer :: i
    integer :: ibegin, iend
    string = str
    ibegin = 1
    if ( present(istart) ) then
      ibegin = max(ibegin,istart)
    end if
    iend = len_trim(str)
    if ( present(istop) ) then
      iend = min(iend,istop)
    end if
    do i = ibegin, iend
      select case (str(i:i))
        case ('A':'Z')
          string(i:i) = char(iachar(str(i:i))+32)
        case default
      end select
    end do
  end function lower

  elemental pure function upper(str,istart,istop) result (string)
    implicit none
    character(*), intent(in) :: str
    character(len(str)) :: string
    integer, intent(in), optional :: istart, istop
    integer :: i
    integer :: ibegin, iend
    string = str
    ibegin = 1
    if ( present(istart) ) then
      ibegin = max(ibegin,istart)
    end if
    iend = len_trim(str)
    if ( present(istop) ) then
      iend = min(iend,istop)
    end if
    do i = ibegin, iend
      select case (str(i:i))
        case ('a':'z')
          string(i:i) = char(iachar(str(i:i))-32)
        case default
      end select
    end do
  end function upper

  subroutine splitstr(input_line,array,delimiters,order,nulls)
    implicit none
    character(len=*), intent(in) :: input_line
    character(len=*), optional, intent(in) :: delimiters
    character(len=*), optional, intent(in) :: order
    character(len=*), optional, intent(in) :: nulls
    character(len=:), allocatable, intent(out) :: array(:)

    integer :: n
    integer, allocatable, dimension(:) :: ibegin
    integer, allocatable, dimension(:) :: iterm
    character(len=:), allocatable :: dlim
    character(len=:), allocatable :: ordr
    character(len=:), allocatable  :: nlls
    integer :: ii, iiii
    integer :: icount
    integer :: ilen
    integer :: i, j, k
    integer :: icol
    integer :: idlim
    integer :: ifound
    integer :: inotnull
    integer :: ireturn
    integer :: imax

    if ( present(delimiters) ) then
      if ( delimiters /= '' ) then
        dlim = delimiters
      else
        dlim = ' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)
      end if
    else
      dlim = ' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)
    end if
    idlim = len(dlim)

    if ( present(order) ) then
      ordr = lower(adjustl(order))
    else
      ordr='sequential'
    end if
    if ( present(nulls) ) then
      nlls = lower(adjustl(nulls))
    else
      nlls = 'ignore'
    end if

    n = len(input_line)+1
    allocate(ibegin(n))
    allocate(iterm(n))
    ibegin(:) = 1
    iterm(:) = 1

    ilen = len(input_line)
    icount = 0
    inotnull = 0
    imax = 0

    select case (ilen)
      case (0)
      case default
        icol = 1
        parseloop: do k = 1, ilen, 1
          ibegin(k) = icol
          if ( index(dlim(1:idlim),input_line(icol:icol)) == 0 ) then
            iterm(k) = ilen
            do i = 1, idlim
              ifound = index(input_line(ibegin(k):ilen),dlim(i:i))
              if ( ifound > 0 ) then
                iterm(k) = min(iterm(k),ifound+ibegin(k)-2)
              end if
            end do
            icol = iterm(k)+2
            inotnull = inotnull+1
          else
            iterm(k) = icol-1
            icol = icol+1
          end if
          imax = max(imax,iterm(k)-ibegin(k)+1)
          icount = k
          if ( icol > ilen ) then
            exit parseloop
          end if
        end do parseloop
    end select

    select case ( trim(adjustl(nlls)) )
      case ('ignore', '', 'ignoreend')
        ireturn = inotnull
      case default
        ireturn = icount
    end select
    allocate(character(len=imax) :: array(ireturn))
    select case ( trim(adjustl(ordr)) )
      case ('reverse', 'right')
        ii = ireturn
        iiii = -1
      case default
        ii = 1
        iiii = 1
    end select

    do j = 1, icount
      if ( iterm(j) < ibegin(j) ) then
        select case ( trim(adjustl(nlls)) )
          case ('ignore', '', 'ignoreend')
          case default
            array(ii) = ' '
            ii = ii+iiii
        end select
      else
        array(ii) = input_line(ibegin(j):iterm(j))
        ii = ii+iiii
      end if
    end do
  end subroutine splitstr

  function basename(path,suffix) result(base)
    implicit none
    character(*), intent(In) :: path
    logical, intent(in), optional :: suffix
    character(:), allocatable :: base
    character(:), allocatable :: file_parts(:)
    logical :: with_suffix

    if ( .not. present(suffix) ) then
      with_suffix = .true.
    else
      with_suffix = suffix
    end if
    call splitstr(path,file_parts,delimiters='\/')
    if ( size(file_parts) > 0 ) then
      base = trim(file_parts(size(file_parts)))
    else
      base = ''
    end if
    if ( .not. with_suffix ) then
      call splitstr(base,file_parts,delimiters='.')
      if ( size(file_parts) >= 2 ) then
        base = trim(file_parts(size(file_parts)-1))
      end if
    end if
  end function basename

end module mod_posix

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
