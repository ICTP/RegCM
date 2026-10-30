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

! Unit test of Share/mod_posix.F90: the host name is the one of the
! system, without the C string terminator; dirlist lists the files of a
! directory created in the current directory, and leaves its result not
! associated for a directory that does not exist.

program test_posix
  use mod_posix
  use share_test_utils
  implicit none
  character(len=64) :: hname, sysname
  type(direntry), pointer, dimension(:) :: de => null()
  integer :: i, ios, nfound
  logical :: la, lb

  call hostname(hname)
  call check(len_trim(hname) > 0 .and. index(hname,char(0)) == 0, &
    'host name without the C terminator')
  call execute_command_line('hostname > test_posix.host',exitstat=ios)
  if ( ios == 0 ) then
    open(10,file='test_posix.host',status='old',iostat=ios)
    if ( ios == 0 ) then
      read(10,'(a)',iostat=ios) sysname
      close(10,status='delete')
    end if
  end if
  if ( ios == 0 ) then
    call check(trim(hname) == trim(sysname),'host name of the system')
  end if

  call execute_command_line('rm -rf test_posix_dir && '// &
    'mkdir test_posix_dir && touch test_posix_dir/a.nc '// &
    'test_posix_dir/file_with_a_long_name.txt')
  call dirlist('test_posix_dir',de)
  call check(associated(de),'dirlist of an existing directory')
  if ( associated(de) ) then
    la = .false.
    lb = .false.
    nfound = 0
    do i = 1, size(de)
      if ( de(i)%ename == 'a.nc' ) la = .true.
      if ( de(i)%ename == 'file_with_a_long_name.txt' ) lb = .true.
      if ( de(i)%ename /= '.' .and. de(i)%ename /= '..' ) &
        nfound = nfound + 1
    end do
    call check(la .and. lb .and. nfound == 2,'dirlist file names')
    deallocate(de)
  end if
  call execute_command_line('rm -rf test_posix_dir')
  call dirlist('test_posix_does_not_exist',de)
  call check(.not. associated(de),'dirlist of a missing directory')

  call test_end('test_posix')
end program test_posix
