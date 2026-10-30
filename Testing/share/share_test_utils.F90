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

! Helpers shared by the unit tests of the Share library: every check prints
! PASS or FAIL, and test_end stops with a non zero exit status if any check
! failed, so that "make check" reports the test as failed.

module share_test_utils
  use mod_intkinds
  use mod_realkinds
  implicit none
  private
  public :: check, test_end, fixed_seed
  integer(ik4) :: nfail = 0
  integer(ik4) :: npass = 0
  contains

  subroutine check(ok,what)
    implicit none
    logical, intent(in) :: ok
    character(len=*), intent(in) :: what
    if ( ok ) then
      npass = npass + 1
      write(*,'(a)') 'PASS: '//what
    else
      nfail = nfail + 1
      write(*,'(a)') 'FAIL: '//what
    end if
  end subroutine check

  subroutine test_end(name)
    implicit none
    character(len=*), intent(in) :: name
    write(*,'(a,i0,a,i0,a)') name//': ', npass, ' passed, ', &
      nfail, ' failed'
    if ( nfail > 0 ) stop 1
  end subroutine test_end

  ! Same random sequence at every run.
  subroutine fixed_seed
    implicit none
    integer :: n, i
    integer, allocatable, dimension(:) :: s
    call random_seed(size=n)
    allocate(s(n))
    s = [(12345 + 7*i, i = 1, n)]
    call random_seed(put=s)
  end subroutine fixed_seed

end module share_test_utils

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort
