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

module mod_mpmessage

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message

  implicit none

  private

  public :: setup_mesg, die, aline, say, note, cry, fatal, checkalloc
  public :: vprntv, vprntm, iprntv

  interface vprntm
    module procedure vprntm_r8
    module procedure vprntm_r4
  end interface

  interface vprntv
    module procedure vprntv_r8
    module procedure vprntv_r4
  end interface vprntv

  contains
  !
  ! Printout helper
  !
  subroutine vprntv_r8(a,n,nam)
    implicit none
    integer(ik4), intent(in) :: n
    character(len=*), intent(in) :: nam
    real(rk8), intent(in), dimension(n) :: a
    integer(ik4) :: k, nn, kk
    integer(ik4), parameter :: npl = 7
    character(len=32) :: myf
    write(stdout,'(2x,a,a,a)') '# ',nam,' #'
    write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
    if ( n > npl ) then
      nn = n/npl
      do k = 1, nn
        kk = (k-1)*npl+1
        write(stdout,myf) a(kk:kk+npl-1)
      end do
      kk = n-(nn*npl)+1
    else
      kk = n
    end if
    if ( kk > 0 ) then
      write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
      write(stdout,myf) a(n-kk+1:)
    end if
  end subroutine vprntv_r8

  subroutine vprntv_r4(a,n,nam)
    implicit none
    integer(ik4), intent(in) :: n
    character(len=*), intent(in) :: nam
    real(rk4), intent(in), dimension(n) :: a
    integer(ik4) :: k, nn, kk
    integer(ik4), parameter :: npl = 7
    character(len=32) :: myf
    write(stdout,'(2x,a,a,a)') '# ',nam,' #'
    write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
    if ( n > npl ) then
      nn = n/npl
      do k = 1, nn
        kk = (k-1)*npl+1
        write(stdout,myf) a(kk:kk+npl-1)
      end do
      kk = n-(nn*npl)+1
    else
      kk = n
    end if
    if ( kk > 0 ) then
      write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
      write(stdout,myf) a(n-kk+1:)
    end if
  end subroutine vprntv_r4

  subroutine iprntv(a,n,nam)
    implicit none
    integer(ik4), intent(in) :: n
    character(len=*), intent(in) :: nam
    integer(ik4), intent(in), dimension(n) :: a
    integer(ik4) :: k, nn, kk
    integer(ik4), parameter :: npl = 7
    character(len=32) :: myf
    write(stdout,'(2x,a,a,a)') '# ',nam,' #'
    write(myf,'(a,i1,a)') '(2x,',npl,'i8)'
    if ( n > npl ) then
      nn = n/npl
      do k = 1, nn
        kk = (k-1)*npl+1
        write(stdout,myf) a(kk:kk+npl-1)
      end do
      kk = n-(nn*npl)+1
    else
      kk = n
    end if
    if ( kk > 0 ) then
      write(myf,'(a,i1,a)') '(2x,',kk,'i8)'
      write(stdout,myf) a(n-kk+1:)
    end if
  end subroutine iprntv
!
  subroutine vprntm_r8(a,n1,n2,nam)
    implicit none
    integer(ik4), intent (in) :: n1, n2
    character(len=*), intent (in) :: nam
    real(rk8), intent (in), dimension(n1,n2) :: a
    integer(ik4) :: k1, k2, nn, kk
    integer(ik4), parameter :: npl = 7
    character(len=32) :: myf
    write(stdout,'(2x,a,a,a)') '# ',nam,' #'
    nn = n1/npl
    do k2 = 1, n2
      write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
      write(stdout,*) '## Row ',k2
      if ( n1 > npl ) then
        do k1 = 1, nn
          kk = (k1-1)*npl+1
          write(stdout,myf) a(kk:kk+npl-1,k2)
        end do
        kk = n1-(nn*npl)+1
      else
        kk = n1
      end if
      if ( kk > 0 ) then
        write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
        write(stdout,myf) a(n1-kk+1:,k2)
      end if
      write(stdout,*) '## '
    end do
  end subroutine vprntm_r8

  subroutine vprntm_r4(a,n1,n2,nam)
    implicit none
    integer(ik4), intent (in) :: n1, n2
    character(len=*), intent (in) :: nam
    real(rk4), intent (in), dimension(n1,n2) :: a
    integer(ik4) :: k1, k2, nn, kk
    integer(ik4), parameter :: npl = 7
    character(len=32) :: myf
    write(stdout,'(2x,a,a,a)') '# ',nam,' #'
    nn = n1/npl
    do k2 = 1, n2
      write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
      write(stdout,*) '## Row ',k2
      if ( n1 > npl ) then
        do k1 = 1, nn
          kk = (k1-1)*npl+1
          write(stdout,myf) a(kk:kk+npl-1,k2)
        end do
        kk = n1-(nn*npl)+1
      else
        kk = n1
      end if
      if ( kk > 0 ) then
        write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
        write(stdout,myf) a(n1-kk+1:,k2)
      end if
      write(stdout,*) '## '
    end do
  end subroutine vprntm_r4
!
end module mod_mpmessage
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
