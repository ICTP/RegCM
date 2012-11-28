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

module mod_mpmessage

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message

  private
!
  public :: setup_mesg , die , aline , say , note , cry , fatal , checkalloc
  public :: vprntv , vprntm

  contains
  !
  ! Printout helper
  !
  subroutine vprntv(a,n,nam)
    implicit none
    integer(ik4) , intent(in) :: n
    character(len=*) , intent(in) :: nam
    real(rk8) , intent(in) , dimension(n) :: a
    integer(ik4) :: k , nn , kk
    integer(ik4) , parameter :: npl = 7
    character(len=32) :: myf
    write(stderr,'(2x,a,a,a)') '# ',nam,' #'
    write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
    if ( n > npl ) then
      nn = n/npl
      do k = 1 , nn
        kk = (k-1)*npl+1
        write(stderr,myf) a(kk:kk+npl-1)
      end do
      kk = n-(nn*npl)+1
    else
      kk = n
    end if
    if ( kk > 0 ) then
      write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
      write(stderr,myf) a(n-kk+1:)
    end if
  end subroutine vprntv
!
  subroutine vprntm(a,n1,n2,nam)
    implicit none
    integer(ik4) , intent (in) :: n1 , n2
    character(len=*) , intent (in) :: nam
    real(rk8) , intent (in) , dimension(n1,n2) :: a
    integer(ik4) :: k1 , k2 , nn , kk
    integer(ik4) , parameter :: npl = 7
    character(len=32) :: myf
    write(stderr,'(2x,a,a,a)') '# ',nam,' #'
    nn = n1/npl
    do k2 = 1 , n2
      write(myf,'(a,i1,a)') '(2x,',npl,'g11.3)'
      write(stderr,*) '## Row ',k2
      if ( n1 > npl ) then
        do k1 = 1 , nn
          kk = (k1-1)*npl+1
          write(stderr,myf) a(kk:kk+npl-1,k2)
        end do
        kk = n1-(nn*npl)+1
      else
        kk = n1
      end if
      if ( kk > 0 ) then
        write(myf,'(a,i1,a)') '(2x,',kk,'g11.3)'
        write(stderr,myf) a(n1-kk+1:,k2)
      end if
      write(stderr,*) '## '
    end do
  end subroutine vprntm
!
end module mod_mpmessage
