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

module mod_template

  use mod_basicmodule

  private

  real(dp) :: arealvar
  integer :: integervar
  real(dp) , dimension(4) :: fixed_dims_array
  real(sp) , allocatable , dimension(:,:) :: allocatable_2d_array

  public :: allocatable_2d_array    ! Data accessed from outside
  public :: init_mod_basicmodule    ! Allocate space and init data
  public :: release_mod_basicmodule ! Release allocated resources and cleanup
  public :: meaningful_subroutine   ! A subroutine doing some work

  contains

  subroutine init_mod_basicmodule(lfag)
    implicit none

    ! Dummy arguments are to be specified of theyr intent

    logical , intent(in) :: lflag   ! Flag to control behaviour

    ! Local variables

    integer :: istat

    !-----------------------------------------------------------

    if ( lflag ) then
      allocate(allocatable_2d_array(iy,jx), stat=istat)
    end if

  end subroutine init_mod_basicmodule

  subroutine release_mod_basicmodule
    implicit none

    !-----------------------------------------------------------

    if ( allocated(allocatable_2d_array) ) then
      deallocate(allocatable_2d_array)
    end if

  end subroutine release_mod_basicmodule

  real(dp) function internal_function(a)
    implicit none
    real , intent(in) :: a

    !-----------------------------------------------------------

    ! Constants are to be declared in double precision
    internal_function = dlog(a) * 3.0D0 + 4.0D0

  end function internal_function

  subroutine meaningful_subroutine
    implicit none

    integer :: i , j

    !-----------------------------------------------------------

    do i = 1 , iy
      do j = 1 , jx
        ! Cast of variables need to be specified
        allocatable_2d_array(i,j) = internal_function(dble(i))
      end do
    end do

  end subroutine meaningful_subroutine

end module mod_template
