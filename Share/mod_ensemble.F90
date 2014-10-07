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

module mod_ensemble
!
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_constants

  implicit none
!
!------------------------------------------------------------------------------
!
! Routines to implement the ensembling method of
! O'Brien, Sloan, and Snyder (2010), Climate Dynamics
!    (DOI: 10.1007/s00382-010-0900-5)
!
!------------------------------------------------------------------------------
!
  private

  interface randify
    module procedure randify3D
    module procedure randify2D
  end interface randify

  public :: randify

  integer(ik4) , dimension(:) , pointer :: seed

  contains
!
! Takes a 3D variable (dVariable3D) and adds random noise to all of its
! values.  At point (i,j,k) in dVariable3D, the random noise will be
! within the half-open interval:
!   [-dFrac*dVariable3D(i,j,k) , +dFrac*dVariable3D(i,j,k))
!
! dVariable3D  -  A 3D floating point array
!
! dFrac        -  Maximum fraction by which to vary any value in dVariable3D
!
  subroutine randify3D(dVariable3D,dFrac,imax,jmax,kmax)
    implicit none
    integer(ik4) , intent(in) :: imax , jmax , kmax
    real(rk8) , dimension(imax,kmax,jmax) , intent(inout) :: dVariable3D
    real(rk8) , intent(in) :: dFrac

    real(rk8) , dimension(imax,kmax,jmax) :: dChange3D , dRand3D
    integer(ik4) :: i
    integer(ik4) :: nseed
    real(rk4) :: cputime

    ! initialize the random number generator with the current clock time

    if ( .not. associated(seed) ) then

      ! get the size of the seed array

      call random_seed(size = nseed)

      ! allocate a new seed array

      call getmem1d(seed,1,nseed,'randify2D:seed')

      ! Get the system time

      call cpu_time(cputime)

      ! TAO:  The odd syntax for this line comes from GNU documentation. I don't
      ! understand why 37 is used as opposed to any other number.

      seed = int(cputime) + 37*(/(i-1,i=1,nseed)/)

      ! Set the seed for the random number generator.  This makes it so that we
      ! get a pseudo-random sequence of numbers

      call random_seed(put = seed)

    end if

    ! Figure out how much to tweak the variable by (at most)

    dChange3D = dVariable3D*dFrac

    ! Generate a random number within +/- this range

    call random_number(dRand3D)
    dRand3D = d_two*dChange3D*(dRand3D - d_half)

    ! Add this random number to the variable

    dVariable3D = dVariable3D + real(dRand3D)

  end subroutine randify3D

! Takes a 2D variable (dVariable2D) and adds random noise to all of its
! values.  At point (i,j,k) in dVariable2D, the random noise will be
! within the half-open interval:
!   [-dFrac*dVariable2D(i,j,k) , +dFrac*dVariable2D(i,j,k))
!
! dVariable2D  -  A 2D floating point array
! dFrac        -  Maximum fraction by which to vary any value in dVariable2D
!
  subroutine randify2D(dVariable2D,dFrac,imax,jmax)
    implicit none
    integer(ik4) , intent(in) :: imax , jmax
    real(rk8) , dimension(imax,jmax) , intent(inout) :: dVariable2D
    real(rk8) , intent(in) :: dFrac

    real(8) , dimension(imax,jmax) :: dRand2D , dChange2D
    integer(ik4) :: i
    integer(ik4) :: nseed
    real(rk4) :: cputime

    ! initialize the random number generator with the current clock time

    if ( .not. associated(seed) ) then

      ! get the size of the seed array

      call random_seed(size = nseed)

      ! allocate a new seed array

      call getmem1d(seed,1,nseed,'randify2D:seed')

      ! Get the system time

      call cpu_time(cputime)

      ! TAO:  The odd syntax for this line comes from GNU documentation. I don't
      ! understand why 37 is used as opposed to any other number.

      seed = int(cputime) + 37*(/(i-1,i=1,nseed)/)

      ! Set the seed for the random number generator.  This makes it so that we
      ! get a pseudo-random sequence of numbers

      call random_seed(put = seed)

    end if

    ! Figure out how much to tweak the variable by (at most)

    dChange2D = dVariable2D*dFrac

    ! Generate a random number within +/- this range

    call random_number(dRand2D)
    dRand2D = d_two*dChange2D*(dRand2D - d_half)

    ! Add this random number to the variable

    dVariable2D = dVariable2D + real(dRand2D)

  end subroutine randify2D

end module mod_ensemble
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
