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

module mod_gridfuncs

  use mod_runparams
  use mod_main , only : atmstate
  use mod_constants , only : d_rfour
  use mpi

  private

  public :: uvcross2dot

  contains
!
! Takes an atmstate variable with u and v on the cross grid (the
! same grid as t, qv, qc, etc.) and interpolates the u and v to
! the dot grid.  This routine sheilds the user of the function
! from the need to worry about the details of the domain
! decomposition.  
!
! Written by Travis A. O'Brien 01/04/11.
!
! type(atmstate),intent(in) :: invar 
!                              An atmstate variable (see mod_main)
!                              that contains the u and v variables
!                              that need to be interpolated.  This
!                              variable is not modified by this
!                              routine.
!
! type(atmstate),intent(inout) :: outvar 
!                              An atmstate variable (see mod_main)
!                              that contains the u and v variables
!                              that will be overwritten by the
!                              interpolation of invar%u and invar%v.
!                              Only u and v are modified in this
!                              routine (t, qv, qc, and tke should
!                              remain unchanged).
!
  subroutine uvcross2dot(invar,outvar)
    implicit none
    type(atmstate) , intent(inout) :: invar
    type(atmstate) , intent(inout) :: outvar
    integer :: ib , ie , jb , je , i , j
    integer :: isendcount , ierr

    ! TODO:  It might make sense to encapsulate the following code
    ! in to a standard routine, since this boundary sending code is
    ! ubiquitous throughout the RegCM code and it is domain
    ! decomposition-dependent.

    ! Send the right-edge of the u/v tendencies to the left
    ! edge of the next process's u/v tendencies (so that
    ! invar%u(i,k,0) holds invar%u(i,k,jxp) of the parallel
    ! chunk next door)

    isendcount = iy*kz
    call mpi_sendrecv(invar%u(:,:,jxp),isendcount,mpi_real8,ieast,30, &
                      invar%u(:,:,0),isendcount,mpi_real8,iwest,30,   &
                      mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_sendrecv(invar%v(:,:,jxp),isendcount,mpi_real8,ieast,31, &
                      invar%v(:,:,0),isendcount,mpi_real8,iwest,31,   &
                      mpi_comm_world,mpi_status_ignore,ierr)

    ! Set j-loop boundaries
    jb = jbegin
    je = jendx
    ! Set i-loop boundaries
    ib = 2
    ie = iym1

    !
    !     x     x     x     x     x     x
    !
    !        o     o     o     o     o 
    !         (i-1,j-1)     (i,j-1)            
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !        o     o  |  o  |  o     o 
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !
    !        o     o     o     o     o 
    !
    !     x     x     x     x     x     x
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do j = jb , je
      do i = ib , ie
        outvar%u(i,:,j) =  outvar%u(i,:,j) +             &
          d_rfour*(invar%u(i,:,j) + invar%u(i,:,j-1) +   &
                   invar%u(i-1,:,j) + invar%u(i-1,:,j-1))
        outvar%v(i,:,j) =  outvar%v(i,:,j) +             &
          d_rfour*(invar%v(i,:,j) + invar%v(i,:,j-1) +   &
                   invar%v(i-1,:,j) + invar%v(i-1,:,j-1))
      end do
    end do
  end subroutine uvcross2dot

end module mod_gridfuncs
