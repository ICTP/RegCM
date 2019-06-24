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
!

#ifdef PNETCDF
subroutine myabort
  use mod_stdio
  use mpi
  implicit none
  integer :: ierr
  write(stderr,*) ' Execution terminated because of runtime error'
  call mpi_abort(mpi_comm_self,1,ierr)
end subroutine myabort
#else
subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort
#endif

program clmbc

#ifdef CLM45

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!  CLMBCC reads ERA5 surface hourly fields and creates hourly input    !
!   used by clmsa program in the RegCM to run the CLM4.5 component of  !
!   the RegCM as stand-alone program decoupled from the atmospheric    !
!   component.                                                         !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_header
  use mod_stdio
  use mod_memutil
  use mod_date
  use mod_grid
  use mod_date
  use mod_write
  use mod_era5
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  integer(ik4) :: nnn
  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdiff , tbdy
  integer(ik4) :: nsteps
  integer(ik4) :: ierr
  character(len=256) :: namelistfile , prgname

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  call header('clmbc')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr /= 0 ) then
    write ( stderr, * ) 'Parameter initialization not completed'
    write ( stderr, * ) 'Usage : '
    write ( stderr, * ) '          ', trim(prgname), ' regcm.in'
    write ( stderr, * ) ' '
    write ( stderr, * ) 'Check argument and namelist syntax'
    call die('clmbc','Check argument and namelist syntax',1)
  end if

  call memory_init

  call init_hgrid(jx,iy,kz)
  call init_houtput

  tdiff = globidate2-globidate1
  tbdy = rcm_time_interval(1,uhrs)
  nsteps = nint(tohours(tdiff)) + 1

  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (stdout,*) 'NSTEPS     : ' , nsteps

  idate = globidate1
  iodate = idate

  call init_era5h

  call newhfile(idate)

  do nnn = 1 , nsteps

    if (.not. lsamemonth(idate, iodate) ) then
      call newhfile(monfirst(idate))
    end if

    call get_era5h(idate)

    call writehf(idate)

    iodate = idate
    idate = idate + tbdy

  end do

  call close_output

  call conclude_era5h

  call dispose_output
  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed CLMBC'
#else
  write(0,*) 'This programs is enabled only if CLM45 is compiled in.'
#endif

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

end program clmbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
