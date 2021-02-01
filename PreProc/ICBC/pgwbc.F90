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

program pgwbc
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
  use mod_projections
  use mod_pgw
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  integer(ik4) :: nnn
  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdiff , tbdy
  integer(ik4) :: nsteps
  integer(ik4) :: ierr
  character(len=256) :: namelistfile, prgname, infilename
  type(anyprojparams) :: pjpara

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  call header('pgwbc')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call get_command_argument(2,value=infilename)
  call initparam(namelistfile, ierr)
  if ( idynamic == 2 ) then
    write(stdout, *) 'Using non hydrostatic parameters'
    write(stdout, '(a,f10.2)') ' base_state_pressure    = ', base_state_pressure
    write(stdout, '(a,f10.2)') ' logp_lrate             = ', logp_lrate
  end if

  if ( idynamic == 3 ) then
    write(stdout, *) 'Using Moloch non-hydrostatic dynamical core'
  end if

  if ( ierr /= 0 ) then
    write ( stderr, * ) 'Parameter initialization not completed'
    write ( stderr, * ) 'Usage : '
    write ( stderr, * ) '          ', trim(prgname), ' regcm.in input.nc'
    write ( stderr, * ) ' '
    write ( stderr, * ) 'Check argument and namelist syntax'
    call die('pgwbc','Check argument and namelist syntax',1)
  end if
!
  call memory_init

  call init_grid(jx,iy,kz)
  call init_output

  pjpara%pcode = iproj
  pjpara%ds = ds*1000.0_rk8
  pjpara%clat = clat
  pjpara%clon = clon
  pjpara%plat = plat
  pjpara%plon = plon
  pjpara%trlat1 = truelatl
  pjpara%trlat2 = truelath
  pjpara%nlon = jx
  pjpara%nlat = iy
  pjpara%rotparam = .true.

  if ( idynamic == 3 ) then
    pjpara%staggerx = .true.
    pjpara%staggery = .false.
    call pju%initialize(pjpara)
    pjpara%staggerx = .false.
    pjpara%staggery = .true.
    call pjv%initialize(pjpara)
  else
    pjpara%staggerx = .true.
    pjpara%staggery = .true.
    call pjd%initialize(pjpara)
  end if

  tdiff = globidate2-globidate1
  tbdy = rcm_time_interval(ibdyfrq,uhrs)
  nsteps = nint(tohours(tdiff))/ibdyfrq + 1

  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (stdout,*) 'NSTEPS     : ' , nsteps

  idate = globidate1
  iodate = idate

  call init_pgw(infilename)
  call newfile(idate)

  do nnn = 1 , nsteps
    if (.not. lsamemonth(idate, iodate) ) then
      call newfile(monfirst(idate))
    end if
    call get_pgw(idate)
    call writef(idate)
    iodate = idate
    idate = idate + tbdy
  end do

  call close_output
  call conclude_pgw

  call pju%destruct( )
  call pjv%destruct( )
  call pjd%destruct( )
  call dispose_output
  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed PGWBC'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

end program pgwbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
