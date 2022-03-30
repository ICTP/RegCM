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

program sst

  use mod_intkinds
  use mod_realkinds
  use mod_header
  use mod_dynparam
  use mod_memutil
  use mod_stdio
  use mod_date
  use mod_message
  use mod_sst_grid
  use mod_sst_1deg
  use mod_sst_ersst
  use mod_sst_gnmnc
  use mod_sst_gndnc
  use mod_sst_gnhnc
  use mod_sst_cmip6
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  integer(ik4) :: ierr
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  ! call and print header
  call header('sst')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Usage : '
    write (stderr,*) '          ', trim(prgname), ' regcm.in'
    write (stderr,*) ' '
    call die('sst','Check argument and namelist syntax.',1)
  end if

  call memory_init

  call init_grid
  terfile=trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
  call read_domain_info(terfile)
  call setup_outvars

  if ( ssttyp == 'CMIP6' ) then
    call cmip6_sst
  else if ( ssttyp == 'GISST' .or. ssttyp == 'OISST' .or.  &
       ssttyp == 'OI_NC' .or. ssttyp == 'OI2ST' .or.       &
       ssttyp == 'OI_WK' .or. ssttyp == 'OI2WK' ) then
    call sst_1deg
  else if ( ssttyp == 'TMIST' ) then
    call sst_gndnc
  else if ( ssttyp == 'ERSST' .or. ssttyp == 'ERSKT' ) then
    call sst_ersst
  else if ( ssttyp == 'CCSST' .or. ssttyp == 'CAM4N' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'CA_' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'HA_' ) then
    if (ical /= y360 ) then
      write(stderr,*) ssttyp//' calendar should be set to 360_day'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'CS_' ) then
    if (ical /= noleap ) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'MI_' ) then
    if (ical /= noleap ) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'EC_' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'IP_' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'GF_' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'CN_' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'CC_' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:2) == 'MP' .or. &
            ssttyp(1:3) == 'ECC' .or. &
            ssttyp == 'E5_A2' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else if ( ssttyp(1:3) == 'NO_' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gndnc
  else if ( ssttyp == 'EIXXX' .or. ssttyp == 'CCSM3' ) then
    if ( ical /= noleap ) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else if ( ssttyp == 'ERA5D' .or. ssttyp(1:3) == 'EID' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gndnc
  else if ( ssttyp(1:3) == 'EIN' .or. ssttyp(1:4) == 'ERA5' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else if ( ssttyp(1:3) == 'LGM' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else if ( ssttyp == 'JRA55' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp(1:3) == 'CFS' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else
    call die('sst', 'Unknown SSTTYP '//ssttyp//' specified in '// &
              trim(namelistfile)//'.',1)
  end if

  call close_sstfile

  call memory_destroy

  call finaltime(0)
  write (stdout,*) 'Successfully generated SST'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

end program sst

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
