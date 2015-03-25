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

subroutine myabort
  implicit none
  call abort
end subroutine myabort

program sst

  use mod_intkinds
  use mod_realkinds
  use mod_header
  use mod_dynparam
  use mod_memutil
  use mod_stdio
  use mod_message
  use mod_sst_grid
  use mod_sst_1deg
  use mod_sst_eh5om
  use mod_sst_ersst
  use mod_sst_fvgcm
  use mod_sst_gnmnc
  use mod_sst_gnhnc

  implicit none

  integer(ik4) :: ierr
  character(len=256) :: namelistfile , prgname
  character(len=256) :: terfile

!     call and print header
  call header('sst')
!
!     Read input global namelist
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

  if ( ssttyp == 'GISST' .or. ssttyp == 'OISST' .or.       &
       ssttyp == 'OI_NC' .or. ssttyp == 'OI2ST' .or.       &
       ssttyp == 'OI_WK' .or. ssttyp == 'OI2WK' ) then
    call sst_1deg
  else if ( ssttyp(1:2) == 'EH' ) then
    call sst_eh5om
  else if ( ssttyp == 'ERSST' .or. ssttyp == 'ERSKT' .or. &
            ssttyp(1:3) == 'EIN' ) then
    call sst_ersst
  else if ( ssttyp == 'FV_A2' .or.  ssttyp == 'FV_B2' ) then
    call sst_fvgcm
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
  else if ( ssttyp(1:3) == 'MP_' .or. ssttyp == 'E5_A2' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnhnc
  else if ( ssttyp == 'EIXXX' .or. ssttyp == 'CCSM3' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
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

  if (debug_level > 2) then
  end if

  call finaltime(0)
  write (stdout,*) 'Successfully generated SST'

end program sst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
