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

  integer :: ierr
  character(256) :: namelistfile , prgname
  character(256) :: terfile

!     call and print header
  call header('sst')
!
!     Read input global namelist
!
  call getarg(0, prgname)
  call getarg(1, namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Usage : '
    write (stderr,*) '          ', trim(prgname), ' regcm.in'
    write (stderr,*) ' '
    call die('sst','Check argument and namelist syntax.',1)
  end if

  if (debug_level > 2) then
  end if

  call memory_init

  call init_grid
  terfile=trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
  call read_domain_info(terfile)

  if ( ssttyp == 'GISST' .or. ssttyp == 'OISST' .or.       &
       ssttyp == 'OI_NC' .or. ssttyp == 'OI2ST' .or.       &
       ssttyp == 'OI_WK' .or. ssttyp == 'OI2WK' ) then
    call sst_1deg
  else if ( ssttyp == 'EH5RF' .or. ssttyp == 'EH5A2' .or.  &
            ssttyp == 'EH5B1' .or. ssttyp == 'EHA1B' ) then
    call sst_eh5om
  else if ( ssttyp == 'ERSST' .or. ssttyp == 'ERSKT' ) then
    call sst_ersst
  else if ( ssttyp == 'FV_RF' .or. ssttyp == 'FV_A2' .or.  &
            ssttyp == 'FV_B2' ) then
    call sst_fvgcm
  else if ( ssttyp == 'CCSST' .or. ssttyp == 'CAM4N' .or. &
            ssttyp == 'CA_RF' .or. ssttyp == 'CA_26' .or. &
            ssttyp == 'CA_45' .or. ssttyp == 'CA_85') then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'HA_RF' .or. ssttyp == 'HA_26' .or. &
            ssttyp == 'HA_45' .or. ssttyp == 'HA_85' ) then
    if (ical /= y360 ) then
      write(stderr,*) ssttyp//' calendar should be set to 360_day'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'CS_RF' .or. ssttyp == 'CS_26' .or. &
            ssttyp == 'CS_45' .or. ssttyp == 'CS_85' ) then
    if (ical /= noleap ) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'EC_RF' .or. ssttyp == 'EC_45' .or. &
            ssttyp == 'EC_85' ) then
    call sst_gnmnc
  else if ( ssttyp == 'IP_RF' .or. ssttyp == 'IP_45' .or. &
            ssttyp == 'IP_85' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'GF_RF' .or. ssttyp == 'GF_45' .or. &
            ssttyp == 'GF_85' ) then
    if (ical /= noleap) then
      write(stderr,*) ssttyp//' calendar should be set to noleap'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'CN_RF' .or. ssttyp == 'CN_45' .or. &
            ssttyp == 'CN_85' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
      call die('sst','Calendar mismatch',1)
    end if
    call sst_gnmnc
  else if ( ssttyp == 'MP_RF' .or. ssttyp == 'MP_45' .or. &
            ssttyp == 'MP_85' ) then
    if (ical /= gregorian) then
      write(stderr,*) ssttyp//' calendar should be set to gregorian'
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
