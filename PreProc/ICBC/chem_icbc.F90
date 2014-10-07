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

program chem_icbc

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date
  use mod_grid
  use mod_wrtoxd
  use mod_header
  use mod_ch_icbc
  use mod_ch_icbc_clim
  use mod_ox_icbc
  use mod_ae_icbc
  use mod_memutil
!
  implicit none
!
!  Local Variables
!
  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdif , tbdy
  integer(ik4) :: nnn , nsteps
  integer(ik4) :: ierr
  character(len=256) :: namelistfile , prgname
  character(len=8)   :: chemsimtype
  integer(ik4) :: ichremlsc , ichremcvc , ichdrdepo , ichcumtra , &
        ichsolver , idirect , ichdustemd , ichdiag , ichsursrc ,  &
        iindirect , ichebdy , ichjphcld , ichbion
  real(rk8) :: rdstemfac
  logical :: dochem , dooxcl , doaero
  data dochem /.false./
  data dooxcl /.false./
  data doaero /.false./
!
  namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo , &
     ichcumtra , ichsolver , idirect , ichdustemd , ichdiag , iindirect ,&
     ichsursrc , ichebdy , rdstemfac , ichjphcld , ichbion

  call header('chem_icbc')
!
! Read input global namelist
!
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Usage : '
    write (stderr,*) '          ', trim(prgname), ' regcm.in'
    write (stderr,*) ' '
    write (stderr,*) 'Check argument and namelist syntax'
    stop
  end if
!
  open(ipunit, file=namelistfile, status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Cannot open ',trim(namelistfile)
    stop
  end if
  read(ipunit, chemparam, iostat=ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Cannot read namelist stanza: chemparam'
    write (stderr,*) 'Assuming nothing to do for this experiment (ichem = 0)'
    stop
  end if
  close(ipunit)
!
  select case (chemsimtype)
    case ( 'CBMZ' )
      dochem = .true.
    case ( 'DUST', 'SSLT', 'DUSS' )
      doaero = .true.
    case ( 'CARB' , 'SULF' , 'SUCA' , 'AERO' )
      doaero = .true.
      dooxcl = .true.
    case ( 'DCCB' )
      dochem = .true.
      doaero = .true.
      dooxcl = .true.
    case default
      write (stderr,*) 'Unknown chemsimtype'
      write (stderr,*) 'Assuming nothing to do for this experiment'
      call finaltime(0)
      write(stdout,*) 'Successfully completed CHEM ICBC'
      stop
  end select

  call memory_init
!
  call init_grid(jx,iy,kz)
  call init_outoxd(chemsimtype)
!
  tdif = globidate2-globidate1
  tbdy = rcm_time_interval(ibdyfrq,uhrs)
  nsteps = idnint(tohours(tdif))/ibdyfrq + 1
!
  write (*,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (*,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (*,*) 'NSTEPS     : ' , nsteps
  write (*,*) 'chemtyp     : ' , chemtyp

  idate = globidate1
  iodate = idate

  if ( dochem ) call newfile_ch_icbc(idate)
  if ( dooxcl ) call newfile_ox_icbc(idate)
  if ( doaero ) call newfile_ae_icbc(idate)

  if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call header_ch_icbc(idate)
  if ( dochem .and. chemtyp .eq. 'MZCLM' ) call header_ch_icbc_clim(idate)
  if ( dooxcl ) call header_ox_icbc
  if ( doaero ) call header_ae_icbc(idate)

  do nnn = 1 , nsteps
   if (.not. lsamemonth(idate, iodate) ) then
     if ( dochem ) call newfile_ch_icbc(monfirst(idate))
     if ( dooxcl ) call newfile_ox_icbc(monfirst(idate))
     if ( doaero ) call newfile_ae_icbc(monfirst(idate))
   end if
   if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call get_ch_icbc(idate)
   if ( dochem .and. chemtyp .eq. 'MZCLM' ) call get_ch_icbc_clim(idate)
   if ( dooxcl ) call get_ox_icbc(idate)
   if ( doaero ) call get_ae_icbc(idate)
   iodate = idate
   idate = idate + tbdy
  end do

  call close_outoxd
  if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call close_ch_icbc
  if ( dochem .and. chemtyp .eq. 'MZCLM' ) call close_ch_icbc_clim
  if ( dooxcl ) call close_ox_icbc
  if ( doaero ) call close_ae_icbc

  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed CHEM ICBC'
end program chem_icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
