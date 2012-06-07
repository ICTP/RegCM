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

  use mod_dynparam
  use mod_date
  use mod_grid
  use mod_wrtoxd
  use mod_header
  use mod_ch_icbc
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
  integer :: nnn , nsteps
  integer :: ierr
  character(256) :: namelistfile , prgname
  character(len=8)   :: chemsimtype
  integer :: ichremlsc , ichremcvc , ichdrdepo , ichcumtra , &
             ichsolver , idirect , ichdustemd
  logical :: dochem , dooxcl , doaero
  data dochem /.false./
  data dooxcl /.false./
  data doaero /.false./
!
  namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo ,  &
    ichcumtra , ichsolver , idirect , ichdustemd
!
  call header('chem_icbc')
!
! Read input global namelist
!
  call getarg(0, prgname)
  call getarg(1, namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write ( 6, * ) 'Parameter initialization not completed'
    write ( 6, * ) 'Usage : '
    write ( 6, * ) '          ', trim(prgname), ' regcm.in'
    write ( 6, * ) ' '
    write ( 6, * ) 'Check argument and namelist syntax'
    stop
  end if
!
  read(ipunit, chemparam, err=101)
!
  select case (chemsimtype)
    case ( 'CBMZ' )
      dochem = .true.
    case ( 'DUST', 'SSLT' )
      doaero = .true.
    case ( 'CARB' , 'SULF' , 'SUCA' , 'AERO' )
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

  idate = globidate1
  iodate = idate

  if ( dochem ) call newfile_ch_icbc(idate)
  if ( dooxcl ) call newfile_ox_icbc(idate)
  if ( doaero ) call newfile_ae_icbc(idate)

  if ( dochem ) call header_ch_icbc(idate)
  if ( dooxcl ) call header_ox_icbc
  if ( doaero ) call header_ae_icbc(idate)

  do nnn = 1 , nsteps
   if (.not. lsamemonth(idate, iodate) ) then
     if ( dochem ) call newfile_ch_icbc(monfirst(idate))
     if ( dooxcl ) call newfile_ox_icbc(monfirst(idate))
     if ( doaero ) call newfile_ae_icbc(monfirst(idate))
   end if
   if ( dochem ) call get_ch_icbc(idate)
   if ( dooxcl ) call get_ox_icbc(idate)
   if ( doaero ) call get_ae_icbc(idate)
   iodate = idate
   idate = idate + tbdy
  end do

  call close_outoxd
  if ( dochem ) call close_ch_icbc
  if ( dooxcl ) call close_ox_icbc
  if ( doaero ) call close_ae_icbc

  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed CHEM ICBC'
  stop

  101 write (stderr,*) 'Cannot read namelist stanza: chemparam'
  write (stderr,*) 'Assuming nothing to do for this experiment (ichem = 0)'
  call finaltime(0)
  write(stdout,*) 'Successfully completed CHEM ICBC'

end program chem_icbc
