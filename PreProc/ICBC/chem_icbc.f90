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
  use mod_ch_oxcl
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
  call memory_init
!  
  call init_grid(jx,iy,kz)
  call init_outoxd
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

  call newfile_ch_oxcl(idate)
  call newfile_ch_icbc(idate)

  call header_ch_oxcl
  call header_ch_icbc

  do nnn = 1 , nsteps
   if (.not. lsamemonth(idate, iodate) ) then
     call newfile_ch_icbc(monfirst(idate))
     call newfile_ch_oxcl(monfirst(idate))
   end if
   call get_ch_icbc(idate)
   call get_ch_oxcl(idate)
   iodate = idate
   idate = idate + tbdy
  end do

  call close_outoxd
  call close_ch_oxcl
  call close_ch_icbc

  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed CHEM ICBC'

end program chem_icbc
