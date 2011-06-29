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

program chem_icbc

  use mod_dynparam
  use mod_date
  use mod_grid
  use mod_wrtoxd
  use mod_header
  use mod_ch_icbc
  use mod_ch_oxcl

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
  call init_grid(iy,jx,kz)
  call init_outoxd
!
  tdif = globidate2-globidate1
  tbdy = rcm_time_interval(ibdyfrq,uhrs)
  nsteps = idnint(tdif%hours())/ibdyfrq + 1
!
  write (*,*) 'GLOBIDATE1 : ' , globidate1%tostring()
  write (*,*) 'GLOBIDATE2 : ' , globidate2%tostring()
  write (*,*) 'NSTEPS     : ' , nsteps

  idate = globidate1
  iodate = idate

  call newfile_ch_oxcl(idate)
  call newfile_ch_icbc(idate)

  call headermozart_ch_oxcl
  call headermozart_ch_icbc

  do nnn = 1 , nsteps
   if (.not. lsamemonth(idate, iodate) ) then
     call newfile_ch_icbc(idate)
     call newfile_ch_oxcl(idate)
   end if
   call getmozart_ch_icbc(idate)
   call getmozart_ch_oxcl(idate)
   iodate = idate
   idate = idate + tbdy
  end do

  call free_grid
  call free_outoxd
  call freemozart_ch_icbc
  call freemozart_ch_oxcl

end program chem_icbc
