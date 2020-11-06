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
  use mod_ch_fnest
  use mod_memutil
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdif , tbdy
  integer(ik4) :: nnn , nsteps
  integer(ik4) :: ierr , ipunit
  character(len=256) :: namelistfile , prgname
  character(len=256) :: cdir , cname
  character(len=8)   :: chemsimtype
  integer(ik4) :: ichremlsc , ichremcvc , ichdrdepo , ichcumtra , &
        ichsolver , idirect , ichdustemd , ichdiag , ichsursrc ,  &
        iindirect , ichebdy , ichjphcld , ichbion , ismoke , ichlinox , &
        isnowdark , ichecold
  logical :: carb_aging_control
  integer(ik4) :: ichem , iclimaaer
  integer(ik4) ibltyp , iboudy , isladvec , iqmsl , icup_lnd , icup_ocn , &
    ipgf , iemiss , lakemod , ipptls , iocnflx , iocncpl , iwavcpl ,      &
    iocnrough , iocnzoq , idcsst , iseaice , idesseas , iconvlwp ,        &
    icldmstrat , icldfrac , irrtm , iclimao3 , isolconst , icumcloud ,    &
    islab_ocean , itweak , ghg_year_const , idiffu , icopcpl ,            &
    iwhitecap , ifixsolar , year_offset , ipcpcool , ichdustparam
  real(rkx) :: temp_tend_maxval , wind_tend_maxval , fixedsolarval
  character(len=8) :: scenario
  real(rkx) :: rdstemfac , rocemfac
  logical :: dochem , dooxcl , doaero
  data dochem /.false./
  data dooxcl /.false./
  data doaero /.false./
  data ichem /0/
  data iclimaaer /0/

  namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo , &
    ichcumtra , ichsolver , idirect , ichdustemd , ichdiag , iindirect , &
    ichsursrc , ichebdy , rdstemfac , rocemfac , ichjphcld , ichbion ,   &
    ismoke , ichlinox , isnowdark , carb_aging_control , ichecold ,      &
    ichdustparam

  namelist /physicsparam/ ibltyp , iboudy , isladvec , iqmsl ,         &
    icup_lnd , icup_ocn , ipgf , iemiss , lakemod , ipptls , idiffu ,  &
    iocnflx , iocncpl , iwavcpl , icopcpl , iocnrough , iocnzoq ,      &
    ichem ,  scenario ,  idcsst , ipcpcool , iwhitecap , iseaice ,     &
    idesseas , iconvlwp , icldmstrat , icldfrac , irrtm , iclimao3 ,   &
    iclimaaer , isolconst , icumcloud , islab_ocean , itweak ,         &
    temp_tend_maxval , wind_tend_maxval , ghg_year_const , ifixsolar , &
    fixedsolarval , year_offset

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

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

  ! Read also chemparam
  open(newunit=ipunit, file=namelistfile, status='old', &
       action='read', iostat=ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Cannot open ',trim(namelistfile)
    stop
  end if
  read(ipunit, physicsparam, iostat=ierr)
  if ( ierr /= 0 ) then
    write (stderr,*) 'Error reading physicsparam namelist'
    write (stderr,*) 'Check namelist content and syntax'
    stop
  end if
  if ( ichem == 1 .or. iclimaaer == 1 ) then
    rewind(ipunit)
    read(ipunit, chemparam, iostat=ierr)
    if ( ierr /= 0 ) then
      write (stderr,*) 'Cannot read namelist stanza: chemparam'
      write (stderr,*) 'Assuming nothing to do for this experiment'
      stop
    end if
    close(ipunit)
  else
    close(ipunit)
    write (stderr,*) 'Assuming nothing to do for this experiment'
    write (stderr,*) 'Both ichem and iclimaaer equal to zero.'
    stop
  end if

  if ( chemtyp .eq. 'FNEST' ) then
    call init_fnestparam(namelistfile,cdir,cname)
  end if

  if ( iclimaaer == 1 ) then
    doaero = .true.
    chemsimtype = 'AERO'
  else
    select case (chemsimtype)
      case ( 'CBMZ' )
        dochem = .true.
      case ( 'DUST', 'DU12', 'SSLT', 'DUSS' )
        doaero = .true.
      case ( 'CARB' , 'SULF' , 'SUCA' , 'AERO' , 'SUCE' )
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
  end if

  call memory_init

  call init_grid(jx,iy,kz)
  call init_outoxd(chemsimtype)

  tdif = globidate2-globidate1
  tbdy = rcm_time_interval(ibdyfrq,uhrs)
  nsteps = nint(tohours(tdif))/ibdyfrq + 1

  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (stdout,*) 'NSTEPS     : ' , nsteps
  if ( dochem ) then
    write (stdout,*) 'chemtyp     : ' , chemtyp
  else
    write (stdout,*) 'Aerosol simulation.'
  end if

  idate = globidate1
  iodate = idate

  if ( dochem ) call newfile_ch_icbc(idate)
  if ( dooxcl ) call newfile_ox_icbc(idate)
  if ( doaero ) call newfile_ae_icbc(idate)

  if ( chemtyp .eq. 'FNEST' ) then
    call init_fnest(idate,cdir,cname,dochem,dooxcl,doaero)
  else
    if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call init_ch_icbc(idate)
    if ( dochem .and. chemtyp .eq. 'MZCLM' ) call init_ch_icbc_clim(idate)
    if ( dooxcl ) call init_ox_icbc(idate)
    if ( doaero ) call init_ae_icbc(idate)
  end if

  do nnn = 1 , nsteps
   if (.not. lsamemonth(idate, iodate) ) then
     if ( dochem ) call newfile_ch_icbc(monfirst(idate))
     if ( dooxcl ) call newfile_ox_icbc(monfirst(idate))
     if ( doaero ) call newfile_ae_icbc(monfirst(idate))
   end if
   if ( chemtyp .eq. 'FNEST' ) then
     call get_fnest(idate)
   else
     if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call get_ch_icbc(idate)
     if ( dochem .and. chemtyp .eq. 'MZCLM' ) call get_ch_icbc_clim(idate)
     if ( dooxcl ) call get_ox_icbc(idate)
     if ( doaero ) call get_ae_icbc(idate)
   end if
   iodate = idate
   idate = idate + tbdy
  end do

  call close_outoxd
  if ( chemtyp .eq. 'FNEST' ) then
    call close_fnest
  else
    if ( dochem .and. chemtyp .eq. 'MZ6HR' ) call close_ch_icbc
    if ( dochem .and. chemtyp .eq. 'MZCLM' ) call close_ch_icbc_clim
    if ( dooxcl ) call close_ox_icbc
    if ( doaero ) call close_ae_icbc
  end if

  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed CHEM ICBC'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

end program chem_icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
