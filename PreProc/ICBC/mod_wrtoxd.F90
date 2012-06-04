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

module mod_wrtoxd

  use netcdf
  use mod_dynparam
  use mod_date
  use mod_grid
  use mod_realkinds
  use mod_memutil
  use mod_message
  use mod_stdio
  use mod_nchelper

  private

  public :: chv4 , oxv4 , aev4
  public :: nchsp , noxsp , naesp
  public :: chspec , oxspec , aespec

  public :: init_outoxd , close_outoxd
  public :: newfile_ch_icbc , newfile_ox_icbc , newfile_ae_icbc
  public :: write_ch_icbc , write_ox_icbc , write_ae_icbc

  integer :: ncoutch , ncoutox , ncoutae
  character(256) :: ofname
  type(rcm_time_and_date) , save :: irefdate
  integer :: itimech , itimeox , itimeae
  integer , dimension(4) :: idims
  integer :: istatus

  integer , parameter :: nchsp = 25
  integer , parameter :: noxsp = 5
  integer :: naesp

  integer , dimension(nchsp+1) :: ichvar
  integer , dimension(noxsp+1) :: ioxvar
  integer , pointer , dimension(:) :: iaevar

  character(len=8) , dimension(nchsp) :: chspec
  character(len=8) , dimension(noxsp) :: oxspec

  character(len=8) , pointer , dimension(:) :: aespec

  character(len=8) , target , dimension(4) :: aedust
  character(len=8) , target , dimension(4) :: aesslt
  character(len=8) , target , dimension(5) :: aecarb
  character(len=8) , target , dimension(1) :: aesulf
  character(len=8) , target , dimension(6) :: aesuca
  character(len=8) , target , dimension(14) :: aeaero

  real(sp) , pointer , dimension(:,:,:,:) :: chv4
  real(sp) , pointer , dimension(:,:,:,:) :: oxv4
  real(sp) , pointer , dimension(:,:,:,:) :: aev4
  
  character(len=128) :: buffer

  data oxspec / 'OH' , 'HO2' , 'O3' , 'NO3' , 'H2O2' /
  data chspec / 'O3' , 'NO' , 'NO2' , 'HNO3' , 'N2O5' , 'H2O2' , 'CH4' , &
                'CO' , 'CH2O' , 'CH3OH' , 'C2H5OH' , 'C2H4' , 'C2H6' ,   &
                'CH3CHO' , 'CH3COCH3' , 'BIGENE' , 'BIGALK' , 'C3H6' ,   &
                'C3H8' , 'ISOP' , 'TOLUENE' , 'PAN' , 'SO2' , 'SO4' , 'DMS' /
  data aedust / 'DST01', 'DST02', 'DST03', 'DST04' /
  data aesslt / 'SSLT01' , 'SSLT02', 'SSLT03', 'SSLT04' /
  data aecarb / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' /
  data aesulf / 'SO4' /
  data aesuca / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' , 'SO4' /
  data aeaero / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' , 'SO4' , &
                'SSLT01' , 'SSLT02', 'SSLT03', 'SSLT04' ,       &
                'DST01', 'DST02', 'DST03', 'DST04' /

  integer :: ioc2 , isoa
  integer :: isslt1 , isslt2 , isslt3 , isslt4
  data ncoutch /-1/
  data ncoutox /-1/
  data ncoutae /-1/

  logical :: sum_soa_to_oc2
  logical :: sum_sslt_bins

  data sum_soa_to_oc2 /.false./
  data sum_sslt_bins  /.false./

  contains

  subroutine init_outoxd(chemsimtype)
    implicit none
    character(len=8) , intent(in) :: chemsimtype
    logical :: doaero , dochem , dooxcl
    integer :: i
    data doaero /.false./
    data dochem /.false./
    data dooxcl /.false./
    select case ( chemsimtype ) 
      case ( 'DUST' )
        naesp = 4
        aespec => aedust
        doaero = .true.
      case ( 'SSLT' )
        naesp = 4
        aespec => aesslt
        doaero = .true.
        sum_sslt_bins = .true.
      case ( 'CARB' )
        naesp = 5
        aespec => aecarb
        doaero = .true.
        sum_soa_to_oc2 = .true.
      case ( 'SULF' )
        naesp = 1
        aespec => aesulf
        doaero = .true.
      case ( 'SUCA' )
        naesp = 6
        aespec => aesuca
        doaero = .true.
        sum_soa_to_oc2 = .true.
      case ( 'AERO' )
        naesp = 14
        aespec => aeaero
        doaero = .true.
        sum_sslt_bins = .true.
        sum_soa_to_oc2 = .true.
      case ( 'CBMZ' )
        dochem = .true.
      case default
        call die('init_outoxd','Unknown chemsimtype')
    end select
    if ( doaero ) then
      call getmem4d(aev4,1,jx,1,iy,1,kz,1,naesp,'mod_wrtoxd:aev4')
      call getmem1d(iaevar,1,naesp+1,'mod_wrtoxd:iaevar')
    end if
    if ( dochem ) call getmem4d(chv4,1,jx,1,iy,1,kz,1,nchsp,'mod_wrtoxd:chv4')
    if ( dooxcl ) call getmem4d(oxv4,1,jx,1,iy,1,kz,1,noxsp,'mod_wrtoxd:oxv4')
    if ( sum_soa_to_oc2 ) then
      ioc2 = -1
      isoa = -1
      do i = 1 , naesp
        if ( aespec(i) == 'SOA' ) isoa = i
        if ( aespec(i) == 'OC2' ) ioc2 = i
      end do
      if ( isoa < 0 .or. ioc2 < 0 ) then
        call fatal(__FILE__,__LINE__,'Logical error: Search SOA error')
      end if
    end if
    if ( sum_sslt_bins ) then
      isslt1 = -1
      isslt2 = -1
      isslt3 = -1
      isslt4 = -1
      do i = 1 , naesp
        if ( aespec(i) == 'SSLT01' ) isslt1 = i
        if ( aespec(i) == 'SSLT02' ) isslt2 = i
        if ( aespec(i) == 'SSLT03' ) isslt3 = i
        if ( aespec(i) == 'SSLT04' ) isslt4 = i
      end do
      if ( isslt1 < 0 .or. isslt2 < 0 .or. isslt3 < 0 .or. isslt4 < 0 ) then
        call fatal(__FILE__,__LINE__,'Logical error: Search SSLT error.')
      end if
    end if
  end subroutine init_outoxd

  subroutine close_outoxd
    implicit none
    if (ncoutch > 0) then
      istatus = nf90_close(ncoutch)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close CH file')
    end if
    if (ncoutox > 0) then
      istatus = nf90_close(ncoutox)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close OX file')
    end if
    if (ncoutae > 0) then
      istatus = nf90_close(ncoutae)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close AE file')
    end if
  end subroutine close_outoxd

  subroutine newfile_ch_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: i , ipnt , istatus
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ihvar
    integer , dimension(2) :: illvar
    real(sp) , pointer , dimension(:) :: xjx , yiy
    character(64) :: csdate
    real(sp) :: hptop

    if (ncoutch > 0) then
      istatus = nf90_close(ncoutch)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close CH file')
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_CHBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itimech = 1

    csdate = tochar(idate1)

    call createfile_withname(ofname,ncoutch)
    call add_common_global_params(ncoutch,'chem_icbc')

    ipnt = 1
    call define_basic_dimensions(ncoutch,jx,iy,kz,ipnt,idims)
    call add_dimension(ncoutch,'time',nf90_unlimited,ipnt,idims)
!
    call define_horizontal_coord(ncoutch,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncoutch,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncoutch,idims,ipnt,illvar)

    istatus = nf90_def_var(ncoutch, 'time', nf90_double, idims(4:4), ichvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncoutch, ichvar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncoutch, ichvar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , nchsp
      istatus = nf90_def_var(ncoutch, chspec(i), nf90_float, idims, ichvar(i+1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding variable '//trim(chspec(i)))
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncoutch, ichvar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//trim(chspec(i)))
#endif
      buffer = 'mass_fraction_of_'//trim(chspec(i))//'_in_air'
      istatus = nf90_put_att(ncoutch, ichvar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(chspec(i)))
      buffer = trim(chspec(i))//' Volume Mixing Ratio'
      istatus = nf90_put_att(ncoutch, ichvar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(chspec(i)))
      istatus = nf90_put_att(ncoutch, ichvar(i+1), 'units', 'kg kg-1')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(chspec(i)))
      istatus = nf90_put_att(ncoutch, ichvar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(chspec(i)))
    end do
    istatus = nf90_enddef(ncoutch)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
    hptop = real(ptop)*10.0
    call write_vertical_coord(ncoutch,sigma2,hptop,izvar)
    call write_horizontal_coord(ncoutch,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncoutch,'xlat',xlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncoutch,'xlon',xlon,ipnt,illvar,no_transpose)

99001 format (a,a,a,a,i10,a)

  end subroutine newfile_ch_icbc

  subroutine newfile_ox_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: i , ipnt , istatus
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ihvar
    integer , dimension(2) :: illvar
    real(sp) , pointer , dimension(:) :: xjx , yiy
    character(64) :: csdate
    real(sp) :: hptop

    if (ncoutox > 0) then
      istatus = nf90_close(ncoutox)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close OX file')
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_OXCL.', toint10(idate1), '.nc'

    irefdate = idate1
    itimeox = 1

    csdate = tochar(idate1)

    call createfile_withname(ofname,ncoutox)
    call add_common_global_params(ncoutox,'oxcl_icbc')

    ipnt = 1
    call define_basic_dimensions(ncoutox,jx,iy,kz,ipnt,idims)
    call add_dimension(ncoutox,'time',nf90_unlimited,ipnt,idims)
!
    call define_horizontal_coord(ncoutox,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncoutox,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncoutox,idims,ipnt,illvar)

    istatus = nf90_def_var(ncoutox, 'time', nf90_double, idims(4:4), ioxvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncoutox, ioxvar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncoutox, ioxvar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , noxsp
      istatus = nf90_def_var(ncoutox, oxspec(i), nf90_float, idims, ioxvar(i+1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding variable '//trim(oxspec(i)))
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncoutox, ioxvar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//trim(oxspec(i)))
#endif
      buffer = 'mass_fraction_of_'//trim(oxspec(i))//'_in_air'
      istatus = nf90_put_att(ncoutox, ioxvar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(oxspec(i)))
      buffer = trim(oxspec(i))//' Volume Mixing Ratio'
      istatus = nf90_put_att(ncoutox, ioxvar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(oxspec(i)))
      istatus = nf90_put_att(ncoutox, ioxvar(i+1), 'units', 'kg kg-1')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(oxspec(i)))
      istatus = nf90_put_att(ncoutox, ioxvar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(oxspec(i)))
    end do
    istatus = nf90_enddef(ncoutox)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
    hptop = real(ptop)*10.0
    call write_vertical_coord(ncoutox,sigma2,hptop,izvar)
    call write_horizontal_coord(ncoutox,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncoutox,'xlat',xlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncoutox,'xlon',xlon,ipnt,illvar,no_transpose)

99001 format (a,a,a,a,i10,a)

  end subroutine newfile_ox_icbc

  subroutine newfile_ae_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: i , ipnt , istatus
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ihvar
    integer , dimension(2) :: illvar
    real(sp) , pointer , dimension(:) :: xjx , yiy
    character(len=64) :: csdate
    character(len=6) :: dustname
    real(sp) :: hptop

    if (ncoutae > 0) then
      istatus = nf90_close(ncoutae)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close AE file')
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_AEBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itimeae = 1

    csdate = tochar(idate1)

    call createfile_withname(ofname,ncoutae)
    call add_common_global_params(ncoutae,'aerc_icbc')

    ipnt = 1
    call define_basic_dimensions(ncoutae,jx,iy,kz,ipnt,idims)
    call add_dimension(ncoutae,'time',nf90_unlimited,ipnt,idims)
!
    call define_horizontal_coord(ncoutae,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncoutae,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncoutae,idims,ipnt,illvar)

    istatus = nf90_def_var(ncoutae, 'time', nf90_double, idims(4:4), iaevar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncoutae, iaevar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncoutae, iaevar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , naesp
      if ( aespec(i) == 'SOA' ) cycle
      if ( aespec(i) == 'SSLT03' ) cycle
      if ( aespec(i) == 'SSLT04' ) cycle
      if ( aespec(i)(1:3) == 'DST' ) then
        dustname = 'DUST'//aespec(i)(4:5)
        istatus = nf90_def_var(ncoutae,dustname,nf90_float,idims,iaevar(i+1))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error adding variable '//trim(aespec(i)))
      else
        istatus = nf90_def_var(ncoutae,aespec(i),nf90_float,idims,iaevar(i+1))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error adding variable '//trim(aespec(i)))
      end if
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncoutae, iaevar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//trim(aespec(i)))
#endif
      buffer = 'mass_fraction_of_'//trim(aespec(i))//'_in_air'
      istatus = nf90_put_att(ncoutae, iaevar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(aespec(i)))
      buffer = trim(aespec(i))//' Volume Mixing Ratio'
      istatus = nf90_put_att(ncoutae, iaevar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(aespec(i)))
      istatus = nf90_put_att(ncoutae, iaevar(i+1), 'units', 'kg kg-1')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(aespec(i)))
      istatus = nf90_put_att(ncoutae, iaevar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(aespec(i)))
    end do
    istatus = nf90_enddef(ncoutae)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
    hptop = real(ptop)*10.0
    call write_vertical_coord(ncoutae,sigma2,hptop,izvar)
    call write_horizontal_coord(ncoutae,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncoutae,'xlat',xlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncoutae,'xlon',xlon,ipnt,illvar,no_transpose)

99001 format (a,a,a,a,i10,a)

  end subroutine newfile_ae_icbc

  subroutine write_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: i , istatus
    integer , dimension(1) :: istart1 , icount1
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: xdate
    type(rcm_time_interval) :: tdif
!
    istart1(1) = itimech
    icount1(1) = 1
    tdif = idate - irefdate
    xdate(1) = tohours(tdif)
    istatus = nf90_put_var(ncoutch, ichvar(1), xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')

    istart(4) = itimech
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx

    do i = 1 , nchsp
      istatus = nf90_put_var(ncoutch,ichvar(i+1),chv4(:,:,:,i),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error variable '//chspec(i)//' write')
    end do

    write (stdout ,*) 'Write ch_icbc : ', tochar(idate)

    itimech = itimech + 1
!
  end subroutine write_ch_icbc

  subroutine write_ox_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: i , istatus
    integer , dimension(1) :: istart1 , icount1
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: xdate
    type(rcm_time_interval) :: tdif
!
    istart1(1) = itimeox
    icount1(1) = 1
    tdif = idate - irefdate
    xdate(1) = tohours(tdif)
    istatus = nf90_put_var(ncoutox, ioxvar(1), xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')

    istart(4) = itimeox
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx

    do i = 1 , noxsp
      istatus = nf90_put_var(ncoutox,ioxvar(i+1),oxv4(:,:,:,i),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error variable '//oxspec(i)//' write')
    end do

    write (stdout ,*) 'Write ox_icbc : ', tochar(idate)

    itimeox = itimeox + 1
!
  end subroutine write_ox_icbc
!
  subroutine write_ae_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: i , istatus
    integer , dimension(1) :: istart1 , icount1
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: xdate
    type(rcm_time_interval) :: tdif
!
    istart1(1) = itimeae
    icount1(1) = 1
    tdif = idate - irefdate
    xdate(1) = tohours(tdif)
    istatus = nf90_put_var(ncoutae, iaevar(1), xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')

    istart(4) = itimeae
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx

    if ( sum_sslt_bins ) then
      aev4(:,:,:,isslt1) = aev4(:,:,:,isslt1) + aev4(:,:,:,isslt2)
      aev4(:,:,:,isslt2) = aev4(:,:,:,isslt3) + aev4(:,:,:,isslt4)
    end if

    if ( sum_soa_to_oc2 ) then
      aev4(:,:,:,ioc2) = aev4(:,:,:,ioc2) + aev4(:,:,:,isoa)
    end if

    do i = 1 , naesp
      if ( aespec(i) == 'SOA' ) cycle
      if ( aespec(i) == 'SSLT03' ) cycle
      if ( aespec(i) == 'SSLT04' ) cycle
      istatus = nf90_put_var(ncoutae,iaevar(i+1),aev4(:,:,:,i),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error variable '//aespec(i)//' write')
    end do

    write (stdout ,*) 'Write ae_icbc : ', tochar(idate)

    itimeae = itimeae + 1
!
  end subroutine write_ae_icbc

end module mod_wrtoxd
