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

  public :: chv4 , oxv4
  public :: nchsp , noxsp
  public :: chspec , oxspec

  public :: init_outoxd , close_outoxd , newfile_ch_icbc , &
            newfile_ch_oxcl , write_ch_icbc , write_ch_oxcl

  integer :: ncout , ncoutox
  character(256) :: ofname
  type(rcm_time_and_date) , save :: irefdate
  integer :: itimech , itimeox
  integer , dimension(4) :: idims
  integer :: istatus

  integer , parameter :: nchsp = 25
  integer , parameter :: noxsp = 5

  integer , dimension(nchsp+1) :: ichvar
  integer , dimension(noxsp+1) :: ioxvar

  character(len=8) , dimension(nchsp) :: chspec
  character(len=8) , dimension(noxsp) :: oxspec

  real(sp) , pointer , dimension(:,:,:,:) :: chv4
  real(sp) , pointer , dimension(:,:,:,:) :: oxv4
  
  character(len=128) :: buffer

  data oxspec / 'OH' , 'HO2' , 'O3' , 'NO3' , 'H2O2' /
  data chspec / 'O3' , 'NO' , 'NO2' , 'HNO3' , 'N2O5' , 'H2O2' , 'CH4' , &
                'CO' , 'CH2O' , 'CH3OH' , 'C2H5OH' , 'C2H4' , 'C2H6' ,   &
                'CH3CHO' , 'CH3COCH3' , 'BIGENE' , 'BIGALK' , 'C3H6' ,   &
                'C3H8' , 'ISOP' , 'TOLUENE' , 'PAN' , 'SO2' , 'SO4' , 'DMS' /

  data ncout   /-1/
  data ncoutox /-1/

  contains

  subroutine init_outoxd
    implicit none
    call getmem4d(chv4,1,jx,1,iy,1,kz,1,nchsp,'mod_wrtoxd:chv4')
    call getmem4d(oxv4,1,jx,1,iy,1,kz,1,noxsp,'mod_wrtoxd:oxv4')
  end subroutine init_outoxd

  subroutine close_outoxd
    implicit none
    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
    end if
    if (ncoutox > 0) then
      istatus = nf90_close(ncoutox)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
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

    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_CHBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itimech = 1

    csdate = tochar(idate1)

    call createfile_withname(ofname,ncout)
    call add_common_global_params(ncout,'chem_icbc')

    ipnt = 1
    call define_basic_dimensions(ncout,jx,iy,kz,ipnt,idims)
    call add_dimension(ncout,'time',nf90_unlimited,ipnt,idims)
!
    call define_horizontal_coord(ncout,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncout,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncout,idims,ipnt,illvar)

    istatus = nf90_def_var(ncout, 'time', nf90_double, idims(4:4), ichvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncout, ichvar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncout, ichvar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , nchsp
      istatus = nf90_def_var(ncout, chspec(i), nf90_float, idims, ichvar(i+1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding variable '//trim(chspec(i)))
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncout, ichvar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//trim(chspec(i)))
#endif
      buffer = 'mass_fraction_of_'//trim(chspec(i))//'_in_air'
      istatus = nf90_put_att(ncout, ichvar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(chspec(i)))
      buffer = trim(chspec(i))//' Volume Mixing Ratio'
      istatus = nf90_put_att(ncout, ichvar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(chspec(i)))
      istatus = nf90_put_att(ncout, ichvar(i+1), 'units', 'kg kg-1')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(chspec(i)))
      istatus = nf90_put_att(ncout, ichvar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(chspec(i)))
    end do
    istatus = nf90_enddef(ncout)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
    hptop = real(ptop)*10.0
    call write_vertical_coord(ncout,sigma2,hptop,izvar)
    call write_horizontal_coord(ncout,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncout,'xlat',xlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'xlon',xlon,ipnt,illvar,no_transpose)

99001 format (a,a,a,a,i10,a)

  end subroutine newfile_ch_icbc

  subroutine newfile_ch_oxcl(idate1)
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
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
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

  end subroutine newfile_ch_oxcl

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
    istatus = nf90_put_var(ncout, ichvar(1), xdate, istart1, icount1)
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
      istatus = nf90_put_var(ncout, ichvar(i+1), chv4(:,:,:,i), istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error variable '//chspec(i)//' write')
    end do

    write (stdout ,*) 'Write ch_icbc : ', tochar(idate)

    itimech = itimech + 1
!
  end subroutine write_ch_icbc

  subroutine write_ch_oxcl(idate)
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
      istatus = nf90_put_var(ncoutox, ioxvar(i+1), &
                             oxv4(:,:,:,i), istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error variable '//oxspec(i)//' write')
    end do

    write (stdout ,*) 'Write ch_oxcl : ', tochar(idate)

    itimeox = itimeox + 1
!
  end subroutine write_ch_oxcl
!
end module mod_wrtoxd
