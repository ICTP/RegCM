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

module mod_sst_grid

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_domain
  use netcdf

  private

  integer(ik4) :: ncid
  integer(ik4) , dimension(4) :: idims
  integer(ik4) , dimension(4) :: ivar
  type (rcm_time_and_date) , save :: refdate
  integer(ik4) :: itime

  real(rk4) , public , pointer , dimension(:,:) :: sstmm , icemm
  real(rk4) , public , pointer , dimension(:,:) :: xlat , xlon
  real(rk4) , pointer , dimension(:,:) :: topo , mask
  real(rk4) , pointer , dimension(:) :: sigma

  public :: init_grid , read_domain_info , open_sstfile , &
            close_sstfile , writerec

  contains

  subroutine init_grid
    implicit none
    call getmem2d(sstmm,1,iy,1,jx,'mod_sst_grid:sstmm')
    call getmem2d(icemm,1,iy,1,jx,'mod_sst_grid:icemm')
    call getmem2d(xlat,1,iy,1,jx,'mod_sst_grid:xlat')
    call getmem2d(xlon,1,iy,1,jx,'mod_sst_grid:xlon')
    call getmem2d(topo,1,iy,1,jx,'mod_sst_grid:topo')
    call getmem2d(mask,1,iy,1,jx,'mod_sst_grid:mask')
    call getmem1d(sigma,1,kzp1,'mod_sst_grid:sigma')
  end subroutine init_grid

  subroutine read_domain_info(terfile)
    implicit none
    character(256) , intent(in) :: terfile
    integer(ik4) :: istatus , incin
    call openfile_withname(terfile,incin)
    call read_domain(incin,sigma,xlat,xlon,ht=topo,mask=mask,ltrans=.true.)
    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error closing file '//trim(terfile))
  end subroutine read_domain_info

  subroutine open_sstfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer(ik4) :: istatus
    character(256) :: sstname
    character(64) :: csdate
    real(rk4) :: hptop
    integer(ik4) , dimension(2) :: izvar
    integer(ik4) , dimension(2) :: ihvar
    integer(ik4) , dimension(4) :: illvar
    integer(ik4) , dimension(3) :: ixdims
    real(rk4) , pointer , dimension(:) :: yiy
    real(rk4) , pointer , dimension(:) :: xjx
    integer(ik4) :: ipnt

    refdate = idate1
    csdate = tochar(refdate)
    itime = 1

    sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'

    call createfile_withname(sstname,ncid)
    call add_common_global_params(ncid,'sst',.false.)

    istatus = nf90_put_att(ncid, nf90_global, 'sst_source', ssttyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global sst_source')
!
    ipnt = 1
    call define_basic_dimensions(ncid,jx,iy,kz+1,ipnt,idims)
    call add_dimension(ncid,'time',nf90_unlimited,ipnt,idims)
    ixdims(1) = idims(1)
    ixdims(2) = idims(2)
    ixdims(3) = idims(4)

    call define_horizontal_coord(ncid,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncid,idims,izvar)
!
    ipnt = 1
    call define_cross_geolocation_coord(ncid,idims,ipnt,illvar)
    call define_topo_and_mask(ncid,idims,ipnt,illvar)

    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(4:4), ivar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncid, ivar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncid, ivar(1), 'calendar', calstr(refdate%calendar))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')
    istatus = nf90_def_var(ncid, 'sst', nf90_float, ixdims, ivar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sst')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncid, ivar(2), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting deflate on sst')
#endif
    istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                           'sea_surface_temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sst standard_name')
    istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                           'Sea Surface Temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sst long_name')
    istatus = nf90_put_att(ncid, ivar(2), 'units', 'K')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sst units')
    istatus = nf90_put_att(ncid, ivar(2), '_FillValue', -9999.0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sst _FillValue')
    istatus = nf90_put_att(ncid, ivar(2), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sst coordinates')
    if ( ssttyp == 'OI2ST' .or. ssttyp == 'OI2WK' ) then
      istatus = nf90_def_var(ncid, 'ice', nf90_float, ixdims, ivar(3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ice')
      istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
                             'sea_ice_thickness')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding ice standard_name')
      istatus = nf90_put_att(ncid, ivar(3), 'long_name', 'Sea ice depth')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding ice long_name')
      istatus = nf90_put_att(ncid, ivar(3), 'units', 'mm')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding ice units')
      istatus = nf90_put_att(ncid, ivar(3), '_FillValue', -9999.0)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding ice _FillValue')
      istatus = nf90_put_att(ncid, ivar(3), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding ice coordinates')
    end if
!
    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error End Definitions NetCDF output')
!
    hptop = real(ptop*10.0D0)
    call write_vertical_coord(ncid,sigma,hptop,izvar)
    call write_horizontal_coord(ncid,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncid,'xlat',xlat,ipnt,illvar,do_transpose)
    call write_var2d_static(ncid,'xlon',xlon,ipnt,illvar,do_transpose)
    call write_var2d_static(ncid,'topo',topo,ipnt,illvar,do_transpose)
    call write_var2d_static(ncid,'mask',mask,ipnt,illvar,do_transpose)

  end subroutine open_sstfile

  subroutine close_sstfile
    implicit none
    integer(ik4) :: istatus
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error closing output file')
  end subroutine close_sstfile

  subroutine writerec(idate,lice)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: lice
    integer(ik4) :: istatus
    integer(ik4) , dimension(1) :: istart1 , icount1
    integer(ik4) , dimension(3) :: istart , icount
    real(rk8) , dimension(1) :: xdate
    type(rcm_time_interval) :: tdiff

    istart(3) = itime
    istart(2) = 1
    istart(1) = 1
    istart1(1) = itime
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    icount1(1) = 1
    tdiff = idate - refdate
    xdate(1) = tohours(tdiff)
    istatus = nf90_put_var(ncid, ivar(1), xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')
    istatus = nf90_put_var(ncid, ivar(2), transpose(sstmm), istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable sst write')
    if (lice) then
      istatus = nf90_put_var(ncid, ivar(3), transpose(icemm), istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error variable ice write')
    end if
    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error sync output file')
    end if
    itime = itime + 1
  end subroutine writerec

end module mod_sst_grid
