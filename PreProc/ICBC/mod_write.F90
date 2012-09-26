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

module mod_write

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_grid
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_ensemble
  use netcdf

  private

  integer(ik4) :: ncout
  character(256) :: ofname
  type(rcm_time_and_date) , save :: irefdate
  integer(ik4) :: itime
  integer(ik4) , dimension(4) :: idims
  integer(ik4) , dimension(8) :: ivar

  real(rk4) , pointer , dimension(:,:) :: ps4 , ts4
  real(rk4) , pointer , dimension(:,:,:) :: h4 , q4
  real(rk4) , pointer , dimension(:,:,:) :: t4 , u4 , v4
  real(rk4) , pointer , dimension(:) :: yiy
  real(rk4) , pointer , dimension(:) :: xjx

  public :: ps4 , ts4 , h4 , q4 , t4 , u4 , v4
  public :: init_output , close_output , newfile , writef

  data ncout /-1/

  contains

  subroutine init_output
  implicit none
    call getmem2d(ps4,1,jx,1,iy,'mod_write:ps4')
    call getmem2d(ts4,1,jx,1,iy,'mod_write:ts4')
    call getmem3d(h4,1,jx,1,iy,1,kz,'mod_write:h4')
    call getmem3d(q4,1,jx,1,iy,1,kz,'mod_write:q4')
    call getmem3d(t4,1,jx,1,iy,1,kz,'mod_write:t4')
    call getmem3d(u4,1,jx,1,iy,1,kz,'mod_write:u4')
    call getmem3d(v4,1,jx,1,iy,1,kz,'mod_write:v4')
    call getmem1d(yiy,1,iy,'mod_write:yiy')
    call getmem1d(xjx,1,jx,'mod_write:xjx')
  end subroutine init_output

  subroutine close_output
    implicit none
    integer(ik4) :: istatus
    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error closing file '//trim(ofname))
    end if
  end subroutine close_output

  subroutine newfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer(ik4) :: ipnt , istatus
    integer(ik4) , dimension(2) :: izvar
    integer(ik4) , dimension(2) :: ihvar
    integer(ik4) , dimension(3) :: ihdims
    integer(ik4) , dimension(6) :: illvar
    character(64) :: csdate
    real(rk4) :: hptop

    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error closing file '//trim(ofname))
    end if

    write (ofname,99001) trim(dirglob), pthsep, &
           trim(domname), '_ICBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itime = 1

    csdate = tochar(idate1)

    call createfile_withname(ofname,ncout)
    call add_common_global_params(ncout,'icbc')

    istatus = nf90_put_att(ncout, nf90_global, 'global_sst_source', ssttyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global sst_source')
    istatus = nf90_put_att(ncout, nf90_global, 'global_atm_source', dattyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global atm_source')

    ipnt = 1
    call define_basic_dimensions(ncout,jx,iy,kz,ipnt,idims)
    call add_dimension(ncout,'time',nf90_unlimited,ipnt,idims)

    ihdims(1) = idims(1)
    ihdims(2) = idims(2)
    ihdims(3) = idims(4)
!
    call define_horizontal_coord(ncout,jx,iy,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncout,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncout,idims,ipnt,illvar)
    call define_dot_geolocation_coord(ncout,idims,ipnt,illvar)
    call define_topo_and_mask(ncout,idims,ipnt,illvar)

    istatus = nf90_def_var(ncout, 'time', nf90_double, idims(4:4), ivar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncout, ivar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncout, ivar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    istatus = nf90_def_var(ncout, 'ps', nf90_float, ihdims, ivar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ps')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(2), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on ps')
#endif
    istatus = nf90_put_att(ncout, ivar(2), 'standard_name', &
                           'surface_air_pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ps standard_name')
    istatus = nf90_put_att(ncout, ivar(2), 'long_name', 'Surface pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ps long_name')
    istatus = nf90_put_att(ncout, ivar(2), 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ps units')
    istatus = nf90_put_att(ncout, ivar(2), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ps coordinates')
    istatus = nf90_def_var(ncout, 'ts', nf90_float, ihdims, ivar(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ts')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(3), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on ts')
#endif
    istatus = nf90_put_att(ncout, ivar(3), 'standard_name', &
                           'surface_temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ts standard_name')
    istatus = nf90_put_att(ncout, ivar(3), 'long_name', 'Surface Temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ts long_name')
    istatus = nf90_put_att(ncout, ivar(3), 'units', 'K')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ts units')
    istatus = nf90_put_att(ncout, ivar(3), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ts coordinates')
    istatus = nf90_def_var(ncout, 'u', nf90_float, idims, ivar(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable u')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(4), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on u')
#endif
    istatus = nf90_put_att(ncout, ivar(4), 'standard_name', 'eastward_wind')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding u standard_name')
    istatus = nf90_put_att(ncout, ivar(4), 'long_name',     &
                           'U component (westerly) of wind')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding u long_name')
    istatus = nf90_put_att(ncout, ivar(4), 'units', 'm s-1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding u units')
    istatus = nf90_put_att(ncout, ivar(4), 'coordinates', 'dlon dlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding u coordinates')
    istatus = nf90_def_var(ncout, 'v', nf90_float, idims, ivar(5))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable v')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(5), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on v')
#endif
    istatus = nf90_put_att(ncout, ivar(5), 'standard_name', 'northward_wind')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding v standard_name')
    istatus = nf90_put_att(ncout, ivar(5), 'long_name',     &
                           'V component (southerly) of wind')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding v long_name')
    istatus = nf90_put_att(ncout, ivar(5), 'units', 'm s-1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding v units')
    istatus = nf90_put_att(ncout, ivar(5), 'coordinates', 'dlon dlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding v coordinates')
    istatus = nf90_def_var(ncout, 't', nf90_float, idims, ivar(6))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable t')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(6), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on t')
#endif
    istatus = nf90_put_att(ncout, ivar(6), 'standard_name', &
                           'air_temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding t standard_name')
    istatus = nf90_put_att(ncout, ivar(6), 'long_name', 'Temperature')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding t long_name')
    istatus = nf90_put_att(ncout, ivar(6), 'units', 'K')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding t units')
    istatus = nf90_put_att(ncout, ivar(6), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding t coordinates')
    istatus = nf90_def_var(ncout, 'qv', nf90_float, idims, ivar(7))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable qv')
#ifdef NETCDF4_HDF5
    istatus = nf90_def_var_deflate(ncout, ivar(7), 1, 1, 9)
    call checkncerr(istatus,__FILE__,__LINE__,'Error setting compression on qv')
#endif
    istatus = nf90_put_att(ncout, ivar(7), 'standard_name', &
                           'humidity_mixing_ratio')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding qv standard_name')
    istatus = nf90_put_att(ncout, ivar(7), 'long_name',     &
                           'Water vapor mixing ratio')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding qv long_name')
    istatus = nf90_put_att(ncout, ivar(7), 'units', 'kg kg-1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding qv units')
    istatus = nf90_put_att(ncout, ivar(7), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding qv coordinates')
!
    istatus = nf90_enddef(ncout)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
!
    hptop = real(ptop*10.0D0)
    call write_vertical_coord(ncout,sigma2,hptop,izvar)
    call write_horizontal_coord(ncout,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncout,'xlat',xlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'xlon',xlon,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'dlat',dlat,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'dlon',dlon,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'topo',topogm,ipnt,illvar,no_transpose)
    call write_var2d_static(ncout,'mask',mask,ipnt,illvar,no_transpose)

    istatus = nf90_sync(ncout)
    call checkncerr(istatus,__FILE__,__LINE__,'Error file sync')

99001 format (a,a,a,a,i10,a)

  end subroutine newfile

  subroutine writef(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdiff
    integer(ik4) :: istatus
    integer(ik4) , dimension(1) :: istart1 , icount1
    integer(ik4) , dimension(4) :: istart , icount
    real(rk8) , dimension(1) :: xdate
!
    if ( ensemble_run ) then
      write(stdout,*) 'Appling perturbation to input dataset:'
      if ( lperturb_ts ) then
        write(stdout,'(a,f7.2,a)') 'TS with value ',perturb_frac_ts*d_100,'%'
        call randify(ts4,perturb_frac_ts,jx,iy)
      end if
      if ( lperturb_ps ) then
        write(stdout,'(a,f7.2,a)') 'PS with value ',perturb_frac_ps*d_100,'%'
        call randify(ps4,perturb_frac_ps,jx,iy)
      end if
      if ( lperturb_t ) then
        write(stdout,'(a,f7.2,a)') 'T  with value ',perturb_frac_t*d_100,'%'
        call randify(t4,perturb_frac_t,jx,iy,kz)
      end if
      if ( lperturb_q ) then
        write(stdout,'(a,f7.2,a)') 'Q  with value ',perturb_frac_q*d_100,'%'
        call randify(q4,perturb_frac_q,jx,iy,kz)
      end if
      if ( lperturb_u ) then
        write(stdout,'(a,f7.2,a)') 'U  with value ',perturb_frac_u*d_100,'%'
        call randify(u4,perturb_frac_u,jx,iy,kz)
      end if
      if ( lperturb_v ) then
        write(stdout,'(a,f7.2,a)') 'V  with value ',perturb_frac_v*d_100,'%'
        call randify(v4,perturb_frac_v,jx,iy,kz)
      end if
    end if
    istart1(1) = itime
    icount1(1) = 1
    tdiff = idate - irefdate
    xdate(1) = tohours(tdiff)
    istatus = nf90_put_var(ncout, ivar(1), xdate, istart1, icount1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable time write')
    istart(3) = itime
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    ps4 = real((dble(ps4)+ptop)*10.0D0)
    istatus = nf90_put_var(ncout, ivar(2), ps4, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ps write')
    istatus = nf90_put_var(ncout, ivar(3), ts4, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ts write')
    istart(4) = itime
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_put_var(ncout, ivar(4), u4, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable u write')
    istatus = nf90_put_var(ncout, ivar(5), v4, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable v write')
    istatus = nf90_put_var(ncout, ivar(6), t4, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable t write')
    istatus = nf90_put_var(ncout, ivar(7), q4, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable qv write')
    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncout)
      call checkncerr(istatus,__FILE__,__LINE__,'Error sync output file')
    end if
    itime = itime + 1
!
  end subroutine writef
!
end module mod_write
