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

  use mod_realkinds
  use netcdf
  use mod_dynparam
  use mod_grid
  use mod_memutil
  use mod_message

  private

  integer :: ncout
  character(256) :: ofname
  type(rcm_time_and_date) , save :: irefdate
  integer :: itime
  integer , dimension(5) :: idims
  integer , dimension(8) :: ivar

  real(sp) , pointer , dimension(:,:) :: ps4 , ts4
  real(sp) , pointer , dimension(:,:,:) :: h4 , q4
  real(sp) , pointer , dimension(:,:,:) :: t4 , u4 , v4
  real(sp) , pointer , dimension(:,:,:) :: sulfate4
  real(sp) , pointer , dimension(:) :: yiy
  real(sp) , pointer , dimension(:) :: xjx

  public :: ps4 , ts4 , h4 , q4 , t4 , u4 , v4 , sulfate4
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
    integer :: istatus
    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__,('Error closing file '//trim(ofname)))
    end if
  end subroutine close_output

  subroutine newfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: istatus
    integer :: i , j
    integer , dimension(8) :: tvals
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(3) :: illvar
    integer , dimension(4) :: x3ddim
    character(64) :: csdate , cdum
    character(256) :: history
    real(sp) , dimension(2) :: trlat
    real(sp) :: hptop


    if (ncout > 0) then
      istatus = nf90_close(ncout)
      call checkncerr(istatus,__FILE__,__LINE__,('Error closing file '//trim(ofname)))
    end if

    write (ofname,99001) trim(dirglob), pthsep, &
           trim(domname), '_ICBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itime = 1

    csdate = tochar(idate1)

#ifdef NETCDF4_HDF5
    istatus = nf90_create(ofname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), ncout)
#else
    istatus = nf90_create(ofname, nf90_clobber, ncout)
#endif
    call checkncerr(istatus,__FILE__,__LINE__,('Error creating file '//trim(ofname)))

    istatus = nf90_put_att(ncout, nf90_global, 'title',  &
              'ICTP Regional Climatic model V4 ICBC program output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global title')
    istatus = nf90_put_att(ncout, nf90_global, 'institution', 'ICTP')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global institution')
    istatus = nf90_put_att(ncout, nf90_global, 'source', &
               'RegCM Model simulation SST output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global source')
    istatus = nf90_put_att(ncout, nf90_global, 'Conventions', 'CF-1.4')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global Conventions')
    call date_and_time(values = tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
         ' : Created by RegCM icbc program'
    istatus = nf90_put_att(ncout, nf90_global, 'history', history)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global history')
    istatus = nf90_put_att(ncout, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global references')
    istatus = nf90_put_att(ncout, nf90_global, 'experiment', domname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global experiment')
    istatus = nf90_put_att(ncout, nf90_global, 'projection', iproj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global projection')
    istatus = nf90_put_att(ncout, nf90_global, 'grid_size_in_meters', ds*1000.0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global gridsize')
    istatus = nf90_put_att(ncout, nf90_global,   &
                 'latitude_of_projection_origin', clat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clat')
    istatus = nf90_put_att(ncout, nf90_global,   &
                 'longitude_of_projection_origin', clon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clon')
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncout, nf90_global, &
                   'grid_north_pole_latitude', plat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plat')
      istatus = nf90_put_att(ncout, nf90_global, &
                   'grid_north_pole_longitude', plon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plon')
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncout, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global truelat')
    end if
    istatus = nf90_put_att(ncout, nf90_global, 'global_data_source', dattyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global data_source')
    istatus = nf90_put_att(ncout, nf90_global, 'sulfate_data_present', cdum)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global sulfate_present')
    istatus = nf90_def_dim(ncout, 'iy', iy, idims(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension iy')
    istatus = nf90_def_dim(ncout, 'jx', jx, idims(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension jx')
    istatus = nf90_def_dim(ncout, 'time', nf90_unlimited, idims(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension time')
    istatus = nf90_def_dim(ncout, 'kz', kz, idims(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension kz')
    x3ddim(1) = idims(1)
    x3ddim(2) = idims(2)
    x3ddim(3) = idims(4)
    x3ddim(4) = idims(3)
!
    istatus = nf90_def_var(ncout, 'sigma', nf90_float, idims(4), izvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sigma')
    istatus = nf90_put_att(ncout, izvar(1), 'standard_name',  &
                           'atmosphere_sigma_coordinate')      
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma standard_name')
    istatus = nf90_put_att(ncout, izvar(1), 'long_name',      &
                           'Sigma at model layers')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma long_name')
    istatus = nf90_put_att(ncout, izvar(1), 'units', '1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma units')
    istatus = nf90_put_att(ncout, izvar(1), 'axis', 'Z')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma axis')
    istatus = nf90_put_att(ncout, izvar(1), 'positive', 'down')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma positive')
    istatus = nf90_put_att(ncout, izvar(1), 'formula_terms',  &
                           'sigma: sigma ps: ps ptop: ptop')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma formula_terms')
    istatus = nf90_def_var(ncout, 'ptop', nf90_float, varid=izvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ptop')
    istatus = nf90_put_att(ncout, izvar(2), 'standard_name', 'air_pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop standard_name')
    istatus = nf90_put_att(ncout, izvar(2), 'long_name',      &
                           'Pressure at model top')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop long_name')
    istatus = nf90_put_att(ncout, izvar(2), 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop units')
    istatus = nf90_def_var(ncout, 'iy', nf90_float, idims(2), ivvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable iy')
    istatus = nf90_put_att(ncout, ivvar(1), 'standard_name',  &
                           'projection_y_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy standard_name')
    istatus = nf90_put_att(ncout, ivvar(1), 'long_name',      &
                           'y-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy long_name')
    istatus = nf90_put_att(ncout, ivvar(1), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy units')
    istatus = nf90_def_var(ncout, 'jx', nf90_float, idims(1), ivvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable jx')
    istatus = nf90_put_att(ncout, ivvar(2), 'standard_name', &
                           'projection_x_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx standard_name')
    istatus = nf90_put_att(ncout, ivvar(2), 'long_name',    &
                           'x-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx long_name')
    istatus = nf90_put_att(ncout, ivvar(2), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx units')
    istatus = nf90_def_var(ncout, 'xlat', nf90_float, idims(1:2), illvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlat')
    istatus = nf90_put_att(ncout, illvar(1), 'standard_name', 'latitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat standard_name')
    istatus = nf90_put_att(ncout, illvar(1), 'long_name',     &
                           'Latitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat long_name')
    istatus = nf90_put_att(ncout, illvar(1), 'units', 'degrees_north')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat units')
    istatus = nf90_def_var(ncout, 'xlon', nf90_float, idims(1:2), illvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlon')
    istatus = nf90_put_att(ncout, illvar(2), 'standard_name', 'longitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon standard_name')
    istatus = nf90_put_att(ncout, illvar(2), 'long_name',     &
                           'Longitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon long_name')
    istatus = nf90_put_att(ncout, illvar(2), 'units', 'degrees_east')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon units')
    istatus = nf90_def_var(ncout, 'topo', nf90_float, idims(1:2), illvar(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable topo')
    istatus = nf90_put_att(ncout, illvar(3), 'standard_name', &
                           'surface_altitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding topo standard_name')
    istatus = nf90_put_att(ncout, illvar(3), 'long_name',     &
                           'Domain surface elevation')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding topo long_name')
    istatus = nf90_put_att(ncout, illvar(3), 'units', 'm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding topo units')
    istatus = nf90_put_att(ncout, illvar(3), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding topo coordinates')
    istatus = nf90_def_var(ncout, 'time', nf90_double, idims(3:3), ivar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncout, ivar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncout, ivar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')
    istatus = nf90_def_var(ncout, 'ps', nf90_float, idims(1:3), ivar(2))
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
    istatus = nf90_def_var(ncout, 'ts', nf90_float, idims(1:3), ivar(3))
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
    istatus = nf90_def_var(ncout, 'u', nf90_float, x3ddim, ivar(4))
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
    istatus = nf90_put_att(ncout, ivar(4), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding u coordinates')
    istatus = nf90_def_var(ncout, 'v', nf90_float, x3ddim, ivar(5))
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
    istatus = nf90_put_att(ncout, ivar(5), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding v coordinates')
    istatus = nf90_def_var(ncout, 't', nf90_float, x3ddim, ivar(6))
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
    istatus = nf90_def_var(ncout, 'qv', nf90_float, x3ddim, ivar(7))
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
    call checkncerr(istatus,__FILE__,__LINE__,'Error End Definitions NetCDF output')
!
    istatus = nf90_put_var(ncout, izvar(1), sigma2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable sigma write')
    hptop = real(ptop*10.0D0)
    istatus = nf90_put_var(ncout, izvar(2), hptop)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ptop write')
    yiy(1) = -real((dble(iy-1)/2.0D0) * ds)
    xjx(1) = -real((dble(jx-1)/2.0D0) * ds)
    do i = 2 , iy
      yiy(i) = real(dble(yiy(i-1))+ds)
    end do
    do j = 2 , jx
      xjx(j) = real(dble(xjx(j-1))+ds)
    end do
    istatus = nf90_put_var(ncout, ivvar(1), yiy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable iy write')
    istatus = nf90_put_var(ncout, ivvar(2), xjx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable jx write')
    istatus = nf90_put_var(ncout, illvar(1), xlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlat write')
    istatus = nf90_put_var(ncout, illvar(2), xlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlon write')
    istatus = nf90_put_var(ncout, illvar(3), topogm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlon write')

    istatus = nf90_sync(ncout)
    call checkncerr(istatus,__FILE__,__LINE__,'Error file sync')


99001 format (a,a,a,a,i10,a)

  end subroutine newfile

  subroutine writef(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdiff
    integer :: istatus
    integer , dimension(1) :: istart1 , icount1
    integer , dimension(4) :: istart , icount
    real(dp) , dimension(1) :: xdate
!
!
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
