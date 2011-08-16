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

  use mod_dynparam
  use mod_date
  use mod_grid
  use m_realkinds
  use m_die
  use netcdf
  use mod_memutil
  use mod_message

  private

  public :: chv4 , oxv4
  public :: nchsp , noxsp
  public :: chspec , oxspec

  public :: init_outoxd , close_outoxd , newfile_ch_icbc , &
            newfile_ch_oxcl , write_ch_icbc , write_ch_oxcl

  integer :: ncid , ncidox
  character(256) :: ofname
  type(rcm_time_and_date) :: irefdate
  integer :: itimech , itimeox
  integer , dimension(5) :: idims
  integer :: istatus

  integer , parameter :: nchsp = 25
  integer , parameter :: noxsp = 5

  integer , dimension(nchsp+1) :: ichvar
  integer , dimension(noxsp+1) :: ioxvar

  character(len=8) , dimension(nchsp) :: chspec
  character(len=8) , dimension(noxsp) :: oxspec

  real(sp) , pointer , dimension(:,:,:,:) :: chv4
  real(sp) , pointer , dimension(:,:,:,:) :: oxv4
  real(sp) , pointer , dimension(:) :: yiy
  real(sp) , pointer , dimension(:) :: xjx
  real(sp) :: hptop
  
  character(len=128) :: buffer

  data oxspec / 'OH' , 'HO2' , 'O3' , 'NO3' , 'H2O2' /
  data chspec / 'O3' , 'NO' , 'NO2' , 'HNO3' , 'N2O5' , 'H2O2' , 'CH4' , &
                'CO' , 'CH2O' , 'CH3OH' , 'C2H5OH' , 'C2H4' , 'C2H6' ,   &
                'CH3CHO' , 'CH3COCH3' , 'BIGENE' , 'BIGALK' , 'C3H6' ,   &
                'C3H8' , 'ISOP' , 'TOLUENE' , 'PAN' , 'SO2' , 'SO4' , 'DMS' /

  data ncid   /-1/
  data ncidox /-1/

  contains

  subroutine init_outoxd
    implicit none
    integer :: i , j
    call getmem4d(chv4,1,jx,1,iy,1,kz,1,nchsp,'mod_wrtoxd:chv4')
    call getmem4d(oxv4,1,jx,1,iy,1,kz,1,noxsp,'mod_wrtoxd:oxv4')
    call getmem1d(yiy,1,iy,'mod_wrtoxd:yiy')
    call getmem1d(xjx,1,jx,'mod_wrtoxd:xjx')
    yiy(1) = -(real(iy-1)/2.0) * real(ds)
    xjx(1) = -(real(jx-1)/2.0) * real(ds)
    do i = 2 , iy
      yiy(i) = yiy(i-1)+real(ds)
    end do
    do j = 2 , jx
      xjx(j) = xjx(j-1)+real(ds)
    end do
    hptop = real(ptop)*10.0
  end subroutine init_outoxd

  subroutine close_outoxd
    implicit none
    if (ncid > 0) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
    end if
  end subroutine close_outoxd

  subroutine newfile_ch_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: i , istatus
    integer , dimension(8) :: tvals
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(2) :: illvar
    integer , dimension(4) :: x3ddim
    character(64) :: csdate
    character(256) :: history
    real(sp) , dimension(2) :: trlat

    if (ncid > 0) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(ofname))
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_CHBC.', toint10(idate1), '.nc'

    irefdate = idate1
    itimech = 1

    csdate = tochar(idate1)

#ifdef NETCDF4_HDF5
    istatus = nf90_create(ofname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid)
#else
    istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
    call checkncerr(istatus,__FILE__,__LINE__,'Error create '//trim(ofname))

    istatus = nf90_put_att(ncid, nf90_global, 'title',  &
          'ICTP Regional Climatic model V4 oxidant program output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global title')
    istatus = nf90_put_att(ncid, nf90_global, 'institution', 'ICTP')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global institution')
    istatus = nf90_put_att(ncid, nf90_global, 'source', &
               'RegCM Model simulation CHEM_ICBC output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global source')
    istatus = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.4')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global Conventions')
    call date_and_time(values=tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
         ' : Created by RegCM chem_icbc program'
    istatus = nf90_put_att(ncid, nf90_global, 'history', history)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global history')
    istatus = nf90_put_att(ncid, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global references')
    istatus = nf90_put_att(ncid, nf90_global, 'experiment', domname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global experiment')
    istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global projection')
    istatus = nf90_put_att(ncid, nf90_global, 'grid_size_in_meters', ds*1000.0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global gridsize')
    istatus = nf90_put_att(ncid, nf90_global, 'latitude_of_projection_origin', clat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clat')
    istatus = nf90_put_att(ncid, nf90_global, 'longitude_of_projection_origin', clon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clon')
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, nf90_global, 'grid_north_pole_latitude', plat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plat')
      istatus = nf90_put_att(ncid, nf90_global, 'grid_north_pole_longitude', plon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plon')
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global truelat')
    end if
    istatus = nf90_put_att(ncid, nf90_global, 'global_data_source', dattyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global data_source')
    istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension iy')
    istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension jx')
    istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension time')
    istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension kz')
    x3ddim(1) = idims(1)
    x3ddim(2) = idims(2)
    x3ddim(3) = idims(4)
    x3ddim(4) = idims(3)
!
    istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4), izvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sigma')
    istatus = nf90_put_att(ncid, izvar(1), 'standard_name',  &
                           'atmosphere_sigma_coordinate')      
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma standard_name')
    istatus = nf90_put_att(ncid, izvar(1), 'long_name', 'Sigma at model layers')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma long_name')
    istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma units')
    istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma axis')
    istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma positive')
    istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                           'sigma: sigma ps: ps ptop: ptop')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma formula_terms')
    istatus = nf90_def_var(ncid, 'ptop', nf90_float, varid=izvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ptop')
    istatus = nf90_put_att(ncid, izvar(2), 'standard_name', 'air_pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop standard_name')
    istatus = nf90_put_att(ncid, izvar(2), 'long_name', 'Pressure at model top')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop long_name')
    istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop units')
    istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ivvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable iy')
    istatus = nf90_put_att(ncid, ivvar(1), 'standard_name', 'projection_y_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy standard_name')
    istatus = nf90_put_att(ncid, ivvar(1), 'long_name',      &
                           'y-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy long_name')
    istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy units')
    istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ivvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable jx')
    istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', 'projection_x_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx standard_name')
    istatus = nf90_put_att(ncid, ivvar(2), 'long_name',    &
                           'x-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx long_name')
    istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx units')
    istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), illvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlat')
    istatus = nf90_put_att(ncid, illvar(1), 'standard_name', 'latitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat standard_name')
    istatus = nf90_put_att(ncid, illvar(1), 'long_name',  'Latitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat long_name')
    istatus = nf90_put_att(ncid, illvar(1), 'units', 'degrees_north')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat units')
    istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), illvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlon')
    istatus = nf90_put_att(ncid, illvar(2), 'standard_name', 'longitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon standard_name')
    istatus = nf90_put_att(ncid, illvar(2), 'long_name', 'Longitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon long_name')
    istatus = nf90_put_att(ncid, illvar(2), 'units', 'degrees_east')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon units')
    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), ichvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncid, ichvar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncid, ichvar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , nchsp
      istatus = nf90_def_var(ncid, chspec(i), nf90_float, x3ddim, ichvar(i+1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding variable '//trim(chspec(i)))
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncid, ichvar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//trim(chspec(i)))
#endif
      buffer = 'mass_fraction_of_'//trim(chspec(i))//'_in_air'
      istatus = nf90_put_att(ncid, ichvar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(chspec(i)))
      buffer = trim(chspec(i))//' Volume Mixing Ratio'
      istatus = nf90_put_att(ncid, ichvar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(chspec(i)))
      istatus = nf90_put_att(ncid, ichvar(i+1), 'units', 'kg kg-1')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(chspec(i)))
      istatus = nf90_put_att(ncid, ichvar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(chspec(i)))
    end do

    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error End Definitions NetCDF output')
!
    istatus = nf90_put_var(ncid, izvar(1), sigma2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable sigma write')
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ptop write')
    istatus = nf90_put_var(ncid, ivvar(1), yiy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable iy write')
    istatus = nf90_put_var(ncid, ivvar(2), xjx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable jx write')
    istatus = nf90_put_var(ncid, illvar(1), xlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlat write')
    istatus = nf90_put_var(ncid, illvar(2), xlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlon write')

99001 format (a,a,a,a,i10,a)

  end subroutine newfile_ch_icbc

  subroutine newfile_ch_oxcl(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: i , istatus
    integer , dimension(8) :: tvals
    integer , dimension(2) :: izvar
    integer , dimension(2) :: ivvar
    integer , dimension(2) :: illvar
    integer , dimension(4) :: x3ddim
    character(64) :: csdate
    character(256) :: history
    real(sp) , dimension(2) :: trlat
    real(sp) :: hptop

    if (ncidox > 0) then
      istatus = nf90_close(ncidox)
      call checkncerr(istatus,__FILE__,__LINE__,('Error closing file '//trim(ofname)))
    end if

    write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
                '_OXCL.', toint10(idate1), '.nc'

    irefdate = idate1
    itimeox = 1

    csdate = tochar(idate1)

#ifdef NETCDF4_HDF5
    istatus = nf90_create(ofname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), ncidox)
#else
    istatus = nf90_create(ofname, nf90_clobber, ncidox)
#endif
    call checkncerr(istatus,__FILE__,__LINE__,('Error creating file '//trim(ofname)))

    istatus = nf90_put_att(ncidox, nf90_global, 'title',  &
          'ICTP Regional Climatic model V4 chem_icbc program output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global title')
    istatus = nf90_put_att(ncidox, nf90_global, 'institution', 'ICTP')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global institution')
    istatus = nf90_put_att(ncidox, nf90_global, 'source', &
               'RegCM Model simulation SST output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global source')
    istatus = nf90_put_att(ncidox, nf90_global, 'Conventions', 'CF-1.4')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global Conventions')
    call date_and_time(values=tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
         ' : Created by RegCM oxidant program'
    istatus = nf90_put_att(ncidox, nf90_global, 'history', history)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global history')
    istatus = nf90_put_att(ncidox, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global references')
    istatus = nf90_put_att(ncidox, nf90_global, 'experiment', domname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global experiment')
    istatus = nf90_put_att(ncidox, nf90_global, 'projection', iproj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global projection')
    istatus = nf90_put_att(ncidox, nf90_global, 'grid_size_in_meters', ds*1000.0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global gridsize')
    istatus = nf90_put_att(ncidox, nf90_global, 'latitude_of_projection_origin', clat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clat')
    istatus = nf90_put_att(ncidox, nf90_global, 'longitude_of_projection_origin', clon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clon')
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncidox, nf90_global, 'grid_north_pole_latitude', plat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plat')
      istatus = nf90_put_att(ncidox, nf90_global, 'grid_north_pole_longitude', plon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plon')
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncidox, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global truelat')
    end if
    istatus = nf90_put_att(ncidox, nf90_global, 'global_data_source', dattyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global data_source')
    istatus = nf90_def_dim(ncidox, 'iy', iy, idims(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension iy')
    istatus = nf90_def_dim(ncidox, 'jx', jx, idims(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension jx')
    istatus = nf90_def_dim(ncidox, 'time', nf90_unlimited, idims(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension time')
    istatus = nf90_def_dim(ncidox, 'kz', kz, idims(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension kz')
    x3ddim(1) = idims(1)
    x3ddim(2) = idims(2)
    x3ddim(3) = idims(4)
    x3ddim(4) = idims(3)
!
    istatus = nf90_def_var(ncidox, 'sigma', nf90_float, idims(4), izvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sigma')
    istatus = nf90_put_att(ncidox, izvar(1), 'standard_name', &
                           'atmosphere_sigma_coordinate')      
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma standard_name')
    istatus = nf90_put_att(ncidox, izvar(1), 'long_name', 'Sigma at model layers')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma long_name')
    istatus = nf90_put_att(ncidox, izvar(1), 'units', '1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma units')
    istatus = nf90_put_att(ncidox, izvar(1), 'axis', 'Z')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma axis')
    istatus = nf90_put_att(ncidox, izvar(1), 'positive', 'down')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma positive')
    istatus = nf90_put_att(ncidox, izvar(1), 'formula_terms',  &
                           'sigma: sigma ps: ps ptop: ptop')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma formula_terms')
    istatus = nf90_def_var(ncidox, 'ptop', nf90_float, varid=izvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable ptop')
    istatus = nf90_put_att(ncidox, izvar(2), 'standard_name', 'air_pressure')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop standard_name')
    istatus = nf90_put_att(ncidox, izvar(2), 'long_name', 'Pressure at model top')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop long_name')
    istatus = nf90_put_att(ncidox, izvar(2), 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop units')
    istatus = nf90_def_var(ncidox, 'iy', nf90_float, idims(2), ivvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable iy')
    istatus = nf90_put_att(ncidox, ivvar(1), 'standard_name',  &
                           'projection_y_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy standard_name')
    istatus = nf90_put_att(ncidox, ivvar(1), 'long_name',      &
                           'y-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy long_name')
    istatus = nf90_put_att(ncidox, ivvar(1), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy units')
    istatus = nf90_def_var(ncidox, 'jx', nf90_float, idims(1), ivvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable jx')
    istatus = nf90_put_att(ncidox, ivvar(2), 'standard_name', &
                           'projection_x_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx standard_name')
    istatus = nf90_put_att(ncidox, ivvar(2), 'long_name',     &
                           'x-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx long_name')
    istatus = nf90_put_att(ncidox, ivvar(2), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx units')
    istatus = nf90_def_var(ncidox, 'xlat', nf90_float, idims(1:2), illvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlat')
    istatus = nf90_put_att(ncidox, illvar(1), 'standard_name', 'latitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat standard_name')
    istatus = nf90_put_att(ncidox, illvar(1), 'long_name', 'Latitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat long_name')
    istatus = nf90_put_att(ncidox, illvar(1), 'units', 'degrees_north')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat units')
    istatus = nf90_def_var(ncidox, 'xlon', nf90_float, idims(1:2), illvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlon')
    istatus = nf90_put_att(ncidox, illvar(2), 'standard_name', 'longitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon standard_name')
    istatus = nf90_put_att(ncidox, illvar(2), 'long_name', 'Longitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon long_name')
    istatus = nf90_put_att(ncidox, illvar(2), 'units', 'degrees_east')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon units')
    istatus = nf90_def_var(ncidox, 'time', nf90_double, idims(3:3), ioxvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncidox, ioxvar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncidox, ioxvar(1), 'calendar', calendar)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')

    do i = 1 , noxsp
      istatus = nf90_def_var(ncidox, oxspec(i), nf90_float, x3ddim, ioxvar(i+1))
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable '//oxspec(i))
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncidox, ioxvar(i+1), 1, 1, 9)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error setting compression on '//oxspec(i))
#endif
      buffer = 'mole_concentration_of_'//trim(oxspec(i))//'_in_air'
      istatus = nf90_put_att(ncidox, ioxvar(i+1), 'standard_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding standard_name on '//trim(oxspec(i)))
      buffer = trim(oxspec(i))//' molarity'
      istatus = nf90_put_att(ncidox, ioxvar(i+1), 'long_name', buffer)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding long_name on '//trim(oxspec(i)))
      istatus = nf90_put_att(ncidox, ioxvar(i+1), 'units', 'mol cm-3')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding units on '//trim(oxspec(i)))
      istatus = nf90_put_att(ncidox, ioxvar(i+1), 'coordinates', 'xlon xlat')
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error adding coordinates on '//trim(oxspec(i)))
    end do

    istatus = nf90_enddef(ncidox)
    call checkncerr(istatus,__FILE__,__LINE__,'Error End Definitions NetCDF output')
!
    istatus = nf90_put_var(ncidox, izvar(1), sigma2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable sigma write')
    istatus = nf90_put_var(ncidox, izvar(2), hptop)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ptop write')
    istatus = nf90_put_var(ncidox, ivvar(1), yiy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable iy write')
    istatus = nf90_put_var(ncidox, ivvar(2), xjx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable jx write')
    istatus = nf90_put_var(ncidox, illvar(1), xlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlat write')
    istatus = nf90_put_var(ncidox, illvar(2), xlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlon write')

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
    istatus = nf90_put_var(ncid, ichvar(1), xdate, istart1, icount1)
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
      istatus = nf90_put_var(ncid, ichvar(i+1), chv4(:,:,:,i), istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error variable '//chspec(i)//' write')
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
    istatus = nf90_put_var(ncidox, ioxvar(1), xdate, istart1, icount1)
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
      istatus = nf90_put_var(ncidox, ioxvar(i+1), oxv4(:,:,:,i), istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error variable '//oxspec(i)//' write')
    end do

    write (stdout ,*) 'Write ch_oxcl : ', tochar(idate)

    itimeox = itimeox + 1
!
  end subroutine write_ch_oxcl
!
end module mod_wrtoxd
