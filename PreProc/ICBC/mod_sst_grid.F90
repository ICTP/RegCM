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

  use mod_realkinds
  use mod_stdio
  use netcdf
  use mod_dynparam
  use mod_memutil
  use mod_message

  private

  integer :: ncid
  integer , dimension(4) :: idims
  integer , dimension(4) :: ivar
  type (rcm_time_and_date) , save :: refdate
  integer :: itime

  real(sp) , public , pointer , dimension(:,:) :: lu , sstmm , icemm ,   &
                                            xlat , xlon , finmat
  real(sp) , pointer , dimension(:) :: sigma
  real(sp) , pointer , dimension(:) :: yiy
  real(sp) , pointer , dimension(:) :: xjx

  public :: init_grid , free_grid , read_domain , open_sstfile , &
            close_sstfile , writerec

  contains

  subroutine init_grid
    implicit none
    call getmem2d(lu,1,iy,1,jx,'mod_sst_grid:lu')
    call getmem2d(sstmm,1,iy,1,jx,'mod_sst_grid:sstmm')
    call getmem2d(icemm,1,iy,1,jx,'mod_sst_grid:icemm')
    call getmem2d(xlat,1,iy,1,jx,'mod_sst_grid:xlat')
    call getmem2d(xlon,1,iy,1,jx,'mod_sst_grid:xlon')
    call getmem2d(finmat,1,jx,1,iy,'mod_sst_grid:finmat')
    call getmem1d(sigma,1,kzp1,'mod_sst_grid:sigma')
    call getmem1d(yiy,1,iy,'mod_sst_grid:yiy')
    call getmem1d(xjx,1,jx,'mod_sst_grid:xjx')
  end subroutine init_grid

  subroutine read_domain(terfile)
    implicit none
    character(256) :: terfile
    intent(in) :: terfile

    integer :: iyy , jxx
    integer :: istatus , incin , idimid , ivarid

    istatus = nf90_open(terfile, nf90_nowrite, incin)
    call checkncerr(istatus,__FILE__,__LINE__,'Error Opening file '//trim(terfile))

    istatus = nf90_inq_dimid(incin, "iy", idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension iy missing')
    istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension iy read error')
    istatus = nf90_inq_dimid(incin, "jx", idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension jx missing')
    istatus = nf90_inquire_dimension(incin, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension jx read error')
    if ( iyy /= iy .or. jxx /= jx ) then
      write (stderr,*) 'IMPROPER DIMENSION SPECIFICATION'
      write (stderr,*) '  namelist   : ' , iy , jx
      write (stderr,*) '  DOMAIN     : ' , iyy , jxx
      call die('read_domain','Dimensions mismatch',1)
    end if

    istatus = nf90_inq_varid(incin, "sigma", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable sigma missing')
    istatus = nf90_get_var(incin, ivarid, sigma)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable sigma read error')
    istatus = nf90_inq_varid(incin, "landuse", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable landuse missing')
    istatus = nf90_get_var(incin, ivarid, finmat)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable landuse read error')
    lu = transpose(finmat)
    istatus = nf90_inq_varid(incin, "xlat", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable xlat missing')
    istatus = nf90_get_var(incin, ivarid, finmat)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable xlat read error')
    xlat = transpose(finmat)
    istatus = nf90_inq_varid(incin, "xlon", ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable xlon missing')
    istatus = nf90_get_var(incin, ivarid, finmat)
    call checkncerr(istatus,__FILE__,__LINE__,'Variable xlon read error')
    xlon = transpose(finmat)
    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__,('Error closing file '//trim(terfile)))

  end subroutine read_domain

  subroutine open_sstfile(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    integer :: istatus
    character(256) :: sstname , history
    character(64) :: csdate
    real(sp) , dimension(2) :: trlat
    real(sp) :: hptop
    integer , dimension(2) :: ivvar
    integer , dimension(2) :: illvar
    integer , dimension(2) :: izvar
    integer , dimension(8) :: tvals
    integer :: i , j

    refdate = idate1
    csdate = tochar(refdate)
    itime = 1

    sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
#ifdef NETCDF4_HDF5
    istatus = nf90_create(sstname, &
              ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid)
#else
    istatus = nf90_create(sstname, nf90_clobber, ncid)
#endif
    call checkncerr(istatus,__FILE__,__LINE__,('Error creating file '//trim(sstname)))
    istatus = nf90_put_att(ncid, nf90_global, 'title',  &
               'ICTP Regional Climatic model V4 SST program output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global title')
    istatus = nf90_put_att(ncid, nf90_global, 'institution','ICTP')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global institution')
    istatus = nf90_put_att(ncid, nf90_global, 'source', &
               'RegCM Model simulation SST output')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global source')
    istatus = nf90_put_att(ncid, nf90_global, 'Conventions','CF-1.4')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global Conventions')
    call date_and_time(values = tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
         ' : Created by RegCM sst program'
    istatus = nf90_put_att(ncid, nf90_global, 'history', history)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global history')

    istatus = nf90_put_att(ncid, nf90_global, 'references', &
               'http://eforge.escience-lab.org/gf/project/regcm')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global references')
    istatus = nf90_put_att(ncid, nf90_global, 'experiment', domname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global experiment')
    istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global projection')
    istatus = nf90_put_att(ncid, nf90_global,   &
                 'grid_size_in_meters', ds*1000.0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global gridsize')
    istatus = nf90_put_att(ncid, nf90_global,   &
                 'latitude_of_projection_origin', clat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clat')
    istatus = nf90_put_att(ncid, nf90_global,   &
                 'longitude_of_projection_origin', clon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global clon')
    if (iproj == 'ROTMER') then
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_north_pole_latitude', plat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plat')
      istatus = nf90_put_att(ncid, nf90_global, &
                   'grid_north_pole_longitude', plon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global plon')
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      istatus = nf90_put_att(ncid, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error adding global truelat')
    end if
    istatus = nf90_put_att(ncid, nf90_global, 'sst_source', ssttyp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding global sst_source')
!
    istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension iy')
    istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension jx')
    istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension time')
    istatus = nf90_def_dim(ncid, 'kz', kz+1, idims(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error creating dimension kz')
!
    istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4), izvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sigma')
    istatus = nf90_put_att(ncid, izvar(1), 'standard_name',  &
                           'atmosphere_sigma_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding sigma standard_name')
    istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                           'Sigma at model layer midpoints')
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
    istatus = nf90_put_att(ncid, izvar(2), 'long_name',      &
                           'Pressure at model top')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop long_name')
    istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding ptop units')
    istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ivvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable iy')
    istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                           'projection_y_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy standard_name')
    istatus = nf90_put_att(ncid, ivvar(1), 'long_name',      &
                           'y-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy long_name')
    istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding iy uits')
    istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ivvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable jx')
    istatus = nf90_put_att(ncid, ivvar(2), 'standard_name',  &
                           'projection_x_coordinate')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx standard_name')
    istatus = nf90_put_att(ncid, ivvar(2), 'long_name',      &
                           'x-coordinate in Cartesian system')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx long_name')
    istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding jx units')
    istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2), illvar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlat')
    istatus = nf90_put_att(ncid, illvar(1), 'standard_name', 'latitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat standard_name')
    istatus = nf90_put_att(ncid, illvar(1), 'long_name',     &
                           'Latitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat long_name')
    istatus = nf90_put_att(ncid, illvar(1), 'units', 'degrees_north')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlat units')
    istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2), illvar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable xlon')
    istatus = nf90_put_att(ncid, illvar(2), 'standard_name', 'longitude')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon standard_name')
    istatus = nf90_put_att(ncid, illvar(2), 'long_name',     &
                           'Longitude at cross points')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon long_name')
    istatus = nf90_put_att(ncid, illvar(2), 'units', 'degrees_east')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding xlon units')
    istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), ivar(1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable time')
    istatus = nf90_put_att(ncid, ivar(1), 'units', 'hours since '//csdate)
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time units')
    istatus = nf90_put_att(ncid, ivar(1), 'calendar', calstr(refdate%calendar))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding time calendar')
    istatus = nf90_def_var(ncid, 'sst', nf90_float, idims(1:3), ivar(2))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable sst')
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
      istatus = nf90_def_var(ncid, 'ice', nf90_float, idims(1:3), ivar(3))
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
    istatus = nf90_def_var(ncid, 'landuse', nf90_float,idims(1:3), ivar(4))
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding variable landuse')
    istatus = nf90_put_att(ncid, ivar(4), 'legend',               &
              '1  => Crop/mixed farming'//char(10)//              &
              '2  => Short grass'//char(10)//                     &
              '3  => Evergreen needleleaf tree'//char(10)//       &
              '4  => Deciduous needleleaf tree'//char(10)//       &
              '5  => Deciduous broadleaf tree'//char(10)//        &
              '6  => Evergreen broadleaf tree'//char(10)//        &
              '7  => Tall grass'//char(10)//                      &
              '8  => Desert'//char(10)//                          &
              '9  => Tundra'//char(10)//                          &
              '10 => Irrigated Crop'//char(10)//                  &
              '11 => Semi-desert'//char(10)//                     &
              '12 => Ice cap/glacier'//char(10)//                 &
              '13 => Bog or marsh'//char(10)//                    &
              '14 => Inland water'//char(10)//                    &
              '15 => Ocean'//char(10)//                           &
              '16 => Evergreen shrub'//char(10)//                 &
              '17 => Deciduous shrub'//char(10)//                 &
              '18 => Mixed Woodland'//char(10)//                  &
              '19 => Forest/Field mosaic'//char(10)//             &
              '20 => Water and Land mixture'//char(10)//          &
              '21 => Urban'//char(10)//                           &
              '22 => Sub-Urban')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding landuse legend')
    istatus = nf90_put_att(ncid, ivar(4), 'standard_name', 'land_type')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding landuse standard_name')
    istatus = nf90_put_att(ncid, ivar(4), 'long_name',            &
                     'Landuse category as defined in BATS1E')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding landuse long_name')
    istatus = nf90_put_att(ncid, ivar(4), 'units', '1')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding landuse units')
    istatus = nf90_put_att(ncid, ivar(4), 'coordinates', 'xlon xlat')
    call checkncerr(istatus,__FILE__,__LINE__,'Error adding landuse coordinates')
!
    istatus = nf90_enddef(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error End Definitions NetCDF output')
!
    istatus = nf90_put_var(ncid, izvar(1), sigma)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable sigma write')
    hptop = real(ptop*10.0D0)
    istatus = nf90_put_var(ncid, izvar(2), hptop)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable ptop write')
    yiy(1) = -real((dble(iy-1)/2.0D0) * ds)
    xjx(1) = -real((dble(jx-1)/2.0D0) * ds)
    do i = 2 , iy
      yiy(i) = real(dble(yiy(i-1))+ds)
    end do
    do j = 2 , jx
      xjx(j) = real(dble(xjx(j-1))+ds)
    end do
    istatus = nf90_put_var(ncid, ivvar(1), yiy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable iy write')
    istatus = nf90_put_var(ncid, ivvar(2), xjx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable jx write')
    istatus = nf90_put_var(ncid, illvar(1), transpose(xlat))
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlat write')
    istatus = nf90_put_var(ncid, illvar(2), transpose(xlon))
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable xlon write')

  end subroutine open_sstfile

  subroutine close_sstfile
    implicit none
    integer :: istatus
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error closing output file')
  end subroutine close_sstfile

  subroutine writerec(idate,lice)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: lice
    integer :: istatus
    integer , dimension(1) :: istart1 , icount1
    integer , dimension(3) :: istart , icount
    real(dp) , dimension(1) :: xdate
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
    istatus = nf90_put_var(ncid, ivar(4), transpose(lu), istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error variable landuse write')
    if ( debug_level > 2 ) then
      istatus = nf90_sync(ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error sync output file')
    end if
    itime = itime + 1
  end subroutine writerec

end module mod_sst_grid
