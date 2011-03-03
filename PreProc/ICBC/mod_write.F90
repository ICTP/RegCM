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

      use mod_dynparam
      use m_die
      use m_zeit
      use m_mall
      use m_realkinds

      private

      integer :: ncout
      character(256) :: ofname
      integer :: irefdate
      integer :: itime
      integer , dimension(5) :: idims
      integer , dimension(8) :: ivar

      real(sp) , allocatable , dimension(:,:) :: ps4 , ts4
      real(sp) , allocatable , dimension(:,:,:) :: h4 , q4
      real(sp) , allocatable , dimension(:,:,:) :: t4 , u4 , v4
      real(sp) , allocatable , dimension(:,:,:) :: sulfate4
      real(sp) , allocatable , dimension(:) :: yiy
      real(sp) , allocatable , dimension(:) :: xjx

      public :: ps4 , ts4 , h4 , q4 , t4 , u4 , v4 , sulfate4
      public :: init_output , free_output , newfile , writef

      data ncout /-1/

      contains

      subroutine init_output
      implicit none
        integer :: ierr
        allocate(ps4(jx,iy), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate ps4',ierr)
        call mall_mci(ps4,'mod_write')
        allocate(ts4(jx,iy), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate ts4',ierr)
        call mall_mci(ts4,'mod_write')
        allocate(h4(jx,iy,kz), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate h4',ierr)
        call mall_mci(h4,'mod_write')
        allocate(q4(jx,iy,kz), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate q4',ierr)
        call mall_mci(q4,'mod_write')
        allocate(t4(jx,iy,kz), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate t4',ierr)
        call mall_mci(t4,'mod_write')
        allocate(u4(jx,iy,kz), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate u4',ierr)
        call mall_mci(u4,'mod_write')
        allocate(v4(jx,iy,kz), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate v4',ierr)
        call mall_mci(v4,'mod_write')
        if ( dattyp=='EH5OM' .and. ehso4) then
          allocate(sulfate4(jx,iy,kz), stat=ierr)
          if (ierr /= 0) call die('init_output','allocate sulfate',ierr)
          call mall_mci(sulfate4,'mod_write')
        end if
        allocate(yiy(iy), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate yiy',ierr)
        call mall_mci(yiy,'mod_write')
        allocate(xjx(jx), stat=ierr)
        if (ierr /= 0) call die('init_output','allocate xjx',ierr)
        call mall_mci(xjx,'mod_write')
      end subroutine init_output

      subroutine free_output
        use netcdf
        implicit none
        integer :: istatus
        call mall_mco(ps4,'mod_write')
        deallocate(ps4)
        call mall_mco(ts4,'mod_write')
        deallocate(ts4)
        call mall_mco(h4,'mod_write')
        deallocate(h4)
        call mall_mco(q4,'mod_write')
        deallocate(q4)
        call mall_mco(t4,'mod_write')
        deallocate(t4)
        call mall_mco(u4,'mod_write')
        deallocate(u4)
        call mall_mco(v4,'mod_write')
        deallocate(v4)
        if ( dattyp=='EH5OM' .and. ehso4) then
          call mall_mco(sulfate4,'mod_write')
          deallocate(sulfate4)
        end if
        deallocate(yiy)
        call mall_mco(yiy,'mod_write')
        deallocate(xjx)
        call mall_mco(xjx,'mod_write')
        if (ncout > 0) then
          istatus = nf90_close(ncout)
          call check_ok(istatus,('Error closing file '//trim(ofname)))
        end if
      end subroutine free_output

      subroutine newfile(idate1)
        use mod_date
        use mod_grid , only : xlat , xlon , topogm , sigma2
        use netcdf
        implicit none
        integer , intent(in) :: idate1
        integer :: istatus
        integer :: iyy , im , id , ih , i , j
        integer , dimension(8) :: tvals
        integer , dimension(2) :: izvar
        integer , dimension(2) :: ivvar
        integer , dimension(3) :: illvar
        integer , dimension(4) :: x3ddim
        character(64) :: csdate , cdum
        character(256) :: history
        real(sp) , dimension(2) :: trlat
        real(sp) :: hptop

        call zeit_ci('newfile')

        if (ncout > 0) then
          istatus = nf90_close(ncout)
          call check_ok(istatus,('Error closing file '//trim(ofname)))
        end if

        write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
              &     '_ICBC.', idate1, '.nc'

        irefdate = idate1
        itime = 1

        call split_idate(idate1,iyy,im,id,ih)
        write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
                & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

#ifdef NETCDF4_HDF5
        istatus = nf90_create(ofname, &
                  ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), &
                  ncout)
#else
        istatus = nf90_create(ofname, nf90_clobber, ncout)
#endif
        call check_ok(istatus,('Error creating file '//trim(ofname)))

        istatus = nf90_put_att(ncout, nf90_global, 'title',  &
                & 'ICTP Regional Climatic model V4 ICBC program output')
        call check_ok(istatus,'Error adding global title')
        istatus = nf90_put_att(ncout, nf90_global, 'institution', &
                 & 'ICTP')
        call check_ok(istatus,'Error adding global institution')
        istatus = nf90_put_att(ncout, nf90_global, 'source', &
                 & 'RegCM Model simulation SST output')
        call check_ok(istatus,'Error adding global source')
        istatus = nf90_put_att(ncout, nf90_global, 'Conventions', &
                 & 'CF-1.4')
        call check_ok(istatus,'Error adding global Conventions')
        call date_and_time(values=tvals)
        write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
             ' : Created by RegCM icbc program'
        istatus = nf90_put_att(ncout, nf90_global, 'history', history)
        call check_ok(istatus,'Error adding global history')
        istatus = nf90_put_att(ncout, nf90_global, 'references', &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        call check_ok(istatus,'Error adding global references')
        istatus = nf90_put_att(ncout, nf90_global, 'experiment', &
                 & domname)
        call check_ok(istatus,'Error adding global experiment')
        istatus = nf90_put_att(ncout, nf90_global, 'projection', iproj)
        call check_ok(istatus,'Error adding global projection')
        istatus = nf90_put_att(ncout, nf90_global,   &
                 &   'grid_size_in_meters', ds*1000.0)
        call check_ok(istatus,'Error adding global gridsize')
        istatus = nf90_put_att(ncout, nf90_global,   &
                 &   'latitude_of_projection_origin', clat)
        call check_ok(istatus,'Error adding global clat')
        istatus = nf90_put_att(ncout, nf90_global,   &
                 &   'longitude_of_projection_origin', clon)
        call check_ok(istatus,'Error adding global clon')
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(ncout, nf90_global, &
                   &   'grid_north_pole_latitude', plat)
          call check_ok(istatus,'Error adding global plat')
          istatus = nf90_put_att(ncout, nf90_global, &
                   &   'grid_north_pole_longitude', plon)
          call check_ok(istatus,'Error adding global plon')
        else if (iproj == 'LAMCON') then
          trlat(1) = truelatl
          trlat(2) = truelath
          istatus = nf90_put_att(ncout, nf90_global, &
                   &   'standard_parallel', trlat)
          call check_ok(istatus,'Error adding global truelat')
        end if
        istatus = nf90_put_att(ncout, nf90_global,  &
                           &   'global_data_source', dattyp)
        call check_ok(istatus,'Error adding global data_source')
        if (dattyp == 'EH5OM') then
          if (ehso4) then
            cdum = 'Yes'
          else
            cdum = 'No'
          end if
          istatus = nf90_put_att(ncout, nf90_global,  &
                             &   'sulfate_data_present', cdum)
          call check_ok(istatus,'Error adding global sulfate_present')
        end if
        istatus = nf90_def_dim(ncout, 'iy', iy, idims(2))
        call check_ok(istatus,'Error creating dimension iy')
        istatus = nf90_def_dim(ncout, 'jx', jx, idims(1))
        call check_ok(istatus,'Error creating dimension jx')
        istatus = nf90_def_dim(ncout, 'time', nf90_unlimited, idims(3))
        call check_ok(istatus,'Error creating dimension time')
        istatus = nf90_def_dim(ncout, 'kz', kz, idims(4))
        call check_ok(istatus,'Error creating dimension kz')
        x3ddim(1) = idims(1)
        x3ddim(2) = idims(2)
        x3ddim(3) = idims(4)
        x3ddim(4) = idims(3)
!
        istatus = nf90_def_var(ncout, 'sigma', nf90_float, idims(4),   &
                            &  izvar(1))
        call check_ok(istatus,'Error adding variable sigma')
        istatus = nf90_put_att(ncout, izvar(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')      
        call check_ok(istatus,'Error adding sigma standard_name')
        istatus = nf90_put_att(ncout, izvar(1), 'long_name',      &
                            &  'Sigma at model layers')
        call check_ok(istatus,'Error adding sigma long_name')
        istatus = nf90_put_att(ncout, izvar(1), 'units', '1')
        call check_ok(istatus,'Error adding sigma units')
        istatus = nf90_put_att(ncout, izvar(1), 'axis', 'Z')
        call check_ok(istatus,'Error adding sigma axis')
        istatus = nf90_put_att(ncout, izvar(1), 'positive', 'down')
        call check_ok(istatus,'Error adding sigma positive')
        istatus = nf90_put_att(ncout, izvar(1), 'formula_terms',  &
                     &         'sigma: sigma ps: ps ptop: ptop')
        call check_ok(istatus,'Error adding sigma formula_terms')
        istatus = nf90_def_var(ncout, 'ptop', nf90_float,         &
                           &   varid=izvar(2))
        call check_ok(istatus,'Error adding variable ptop')
        istatus = nf90_put_att(ncout, izvar(2), 'standard_name',  &
                            &  'air_pressure')
        call check_ok(istatus,'Error adding ptop standard_name')
        istatus = nf90_put_att(ncout, izvar(2), 'long_name',      &
                            &  'Pressure at model top')
        call check_ok(istatus,'Error adding ptop long_name')
        istatus = nf90_put_att(ncout, izvar(2), 'units', 'hPa')
        call check_ok(istatus,'Error adding ptop units')
        istatus = nf90_def_var(ncout, 'iy', nf90_float, idims(2), &
                            &  ivvar(1))
        call check_ok(istatus,'Error adding variable iy')
        istatus = nf90_put_att(ncout, ivvar(1), 'standard_name',  &
                            &  'projection_y_coordinate')
        call check_ok(istatus,'Error adding iy standard_name')
        istatus = nf90_put_att(ncout, ivvar(1), 'long_name',      &
                            &  'y-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding iy long_name')
        istatus = nf90_put_att(ncout, ivvar(1), 'units', 'km')
        call check_ok(istatus,'Error adding iy units')
        istatus = nf90_def_var(ncout, 'jx', nf90_float, idims(1), &
                            &  ivvar(2))
        call check_ok(istatus,'Error adding variable jx')
        istatus = nf90_put_att(ncout, ivvar(2), 'standard_name', &
                            &  'projection_x_coordinate')
        call check_ok(istatus,'Error adding jx standard_name')
        istatus = nf90_put_att(ncout, ivvar(2), 'long_name',    &
                            &  'x-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding jx long_name')
        istatus = nf90_put_att(ncout, ivvar(2), 'units', 'km')
        call check_ok(istatus,'Error adding jx units')
        istatus = nf90_def_var(ncout, 'xlat', nf90_float, idims(1:2),  &
                            &  illvar(1))
        call check_ok(istatus,'Error adding variable xlat')
        istatus = nf90_put_att(ncout, illvar(1), 'standard_name', &
                            &  'latitude')
        call check_ok(istatus,'Error adding xlat standard_name')
        istatus = nf90_put_att(ncout, illvar(1), 'long_name',     &
                            &  'Latitude at cross points')
        call check_ok(istatus,'Error adding xlat long_name')
        istatus = nf90_put_att(ncout, illvar(1), 'units',         &
                            &  'degrees_north')
        call check_ok(istatus,'Error adding xlat units')
        istatus = nf90_def_var(ncout, 'xlon', nf90_float, idims(1:2),  &
                            &  illvar(2))
        call check_ok(istatus,'Error adding variable xlon')
        istatus = nf90_put_att(ncout, illvar(2), 'standard_name', &
                            &  'longitude')
        call check_ok(istatus,'Error adding xlon standard_name')
        istatus = nf90_put_att(ncout, illvar(2), 'long_name',     &
                            &  'Longitude at cross points')
        call check_ok(istatus,'Error adding xlon long_name')
        istatus = nf90_put_att(ncout, illvar(2), 'units',         &
                            &  'degrees_east')
        call check_ok(istatus,'Error adding xlon units')
        istatus = nf90_def_var(ncout, 'topo', nf90_float, idims(1:2),  &
                            &  illvar(3))
        call check_ok(istatus,'Error adding variable topo')
        istatus = nf90_put_att(ncout, illvar(3), 'standard_name', &
                            &  'surface_altitude')
        call check_ok(istatus,'Error adding topo standard_name')
        istatus = nf90_put_att(ncout, illvar(3), 'long_name',     &
                            &  'Domain surface elevation')
        call check_ok(istatus,'Error adding topo long_name')
        istatus = nf90_put_att(ncout, illvar(3), 'units',         &
                            &  'm')
        call check_ok(istatus,'Error adding topo units')
        istatus = nf90_put_att(ncout, illvar(3), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding topo coordinates')
        istatus = nf90_def_var(ncout, 'time', nf90_double, idims(3:3),  &
                            &  ivar(1))
        call check_ok(istatus,'Error adding variable time')
        istatus = nf90_put_att(ncout, ivar(1), 'units', &
                       &   'hours since '//csdate)
        call check_ok(istatus,'Error adding time units')
        istatus = nf90_def_var(ncout, 'ps', nf90_float, idims(1:3),  &
                            &  ivar(2))
        call check_ok(istatus,'Error adding variable ps')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(2), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on ps')
#endif
        istatus = nf90_put_att(ncout, ivar(2), 'standard_name', &
                            &  'surface_air_pressure')
        call check_ok(istatus,'Error adding ps standard_name')
        istatus = nf90_put_att(ncout, ivar(2), 'long_name',     &
                            &  'Surface pressure')
        call check_ok(istatus,'Error adding ps long_name')
        istatus = nf90_put_att(ncout, ivar(2), 'units', 'hPa')
        call check_ok(istatus,'Error adding ps units')
        istatus = nf90_put_att(ncout, ivar(2), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding ps coordinates')
        istatus = nf90_def_var(ncout, 'ts', nf90_float, idims(1:3),  &
                            &  ivar(3))
        call check_ok(istatus,'Error adding variable ts')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(3), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on ts')
#endif
        istatus = nf90_put_att(ncout, ivar(3), 'standard_name', &
                            &  'surface_temperature')
        call check_ok(istatus,'Error adding ts standard_name')
        istatus = nf90_put_att(ncout, ivar(3), 'long_name',     &
                            &  'Surface Temperature')
        call check_ok(istatus,'Error adding ts long_name')
        istatus = nf90_put_att(ncout, ivar(3), 'units', 'K')
        call check_ok(istatus,'Error adding ts units')
        istatus = nf90_put_att(ncout, ivar(3), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding ts coordinates')
        istatus = nf90_def_var(ncout, 'u', nf90_float, x3ddim,  &
                            &  ivar(4))
        call check_ok(istatus,'Error adding variable u')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(4), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on u')
#endif
        istatus = nf90_put_att(ncout, ivar(4), 'standard_name', &
                            &  'eastward_wind')
        call check_ok(istatus,'Error adding u standard_name')
        istatus = nf90_put_att(ncout, ivar(4), 'long_name',     &
                            &  'U component (westerly) of wind')
        call check_ok(istatus,'Error adding u long_name')
        istatus = nf90_put_att(ncout, ivar(4), 'units', 'm s-1')
        call check_ok(istatus,'Error adding u units')
        istatus = nf90_put_att(ncout, ivar(4), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding u coordinates')
        istatus = nf90_def_var(ncout, 'v', nf90_float, x3ddim,  &
                            &  ivar(5))
        call check_ok(istatus,'Error adding variable v')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(5), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on v')
#endif
        istatus = nf90_put_att(ncout, ivar(5), 'standard_name', &
                            &  'northward_wind')
        call check_ok(istatus,'Error adding v standard_name')
        istatus = nf90_put_att(ncout, ivar(5), 'long_name',     &
                            &  'V component (southerly) of wind')
        call check_ok(istatus,'Error adding v long_name')
        istatus = nf90_put_att(ncout, ivar(5), 'units', 'm s-1')
        call check_ok(istatus,'Error adding v units')
        istatus = nf90_put_att(ncout, ivar(5), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding v coordinates')
        istatus = nf90_def_var(ncout, 't', nf90_float, x3ddim,  &
                            &  ivar(6))
        call check_ok(istatus,'Error adding variable t')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(6), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on t')
#endif
        istatus = nf90_put_att(ncout, ivar(6), 'standard_name', &
                            &  'air_temperature')
        call check_ok(istatus,'Error adding t standard_name')
        istatus = nf90_put_att(ncout, ivar(6), 'long_name',     &
                            &  'Temperature')
        call check_ok(istatus,'Error adding t long_name')
        istatus = nf90_put_att(ncout, ivar(6), 'units', 'K')
        call check_ok(istatus,'Error adding t units')
        istatus = nf90_put_att(ncout, ivar(6), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding t coordinates')
        istatus = nf90_def_var(ncout, 'qv', nf90_float, x3ddim,  &
                            &  ivar(7))
        call check_ok(istatus,'Error adding variable qv')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncout, ivar(7), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on qv')
#endif
        istatus = nf90_put_att(ncout, ivar(7), 'standard_name', &
                            &  'humidity_mixing_ratio')
        call check_ok(istatus,'Error adding qv standard_name')
        istatus = nf90_put_att(ncout, ivar(7), 'long_name',     &
                            &  'Water vapor mixing ratio')
        call check_ok(istatus,'Error adding qv long_name')
        istatus = nf90_put_att(ncout, ivar(7), 'units', 'kg kg-1')
        call check_ok(istatus,'Error adding qv units')
        istatus = nf90_put_att(ncout, ivar(7), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding qv coordinates')
        if ( dattyp=='EH5OM' .and. ehso4) then
          istatus = nf90_def_var(ncout, 'so4', nf90_float, x3ddim, &
                              &  ivar(8))
          call check_ok(istatus,'Error adding variable so4')
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var_deflate(ncout, ivar(8), 1, 1, 9)
          call check_ok(istatus,'Error setting compression on qv')
#endif
          istatus = nf90_put_att(ncout, ivar(8), 'standard_name', &
                              &  'atmosphere_sulfate_content')
          call check_ok(istatus,'Error adding so4 standard_name')
          istatus = nf90_put_att(ncout, ivar(8), 'long_name',     &
                              &  'Sulfate')
          call check_ok(istatus,'Error adding so4 long_name')
          istatus = nf90_put_att(ncout, ivar(8), 'units', 'kg m-2')
          call check_ok(istatus,'Error adding so4 units')
          istatus = nf90_put_att(ncout, ivar(8), 'coordinates', &
                              &  'xlon xlat')
          call check_ok(istatus,'Error adding so4 coordinates')
        end if
!
        istatus = nf90_enddef(ncout)
        call check_ok(istatus,'Error End Definitions NetCDF output')
!
        istatus = nf90_put_var(ncout, izvar(1), sigma2)
        call check_ok(istatus,'Error variable sigma write')
        hptop = ptop*10.0
        istatus = nf90_put_var(ncout, izvar(2), hptop)
        call check_ok(istatus,'Error variable ptop write')
        yiy(1) = -(dble(iy-1)/2.0) * ds
        xjx(1) = -(dble(jx-1)/2.0) * ds
        do i = 2 , iy
          yiy(i) = yiy(i-1)+ds
        end do
        do j = 2 , jx
          xjx(j) = xjx(j-1)+ds
        end do
        istatus = nf90_put_var(ncout, ivvar(1), yiy)
        call check_ok(istatus,'Error variable iy write')
        istatus = nf90_put_var(ncout, ivvar(2), xjx)
        call check_ok(istatus,'Error variable jx write')
        istatus = nf90_put_var(ncout, illvar(1), xlat)
        call check_ok(istatus,'Error variable xlat write')
        istatus = nf90_put_var(ncout, illvar(2), xlon)
        call check_ok(istatus,'Error variable xlon write')
        istatus = nf90_put_var(ncout, illvar(3), topogm)
        call check_ok(istatus,'Error variable xlon write')

        istatus = nf90_sync(ncout)
        call check_ok(istatus,'Error file sync')

        call zeit_co('newfile')

99001 format (a,a,a,a,i10,a)

      end subroutine newfile

      subroutine writef(idate)
        use netcdf
        use mod_date
        implicit none
        integer , intent(in) :: idate
        integer :: istatus
        integer , dimension(1) :: istart1 , icount1
        integer , dimension(4) :: istart , icount
        real(dp) , dimension(1) :: xdate
!
        call zeit_ci('writef')
!
        istart1(1) = itime
        icount1(1) = 1
        xdate(1) = dble(idatediff(idate,irefdate))
        istatus = nf90_put_var(ncout, ivar(1), xdate, istart1, icount1)
        call check_ok(istatus,'Error variable time write')
        istart(3) = itime
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = iy
        icount(1) = jx
        ps4 = (ps4+ptop)*10.0
        istatus = nf90_put_var(ncout, ivar(2), ps4, istart(1:3), &
                               icount(1:3))
        call check_ok(istatus,'Error variable ps write')
        istatus = nf90_put_var(ncout, ivar(3), ts4, istart(1:3), &
                               icount(1:3))
        call check_ok(istatus,'Error variable ts write')
        istart(4) = itime
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = kz
        icount(2) = iy
        icount(1) = jx
        istatus = nf90_put_var(ncout, ivar(4), u4, istart, icount)
        call check_ok(istatus,'Error variable u write')
        istatus = nf90_put_var(ncout, ivar(5), v4, istart, icount)
        call check_ok(istatus,'Error variable v write')
        istatus = nf90_put_var(ncout, ivar(6), t4, istart, icount)
        call check_ok(istatus,'Error variable t write')
        istatus = nf90_put_var(ncout, ivar(7), q4, istart, icount)
        call check_ok(istatus,'Error variable qv write')
        if ( dattyp=='EH5OM' .and. ehso4) then
          istatus = nf90_put_var(ncout, ivar(8), sulfate4,  &
                            &    istart, icount)
          call check_ok(istatus,'Error variable so4 write')
        end if
        itime = itime + 1
!
        call zeit_co('writef')
      end subroutine writef
!
      subroutine check_ok(ierr,message)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        character(*) :: message
        if (ierr /= nf90_noerr) then
          call die('writef',message,1,nf90_strerror(ierr),ierr)
        end if
      end subroutine check_ok
!
      end module mod_write
