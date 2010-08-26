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

      contains

      subroutine setup(iy,jx,ntypec,iproj,ds,clat,clon)
      use mod_block
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds
      integer :: iy , jx , ntypec
      character(6) :: iproj
      intent (in) clat , clon , ds , iproj , iy , jx , ntypec
!
      write (6,*) ' '
      write (6,*) 'Doing Domain Setup with following parameters'
      write (6,*) ' '
      write (6,*) 'ntypec = ' , ntypec
      write (6,*) 'iy     = ' , iy
      write (6,*) 'jx     = ' , jx
      write (6,*) 'ds     = ' , ds
      write (6,*) 'clat   = ' , clat
      write (6,*) 'clon   = ' , clon
      write (6,*) 'rin    = ' , rin
      write (6,*) 'iproj  = ' , iproj
      write (6,*) ' '
!
      rin = 1.5          ! 1.5 rad of influence-coarse mesh
!
      dsinm = ds*1000.
!
      nnc = nint(60./float(ntypec))
      xnc = float(ntypec)/60.
      print * , '***** Terrain resolution (min): ' , xnc*60.
!
      end subroutine setup
!
      subroutine write_domain(lsub)
        use netcdf
        use mod_dynparam
        use mod_block
        use mod_maps
        implicit none
        logical , intent (in) :: lsub

        integer :: istatus , i , j
        integer :: incout
        integer , dimension(4) :: idims
        integer , dimension(3) :: istart
        integer , dimension(3) :: icount
        integer , dimension(12) :: ivar
        integer , dimension(2) :: itvar
        integer , dimension(2) :: ivdim
        integer , dimension(2) :: izdim
        integer , dimension(8) :: tvals
        character(256) :: fname , history
        character(3) :: cnsg
        real(4) , dimension(2) :: trlat
        real(4) :: hptop
        real(4) , allocatable , dimension(:) :: yiy
        real(4) , allocatable , dimension(:) :: xjx

        trlat(1) = truelatl
        trlat(2) = truelath

        if (lsub) then
          allocate(yiy(iysg))
          allocate(xjx(jxsg))
          yiy(1) = -(dble(iysg-1)/2.0) * ds
          xjx(1) = -(dble(jxsg-1)/2.0) * ds
          do i = 2 , iysg
            yiy(i) = yiy(i-1)+ds
          end do
          do j = 2 , jxsg
            xjx(j) = xjx(j-1)+ds
          end do
        else
          allocate(yiy(iy))
          allocate(xjx(jx))
          yiy(1) = -(dble(iy-1)/2.0) * ds
          xjx(1) = -(dble(jx-1)/2.0) * ds
          do i = 2 , iy
            yiy(i) = yiy(i-1)+ds
          end do
          do j = 2 , jx
            xjx(j) = xjx(j-1)+ds
          end do
        end if

        if (lsub) then
          write (cnsg, '(i0.3)') nsg
          fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//      &
                    & cnsg//'.nc'
        else
          fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_create(fname, ior(nf90_clobber,nf90_hdf5),       &
                         &    incout)
#else
        istatus = nf90_create(fname, nf90_clobber, incout)
#endif
        call check_ok(istatus, &
                & 'Error creating NetCDF output '//trim(fname))

        istatus = nf90_put_att(incout, nf90_global, 'title',            &
                 & 'ICTP Regional Climatic model V4 domain')
        call check_ok(istatus,'Error adding global title')
        istatus = nf90_put_att(incout, nf90_global, 'institution',      &
                 & 'ICTP')
        call check_ok(istatus,'Error adding global institution')
        istatus = nf90_put_att(incout, nf90_global, 'source',           &
                 & 'RegCM Model simulation Terrain output')
        call check_ok(istatus,'Error adding global source')
        istatus = nf90_put_att(incout, nf90_global, 'Conventions',      &
                 & 'CF-1.4')
        call check_ok(istatus,'Error adding global Conventions')
        call date_and_time(values=tvals)
        write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
             ' : Created by RegCM terrain'
        istatus = nf90_put_att(incout, nf90_global, 'history', history)
        call check_ok(istatus,'Error adding global history')
        istatus = nf90_put_att(incout, nf90_global, 'references',       &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        call check_ok(istatus,'Error adding global references')
        istatus = nf90_put_att(incout, nf90_global, 'experiment',       &
                 & domname)
        call check_ok(istatus,'Error adding global experiment')
        istatus = nf90_put_att(incout, nf90_global, 'projection', iproj)
        call check_ok(istatus,'Error adding global projection')
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'grid_size_in_meters', ds*1000.0)
        call check_ok(istatus,'Error adding global gridsize')
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'latitude_of_projection_origin', clat)
        call check_ok(istatus,'Error adding global clat')
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'longitude_of_projection_origin', clon)
        call check_ok(istatus,'Error adding global clon')
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'latitude_of_projection_pole', plat)
          call check_ok(istatus,'Error adding global plat')
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'longitude_of_projection_pole', plon)
          call check_ok(istatus,'Error adding global plon')
        else if (iproj == 'LAMCON') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'standard_parallel', trlat)
          call check_ok(istatus,'Error adding global truelat')
        end if
        if (ifanal) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Cressman type objective analysis')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Overlapping parabolic 16 points')
        end if
        call check_ok(istatus,'Error adding global data_interpolation')
        if (smthbdy) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'No')
        end if
        call check_ok(istatus,'Error adding global boundary_smoothing')
        if (lakadj) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'No')
        end if
        call check_ok(istatus,'Error adding global great_lakes_adjust')
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'minimum_h2o_pct_for_water', h2opct)
        call check_ok(istatus,'Error adding global min_h2o_pct_for_wat')

        if (lsub) then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec_s)
          call check_ok(istatus,'Error adding global ntypec_s')
          if (fudge_lnd_s) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          call check_ok(istatus,'Error adding global landuse_fudging_s')
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex_s) then
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'No')
            end if
            call check_ok(istatus,'Error adding global texture_fudge_s')
          end if
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec)
          call check_ok(istatus,'Error adding global ntypec')
          if (fudge_lnd) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          call check_ok(istatus,'Error adding global landuse_fudging')
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex) then
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'No')
            end if
            call check_ok(istatus,'Error adding global texture_fudging')
          end if
        end if

        istatus = nf90_put_att(incout, nf90_global, 'grid_factor', xn)
        call check_ok(istatus,'Error adding global grid_factor')

        if (lsub) then
          istatus = nf90_def_dim(incout, 'iy', iysg, idims(2))
          call check_ok(istatus,'Error adding dimension iy')
          istatus = nf90_def_dim(incout, 'jx', jxsg, idims(1))
          call check_ok(istatus,'Error adding dimension jx')
        else
          istatus = nf90_def_dim(incout, 'iy', iy, idims(2))
          call check_ok(istatus,'Error adding dimension iy')
          istatus = nf90_def_dim(incout, 'jx', jx, idims(1))
          call check_ok(istatus,'Error adding dimension jx')
        end if
        if ( aertyp(7:7)=='1' ) then
          istatus = nf90_def_dim(incout, 'ntex', ntex, idims(3))
          call check_ok(istatus,'Error adding dimension NVEG')
        end if
        istatus = nf90_def_dim(incout, 'kz', kz+1, idims(4))
        call check_ok(istatus,'Error adding dimension kz')

        istatus = nf90_def_var(incout, 'sigma', nf90_float, idims(4),   &
                            &  izdim(1))
        call check_ok(istatus,'Error adding variable sigma')
        istatus = nf90_put_att(incout, izdim(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')
        call check_ok(istatus,'Error adding sigma standard_name')
        istatus = nf90_put_att(incout, izdim(1), 'long_name',           &
                            &  'Sigma at model layer midpoints')
        call check_ok(istatus,'Error adding sigma long_name')
        istatus = nf90_put_att(incout, izdim(1), 'units', '1')
        call check_ok(istatus,'Error adding sigma units')
        istatus = nf90_put_att(incout, izdim(1), 'axis', 'Z')
        call check_ok(istatus,'Error adding sigma axis')
        istatus = nf90_put_att(incout, izdim(1), 'positive', 'down')
        call check_ok(istatus,'Error adding sigma positive')
        istatus = nf90_put_att(incout, izdim(1), 'formula_terms',       &
                     &         'sigma: sigma ps: ps ptop: ptop')
        call check_ok(istatus,'Error adding sigma formula_terms')
        istatus = nf90_def_var(incout, 'ptop', nf90_float,              &
                           &   varid=izdim(2))
        call check_ok(istatus,'Error adding variable ptop')
        istatus = nf90_put_att(incout, izdim(2), 'standard_name',       &
                            &  'air_pressure')
        call check_ok(istatus,'Error adding ptop standard_name')
        istatus = nf90_put_att(incout, izdim(2), 'long_name',           &
                            &  'Pressure at model top')
        call check_ok(istatus,'Error adding ptop long_name')
        istatus = nf90_put_att(incout, izdim(2), 'units', 'hPa')
        call check_ok(istatus,'Error adding ptop units')
        istatus = nf90_def_var(incout, 'iy', nf90_float, idims(2),      &
                            &  ivdim(1))
        call check_ok(istatus,'Error adding variable iy')
        istatus = nf90_put_att(incout, ivdim(1), 'standard_name',       &
                            &  'projection_y_coordinate')
        call check_ok(istatus,'Error adding iy standard_name')
        istatus = nf90_put_att(incout, ivdim(1), 'long_name',           &
                            &  'y-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding iy long_name')
        istatus = nf90_put_att(incout, ivdim(1), 'units', 'km')
        call check_ok(istatus,'Error adding iy units')
        istatus = nf90_def_var(incout, 'jx', nf90_float, idims(1),      &
                            &  ivdim(2))
        call check_ok(istatus,'Error adding variable jx')
        istatus = nf90_put_att(incout, ivdim(2), 'standard_name',       &
                            &  'projection_x_coordinate')
        call check_ok(istatus,'Error adding jx standard_name')
        istatus = nf90_put_att(incout, ivdim(2), 'long_name',           &
                            &  'x-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding jx long_name')
        istatus = nf90_put_att(incout, ivdim(2), 'units', 'km')
        call check_ok(istatus,'Error adding jx units')

        ! XLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1))
#endif
        call check_ok(istatus,'Error adding variable xlat')
        istatus = nf90_put_att(incout, ivar(1), 'standard_name',        &
                            &  'latitude')
        call check_ok(istatus,'Error adding xlat standard_name')
        istatus = nf90_put_att(incout, ivar(1), 'long_name',            &
                            &  'Latitude at cross points')
        call check_ok(istatus,'Error adding xlat long_name')
        istatus = nf90_put_att(incout, ivar(1), 'units',                &
                            &  'degrees_north')
        call check_ok(istatus,'Error adding xlat units')
        ! XLAT

        ! XLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2))
#endif
        call check_ok(istatus,'Error adding variable xlon')
        istatus = nf90_put_att(incout, ivar(2), 'standard_name',        &
                            &  'longitude')
        call check_ok(istatus,'Error adding xlon standard_name')
        istatus = nf90_put_att(incout, ivar(2), 'long_name',            &
                            &  'Longitude at cross points')
        call check_ok(istatus,'Error adding xlon long_name')
        istatus = nf90_put_att(incout, ivar(2), 'units',                &
                            &  'degrees_east')
        call check_ok(istatus,'Error adding xlon units')
        ! XLON

        ! DLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3))
#endif
        call check_ok(istatus,'Error adding variable dlat')
        istatus = nf90_put_att(incout, ivar(3), 'standard_name',        &
                            &  'latitude')
        call check_ok(istatus,'Error adding dlat standard_name')
        istatus = nf90_put_att(incout, ivar(3), 'long_name',            &
                            &  'Latitude at dot points')
        call check_ok(istatus,'Error adding dlat long_name')
        istatus = nf90_put_att(incout, ivar(3), 'units',                &
                            &  'degrees_north')
        call check_ok(istatus,'Error adding dlat units')
        ! DLAT

        ! DLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4))
#endif
        call check_ok(istatus,'Error adding variable dlon')
        istatus = nf90_put_att(incout, ivar(4), 'standard_name',        &
                            &  'longitude')
        call check_ok(istatus,'Error adding dlon standard_name')
        istatus = nf90_put_att(incout, ivar(4), 'long_name',            &
                            &  'Longitude at dot points')
        call check_ok(istatus,'Error adding dlon long_name')
        istatus = nf90_put_att(incout, ivar(4), 'units',                &
                            &  'degrees_east')
        call check_ok(istatus,'Error adding dlon units')
        ! DLON

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5))
#endif
        call check_ok(istatus,'Error adding variable topo')
        istatus = nf90_put_att(incout, ivar(5), 'standard_name',        &
                            &  'surface_altitude')
        call check_ok(istatus,'Error adding topo standard_name')
        istatus = nf90_put_att(incout, ivar(5), 'long_name',            &
                            &  'Domain surface elevation')
        call check_ok(istatus,'Error adding topo long_name')
        istatus = nf90_put_att(incout, ivar(5), 'units', 'm')
        call check_ok(istatus,'Error adding topo units')
        istatus = nf90_put_att(incout, ivar(5), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding topo coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6))
#endif
        call check_ok(istatus,'Error adding variable htsd')
        istatus = nf90_put_att(incout, ivar(6), 'cell_method',          &
                            &  'area: standard_deviation')
        call check_ok(istatus,'Error adding htsd cell_method')
        istatus = nf90_put_att(incout, ivar(6), 'standard_name',        &
                            &  'surface_altitude')
        call check_ok(istatus,'Error adding htsd standard_name')
        istatus = nf90_put_att(incout, ivar(6), 'long_name',            &
                      &  'Domain surface elevation standard deviation')
        call check_ok(istatus,'Error adding htsd long_name')
        istatus = nf90_put_att(incout, ivar(6), 'units', 'm')
        call check_ok(istatus,'Error adding htsd units')
        istatus = nf90_put_att(incout, ivar(6), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding htsd coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7))
#endif
        call check_ok(istatus,'Error adding variable landuse')
        istatus = nf90_put_att(incout, ivar(7), 'legend',               &
                & '1  => Crop/mixed farming'//char(10)//                &
                & '2  => Short grass'//char(10)//                       &
                & '3  => Evergreen needleleaf tree'//char(10)//         &
                & '4  => Deciduous needleleaf tree'//char(10)//         &
                & '5  => Deciduous broadleaf tree'//char(10)//          &
                & '6  => Evergreen broadleaf tree'//char(10)//          &
                & '7  => Tall grass'//char(10)//                        &
                & '8  => Desert'//char(10)//                            &
                & '9  => Tundra'//char(10)//                            &
                & '10 => Irrigated Crop'//char(10)//                    &
                & '11 => Semi-desert'//char(10)//                       &
                & '12 => Ice cap/glacier'//char(10)//                   &
                & '13 => Bog or marsh'//char(10)//                      &
                & '14 => Inland water'//char(10)//                      &
                & '15 => Ocean'//char(10)//                             &
                & '16 => Evergreen shrub'//char(10)//                   &
                & '17 => Deciduous shrub'//char(10)//                   &
                & '18 => Mixed Woodland'//char(10)//                    &
                & '19 => Forest/Field mosaic'//char(10)//               &
                & '20 => Water and Land mixture')
        call check_ok(istatus,'Error adding landuse legend')
        istatus = nf90_put_att(incout, ivar(7), 'standard_name',        &
                            &  'land_type')
        call check_ok(istatus,'Error adding landuse standard_name')
        istatus = nf90_put_att(incout, ivar(7), 'long_name',            &
                      &  'Landuse category as defined in BATS1E')
        call check_ok(istatus,'Error adding landuse long_name')
        istatus = nf90_put_att(incout, ivar(7), 'units', '1')
        call check_ok(istatus,'Error adding landuse units')
        istatus = nf90_put_att(incout, ivar(7), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding landuse coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8))
#endif
        call check_ok(istatus,'Error adding variable xmap')
        istatus = nf90_put_att(incout, ivar(8), 'standard_name',        &
                            &  'map_factor')
        call check_ok(istatus,'Error adding xmap standard_name')
        istatus = nf90_put_att(incout, ivar(8), 'long_name',            &
                            &  'Map factor in domain cross points')
        call check_ok(istatus,'Error adding xmap long_name')
        istatus = nf90_put_att(incout, ivar(8), 'units', '1')
        call check_ok(istatus,'Error adding xmap units')
        istatus = nf90_put_att(incout, ivar(8), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding xmap coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9))
#endif
        call check_ok(istatus,'Error adding variable dmap')
        istatus = nf90_put_att(incout, ivar(9), 'standard_name',        &
                            &  'map_factor')
        call check_ok(istatus,'Error adding dmap standard_name')
        istatus = nf90_put_att(incout, ivar(9), 'long_name',            &
                            &  'Map factor in domain dot points')
        call check_ok(istatus,'Error adding dmap long_name')
        istatus = nf90_put_att(incout, ivar(9), 'units', '1')
        call check_ok(istatus,'Error adding dmap units')
        istatus = nf90_put_att(incout, ivar(9), 'coordinates',          &
                            &  'dlon dlat')
        call check_ok(istatus,'Error adding dmap coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10))
#endif
        call check_ok(istatus,'Error adding variable coriol')
        istatus = nf90_put_att(incout, ivar(10), 'standard_name',       &
                            &  'coriolis_parameter')
        call check_ok(istatus,'Error adding coriol standard_name')
        istatus = nf90_put_att(incout, ivar(10), 'long_name',           &
                            &  'Coriolis force parameter')
        call check_ok(istatus,'Error adding coriol long_name')
        istatus = nf90_put_att(incout, ivar(10), 'units', 's-1')
        call check_ok(istatus,'Error adding coriol units')
        istatus = nf90_put_att(incout, ivar(10), 'coordinates',         &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding coriol coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11))
#endif
        call check_ok(istatus,'Error adding variable snowam')
        istatus = nf90_put_att(incout, ivar(11), 'standard_name',       &
                            &  'snowfall_amount')
        call check_ok(istatus,'Error adding snowam standard_name')
        istatus = nf90_put_att(incout, ivar(11), 'long_name',           &
                            &  'Snow initial amount in mm')
        call check_ok(istatus,'Error adding snowam long_name')
        istatus = nf90_put_att(incout, ivar(11), 'units', 'kg m-2')
        call check_ok(istatus,'Error adding snowam units')
        istatus = nf90_put_att(incout, ivar(11), 'coordinates',         &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding snowam coordinates')

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12))
#endif
        call check_ok(istatus,'Error adding variable mask')
        istatus = nf90_put_att(incout, ivar(12), 'standard_name',       &
                            &  'land_binary_mask')
        call check_ok(istatus,'Error adding mask standard_name')
        istatus = nf90_put_att(incout, ivar(12), 'long_name',           &
                            &  'Land Sea mask')
        call check_ok(istatus,'Error adding mask long_name')
        istatus = nf90_put_att(incout, ivar(12), 'units', '1')
        call check_ok(istatus,'Error adding mask units')
        istatus = nf90_put_att(incout, ivar(12), 'coordinates',         &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding mask coordinates')

        if ( aertyp(7:7)=='1' ) then
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1))
#endif
          call check_ok(istatus,'Error adding variable texture')
          istatus = nf90_put_att(incout, itvar(1), 'legend',            &
                & '1  => Sand'//char(10)//                              &
                & '2  => Loamy Sand'//char(10)//                        &
                & '3  => Sandy Loam'//char(10)//                        &
                & '4  => Silt Loam'//char(10)//                         &
                & '5  => Silt'//char(10)//                              &
                & '6  => Loam'//char(10)//                              &
                & '7  => Sandy Clay Loam'//char(10)//                   &
                & '8  => Silty Clay Loam'//char(10)//                   &
                & '9  => Clay Loam'//char(10)//                         &
                & '10 => Sandy Clay'//char(10)//                        &
                & '11 => Silty Clay'//char(10)//                        &
                & '12 => Clay'//char(10)//                              &
                & '13 => OM'//char(10)//                                &
                & '14 => Water'//char(10)//                             &
                & '15 => Bedrock'//char(10)//                           &
                & '16 => Other'//char(10)//                             &
                & '17 => No data')
          call check_ok(istatus,'Error adding texture legend')
          istatus = nf90_put_att(incout, itvar(1), 'standard_name',     &
                              &  'soil_type')
          call check_ok(istatus,'Error adding texture standard_name')
          istatus = nf90_put_att(incout, itvar(1), 'long_name',         &
                              &  'Texture dominant category')
          call check_ok(istatus,'Error adding texture long_name')
          istatus = nf90_put_att(incout, itvar(1), 'units', '1')
          call check_ok(istatus,'Error adding texture units')
          istatus = nf90_put_att(incout, itvar(1), 'coordinates',       &
                              &  'xlon xlat')
          call check_ok(istatus,'Error adding texture coordinates')

#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims(1:3), itvar(2), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims(1:3), itvar(2))
#endif
          call check_ok(istatus,'Error adding variable texture_frac')
          istatus = nf90_put_att(incout, itvar(2), 'standard_name',     &
                              &  'soil_type_fraction')
          call check_ok(istatus,'Error adding text_frac standard_name')
          istatus = nf90_put_att(incout, itvar(2), 'long_name',         &
                              &  'Texture category fraction')
          call check_ok(istatus,'Error adding text_frac long_name')
          istatus = nf90_put_att(incout, itvar(2), 'units', '1')
          call check_ok(istatus,'Error adding text_frac units')
          istatus = nf90_put_att(incout, itvar(2), 'coordinates',       &
                              &  'xlon xlat')
          call check_ok(istatus,'Error adding text_frac coordinates')
        end if
!
!-----------------------------------------------------------------------
!
        istatus = nf90_enddef(incout)
        call check_ok(istatus,'Error End Definitions NetCDF output')
!
!-----------------------------------------------------------------------
!
        istatus = nf90_put_var(incout, izdim(1), sigma)
        call check_ok(istatus,'Error variable sigma write')
        hptop = ptop * 10.0
        istatus = nf90_put_var(incout, izdim(2), hptop)
        call check_ok(istatus,'Error variable ptop write')
        istatus = nf90_put_var(incout, ivdim(1), yiy)
        call check_ok(istatus,'Error variable iy write')
        istatus = nf90_put_var(incout, ivdim(2), xjx)
        call check_ok(istatus,'Error variable jx write')

        if (lsub) then
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat_s))
        else
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat))
        end if
        call check_ok(istatus,'Error variable xlat write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon_s))
        else
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon))
        end if
        call check_ok(istatus,'Error variable xlon write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat_s))
        else
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat))
        end if
        call check_ok(istatus,'Error variable dlat write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon_s))
        else
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon))
        end if
        call check_ok(istatus,'Error variable dlon write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid))
        end if
        call check_ok(istatus,'Error variable topo write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid))
        end if
        call check_ok(istatus,'Error variable htsd write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout_s))
        else
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout))
        end if
        call check_ok(istatus,'Error variable landuse write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap_s))
        else
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap))
        end if
        call check_ok(istatus,'Error variable xmap write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap_s))
        else
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap))
        end if
        call check_ok(istatus,'Error variable dmap write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol_s))
        else
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol))
        end if
        call check_ok(istatus,'Error variable coriol write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam_s))
        else
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam))
        end if
        call check_ok(istatus,'Error variable snowam write')
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(12), transpose(mask_s))
        else
          istatus = nf90_put_var(incout, ivar(12), transpose(mask))
        end if
        call check_ok(istatus,'Error variable mask write')

        if ( aertyp(7:7)=='1' ) then
          if (lsub) then
            istatus = nf90_put_var(incout, itvar(1), transpose(texout_s))
          else
            istatus = nf90_put_var(incout, itvar(1), transpose(texout))
          end if
          call check_ok(istatus,'Error variable texture write')
          istart(1) = 1
          istart(2) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          do i = 1 , ntex
            istart(3) = i
            if (lsub) then
              istatus = nf90_put_var(incout, itvar(2),                  &
                         & transpose(frac_tex_s(:,:,i)),istart,icount)
            else
              istatus = nf90_put_var(incout, itvar(2),                  &
                         & transpose(frac_tex(:,:,i)),istart,icount)
            end if
            call check_ok(istatus,'Error variable texture_frac write')
          end do
        end if

        istatus = nf90_close(incout)
        call check_ok(istatus, &
                 & ('Error closing NetCDF output '//trim(fname)))

        deallocate(yiy)
        deallocate(xjx)

      end subroutine write_domain

      subroutine check_ok(ierr,message)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        character(*) :: message
        if (ierr /= nf90_noerr) then 
          write (6,*) message
          write (6,*) nf90_strerror(ierr)
          stop
        end if
      end subroutine check_ok

      end module mod_write
