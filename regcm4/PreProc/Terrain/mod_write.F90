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
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating NetCDF output ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_put_att(incout, nf90_global, 'title',            &
                 & 'ICTP Regional Climatic model V4 domain')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global title'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'institution',      &
                 & 'ICTP')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global institution'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'source',           &
                 & 'RegCM Model simulation Terrain output')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global source'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'Conventions',      &
                 & 'CF-1.4')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global Conventions'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        call date_and_time(values=tvals)
        write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
             ' : Created by RegCM terrain'
        istatus = nf90_put_att(incout, nf90_global, 'history', history)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global history'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'references',       &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global references'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'experiment',       &
                 & domname)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global experiment'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'projection', iproj)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global projection'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'grid_size_in_meters', ds*1000.0)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global gridsize'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'latitude_of_projection_origin', clat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clat'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'longitude_of_projection_origin', clon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clon'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'latitude_of_projection_pole', plat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'longitude_of_projection_pole', plon)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plon'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        else if (iproj == 'LAMCON') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'standard_parallel', trlat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global truelat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        if (ifanal) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Cressman type objective analysis')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Overlapping parabolic 16 points')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global data_interpolation'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (smthbdy) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'No')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global boundary_smoothing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lakadj) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'No')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global great_lakes_adjustment'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'minimum_h2o_pct_for_water', h2opct)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global minimum_h2o_pct_for_water'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (lsub) then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec_s)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global ntypec_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if (fudge_lnd_s) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global landuse_fudging_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex_s) then
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'No')
            end if
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error adding global texture_fudging_s'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global ntypec_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if (fudge_lnd) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global landuse_fudging'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex) then
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'No')
            end if
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error adding global texture_fudging'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
        end if

        istatus = nf90_put_att(incout, nf90_global, 'grid_factor', xn)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global grid_factor attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (lsub) then
          istatus = nf90_def_dim(incout, 'iy', iysg, idims(2))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension iy'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_def_dim(incout, 'jx', jxsg, idims(1))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension jx'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        else
          istatus = nf90_def_dim(incout, 'iy', iy, idims(2))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension iy'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_def_dim(incout, 'jx', jx, idims(1))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension jx'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        if ( aertyp(7:7)=='1' ) then
          istatus = nf90_def_dim(incout, 'ntex', ntex, idims(3))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension NVEG'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        istatus = nf90_def_dim(incout, 'kz', kz+1, idims(4))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension kz'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_def_var(incout, 'sigma', nf90_float, idims(4),   &
                            &  izdim(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'long_name',           &
                            &  'Sigma at model layer midpoints')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'axis', 'Z')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma axis attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'positive', 'down')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma positive attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(1), 'formula_terms',       &
                     &         'sigma: sigma ps: ps ptop: ptop')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma formula_terms attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(incout, 'ptop', nf90_float,              &
                           &   varid=izdim(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(2), 'standard_name',       &
                            &  'air_pressure')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(2), 'long_name',           &
                            &  'Pressure at model top')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, izdim(2), 'units', 'hPa')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(incout, 'iy', nf90_float, idims(2),      &
                            &  ivdim(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'standard_name',       &
                            &  'projection_y_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'long_name',           &
                            &  'y-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(incout, 'jx', nf90_float, idims(1),      &
                            &  ivdim(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'standard_name',       &
                            &  'projection_x_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'long_name',           &
                            &  'x-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        ! XLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'standard_name',        &
                            &  'latitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'long_name',            &
                            &  'Latitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'units',                &
                            &  'degrees_north')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! XLAT

        ! XLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'standard_name',        &
                            &  'longitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'long_name',            &
                            &  'Longitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'units',                &
                            &  'degrees_east')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! XLON

        ! DLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'standard_name',        &
                            &  'latitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'long_name',            &
                            &  'Latitude at dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'units',                &
                            &  'degrees_north')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! DLAT

        ! DLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'standard_name',        &
                            &  'longitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'long_name',            &
                            &  'Longitude at dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'units',                &
                            &  'degrees_east')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! DLON

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'standard_name',        &
                            &  'surface_altitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'long_name',            &
                            &  'Domain surface elevation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'units', 'm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'cell_method',          &
                            &  'area: standard_deviation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd cell_method attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'standard_name',        &
                            &  'surface_altitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'long_name',            &
                      &  'Domain surface elevation standard deviation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'units', 'm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
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
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse legend attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'standard_name',        &
                            &  'land_type')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'long_name',            &
                      &  'Landuse category as defined in BATS1E')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'standard_name',        &
                            &  'map_factor')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'long_name',            &
                            &  'Map factor in domain cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'standard_name',        &
                            &  'map_factor')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'long_name',            &
                            &  'Map factor in domain dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'coordinates',          &
                            &  'dlon dlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'standard_name',       &
                            &  'coriolis_parameter')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'long_name',           &
                            &  'Coriolis force parameter')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'units', 's-1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'standard_name',       &
                            &  'snowfall_amount')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'long_name',           &
                            &  'Snow initial amount in mm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'units', 'kg m-2')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'standard_name',       &
                            &  'land_binary_mask')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'long_name',           &
                            &  'Land Sea mask')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if ( aertyp(7:7)=='1' ) then
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1))
#endif
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture def in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
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
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture legend attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'standard_name',     &
                              &  'soil_type')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture standard_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'long_name',         &
                              &  'Texture dominant category')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture long_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'units', '1')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture units attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'coordinates',       &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if

#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims(1:3), itvar(2), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims(1:3), itvar(2))
#endif
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction def in NetCDF'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'standard_name',     &
                              &  'soil_type_fraction')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction standard_name'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'long_name',         &
                              &  'Texture category fraction')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction long_name'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'units', '1')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction units'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'coordinates',       &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
!
!-----------------------------------------------------------------------
!
        istatus = nf90_enddef(incout)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error End Definitions NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!-----------------------------------------------------------------------
!
        istatus = nf90_put_var(incout, izdim(1), sigma)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        hptop = ptop * 10.0
        istatus = nf90_put_var(incout, izdim(2), hptop)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(incout, ivdim(1), yiy)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(incout, ivdim(2), xjx)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (lsub) then
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat_s))
        else
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon_s))
        else
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat_s))
        else
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon_s))
        else
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout_s))
        else
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap_s))
        else
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap_s))
        else
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol_s))
        else
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam_s))
        else
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(12), transpose(mask_s))
        else
          istatus = nf90_put_var(incout, ivar(12), transpose(mask))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if ( aertyp(7:7)=='1' ) then
          if (lsub) then
            istatus = nf90_put_var(incout, itvar(1), transpose(texout_s))
          else
            istatus = nf90_put_var(incout, itvar(1), transpose(texout))
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture write in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
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
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error Variable texture_frac write in NetCDF'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end do
        end if

        istatus = nf90_close(incout)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error closing NetCDF output ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        deallocate(yiy)
        deallocate(xjx)

      end subroutine write_domain

      end module mod_write
