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

      implicit none

      integer , private :: ncid
      character(256) , private :: ofname
      integer , private :: irefdate
      integer , private :: itime
      integer , dimension(5) , private :: idims
      integer , dimension(8) , private :: ivar

      real(4) , allocatable , dimension(:,:) :: ps4 , ts4
      real(4) , allocatable , dimension(:,:,:) :: c4 , h4 , q4
      real(4) , allocatable , dimension(:,:,:) :: t4 , u4 , v4
      real(4) , allocatable , dimension(:,:,:) :: sulfate4

      data ncid /-1/

      contains

      subroutine init_output
      implicit none
        allocate(ps4(jx,iy))
        allocate(ts4(jx,iy))
        allocate(c4(jx,iy,kz))
        allocate(h4(jx,iy,kz))
        allocate(q4(jx,iy,kz))
        allocate(t4(jx,iy,kz))
        allocate(u4(jx,iy,kz))
        allocate(v4(jx,iy,kz))
        if ( dattyp=='EH5OM' .and. ehso4) then
          allocate(sulfate4(jx,iy,kz))
        end if
      end subroutine init_output

      subroutine free_output
        use netcdf
        implicit none
        integer :: istatus
        deallocate(ps4)
        deallocate(ts4)
        deallocate(c4)
        deallocate(h4)
        deallocate(q4)
        deallocate(t4)
        deallocate(u4)
        deallocate(v4)
        if ( dattyp=='EH5OM' .and. ehso4) then
          deallocate(sulfate4)
        end if
        if (ncid > 0) then
          istatus = nf90_close(ncid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error closing file ', trim(ofname)
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
      end subroutine free_output

      subroutine newfile(idate1)
        use mod_date
        use mod_grid , only : xlat , xlon , sigma2
        use netcdf
        implicit none
        integer , intent(in) :: idate1
        integer :: istatus
        integer :: iyy , im , id , ih , i , j
        integer , dimension(8) :: tvals
        integer , dimension(2) :: izvar
        integer , dimension(2) :: ivvar
        integer , dimension(2) :: illvar
        integer , dimension(4) :: x3ddim
        real(4) , allocatable , dimension(:) :: yiy
        real(4) , allocatable , dimension(:) :: xjx
        character(64) :: csdate , cdum
        character(256) :: history
        real(4) , dimension(2) :: trlat
        real(4) :: hptop

        if (ncid > 0) then
          istatus = nf90_close(ncid)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error closing file ', trim(ofname)
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if

        write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
              &     '_ICBC.', idate1, '.nc'

        irefdate = idate1
        itime = 1

        call split_idate(idate1,iyy,im,id,ih)
        write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
                & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

#ifdef NETCDF4_HDF5
        istatus = nf90_create(ofname, ior(nf90_clobber,nf90_hdf5), &
                              ncid)
#else
        istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating NetCDF output ', trim(ofname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_put_att(ncid, nf90_global, 'title',  &
                & 'ICTP Regional Climatic model V4 ICBC program output')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global title'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'institution', &
                 & 'ICTP')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global institution'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'source', &
                 & 'RegCM Model simulation SST output')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global source'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
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
             ' : Created by RegCM sst program'
        istatus = nf90_put_att(ncid, nf90_global, 'history', history)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global history'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'references', &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global references'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'experiment', &
                 & domname)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global experiment'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global projection'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'grid_size_in_meters', ds*1000.0)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global gridsize'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'latitude_of_projection_origin', clat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clat'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'longitude_of_projection_origin', clon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clon'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'latitude_of_projection_pole', plat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'longitude_of_projection_pole', plon)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plon'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        else if (iproj == 'LAMCON') then
          trlat(1) = truelatl
          trlat(2) = truelath
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'standard_parallel', trlat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global truelat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        istatus = nf90_put_att(ncid, nf90_global,  &
                           &   'global_data_source', dattyp)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global_data_source'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (dattyp == 'EH5OM') then
          if (ehso4) then
            cdum = 'Yes'
          else
            cdum = 'No'
          end if
          istatus = nf90_put_att(ncid, nf90_global,  &
                             &   'sulfate_data_present', cdum)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding sulfate_data_present'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension iy'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension jx'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension time'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension kz'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        x3ddim(1) = idims(1)
        x3ddim(2) = idims(2)
        x3ddim(3) = idims(4)
        x3ddim(4) = idims(3)
!
        istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4),   &
                            &  izvar(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')      
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                            &  'Sigma at model layers')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma axis attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma positive attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                     &         'sigma: sigma ps: ps ptop: ptop')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma formula_terms attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'ptop', nf90_float,         &
                           &   varid=izvar(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(2), 'standard_name',  &
                            &  'air_pressure')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(2), 'long_name',      &
                            &  'Pressure at model top')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), &
                            &  ivvar(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                            &  'projection_y_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(1), 'long_name',      &
                            &  'y-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), &
                            &  ivvar(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                            &  'projection_x_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(2), 'long_name',    &
                            &  'x-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2),  &
                            &  illvar(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(1), 'standard_name', &
                            &  'latitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(1), 'long_name',     &
                            &  'Latitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(1), 'units',         &
                            &  'degrees_north')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2),  &
                            &  illvar(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(2), 'standard_name', &
                            &  'longitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(2), 'long_name',     &
                            &  'Longitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, illvar(2), 'units',         &
                            &  'degrees_east')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3),  &
                            &  ivar(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable time definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(1), 'units', &
                       &   'hours since '//csdate)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable time units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3),  &
                            &  ivar(2), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3),  &
                            &  ivar(2))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                            &  'surface_air_pressure')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                            &  'Surface pressure')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'units', 'hPa')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 'ts', nf90_float, idims(1:3),  &
                            &  ivar(3), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 'ts', nf90_float, idims(1:3),  &
                            &  ivar(3))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
                            &  'surface_temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(3), 'long_name',     &
                            &  'Surface Temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(3), 'units', 'K')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(3), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 'u', nf90_float, x3ddim,  &
                            &  ivar(4), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 'u', nf90_float, x3ddim,  &
                            &  ivar(4))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'standard_name', &
                            &  'eastward_wind')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'long_name',     &
                            &  'U component (westerly) of wind')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'units', 'm s-1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 'v', nf90_float, x3ddim,  &
                            &  ivar(5), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 'v', nf90_float, x3ddim,  &
                            &  ivar(5))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(5), 'standard_name', &
                            &  'northward_wind')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(5), 'long_name',     &
                            &  'V component (southerly) of wind')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(5), 'units', 'm s-1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(5), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 't', nf90_float, x3ddim,  &
                            &  ivar(6), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 't', nf90_float, x3ddim,  &
                            &  ivar(6))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(6), 'standard_name', &
                            &  'air_temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(6), 'long_name',     &
                            &  'Temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(6), 'units', 'K')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(6), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(ncid, 'qv', nf90_float, x3ddim,  &
                            &  ivar(7), deflate_level=9)
#else
        istatus = nf90_def_var(ncid, 'qv', nf90_float, x3ddim,  &
                            &  ivar(7))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(7), 'standard_name', &
                            &  'humidity_mixing_ratio')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(7), 'long_name',     &
                            &  'Water vapor mixing ratio')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(7), 'units', 'kg kg-1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(7), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if ( dattyp=='EH5OM' .and. ehso4) then
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(ncid, 'so4', nf90_float, x3ddim, &
                              &  ivar(8), deflate_level=9)
#else
          istatus = nf90_def_var(ncid, 'so4', nf90_float, x3ddim, &
                              &  ivar(8))
#endif
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 definition ', &
                      & 'in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(8), 'standard_name', &
                              &  'atmosphere_sulfate_content')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 standard_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(8), 'long_name',     &
                              &  'Sulfate')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 long_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(8), 'units', 'kg m-2')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 units attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(8), 'coordinates', &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
!
        istatus = nf90_enddef(ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error End Definitions NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
        istatus = nf90_put_var(ncid, izvar(1), sigma2)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        hptop = ptop*10.0
        istatus = nf90_put_var(ncid, izvar(2), hptop)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ptop write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
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
        istatus = nf90_put_var(ncid, ivvar(1), yiy)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivvar(2), xjx)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        deallocate(yiy)
        deallocate(xjx)
        istatus = nf90_put_var(ncid, illvar(1), xlat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, illvar(2), xlon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

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
        real(8) , dimension(1) :: xdate
!
        istart1(1) = itime
        icount1(1) = 1
        xdate(1) = dble(idatediff(idate,irefdate))
        istatus = nf90_put_var(ncid, ivar(1), xdate, istart1, icount1)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable time write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istart(3) = itime
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = iy
        icount(1) = jx
        ps4 = (ps4+ptop)*10.0
        istatus = nf90_put_var(ncid, ivar(2), ps4, istart(1:3), &
                               icount(1:3))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ps write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivar(3), ts4, istart(1:3), &
                               icount(1:3))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ts write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istart(4) = itime
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = kz
        icount(2) = iy
        icount(1) = jx
        istatus = nf90_put_var(ncid, ivar(4), u4, istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable u write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivar(5), v4, istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable v write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivar(6), t4, istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable t write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivar(7), q4, istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable qv write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if ( dattyp=='EH5OM' .and. ehso4) then
          istatus = nf90_put_var(ncid, ivar(8), sulfate4,  &
                            &    istart, icount)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable so4 write in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        itime = itime + 1
!
      end subroutine writef
!
      end module mod_write
