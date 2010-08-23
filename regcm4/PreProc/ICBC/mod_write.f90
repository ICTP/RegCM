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

      integer :: noutrec

      real(4) , allocatable , dimension(:,:) :: ps4 , ts4
      real(4) , allocatable , dimension(:,:,:) :: c4 , h4 , q4
      real(4) , allocatable , dimension(:,:,:) :: t4 , u4 , v4
      real(4) , allocatable , dimension(:,:,:) :: sulfate4

      data noutrec /0/
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

        istatus = nf90_create(ofname, nf90_clobber, ncid)
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
        istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3),  &
                            &  ivar(2))
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
        istatus = nf90_def_var(ncid, 'ts', nf90_float, idims(1:3),  &
                            &  ivar(3))
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
        istatus = nf90_def_var(ncid, 'u', nf90_float, x3ddim,  &
                            &  ivar(4))
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
        istatus = nf90_def_var(ncid, 'v', nf90_float, x3ddim,  &
                            &  ivar(5))
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
        istatus = nf90_def_var(ncid, 't', nf90_float, x3ddim,  &
                            &  ivar(6))
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
        istatus = nf90_def_var(ncid, 'qv', nf90_float, x3ddim,  &
                            &  ivar(7))
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
          istatus = nf90_def_var(ncid, 'so4', nf90_float, x3ddim, &
                              &  ivar(8))
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
        integer :: i , j , k
        integer :: istatus
        integer , dimension(1) :: istart1 , icount1
        integer , dimension(4) :: istart , icount
        real(8) , dimension(1) :: xdate
!
        noutrec = noutrec + 1
        write (64,rec=noutrec) idate , jx , iy , kz
        do k = kz , 1 , -1
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((u4(i,j,k),i=1,jx),j=1,iy)
        end do
        do k = kz , 1 , -1
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((v4(i,j,k),i=1,jx),j=1,iy)
        end do
        do k = kz , 1 , -1
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((t4(i,j,k),i=1,jx),j=1,iy)
        end do
        do k = kz , 1 , -1
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((q4(i,j,k),i=1,jx),j=1,iy)
        end do
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jx),j=1,iy)
        noutrec = noutrec + 1
        write (64,rec=noutrec) ts4
        if ( dattyp=='EH5OM' .and. ehso4) then
          do k = kz , 1 , -1
            noutrec = noutrec + 1
            write (64,rec=noutrec) ((sulfate4(j,i,k),j=1,jx),i=1,iy)
          end do
        end if
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

      subroutine gradsctl(finame,idate,inumber)
      use mod_grid
      implicit none
!
! Dummy arguments
!
      character(*) :: finame
      integer :: idate , inumber
      intent (in) finame , idate , inumber
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , j , k , month , nday , nhour , nx , ny , nyear
!
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      open (71,file=trim(finame)//'.CTL',status='replace')
      write (71,'(a,a,a,i10)') 'dset ^',trim(domname),'_ICBC.',idate
      write (71,'(a)') 'title ICBC fields for RegCM domain'
      if ( ibigend==1 ) then
        write (71,'(a)') 'options big_endian'
      else
        write (71,'(a)') 'options little_endian'
      end if
      write (71,'(a)') 'undef -9999.'
      if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
        do j = 1 , jx
          if ( xlat(j,1)<alatmin ) alatmin = xlat(j,1)
          if ( xlat(j,iy)>alatmax ) alatmax = xlat(j,iy)
        end do
        do i = 1 , iy
          do j = 1 , jx
            if ( clon>=0.0 ) then
              if ( xlon(j,i)>=0.0 ) then
                alonmin = amin1(alonmin,xlon(j,i))
                alonmax = amax1(alonmax,xlon(j,i))
              else if ( abs(clon-xlon(j,i))<abs(clon-(xlon(j,i)+360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(j,i))
                alonmax = amax1(alonmax,xlon(j,i))
              else
                alonmin = amin1(alonmin,xlon(j,i)+360.)
                alonmax = amax1(alonmax,xlon(j,i)+360.)
              end if
            else if ( xlon(j,i)<0.0 ) then
              alonmin = amin1(alonmin,xlon(j,i))
              alonmax = amax1(alonmax,xlon(j,i))
            else if ( abs(clon-xlon(j,i))<abs(clon-(xlon(j,i)-360.)) )  &
                    & then
              alonmin = amin1(alonmin,xlon(j,i))
              alonmax = amax1(alonmax,xlon(j,i))
            else
              alonmin = amin1(alonmin,xlon(j,i)-360.)
              alonmax = amax1(alonmax,xlon(j,i)-360.)
            end if
          end do
        end do
        delx = ds*1000.0
        rlatinc = delx*0.001/111./2.
        rloninc = delx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
        centerj = jx/2.
        centeri = iy/2.
      end if
      if ( iproj=='LAMCON' ) then       ! Lambert projection
        write (71,99001) jx , iy , clat , clon , centerj , centeri ,    &
                       & truelatl , truelath , clon , delx , delx
        write (71,99002) nx + 2 , alonmin - rloninc , rloninc
        write (71,99003) ny + 2 , alatmin - rlatinc , rlatinc
      else if ( iproj=='POLSTR' ) then  !
      else if ( iproj=='NORMER' ) then
        write (71,99004) jx , xlon(1,1) , xlon(2,1) - xlon(1,1)
        write (71,99005) iy
        write (71,99006) (xlat(1,i),i=1,iy)
      else if ( iproj=='ROTMER' ) then
        write (*,*) 'Note that rotated Mercartor (ROTMER)' ,            &
                   &' projections are not supported by GrADS.'
        write (*,*) '  Although not exact, the eta.u projection' ,      &
                   &' in GrADS is somewhat similar.'
        write (*,*) ' FERRET, however, does support this projection.'
        write (71,99007) jx , iy , plon , plat , delx/111000. ,         &
                       & delx/111000.*.95238
        write (71,99002) nx + 2 , alonmin - rloninc , rloninc
        write (71,99003) ny + 2 , alatmin - rlatinc , rlatinc
      else
        write (*,*) 'Are you sure your map projection is correct ?'
        stop
      end if
      write (71,99008) kz , ((1013.25-ptop*10.)*sigma2(k)+ptop*10.,k=kz,&
                     & 1,-1)
      nyear = idate/1000000
      month = (idate-nyear*1000000)/10000
      nday = (idate-nyear*1000000-month*10000)/100
      nhour = mod(idate,100)
      write (71,99009) inumber , nhour , cday(nday) , cmonth(month) ,    &
                     & nyear
      if ( dattyp=='EH5OM' ) then
        if ( ehso4 ) then
          write (71,99010) 8
        else
          write (71,99010) 7
        end if
      else
        write (71,99010) 7
      end if
      write (71,'(a)') 'date 0 99 header information'
      if ( iproj=='LAMCON' ) then       ! Lambert projection
        write (71,99013) 'u   ' , kz , 'westerly wind    '
        write (71,99014) 'v   ' , kz , 'southerly wind   '
      else
        write (71,99012) 'u   ' , kz , 'westerly wind    '
        write (71,99012) 'v   ' , kz , 'southerly wind   '
      end if
      write (71,99012) 't   ' , kz , 'air temperature  '
      write (71,99012) 'q   ' , kz , 'specific moisture'
      write (71,99015) 'px  ' , 'surface pressure           '
      write (71,99015) 'ts  ' , 'surface air temperature    '
      if ( dattyp=='EH5OM' .and. ehso4 ) write (71,99012) 'so4 ' , kz , &
          &'sulfate amount   '
      write (71,'(a)') 'endvars'
      close (71)
99001 format ('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i2,' levels ',30F7.2)
99009 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' 6hr')
99010 format ('vars ',i1)
99012 format (a4,i2,' 0 ',a17)
99013 format (a4,i2,' 33,100 ',a17)
99014 format (a4,i2,' 34,100 ',a17)
99015 format (a4,'0 99 ',a26)
!
      end subroutine gradsctl
!
!-----------------------------------------------------------------------
!
      subroutine fexist(filnam)
      implicit none
!
! Dummy arguments
!
      character(*) :: filnam
      intent (inout) filnam
!
! Local variables
!
      logical :: there
      character(1) :: yesno
 
 100  continue
      inquire (file=filnam,exist=there)
      if ( there ) then
 150    continue
        print * , ' '
        print * , ' '
        print * , '**************************************************'
        print * , 'FILE ALREADY EXISTS:  ' , filnam
        print * , 'Do you want to overwrite the existing file? [y/n/q]'
        read (*,*) yesno
        if ( yesno=='y' ) then
          return
        else if ( yesno=='n' ) then
          print * , 'ENTER NEW FILE NAME'
          read (*,'(a)') filnam
          go to 100
        else if ( yesno=='q' ) then
          stop 999
        else
          go to 150
        end if
      end if
 
      end subroutine fexist

      end module mod_write
