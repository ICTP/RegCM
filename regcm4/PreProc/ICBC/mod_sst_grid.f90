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
      use netcdf
      use mod_dynparam
      implicit none

      integer , private :: ncid
      integer , dimension(4) , private :: idims
      integer , dimension(4) , private :: ivar
      integer , private :: irefdate
      integer , private :: itime

      real(4) , allocatable , dimension(:,:) :: lu , sstmm , icemm ,    &
                                  &             xlat , xlon , finmat
      real(4) , allocatable , dimension(:) :: sigma

      contains

      subroutine init_grid
        implicit none
        allocate(lu(iy,jx))
        allocate(sstmm(iy,jx))
        allocate(icemm(iy,jx))
        allocate(xlat(iy,jx))
        allocate(xlon(iy,jx))
        allocate(finmat(jx,iy))
        allocate(sigma(kzp1))
      end subroutine init_grid

      subroutine free_grid
        implicit none
        deallocate(lu)
        deallocate(sstmm)
        deallocate(icemm)
        deallocate(xlat)
        deallocate(xlon)
        deallocate(finmat)
        deallocate(sigma)
      end subroutine free_grid

      subroutine read_domain(terfile)
        implicit none
        character(256) :: terfile
        intent(in) :: terfile

        integer :: iyy , jxx
        integer :: istatus , incin , idimid , ivarid

        istatus = nf90_open(terfile, nf90_nowrite, incin)
        if ( istatus /= nf90_noerr) then
          write (6,*) 'Error Opening Domain file ', trim(terfile)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_inq_dimid(incin, "iy", idimid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Dimension iy missing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dimension iy'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_dimid(incin, "jx", idimid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Dimension jx missing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inquire_dimension(incin, idimid, len=jxx)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error dimension jx'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if ( iyy/=iy .or. jxx/=jx ) then
          print * , 'IMPROPER DIMENSION SPECIFICATION'
          print * , '  namelist   : ' , iy , jx
          print * , '  DOMAIN     : ' , iyy , jxx
          stop 'Dimensions mismatch'
        end if

        istatus = nf90_inq_varid(incin, "sigma", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error sigma variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
         istatus = nf90_get_var(incin, ivarid, sigma)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading sigma variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_inq_varid(incin, "landuse", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error landuse variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
         istatus = nf90_get_var(incin, ivarid, finmat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading landuse variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        lu = transpose(finmat)
        istatus = nf90_inq_varid(incin, "xlat", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error xlat variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
         end if
        istatus = nf90_get_var(incin, ivarid, finmat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading xlat variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        xlat = transpose(finmat)
        istatus = nf90_inq_varid(incin, "xlon", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error xlon variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_get_var(incin, ivarid, finmat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading xlon variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        xlon = transpose(finmat)
        istatus = nf90_close(incin)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error closing Domain file ', trim(terfile)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

      end subroutine read_domain

      subroutine open_sstfile(idate1)
        use mod_date , only : split_idate
        implicit none
        integer , intent(in) :: idate1
        integer :: istatus
        character(256) :: sstname , history
        character(64) :: csdate
        real(4) , dimension(2) :: trlat
        real(4) , allocatable , dimension(:) :: yiy
        real(4) , allocatable , dimension(:) :: xjx
        integer , dimension(2) :: ivvar
        integer , dimension(2) :: illvar
        integer , dimension(2) :: izvar
        integer , dimension(8) :: tvals
        integer :: iyy , im , id , ih , i , j

        irefdate = idate1
        itime = 1

        call split_idate(idate1,iyy,im,id,ih)
        write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
                & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

        sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
        istatus = nf90_create(sstname, nf90_clobber, ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating NetCDF output ', trim(sstname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, nf90_global, 'title',  &
                 & 'ICTP Regional Climatic model V4 SST program output')
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
                           &   'sst_source', ssttyp)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global data_interpolation'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!
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
        istatus = nf90_def_dim(ncid, 'kz', kz+1, idims(4))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating dimension kz'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!
!
        istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4),   &
                            &  izvar(1))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                            &  'Sigma at model layer midpoints')
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
        istatus = nf90_def_var(ncid, 'sst', nf90_float, idims(1:3),  &
                            &  ivar(2))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                            &  'sea_surface_temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                            &  'Sea Surface Temperature')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'units', 'K')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), '_FillValue', -9999.0)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst _FillValue attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(2), 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          istatus = nf90_def_var(ncid, 'ice', nf90_float, idims(1:3),  &
                              &  ivar(3))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice definition in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
                              &  'sea_ice_thickness')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice standard_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(3), 'long_name',     &
                              &  'Sea ice depth')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice long_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(3), 'units', 'mm')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice units attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(3), '_FillValue', -9999.0)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice _FillValue attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(ncid, ivar(3), 'coordinates', &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        istatus = nf90_def_var(ncid, 'landuse', nf90_float,idims(1:3),&
                            &  ivar(4))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'legend',               &
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
        istatus = nf90_put_att(ncid, ivar(4), 'standard_name',        &
                            &  'land_type')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'long_name',            &
                      &  'Landuse category as defined in BATS1E')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar(4), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!
!
        istatus = nf90_enddef(ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error End Definitions NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!
!
        istatus = nf90_put_var(ncid, izvar(1), sigma)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sigma write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, izvar(2), ptop)
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
        istatus = nf90_put_var(ncid, illvar(1), transpose(xlat))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, illvar(2), transpose(xlon))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

      end subroutine open_sstfile

      subroutine close_sstfile
        implicit none
        integer :: istatus
        istatus = nf90_close(ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Closing output sst file'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
      end subroutine close_sstfile

      subroutine writerec(idate,lice)
        use mod_date , only : idatediff
        implicit none
        integer , intent(in) :: idate
        logical , intent(in) , optional :: lice
        integer :: istatus
        logical :: ladd
        integer , dimension(1) :: istart1 , icount1
        integer , dimension(3) :: istart , icount
        real(8) , dimension(1) :: xdate
        ladd = .false.
        if (present(lice)) then
          ladd = lice
        end if
        istart(3) = itime
        istart(2) = 1
        istart(1) = 1
        istart1(1) = itime
        icount(3) = 1
        icount(2) = iy
        icount(1) = jx
        icount1(1) = 1
        xdate(1) = dble(idatediff(idate,irefdate))
        istatus = nf90_put_var(ncid, ivar(1), xdate, istart1, icount1)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable time write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(ncid, ivar(2), transpose(sstmm), &
                               istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable sst write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (ladd) then
          istatus = nf90_put_var(ncid, ivar(3), transpose(icemm), &
                                 istart, icount)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable ice write in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        istatus = nf90_put_var(ncid, ivar(4), transpose(lu), &
                               istart, icount)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        itime = itime + 1
      end subroutine writerec

      end module mod_sst_grid
