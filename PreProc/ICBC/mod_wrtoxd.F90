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

      private

      public :: oh4 , ho24 , o34 , no34 , h2o24
      public :: init_outoxd , free_outoxd , newfile , writeox

      integer :: ncid
      character(256) :: ofname
      integer :: irefdate
      integer :: itime
      integer , dimension(5) :: idims
      integer , dimension(6) :: ivar
      integer :: istatus

      real(sp) , allocatable , dimension(:,:,:) :: oh4 , ho24 , o34 , &
                                                   no34 , h2o24

      data ncid /-1/

      contains

      subroutine init_outoxd
        implicit none
        allocate(oh4(jx,iy,kz))
        allocate(ho24(jx,iy,kz))
        allocate(o34(jx,iy,kz))
        allocate(no34(jx,iy,kz))
        allocate(h2o24(jx,iy,kz))
      end subroutine init_outoxd

      subroutine free_outoxd
        use netcdf
        implicit none
        deallocate(oh4)
        deallocate(ho24)
        deallocate(o34)
        deallocate(no34)
        deallocate(h2o24)
        if (ncid > 0) then
          istatus = nf90_close(ncid)
          call check_ok(istatus,('Error closing file '//trim(ofname)))
        end if
      end subroutine free_outoxd

      subroutine newfile(idate1)
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
        real(sp) , allocatable , dimension(:) :: yiy
        real(sp) , allocatable , dimension(:) :: xjx
        character(64) :: csdate
        character(256) :: history
        real(sp) , dimension(2) :: trlat
        real(sp) :: hptop

        if (ncid > 0) then
          istatus = nf90_close(ncid)
          call check_ok(istatus,('Error closing file '//trim(ofname)))
        end if

        write (ofname,99001) trim(dirglob), pthsep, trim(domname),    &
              &     '_OXBC.', idate1, '.nc'

        irefdate = idate1
        itime = 1

        call split_idate(idate1,iyy,im,id,ih)
        write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
                & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

#ifdef NETCDF4_HDF5
        istatus = nf90_create(ofname, &
                  ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), &
                  ncid)
#else
        istatus = nf90_create(ofname, nf90_clobber, ncid)
#endif
        call check_ok(istatus,('Error creating file '//trim(ofname)))

        istatus = nf90_put_att(ncid, nf90_global, 'title',  &
            & 'ICTP Regional Climatic model V4 oxidant program output')
        call check_ok(istatus,'Error adding global title')
        istatus = nf90_put_att(ncid, nf90_global, 'institution', &
                 & 'ICTP')
        call check_ok(istatus,'Error adding global institution')
        istatus = nf90_put_att(ncid, nf90_global, 'source', &
                 & 'RegCM Model simulation SST output')
        call check_ok(istatus,'Error adding global source')
        istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
                 & 'CF-1.4')
        call check_ok(istatus,'Error adding global Conventions')
        call date_and_time(values=tvals)
        write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
             ' : Created by RegCM oxidant program'
        istatus = nf90_put_att(ncid, nf90_global, 'history', history)
        call check_ok(istatus,'Error adding global history')
        istatus = nf90_put_att(ncid, nf90_global, 'references', &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        call check_ok(istatus,'Error adding global references')
        istatus = nf90_put_att(ncid, nf90_global, 'experiment', &
                 & domname)
        call check_ok(istatus,'Error adding global experiment')
        istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
        call check_ok(istatus,'Error adding global projection')
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'grid_size_in_meters', ds*1000.0)
        call check_ok(istatus,'Error adding global gridsize')
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'latitude_of_projection_origin', clat)
        call check_ok(istatus,'Error adding global clat')
        istatus = nf90_put_att(ncid, nf90_global,   &
                 &   'longitude_of_projection_origin', clon)
        call check_ok(istatus,'Error adding global clon')
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'grid_north_pole_latitude', plat)
          call check_ok(istatus,'Error adding global plat')
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'grid_north_pole_longitude', plon)
          call check_ok(istatus,'Error adding global plon')
        else if (iproj == 'LAMCON') then
          trlat(1) = truelatl
          trlat(2) = truelath
          istatus = nf90_put_att(ncid, nf90_global, &
                   &   'standard_parallel', trlat)
          call check_ok(istatus,'Error adding global truelat')
        end if
        istatus = nf90_put_att(ncid, nf90_global,  &
                           &   'global_data_source', dattyp)
        call check_ok(istatus,'Error adding global data_source')
        istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
        call check_ok(istatus,'Error creating dimension iy')
        istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
        call check_ok(istatus,'Error creating dimension jx')
        istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
        call check_ok(istatus,'Error creating dimension time')
        istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
        call check_ok(istatus,'Error creating dimension kz')
        x3ddim(1) = idims(1)
        x3ddim(2) = idims(2)
        x3ddim(3) = idims(4)
        x3ddim(4) = idims(3)
!
        istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4),   &
                            &  izvar(1))
        call check_ok(istatus,'Error adding variable sigma')
        istatus = nf90_put_att(ncid, izvar(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')      
        call check_ok(istatus,'Error adding sigma standard_name')
        istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                            &  'Sigma at model layers')
        call check_ok(istatus,'Error adding sigma long_name')
        istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
        call check_ok(istatus,'Error adding sigma units')
        istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
        call check_ok(istatus,'Error adding sigma axis')
        istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
        call check_ok(istatus,'Error adding sigma positive')
        istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                     &         'sigma: sigma ps: ps ptop: ptop')
        call check_ok(istatus,'Error adding sigma formula_terms')
        istatus = nf90_def_var(ncid, 'ptop', nf90_float,         &
                           &   varid=izvar(2))
        call check_ok(istatus,'Error adding variable ptop')
        istatus = nf90_put_att(ncid, izvar(2), 'standard_name',  &
                            &  'air_pressure')
        call check_ok(istatus,'Error adding ptop standard_name')
        istatus = nf90_put_att(ncid, izvar(2), 'long_name',      &
                            &  'Pressure at model top')
        call check_ok(istatus,'Error adding ptop long_name')
        istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
        call check_ok(istatus,'Error adding ptop units')
        istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), &
                            &  ivvar(1))
        call check_ok(istatus,'Error adding variable iy')
        istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                            &  'projection_y_coordinate')
        call check_ok(istatus,'Error adding iy standard_name')
        istatus = nf90_put_att(ncid, ivvar(1), 'long_name',      &
                            &  'y-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding iy long_name')
        istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
        call check_ok(istatus,'Error adding iy units')
        istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), &
                            &  ivvar(2))
        call check_ok(istatus,'Error adding variable jx')
        istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                            &  'projection_x_coordinate')
        call check_ok(istatus,'Error adding jx standard_name')
        istatus = nf90_put_att(ncid, ivvar(2), 'long_name',    &
                            &  'x-coordinate in Cartesian system')
        call check_ok(istatus,'Error adding jx long_name')
        istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
        call check_ok(istatus,'Error adding jx units')
        istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2),  &
                            &  illvar(1))
        call check_ok(istatus,'Error adding variable xlat')
        istatus = nf90_put_att(ncid, illvar(1), 'standard_name', &
                            &  'latitude')
        call check_ok(istatus,'Error adding xlat standard_name')
        istatus = nf90_put_att(ncid, illvar(1), 'long_name',     &
                            &  'Latitude at cross points')
        call check_ok(istatus,'Error adding xlat long_name')
        istatus = nf90_put_att(ncid, illvar(1), 'units',         &
                            &  'degrees_north')
        call check_ok(istatus,'Error adding xlat units')
        istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2),  &
                            &  illvar(2))
        call check_ok(istatus,'Error adding variable xlon')
        istatus = nf90_put_att(ncid, illvar(2), 'standard_name', &
                            &  'longitude')
        call check_ok(istatus,'Error adding xlon standard_name')
        istatus = nf90_put_att(ncid, illvar(2), 'long_name',     &
                            &  'Longitude at cross points')
        call check_ok(istatus,'Error adding xlon long_name')
        istatus = nf90_put_att(ncid, illvar(2), 'units',         &
                            &  'degrees_east')
        call check_ok(istatus,'Error adding xlon units')
        istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3),  &
                            &  ivar(1))
        call check_ok(istatus,'Error adding variable time')
        istatus = nf90_put_att(ncid, ivar(1), 'units', &
                       &   'hours since '//csdate)
        call check_ok(istatus,'Error adding time units')
        istatus = nf90_def_var(ncid, 'oh', nf90_float, x3ddim,  &
                            &  ivar(2))
        call check_ok(istatus,'Error adding variable oh')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncid, ivar(2), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on oh')
#endif
        istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                            &  'hydroxide_molecular_density')
        call check_ok(istatus,'Error adding oh standard_name')
        istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                            & 'Hydroxide molecular density')
        call check_ok(istatus,'Error adding oh long_name')
        istatus = nf90_put_att(ncid, ivar(2), 'units', 'moleculs cm-3')
        call check_ok(istatus,'Error adding oh units')
        istatus = nf90_put_att(ncid, ivar(2), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding oh coordinates')
        istatus = nf90_def_var(ncid, 'ho2', nf90_float, x3ddim,  &
                            &  ivar(3))
        call check_ok(istatus,'Error adding variable ho2')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncid, ivar(3), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on ho2')
#endif
        istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
                          & 'hydrogen_dioxide_molecular_density')
        call check_ok(istatus,'Error adding ho2 standard_name')
        istatus = nf90_put_att(ncid, ivar(3), 'long_name',     &
                          & 'Hydrogen dioxide molecular density')
        call check_ok(istatus,'Error adding ho2 long_name')
        istatus = nf90_put_att(ncid, ivar(3), 'units', 'moleculs cm-3')
        call check_ok(istatus,'Error adding ho2 units')
        istatus = nf90_put_att(ncid, ivar(3), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding ho2 coordinates')
        istatus = nf90_def_var(ncid, 'o3', nf90_float, x3ddim,  &
                            &  ivar(4))
        call check_ok(istatus,'Error adding variable o3')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncid, ivar(4), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on o3')
#endif
        istatus = nf90_put_att(ncid, ivar(4), 'standard_name', &
                            &  'ozone_mixing_ratio')
        call check_ok(istatus,'Error adding o3 standard_name')
        istatus = nf90_put_att(ncid, ivar(4), 'long_name',     &
                            &  'Ozone mixing ratio')
        call check_ok(istatus,'Error adding o3 long_name')
        istatus = nf90_put_att(ncid, ivar(4), 'units', 'm^3/m^3')
        call check_ok(istatus,'Error adding o3 units')
        istatus = nf90_put_att(ncid, ivar(4), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding o3 coordinates')
        istatus = nf90_def_var(ncid, 'no3', nf90_float, x3ddim,  &
                            &  ivar(5))
        call check_ok(istatus,'Error adding variable no3')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncid, ivar(5), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on no3')
#endif
        istatus = nf90_put_att(ncid, ivar(5), 'standard_name', &
                            &  'nitrate_molecular_density')
        call check_ok(istatus,'Error adding no3 standard_name')
        istatus = nf90_put_att(ncid, ivar(5), 'long_name',     &
                            &  'Nitrate molecular density')
        call check_ok(istatus,'Error adding no3 long_name')
        istatus = nf90_put_att(ncid, ivar(5), 'units', 'moleculs cm^3')
        call check_ok(istatus,'Error adding no3 units')
        istatus = nf90_put_att(ncid, ivar(5), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding no3 coordinates')
        istatus = nf90_def_var(ncid, 'h2o2', nf90_float, x3ddim, &
                              &  ivar(6))
        call check_ok(istatus,'Error adding variable h2o2')
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var_deflate(ncid, ivar(6), 1, 1, 9)
        call check_ok(istatus,'Error setting compression on h2o2')
#endif
        istatus = nf90_put_att(ncid, ivar(6), 'standard_name', &
                           &  'hydrogen_peroxide_mixing_ratio')
        call check_ok(istatus,'Error adding h2o2 standard_name')
        istatus = nf90_put_att(ncid, ivar(6), 'long_name',     &
                           &  'Hydrogen peroxide mixing ratio')
        call check_ok(istatus,'Error adding h2o2 long_name')
        istatus = nf90_put_att(ncid, ivar(6), 'units', 'm^3/m^3')
        call check_ok(istatus,'Error adding h2o2 units')
        istatus = nf90_put_att(ncid, ivar(6), 'coordinates', &
                           &  'xlon xlat')
        call check_ok(istatus,'Error adding h2o2 coordinates')
!
        istatus = nf90_enddef(ncid)
        call check_ok(istatus,'Error End Definitions NetCDF output')
!
        istatus = nf90_put_var(ncid, izvar(1), sigma2)
        call check_ok(istatus,'Error variable sigma write')
        hptop = ptop*10.0
        istatus = nf90_put_var(ncid, izvar(2), hptop)
        call check_ok(istatus,'Error variable ptop write')
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
        call check_ok(istatus,'Error variable iy write')
        istatus = nf90_put_var(ncid, ivvar(2), xjx)
        call check_ok(istatus,'Error variable jx write')
        deallocate(yiy)
        deallocate(xjx)
        istatus = nf90_put_var(ncid, illvar(1), xlat)
        call check_ok(istatus,'Error variable xlat write')
        istatus = nf90_put_var(ncid, illvar(2), xlon)
        call check_ok(istatus,'Error variable xlon write')

99001 format (a,a,a,a,i10,a)

      end subroutine newfile

      subroutine writeox(idate)
        use netcdf
        implicit none
        integer , intent(in) :: idate
        integer :: istatus
        integer , dimension(1) :: istart1 , icount1
        integer , dimension(4) :: istart , icount
        real(dp) , dimension(1) :: xdate
!
        istart1(1) = itime
        icount1(1) = 1
        xdate(1) = dble(idatediff(idate,irefdate))
        istatus = nf90_put_var(ncid, ivar(1), xdate, istart1, icount1)
        call check_ok(istatus,'Error variable time write')
        istart(4) = itime
        istart(3) = 1
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = kz
        icount(2) = iy
        icount(1) = jx
        istatus = nf90_put_var(ncid, ivar(2), oh4, istart, icount)
        call check_ok(istatus,'Error variable oh write')
        istatus = nf90_put_var(ncid, ivar(3), ho24, istart, icount)
        call check_ok(istatus,'Error variable ho2 write')
        istatus = nf90_put_var(ncid, ivar(4), o34, istart, icount)
        call check_ok(istatus,'Error variable o3 write')
        istatus = nf90_put_var(ncid, ivar(5), no34, istart, icount)
        call check_ok(istatus,'Error variable no3 write')
        istatus = nf90_put_var(ncid, ivar(6), h2o24, istart, icount)
        call check_ok(istatus,'Error variable h2o2 write')
        itime = itime + 1
!
      end subroutine writeox
!
      subroutine check_ok(ierr,message)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        character(*) :: message
        if (ierr /= nf90_noerr) then
          call die('writeox',message,1,nf90_strerror(ierr),ierr)
        end if
      end subroutine check_ok
!
      end module mod_wrtoxd
