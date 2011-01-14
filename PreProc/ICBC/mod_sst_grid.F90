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
      use m_realkinds
      use m_die
      use m_zeit
      use m_mall
      use netcdf
      use mod_dynparam

      integer , private :: ncid
      integer , dimension(4) , private :: idims
      integer , dimension(4) , private :: ivar
      integer , private :: irefdate
      integer , private :: itime

      real(sp) , allocatable , dimension(:,:) :: lu , sstmm , icemm ,   &
                                  &             xlat , xlon , finmat
      real(sp) , allocatable , dimension(:) :: sigma
      real(sp) , allocatable , dimension(:) :: yiy
      real(sp) , allocatable , dimension(:) :: xjx

      contains

      subroutine init_grid
        implicit none
        integer :: ierr
        allocate(lu(iy,jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate lu',ierr)
        call mall_mci(lu,'mod_sst_grid')
        allocate(sstmm(iy,jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate sstmm',ierr)
        call mall_mci(sstmm,'mod_sst_grid')
        allocate(icemm(iy,jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate icemm',ierr)
        call mall_mci(icemm,'mod_sst_grid')
        allocate(xlat(iy,jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate xlat',ierr)
        call mall_mci(xlat,'mod_sst_grid')
        allocate(xlon(iy,jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate xlon',ierr)
        call mall_mci(xlon,'mod_sst_grid')
        allocate(finmat(jx,iy), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate finmat',ierr)
        call mall_mci(finmat,'mod_sst_grid')
        allocate(sigma(kzp1), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate sigma',ierr)
        call mall_mci(sigma,'mod_sst_grid')
        allocate(yiy(iy), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate yiy',ierr)
        call mall_mci(yiy,'mod_sst_grid')
        allocate(xjx(jx), stat=ierr)
        if (ierr /= 0) call die('init_grid','allocate xjx',ierr)
        call mall_mci(xjx,'mod_sst_grid')
      end subroutine init_grid

      subroutine free_grid
        implicit none
        call mall_mco(lu,'mod_sst_grid')
        deallocate(lu)
        call mall_mco(sstmm,'mod_sst_grid')
        deallocate(sstmm)
        call mall_mco(icemm,'mod_sst_grid')
        deallocate(icemm)
        call mall_mco(xlat,'mod_sst_grid')
        deallocate(xlat)
        call mall_mco(xlon,'mod_sst_grid')
        deallocate(xlon)
        call mall_mco(finmat,'mod_sst_grid')
        deallocate(finmat)
        call mall_mco(sigma,'mod_sst_grid')
        deallocate(sigma)
        call mall_mco(yiy,'mod_sst_grid')
        deallocate(yiy)
        call mall_mco(xjx,'mod_sst_grid')
        deallocate(xjx)
      end subroutine free_grid

      subroutine read_domain(terfile)
        implicit none
        character(256) :: terfile
        intent(in) :: terfile

        integer :: iyy , jxx
        integer :: istatus , incin , idimid , ivarid

        call zeit_ci('read_domain')
        istatus = nf90_open(terfile, nf90_nowrite, incin)
        call check_ok(istatus, &
             &        'Error Opening Domain file '//trim(terfile))

        istatus = nf90_inq_dimid(incin, "iy", idimid)
        call check_ok(istatus,'Dimension iy missing')
        istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
        call check_ok(istatus,'Dimension iy read error')
        istatus = nf90_inq_dimid(incin, "jx", idimid)
        call check_ok(istatus,'Dimension jx missing')
        istatus = nf90_inquire_dimension(incin, idimid, len=jxx)
        call check_ok(istatus,'Dimension jx read error')
        if ( iyy/=iy .or. jxx/=jx ) then
          write (stderr,*) 'IMPROPER DIMENSION SPECIFICATION'
          write (stderr,*) '  namelist   : ' , iy , jx
          write (stderr,*) '  DOMAIN     : ' , iyy , jxx
          call die('read_domain','Dimensions mismatch',1)
        end if

        istatus = nf90_inq_varid(incin, "sigma", ivarid)
        call check_ok(istatus,'Variable sigma missing')
        istatus = nf90_get_var(incin, ivarid, sigma)
        call check_ok(istatus,'Variable sigma read error')
        istatus = nf90_inq_varid(incin, "landuse", ivarid)
        call check_ok(istatus,'Variable landuse missing')
        istatus = nf90_get_var(incin, ivarid, finmat)
        call check_ok(istatus,'Variable landuse read error')
        lu = transpose(finmat)
        istatus = nf90_inq_varid(incin, "xlat", ivarid)
        call check_ok(istatus,'Variable xlat missing')
        istatus = nf90_get_var(incin, ivarid, finmat)
        call check_ok(istatus,'Variable xlat read error')
        xlat = transpose(finmat)
        istatus = nf90_inq_varid(incin, "xlon", ivarid)
        call check_ok(istatus,'Variable xlon missing')
        istatus = nf90_get_var(incin, ivarid, finmat)
        call check_ok(istatus,'Variable xlon read error')
        xlon = transpose(finmat)
        istatus = nf90_close(incin)
        call check_ok(istatus, &
             &        ('Error closing Domain file '//trim(terfile)))
        call zeit_co('read_domain')

      end subroutine read_domain

      subroutine open_sstfile(idate1)
        use mod_date , only : split_idate
        implicit none
        integer , intent(in) :: idate1
        integer :: istatus
        character(256) :: sstname , history
        character(64) :: csdate
        real(sp) , dimension(2) :: trlat
        real(sp) :: hptop
        integer , dimension(2) :: ivvar
        integer , dimension(2) :: illvar
        integer , dimension(2) :: izvar
        integer , dimension(8) :: tvals
        integer :: iyy , im , id , ih , i , j

        call zeit_ci('newfile')
        irefdate = idate1
        itime = 1

        call split_idate(idate1,iyy,im,id,ih)
        write (csdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
                & iyy,'-',im,'-',id,' ',ih,':00:00 UTC'

        sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
#ifdef NETCDF4_HDF5
        istatus = nf90_create(sstname, &
                  ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model), &
                  ncid)
#else
        istatus = nf90_create(sstname, nf90_clobber, ncid)
#endif
        call check_ok(istatus, &
                      ('Error creating NetCDF output '//trim(sstname)))
        istatus = nf90_put_att(ncid, nf90_global, 'title',  &
                 & 'ICTP Regional Climatic model V4 SST program output')
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
             ' : Created by RegCM sst program'
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
                           &   'sst_source', ssttyp)
        call check_ok(istatus,'Error adding global sst_source')
!
        istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
        call check_ok(istatus,'Error creating dimension iy')
        istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
        call check_ok(istatus,'Error creating dimension jx')
        istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
        call check_ok(istatus,'Error creating dimension time')
        istatus = nf90_def_dim(ncid, 'kz', kz+1, idims(4))
        call check_ok(istatus,'Error creating dimension kz')
!
        istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4),   &
                            &  izvar(1))
        call check_ok(istatus,'Error adding variable sigma')
        istatus = nf90_put_att(ncid, izvar(1), 'standard_name',       &
                            &  'atmosphere_sigma_coordinate')
        call check_ok(istatus,'Error adding sigma standard_name')
        istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                            &  'Sigma at model layer midpoints')
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
        call check_ok(istatus,'Error adding iy uits')
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
        istatus = nf90_def_var(ncid, 'sst', nf90_float, idims(1:3),  &
                            &  ivar(2))
        call check_ok(istatus,'Error adding variable sst')
        istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                            &  'sea_surface_temperature')
        call check_ok(istatus,'Error adding sst standard_name')
        istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                            &  'Sea Surface Temperature')
        call check_ok(istatus,'Error adding sst long_name')
        istatus = nf90_put_att(ncid, ivar(2), 'units', 'K')
        call check_ok(istatus,'Error adding sst units')
        istatus = nf90_put_att(ncid, ivar(2), '_FillValue', -9999.0)
        call check_ok(istatus,'Error adding sst _FillValue')
        istatus = nf90_put_att(ncid, ivar(2), 'coordinates', &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding sst coordinates')
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          istatus = nf90_def_var(ncid, 'ice', nf90_float, idims(1:3),  &
                              &  ivar(3))
          call check_ok(istatus,'Error adding variable ice')
          istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
                              &  'sea_ice_thickness')
          call check_ok(istatus,'Error adding ice standard_name')
          istatus = nf90_put_att(ncid, ivar(3), 'long_name',     &
                              &  'Sea ice depth')
          call check_ok(istatus,'Error adding ice long_name')
          istatus = nf90_put_att(ncid, ivar(3), 'units', 'mm')
          call check_ok(istatus,'Error adding ice units')
          istatus = nf90_put_att(ncid, ivar(3), '_FillValue', -9999.0)
          call check_ok(istatus,'Error adding ice _FillValue')
          istatus = nf90_put_att(ncid, ivar(3), 'coordinates', &
                              &  'xlon xlat')
          call check_ok(istatus,'Error adding ice coordinates')
        end if
        istatus = nf90_def_var(ncid, 'landuse', nf90_float,idims(1:3),&
                            &  ivar(4))
        call check_ok(istatus,'Error adding variable landuse')
        istatus = nf90_put_att(ncid, ivar(4), 'legend',               &
                & '1  => Crop/mixed farming'//char(10)//              &
                & '2  => Short grass'//char(10)//                     &
                & '3  => Evergreen needleleaf tree'//char(10)//       &
                & '4  => Deciduous needleleaf tree'//char(10)//       &
                & '5  => Deciduous broadleaf tree'//char(10)//        &
                & '6  => Evergreen broadleaf tree'//char(10)//        &
                & '7  => Tall grass'//char(10)//                      &
                & '8  => Desert'//char(10)//                          &
                & '9  => Tundra'//char(10)//                          &
                & '10 => Irrigated Crop'//char(10)//                  &
                & '11 => Semi-desert'//char(10)//                     &
                & '12 => Ice cap/glacier'//char(10)//                 &
                & '13 => Bog or marsh'//char(10)//                    &
                & '14 => Inland water'//char(10)//                    &
                & '15 => Ocean'//char(10)//                           &
                & '16 => Evergreen shrub'//char(10)//                 &
                & '17 => Deciduous shrub'//char(10)//                 &
                & '18 => Mixed Woodland'//char(10)//                  &
                & '19 => Forest/Field mosaic'//char(10)//             &
                & '20 => Water and Land mixture'//char(10)//          &
                & '21 => Urban'//char(10)//                           &
                & '22 => Sub-Urban')
        call check_ok(istatus,'Error adding landuse legend')
        istatus = nf90_put_att(ncid, ivar(4), 'standard_name',        &
                            &  'land_type')
        call check_ok(istatus,'Error adding landuse standard_name')
        istatus = nf90_put_att(ncid, ivar(4), 'long_name',            &
                      &  'Landuse category as defined in BATS1E')
        call check_ok(istatus,'Error adding landuse long_name')
        istatus = nf90_put_att(ncid, ivar(4), 'units', '1')
        call check_ok(istatus,'Error adding landuse units')
        istatus = nf90_put_att(ncid, ivar(4), 'coordinates',          &
                            &  'xlon xlat')
        call check_ok(istatus,'Error adding landuse coordinates')
!
        istatus = nf90_enddef(ncid)
        call check_ok(istatus,'Error End Definitions NetCDF output')
!
        istatus = nf90_put_var(ncid, izvar(1), sigma)
        call check_ok(istatus,'Error variable sigma write')
        hptop = ptop*10.0
        istatus = nf90_put_var(ncid, izvar(2), hptop)
        call check_ok(istatus,'Error variable ptop write')
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
        istatus = nf90_put_var(ncid, illvar(1), transpose(xlat))
        call check_ok(istatus,'Error variable xlat write')
        istatus = nf90_put_var(ncid, illvar(2), transpose(xlon))
        call check_ok(istatus,'Error variable xlon write')
        call zeit_co('newfile')

      end subroutine open_sstfile

      subroutine close_sstfile
        implicit none
        integer :: istatus
        istatus = nf90_close(ncid)
        call check_ok(istatus,'Error Closing output sst file')
      end subroutine close_sstfile

      subroutine writerec(idate,lice)
        use mod_date , only : idatediff
        implicit none
        integer , intent(in) :: idate
        logical , intent(in) :: lice
        integer :: istatus
        integer , dimension(1) :: istart1 , icount1
        integer , dimension(3) :: istart , icount
        real(dp) , dimension(1) :: xdate

        call zeit_ci('writerec')
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
        call check_ok(istatus,'Error variable time write')
        istatus = nf90_put_var(ncid, ivar(2), transpose(sstmm), &
                               istart, icount)
        call check_ok(istatus,'Error variable sst write')
        if (lice) then
          istatus = nf90_put_var(ncid, ivar(3), transpose(icemm), &
                                 istart, icount)
          call check_ok(istatus,'Error variable ice write')
        end if
        istatus = nf90_put_var(ncid, ivar(4), transpose(lu), &
                               istart, icount)
        call check_ok(istatus,'Error variable landuse write')
        itime = itime + 1
        call zeit_co('writerec')
      end subroutine writerec

      subroutine check_ok(ierr,message)
        use netcdf
        implicit none
        integer , intent(in) :: ierr
        character(*) :: message
        if (ierr /= nf90_noerr) then 
          call die('sst_grid',message,ierr,nf90_strerror(ierr),0)
        end if
      end subroutine check_ok

      end module mod_sst_grid
