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

      program aerosol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! EDGAR                                                              c
!                                                                    c
! LIOUSSE96                                                          c
!                                                                    c
! BOND                                                               c
!                                                                    c
! GEIA                                                               c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_date
      use netcdf

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      logical :: there
      integer :: i , j , nrec , ierr
      integer :: iyy , jxx , kzz
      integer :: istatus , incin , ncid , idimid , ivarid
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      real(4) , dimension(ilon,jlat) :: aer2
      real(4) , allocatable , dimension(:,:) :: aermm , xlat , xlon
      real(4) , allocatable , dimension(:,:) :: finmat
      real(4) , allocatable , dimension(:) :: sigma
      character(256) :: namelistfile, prgname , terfile , aerofile
      character(64) :: history , csdate , aerdesc
      integer , dimension(4) :: idims
      integer , dimension(7) :: ivar
      integer :: irefdate , imondate , year , month , day , hour , imon
      real(4) , dimension(2) :: trlat
      real(4) :: hptop
      real(4) , allocatable , dimension(:) :: yiy
      real(4) , allocatable , dimension(:) :: xjx
      integer , dimension(2) :: ivvar
      integer , dimension(2) :: illvar
      integer , dimension(2) :: izvar
      integer , dimension(8) :: tvals
      integer , dimension(1) :: istart1 , icount1
      integer , dimension(3) :: istart , icount
      real(8) , dimension(1) :: xdate
!
!     Read input global namelist
!
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr/=0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
!
      allocate(aermm(iy,jx))
      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))
      allocate(finmat(jx,iy))

      inquire (file=trim(inpglob)//'/AERGLOB/AEROSOL.dat',exist=there)
      if ( .not.there ) print * , 'AEROSOL.dat is not available' ,      &
                             &' under ',trim(inpglob),'/AERGLOB/'
      open (11,file=trim(inpglob)//'/AERGLOB/AEROSOL.dat',              &
          & form='unformatted',recl=ilon*jlat*ibyte,access='direct',    &
          & status='old',err=100)

      aerofile = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
      istatus = nf90_create(aerofile, nf90_clobber, ncid)
      call check_ok(istatus, &
                &   ('Error creating NetCDF output '//trim(aerofile)))

      terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      istatus = nf90_open(terfile, nf90_nowrite, incin)
      call check_ok(istatus, &
                &   ('Error Opening Domain '//trim(terfile)))

      istatus = nf90_inq_dimid(incin, "iy", idimid)
      call check_ok(istatus, 'Dimension iy missing')
      istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
      call check_ok(istatus, 'Dimension iy read error')
      istatus = nf90_inq_dimid(incin, "jx", idimid)
      call check_ok(istatus, 'Dimension jx missing')
      istatus = nf90_inquire_dimension(incin, idimid, len=jxx)
      call check_ok(istatus, 'Dimension jx read error')
      istatus = nf90_inq_dimid(incin, "kz", idimid)
      call check_ok(istatus, 'Dimension kz missing')
      istatus = nf90_inquire_dimension(incin, idimid, len=kzz)
      call check_ok(istatus, 'Dimension kz read error')
      if ( iyy/=iy .or. jxx/=jx .or. kzz/=kz+1) then
        print * , 'IMPROPER DIMENSION SPECIFICATION'
        print * , '  namelist   : ' , iy , jx , kz
        print * , '  DOMAIN     : ' , iyy , jxx , kzz-1
        stop 'Dimensions mismatch'
      end if
      allocate(sigma(kzp1))
      istatus = nf90_inq_varid(incin, "sigma", ivarid)
      call check_ok(istatus, 'Variable sigma missing')
      istatus = nf90_get_var(incin, ivarid, sigma)
      call check_ok(istatus, 'Variable sigma read error')
      istatus = nf90_inq_varid(incin, "xlat", ivarid)
      call check_ok(istatus, 'Variable xlat missing')
      istatus = nf90_get_var(incin, ivarid, finmat)
      call check_ok(istatus, 'Variable xlat read error')
      xlat = transpose(finmat)
      istatus = nf90_inq_varid(incin, "xlon", ivarid)
      call check_ok(istatus, 'Variable xlon missing')
      istatus = nf90_get_var(incin, ivarid, finmat)
      call check_ok(istatus, 'Variable xlon read error')
      xlon = transpose(finmat)
      istatus = nf90_close(incin)
      call check_ok(istatus, &
                 &  ('Error closing Domain file '//trim(terfile)))
      deallocate(finmat)
!
      istatus = nf90_put_att(ncid, nf90_global, 'title',  &
           & 'ICTP Regional Climatic model V4 Aerosol program output')
      call check_ok(istatus, 'Error adding global title')
      istatus = nf90_put_att(ncid, nf90_global, 'institution', &
               & 'ICTP')
      call check_ok(istatus, 'Error adding global institution')
      istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
               & 'CF-1.4')
      call check_ok(istatus, 'Error adding global Conventions')
      call date_and_time(values=tvals)
      write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
           tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
           tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
           ' : Created by RegCM aerosol program'
      istatus = nf90_put_att(ncid, nf90_global, 'history', history)
      call check_ok(istatus, 'Error adding global history')

      istatus = nf90_put_att(ncid, nf90_global, 'references', &
               & 'http://eforge.escience-lab.org/gf/project/regcm')
      call check_ok(istatus, 'Error adding global references')
      istatus = nf90_put_att(ncid, nf90_global, 'experiment', &
               & domname)
      call check_ok(istatus, 'Error adding global experiment')
      istatus = nf90_put_att(ncid, nf90_global, 'projection', iproj)
      call check_ok(istatus, 'Error adding global projection')

      istatus = nf90_put_att(ncid, nf90_global,   &
               &   'grid_size_in_meters', ds*1000.0)
      call check_ok(istatus, 'Error adding global gridsize')
      istatus = nf90_put_att(ncid, nf90_global,   &
               &   'latitude_of_projection_origin', clat)
      call check_ok(istatus, 'Error adding global clat')
      istatus = nf90_put_att(ncid, nf90_global,   &
               &   'longitude_of_projection_origin', clon)
      call check_ok(istatus, 'Error adding global clon')
      if (iproj == 'ROTMER') then
        istatus = nf90_put_att(ncid, nf90_global, &
                 &   'latitude_of_projection_pole', plat)
        call check_ok(istatus, 'Error adding global plat')
        istatus = nf90_put_att(ncid, nf90_global, &
                 &   'longitude_of_projection_pole', plon)
        call check_ok(istatus, 'Error adding global plon')
      else if (iproj == 'LAMCON') then
        trlat(1) = truelatl
        trlat(2) = truelath
        istatus = nf90_put_att(ncid, nf90_global, &
                 &   'standard_parallel', trlat)
        call check_ok(istatus, 'Error adding global truelat')
      end if
      if (aertyp == 'AER00D0') then
        aerdesc = 'Neither aerosol, nor dust used'
      else if (aertyp == 'AER01D0') then
        aerdesc = 'Biomass burning, SO2 + BC + OC, no dust.'
      else if (aertyp == 'AER10D0') then
        aerdesc = 'Fossil fuel burning, SO2 + BC + OC, no dust'
      else if (aertyp == 'AER11D0') then
        aerdesc = 'Fossil fuel + Biomass burning, SO2 + BC'// &
                  & ' + OC, no dust'
      else if (aertyp == 'AER00D1') then
        aerdesc = 'Dust only'
      else if (aertyp == 'AER01D1') then
        aerdesc = 'Biomass, SO2 + BC + OC, with dust'
      else if (aertyp == 'AER10D1') then
        aerdesc = 'Fossil fuel, SO2 + BC + OC, with dust'
      else if (aertyp == 'AER11D1') then
        aerdesc = 'Fossil fuel + Biomass burning, SO2 + BC'// &
                  & ' + OC, with dust'
      end if

      istatus = nf90_put_att(ncid, nf90_global,  &
                         &   'aerosol_type', aerdesc)
      call check_ok(istatus, 'Error adding global aerosol_type')
      istatus = nf90_def_dim(ncid, 'iy', iy, idims(2))
      call check_ok(istatus, 'Error creating dimension iy')
      istatus = nf90_def_dim(ncid, 'jx', jx, idims(1))
      call check_ok(istatus, 'Error creating dimension jx')
      istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
      call check_ok(istatus, 'Error creating dimension time')
      istatus = nf90_def_dim(ncid, 'kz', kz+1, idims(4))
      call check_ok(istatus, 'Error creating dimension kz')
      istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4),   &
                          &  izvar(1))
      call check_ok(istatus, 'Error adding variable sigma')
      istatus = nf90_put_att(ncid, izvar(1), 'standard_name',       &
                          &  'atmosphere_sigma_coordinate')    
      call check_ok(istatus, 'Error adding sigma standard_name')
      istatus = nf90_put_att(ncid, izvar(1), 'long_name',      &
                          &  'Sigma at model layer midpoints')
      call check_ok(istatus, 'Error adding sigma long_name')
      istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
      call check_ok(istatus, 'Error adding sigma units')
      istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
      call check_ok(istatus, 'Error adding sigma axis')
      istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
      call check_ok(istatus, 'Error adding sigma positive')
      istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                   &         'sigma: sigma ps: ps ptop: ptop')
      call check_ok(istatus, 'Error adding sigma formula_terms')
      istatus = nf90_def_var(ncid, 'ptop', nf90_float,         &
                         &   varid=izvar(2))
      call check_ok(istatus, 'Error adding variable ptop')
      istatus = nf90_put_att(ncid, izvar(2), 'standard_name',  &
                          &  'air_pressure')
      call check_ok(istatus, 'Error adding ptop standard_name')
      istatus = nf90_put_att(ncid, izvar(2), 'long_name',      &
                          &  'Pressure at model top')
      call check_ok(istatus, 'Error adding ptop long_name')
      istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
      call check_ok(istatus, 'Error adding ptop units')
      istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), &
                          &  ivvar(1))
      call check_ok(istatus, 'Error adding variable iy')
      istatus = nf90_put_att(ncid, ivvar(1), 'standard_name',  &
                          &  'projection_y_coordinate')
      call check_ok(istatus, 'Error adding iy standard_name')
      istatus = nf90_put_att(ncid, ivvar(1), 'long_name',      &
                          &  'y-coordinate in Cartesian system')
      call check_ok(istatus, 'Error adding iy long_name')
      istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
      call check_ok(istatus, 'Error adding iy units')
      istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), &
                          &  ivvar(2))
      call check_ok(istatus, 'Error adding variable jx')
      istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                          &  'projection_x_coordinate')
      call check_ok(istatus, 'Error adding jx standard_name')
      istatus = nf90_put_att(ncid, ivvar(2), 'long_name',    &
                          &  'x-coordinate in Cartesian system')
      call check_ok(istatus, 'Error adding jx long_name')
      istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
      call check_ok(istatus, 'Error adding jx units')
      istatus = nf90_def_var(ncid, 'xlat', nf90_float, idims(1:2),  &
                          &  illvar(1))
      call check_ok(istatus, 'Error adding variable xlat')
      istatus = nf90_put_att(ncid, illvar(1), 'standard_name', &
                          &  'latitude')
      call check_ok(istatus, 'Error adding xlat standard_name')
      istatus = nf90_put_att(ncid, illvar(1), 'long_name',     &
                          &  'Latitude at cross points')
      call check_ok(istatus, 'Error adding xlat long_name')
      istatus = nf90_put_att(ncid, illvar(1), 'units',         &
                          &  'degrees_north')
      call check_ok(istatus, 'Error adding xlat units')
      istatus = nf90_def_var(ncid, 'xlon', nf90_float, idims(1:2),  &
                          &  illvar(2))
      call check_ok(istatus, 'Error adding variable xlon')
      istatus = nf90_put_att(ncid, illvar(2), 'standard_name', &
                          &  'longitude')
      call check_ok(istatus, 'Error adding xlon standard_name')
      istatus = nf90_put_att(ncid, illvar(2), 'long_name',     &
                          &  'Longitude at cross points')
      call check_ok(istatus, 'Error adding xlon long_name')
      istatus = nf90_put_att(ncid, illvar(2), 'units',         &
                          &  'degrees_east')
      call check_ok(istatus, 'Error adding xlon units')
      istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3),  &
                          &  ivar(1))
      call check_ok(istatus, 'Error adding variable time')

      irefdate = globidate1
      call split_idate(irefdate, year, month, day, hour)
      write (csdate, '(i0.4,a)') year,'-01-16 00:00:00'
      irefdate = mkidate(year, 1, 16, 0)
      istatus = nf90_put_att(ncid, ivar(1), 'units', &
                     &   'hours since '//csdate)
      call check_ok(istatus, 'Error adding time units')
      istatus = nf90_def_var(ncid, 'so2', nf90_float, idims(1:2),  &
                          &  ivar(2))
      call check_ok(istatus, 'Error adding variable so2')
      istatus = nf90_put_att(ncid, ivar(2), 'standard_name', &
                          &  'atmosphere_mass_content_of_sulfur_dioxide')
      call check_ok(istatus, 'Error adding so2 standard_name')
      istatus = nf90_put_att(ncid, ivar(2), 'long_name',     &
                          &  'Anthropogenic SO2 emission, EDGAR')
      call check_ok(istatus, 'Error adding so2 long_name')
      istatus = nf90_put_att(ncid, ivar(2), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding so2 units')
      istatus = nf90_put_att(ncid, ivar(2), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding so2 coordinates')
      istatus = nf90_def_var(ncid, 'bc', nf90_float, idims(1:2),  &
                          &  ivar(3))
      call check_ok(istatus, 'Error adding variable bc')
      istatus = nf90_put_att(ncid, ivar(3), 'standard_name', &
             &  'atmosphere_mass_content_of_black_carbon_dry_aerosol')
      call check_ok(istatus, 'Error adding bc standard_name')
      istatus = nf90_put_att(ncid, ivar(3), 'long_name',     &
                          &  'Anthropogenic Black Carbon (BC), EDGAR')
      call check_ok(istatus, 'Error adding bc long_name')
      istatus = nf90_put_att(ncid, ivar(3), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding bc units')
      istatus = nf90_put_att(ncid, ivar(3), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding bc coordinates')
      istatus = nf90_def_var(ncid, 'oc', nf90_float, idims(1:2),  &
                          &  ivar(4))
      call check_ok(istatus, 'Error adding variable oc')
      istatus = nf90_put_att(ncid, ivar(4), 'standard_name', &
             &  'atmosphere_mass_content_of_organic_carbon_dry_aerosol')
      call check_ok(istatus, 'Error adding oc standard_name')
      istatus = nf90_put_att(ncid, ivar(4), 'long_name',     &
                          &  'Anthropogenic Organic Carbon (OC), EDGAR')
      call check_ok(istatus, 'Error adding oc long_name')
      istatus = nf90_put_att(ncid, ivar(4), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding oc units')
      istatus = nf90_put_att(ncid, ivar(4), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding oc coordinates')
      istatus = nf90_def_var(ncid, 'so2_monthly', nf90_float,           &
                          &  idims(1:3), ivar(5))
      call check_ok(istatus, 'Error adding variable so2_monthly')
      istatus = nf90_put_att(ncid, ivar(5), 'standard_name', &
                         &  'atmosphere_mass_content_of_sulfur_dioxide')
      call check_ok(istatus, 'Error adding so2_monthly standard_name')
      istatus = nf90_put_att(ncid, ivar(5), 'long_name',     &
                         &  'Anthropogenic SO2 emission monthly, EDGAR')
      call check_ok(istatus, 'Error adding so2_monthly long_name')
      istatus = nf90_put_att(ncid, ivar(5), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding so2_monthly units')
      istatus = nf90_put_att(ncid, ivar(5), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding so2_monthly coordinates')
      istatus = nf90_def_var(ncid, 'bc_monthly', nf90_float, idims(1:3),&
                          &  ivar(6))
      call check_ok(istatus, 'Error adding variable bc_monthly')
      istatus = nf90_put_att(ncid, ivar(6), 'standard_name', &
             &  'atmosphere_mass_content_of_black_carbon_dry_aerosol')
      call check_ok(istatus, 'Error adding bc_monthly standard_name')
      istatus = nf90_put_att(ncid, ivar(6), 'long_name',     &
                     &  'Anthropogenic Black Carbon (BC), LIOUSSE')
      call check_ok(istatus, 'Error adding bc_monthly long_name')
      istatus = nf90_put_att(ncid, ivar(6), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding bc_monthly units')
      istatus = nf90_put_att(ncid, ivar(6), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding bc_monthly coordinates')
      istatus = nf90_def_var(ncid, 'oc_monthly', nf90_float, idims(1:3),&
                          &  ivar(7))
      call check_ok(istatus, 'Error adding variable oc_monthly')
      istatus = nf90_put_att(ncid, ivar(7), 'standard_name', &
             &  'atmosphere_mass_content_of_organic_carbon_dry_aerosol')
      call check_ok(istatus, 'Error adding oc_monthly standard_name')
      istatus = nf90_put_att(ncid, ivar(7), 'long_name',     &
                       &  'Anthropogenic Organic Carbon (OC), LIOUSSE')
      call check_ok(istatus, 'Error adding oc_monthly long_name')
      istatus = nf90_put_att(ncid, ivar(7), 'units', 'kg m-2')
      call check_ok(istatus, 'Error adding oc_monthly units')
      istatus = nf90_put_att(ncid, ivar(7), 'coordinates', &
                          &  'xlon xlat')
      call check_ok(istatus, 'Error adding oc_monthly coordinates')
!
      istatus = nf90_enddef(ncid)
      call check_ok(istatus, 'Error End Definitions NetCDF output')
!
      istatus = nf90_put_var(ncid, izvar(1), sigma)
      call check_ok(istatus, 'Error variable sigma write')
      deallocate(sigma)
      hptop = ptop * 10.0
      istatus = nf90_put_var(ncid, izvar(2), hptop)
      call check_ok(istatus, 'Error variable ptop write')
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
      call check_ok(istatus, 'Error variable iy write')
      istatus = nf90_put_var(ncid, ivvar(2), xjx)
      call check_ok(istatus, 'Error variable jx write')
      deallocate(yiy)
      deallocate(xjx)
      istatus = nf90_put_var(ncid, illvar(1), transpose(xlat))
      call check_ok(istatus, 'Error variable xlat write')
      istatus = nf90_put_var(ncid, illvar(2), transpose(xlon))
      call check_ok(istatus, 'Error variable xlon write')
!
!     ******    SET UP LONGITUDES AND LATITUDES FOR AEROSOL DATA
!
      do i = 1 , ilon
        loni(i) = -179.5 + float(i-1)
      end do
      do j = 1 , jlat
        lati(j) = -89.5 + 1.*float(j-1)
      end do
! 
!     ****** ALL AEROSOL DATA, 1 Deg data, Climate value
! 
      do nrec = 1 , 3
        read (11,rec=nrec) aer2
        call bilinx(aer2,aermm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        istatus = nf90_put_var(ncid, ivar(1+nrec), transpose(aermm))
        call check_ok(istatus, 'Error variable write')
      end do

      imondate = irefdate
      do imon = 1 , 12
        istart1(1) = imon
        icount1(1) = 1
        xdate(1) = dble(idatediff(imondate,irefdate))
        istatus = nf90_put_var(ncid, ivar(1), xdate, istart1, icount1)
        call check_ok(istatus, 'Error variable time write')
        imondate = inextmon(imondate)
      end do

      nrec = 4
      do i = 1 , 3
        do imon = 1, 12
          read (11,rec=nrec) aer2
          nrec = nrec + 1
          call bilinx(aer2,aermm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          istart(3) = imon
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          istatus = nf90_put_var(ncid, ivar(4+i), transpose(aermm),     &
                                 istart, icount)
          call check_ok(istatus, 'Error monthly variable write')
        end do
      end do
 
      deallocate(aermm)
      deallocate(xlat)
      deallocate(xlon)

      istatus = nf90_close(ncid)
      call check_ok(istatus, 'Error Closing output sst file')

      print *, 'Successfully built aerosol data for domain'
      stop

 100  continue
      print * , 'ERROR OPENING AEROSOL FILE'
      stop

      contains
!
      subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,iy,jx, &
                      & nflds)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nflds , nlati , nloni
      real(4) , dimension(nloni,nlati,nflds) :: fin
      real(4) , dimension(nlati) :: lati
      real(4) , dimension(iy,jx) :: lato , lono
      real(4) , dimension(nloni) :: loni
      real(4) , dimension(iy,jx,nflds) :: fout
      intent (in) fin , iy , jx , lati , lato , loni , lono , nflds ,   &
                & nlati , nloni
      intent (out) fout
!
! Local variables
!
      real(4) :: bas , lon180 , p , q , xsum , xind , yind
      integer :: i , ip , ipp1 , j , jq , jqp1 , l
!
!
!     PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
!     BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF
!     GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
!     GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
!     TRAPPED POINT.. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES
!     IN BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
!     INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
!     THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!     INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
!     THE CROSS POINTS IN THE MM4 MODEL.
!
!     IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!     OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL
!     GRID. LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     P.........EAST-WEST WEIGHTING FACTOR.
!     Q.........NORTH-SOUTH WEIGHTING FACTOR.
!     IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!     IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID
 
!     POINT.
 
      do j = 1 , jx
        do i = 1 , iy
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jq = max0(jq,1)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon180 = lono(i,j)
          if ( lono(i,j)<-180. ) lon180 = lono(i,j) + 360.
          if ( lono(i,j)>180. ) lon180 = lono(i,j) - 360.
          xind = (((lon180-loni(1))/(loni(nloni)-loni(1)))              &
               & *float(nloni-1)) + 1.
          ip = int(xind)
          ip = max0(ip,1)
          ipp1 = min0(ip+1,nloni)
          p = xind - ip
 
          do l = 1 , nflds
            xsum = 0.0
            bas = 0.0
            if ( fin(ip,jq,l)<-9990.0 .and. fin(ipp1,jq,l)<-9990.0 .and.&
               & fin(ipp1,jqp1,l)<-9990.0 .and. fin(ip,jqp1,l)<-9990.0 )&
               & then
              fout(i,j,l) = -9999.
            else
              if ( fin(ip,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*(1.-p)*fin(ip,jq,l)
                bas = bas + (1.-q)*(1.-p)
              end if
              if ( fin(ipp1,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*p*fin(ipp1,jq,l)
                bas = bas + (1.-q)*p
              end if
              if ( fin(ipp1,jqp1,l)>-9990.0 ) then
                xsum = xsum + q*p*fin(ipp1,jqp1,l)
                bas = bas + q*p
              end if
              if ( fin(ip,jqp1,l)>-9990.0 ) then
                xsum = xsum + q*(1.-p)*fin(ip,jqp1,l)
                bas = bas + q*(1.-p)
              end if
              fout(i,j,l) = xsum/bas
            end if
          end do
        end do
 
      end do
 
      end subroutine bilinx

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

      end program aerosol
