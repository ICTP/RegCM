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

      program clmproc
 
      use netcdf
      use mod_nclib
      use mod_dynparam
      use mod_read_domain
      use mod_param_clm
      use mod_date
      use mod_clm3grid

      implicit none
!
      real(4) , parameter :: vmisdat=-9999
      integer , parameter :: ndim = 3
      logical , parameter :: bvoc = .false.
!
! Local variables
!
      integer :: istatus , ncid , idum
      integer , dimension(4) :: idims
      integer , dimension(4) :: ivdims
      integer :: irefdate , imondate , ldim , ivar
      real(4) , dimension(2) :: trlat
      real(4) , allocatable , dimension(:) :: yiy
      real(4) , allocatable , dimension(:) :: xjx
      integer , dimension(8) :: tvals
      real(4) :: hptop , xmiss
      real(8) , dimension(1) :: xdate
      integer , dimension(2) :: ivvar
      integer , dimension(2) :: illvar
      integer , dimension(2) :: izvar
      integer , dimension(4) :: icount , istart
      integer , dimension(1) :: istart1 , icount1
      integer , dimension(3) :: iadim
      character(64) , dimension(nfld) :: lnam
      character(64) :: cdum
      logical :: there
      character(64) , dimension(nfld) :: units
      real(4) , dimension(3) :: varmax , varmin
      real(8) :: xhr
      real(4) :: offset , xscale , xlatmin , xlatmax , xlonmin , xlonmax
      real(4) :: perr , pmax
      real(4) , allocatable , dimension(:) :: glat , glon , zlat ,      &
               &                           zlev , zlon
      real(4) , allocatable , dimension(:,:) :: mpu
      real(4) , allocatable , dimension(:,:,:) :: regxyz
      real(4) , allocatable , dimension(:,:,:,:) :: regyxzt , zoom ,    &
                                         & dumw
      real(4) , allocatable , dimension(:,:) :: landmask , sandclay
      integer :: ipathdiv , ierr
      integer :: i , iz , it , j , k , l , kmax
      integer :: jotyp , idin , idout , ifield , ifld , imap
      integer :: idatex , iyr , imo , idy , ihr , julnc
      character(256) :: namelistfile , prgname
      character(256) :: inpfile , terfile , checkfile
      character(256) :: outfil_nc
      character(64) :: history , csdate , cldim
      integer , dimension(8) :: ilevs
!
      data ilevs /-1,-1,-1,-1,-1,-1,-1,-1/
      data xmiss /-9999.0/
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
      if ( nsg/=1 ) then
        write ( 6,* ) 'CLM does not work with subgridding enable.'
        write ( 6,* ) 'Please set nsg=1 in regcm.in'
        stop
      end if

      call allocate_domain
!
!     ** Get latitudes and longitudes from DOMAIN file
!
      terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      call read_domain(terfile)

      if ( clatx/=clat .or. clonx/=clon) then
        print * , 'DOMAIN file is inconsistent with regcm.in'
        print * , '  namelist       :  clat=' , clat , ' clon=' , clon
        print * , '  DOMAIN file    :  clat=' , clatx , ' clon=' , clonx
        stop 782
      end if
      if ( iprojx/=iproj ) then
        print * , 'DOMAIN file is inconsistent with regcm.in'
        print * , '  namelist       : iproj=' , iproj
        print * , '  DOMAIN file    : iproj=' , iprojx
        stop 783
      end if
 
!     ** Set output variables
      jotyp = 2
      xscale = 1.
      offset = 0.
 
!     ** Open Output checkfile in NetCDF format

      checkfile = trim(dirglob)//pthsep//trim(domname)//'_CLM3.nc'
      istatus = nf90_create(checkfile, nf90_clobber, ncid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error creating NetCDF output ', trim(checkfile)
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      istatus = nf90_put_att(ncid, nf90_global, 'title',  &
           & 'ICTP Regional Climatic model V4 clm2rcm program output')
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
      istatus = nf90_put_att(ncid, nf90_global, 'Conventions', &
               & 'None')
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error adding global Conventions'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      call date_and_time(values=tvals)
      write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
           tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
           tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
           ' : Created by RegCM aerosol program'
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
                          &  ivar)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error Variable time definition in NetCDF output'
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      irefdate = globidate1
      call split_idate(irefdate, iyr , imo , idy , ihr)
      irefdate = imonmiddle(mkidate(iyr, 1, 1, 0))
      call split_idate(irefdate, iyr , imo , idy , ihr)
      write (csdate, '(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') iyr, '-', imo, &
           & '-', idy, ' ', ihr, ':00:00 UTC'
      istatus = nf90_put_att(ncid, ivar, 'units', &
                     &   'hours since '//csdate)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error Variable time units attribute'
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      istatus = nf90_enddef(ncid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error End Definitions NetCDF output'
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      istatus = nf90_put_var(ncid, izvar(1), sigx)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error Variable sigma write in NetCDF output'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      hptop = ptop * 10.0
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

      ! Pre-allocate a twelve month for monthly vars
      imondate = irefdate
      do it = 1 , 12
        istart1(1) = it
        icount1(1) = 1
        xdate(1) = dble(idatediff(imondate,irefdate))
        istatus = nf90_put_var(ncid, ivar, xdate, istart1, icount1)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable time write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        imondate = inextmon(imondate)
        imondate = imonmiddle(imondate)
      end do

!     ** determine which files to create (emission factor map or not)
      call comp(ifield,bvoc)
 
!     ** Loop over the fields defined in clm.param
      do ifld = 1 , ifield
 
!       ** Open input and output files
!       **   Some files have more than one required variable. 
!       Therefore, **   the output file should only be opened once.
        inpfile = trim(inpglob)//infil(ifld)
        inquire (file=inpfile,exist=there)
        if ( .not.there ) then
          print * , 'CLM Input file does not exist: ', trim(inpfile)
          stop 'NON-EXISTENT FILE'
        end if
        if ( ifld==ipft .or. ifld==ilai .or. ifld==ilak .or.            &
           & ifld==iglc .or. ifld==iurb .or. ifld==isnd .or.            &
           & ifld==icol .or. ifld==ioro .or. ifld==iiso .or.            &
           & ifld==ifma .or. ifld==imbo .or. ifld==ibpin .or.           &
           & ifld==iapin ) then
!         ************************ CHANGED LINE ABOVE to include iiso
!         ************************
          print * , 'OPENING Input NetCDF FILE: ' , trim(inpfile)
          ierr = nf90_open(inpfile,nf90_nowrite,idin)
          if ( ierr/=nf90_noerr ) then
            write (6,*) 'Cannot open input file ', trim(inpfile)
            stop 'INPUT NOT READY'
          end if
          ipathdiv = scan(inpfile, pthsep, .true.)
          if ( ipathdiv/=0 ) then
            outfil_nc = trim(dirglob)//pthsep//trim(domname)//          &
                   &  '_RCM'//inpfile(ipathdiv+7:)
          else
            outfil_nc = trim(dirglob)//pthsep//trim(domname)//          &
                   & '_RCM'//inpfile(7:)
          endif
!         CALL FEXIST(outfil_nc)
          print * , 'OPENING Output NetCDF FILE: ' , trim(outfil_nc)
          call rcrecdf(outfil_nc,idout,varmin,varmax,3,ierr)
        end if
 
!       ** Setup RegCM3 grid variables
        call param(jx,iy,nlev(ifld),xlat,xlon,varmin,varmax,            &
                 & xlat1d,xlon1d,xlonmin,xlonmax,xlatmin,xlatmax,       &
                 & iadim,ndim)
 
!       ** Setup CLM3 grid variables, including subgrid region
        allocate(glon(nlon(ifld)))
        allocate(glat(nlat(ifld)))
 
        call clm3grid1(nlon(ifld),nlat(ifld),nlev(ifld),ntim(ifld),     &
                     & glon1(ifld),glon2(ifld),glat1(ifld),glat2(ifld), &
                     & xlonmin,xlonmax,xlatmin,xlatmax,glon,glat,istart,&
                     & icount)
 
        if ( ifld==isnd .or. ifld==icly ) then
          istart(4) = 1
          icount(4) = 1
        end if
 
        allocate(zoom(icount(1),icount(2),icount(3),icount(4)))
        allocate(zlon(icount(1)))
        allocate(zlat(icount(2)))
        allocate(zlev(icount(3)))
        allocate(landmask(icount(1),icount(2)))
!
        call clm3grid2(nlon(ifld),nlat(ifld),glon,glat,istart,          &
                     & icount,zlon,zlat,zlev)
!
        print *, 'Reading variables from input file'
!
!       ** Read in the variables.
!       In some cases, special reads need to be performed:
!       - Sand/clay fractions need a mapping variable
!       - Lakes, wetlands, soil color, and orography need a
!       180 degree longitiude shift.
!       - Soil color and Orography do not have landmasks (must be made)
        if ( ifld==isnd .or. ifld==icly ) then
          allocate(sandclay(ntim(ifld),nlev(ifld)))
          allocate(mpu(icount(1),icount(2)))
          call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),1,      &
                       & ntim(ifld),1,nlev(ifld),1,1,1,1,sandclay)
          print *, 'Read ', trim(lnam(ifld))
          call readcdfr4(idin,vnam_lm,cdum,cdum,istart(1),              &
                       & icount(1),istart(2),icount(2),1,1,1,1,landmask)
          call readcdfr4(idin,vnam_st,cdum,cdum,istart(1),              &
                       & icount(1),istart(2),icount(2),1,1,1,1,mpu)
          do j = 1 , icount(2)
            do i = 1 , icount(1)
              imap = nint(mpu(i,j))
              do k = 1 , icount(3)
                if ( imap>0 .and. landmask(i,j)>0.5 ) then
                  zoom(i,j,k,1) = sandclay(imap,k)
                else
                  zoom(i,j,k,1) = vmisdat
                end if
              end do
            end do
          end do
          deallocate(mpu)
          deallocate(sandclay)
          do k = 1 , icount(3)
            zlev(k) = glev_st(k)
          end do
          ntim(ifld) = 1
 
        else if ( ifld==iiso .or. ifld==ibpin .or. ifld==iapin .or.     &
                & ifld==imbo ) then
          call readcdfr4_iso(idin,vnam(ifld),lnam(ifld),units(ifld),    &
                           & istart(1),icount(1),istart(2),icount(2),   &
                           & istart(3),icount(3),istart(4),icount(4),   &
                           & zoom)
          print *, 'Read ', trim(lnam(ifld))
        else
          if ( ifld/=icol ) then
            call readcdfr4(idin,vnam_lm,lnam(ifld),units(ifld),         &
                 & istart(1),icount(1),istart(2),icount(2),1,1,1,1,     &
                 & landmask)
            print *, 'Read ', trim(lnam(ifld))
          end if
          call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),        &
                       & istart(1),icount(1),istart(2),icount(2),       &
                       & istart(3),icount(3),istart(4),icount(4),zoom)
          print *, 'Read ', trim(lnam(ifld))
        end if
 
        if ( ifld==icol .or. ifld==iiso .or. ifld==ibpin .or.           &
           & ifld==imbo .or. ifld==iapin ) then
          print *, 'Adjusting landmask'
          do j = 1 , icount(2)
            do i = 1 , icount(1)
              if ( zoom(i,j,1,1)>vmin(ifld) ) then
                landmask(i,j) = 1.0
              else
                landmask(i,j) = vmisdat
              end if
            end do
          end do
        end if

        print * , 'READ/WRITE: ' , vnam(ifld) , lnam(ifld) , units(ifld)
 
!       ** Set the non-land values to missing for interpolation purposes

        if ( ifld/=ioro ) call maskme(landmask,zoom,vmisdat,icount(1),  &
                                    & icount(2),icount(3),icount(4))
 
!       ** Interpolate data to RegCM3 grid

        allocate(regyxzt(iy,jx,nlev(ifld),ntim(ifld)))
        allocate(regxyz(jx,iy,nlev(ifld)))

        call bilinx4d(zoom,zlon,zlat,icount(1),icount(2),regyxzt,xlon,  &
                    & xlat,iy,jx,icount(3),icount(4),vmin(ifld),vmisdat)
 
!       ** Write the interpolated data to NetCDF for CLM and checkfile

        do l = 1 , ntim(ifld)
          idatex = 1900000000 + l*10000 + 1500
          call julian(idatex,julnc,iyr,imo,idy,ihr,xhr)
          if ( ifld==ipft ) then
            do j = 1 , iy
              do i = 1 , jx
                if ( regyxzt(j,i,1,l)>-99. ) then
                  do k = 1 , nlev(ifld)
                    regyxzt(j,i,k,l) = nint(regyxzt(j,i,k,l))
                  end do
                  perr = 100.
                  kmax = -1
                  pmax = -99.
                  do k = 1 , nlev(ifld)
                    perr = perr - regyxzt(j,i,k,l)
                    if ( regyxzt(j,i,k,l)>pmax ) then
                      pmax = regyxzt(j,i,k,l)
                      kmax = k
                    end if
                  end do
                  regyxzt(j,i,kmax,l) = regyxzt(j,i,kmax,l) + perr
!                 print*,i,j,perr,pmax,regyxzt(j,i,kmax,l)
                end if
              end do
            end do
          end if
          do k = 1 , nlev(ifld)
            do j = 1 , iy
              do i = 1 , jx
                if ( ifld==icol ) regyxzt(j,i,k,l)                      &
                   & = float(nint(regyxzt(j,i,k,l)))
                if ( regyxzt(j,i,k,l)>vmin(ifld) ) then
                  regxyz(i,j,k) = regyxzt(j,i,k,l)
                else
                  regxyz(i,j,k) = 0
                end if
              end do
            end do
          end do
 
          call writecdf(idout,vnam(ifld),regxyz,jx,iy,nlev(ifld),iadim, &
                      & xhr,lnam(ifld),units(ifld),xscale,offset,varmin,&
                      & varmax,xlat1d,xlon1d,zlev,0,vmisdat,jotyp)
        end do
 
        istatus = nf90_redef(ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Redef in file ', trim(checkfile)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (nlev(ifld) > 1) then
          ldim = -1
          do iz = 1 , 8
            if (ilevs(iz) < 0) then
              ldim = -iz
              exit
            end if
            if (nlev(ifld) == ilevs(iz)) then
              ldim = iz
              exit
            end if
          end do
          if (ldim < 0) then
            ldim = -ldim
            ilevs(ldim) = nlev(ifld)
            write (cldim,'(a,i0.3)') 'level_', nlev(ifld)
            istatus = nf90_def_dim(ncid, cldim, nlev(ifld), idum)
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error creating dimension ', trim(cldim)
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
        end if

        if (vnam(ifld)(1:8) == 'MONTHLY_') then
          vnam(ifld) = vnam(ifld)(9:)
        end if

        if (ntim(ifld) > 1) then
          if (nlev(ifld) > 1) then
            ivdims(1) = idims(1)
            ivdims(2) = idims(2)
            ivdims(3) = idims(4)+ldim
            ivdims(4) = idims(3)
            istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, &
                                   ivdims, ivar)
          else
            ivdims(1) = idims(1)
            ivdims(2) = idims(2)
            ivdims(3) = idims(3)
            istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, &
                                   ivdims(1:3), ivar)
          end if
        else
          if (nlev(ifld) > 1) then
            ivdims(1) = idims(1)
            ivdims(2) = idims(2)
            ivdims(3) = idims(4)+ldim
            istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, &
                               &   ivdims(1:3), ivar)
          else
            ivdims(1) = idims(1)
            ivdims(2) = idims(2)
            istatus = nf90_def_var(ncid, vnam(ifld), nf90_float, &
                               &   ivdims(1:2), ivar)
          end if
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating variable ', trim(vnam(ifld))
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        idum = len_trim(lnam(ifld))
        istatus = nf90_put_att(ncid, ivar, 'long_name', &
                             & lnam(ifld)(1:idum))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ', trim(vnam(ifld)), &
                      ' long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        idum = len_trim(units(ifld))
        istatus = nf90_put_att(ncid, ivar, 'units', & 
                             & units(ifld)(1:idum))
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ', trim(vnam(ifld)), &
                      ' units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar, '_FillValue', xmiss)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ', trim(vnam(ifld)), &
                      ' _FillValue attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(ncid, ivar, 'coordinates', &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ', trim(vnam(ifld)), &
                      ' coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_enddef(ncid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error End Definitions NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        allocate(dumw(jx,iy,nlev(ifld),ntim(ifld)))
        do j = 1 , iy
          do i = 1 , jx
             dumw(i,j,:,:) = regyxzt(j,i,:,:)
          end do
        end do
        istatus = nf90_put_var(ncid, ivar, dumw)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ', trim(vnam(ifld)), &
                      ' write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        deallocate(dumw)

!       ** Deallocate variables for next CLM3 field
        deallocate(glon)
        deallocate(glat)
        deallocate(zoom)
        deallocate(zlon)
        deallocate(zlat)
        deallocate(zlev)
        deallocate(landmask)
        deallocate(regyxzt)
        deallocate(regxyz)
 
      end do  ! End nfld loop

      call free_domain

!     ** Close files

      call clscdf(idin,ierr)
      call clscdf(idout,ierr)

      istatus = nf90_close(ncid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error Closing output file ', trim(checkfile)
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      print *, 'Successfully completed CLM preprocessing.'
 
      contains

      subroutine julian(idate,julnc,iyr,imo,idy,ihr,xhr)
 
      implicit none
!
! Dummy arguments
!
      integer :: idate , idy , ihr , imo , iyr , julnc
      real(8) :: xhr
      intent (out) xhr
      intent (inout) idate , idy , ihr , imo , iyr , julnc
!
! Local variables
!
      integer :: ileap , iyrm1 , j , julday
      integer , dimension(12) :: jprev , lenmon
!
      data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 ,&
         & 31/
 
      iyr = idate/1000000
      imo = idate/10000 - iyr*100
      idy = idate/100 - iyr*10000 - imo*100
      ihr = idate - idate/100*100
      ileap = mod(iyr,4)
      if ( ileap==0 ) then
        lenmon(2) = 29
      else
        lenmon(2) = 28
      end if
 
      if ( ihr>23 ) then
        ihr = ihr - 24
        idy = idy + 1
      end if
      if ( idy>lenmon(imo) ) then
        idy = 1
        imo = imo + 1
      end if
      if ( imo>12 ) then
        imo = 1
        iyr = iyr + 1
      end if
      idate = iyr*1000000 + imo*10000 + idy*100 + ihr
 
      iyrm1 = iyr - 1
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + lenmon(j-1)
      end do
 
      julday = idy + jprev(imo) - 1
 
      julnc = ((iyr-1900)*365+julday+int((iyrm1-1900)/4))*24 + ihr
      xhr = float(julnc)
 
      end subroutine julian

      end program clmproc
