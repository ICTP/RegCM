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

#ifdef PNETCDF
subroutine myabort
  use mod_stdio
  use mpi
  implicit none
  integer :: ierr
  write(stderr,*) ' Execution terminated because of runtime error'
  call mpi_abort(mpi_comm_self,1,ierr)
end subroutine myabort
#else
subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort
#endif

program terrain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   TERRAIN is the first component of the REGional Climate Modeling    !
!   (RegCM) system version 4.0 and used to access archived terrain     !
!   height and landuse charactistics data at regular latitude-         !
!   longititude intervals and interpolate to the mesoscale grid for    !
!   a specified map projection.                                        !
!                                                                      !
!                                     PWC group, Abdus Salam ICTP      !
!                                                   May. 27, 2006      !
!                                                                      !
!   The authors wish to acknowledge the authors of NCAR MM4/5          !
!   terrain codes, their works is our code base.                       !
!                                                                      !
!       MM4 terrain code: A. Mcnab, T. Tarbell and N. Seaman           !
!                         C. Larkin and N. Seaman  (before 1986)       !
!       MM5 terrain code: Yong-Run Guo and Sue Chen  10/21/1993        !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  In the present version of RegCM, all the variables, which directly
!  related to map factors, are just calculated and stored in TERRAIN
!  for later use.
!
!  This program reads terrain height data from :
!
!      GTOPO or GMTED 30s DEM in NetCDF format
!      GLCC V2 BATS           in NetCDF format
!      SOIL ZOBLER            in NetCDF format
!      ETOPO BATHYM           in NetCDF format
!
!  and analyzes heights and landuse values to a given grid.
!
!---------------------------------------------------------------------
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date
  use mod_constants
  use mod_maps
  use mod_smooth
  use mod_maputils
  use mod_intldtr
  use mod_fudge
  use mod_rdldtr
  use mod_write
  use mod_header
  use mod_stdio
  use mod_message
  use mod_memutil
  use mod_moist
  use mod_sigma
  use mod_nhinterp
  use mod_earth
  use mod_stdatm
  use mod_zita
#ifdef PNETCDF
  use mpi
#endif

  implicit none
  character(len=256) :: char_lnd , char_tex , char_lak
  character(len=256) :: namelistfile , prgname , outname
  integer(ik4) :: i , j , k , ierr , ism , mmx
  integer(ik4) :: year , month , day , hour
  logical :: ibndry
  real(rkx) :: clong , dsnsg
  integer(ik4) :: ntypec , ntypec_s
  real(rkx) , allocatable , dimension(:,:) :: tmptex
  real(rkx) , pointer , dimension(:,:) :: values
  real(rkx) :: psig , psig1 , zsig , pstar , tswap , tsig
  real(rkx) :: ts0 , ptrop
  !type(globalfile) :: gfile
  character(len=*) , parameter :: f99001 = '(a,a,a,a,i0.3)'
  character(len=*) , parameter :: f99002 = '(a,a,a,a)'
  logical :: laround0 = .true.

  data ibndry /.false./

  ! You should nor modify those, but it may help with "difficult" domains
  ! to obtain better representation of topography and landuse

  ! RESAMPLING METHODS if lresamp === true
  ! 1 : mean of values
  ! 2 : median of values
  ! 3 : most present for classes
  ! 4 : mean of central 50% after median filtering

  integer(ik4) , parameter :: topo_resamp_method = 4
  integer(ik4) , parameter :: class_resamp_method = 3
  integer(ik4) , parameter :: moisture_resamp_method = 1

  ! INTERPOLATION METHODS
  ! 1 : bilinear interpolation
  ! 2 : bicubic interpolation
  ! 3 : nearest neighbour interpolation
  ! 4 : most represented class
  ! 5 : percent of most represented classes
  ! 6 : distance weigth with roidem radius in ds units

  integer(ik4) , parameter :: topo_interp_method = 6
  integer(ik4) , parameter :: class_interp_method = 4
  integer(ik4) , parameter :: percent_interp_method = 5
  integer(ik4) , parameter :: moist_interp_method = 1

  type(regcm_projection) :: pjx , pjd , pju , pjv
  type(anyprojparams) :: pjpara

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  call header(1)
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( ierr/=0 ) then
    write (stderr,*) 'Parameter initialization not completed'
    write (stderr,*) 'Usage : '
    write (stderr,*) '          ', trim(prgname), ' regcm.in'
    write (stderr,*) ' '
    write (stderr,*) 'Check argument and namelist syntax'
    call die('terrain')
  end if

  call memory_init
  call split_idate(globidate1,year,month,day,hour)

  call prepare_grid(jx,iy,kz,ntex,num_soil_layers,idynamic)
  if ( nsg>1 ) then
    call prepare_subgrid(jxsg,iysg,kz,ntex,num_soil_layers,idynamic)
  end if
  call setup_outvars

  if ( idynamic < 3 ) then
    call init_sigma(kz,dsmax,dsmin)
    sigma(:) = sigma_coordinate(:)
  else if ( idynamic == 3 ) then
    call model_zitaf(zita)
    sigma = sigmazita(zita)
    ak = md_ak(zita)
    bk = md_bk(zita)
  else
    write(stderr, *) 'UNKNOWN DYNAMICAL CORE'
  end if

  if ( idynamic == 2 ) then
    write(stdout, *) 'Using non hydrostatic parameters'
    write(stdout, '(a,f10.2)') ' base_state_pressure    = ', base_state_pressure
    write(stdout, '(a,f10.2)') ' logp_lrate             = ', logp_lrate
  end if

  if ( idynamic == 3 ) then
    write(stdout, *) 'Using non hydrostatic Moloch dynamic'
  else
    if ( iproj == 'ROTLLR' ) then
    write (stderr,*) 'Temporary fix:'
    write (stderr,*) 'Rotated Longitude/Latitude projection supported only&
                    & by MOLOCH'
    write (stderr,*) 'Either change iproj or use idynamic=3 in &
                    &coreparam (MOLOCH)'
    call die('terrain')
    end if
  end if

  clong = clon
  if ( clong>180. ) clong = clong - 360.
  if ( clong<=-180. ) clong = clong + 360.

  pjpara%pcode = iproj
  pjpara%clat = clat
  pjpara%clon = clon
  pjpara%plat = plat
  pjpara%plon = plon
  pjpara%trlat1 = truelatl
  pjpara%trlat2 = truelath
  pjpara%rotparam = .false.

  if ( nsg > 1 ) then
    write (stdout,*) ''
    write (stdout,*) 'Doing Horizontal Subgrid with following parameters'
    write (stdout,*) 'iy     = ' , iysg
    write (stdout,*) 'jx     = ' , jxsg
    write (stdout,*) 'kz     = ' , kz
    write (stdout,*) 'ds     = ' , ds/nsg
    write (stdout,*) 'clat   = ' , clat
    write (stdout,*) 'clon   = ' , clong
    write (stdout,*) 'iproj  = ' , iproj
    dsnsg = (ds/real(nsg,rkx))
    write(stdout,*) 'Subgrid setup done'

    pjpara%ds = dsnsg*1000.0_rk8
    pjpara%nlon = jxsg
    pjpara%nlat = iysg
    pjpara%staggerx = .false.
    pjpara%staggery = .false.
    call getcoord(pjpara,xlon_s,xlat_s,pjx,jxsg,iysg)
    if ( idynamic == 3 ) then
      pjpara%staggerx = .true.
      pjpara%staggery = .false.
      call getcoord(pjpara,ulon_s,ulat_s,pju,jxsg,iysg)
      pjpara%staggerx = .false.
      pjpara%staggery = .true.
      call getcoord(pjpara,vlon_s,vlat_s,pjv,jxsg,iysg)
      pjpara%staggerx = .true.
      pjpara%staggery = .true.
      call getcoord(pjpara,dlon_s,dlat_s,pjd,jxsg,iysg)
    else
      pjpara%staggerx = .true.
      pjpara%staggery = .true.
      call getcoord(pjpara,dlon_s,dlat_s,pjd,jxsg,iysg)
    end if

    if ( idynamic == 3 ) then
      call mappar(pjx,xlat_s,xlon_s,xmap_s)
      call mappar(pju,ulat_s,ulon_s,umap_s)
      call mappar(pjv,vlat_s,vlon_s,vmap_s)
      call corpar(xlat_s,coriol_s)
    else
      call mappar(pjx,xlat_s,xlon_s,xmap_s)
      call mappar(pjd,dlat_s,dlon_s,dmap_s)
      call corpar(dlat_s,coriol_s)
    end if

    write(stdout,*) 'Subgrid Geo mapping done'
    !
    ! reduce the search area for the domain
    !
    call mxmnll(jxsg,iysg,xlon_s,xlat_s,i_band)
    write(stdout,*) 'Determined Subgrid coordinate range'

    if ( lresamp ) then
      ntypec_s = int(dsnsg/6.0_rkx)
      if ( ntypec_s > 0 ) then
        do while ( mod(3600,ntypec_s*60) /= 0 .and. ntypec_s > 1 )
          ntypec_s = ntypec_s - 1
        end do
        ntypec_s = max(ntypec_s, 0)
      end if
    else
      ntypec_s = 0
    end if
    if ( ntypec_s > 0 ) then
      write(stdout,*) 'Using resampling at ', ntypec_s, ' minutes.'
    else
      write(stdout,*) 'No resampling used.'
    end if
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                     pthsep//trim(tersrc)//'_DEM_30s.nc','z',       &
                     ntypec_s,topo_resamp_method,i_band,            &
                     xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                     nlatin,nlonin,values)
    write(stdout,*)'Static DEM data successfully read in'
    call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,htgrid_s, &
                values,topo_interp_method,rdem=roidem)
    call relmem2d(values)
    write(stdout,*)'Interpolated DEM on SUBGRID'

    call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                     pthsep//'GLCC_BATS_30s.nc','landcover',        &
                     ntypec_s,class_resamp_method,i_band,           &
                     xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                     nlatin,nlonin,values)
    write(stdout,*)'Static landcover BATS data successfully read in'
    call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,lndout_s,values, &
                class_interp_method,ibnty=1,h2opct=h2opct,rdem=roidem)
    call filter1plakes(jxsg,iysg,lndout_s)
    call relmem2d(values)
    write(stdout,*)'Interpolated landcover on SUBGRID'

    if ( lsmoist ) then
      write(stdout, *) 'Reading ', trim(smsrc)//'-SOILMOISTURE.nc'
      if ( smsrc(1:3) == 'CPC' ) then
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                         pthsep//trim(smsrc)//'-SOILMOISTURE.nc',       &
                         'soilw',30,moisture_resamp_method,i_band,      &
                         xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                         nlatin,nlonin,values,month)
        where ( values > d_zero )
          ! The data are mm content in a column of 1.6 m
          values = values * 0.0016_rkx
        end where
      else if ( smsrc(1:6) == 'ERA20C' ) then
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                         pthsep//trim(smsrc)//'-SOILMOISTURE.nc',       &
                         'swvl1',30,moisture_resamp_method,i_band,      &
                         xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                         nlatin,nlonin,values,month)
      else
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                         pthsep//trim(smsrc)//'-SOILMOISTURE.nc',       &
                         'sm',15,moisture_resamp_method,i_band,         &
                         xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                         nlatin,nlonin,values,month)
        if ( smsrc(1:6) == 'ESACCI' ) then
          where ( values > d_zero )
            ! This is the scale factor from original data in the file
            values = values * 0.0001_rkx
          end where
        end if
      end if
      write(stdout,*)'Satellite soil moisture data successfully read in'
      mmx = (2*minval(shape(values))/2+1)**2
      laround0 = .true.
      do i = 1 , nlatin
        do j = 1 , nlonin
          if ( values(j,i) < d_zero ) then
            call findaround(values,i,j,nlatin,nlonin,mmx)
          end if
        end do
      end do
      call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,smoist_s,values, &
                  moist_interp_method)
      call relmem2d(values)
      write(stdout,*)'Interpolated soil moisture on SUBGRID'
    end if
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                     pthsep//'GLZB_SOIL_30s.nc','soiltype',         &
                     ntypec_s,class_resamp_method,i_band,           &
                     xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                     nlatin,nlonin,values)
    write(stdout,*)'Static texture data successfully read in'
    call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,texout_s,values, &
                class_interp_method,ibnty=2,h2opct=h2opct,rdem=roidem)
!$OMP PARALLEL DO
    do i = 1 , ntex
      call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,frac_tex_s(:,:,i),values, &
                  percent_interp_method,ival=i,rdem=roidem)
    end do
!$OMP END PARALLEL DO
    call relmem2d(values)
    write(stdout,*)'Interpolated texture on SUBGRID'

    if ( lakedpth ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'//              &
                       pthsep//'ETOPO_BTM_30s.nc','z',                &
                       ntypec_s,topo_resamp_method,i_band,            &
                       xlat_s,xlon_s,grdlnma,grdlnmn,grdltma,grdltmn, &
                       nlatin,nlonin,values)
      write(stdout,*)'Static bathymetry data successfully read in'
      call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,dpth_s,values, &
                  topo_interp_method,rdem=roidem)
      call relmem2d(values)
      write(stdout,*)'Interpolated bathymetry on SUBGRID'
    end if

    do i = 1 , iysg
      do j = 1 , jxsg
        snowam_s(j,i) = 0.0
        rmoist_s(j,i,:) = -1.0
      end do
    end do

    write (char_lnd,f99001) trim(dirter), pthsep, trim(domname), &
           '_LANDUSE' , nsg
    call lndfudge(fudge_lnd_s,lndout_s,jxsg,iysg,trim(char_lnd))
    write (char_tex,f99001) trim(dirter), pthsep, trim(domname), &
           '_TEXTURE' , nsg
    allocate(tmptex(jxsg,iysg))
    tmptex(:,:) = texout_s(:,:)
    call texfudge(fudge_tex_s,texout_s,lndout_s,jxsg,iysg,trim(char_tex))
    do i = 1 , iysg
      do j = 1 , jxsg
        ! Swap percentages of the old class and the new requested
        if ( texout_s(j,i) /= tmptex(j,i) ) then
          tswap = frac_tex_s(j,i,int(tmptex(j,i)))
          frac_tex_s(j,i,int(tmptex(j,i))) = &
                  frac_tex_s(j,i,int(texout_s(j,i)))
          frac_tex_s(j,i,int(texout_s(j,i))) = tswap
        end if
      end do
    end do
    deallocate(tmptex)
    write(stdout,*) 'Fudging data (if requested) succeeded'

  end if
  !
  ! set up the parameters and constants
  !
  write (stdout,*) ''
  write (stdout,*) 'Doing Horizontal Grid with following parameters'
  write (stdout,*) 'iy     = ' , iy
  write (stdout,*) 'jx     = ' , jx
  write (stdout,*) 'kz     = ' , kz
  write (stdout,*) 'ds     = ' , ds
  write (stdout,*) 'clat   = ' , clat
  write (stdout,*) 'clon   = ' , clong
  write (stdout,*) 'iproj  = ' , iproj
  write(stdout,*) 'Grid setup done'
  !
  ! calling the map projection subroutine
  !
  pjpara%ds = ds*1000.0_rk8
  pjpara%nlon = jx
  pjpara%nlat = iy
  pjpara%staggerx = .false.
  pjpara%staggery = .false.
  call getcoord(pjpara,xlon,xlat,pjx,jx,iy)
  if ( idynamic == 3 ) then
    pjpara%staggerx = .true.
    pjpara%staggery = .false.
    call getcoord(pjpara,ulon,ulat,pju,jx,iy)
    pjpara%staggerx = .false.
    pjpara%staggery = .true.
    call getcoord(pjpara,vlon,vlat,pjv,jx,iy)
    pjpara%staggerx = .true.
    pjpara%staggery = .true.
    call getcoord(pjpara,dlon,dlat,pjd,jx,iy)
  else
    pjpara%staggerx = .true.
    pjpara%staggery = .true.
    call getcoord(pjpara,dlon,dlat,pjd,jx,iy)
  end if

  if ( idynamic == 3 ) then
    call mappar(pjx,xlat,xlon,xmap)
    call mappar(pju,ulat,ulon,umap)
    call mappar(pjv,vlat,vlon,vmap)
    call corpar(xlat,coriol)
  else
    call mappar(pjx,xlat,xlon,xmap)
    call mappar(pjd,dlat,dlon,dmap)
    call corpar(dlat,coriol)
  end if

  write(stdout,*) 'Geo mapping done'

  !
  ! reduce the search area for the domain
  !
  call mxmnll(jx,iy,xlon,xlat,i_band)
  write(stdout,*)'Determined Grid coordinate range'

  if ( lresamp ) then
    ntypec = int(ds/6.0_rkx)
    if ( ntypec > 0 ) then
      do while ( mod(3600,ntypec*60) /= 0 .and. ntypec > 1 )
        ntypec = ntypec - 1
      end do
      ntypec = max(ntypec, 0)
    end if
  else
    ntypec = 0
  end if
  if ( ntypec > 0 ) then
    write(stdout,*) 'Using resampling at ', ntypec, ' minutes.'
  else
    write(stdout,*) 'No resampling used.'
  end if

  !call gfopen(gfile,trim(inpter)//pthsep//'SURFACE'//        &
  !            pthsep//trim(tersrc)//'_DEM_30s.nc',xlat,xlon, &
  !            ds,roidem,i_band)
  !call gfread(gfile,'z',htgrid,d_zero)
  !call gfclose(gfile)
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                   pthsep//trim(tersrc)//'_DEM_30s.nc','z',   &
                   ntypec,topo_resamp_method,i_band,          &
                   xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                   nlatin,nlonin,values)
  write(stdout,*)'Static DEM data successfully read in'
  call interp(ds,jx,iy,xlat,xlon,htgrid,values, &
              topo_interp_method,rdem=roidem)
  call relmem2d(values)
  write(stdout,*)'Interpolated DEM on model GRID'

  !call gfopen(gfile,trim(inpter)//pthsep//'SURFACE'//  &
  !            pthsep//'GLCC_BATS_30s.nc', xlat,xlon,   &
  !            ds,roidem,i_band)
  !call gfread(gfile,'landcover',lndout,15,h2opct,d_zero)
  !call gfclose(gfile)
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                   pthsep//'GLCC_BATS_30s.nc','landcover',    &
                   ntypec,class_resamp_method,i_band,         &
                   xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                   nlatin,nlonin,values)
  write(stdout,*)'Static landcover BATS data successfully read in'
  call interp(ds,jx,iy,xlat,xlon,lndout,values, &
              class_interp_method,ibnty=1,h2opct=h2opct,rdem=roidem)
  call filter1plakes(jx,iy,lndout)
  call relmem2d(values)
  write(stdout,*)'Interpolated landcover on model GRID'

  if ( lsmoist ) then
    write(stdout, *) 'Reading ', trim(smsrc)//'-SOILMOISTURE.nc'
    if ( smsrc(1:3) == 'CPC' ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                       pthsep//trim(smsrc)//'-SOILMOISTURE.nc',   &
                       'soilw',30,moisture_resamp_method,i_band,  &
                       xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                       nlatin,nlonin,values,month)
      where ( values > d_zero )
        ! The data are mm content in a column of 1.6 m
        values = values * 0.0016_rkx
      end where
    else if ( smsrc(1:6) == 'ERA20C' ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                       pthsep//trim(smsrc)//'-SOILMOISTURE.nc',   &
                       'swvl1',30,moisture_resamp_method,i_band,  &
                       xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                       nlatin,nlonin,values,month)
    else
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                       pthsep//trim(smsrc)//'-SOILMOISTURE.nc',   &
                       'sm',15,moisture_resamp_method,i_band,     &
                       xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                       nlatin,nlonin,values,month)
      if ( smsrc(1:6) == 'ESACCI' ) then
        where ( values > d_zero )
          ! This is the scale factor from original data in the file
          values = values * 0.0001_rkx
        end where
      end if
    end if
    write(stdout,*)'Satellite soil moisture data successfully read in'
    mmx = (2*minval(shape(values))/2+1)**2
    laround0 = .true.
    do i = 1 , nlatin
      do j = 1 , nlonin
        if ( values(j,i) < d_zero ) then
          call findaround(values,i,j,nlatin,nlonin,mmx)
        end if
      end do
    end do
    call interp(ds,jx,iy,xlat,xlon,smoist,values,moist_interp_method)
    call relmem2d(values)
    write(stdout,*)'Interpolated soil moisture on model GRID'
  end if
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'//          &
                   pthsep//'GLZB_SOIL_30s.nc','soiltype',     &
                   ntypec,class_resamp_method,i_band,         &
                   xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                   nlatin,nlonin,values)
  write(stdout,*)'Static texture data successfully read in'
  call interp(ds,jx,iy,xlat,xlon,texout,values, &
              class_interp_method,ibnty=2,h2opct=h2opct,rdem=roidem)
!$OMP PARALLEL DO
  do i = 1 , ntex
    call interp(ds,jx,iy,xlat,xlon,frac_tex(:,:,i),values, &
                percent_interp_method,ival=i,rdem=roidem)
  end do
!$OMP END PARALLEL DO
  call relmem2d(values)
  write(stdout,*)'Interpolated texture on model GRID'

  if ( lakedpth ) then
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'//        &
                   pthsep//'ETOPO_BTM_30s.nc','z',            &
                   ntypec,topo_resamp_method,i_band,          &
                   xlat,xlon,grdlnma,grdlnmn,grdltma,grdltmn, &
                   nlatin,nlonin,values)
    write(stdout,*)'Static bathymetry data successfully read in'
    call interp(ds,jx,iy,xlat,xlon,dpth,values, &
                topo_interp_method,rdem=roidem)
    call relmem2d(values)
    write(stdout,*)'Interpolated bathymetry on model GRID'
  end if

  do i = 1 , iy
    do j = 1 , jx
      snowam(j,i) = 0.0
      rmoist(j,i,:) = -1.0
    end do
  end do

  write (char_lnd,f99002) trim(dirter), pthsep, trim(domname),'_LANDUSE'
  call lndfudge(fudge_lnd,lndout,jx,iy,trim(char_lnd))
  write (char_tex,f99002) trim(dirter), pthsep, trim(domname),'_TEXTURE'
  allocate(tmptex(jx,iy))
  tmptex(:,:) = texout(:,:)
  call texfudge(fudge_tex,texout,lndout,jx,iy,trim(char_tex))
  do i = 1 , iy
    do j = 1 , jx
      ! Swap percentages of the old class and the new requested
      if ( texout(j,i) /= tmptex(j,i) ) then
        tswap = frac_tex(j,i,int(tmptex(j,i)))
        frac_tex(j,i,int(tmptex(j,i))) = frac_tex(j,i,int(texout(j,i)))
        frac_tex(j,i,int(texout(j,i))) = tswap
      end if
    end do
  end do
  deallocate(tmptex)

#ifdef CLM45
  if ( lclm45lake ) then
    where ( lndout > 13.5_rkx .and. lndout < 14.5_rkx )
      lndout = 8.0_rkx
    end where
  end if
#endif
  where ( lndout > 14.5_rkx .and. lndout < 15.5_rkx .and. htgrid < 0.0_rkx )
    htgrid = 0.0_rkx
  end where
  where ( lndout > 13.5_rkx .and. lndout < 15.5_rkx )
    mask = 0.0_rkx
  elsewhere
    mask = 2.0_rkx
  end where
  if (lakedpth) then
    where ( mask > 1.0 )
      dpth = 0.0_rkx
    end where
    where (mask < 1.0_rkx .and. dpth < 2.0_rkx)
      dpth = 2.0_rkx
    end where
  end if
  if ( lsmoist ) then
    where ( mask < 1.0_rkx )
      smoist = smissval
    end where
  else
    smoist = smissval
  end if

  ! Preliminary heavy smoothing of boundaries
  if ( smthbdy ) call smthtr(htgrid,jx,iy,nspgx)

  ! Grell smoothing to eliminate 2 delx wave
  call smtdsmt(htgrid,jx,iy)

  ! Smoothing using 1-2-1 smoother
  do ism = 1 , ismthlev
    call smth121(htgrid,jx,iy)
  end do

  if ( .not. h2ohgt ) then
    where ( lndout > 14.5_rkx .and. &
            lndout < 15.5_rkx .and. &
            htgrid > 0.0_rkx )
      htgrid = 0.0_rkx
    end where
  end if

  if ( ibndry ) then
    do j = 1 , jx
      htgrid(j,1) = htgrid(j,2)
      htgrid(j,iy-1) = htgrid(j,iy-2)
      htgrid(j,iy) = htgrid(j,iy-1)
      lndout(j,1) = lndout(j,2)
      lndout(j,iy-1) = lndout(j,iy-2)
      lndout(j,iy) = lndout(j,iy-1)
      mask(j,1) = mask(j,2)
      mask(j,iy-1) = mask(j,iy-2)
      mask(j,iy) = mask(j,iy-1)
      texout(j,1) = texout(j,2)
      texout(j,iy-1) = texout(j,iy-2)
      texout(j,iy) = texout(j,iy-1)
      do k = 1 , ntex
        frac_tex(j,1,k) = frac_tex(j,2,k)
        frac_tex(j,iy-1,k) = frac_tex(j,iy-2,k)
        frac_tex(j,iy,k) = frac_tex(j,iy-1,k)
      end do
    end do
    if ( i_band /= 1 ) then
      do i = 2 , iy-1
        htgrid(1,i) = htgrid(2,i)
        htgrid(jx-1,i) = htgrid(jx-2,i)
        htgrid(jx,i) = htgrid(jx-1,i)
        lndout(1,i) = lndout(2,i)
        lndout(jx-1,i) = lndout(jx-2,i)
        lndout(jx,i) = lndout(jx-1,i)
        mask(1,i) = mask(2,i)
        mask(jx-1,i) = mask(jx-2,i)
        mask(jx,i) = mask(jx-1,i)
        texout(1,i) = texout(2,i)
        texout(jx-1,i) = texout(jx-2,i)
        texout(jx,i) = texout(jx-1,i)
        do k = 1 , ntex
          frac_tex(1,i,k) = frac_tex(2,i,k)
          frac_tex(jx-1,i,k) = frac_tex(jx-2,i,k)
          frac_tex(jx,i,k) = frac_tex(jx-1,i,k)
        end do
      end do
    end if
  end if

  if ( lakedpth ) then
    write (char_lak,f99002) trim(dirter), pthsep, trim(domname), &
             '_LAK'
    call lakfudge(fudge_lak,dpth,lndout,jx,iy,trim(char_lak))
  end if
  write(stdout,*) 'Fudging data (if requested) succeeded'

  if ( nsg > 1 ) then
#ifdef CLM45
    if ( lclm45lake ) then
      where ( lndout_s > 13.5_rkx .and. lndout_s < 14.5_rkx )
        lndout_s = 8.0_rkx
      end where
    end if
#endif
    where ( lndout_s > 14.5_rkx .and. &
            lndout_s < 15.5_rkx .and. &
            htgrid_s < 0.0_rkx )
      htgrid_s = 0.0_rkx
    end where
    where ( lndout_s > 13.5_rkx .and. lndout_s < 15.5_rkx )
      mask_s = 0.0_rkx
    elsewhere
      mask_s = 2.0_rkx
    end where
    if (lakedpth) then
      where ( mask_s > 1.0_rkx )
        dpth_s = 0.0_rkx
      end where
      where ( mask_s < 1.0_rkx .and. dpth_s < 2.0_rkx )
        dpth_s = 2.0_rkx
      end where
    end if

    ! Grell smoothing to eliminate 2 delx wave
    call smtdsmt(htgrid_s,jxsg,iysg)

    ! Smoothing using 1-2-1 smoother
    do ism = 1 , ismthlev
      call smth121(htgrid_s,jxsg,iysg)
    end do

    if ( .not. h2ohgt ) then
      where ( lndout_s > 14.5_rkx .and. &
              lndout_s < 15.5_rkx .and. &
              htgrid_s > 0.0_rkx)
        htgrid_s = 0.0_rkx
      end where
    end if
    if ( lakedpth ) then
      write (char_lak,f99001) trim(dirter), pthsep, trim(domname), &
               '_LAK', nsg
      call lakfudge(fudge_lak_s,dpth_s,lndout_s,jxsg,iysg, &
                    trim(char_lak))
    end if
    if ( lsmoist ) then
      where ( mask_s < 1.0_rkx )
        smoist_s = smissval
      end where
    else
      smoist_s = smissval
    end if

    if ( ibndry ) then
      do j = 1 , jxsg
        htgrid_s(j,1) = htgrid_s(j,2)
        htgrid_s(j,iysg-1) = htgrid_s(j,iysg-2)
        htgrid_s(j,iysg) = htgrid_s(j,iysg-1)
        lndout_s(j,1) = lndout_s(j,2)
        lndout_s(j,iysg-1) = lndout_s(j,iysg-2)
        lndout_s(j,iysg) = lndout_s(j,iysg-1)
        mask_s(j,1) = mask_s(j,2)
        mask_s(j,iysg-1) = mask_s(j,iysg-2)
        mask_s(j,iysg) = mask_s(j,iysg-1)
        texout_s(j,1) = texout_s(j,2)
        texout_s(j,iysg-1) = texout_s(j,iysg-2)
        texout_s(j,iysg) = texout_s(j,iysg-1)
        do k = 1 , ntex
          frac_tex_s(j,1,k) = frac_tex_s(j,2,k)
          frac_tex_s(j,iysg-1,k) = frac_tex_s(j,iysg-2,k)
          frac_tex_s(j,iysg,k) = frac_tex_s(j,iysg-1,k)
        end do
      end do
      if ( i_band /= 1 ) then
        do i = 2 , iysg-1
          htgrid_s(1,i) = htgrid_s(2,i)
          htgrid_s(jxsg-1,i) = htgrid_s(jxsg-2,i)
          htgrid_s(jxsg,i) = htgrid_s(jxsg-1,i)
          lndout_s(1,i) = lndout_s(2,i)
          lndout_s(jxsg-1,i) = lndout_s(jxsg-2,i)
          lndout_s(jxsg,i) = lndout_s(jxsg-1,i)
          mask_s(1,i) = mask_s(2,i)
          mask_s(jxsg-1,i) = mask_s(jxsg-2,i)
          mask_s(jxsg,i) = mask_s(jxsg-1,i)
          texout_s(1,i) = texout_s(2,i)
          texout_s(jxsg-1,i) = texout_s(jxsg-2,i)
          texout_s(jxsg,i) = texout_s(jxsg-1,i)
          do k = 1 , ntex
            frac_tex_s(1,i,k) = frac_tex_s(2,i,k)
            frac_tex_s(jxsg-1,i,k) = frac_tex_s(jxsg-2,i,k)
            frac_tex_s(jxsg,i,k) = frac_tex_s(jxsg-1,i,k)
          end do
        end do
      end if
    end if

    if ( idynamic == 2 ) then
      ts0 = base_state_temperature(1,iysg,1,jxsg,xlat_s)
      call nhsetup(ptop,base_state_pressure,logp_lrate,ts0)
      call nhbase(1,iysg,1,jxsg,kz+1,sigma,htgrid_s, &
                  ps0_s,pr0_s,t0_s,rho0_s,z0_s)
    end if

    if ( idynamic == 3 ) then
      ! Here we compute zeta on full model levels.
      do k = 1 , kzp1
        do i = 1 , iysg
          do j = 1 , jxsg
            zeta_s(j,i,k) = md_zeta_h(zita(k),htgrid_s(j,i))
            fmz_s(j,i,k) = md_fmz_h(zita(k),htgrid_s(j,i))
          end do
        end do
      end do
    end if

    write (outname,'(a,i0.3,a)') &
       trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
    call write_domain(outname,.true.,fudge_lnd_s,fudge_tex_s,fudge_lak_s,   &
                      ntypec_s,sigma,xlat_s,xlon_s,dlat_s,dlon_s,ulat_s,    &
                      ulon_s,vlat_s,vlon_s,xmap_s,dmap_s,umap_s,vmap_s,     &
                      coriol_s,mask_s,htgrid_s,lndout_s,snowam_s,smoist_s,  &
                      rmoist_s,dpth_s,texout_s,frac_tex_s,ps0_s,pr0_s,t0_s, &
                      rho0_s,z0_s,ts0,zeta_s,fmz_s)
    write(stdout,*) 'Subgrid data written to output file'
  end if

  call read_moist(moist_filename,rmoist,snowam,jx,iy,num_soil_layers,lrmoist)

  if ( idynamic == 1 ) then
    ! Write the levels out to the screen
    write (stdout,*) 'Vertical Grid Description (T estimated)'
    write (stdout,*) ''
    write (stdout,*) '--------------------------------------------------'
    write (stdout,*) 'k        sigma       p(mb)           z(m)     T(K)'
    write (stdout,*) '--------------------------------------------------'
    pstar = stdpmb - d_10*ptop
    zsig = d_zero
    tsig = stdatm_val(xlat(jx/2,iy/2),stdpmb,istdatm_tempk)
    psig = pstar*sigma(kz+1) + d_10*ptop
    write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2,4x,f6.1)') &
             kz+1, sigma(kz+1), psig, zsig, tsig
    do k = kz, 1, -1
      psig = pstar*sigma(k) + d_10*ptop
      psig1 = pstar*d_half*(sigma(k)+sigma(k+1)) + d_10*ptop
      tsig = stdatm_val(xlat(jx/2,iy/2),psig1,istdatm_tempk)
      psig1 = pstar*sigma(k+1) + d_10*ptop
      zsig = zsig + rovg*tsig*log(psig1/psig)
      write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2,4x,f6.1)') &
             k, sigma(k), psig, zsig, tsig
    end do
  else if ( idynamic == 2 ) then
    ts0 = base_state_temperature(1,iy,1,jx,xlat)
    call nhsetup(ptop,base_state_pressure,logp_lrate,ts0)
    call nhbase(1,iy,1,jx,kz+1,sigma,htgrid,ps0,pr0,t0,rho0,z0)
    write (stdout,*) 'Vertical Grid Description (mean over domain)'
    write (stdout,*) ''
    write (stdout,*) '--------------------------------------------------'
    write (stdout,*) 'k        sigma       p(mb)           z(m)     T(K)'
    write (stdout,*) '--------------------------------------------------'
    do k = kzp1, 1, -1
      write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2,4x,f6.1)') k, sigma(k), &
        d_r100*sum(pr0(:,:,k))/real(jx*iy,rkx), &
        sum(z0(:,:,k))/real(jx*iy,rkx), sum(t0(:,:,k))/real(jx*iy,rkx)
    end do
  else if ( idynamic == 3 ) then
    do k = 1 , kzp1
      do i = 1 , iy
        do j = 1 , jx
          zeta(j,i,k) = md_zeta_h(zita(k),htgrid(j,i))
          fmz(j,i,k) = md_fmz_h(zita(k),htgrid(j,i))
        end do
      end do
    end do
    ! Write the levels out to the screen
    write (stdout,*) 'Vertical Grid Description (mean over domain)'
    write (stdout,*) ''
    write (stdout,*) '--------------------------------------------------'
    write (stdout,*) 'k        sigma       p(mb)          h(m)      T(K)'
    write (stdout,*) '--------------------------------------------------'
    ptrop = 17.0e3_rkx - 6.0e3_rkx*cos(sum(xlat)/real(jx*iy,rkx)*degrad)**2
    do k = kzp1, 1, -1
      zsig = sum(zeta(:,:,k))/real(jx*iy,rkx)
      if ( zsig > ptrop ) then
        if ( zsig > ptrop + 1000.0_rkx ) then
          tsig = (tzero - 56.5_rkx) + 0.00045_rkx * (zsig-ptrop+1000.0_rkx)
        else
          tsig = tzero - 56.5_rkx
        end if
      else
        tsig = stdt - lrate * zsig
      end if
      psig = stdp*(d_one - &
          0.0065/stdt*zsig)**((egrav*0.0289644)/(8.31432*0.0065))
      write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2,4x,f6.1)') k, sigma(k), &
        psig * d_r100 , zsig , tsig
    end do
    ptop = psig * d_r100 ! Approximation of top pressure. hPa
  end if

  write (outname,'(a,i0.3,a)') &
     trim(dirter)//pthsep//trim(domname)//'_DOMAIN',0,'.nc'
  call write_domain(outname,.false.,fudge_lnd,fudge_tex,fudge_lak,ntypec, &
                    sigma,xlat,xlon,dlat,dlon,ulat,ulon,vlat,vlon,xmap,   &
                    dmap,umap,vmap,coriol,mask,htgrid,lndout,snowam,      &
                    smoist,rmoist,dpth,texout,frac_tex,ps0,pr0,t0,rho0,   &
                    z0,ts0,zeta,fmz)
  write(stdout,*) 'Grid data written to output file'

  if ( debug_level > 2 ) then
    call viz_init(48,24)
    call viz_clear
    call viz_plot(mask)
    call viz_done
  end if

  call pjx%destruct( )
  call pjd%destruct( )
  call pju%destruct( )
  call pjv%destruct( )

  call memory_destroy

  write(stdout,*)'Successfully completed terrain fields generation'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

  contains

  subroutine findaround(xx,i,j,imax,jmax,mmx)
    implicit none
    integer(ik4) , intent(in) :: i , j , imax , jmax , mmx
    real(rkx) , dimension(:,:) , intent(inout) :: xx
    real(rk8) , dimension (mmx) :: vals
    integer(ik4) :: ii , jj , js , is , ip , il , maxil
    real(rk8) :: countcc
    real(rk8) , save :: mincc , maxcc
    if ( laround0 ) then
      countcc = 1.0_rk8
      maxcc = 0.0_rk8
      mincc = 0.0_rk8
      do ii = 1 , imax
        do jj = 1 , jmax
          if ( xx(jj,ii) > 0.0_rkx ) then
            countcc = countcc + 1.0_rk8
            if ( maxcc < xx(jj,ii) ) maxcc = xx(jj,ii)
            if ( mincc > xx(jj,ii) ) mincc = xx(jj,ii)
          end if
        end do
      end do
      if ( countcc > 0.0_rk8 ) then
        mincc = mincc / countcc
        maxcc = maxcc / countcc
      end if
      laround0 = .false.
    end if
    maxil = minval(shape(xx))/2
    il = 1
    do
      ip = 0
      vals(:) = d_zero
      do ii = i - il , i + il
        do jj = j - il , j + il
          is = ii
          js = jj
          if ( js < 1 ) js = 1-js
          if ( js > jmax ) js = 2*jmax - js
          if ( is < 1 ) is = 1-is
          if ( is > imax ) is = 2*imax - is
          if ( xx(js,is) > d_zero ) then
            ip = ip + 1
            vals(ip) = vals(ip) + xx(js,is)
          end if
        end do
      end do
      if ( ip > 0 ) then
        xx(j,i) = real(sum(vals(1:ip))/real(ip,rk8),rkx)
        exit
      else
        il = il + 1
        if ( il == maxil ) then
          write(stderr,*) 'At point lat = ',xlat(j,i)
          write(stderr,*) '         lon = ',xlon(j,i)
          write(stderr,*) 'Setting to default value '
          xx(j,i) = real((maxcc+mincc) * 0.5_rk8,rkx)
        end if
      end if
    end do
  end subroutine findaround

end program terrain
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
