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

subroutine myabort
  implicit none
  call abort
end subroutine myabort

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

  implicit none
  character(len=256) :: char_lnd , char_tex , char_lak
  character(len=256) :: namelistfile , prgname , outname
  integer(ik4) :: i , j , k , ierr , ism , mmx
  integer(ik4) :: year , month , day , hour
  logical :: ibndry
  real(rkx) :: clong , dsinm , dsnsg
  integer(ik4) :: ntypec , ntypec_s
  real(rkx) , allocatable , dimension(:,:) :: tmptex
  real(rkx) , pointer , dimension(:,:) :: values
  real(rkx) :: psig , zsig , pstar , tswap
  real(rkx) :: ts0
  !type(globalfile) :: gfile
  character(len=*) , parameter :: f99001 = '(a,a,a,a,i0.3)'
  character(len=*) , parameter :: f99002 = '(a,a,a,a)'
  logical :: laround0 = .true.

  data ibndry /.true./

  ! You should nor modify those, but it may help with "difficult" domains
  ! to obtain better representation of topography and landuse

  ! RESAMPLING METHODS if lresamp === true
  ! 1 : mean of values
  ! 2 : median of values
  ! 3 : most present for classes
  ! 4 : mean of central 50% after median filtering

  integer(ik4) :: topo_resamp_method = 4
  integer(ik4) :: class_resamp_method = 3
  integer(ik4) :: moisture_resamp_method = 1

  ! INTERPOLATION METHODS
  ! 1 : bilinear interpolation
  ! 2 : bicubic interpolation
  ! 3 : nearest neighbour interpolation
  ! 4 : most represented class
  ! 5 : percent of most represented classes
  ! 6 : distance weigth with roidem radius in ds units

  integer(ik4) :: topo_interp_method = 6
  integer(ik4) :: class_interp_method = 4
  integer(ik4) :: percent_interp_method = 5
  integer(ik4) :: moist_interp_method = 1

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

  call init_sigma(kz,dsmax,dsmin)
  sigma(:) = sigma_coordinate(:)

  if ( idynamic == 1 ) then
    ! Write the levels out to the screen
    write (stdout,*) 'Vertical Grid Description'
    write (stdout,*) ''
    pstar = stdpmb - d_10*ptop
    zsig = d_zero
    psig = pstar*sigma(kz+1) + d_10*ptop
    write (stdout,*) '-----------------------------------------'
    write (stdout,*) 'k        sigma       p(mb)           z(m)'
    write (stdout,*) '-----------------------------------------'
    write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2)') kz+1, sigma(kz+1), psig, zsig
    do k = kz, 1, -1
      psig = pstar*sigma(k) + d_10*ptop
      zsig = zsig+rgas*stdt*pstar*sigma_delta(k)/(psig*egrav)
      write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2)') k, sigma(k), psig, zsig
    end do
  else if ( idynamic == 2 ) then
    write(stdout, *) 'Using non hydrostatic parameters'
    write(stdout, '(a,f10.2)') ' base_state_pressure    = ', base_state_pressure
    write(stdout, '(a,f10.2)') ' logp_lrate             = ', logp_lrate
  end if

  clong = clon
  if ( clong>180. ) clong = clong - 360.
  if ( clong<=-180. ) clong = clong + 360.

  if ( nsg > 1 ) then

    write (stdout,*) ''
    write (stdout,*) 'Doing Horizontal Subgrid with following parameters'
    write (stdout,*) 'iy     = ' , iysg
    write (stdout,*) 'jx     = ' , jxsg
    write (stdout,*) 'ds     = ' , ds/nsg
    write (stdout,*) 'clat   = ' , clat
    write (stdout,*) 'clon   = ' , clong
    write (stdout,*) 'iproj  = ' , iproj
    dsnsg = (ds/real(nsg,rkx))
    dsinm = dsnsg*d_1000
    write(stdout,*) 'Subgrid setup done'

    if ( iproj=='LAMCON' ) then
      call lambrt(xlon_s,xlat_s,xmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,0,xcone,truelatl,truelath)
      call lambrt(dlon_s,dlat_s,dmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,1,xcone,truelatl,truelath)
    else if ( iproj=='POLSTR' ) then
      call mappol(xlon_s,xlat_s,xmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,0)
      call mappol(dlon_s,dlat_s,dmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,1)
      xcone = d_one
    else if ( iproj=='NORMER' ) then
      call normer(xlon_s,xlat_s,xmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,0)
      call normer(dlon_s,dlat_s,dmap_s,coriol_s,jxsg,iysg,clong,    &
                  clat,dsinm,1)
      xcone = d_zero
    else if ( iproj=='ROTMER' ) then
      call rotmer(xlon_s,xlat_s,xmap_s,coriol_s,jxsg,iysg,clon,     &
                  clat,plon,plat,dsinm,0)
      call rotmer(dlon_s,dlat_s,dmap_s,coriol_s,jxsg,iysg,clon,     &
                  clat,plon,plat,dsinm,1)
      xcone = d_zero
    else
      write (stderr,*) 'iproj = ', iproj
      write (stderr,*) 'Unrecognized or unsupported projection'
      write (stderr,*) 'Set iproj in LAMCON,POLSTR,NORMER,ROTMER'
      call die('terrain')
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
    do i = 1 , ntex
      call interp(dsnsg,jxsg,iysg,xlat_s,xlon_s,frac_tex_s(:,:,i),values, &
                  percent_interp_method,ival=i,rdem=roidem)
    end do
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
  write (stdout,*) 'iy     = ' , iysg
  write (stdout,*) 'jx     = ' , jxsg
  write (stdout,*) 'ds     = ' , ds
  write (stdout,*) 'clat   = ' , clat
  write (stdout,*) 'clon   = ' , clong
  write (stdout,*) 'iproj  = ' , iproj
  dsinm = ds*d_1000
  write(stdout,*) 'Grid setup done'
  !
  ! calling the map projection subroutine
  !
  if ( iproj=='LAMCON' ) then
    call lambrt(xlon,xlat,xmap,coriol,jx,iy,clong,clat,dsinm,0,xcone,  &
                truelatl,truelath)
    call lambrt(dlon,dlat,dmap,coriol,jx,iy,clong,clat,dsinm,1,xcone,  &
                truelatl,truelath)
  else if ( iproj=='POLSTR' ) then
    call mappol(xlon,xlat,xmap,coriol,jx,iy,clong,clat,dsinm,0)
    call mappol(dlon,dlat,dmap,coriol,jx,iy,clong,clat,dsinm,1)
    xcone = d_one
  else if ( iproj=='NORMER' ) then
    call normer(xlon,xlat,xmap,coriol,jx,iy,clong,clat,dsinm,0)
    call normer(dlon,dlat,dmap,coriol,jx,iy,clong,clat,dsinm,1)
    xcone = d_zero
  else if ( iproj=='ROTMER' ) then
    call rotmer(xlon,xlat,xmap,coriol,jx,iy,clong,clat,plon,plat,   &
                dsinm,0)
    call rotmer(dlon,dlat,dmap,coriol,jx,iy,clong,clat,plon,plat,   &
                dsinm,1)
    xcone = d_zero
  else
    write (stderr,*) 'iproj = ', iproj
    write (stderr,*) 'Unrecognized or unsupported projection'
    write (stderr,*) 'Set iproj in LAMCON, POLSTR, NORMER, ROTMER'
    call die('terrain')
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
  do i = 1 , ntex
    call interp(ds,jx,iy,xlat,xlon,frac_tex(:,:,i),values, &
                percent_interp_method,ival=i,rdem=roidem)
  end do
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

  where ( htgrid < 0.0 )
    htgrid = 0.0
  end where
  where ( lndout > 13.5 .and. lndout < 15.5 )
    mask = 0.0
  elsewhere
    mask = 2.0
  end where
  where ( lndout > 14.5 .and. lndout < 15.5 )
    htgrid = 0.0
  end where
  ! Correction for Black Sea
  where ( xlat > 40.0 .and. xlat < 47.2 .and. &
          xlon > 26.5 .and. xlon < 42.0 .and. &
          lndout > 13.5 .and. lndout < 14.5 )
    htgrid = 0.0
  end where
  if (lakedpth) then
    where ( mask > 1.0 )
      dpth = 0.0
    end where
    where (mask < 1.0 .and. dpth < 2.0)
      dpth = 2.0
    end where
  end if
  if ( lsmoist ) then
    where ( mask == 0 )
      smoist = smissval
    end where
  else
    smoist = smissval
  end if

  if ( h2ohgt ) then
    call h2oelev(jx,iy,htgrid,mask)
  end if

  ! preliminary heavy smoothing of boundaries
  if ( smthbdy ) call smthtr(htgrid,jx,iy,nspgx)

  ! grell smoothing to eliminate 2 delx wave (6/90):
  call smtdsmt(htgrid,jx,iy)
  do ism = 2 , ismthlev
    call smth121(htgrid,jx,iy)
  end do

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
    where ( htgrid_s < 0.0 )
      htgrid_s = 0.0
    end where
    where ( lndout_s > 13.5 .and. lndout_s < 15.5 )
      mask_s = 0.0
    elsewhere
      mask_s = 2.0
    end where
    where ( lndout_s > 14.5 .and. lndout_s < 15.5 )
      htgrid_s = 0.0
    end where
    ! Correction for Black Sea
    where ( xlat_s > 40.0 .and. xlat_s < 47.2 .and. &
            xlon_s > 26.5 .and. xlon_s < 42.0 .and. &
            lndout_s > 13.5 .and. lndout_s < 14.5 )
      htgrid_s = 0.0
    end where

    if ( h2ohgt ) then
      call h2oelev(jxsg,iysg,htgrid_s,mask_s)
    end if

    ! grell smoothing to eliminate 2 delx wave (6/90):
    call smtdsmt(htgrid_s,jxsg,iysg)
    do ism = 2 , ismthlev
      call smth121(htgrid_s,jxsg,iysg)
    end do

    if (lakedpth) then
      where ( mask_s > 1.0 )
        dpth_s = 0.0
      end where
      where ( mask_s < 1.0 .and. dpth_s < 2.0 )
        dpth_s = 2.0
      end where
    end if
    if ( lakedpth ) then
      write (char_lak,f99001) trim(dirter), pthsep, trim(domname), &
               '_LAK', nsg
      call lakfudge(fudge_lak_s,dpth_s,lndout_s,jxsg,iysg, &
                    trim(char_lak))
    end if
    if ( lsmoist ) then
      where ( mask_s == 0 )
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
        mask(j,1) = mask(j,2)
        mask(j,iysg-1) = mask(j,iysg-2)
        mask(j,iysg) = mask(j,iysg-1)
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
          mask(1,i) = mask(2,i)
          mask(jxsg-1,i) = mask(jxsg-2,i)
          mask(jxsg,i) = mask(jxsg-1,i)
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

    write (outname,'(a,i0.3,a)') &
       trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
    call write_domain(outname,.true.,fudge_lnd_s,fudge_tex_s,fudge_lak_s, &
                      ntypec_s,sigma,xlat_s,xlon_s,dlat_s,dlon_s,xmap_s,  &
                      dmap_s,coriol_s,mask_s,htgrid_s,lndout_s,snowam_s,  &
                      smoist_s,rmoist_s,dpth_s,texout_s,frac_tex_s,ps0_s, &
                      pr0_s,t0_s,rho0_s,z0_s,ts0)
    write(stdout,*) 'Subgrid data written to output file'
  end if

  call read_moist(moist_filename,rmoist,snowam,jx,iy,num_soil_layers,lrmoist)

  if ( idynamic == 2 ) then
    ts0 = base_state_temperature(1,iy,1,jx,xlat)
    call nhsetup(ptop,base_state_pressure,logp_lrate,ts0)
    call nhbase(1,iy,1,jx,kz+1,sigma,htgrid,ps0,pr0,t0,rho0,z0)
    write (stdout,*) 'Vertical Grid Description'
    write (stdout,*) ''
    write (stdout,*) '-----------------------------------------'
    write (stdout,*) 'k        sigma       p(mb)           z(m)'
    write (stdout,*) '-----------------------------------------'
    do k = kzp1, 1, -1
      write (stdout,'(i3,4x,f8.3,4x,f8.2,4x,f10.2)') k, sigma(k), &
        d_r100*sum(pr0(:,:,k))/real(jx*iy,rkx), sum(z0(:,:,k))/real(jx*iy,rkx)
    end do
  end if

  write (outname,'(a,i0.3,a)') &
     trim(dirter)//pthsep//trim(domname)//'_DOMAIN',0,'.nc'
  call write_domain(outname,.false.,fudge_lnd,fudge_tex,fudge_lak,ntypec, &
                    sigma,xlat,xlon,dlat,dlon,xmap,dmap,coriol,mask,      &
                    htgrid,lndout,snowam,smoist,rmoist,dpth,texout,       &
                    frac_tex,ps0,pr0,t0,rho0,z0,ts0)
  write(stdout,*) 'Grid data written to output file'

  if ( debug_level > 2 ) then
    call viz_init(48,24)
    call viz_clear
    call viz_plot(mask)
    call viz_done
  end if

  call memory_destroy

  write(stdout,*)'Successfully completed terrain fields generation'

  contains

  subroutine findaround(xx,i,j,imax,jmax,mmx)
    implicit none
    integer(ik4) , intent(in) :: i , j , imax , jmax , mmx
    real(rkx) , dimension(:,:) , intent(inout) :: xx
    real(rkx) , dimension (mmx) :: vals
    integer(ik4) :: ii , jj , js , is , ip , il , maxil
    real(rkx) :: countcc
    real(rkx) , save :: mincc , maxcc
    if ( laround0 ) then
      countcc = 1.0_rkx
      maxcc = 0.0_rkx
      mincc = 0.0_rkx
      do ii = 1 , imax
        do jj = 1 , jmax
          if ( xx(jj,ii) > 0.0_rkx ) then
            countcc = countcc + 1.0_rkx
            if ( maxcc < xx(jj,ii) ) maxcc = xx(jj,ii)
            if ( mincc > xx(jj,ii) ) mincc = xx(jj,ii)
          end if
        end do
      end do
      if ( countcc > 0.0_rkx ) then
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
        xx(j,i) = sum(vals(1:ip))/real(ip,rkx)
        exit
      else
        il = il + 1
        if ( il == maxil ) then
          write(stderr,*) 'At point lat = ',xlat(j,i)
          write(stderr,*) '         lon = ',xlon(j,i)
          write(stderr,*) 'Setting to default value '
          xx(j,i) = (maxcc+mincc) * d_half
        end if
      end if
    end do
  end subroutine findaround

end program terrain
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
