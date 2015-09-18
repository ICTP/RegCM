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
  use mod_block
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

  implicit none
  character(len=256) :: char_lnd , char_tex , char_lak
  character(len=256) :: namelistfile , prgname , outname
  integer(ik4) :: i , j , k , ierr , ism , mmx
  integer(ik4) :: year , month , day , hour
  logical :: ibndry
  real(rk8) :: clong , dsinm
  integer(ik4) :: ntypec , ntypec_s
  real(rk8) , allocatable , dimension(:,:) :: tmptex
  real(rk8) :: psig , zsig , pstar , tswap
  integer(ik4) :: ifupr
  real(rk8) :: logp_lrate = 50.0D0
  data ibndry /.true./

  namelist /nonhydroparam/ ifupr , logp_lrate

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
    open(ipunit, file=namelistfile, status='old', &
         action='read', iostat=ierr)
    rewind(ipunit)
    read(ipunit, nml=nonhydroparam, iostat=ierr)
    if ( ierr /= 0 ) then
      write(stdout, *) 'Using default non hydrostatic parameters'
      ierr = 0
    end if
    close(ipunit)
  end if

  clong = clon
  if ( clong>180. ) clong = clong - 360.
  if ( clong<=-180. ) clong = clong + 360.

  if ( nsg>1 ) then

    write (stdout,*) ''
    write (stdout,*) 'Doing Horizontal Subgrid with following parameters'
    write (stdout,*) 'iy     = ' , iysg
    write (stdout,*) 'jx     = ' , jxsg
    write (stdout,*) 'ds     = ' , ds/nsg
    write (stdout,*) 'clat   = ' , clat
    write (stdout,*) 'clon   = ' , clong
    write (stdout,*) 'iproj  = ' , iproj
    dsinm = (ds/dble(nsg))*d_1000
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
!       reduce the search area for the domain
    call mxmnll(jxsg,iysg,xlon_s,xlat_s,i_band)
    write(stdout,*) 'Determined Subgrid coordinate range'
!
    ntypec_s = idnint((ds/dble(nsg)/d_two)*60.0D0/110.0)
    do while ( mod(3600,ntypec_s*60) /= 0 )
      ntypec_s = ntypec_s -1
    end do
    write(stdout,*) 'Using resampling at ', ntypec_s, ' minutes.'
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//trim(tersrc)//'_DEM_30s.nc','z',   &
                     30,ntypec_s,.true.,2)
    write(stdout,*)'Static DEM data successfully read in'
    call interp(jxsg,iysg,xlat_s,xlon_s,htgrid_s, &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec_s,1,lonwrap,lcrosstime)
    call relmem2d(values)
    write(stdout,*)'Interpolated DEM on SUBGRID'
!
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'GLCC_BATS_30s.nc',       &
                     'landcover',30,ntypec_s,.true.,3)
    write(stdout,*)'Static landcover BATS data successfully read in'
    call interp(jxsg,iysg,xlat_s,xlon_s,lndout_s,  &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec_s,4,lonwrap,lcrosstime,     &
                ibnty=1,h2opct=h2opct)
    call filter1plakes(jxsg,iysg,lndout_s)
    call relmem2d(values)
    write(stdout,*)'Interpolated landcover on SUBGRID'
!
    if ( lsmoist ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'ESACCI-SOILMOISTURE.nc', &
                       'sm',900,15,.true.,2,month)
      write(stdout,*)'Satellite soil moisture data successfully read in'
      mmx = (2*minval(shape(values))/2+1)**2
      do i = 1 , nlatin
        do j = 1 , nlonin
          if ( values(j,i) < 0 ) then
            call findaround(values,i,j,nlatin,nlonin,mmx)
          end if
        end do
      end do
      call interp(jxsg,iysg,xlat_s,xlon_s,smoist_s,  &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  15,1,lonwrap,lcrosstime)
      call relmem2d(values)
      write(stdout,*)'Interpolated soil moisture on SUBGRID'
    end if
!
    if ( ltexture ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'GLZB_SOIL_30s.nc',       &
                       'soiltype',30,ntypec_s,.true.,3)
      write(stdout,*)'Static texture data successfully read in'
      call interp(jxsg,iysg,xlat_s,xlon_s,texout_s,   &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec_s,4,lonwrap,lcrosstime,      &
                  ibnty=2,h2opct=h2opct)
      do i = 1 , ntex
        call interp(jxsg,iysg,xlat_s,xlon_s,frac_tex_s(:,:,i), &
                    nlatin,nlonin,grdltmn,grdlnmn,values,      &
                    ntypec_s,5,lonwrap,lcrosstime,ival=i)
      end do
      call relmem2d(values)
      write(stdout,*)'Interpolated texture on SUBGRID'
    end if

    if ( lakedpth ) then
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'ETOPO_BTM_30s.nc',       &
                       'z',30,ntypec_s,.true.,2)
      write(stdout,*)'Static bathymetry data successfully read in'
      call interp(jxsg,iysg,xlat_s,xlon_s,dpth_s,    &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec_s,1,lonwrap,lcrosstime)
      call relmem2d(values)
      write(stdout,*)'Interpolated bathymetry on SUBGRID'
    end if

!     ******           grell smoothing to eliminate 2 delx wave (6/90):

    do ism = 1 , ismthlev
      call smth121(htgrid_s,jxsg,iysg)
    end do

    if ( ibndry ) then
      do j = 2 , jxsg - 1
        htgrid_s(j,1) = htgrid_s(j,2)
        htgrid_s(j,iysg) = htgrid_s(j,iysg-1)
        lndout_s(j,1) = lndout_s(j,2)
        lndout_s(j,iysg) = lndout_s(j,iysg-1)

        if ( ltexture ) then
          texout_s(j,1) = texout_s(j,2)
          texout_s(j,iysg) = texout_s(j,iysg-1)
          do k = 1 , ntex
            frac_tex_s(j,1,k) = frac_tex_s(j,2,k)
            frac_tex_s(j,iysg,k) = frac_tex_s(j,iysg-1,k)
          end do
        end if
      end do
      do i = 1 , iysg
        htgrid_s(1,i) = htgrid_s(2,i)
        htgrid_s(jxsg,i) = htgrid_s(jxsg-1,i)
        lndout_s(1,i) = lndout_s(2,i)
        lndout_s(jxsg,i) = lndout_s(jxsg-1,i)

        if ( ltexture ) then
          texout_s(1,i) = texout_s(2,i)
          texout_s(jxsg,i) = texout_s(jxsg-1,i)
          do k = 1 , ntex
            frac_tex_s(1,i,k) = frac_tex_s(2,i,k)
            frac_tex_s(jxsg,i,k) = frac_tex_s(jxsg-1,i,k)
          end do
        end if
      end do
    end if
    do i = 1 , iysg
      do j = 1 , jxsg
        snowam_s(j,i) = 0.0
        rmoist_s(j,i,:) = -1.0
      end do
    end do

    write (char_lnd,99001) trim(dirter), pthsep, trim(domname), &
           '_LANDUSE' , nsg
    call lndfudge(fudge_lnd_s,lndout_s,htgrid_s,jxsg,iysg, &
                  trim(char_lnd))
    if ( ltexture ) then
      write (char_tex,99001) trim(dirter), pthsep, trim(domname), &
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
    end if
    write(stdout,*) 'Fudging data (if requested) succeeded'

  end if
!
!     set up the parameters and constants
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
!-----calling the map projection subroutine
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
!     reduce the search area for the domain
  call mxmnll(jx,iy,xlon,xlat,i_band)
  write(stdout,*)'Determined Grid coordinate range'
!
  ntypec = idnint((ds/d_two)*60.0D0/110.0)
  do while ( mod(3600,ntypec*60) /= 0 .and. ntypec > 1 )
    ntypec = ntypec - 1
  end do
  write(stdout,*) 'Using resampling at ', ntypec, ' minutes.'
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                   pthsep//trim(tersrc)//'_DEM_30s.nc','z',   &
                   30,ntypec,.true.,2)
  write(stdout,*)'Static DEM data successfully read in'
  call interp(jx,iy,xlat,xlon,htgrid,           &
              nlatin,nlonin,grdltmn,grdlnmn,values, &
              ntypec,1,lonwrap,lcrosstime)
  call relmem2d(values)
  write(stdout,*)'Interpolated DEM on model GRID'
!
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                   pthsep//'GLCC_BATS_30s.nc',       &
                   'landcover',30,ntypec,.true.,3)
  write(stdout,*)'Static landcover BATS data successfully read in'
  call interp(jx,iy,xlat,xlon,lndout,            &
              nlatin,nlonin,grdltmn,grdlnmn,values, &
              ntypec,4,lonwrap,lcrosstime,       &
              ibnty=1,h2opct=h2opct)
  call filter1plakes(jx,iy,lndout)
  call relmem2d(values)
  write(stdout,*)'Interpolated landcover on model GRID'
!
  if ( lsmoist ) then
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'ESACCI-SOILMOISTURE.nc', &
                     'sm',900,15,.true.,2,month)
    write(stdout,*)'Satellite soil moisture data successfully read in'
    mmx = (2*minval(shape(values))/2+1)**2
    do i = 1 , nlatin
      do j = 1 , nlonin
        if ( values(j,i) < 0 ) then
          call findaround(values,i,j,nlatin,nlonin,mmx)
        end if
      end do
    end do
    call interp(jx,iy,xlat,xlon,smoist,               &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                15,1,lonwrap,lcrosstime)
    call relmem2d(values)
    write(stdout,*)'Interpolated soil moisture on model GRID'
  end if
!
  if ( ltexture ) then
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'GLZB_SOIL_30s.nc',       &
                     'soiltype',30,ntypec,.true.,3)
    write(stdout,*)'Static texture data successfully read in'
    call interp(jx,iy,xlat,xlon,texout,             &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec,4,lonwrap,lcrosstime,        &
                ibnty=2,h2opct=h2opct)
    do i = 1 , ntex
      call interp(jx,iy,xlat,xlon,frac_tex(:,:,i),   &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec,5,lonwrap,lcrosstime,ival=i)
    end do
    call relmem2d(values)
    write(stdout,*)'Interpolated texture on model GRID'
  end if

  if ( lakedpth ) then
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'ETOPO_BTM_30s.nc',       &
                     'z',30,ntypec,.true.,2)
    write(stdout,*)'Static bathymetry data successfully read in'
    call interp(jx,iy,xlat,xlon,dpth,              &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec,1,lonwrap,lcrosstime)
    call relmem2d(values)
    write(stdout,*)'Interpolated bathymetry on model GRID'
  end if

!     ******           preliminary heavy smoothing of boundaries
  if ( smthbdy ) call smthtr(htgrid,jx,iy)

!     ******           grell smoothing to eliminate 2 delx wave (6/90):
  do ism = 1 , ismthlev
    call smth121(htgrid,jx,iy)
  end do

  if ( ibndry ) then
    do j = 2 , jx - 1
      htgrid(j,1) = htgrid(j,2)
      htgrid(j,iy) = htgrid(j,iy-1)
      lndout(j,1) = lndout(j,2)
      lndout(j,iy) = lndout(j,iy-1)

      if ( ltexture ) then
        texout(j,1) = texout(j,2)
        texout(j,iy) = texout(j,iy-1)
        do k = 1 , ntex
          frac_tex(j,1,k) = frac_tex(j,2,k)
          frac_tex(j,iy,k) = frac_tex(j,iy-1,k)
        end do
      end if
    end do
    do i = 1 , iy
      htgrid(1,i) = htgrid(2,i)
      htgrid(jx,i) = htgrid(jx-1,i)
      lndout(1,i) = lndout(2,i)
      lndout(jx,i) = lndout(jx-1,i)

      if ( ltexture ) then
        texout(1,i) = texout(2,i)
        texout(jx,i) = texout(jx-1,i)
        do k = 1 , ntex
          frac_tex(1,i,k) = frac_tex(2,i,k)
          frac_tex(jx,i,k) = frac_tex(jx-1,i,k)
        end do
      end if
    end do
  end if

  do i = 1 , iy
    do j = 1 , jx
      snowam(j,i) = 0.0
      rmoist(j,i,:) = -1.0
    end do
  end do

  write (char_lnd,99002) trim(dirter), pthsep, trim(domname),'_LANDUSE'
  call lndfudge(fudge_lnd,lndout,htgrid,jx,iy,trim(char_lnd))

  if ( ltexture ) then
    write (char_tex,99002) trim(dirter), pthsep, trim(domname),'_TEXTURE'
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
  end if

  where ( htgrid < 0.0 )
    htgrid = 0.0
  end where
  if ( .not. h2ohgt ) then
    where ( lndout > 14.5 .and. lndout < 15.5 )
      htgrid = 0.0
    end where
  end if
  where ( lndout > 13.5 .and. lndout < 15.5 )
    mask = 0.0
  elsewhere
    mask = 2.0
    dpth = 0.0
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
    else where
      smoist = smoist * 0.0001
    end where
  else
    smoist = smissval
  end if

  if ( lakedpth ) then
    write (char_lak,99002) trim(dirter), pthsep, trim(domname), &
             '_LAK'
    call lakfudge(fudge_lak,dpth,lndout,jx,iy,trim(char_lak))
  end if
  write(stdout,*) 'Fudging data (if requested) succeeded'

  if ( nsg > 1 ) then
    where ( htgrid_s < 0.0 )
      htgrid_s = 0.0
    end where
    if ( .not. h2ohgt ) then
      where ( lndout_s > 14.5 .and. lndout_s < 15.5 )
        htgrid_s = 0.0
      end where
    end if
    where ( lndout_s > 13.5 .and. lndout_s < 15.5 )
      mask_s = 0.0
    elsewhere
      mask_s = 2.0
    end where
    if (lakedpth) then
      where ( mask_s > 1.0 )
        dpth_s = 0.0
      end where
      where ( mask_s < 1.0 .and. dpth_s < 2.0 )
        dpth_s = 2.0
      end where
    end if
    if ( lakedpth ) then
      write (char_lak,99001) trim(dirter), pthsep, trim(domname), &
               '_LAK', nsg
      call lakfudge(fudge_lak_s,dpth_s,lndout_s,jxsg,iysg, &
                    trim(char_lak))
    end if
    if ( lsmoist ) then
      where ( mask_s == 0 )
        smoist_s = smissval
      else where
        smoist_s = smoist_s * 0.0001
      end where
    else
      smoist_s = smissval
    end if

    if ( idynamic == 2 ) then
      call nhsetup(ptop,stdp,stdt,logp_lrate)
      call nhbase(1,iysg,1,jxsg,kz+1,sigma,htgrid_s,ps0_s,pr0_s,t0_s,rho0_s)
    end if

    write (outname,'(a,i0.3,a)') &
       trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
    call write_domain(outname,.true.,fudge_lnd_s,fudge_tex_s,fudge_lak_s, &
                      ntypec_s,sigma,xlat_s,xlon_s,dlat_s,dlon_s,xmap_s,  &
                      dmap_s,coriol_s,mask_s,htgrid_s,lndout_s,snowam_s,  &
                      smoist_s,rmoist_s,dpth_s,texout_s,frac_tex_s,ps0_s, &
                      pr0_s,t0_s,rho0_s)
    write(stdout,*) 'Subgrid data written to output file'
  end if

  call read_moist(moist_filename,rmoist,snowam,jx,iy,num_soil_layers,lrmoist)

  if ( idynamic == 2 ) then
    call nhsetup(ptop,stdp,stdt,logp_lrate)
    call nhbase(1,iy,1,jx,kz+1,sigma,htgrid,ps0,pr0,t0,rho0)
  end if

  write (outname,'(a,i0.3,a)') &
     trim(dirter)//pthsep//trim(domname)//'_DOMAIN',0,'.nc'
  call write_domain(outname,.false.,fudge_lnd,fudge_tex,fudge_lak,ntypec, &
                    sigma,xlat,xlon,dlat,dlon,xmap,dmap,coriol,mask,      &
                    htgrid,lndout,snowam,smoist,rmoist,dpth,texout,       &
                    frac_tex,ps0,pr0,t0,rho0)
  write(stdout,*) 'Grid data written to output file'

  if ( debug_level > 2 ) then
    call viz_init(48,24)
    call viz_clear
    call viz_plot(mask)
    call viz_done
  end if

  call memory_destroy

  write(stdout,*)'Successfully completed terrain fields generation'
!
99001 format (a,a,a,a,i0.3)
99002 format (a,a,a,a)

  contains

  subroutine findaround(xx,i,j,imax,jmax,mmx)
    implicit none
    integer(ik4) , intent(in) :: i , j , imax , jmax , mmx
    real(rk8) , dimension(:,:) , intent(inout) :: xx
    real(rk8) , dimension (mmx) :: vals
    integer(ik4) :: ii , jj , js , is , ip , il , maxil
    il = 1
    maxil = minval(shape(xx))/2
    do
      ip = 0
      vals(:) = 0.0D0
      do ii = i - il , i + il
        do jj = j - il , j + il
          is = ii
          js = jj
          if ( js < 1 ) js = 1-js
          if ( js > jmax ) js = 2*jmax - js
          if ( is < 1 ) is = 1-is
          if ( is > imax ) is = 2*imax - is
          if ( xx(js,is) > 0.0 ) then
            ip = ip + 1
            vals(ip) = vals(ip) + xx(js,is)
          end if
        end do
      end do
      if ( ip > 0 ) then
        xx(j,i) = sum(vals(:))/real(ip)
        exit
      else
        il = il + 1
        if ( il == maxil ) then
          write(stderr,*) 'At point lat = ',xlat(j,i)
          write(stderr,*) '         lon = ',xlon(j,i)
          call die(__FILE__,'Not finding anything around !',__LINE__)
          exit
        end if
      end if
    end do
  end subroutine findaround

end program terrain
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
