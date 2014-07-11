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
!      GTOPO 30s DEM in NetCDF format
!      GLCC V2 BATS  in NetCDF format
!      SOIL ZOBLER   in NetCDF format
!      ETOPO BATHYM  in NetCDF format
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
  use mod_snow

  implicit none
!
  character(len=256) :: char_lnd , char_tex , char_lak
  character(len=256) :: namelistfile , prgname , outname
  integer(ik4) :: i , j , k , ierr , i0 , j0 , m , n , ism
  logical :: ibndry
  real(rk8) :: clong , hsum , have , dsinm
  integer(ik4) :: ntypec , ntypec_s
!
  data ibndry /.true./
  real(rk8) :: jumpsize , apara , bpara , func , funcprev
  real(rk8) , pointer , dimension(:) :: alph , dsig
  real(rk8) , allocatable , dimension(:,:) :: tmptex
  real(rk8) :: psig , zsig , pstar , tswap
  integer(ik4) :: icount
  integer(ik4) , parameter :: maxiter = 1000000
  real(rk8) , parameter :: conv_crit = 0.00001D0
!
  call header(1)
!
!     Read input global namelist
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
!
  call prepare_grid(jx,iy,kz,ntex)

  if ( nsg>1 ) then
    call prepare_subgrid(jxsg,iysg,ntex)
  end if
  call setup_outvars
!
!  Setup hardcoded sigma levels
!
  if ( kz==14 ) then                      ! RegCM2
    sigma(1) = 0.0D0
    sigma(2) = 0.04D0
    sigma(3) = 0.10D0
    sigma(4) = 0.17D0
    sigma(5) = 0.25D0
    sigma(6) = 0.35D0
    sigma(7) = 0.46D0
    sigma(8) = 0.56D0
    sigma(9) = 0.67D0
    sigma(10) = 0.77D0
    sigma(11) = 0.86D0
    sigma(12) = 0.93D0
    sigma(13) = 0.97D0
    sigma(14) = 0.99D0
    sigma(15) = 1.0D0
  else if ( kz==18 ) then                 ! RegCM3, default
    sigma(1) = 0.0D0
    sigma(2) = 0.05D0
    sigma(3) = 0.10D0
    sigma(4) = 0.16D0
    sigma(5) = 0.23D0
    sigma(6) = 0.31D0
    sigma(7) = 0.39D0
    sigma(8) = 0.47D0
    sigma(9) = 0.55D0
    sigma(10) = 0.63D0
    sigma(11) = 0.71D0
    sigma(12) = 0.78D0
    sigma(13) = 0.84D0
    sigma(14) = 0.89D0
    sigma(15) = 0.93D0
    sigma(16) = 0.96D0
    sigma(17) = 0.98D0
    sigma(18) = 0.99D0
    sigma(19) = 1.0D0
  else if ( kz==23 ) then                 ! MM5V3
    sigma(1) = 0.0D0
    sigma(2) = 0.05D0
    sigma(3) = 0.1D0
    sigma(4) = 0.15D0
    sigma(5) = 0.2D0
    sigma(6) = 0.25D0
    sigma(7) = 0.3D0
    sigma(8) = 0.35D0
    sigma(9) = 0.4D0
    sigma(10) = 0.45D0
    sigma(11) = 0.5D0
    sigma(12) = 0.55D0
    sigma(13) = 0.6D0
    sigma(14) = 0.65D0
    sigma(15) = 0.7D0
    sigma(16) = 0.75D0
    sigma(17) = 0.8D0
    sigma(18) = 0.85D0
    sigma(19) = 0.89D0
    sigma(20) = 0.93D0
    sigma(21) = 0.96D0
    sigma(22) = 0.98D0
    sigma(23) = 0.99D0
    sigma(24) = 1.0D0
  else
    call getmem1d(alph,1,kz,'terrain:alph')
    call getmem1d(dsig,1,kz,'terrain:dsig')
    write (stdout,*) 'Creating a custom set of sigma levels : '
    if ( (dsmax*dble(kz)) < d_one ) then
      write (stderr,*) 'dsmax must be greater than ', d_one/dble(kz)
      write (stderr,*) 'or kz must be less than ', d_one/dsmax
      call die('terrain','Maximum resolution, dsmax, is too low.',1)
    end if
    if ( (dsmin*dble(kz)) >= d_one ) then
      write (stderr,*) 'dsmin must be less than ', d_one/dble(kz)
      write (stderr,*) 'or kz must be greater than ', d_one/dsmax
      call die('terrain','Minimum resolution, dsmin, is too large.',1)
    end if
    ! Do a function minimization to determine the a,b coefficients for the
    ! following equation:
    !
    !        dsig(i) = dsmax*a^(i-1)*b^(0.5*(i-2)*(i-1))
    !
    ! which is derived from the recursive relation:
    !
    !        dsig(i) = a(i)*dsig(i-1) with a(i) = b*a(i-1)
    !
    ! The function provides level spacings between dsmin and dsmax such that
    ! level thicknesses vary approximately exponentially up to the TOA and gives
    ! more levels toward the surface.
    !
    ! Set the initial conditions for the function minimization
    !
    jumpsize = 0.0015D0
    bpara = 0.99573D0
    apara = ((dsmin/dsmax)**(d_one/dble(kz-1)))*(bpara**(-d_half*dble(kz-2)))
    alph(1) = apara/bpara
    dsig(1) = dsmax
    do k = 2 , kz
      alph(k) = bpara*alph(k-1)
      dsig(k) = alph(k)*dsig(k-1)
    end do
    func = sum(dsig)-d_one

    ! Loop through the minimization until the convergence criterion is satisfied
    do icount = 1 , maxiter
      funcprev = func
      bpara = bpara + jumpsize
      if ( bpara < d_zero ) bpara = 1D-10
      apara = ((dsmin/dsmax)**(d_one/dble(kz-1)))*(bpara**(-d_half*dble(kz-2)))
      alph(1) = apara/bpara
      dsig(1) = dsmax
      do k = 2 , kz
        alph(k) = bpara*alph(k-1)
        dsig(k) = alph(k)*dsig(k-1)
      end do
      func = sum(dsig)-d_one
      ! If we overshot 0, then reduce the jump size and reverse course
      if ( func*funcprev < d_zero ) then
        jumpsize = -jumpsize/d_two
      else if ( abs(func) > abs(funcprev) ) then
        jumpsize = -jumpsize
      end if
      ! If we converged, then jump out of the loop
      if ( abs(func) < conv_crit ) then
        write (stdout,*) 'Convergence reached.'
        write (stdout,*) '#', apara, bpara, icount
        write (stdout,*) '#', dsmax, dsmin, sum(dsig)
        exit
      end if
      if ( icount == maxiter-1 ) then
        write (stderr,*) 'Failed to converge.'
        write (stderr,*) 'b,a,jumpsize,func,funcprev:', bpara, apara, &
                         jumpsize,func,funcprev
        call die('terrain','Error setting up custom sigma levels')
      end if
    end do
    sigma(1) = 0.0D0
    pstar = stdpmb - d_10*ptop
    do k = 1 , kz-1
      sigma(k+1) = sigma(k)+dsig(k)
    end do
    sigma(kz+1) = 1.0D0
    ! Write the levels out to the screen
    zsig = d_zero
    psig = pstar*sigma(kz+1) + d_10*ptop
    write (stdout,*) 'k        sigma       p(mb)           z(m)'
    write (stdout,*) '-----------------------------------------'
    write (stdout,5925) kz+1, sigma(kz+1), psig, zsig
    do k = kz, 1, -1
      psig = pstar*sigma(k) + d_10*ptop
      zsig = zsig+rgas*stdt*pstar*dsig(k)/(psig*egrav)
      write (stdout,5925) k, sigma(k), psig, zsig
    end do
5925  format(i3,5x,f8.3,5x,f8.2,5x,f10.2)
    call relmem1d(dsig)
    call relmem1d(alph)
  end if

!---------------------------------------------------------------------
!
  clong = clon
  if ( clong>180. ) clong = clong - 360.
  if ( clong<=-180. ) clong = clong + 360.

  if ( nsg>1 ) then

    write (stdout,*) 'Doing Subgrid with following parameters'
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
                     pthsep//'GTOPO_DEM_30s.nc','z',   &
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
  write (stdout,*) 'Doing Grid with following parameters'
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
                   pthsep//'GTOPO_DEM_30s.nc','z',   &
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

  if ( lakedpth ) then
    write (char_lak,99002) trim(dirter), pthsep, trim(domname), &
             '_LAK'
    call lakfudge(fudge_lak,dpth,lndout,jx,iy,trim(char_lak))
  end if
  write(stdout,*) 'Fudging data (if requested) succeeded'

  if ( nsg > 1 ) then
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

    write (outname,'(a,i0.3,a)') &
       trim(dirter)//pthsep//trim(domname)//'_DOMAIN',nsg,'.nc'
    call write_domain(outname,.true.,fudge_lnd_s,fudge_tex_s,fudge_lak_s,  &
                      ntypec_s,sigma,xlat_s,xlon_s,dlat_s,dlon_s,xmap_s, &
                      dmap_s,coriol_s,mask_s,htgrid_s,lndout_s,snowam_s,   &
                      dpth_s,texout_s,frac_tex_s)
    write(stdout,*) 'Subgrid data written to output file'
  end if

  call read_snow(snowam,jx,iy)

  write (outname,'(a,i0.3,a)') &
     trim(dirter)//pthsep//trim(domname)//'_DOMAIN',0,'.nc'
  call write_domain(outname,.false.,fudge_lnd,fudge_tex,fudge_lak,ntypec, &
                    sigma,xlat,xlon,dlat,dlon,xmap,dmap,coriol,mask,    &
                    htgrid,lndout,snowam,dpth,texout,frac_tex)
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

end program terrain
