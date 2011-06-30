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
  use mod_dynparam
  use mod_constants
  use mod_maps
  use mod_block
  use mod_smooth
  use mod_maputils
  use mod_interp
  use mod_fudge
  use mod_rdldtr
  use mod_write
  use mod_header
  use m_stdio
  use m_die
  use m_zeit
  use mod_memutil
  implicit none
!
  character(256) :: char_lnd , char_tex , char_lak
  character(256) :: namelistfile , prgname
  integer :: i , j , k , ierr , i0 , j0 , m , n
  logical :: ibndry
  real(dp) :: clong , hsum , have , dsinm
!
  data ibndry /.true./
  real(dp) :: jumpsize , apara , bpara , func , funcprev
  real(dp) , pointer , dimension(:) :: alph , dsig
  real(dp) :: psig , zsig , pstar
  integer :: icount
  integer , parameter :: maxiter = 1000000
  real(dp) , parameter :: conv_crit = 0.00001D0
!
  call header(1)
!
!     Read input global namelist
!
  call getarg(0, prgname)
  call getarg(1, namelistfile)
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
!
  if (debug_level > 2) then
    call zeit_ci('terrain')
  end if

  call prepare_grid(iy,jx,kz,ntex)

  if ( nsg>1 ) then
    call prepare_subgrid(iysg,jxsg,ntex)
  end if
!
!  Setup hardcoded sigma levels
!
  if ( kz==14 ) then                      ! RegCM2
    sigma(1) = 0.0
    sigma(2) = 0.04
    sigma(3) = 0.10
    sigma(4) = 0.17
    sigma(5) = 0.25
    sigma(6) = 0.35
    sigma(7) = 0.46
    sigma(8) = 0.56
    sigma(9) = 0.67
    sigma(10) = 0.77
    sigma(11) = 0.86
    sigma(12) = 0.93
    sigma(13) = 0.97
    sigma(14) = 0.99
    sigma(15) = 1.0
  else if ( kz==18 ) then                 ! RegCM3, default
    sigma(1) = 0.0
    sigma(2) = 0.05
    sigma(3) = 0.10
    sigma(4) = 0.16
    sigma(5) = 0.23
    sigma(6) = 0.31
    sigma(7) = 0.39
    sigma(8) = 0.47
    sigma(9) = 0.55
    sigma(10) = 0.63
    sigma(11) = 0.71
    sigma(12) = 0.78
    sigma(13) = 0.84
    sigma(14) = 0.89
    sigma(15) = 0.93
    sigma(16) = 0.96
    sigma(17) = 0.98
    sigma(18) = 0.99
    sigma(19) = 1.0
  else if ( kz==23 ) then                 ! MM5V3
    sigma(1) = 0.0
    sigma(2) = 0.05
    sigma(3) = 0.1
    sigma(4) = 0.15
    sigma(5) = 0.2
    sigma(6) = 0.25
    sigma(7) = 0.3
    sigma(8) = 0.35
    sigma(9) = 0.4
    sigma(10) = 0.45
    sigma(11) = 0.5
    sigma(12) = 0.55
    sigma(13) = 0.6
    sigma(14) = 0.65
    sigma(15) = 0.7
    sigma(16) = 0.75
    sigma(17) = 0.8
    sigma(18) = 0.85
    sigma(19) = 0.89
    sigma(20) = 0.93
    sigma(21) = 0.96
    sigma(22) = 0.98
    sigma(23) = 0.99
    sigma(24) = 1.0
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
    sigma(1) = d_zero
    pstar = stdpmb - d_10*ptop
    do k = 1 , kz-1
      sigma(k+1) = sigma(k)+dsig(k)
    end do
    sigma(kz+1) = d_one
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

    if (debug_level > 2) call zeit_ci('subgrid')

    write (stdout,*) 'Doing Subgrid with following parameters'
    write (stdout,*) 'ntypec = ' , ntypec_s
    write (stdout,*) 'iy     = ' , iysg
    write (stdout,*) 'jx     = ' , jxsg
    write (stdout,*) 'ds     = ' , ds/nsg
    write (stdout,*) 'clat   = ' , clat
    write (stdout,*) 'clon   = ' , clong
    write (stdout,*) 'iproj  = ' , iproj
    dsinm = (ds/dble(nsg))*d_1000
    write(stdout,*) 'Subgrid setup done'

    if ( iproj=='LAMCON' ) then
      call lambrt(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,0,xn,truelatl,truelath)
      call lambrt(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,1,xn,truelatl,truelath)
    else if ( iproj=='POLSTR' ) then
      call mappol(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,0)
      call mappol(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,1)
      xn = d_one
    else if ( iproj=='NORMER' ) then
      call normer(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,0)
      call normer(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                  clat,dsinm,1)
      xn = d_zero
    else if ( iproj=='ROTMER' ) then
      call rotmer(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clon,     &
                  clat,plon,plat,dsinm,0)
      call rotmer(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clon,     &
                  clat,plon,plat,dsinm,1)
      xn = d_zero
    else
      write (stderr,*) 'iproj = ', iproj
      write (stderr,*) 'Unrecognized or unsupported projection'
      write (stderr,*) 'Set iproj in LAMCON,POLSTR,NORMER,ROTMER'
      call die('terrain')
    end if
    write(stdout,*) 'Subgrid Geo mapping done'
!
!       reduce the search area for the domain
    call mxmnll(iysg,jxsg,xlon_s,xlat_s,i_band)
    write(stdout,*) 'Determined Subgrid coordinate range'
!
    if (debug_level > 2) call zeit_ci('TOPO read')
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'GTOPO_DEM_30s.nc','z',   &
                     30,ntypec_s,.true.,0)
    write(stdout,*)'Static DEM data successfully read in'
    call interp(iysg,jxsg,xlat_s,xlon_s,htgrid_s, &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec_s,3,lonwrap,lcrosstime)
    write(stdout,*)'Interpolated DEM on SUBGRID'
    if (debug_level > 2) call zeit_co('TOPO read')
!
    if (debug_level > 2) call zeit_ci('LAND read')
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'GLCC_BATS_30s.nc',       &
                     'landcover',30,ntypec_s,.true.,0)
    write(stdout,*)'Static landcover BATS data successfully read in'
    call interp(iysg,jxsg,xlat_s,xlon_s,lndout_s,  &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec_s,4,lonwrap,lcrosstime,     &
                ibnty=1,h2opct=h2opct)
    call filter1plakes(iysg,jxsg,lndout_s)
    write(stdout,*)'Interpolated landcover on SUBGRID'
    if (debug_level > 2) call zeit_co('LAND read')
!
    if ( aertyp(7:7)=='1' ) then
      if (debug_level > 2) call zeit_ci('SOIL read')
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'GLZB_SOIL_30s.nc',       &
                       'soiltype',30,ntypec_s,.true.,0)
      write(stdout,*)'Static texture data successfully read in'
      call interp(iysg,jxsg,xlat_s,xlon_s,texout_s,   &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec_s,4,lonwrap,lcrosstime,      &
                  ibnty=2,h2opct=h2opct)
      do i = 1 , ntex
        call interp(iysg,jxsg,xlat_s,xlon_s,frac_tex_s(:,:,i), &
                    nlatin,nlonin,grdltmn,grdlnmn,values,      &
                    ntypec_s,5,lonwrap,lcrosstime,ival=i)
      end do
      write(stdout,*)'Interpolated texture on SUBGRID'
      if (debug_level > 2) call zeit_co('SOIL read')
    end if

    if ( lakedpth ) then
      if (debug_level > 2) call zeit_ci('BATH read')
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'ETOPO_BTM_30s.nc',       &
                       'z',30,ntypec_s,.true.,0)
      write(stdout,*)'Static bathymetry data successfully read in'
      call interp(iysg,jxsg,xlat_s,xlon_s,dpth_s,    &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec_s,3,lonwrap,lcrosstime)
      write(stdout,*)'Interpolated bathymetry on SUBGRID'
      if (debug_level > 2) call zeit_co('BATH read')
    end if

!     ******           grell smoothing to eliminate 2 delx wave (6/90):
    call smth121(htgrid_s,iysg,jxsg)
    call smth121(htgrid_s,iysg,jxsg)

    if ( ibndry ) then
      do j = 2 , jxsg - 1
        htgrid_s(1,j) = htgrid_s(2,j)
        htgrid_s(iysg,j) = htgrid_s(iysg-1,j)
        lndout_s(1,j) = lndout_s(2,j)
        lndout_s(iysg,j) = lndout_s(iysg-1,j)
 
        if ( aertyp(7:7)=='1' ) then
          texout_s(1,j) = texout_s(2,j)
          texout_s(iysg,j) = texout_s(iysg-1,j)
          do k = 1 , ntex
            frac_tex_s(1,j,k) = frac_tex_s(2,j,k)
            frac_tex_s(iysg,j,k) = frac_tex_s(iysg-1,j,k)
          end do
        end if
      end do
      do i = 1 , iysg
        htgrid_s(i,1) = htgrid_s(i,2)
        htgrid_s(i,jxsg) = htgrid_s(i,jxsg-1)
        lndout_s(i,1) = lndout_s(i,2)
        lndout_s(i,jxsg) = lndout_s(i,jxsg-1)
 
        if ( aertyp(7:7)=='1' ) then
          texout_s(i,1) = texout_s(i,2)
          texout_s(i,jxsg) = texout_s(i,jxsg-1)
          do k = 1 , ntex
            frac_tex_s(i,1,k) = frac_tex_s(i,2,k)
            frac_tex_s(i,jxsg,k) = frac_tex_s(i,jxsg-1,k)
          end do
        end if
      end do
    end if
    do i = 1 , iysg
      do j = 1 , jxsg
        snowam_s(i,j) = 0.0
      end do
    end do

    write (char_lnd,99001) trim(dirter), pthsep, trim(domname), &
           '_LANDUSE' , nsg
    call lndfudge(fudge_lnd_s,lndout_s,htgrid_s,iysg,jxsg, &
                  trim(char_lnd))
    if ( aertyp(7:7)=='1' ) then
      write (char_tex,99001) trim(dirter), pthsep, trim(domname), &
             '_TEXTURE' , nsg
      call texfudge(fudge_tex_s,texout_s,htgrid_s,iysg,jxsg, &
                    trim(char_tex))
    end if
    write(stdout,*) 'Fudging data (if requested) succeeded'
    if (debug_level > 2) call zeit_co('subgrid')

  end if
!
!     set up the parameters and constants
!
  if (debug_level > 2) call zeit_ci('grid')
  write (stdout,*) 'Doing Grid with following parameters'
  write (stdout,*) 'ntypec = ' , ntypec_s
  write (stdout,*) 'iy     = ' , iysg
  write (stdout,*) 'jx     = ' , jxsg
  write (stdout,*) 'ds     = ' , ds/nsg
  write (stdout,*) 'clat   = ' , clat
  write (stdout,*) 'clon   = ' , clong
  write (stdout,*) 'iproj  = ' , iproj
  dsinm = ds*d_1000
  write(stdout,*) 'Grid setup done'
!
!-----calling the map projection subroutine
  if ( iproj=='LAMCON' ) then
    call lambrt(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0,xn,  &
                truelatl,truelath)
    call lambrt(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1,xn,  &
                truelatl,truelath)
  else if ( iproj=='POLSTR' ) then
    call mappol(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0)
    call mappol(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1)
    xn = d_one
  else if ( iproj=='NORMER' ) then
    call normer(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0)
    call normer(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1)
    xn = d_zero
  else if ( iproj=='ROTMER' ) then
    call rotmer(xlon,xlat,xmap,coriol,iy,jx,clong,clat,plon,plat,   &
                dsinm,0)
    call rotmer(dlon,dlat,dmap,coriol,iy,jx,clong,clat,plon,plat,   &
                dsinm,1)
    xn = d_zero
  else
    write (stderr,*) 'iproj = ', iproj
    write (stderr,*) 'Unrecognized or unsupported projection'
    write (stderr,*) 'Set iproj in LAMCON, POLSTR, NORMER, ROTMER'
    call die('terrain')
  end if
  write(stdout,*) 'Geo mapping done'
!
!     reduce the search area for the domain
  call mxmnll(iy,jx,xlon,xlat,i_band)
  write(stdout,*)'Determined Grid coordinate range'
!
  if (debug_level > 2) call zeit_ci('TOPO read')
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                   pthsep//'GTOPO_DEM_30s.nc','z',   &
                   30,ntypec,.true.,0)
  write(stdout,*)'Static DEM data successfully read in'
  call interp(iy,jx,xlat,xlon,htgrid,           &
              nlatin,nlonin,grdltmn,grdlnmn,values, &
              ntypec,3,lonwrap,lcrosstime)
  write(stdout,*)'Interpolated DEM on model GRID'
  if (debug_level > 2) call zeit_co('TOPO read')
!
  if (debug_level > 2) call zeit_ci('LAND read')
  call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                   pthsep//'GLCC_BATS_30s.nc',       &
                   'landcover',30,ntypec,.true.,0)
  write(stdout,*)'Static landcover BATS data successfully read in'
  call interp(iy,jx,xlat,xlon,lndout,            &
              nlatin,nlonin,grdltmn,grdlnmn,values, &
              ntypec,4,lonwrap,lcrosstime,       &
              ibnty=1,h2opct=h2opct)
  call filter1plakes(iy,jx,lndout)
  write(stdout,*)'Interpolated landcover on model GRID'
  if (debug_level > 2) call zeit_co('LAND read')
!
  if ( aertyp(7:7)=='1' ) then
    if (debug_level > 2) call zeit_ci('SOIL read')
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'GLZB_SOIL_30s.nc',       &
                     'soiltype',30,ntypec,.true.,0)
    write(stdout,*)'Static texture data successfully read in'
    call interp(iy,jx,xlat,xlon,texout,             &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec,4,lonwrap,lcrosstime,        &
                ibnty=2,h2opct=h2opct)
    do i = 1 , ntex
      call interp(iy,jx,xlat,xlon,frac_tex(:,:,i),   &
                  nlatin,nlonin,grdltmn,grdlnmn,values, &
                  ntypec,5,lonwrap,lcrosstime,ival=i)
    end do
    write(stdout,*)'Interpolated texture on model GRID'
    if (debug_level > 2) call zeit_co('SOIL read')
  end if

  if ( lakedpth ) then
    if (debug_level > 2) call zeit_ci('BATH read')
    call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                     pthsep//'ETOPO_BTM_30s.nc',       &
                     'z',30,ntypec,.true.,0)
    write(stdout,*)'Static bathymetry data successfully read in'
    call interp(iy,jx,xlat,xlon,dpth,              &
                nlatin,nlonin,grdltmn,grdlnmn,values, &
                ntypec,3,lonwrap,lcrosstime)
    write(stdout,*)'Interpolated bathymetry on model GRID'
    if (debug_level > 2) call zeit_co('BATH read')
  end if

!     ******           preliminary heavy smoothing of boundaries
  if ( smthbdy ) call smthtr(htgrid,iy,jx)
 
!     ******           grell smoothing to eliminate 2 delx wave (6/90):
  call smth121(htgrid,iy,jx)
  call smth121(htgrid,iy,jx)

  if ( ibndry ) then
    do j = 2 , jx - 1
      htgrid(1,j) = htgrid(2,j)
      htgrid(iy,j) = htgrid(iy-1,j)
      lndout(1,j) = lndout(2,j)
      lndout(iy,j) = lndout(iy-1,j)
 
      if ( aertyp(7:7)=='1' ) then
        texout(1,j) = texout(2,j)
        texout(iy,j) = texout(iy-1,j)
        do k = 1 , ntex
          frac_tex(1,j,k) = frac_tex(2,j,k)
          frac_tex(iy,j,k) = frac_tex(iy-1,j,k)
        end do
      end if
    end do
    do i = 1 , iy
      htgrid(i,1) = htgrid(i,2)
      htgrid(i,jx) = htgrid(i,jx-1)
      lndout(i,1) = lndout(i,2)
      lndout(i,jx) = lndout(i,jx-1)
 
      if ( aertyp(7:7)=='1' ) then
        texout(i,1) = texout(i,2)
        texout(i,jx) = texout(i,jx-1)
        do k = 1 , ntex
          frac_tex(i,1,k) = frac_tex(i,2,k)
          frac_tex(i,jx,k) = frac_tex(i,jx-1,k)
        end do
      end if
    end do
  end if
 
  do i = 1 , iy
    do j = 1 , jx
      snowam(i,j) = 0.0
    end do
  end do

  write (char_lnd,99002) trim(dirter), pthsep, trim(domname), &
           '_LANDUSE'
  call lndfudge(fudge_lnd,lndout,htgrid,iy,jx,trim(char_lnd))

  if ( aertyp(7:7)=='1' ) then
    write (char_tex,99002) trim(dirter), pthsep, trim(domname), &
             '_TEXTURE'
    call texfudge(fudge_tex,texout,htgrid,iy,jx,trim(char_tex))
  end if

  where ( lndout > 14.5 .and. lndout < 15.5 )
    htgrid = 0.0
  end where
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
    call lakfudge(fudge_lak,dpth,lndout,iy,jx,trim(char_lak))
  end if
  write(stdout,*) 'Fudging data (if requested) succeeded'
  if (debug_level > 2) call zeit_co('grid')

  if ( nsg>1 ) then
    do i = 1 , iy
      do j = 1 , jx
        i0 = (i-1)*nsg
        j0 = (j-1)*nsg
        hsum = d_zero
        do m = 1 , nsg
          do n = 1 , nsg
            if ( htgrid(i,j)<0.1 .and. &
                (lndout(i,j)>14.5 .and. lndout(i,j)<15.5) ) then
              htgrid_s(i0+m,j0+n) = 0.0
              lndout_s(i0+m,j0+n) = 15.0
            end if
            hsum = hsum + dble(htgrid_s(i0+m,j0+n))
          end do
        end do
        have = hsum/dble(nnsg)
        do m = 1 , nsg
          do n = 1 , nsg
            htgrid_s(i0+m,j0+n) = htgrid(i,j) + &
                                  (htgrid_s(i0+m,j0+n) - real(have))
          end do
        end do
      end do
    end do

    where ( lndout_s > 14.5 .and. lndout_s < 15.5 )
      htgrid_s = 0.0
    end where
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
      call lakfudge(fudge_lak_s,dpth_s,lndout_s,iysg,jxsg, &
                    trim(char_lak))
    end if

    call write_domain(.true.)
    write(stdout,*) 'Subgrid data written to output file'
  end if

  call write_domain(.false.)
  write(stdout,*) 'Grid data written to output file'

  call memory_destroy

  if (debug_level > 2) then
    call zeit_co('terrain')
    call zeit_flush(stdout)
  end if

  write(stdout,*)'Successfully completed terrain fields generation'
!
99001 format (a,a,a,a,i0.3)
99002 format (a,a,a,a)
end program terrain
