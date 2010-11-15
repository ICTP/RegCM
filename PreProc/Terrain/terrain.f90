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
      use mod_maps
      use mod_block
      use mod_smooth , only : smth121 , smthtr
      use mod_maputils , only : lambrt , mappol , normer , rotmer
      use mod_interp
      use mod_fudge
      use mod_rdldtr , only : read_ncglob
      use mod_write , only : setup , write_domain , dsinm
      use mod_header , only : header
      implicit none
!
! Local variables
!
      character(256) :: char_lnd , char_tex , char_lak
      character(256) :: namelistfile , prgname
      integer :: i , j , k , minsize , ierr , i0 , j0 , m , n
      logical :: ibndry
      real(4) :: clong , hsum , have
!
      data ibndry /.true./
!
      call header(1)
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
!     Preliminary consistency check to avoid I/O format errors
!
      minsize = (kz+1)+16
      if ((iym2*jxm2) .lt. minsize) then
        write (6, *) 'Please increase domain size.'
        write (6, *) 'Minsize (iy-2)*(jx-2) is ', minsize
        call abort
      end if

      call allocate_grid(iy,jx,kz,ntex)
      if ( nsg>1 ) call allocate_subgrid(iysg,jxsg,ntex)
!
!     Setup hardcoded sigma levels

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
        write (*,*) 'You vertical level number is not 14, 18, or 23'
        write (*,*) 'Please set your sigma parameters in OUTPUT'
        stop
      end if

!---------------------------------------------------------------------
!
      clong = clon
      if ( clong>180. ) clong = clong - 360.
      if ( clong<=-180. ) clong = clong + 360.

      if ( nsg>1 ) then

        call setup(iysg,jxsg,ntypec_s,iproj,ds/nsg,clat,clong)
        print * , 'Subgrid setup done'

        if ( iproj=='LAMCON' ) then
          call lambrt(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,0,xn,truelatl,truelath)
          call lambrt(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,1,xn,truelatl,truelath)
        else if ( iproj=='POLSTR' ) then
          call mappol(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,0)
          call mappol(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,1)
          xn = 1.
        else if ( iproj=='NORMER' ) then
          call normer(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,0)
          call normer(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,1)
          xn = 0.
        else if ( iproj=='ROTMER' ) then
          call rotmer(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clon,     &
                    & clat,plon,plat,dsinm,0)
          call rotmer(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clon,     &
                    & clat,plon,plat,dsinm,1)
          xn = 0.
        else
          write (6,*) 'iproj = ', iproj
          write (6,*) 'Unrecognized or unsupported projection'
          write (6,*) 'Set iproj to one in LAMCON,POLSTR,NORMER,ROTMER'
          stop
        end if
        print * , 'Subgrid Geo mapping done'
!
!       reduce the search area for the domain
        call mxmnll(iysg,jxsg,xlon_s,xlat_s,i_band)
        print * , 'Determined Subgrid coordinate range'
!
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                         pthsep//'GTOPO_DEM_30s.nc','z',   &
                         30,ntypec_s,.true.,0,ht)
        print *, 'Static DEM data successfully read in'
        call interp(iysg,jxsg,xlat_s,xlon_s,htgrid_s, &
                    nlatin,nlonin,grdltmn,grdlnmn,ht, &
                    ntypec_s,2,lonwrap,lcrosstime)
        print *, 'Interpolated DEM on SUBGRID'
        deallocate(ht)
!
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                         pthsep//'GLCC_BATS_30s.nc',       &
                         'landcover',30,ntypec_s,.true.,0,lnd)
        print *, 'Static landcover BATS data successfully read in'
        call interp(iysg,jxsg,xlat_s,xlon_s,lndout_s,  &
                    nlatin,nlonin,grdltmn,grdlnmn,lnd, &
                    ntypec_s,4,lonwrap,lcrosstime,     &
                    ibnty=1,h2opct=h2opct)
        call filter1plakes(iysg,jxsg,lndout_s)
        print *, 'Interpolated landcover on SUBGRID'
        deallocate(lnd)
!
        if ( aertyp(7:7)=='1' ) then
          call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                           pthsep//'GLZB_SOIL_30s.nc',       &
                           'soiltype',30,ntypec_s,.true.,0,text)
          print *, 'Static texture data successfully read in'
          call interp(iysg,jxsg,xlat_s,xlon_s,texout_s,   &
                      nlatin,nlonin,grdltmn,grdlnmn,text, &
                      ntypec_s,4,lonwrap,lcrosstime,      &
                      ibnty=2,h2opct=h2opct)
          do i = 1 , ntex
            call interp(iysg,jxsg,xlat_s,xlon_s,frac_tex_s(:,:,i), &
                        nlatin,nlonin,grdltmn,grdlnmn,text,         &
                        ntypec_s,5,lonwrap,lcrosstime,ival=i)
          end do
          print *, 'Interpolated texture on SUBGRID'
          deallocate(text)
        end if

        if ( lakedpth ) then
          call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                           pthsep//'ETOPO_BTM_30s.nc',       &
                           'z',30,ntypec_s,.true.,0,dpt)
          print *, 'Static bathymetry data successfully read in'
          call interp(iysg,jxsg,xlat_s,xlon_s,dpth_s,    &
                      nlatin,nlonin,grdltmn,grdlnmn,dpt, &
                      ntypec_s,2,lonwrap,lcrosstime)
          print *, 'Interpolated bathymetry on SUBGRID'
          deallocate(dpt)
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
           &   '_LANDUSE' , nsg
        call lndfudge(fudge_lnd_s,lndout_s,htgrid_s,iysg,jxsg, &
                    & trim(char_lnd))
        if ( aertyp(7:7)=='1' ) then
          write (char_tex,99001) trim(dirter), pthsep, trim(domname), &
             &   '_TEXTURE' , nsg
          call texfudge(fudge_tex_s,texout_s,htgrid_s,iysg,jxsg, &
                        trim(char_tex))
        end if
        print * , 'Fudging data (if requested) succeeded'

      end if
!
!     set up the parameters and constants
!
      call setup(iy,jx,ntypec,iproj,ds,clat,clong)
      print * , 'Grid setup done'
!
!-----calling the map projection subroutine
      if ( iproj=='LAMCON' ) then
        call lambrt(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0,xn,  &
                  & truelatl,truelath)
        call lambrt(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1,xn,  &
                  & truelatl,truelath)
      else if ( iproj=='POLSTR' ) then
        call mappol(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0)
        call mappol(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1)
        xn = 1.
      else if ( iproj=='NORMER' ) then
        call normer(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0)
        call normer(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1)
        xn = 0.
      else if ( iproj=='ROTMER' ) then
        call rotmer(xlon,xlat,xmap,coriol,iy,jx,clong,clat,plon,plat,   &
                  & dsinm,0)
        call rotmer(dlon,dlat,dmap,coriol,iy,jx,clong,clat,plon,plat,   &
                  & dsinm,1)
        xn = 0.
      else
        write (6,*) 'iproj = ', iproj
        write (6,*) 'Unrecognized or unsupported projection'
        write (6,*) 'Set iproj to one in LAMCON, POLSTR, NORMER, ROTMER'
        stop
      end if
      print * , 'Geo mapping done'
!
!     reduce the search area for the domain
      call mxmnll(iy,jx,xlon,xlat,i_band)
      print *, 'Determined Grid coordinate range'
!
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'GTOPO_DEM_30s.nc','z',   &
                       30,ntypec,.true.,0,ht)
      print *, 'Static DEM data successfully read in'
      call interp(iy,jx,xlat,xlon,htgrid,           &
                  nlatin,nlonin,grdltmn,grdlnmn,ht, &
                  ntypec,2,lonwrap,lcrosstime)
      print *, 'Interpolated DEM on model GRID'
      deallocate(ht)
!
      call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                       pthsep//'GLCC_BATS_30s.nc',       &
                       'landcover',30,ntypec,.true.,0,lnd)
      print *, 'Static landcover BATS data successfully read in'
      call interp(iy,jx,xlat,xlon,lndout,            &
                  nlatin,nlonin,grdltmn,grdlnmn,lnd, &
                  ntypec,4,lonwrap,lcrosstime,       &
                  ibnty=1,h2opct=h2opct)
      call filter1plakes(iy,jx,lndout)
      print *, 'Interpolated landcover on model GRID'
      deallocate(lnd)
!
      if ( aertyp(7:7)=='1' ) then
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                         pthsep//'GLZB_SOIL_30s.nc',       &
                         'soiltype',30,ntypec,.true.,0,text)
        print *, 'Static texture data successfully read in'
        call interp(iy,jx,xlat,xlon,texout,             &
                    nlatin,nlonin,grdltmn,grdlnmn,text, &
                    ntypec,4,lonwrap,lcrosstime,        &
                    ibnty=2,h2opct=h2opct)
        do i = 1 , ntex
          call interp(iy,jx,xlat,xlon,frac_tex(:,:,i),   &
                      nlatin,nlonin,grdltmn,grdlnmn,text, &
                      ntypec,5,lonwrap,lcrosstime,ival=i)
        end do
        print *, 'Interpolated texture on model GRID'
        deallocate(text)
      end if

      if ( lakedpth ) then
        call read_ncglob(trim(inpter)//pthsep//'SURFACE'// &
                         pthsep//'ETOPO_BTM_30s.nc',       &
                         'z',30,ntypec,.true.,0,dpt)
        print *, 'Static bathymetry data successfully read in'
        call interp(iy,jx,xlat,xlon,dpth,              &
                    nlatin,nlonin,grdltmn,grdlnmn,dpt, &
                    ntypec,2,lonwrap,lcrosstime)
        print *, 'Interpolated bathymetry on model GRID'
        deallocate(dpt)
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
           &   '_LANDUSE'
      call lndfudge(fudge_lnd,lndout,htgrid,iy,jx,trim(char_lnd))

      if ( aertyp(7:7)=='1' ) then
        write (char_tex,99002) trim(dirter), pthsep, trim(domname), &
             &   '_TEXTURE'
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
             &   '_LAK'
        call lakfudge(fudge_lak,dpth,lndout,iy,jx,trim(char_lak))
      end if
      print * , 'Fudging data (if requested) succeeded'

      if ( nsg>1 ) then
        do i = 1 , iy
          do j = 1 , jx
            i0 = (i-1)*nsg
            j0 = (j-1)*nsg
            hsum = 0.0
            do m = 1 , nsg
              do n = 1 , nsg
                if ( htgrid(i,j)<0.1 .and. &
                    (lndout(i,j)>14.5 .and. lndout(i,j)<15.5) ) then
                  htgrid_s(i0+m,j0+n) = 0.0
                  lndout_s(i0+m,j0+n) = 15.0
                end if
                hsum = hsum + htgrid_s(i0+m,j0+n)
              end do
            end do
            have = hsum/float(nnsg)
            do m = 1 , nsg
              do n = 1 , nsg
                htgrid_s(i0+m,j0+n) = htgrid(i,j) + &
                                      (htgrid_s(i0+m,j0+n) - have)
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
               &   '_LAK', nsg
          call lakfudge(fudge_lak_s,dpth_s,lndout_s,iysg,jxsg, &
                        trim(char_lak))
        end if

        call write_domain(.true.)
        print * , 'Subgrid data written to output file'
        call free_subgrid
      end if

      call write_domain(.false.)
      print * , 'Grid data written to output file'
      call free_grid

      print *, 'Successfully completed terrain fields generation'
!
99001 format (a,a,a,a,i0.3)
99002 format (a,a,a,a)
      end program terrain
