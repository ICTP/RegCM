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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   TERRAIN is the first component of the REGional Climate Modeling    !
!   (RegCM) system version 3.0 and used to access archived terrain     !
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
!  This program reads terrain height data from the NCAR air force
!  terrain tapes (1 deg., 30', or 5') or terrain and landuse data
!  from the PSU/NCAR combined landuse tapes (1 deg., 30', 10'), or 10'
!  GLCC landuse data (Loveland et al 1999) and analyzes heights and/or
!  landuse values to a given grid.
!
!  terrain will be read from unit 12, landuse from unit 10 or 11,
!  output variables unit 9.
!---------------------------------------------------------------------
      use mod_dynparam
      use mod_smooth , only : smth121 , smthtr
      use mod_projections , only : lambrt , mappol , normer , rotmer
      use mod_interp , only : anal2 , interp
      use mod_fudge , only : lndfudge , texfudge , lakeadj
      use mod_rdldtr , only : rdldtr , rdldtr_nc
      use mod_maps
      use mod_block
      use mod_interfaces
      implicit none
!
! Local variables
!
      integer :: maxiter , maxjter , maxdim
      character(256) :: char_lnd , char_tex
      character(256) :: namelistfile , prgname
      character(256) :: ctlfile_s , datafile_s
      character(256) :: ctlfile , datafile
      integer :: i , j , k , minsize , ierr , i0 , j0 , m , n
      logical :: ibndry
      integer :: nunitc , nunitc_s , ctlunit , ctlunit_s
      real(4) :: clong , dsx , dsx_s , htave , htgrid_a
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

      call allocate_grid(iy,jx,kz,nveg,ntex)
      if ( nsg>1 ) call allocate_subgrid(iysg,jxsg,nveg,ntex)
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
!     iblk = dimension of arrays xobs,yobs,ht,htsd
!     estimate iblk = ihmax*jhmax where ihmax = (xmaxlat-xminlat)/xnc
!     and jhmax = (xmaxlon-xminlon)/xnc.  xnc = 1,0.5,1./6.or 5./60. for
!     1 deg, 30 min, 10 min or 5 min data, respectively.  add 1-2000 to
!     estimate of iblk for safety.  if you underestimate the
!     dimensions the program will abort with a diagnostic message
!     indicating the correct dimensions required .
!     iter,jter = dimensions of array lnd8.
!     within the search region. at present, iter and jter must be
!     equal and >= max(ihmax,jhmax) where ihmax and jhmax are
!     as calculated above.
!
      open (48,status='scratch',form='unformatted')

      dxcen = 0.0
      dycen = 0.0
!
      clong = clon
      if ( clong>180. ) clong = clong - 360.
      if ( clong<=-180. ) clong = clong + 360.
      nunitc  = 109
      ctlunit = 110
      if ( nsg>1 ) then
        nunitc_s  = 119
        ctlunit_s = 120
        write (datafile_s,99001)                                        &
              & trim(dirter), pthsep, trim(domname) , nsg , '.INFO'
        write (ctlfile_s,99001)                                         &
              & trim(dirter), pthsep, trim(domname) , nsg , '.CTL'
        call setup(nunitc_s,ctlunit_s,iysg,jxsg,ntypec_s,iproj,ds/nsg,  &
                 & clat,clong,igrads,ibyte,datafile_s,ctlfile_s)
        if ( iproj=='LAMCON' ) then
          call lambrt(xlon_s,xlat_s,xmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,0,xn,truelatl,truelath)
          call lambrt(dlon_s,dlat_s,dmap_s,coriol_s,iysg,jxsg,clong,    &
                    & clat,dsinm,1,xn,truelatl,truelath)
          write (*,*) 'XN,TRUELATL,TRUELATH = ' , xn , truelatl ,       &
                    & truelath
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
          print * , 'iproj MAP PROJECTION IS NOT AN OPTION'
          stop 999
        end if
        dsx_s = dsinm
        print * , 'after calling MAP PROJECTION, for subgrid'
!
!       reduce the search area for the domain
!       [minlat:maxlat,minlon:maxlon]
        call mxmnll(iysg,jxsg,clong,xlon_s,xlat_s,ntypec_s)
        print * , 'after calling MXMNLL, for subgrid'
!
        maxiter = (xmaxlat-xminlat)/xnc
        maxjter = (xmaxlon-xminlon)/xnc
        maxdim = max(maxiter,maxjter) + 1000
        print *, 'Allocating ' , maxdim
        call allocate_block(maxdim,maxdim)
!
!       read in the terrain & landuse data
        if ( itype_in==1 ) then
          call rdldtr(ntypec_s,nveg,ntex,lsmtyp,aertyp,ibyte)
          print * , 'after calling RDLDTR_s, for subgrid'
        else if ( itype_in==2 ) then
          call rdldtr_nc(ntypec_s,nveg,ntex,lsmtyp,aertyp)
          print * , 'after calling RDLDTR_nc, for subgrid'
        else
          print * , 'Unknown Itype for input'
          stop
        endif
        if ( ifanal ) then
!         convert xobs and yobs from LON and LAT to x and y in mesh
          call xyobsll(iysg,jxsg,iproj,clat,clong,plat,plon,        &
                     & truelath)
          print * , 'after calling XYOBSLL, for subgrid'
!
!         create the terrain height fields
          call anal2(htsdgrid_s,ht2,nobs,iysg,jxsg,corc_s,sumc_s,       &
                   & nsc_s,wtmaxc_s,htsavc_s)
          call anal2(htgrid_s,ht,nobs,iysg,jxsg,corc_s,sumc_s,nsc_s,    &
                   & wtmaxc_s,htsavc_s)
          print * , 'after calling ANAL2, for subgrid'
          do j = 1 , jxsg
            do i = 1 , iysg
              htgrid_s(i,j) = amax1(htgrid_s(i,j)*100.,0.0)
              htsdgrid_s(i,j) = amax1(htsdgrid_s(i,j)*100000.,0.0)
              htsdgrid_s(i,j) = sqrt(amax1(htsdgrid_s(i,j)-htgrid_s(i,j)&
                              & **2,0.0))
            end do
          end do
        else
          call interp(jxsg,iysg,xlat_s,xlon_s,htgrid_s,htsdgrid_s,      &
                    & ntypec_s)
          print * , 'after calling INTERP, for subgrid'
!         print*, '  Note that the terrain standard deviation is'
!         print*, '  underestimated using INTERP. (I dont know why?)'
        end if
!       create surface landuse types
        call surf(xlat_s,xlon_s,lnduse_s,iysg,jxsg,nnc,       &
                & xnc,lndout_s,land_s,nobs,h2opct, &
                & lsmtyp,sanda_s,sandb_s,claya_s,clayb_s,frac_lnd_s,    &
                & nveg,aertyp,intext_s,texout_s,frac_tex_s,ntex)
        print * , 'after calling SURF, for subgrid'
!       **** Adjust the Great Lake Heights to their actual values.
        if ( lakadj ) then
          print * ,                                                     &
               &'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
          call lakeadj(lsmtyp,lnduse_s,htgrid_s,xlat_s,xlon_s,iysg,     &
                     & jxsg)
          print * , 'after calling LAKEADJ, for subgrid'
        end if
        call smth121(htgrid_s,iysg,jxsg,hscr1_s)
        call smth121(htsdgrid_s,iysg,jxsg,hscr1_s)
!       **** Readjust the Great Lake Heights to their actual values
!       again.
        if ( lakadj ) then
          print * ,                                                     &
               &'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
          call lakeadj(lsmtyp,lnduse_s,htgrid_s,xlat_s,xlon_s,iysg,     &
                     & jxsg)
          print * , 'after calling LAKEADJ, for subgrid'
        end if
        ibndry = .true.
        if ( ibndry ) then
          do j = 2 , jxsg - 1
            htgrid_s(1,j) = htgrid_s(2,j)
            htgrid_s(iysg,j) = htgrid_s(iysg-1,j)
            lnduse_s(1,j) = lnduse_s(2,j)
            lnduse_s(iysg,j) = lnduse_s(iysg-1,j)
            lndout_s(1,j) = lndout_s(2,j)
            lndout_s(iysg,j) = lndout_s(iysg-1,j)
 
            if ( lsmtyp=='USGS' ) then
              sanda_s(1,j) = sanda_s(2,j)
              sanda_s(iysg,j) = sanda_s(iysg-1,j)
              sandb_s(1,j) = sandb_s(2,j)
              sandb_s(iysg,j) = sandb_s(iysg-1,j)
              claya_s(1,j) = claya_s(2,j)
              claya_s(iysg,j) = claya_s(iysg-1,j)
              clayb_s(1,j) = clayb_s(2,j)
              clayb_s(iysg,j) = clayb_s(iysg-1,j)
              do k = 1 , nveg
                frac_lnd_s(1,j,k) = frac_lnd_s(2,j,k)
                frac_lnd_s(iysg,j,k) = frac_lnd_s(iysg-1,j,k)
              end do
            end if
            if ( aertyp(7:7)=='1' ) then
              intext_s(1,j) = intext_s(2,j)
              intext_s(iysg,j) = intext_s(iysg-1,j)
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
            lnduse_s(i,1) = lnduse_s(i,2)
            lnduse_s(i,jxsg) = lnduse_s(i,jxsg-1)
            lndout_s(i,1) = lndout_s(i,2)
            lndout_s(i,jxsg) = lndout_s(i,jxsg-1)
 
            if ( lsmtyp=='USGS' ) then
              sanda_s(i,1) = sanda_s(i,2)
              sanda_s(i,jxsg) = sanda_s(i,jxsg-1)
              sandb_s(i,1) = sandb_s(i,2)
              sandb_s(i,jxsg) = sandb_s(i,jxsg-1)
              claya_s(i,1) = claya_s(i,2)
              claya_s(i,jxsg) = claya_s(i,jxsg-1)
              clayb_s(i,1) = clayb_s(i,2)
              clayb_s(i,jxsg) = clayb_s(i,jxsg-1)
              do k = 1 , nveg
                frac_lnd_s(i,1,k) = frac_lnd_s(i,2,k)
                frac_lnd_s(i,jxsg,k) = frac_lnd_s(i,jxsg-1,k)
              end do
            end if
            if ( aertyp(7:7)=='1' ) then
              intext_s(i,1) = intext_s(i,2)
              intext_s(i,jxsg) = intext_s(i,jxsg-1)
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
!       land/sea mask fudging
        write (char_lnd,99003) trim(dirter), pthsep, trim(domname),     &
           &   '_LANDUSE_' , nsg
        write (char_tex,99003) trim(dirter), pthsep, trim(domname),     &
           &   'TEXTURE_' , nsg
        call lndfudge(fudge_lnd_s,ch_s,lndout_s,htgrid_s,iysg,jxsg,     &
                    & lsmtyp,trim(char_lnd))
        if ( aertyp(7:7)=='1' ) call texfudge(fudge_tex_s,ch_s,texout_s,&
           & htgrid_s,iysg,jxsg,trim(char_tex))
        print * , 'after calling FUDGE, for subgrid'
!       output terrestrial fields
!       OUTPUT is used to output also the fraction of each
!       LANDUSE legend and TEXTURE type

        call free_block

      end if
!
      dxcen = 0.0
      dycen = 0.0
!
!     set up the parameters and constants
      write (datafile,99002)                                            &
           & trim(dirter), pthsep, trim(domname) , '.INFO'
      write (ctlfile,99002)                                             &
           & trim(dirter), pthsep, trim(domname) , '.CTL'
      call setup(nunitc,ctlunit,iy,jx,ntypec,iproj,ds,clat,clong,igrads,&
               & ibyte,datafile,ctlfile)
      print * , 'after calling SETUP'
!
!-----calling the map projection subroutine
      if ( iproj=='LAMCON' ) then
        call lambrt(xlon,xlat,xmap,coriol,iy,jx,clong,clat,dsinm,0,xn,  &
                  & truelatl,truelath)
        call lambrt(dlon,dlat,dmap,coriol,iy,jx,clong,clat,dsinm,1,xn,  &
                  & truelatl,truelath)
        write (*,*) 'XN,TRUELATL,TRUELATH = ' , xn , truelatl , truelath
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
        print * , 'iproj MAP PROJECTION IS NOT AN OPTION'
        stop 999
      end if
      dsx = dsinm
      print * , 'after calling MAP PROJECTION'
!
!     reduce the search area for the domain
!     [minlat:maxlat,minlon:maxlon]
      call mxmnll(iy,jx,clong,xlon,xlat,ntypec)
      print * , 'after calling MXMNLL'

      maxiter = (xmaxlat-xminlat)/xnc
      maxjter = (xmaxlon-xminlon)/xnc
      maxdim = max(maxiter,maxjter) + 1000
      print *, 'Allocating ' , maxdim
      call allocate_block(maxdim,maxdim)
!
!     read in the terrain & landuse data
      if ( itype_in==1 ) then
        call rdldtr(ntypec,nveg,ntex,lsmtyp,aertyp,ibyte)
        print * , 'after calling RDLDTR'
      else if (itype_in==2 ) then
        call rdldtr_nc(ntypec,nveg,ntex,lsmtyp,aertyp)
        print * , 'after calling RDLDTR_nc'
      else
      endif
 
!     compute the scaled standard deviation of terrain height.
!     (must be called before XYOBSLL because xobs/yobs are modified)
!     CALL SCALESD
!     print*, 'after calling SCALESD'
 
      if ( ifanal ) then
!       convert xobs and yobs from LON and LAT to x and y in mesh
        call xyobsll(iy,jx,iproj,clat,clong,plat,plon,truelath)
        print * , 'after calling XYOBSLL'
 
!       create the terrain height fields
        call anal2(htsdgrid,ht2,nobs,iy,jx,corc,sumc,nsc,wtmaxc,htsavc)
        call anal2(htgrid,ht,nobs,iy,jx,corc,sumc,nsc,wtmaxc,htsavc)
        print * , 'after calling ANAL2'
        do j = 1 , jx
          do i = 1 , iy
            htgrid(i,j) = amax1(htgrid(i,j)*100.,0.0)
            htsdgrid(i,j) = amax1(htsdgrid(i,j)*100000.,0.0)
            htsdgrid(i,j) = sqrt(amax1(htsdgrid(i,j)-htgrid(i,j)**2,0.0)&
                          & )
          end do
        end do
      else
        call interp(jx,iy,xlat,xlon,htgrid,htsdgrid,ntypec)
!       print*, 'after calling INTERP'
!       print*, '  Note that the terrain standard deviation is'
!       print*, '  underestimated using INTERP. (I dont know why!)'
      end if
 
 
!     create surface landuse types
      call surf(xlat,xlon,lnduse,iy,jx,nnc,xnc,lndout,land,   &
              & nobs,h2opct,lsmtyp,sanda,sandb,    &
              & claya,clayb,frac_lnd,nveg,aertyp,intext,texout,frac_tex,&
              & ntex)
      print * , 'after calling SURF'
 
!     **** Adjust the Great Lake Heights to their actual values.
      if ( lakadj ) then
        print * , 'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
        call lakeadj(lsmtyp,lnduse,htgrid,xlat,xlon,iy,jx)
        print * , 'after calling LAKEADJ'
      end if
 
!     ******           preliminary heavy smoothing of boundaries
      if ( smthbdy ) call smthtr(htgrid,iy,jx)
 
!     ******           grell smoothing to eliminate 2 delx wave (6/90):
      call smth121(htgrid,iy,jx,hscr1)
      call smth121(htsdgrid,iy,jx,hscr1)
 
!     **** Readjust the Great Lake Heights to their actual values again.
      if ( lakadj ) then
        print * , 'CALLING LAKEADJ FOR THE SECOND TIME (after 2dx pass)'
        call lakeadj(lsmtyp,lnduse,htgrid,xlat,xlon,iy,jx)
      end if
 
      ibndry = .true.
      if ( ibndry ) then
        do j = 2 , jx - 1
          htgrid(1,j) = htgrid(2,j)
          htgrid(iy,j) = htgrid(iy-1,j)
          lnduse(1,j) = lnduse(2,j)
          lnduse(iy,j) = lnduse(iy-1,j)
          lndout(1,j) = lndout(2,j)
          lndout(iy,j) = lndout(iy-1,j)
 
          if ( lsmtyp=='USGS' ) then
            sanda(1,j) = sanda(2,j)
            sanda(iy,j) = sanda(iy-1,j)
            sandb(1,j) = sandb(2,j)
            sandb(iy,j) = sandb(iy-1,j)
            claya(1,j) = claya(2,j)
            claya(iy,j) = claya(iy-1,j)
            clayb(1,j) = clayb(2,j)
            clayb(iy,j) = clayb(iy-1,j)
            do k = 1 , nveg
              frac_lnd(1,j,k) = frac_lnd(2,j,k)
              frac_lnd(iy,j,k) = frac_lnd(iy-1,j,k)
            end do
          end if
          if ( aertyp(7:7)=='1' ) then
            intext(1,j) = intext(2,j)
            intext(iy,j) = intext(iy-1,j)
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
          lnduse(i,1) = lnduse(i,2)
          lnduse(i,jx) = lnduse(i,jx-1)
          lndout(i,1) = lndout(i,2)
          lndout(i,jx) = lndout(i,jx-1)
 
          if ( lsmtyp=='USGS' ) then
            sanda(i,1) = sanda(i,2)
            sanda(i,jx) = sanda(i,jx-1)
            sandb(i,1) = sandb(i,2)
            sandb(i,jx) = sandb(i,jx-1)
            claya(i,1) = claya(i,2)
            claya(i,jx) = claya(i,jx-1)
            clayb(i,1) = clayb(i,2)
            clayb(i,jx) = clayb(i,jx-1)
            do k = 1 , nveg
              frac_lnd(i,1,k) = frac_lnd(i,2,k)
              frac_lnd(i,jx,k) = frac_lnd(i,jx-1,k)
            end do
          end if
          if ( aertyp(7:7)=='1' ) then
            intext(i,1) = intext(i,2)
            intext(i,jx) = intext(i,jx-1)
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
!     land/sea mask fudging
      write (char_lnd,99004) trim(dirter), pthsep, trim(domname),       &
           &   '_LANDUSE'
      write (char_tex,99004) trim(dirter), pthsep, trim(domname),       &
           &   '_TEXTURE'
      call lndfudge(fudge_lnd,ch,lndout,htgrid,iy,jx,lsmtyp,            &
         & trim(char_lnd))
      if ( aertyp(7:7)=='1' ) call texfudge(fudge_tex,ch,texout,htgrid, &
         & iy,jx,trim(char_tex))
      print * , 'after calling FUDGE'
!     output terrestrial fields
!     OUTPUT is used to output also the fraction of each
!     LANDUSE legend and TEXTURE type
      call free_block

      if ( nsg>1 ) then
        do i = 1 , iy
          do j = 1 , jx
            i0 = (i-1)*nsg
            j0 = (j-1)*nsg
            htave = 0.0
            do m = 1 , nsg
              do n = 1 , nsg
                if ( lsmtyp=='BATS' ) then
                  if ( htgrid(i,j)<0.1 .and.                            &
                     & (lndout(i,j)>14.5 .and. lndout(i,j)<15.5) ) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 15.
                  end if
                else if ( lsmtyp=='USGS' ) then
                  if ( htgrid(i,j)<0.1 .and. lndout(i,j)>24.5 ) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 25.
                  end if
                else
                end if
                htave = htave + htgrid_s(i0+m,j0+n)
              end do
            end do
            htgrid_a = htave/float(nsg*nsg)
            do m = 1 , nsg
              do n = 1 , nsg
                htgrid_s(i0+m,j0+n) = htgrid_s(i0+m,j0+n) - htgrid_a +  &
                                    & htgrid(i,j)
              end do
            end do
          end do
        end do
        call output(nunitc_s,ctlunit_s,iysg,jxsg,dsx_s,clat,clong,      &
                  & plat,plon,iproj,htgrid_s,htsdgrid_s,lndout_s,       &
                  & xlat_s,xlon_s,dlat_s,dlon_s,xmap_s,dattyp,dmap_s,   &
                  & coriol_s,snowam_s,igrads,ibigend,kz,sigma,mask_s,   &
                  & ptop,nsg,truelatl,truelath,xn,domname,lsmtyp,       &
                  & sanda_s,sandb_s,claya_s,clayb_s,frac_lnd_s,nveg,    &
                  & aertyp,texout_s,frac_tex_s,ntex,.false.)
        print * , 'after calling OUTPUT, for subgrid'
        call free_subgrid
      end if

      call output(nunitc,ctlunit,iy,jx,dsx,clat,clong,plat,plon,iproj,  &
                & htgrid,htsdgrid,lndout,xlat,xlon,dlat,dlon,xmap,      &
                & dattyp,dmap,coriol,snowam,igrads,ibigend,kz,sigma,    &
                & mask,ptop,nsg,truelatl,truelath,xn,domname,           &
                & lsmtyp,sanda,sandb,claya,clayb,frac_lnd,nveg,aertyp,  &
                & texout,frac_tex,ntex,.true.)
      print * , 'after calling OUTPUT'
      call free_grid

      close (48, status='delete')

      print *, 'Successfully completed terrain fields generation'
 
!     stop 9999
!
99001 format (a,a,a,i0.3,a)
99002 format (a,a,a,a)
99003 format (a,a,a,a8,i0.3)
99004 format (a,a,a,a8)
      end program terrain
