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
      use mod_dynparam
      use mod_param_clm

      implicit none
!
      real(4) , parameter :: vmisdat=-9999
      integer , parameter :: ndim = 3
      logical , parameter :: bvoc = .false.
!
! Local variables
!
      real(4) :: clatx , clonx , dsx , grdfacx , offset ,               &
           & perr , platx , plonx , pmax , ptopx , xscale , xlatmax ,   &
           & xlatmin , xlonmax , xlonmin
      integer :: i , ibigendx , idatex , idin , idout , idy , ierr ,    &
               & ifield , ifld , igradsx , ihr , imap , imo , jotyp ,   &
               & irec , iyr , iyy , j , julnc , jxx , kmax , kzz , l
      integer :: k
      integer :: iunout , iunctl
      integer , dimension(4) :: icount , istart
      integer , dimension(3) :: iadim
      character(6) :: iprojx
      character(64) , dimension(nfld) :: lnam
      character(256) :: outfil_ctl , outfil_gr , outfil_nc
      logical :: there
      character(64) , dimension(nfld) :: units
      real(4) , dimension(3) :: varmax , varmin
      character(64) , dimension(nfld) :: vnam_o
      real(8) :: xhr
      real(4) , allocatable , dimension(:) :: glat , glon , zlat ,      &
               &                           zlev , zlon
      real(4) , allocatable , dimension(:,:,:) :: regxyz
      real(4) , allocatable , dimension(:,:,:,:) :: regyxzt , zoom
      real(4) , allocatable , dimension(:,:) :: landmask , regxy ,      &
               &                             sandclay
      real(4) , allocatable , dimension(:,:) :: xlat , xlon
      real(4) , allocatable , dimension(:) :: xlat1d
      real(4) , allocatable , dimension(:) :: xlon1d
      real(4) , allocatable , dimension(:) :: sigx
      integer :: ipathdiv
      character(256) :: namelistfile, prgname
      character(256) :: terfile , inpfile
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

      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))
      allocate(xlat1d(iy))
      allocate(xlon1d(jx))
      allocate(sigx(kz+1))
!
!     ** Get latitudes and longitudes from DOMAIN.INFO
!
      write (terfile,99001)                                             &
        & trim(dirter), pthsep, trim(domname), '.INFO'
      open (unit=10,file=terfile,status='old',form='unformatted',       &
          & recl=iy*jx*ibyte,access='direct')
      read (10,rec=1) iyy , jxx , kzz , dsx , clatx , clonx , platx ,   &
                    & plonx , grdfacx , iprojx , sigx , ptopx ,         &
                    & igradsx , ibigendx
      if ( iyy/=iy .or. jxx/=jx ) then
        print * , 'DOMAIN.INFO is inconsistent with regcm.in'
        print * , '  namelist       :    iy=' , iy , '   jx=' , jx
        print * , '  DOMAIN.INFO    :    iy=' , iyy , '   jx=' , jxx
        stop 780
      end if
      if ( abs(ds*1000.-dsx)>0.01 ) then
        print * , 'DOMAIN.INFO is inconsistent with regcm.in'
        print * , '  namelist       : ds=' , ds*1000.
        print * , '  DOMAIN.INFO    : ds=' , dsx
        stop 781
      end if
      if ( clatx/=clat .or. clonx/=clon .or. platx/=plat .or.           &
         & plonx/=plon ) then
        print * , 'DOMAIN.INFO is inconsistent with regcm.in'
        print * , '  namelist       :  clat=' , clat , ' clon=' , clon
        print * , '  DOMAIN.INFO    :  clat=' , clatx , ' clon=' , clonx
        print * , '  namelist       :  plat=' , plat , ' plon=' , plon
        print * , '  DOMAIN.INFO    :  plat=' , platx , ' plon=' , plonx
        stop 782
      end if
      if ( iprojx/=iproj ) then
        print * , 'DOMAIN.INFO is inconsistent with regcm.in'
        print * , '  namelist       : iproj=' , iproj
        print * , '  DOMAIN.INFO    : iproj=' , iprojx
        stop 783
      end if
      read (10,rec=5) ((xlat(i,j),j=1,jx),i=1,iy)
      read (10,rec=6) ((xlon(i,j),j=1,jx),i=1,iy)
      close (10)
 
!     ** Set output variables
      jotyp = 2
      xscale = 1.
      offset = 0.
      iunout = 101
      iunctl = 102
 
!     ** Open direct access CLM3/RegCM3 output file
      outfil_gr = trim(dirglob)//pthsep//trim(domname)//'_CLM3.INFO'
      call fexist(outfil_gr)
      print *, 'Open ', trim(outfil_gr)
      open (iunout,file=outfil_gr,status='replace',                     &
          & form='unformatted',recl=jx*iy*ibyte,access='direct')
      irec = 1
      write (iunout,rec=irec) iyy , jxx , npft , nsoi , dsx , clatx ,   &
                            & clonx , platx , plonx , grdfacx , iprojx
!abt  added below
!     ** determine which files to create (emission factor map or not)
      call comp(ifield,bvoc)
!abt  added above
 
!     ** Loop over the fields defined in clm.param
!     do ifld=1,nfld
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
          call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),1,      &
                       & ntim(ifld),1,nlev(ifld),1,1,1,1,sandclay)
          print *, 'Read ', trim(lnam(ifld))
          call readcdfr4(idin,vnam_lm,lnam(ifld),units(ifld),istart(1), &
                       & icount(1),istart(2),icount(2),1,1,1,1,landmask)
          print *, 'Read ', trim(lnam(ifld))
          call readcdfr4(idin,vnam_st,lnam(ifld),units(ifld),istart(1), &
                       & icount(1),istart(2),icount(2),1,1,1,1,zoom)
          print *, 'Read ', trim(lnam(ifld))
          do j = 1 , icount(2)
            do i = 1 , icount(1)
              imap = nint(zoom(i,j,1,1))
              do k = 1 , icount(3)
                if ( imap>0 .and. landmask(i,j)>0.5 ) then
                  zoom(i,j,k,1) = sandclay(imap,k)
                else
                  zoom(i,j,k,1) = vmisdat
                end if
              end do
            end do
          end do
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
        allocate(regxy(jx,iy))

        call bilinx4d(zoom,zlon,zlat,icount(1),icount(2),regyxzt,xlon,  &
                    & xlat,iy,jx,icount(3),icount(4),vmin(ifld),vmisdat)
 
!       ** Write the interpolated data to NetCDF and direct-access
!       output
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
                  regxy(i,j) = regyxzt(j,i,k,l)
                else
                  regxyz(i,j,k) = 0
                  regxy(i,j) = 0
                end if
              end do
            end do
            irec = irec + 1
            write (iunout,rec=irec) regxy
          end do
 
          call writecdf(idout,vnam(ifld),regxyz,jx,iy,nlev(ifld),iadim, &
                      & xhr,lnam(ifld),units(ifld),xscale,offset,varmin,&
                      & varmax,xlat1d,xlon1d,zlev,0,vmisdat,jotyp)
        end do
 
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
        deallocate(regxy)
        if ( ifld==isnd .or. ifld==icly ) deallocate(sandclay)
 
      end do  ! End nfld loop

      deallocate(xlat)
      deallocate(xlon)
      deallocate(xlat1d)
      deallocate(xlon1d)
      deallocate(sigx)
 
!     ** Make GrADS CTL file
      outfil_ctl = trim(dirglob)//pthsep//trim(domname)//'_CLM3.CTL'
      open (iunctl,file=outfil_ctl,status='unknown')
!     do ifld=1,nfld
      do ifld = 1 , ifield
        vnam_o(ifld) = vnam(ifld)
      end do
      call makectl(iunctl,domname,jx,iy,npft,ifield,dsx,clatx,clonx,    &
                 & platx,plonx,iprojx,ibigendx,truelatl,truelath,       &
                 & vmisdat,vnam_o,lnam,nlev,ntim,varmin,varmax)
!     ** Close files
      call clscdf(idin,ierr)
      call clscdf(idout,ierr)
      close (iunout)
      close (iunctl)
      print *, 'Successfully completed CLM preprocessing.'
 
99001 format (a,a,a,a)
      end program clmproc
