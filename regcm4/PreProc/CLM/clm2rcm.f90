      program clmproc
 
      use netcdf
      use mod_param
      use mod_param_clm

      implicit none
!
      real(4) , parameter :: vmisdat=-99.e30
      integer , parameter :: ndim = 3
      logical , parameter :: bvoc = .false.
!
! Local variables
!
      real :: clatx , clonx , dsx , grdfac , grdfacx , offset , perr ,  &
           & platx , plonx , pmax , ptopx , xscale , xlatmax , xlatmin ,&
           & xlonmax , xlonmin
      real , allocatable , dimension(:) :: glat , glon , zlat ,         &
               &                           zlev , zlon
      integer :: i , ibigendx , idatex , idin , idout , idy , ierr ,    &
               & ifield , ifld , igradsx , ihr , ilatmax , ilatmin ,    &
               & ilonmax , ilonmin , imap , imo , iotyp , irec , iyr ,  &
               & ixx , j , julnc , jxx , kk , kmax , kzz , l , nf_close
      integer :: iunout , iunctl
      integer , dimension(4) :: icount , istart
      integer , dimension(3) :: iadim
      character(6) :: iprojx
      real , allocatable , dimension(:,:) :: landmask , regxy , sandclay
      character(30) , dimension(nfld) :: lnam
      character(53) :: outfil_ctl , outfil_gr , outfil_nc
      real , allocatable , dimension(:,:,:) :: regxyz
      real , allocatable , dimension(:,:,:,:) :: regyxzt , zoom
      real , dimension(1) :: sig1d
      real , dimension(2) :: sigf
      real , dimension(kz+1) :: sigx
      logical :: there
      character(13) , dimension(nfld) :: units
      real , dimension(3) :: varmax , varmin
      character(15) , dimension(nfld) :: vnam_o
      real(8) :: xhr
      real , dimension(ix,jx) :: xlat , xlon
      real , dimension(ix) :: xlat1d
      real , dimension(jx) :: xlon1d
!
!     ** Get latitudes and longitudes from DOMAIN.INFO
      open (unit=10,file='fort.10',status='old',form='unformatted',     &
          & recl=ix*jx*ibyte,access='direct')
      read (10,rec=1) ixx , jxx , kzz , dsx , clatx , clonx , platx ,   &
                    & plonx , grdfacx , iprojx , sigx , ptopx ,         &
                    & igradsx , ibigendx
      if ( ixx/=ix .or. jxx/=jx ) then
        print * , 'DOMAIN.INFO is inconsistent with domain.param'
        print * , '  domain.param:    ix=' , ix , '   jx=' , jx
        print * , '  DOMAIN.INFO:     ix=' , ixx , '   jx=' , jxx
        stop 780
      end if
      if ( abs(ds*1000.-dsx)>0.01 ) then
        print * , 'DOMAIN.INFO is inconsistent with domain.param'
        print * , '  domain.param: ds=' , ds*1000.
        print * , '  DOMAIN.INFO:  ds=' , dsx
        stop 781
      end if
      if ( clatx/=clat .or. clonx/=clon .or. platx/=plat .or.           &
         & plonx/=plon ) then
        print * , 'DOMAIN.INFO is inconsistent with domain.param'
        print * , '  domain.param:  clat=' , clat , ' clon=' , clon
        print * , '  DOMAIN.INFO:   clat=' , clatx , ' clon=' , clonx
        print * , '  domain.param:  plat=' , plat , ' plon=' , plon
        print * , '  DOMAIN.INFO:   plat=' , platx , ' plon=' , plonx
        stop 782
      end if
      if ( iprojx/=iproj ) then
        print * , 'DOMAIN.INFO is inconsistent with domain.param'
        print * , '  domain.param: iproj=' , iproj
        print * , '  DOMAIN.INFO:  iproj=' , iprojx
        stop 783
      end if
      read (10,rec=5) ((xlat(i,j),j=1,jx),i=1,ix)
      read (10,rec=6) ((xlon(i,j),j=1,jx),i=1,ix)
      close (10)
 
!     ** Set output variables
      iotyp = 2
      xscale = 1.
      offset = 0.
 
!     ** Open direct access CLM3/RegCM3 output file
      outfil_gr = trim(outdir)//trim(outfil)
      call fexist(outfil_gr)
      open (iunout,file=outfil_gr,status='unknown',form='unformatted',  &
          & recl=jx*ix*ibyte,access='direct')
      irec = 1
      write (iunout,rec=irec) ixx , jxx , npft , nsoi , dsx , clatx ,   &
                            & clonx , platx , plonx , grdfac , iprojx
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
        inquire (file=infil(ifld),exist=there)
        if ( .not.there ) then
          print * , 'CLM Input file does not exist: ' , infil(ifld)
          stop 'NON-EXISTANT FILE'
        end if
        if ( ifld==ipft .or. ifld==ilai .or. ifld==ilak .or.            &
           & ifld==iglc .or. ifld==iurb .or. ifld==isnd .or.            &
           & ifld==icol .or. ifld==ioro .or. ifld==iiso .or.            &
           & ifld==ifma .or. ifld==imbo .or. ifld==ibpin .or.           &
           & ifld==iapin ) then
!         ************************ CHANGED LINE ABOVE to include iiso
!         ************************
          if ( ifld==iiso .or. ifld==ibpin .or. ifld==iapin .or.        &
             & ifld==imbo ) then   !****** added if statement to read isoprene differently ABT  ***
            ierr = nf90_open(infil(ifld),nf90_nowrite,idin)
          else
            ierr = nf90_open(infil(ifld),nf90_nowrite,idin)
          end if
!******************************************* ABT ***************************************
          outfil_nc = trim(outdir)//'RCM'//trim(infil(ifld)(7:50))
!         CALL FEXIST(outfil_nc)
          print * , 'OPENING NetCDF FILE: ' , outfil_nc
          call rcrecdf(outfil_nc,idout,varmin,varmax,3,ierr)
        end if
 
!       ** Setup RegCM3 grid variables
        call param(jx,ix,nlev(ifld),ds,clatx,clonx,platx,plonx,xlat,    &
                 & xlon,varmin,varmax,xlat1d,xlon1d,xlonmin,xlonmax,    &
                 & xlatmin,xlatmax,iadim,ndim)
 
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
        call clm3grid2(nlon(ifld),nlat(ifld),ifield,glon,glat,istart,   &
                     & icount,ifld,zlon,zlat,zlev)
 
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
          call readcdfr4(idin,vnam_lm,lnam(ifld),units(ifld),istart(1), &
                       & icount(1),istart(2),icount(2),1,1,1,1,landmask)
          call readcdfr4(idin,vnam_st,lnam(ifld),units(ifld),istart(1), &
                       & icount(1),istart(2),icount(2),1,1,1,1,zoom)
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
 
!****************** ADDED IF STATEMENT TO INCLUDE BVOCS ABT **************
        else if ( ifld==iiso .or. ifld==ibpin .or. ifld==iapin .or.     &
                & ifld==imbo ) then
          call readcdfr4_iso(idin,vnam(ifld),lnam(ifld),units(ifld),    &
                           & istart(1),icount(1),istart(2),icount(2),   &
                           & istart(3),icount(3),istart(4),icount(4),   &
                           & zoom)
!***************** ADDED 'IF' TO INCLUDE BVOCS ABT ************************
 
        else
          if ( ifld/=icol ) call readcdfr4(idin,vnam_lm,lnam(ifld),     &
             & units(ifld),istart(1),icount(1),istart(2),icount(2),1,1, &
             & 1,1,landmask)      !HI-RES SHOULD HAVE LANDMASK FOR ALL EXCEPT SOILCOLOR (abt edit)
          call readcdfr4(idin,vnam(ifld),lnam(ifld),units(ifld),        &
                       & istart(1),icount(1),istart(2),icount(2),       &
                       & istart(3),icount(3),istart(4),icount(4),zoom)
        end if
 
        if ( ifld==icol .or. ifld==iiso .or. ifld==ibpin .or.           &
           & ifld==imbo .or. ifld==iapin ) then         ! ******** ABT added iiso for calculating landmask *******
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
        allocate(regyxzt(ix,jx,nlev(ifld),ntim(ifld)))
        allocate(regxyz(jx,ix,nlev(ifld)))
        allocate(regxy(jx,ix))
        call bilinx4d(zoom,zlon,zlat,icount(1),icount(2),regyxzt,xlon,  &
                    & xlat,ix,jx,icount(3),icount(4),vmin(ifld),vmisdat)
 
!       ** Write the interpolated data to NetCDF and direct-access
!       output
        do l = 1 , ntim(ifld)
          idatex = 1900000000 + l*10000 + 1500
          call julian(idatex,julnc,iyr,imo,idy,ihr,xhr)
          if ( ifld==ipft ) then
            do j = 1 , ix
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
            do j = 1 , ix
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
 
          call writecdf(idout,vnam(ifld),regxyz,jx,ix,nlev(ifld),iadim, &
                      & xhr,lnam(ifld),units(ifld),xscale,offset,varmin,&
                      & varmax,xlat1d,xlon1d,zlev,0,vmisdat,iotyp)
        end do
 
!       ** Deallocate variables for next CLM3 field
        deallocate(glon)
        deallocate(glat)
        deallocate(zoom)
        deallocate(zlon)
        deallocate(zlat)
        deallocate(zlev)
        deallocate(regyxzt)
        deallocate(regxyz)
        if ( ifld==isnd .or. ifld==icly ) deallocate(sandclay)
 
      end do  ! End nfld loop
 
!     ** Make GrADS CTL file
      outfil_ctl = trim(outdir)//trim(outfil)//'.CTL'
      open (iunctl,file=outfil_ctl,status='unknown')
!     do ifld=1,nfld
      do ifld = 1 , ifield
        if ( ifld==ilai .or. ifld==isai .or. ifld==itop .or.            &
           & ifld==ibot ) then
          vnam_o(ifld) = trim(vnam(ifld)(9:20))
        else
          vnam_o(ifld) = trim(vnam(ifld)(1:15))
        end if
      end do
      call makectl(iunctl,outfil,jx,ix,npft,ifield,dsx,clatx,clonx,     &
                 & platx,plonx,iprojx,ibigendx,truelatl,truelath,       &
                 & grdfacx,vmisdat,vnam_o,lnam,nlev,ntim,varmin,varmax)
!     ** Close files
      call clscdf(idin,ierr)
      call clscdf(idout,ierr)
      close (iunout)
      close (iunctl)
 
      end program clmproc
