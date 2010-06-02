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

      module mod_sst_1deg

      contains

      subroutine sst_1deg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! GISST2.3b    UKMO SST (Rayner et al 1996), 1 degree                c
!              from UKMO DATA archive (http://www.badc.rl.ac.uk/)    c
!              and reformed as direct-accessed binary GrADS format   c
!              in file GISST_187101_200209                           c
!              ML= 1 is-179.5; ML= 2 is-178.5; => ML=360 is 179.5E   c
!              NL= 1 is -89.5; NL= 2 is -88.5; => NL=180 is  89.5    c
!              see the GrADS control file for details.               c
!                                                                    c
! OISST        from CAC Optimal Interpolation dataset.               c
!              in the original netCDF format.                        c
!              ML= 1 is   0.5; ML= 2 is   1.5; => ML=360 is 359.5E   c
!              NL= 1 is -89.5; NL= 2 is -88.5; => NL=180 is  89.5    c
!                                                                    c
! OI2ST        both SST and SeaIce in the original netCDF format.    c
!                                                                    c
! OI_WK        weekly OISST in the original netCDF format.           c
!                                                                    c
! OI2WK        weekly OISST and SeaIce in the original netCDF format.c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_date
      use mod_interp , only : bilinx
      use mod_printl

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
      integer , parameter :: idtbc = 6
!
! Local variables
!
      real(4) , dimension(ilon,jlat) :: sst , ice
      integer :: idate , idate0 , nsteps , istep
      integer :: i , idatef , idateo , j , k , ludom , lumax , mrec ,   &
               & nday , nmo , nrec , nyear , nh
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      integer , dimension(25) :: lund
      character(256) :: terfile , sstfile , inpfile
      logical :: there
      real(4) , allocatable , dimension(:,:) :: lu , sstmm , icemm ,    &
                                  &             xlat , xlon
!
      allocate(lu(iy,jx))
      allocate(sstmm(iy,jx))
      allocate(icemm(iy,jx))
      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))
!
      if ( ssttyp=='GISST' ) then
        if ( globidate1<1947121512 .or. globidate2>2002091512 ) then
          print * , 'GISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        open (11,file=trim(inpglob)//'/SST/GISST_194712_200209',        &
             &form='unformatted',recl=ilon*jlat*ibyte,access='direct',  &
            & status='old',err=100)
      else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.               &
            &   ssttyp=='OI2ST' ) then
        if ( globidate1<1981121512 .or. globidate2<1981121512 ) then
          print * , 'OISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        inquire (file=trim(inpglob)//'/SST/sst.mnmean.nc',exist=there)
        if ( .not.there ) print * , 'sst.mnmean.nc is not available' ,  &
                               &' under ',trim(inpglob),'/SST/'
        if ( ssttyp=='OI2ST' ) then
          inquire (file=trim(inpglob)//'/SST/icec.mnmean.nc',           &
              &    exist=there)
          if ( .not. there )                                            &
            & print *, 'icec.mnmean.nc is not available' ,              &
                   &' under ',trim(inpglob),'/SST/'
        end if
      else if ( ssttyp=='OI_WK' .or. ssttyp=='OI2WK' ) then
        if ( globidate1<1981110100 .or. globidate2<1981110106 ) then
          print * , 'OI_WK (or OI2WK) data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        inquire (file=trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc',    &
             &   exist=there)
        if ( .not.there ) print * ,                                     &
                             &'sst.wkmean.1981-1989.nc is not available'&
                            & , ' under ',trim(inpglob),'/SST/'
        inquire (file=trim(inpglob)//'/SST/sst.wkmean.1990-present.nc', &
               & exist=there)
        if ( .not.there ) print * ,                                     &
                          &'sst.wkmean.1990-present.nc is not available'&
                         & , ' under ',trim(inpglob),'/SST/'
        if ( ssttyp=='OI2WK' ) then
          inquire (file=trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc', &
                & exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1981-1989.nc is not available',       &
                &    ' under ',trim(inpglob),'/SST/'
          end if
          inquire (file=trim(inpglob)//                                 &
                 '/SST/icec.wkmean.1990-present.nc',exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1990-present.nc is not available',    &
                &    ' under ',trim(inpglob),'/SST/'
          end if
        end if
      else
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported are GISST OISST OI_NC OI2ST OI_WK OI2WK'
        stop
      end if
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname),       &
             & '_SST.RCM'
      open (21,file=sstfile,form='unformatted',status='replace')
 
!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      write (terfile,99001)                                             &
        & trim(dirter), pthsep, trim(domname), '.INFO'
      open (10,file=terfile,form='unformatted',recl=iy*jx*ibyte,        &
          & access='direct',status='unknown',err=200)
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname),       &
          & '_RCM_SST.dat'
      open (25,file=sstfile,status='unknown',form='unformatted',        &
          & recl=iy*jx*ibyte,access='direct')
      if ( igrads==1 ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname),     &
           &  '_RCM_SST.ctl'
        open (31,file=sstfile,status='replace')
        write (31,'(a,a,a)') 'dset ^',trim(domname),'_RCM_SST.dat'
      end if
      if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname),     &
           &  '_RCM_ICE.dat'
        open (26,file=sstfile,status='unknown',form='unformatted',      &
          & recl=iy*jx*ibyte,access='direct')
        if ( igrads==1 ) then
          write (sstfile,99001) trim(dirglob), pthsep, trim(domname),   &
            & '_RCM_ICE.ctl'
          open ( 32,file=sstfile,status='replace')
          write (32,'(a,a,a)') 'dset ^',trim(domname),'_RCM_ICE.dat'
        end if
      end if
!#####
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!#####
        idateo = iprevmon(globidate1)
        idatef = inextmon(globidate2)
        print * , globidate1 , globidate2 , idateo , idatef
        call gridml1d(xlon,xlat,lu,iy,jx,idateo,idatef,ibyte,ssttyp)
!#####
      else
!#####
        idateo = ifodweek(globidate1)
        idatef = iladweek(globidate2)
        nsteps = iwkdiff(idatef , idateo)
        print * , globidate1 , globidate2 , idateo , idatef , nsteps
        call gridml2(xlon,xlat,lu,iy,jx,idateo,nsteps,ibyte,ssttyp)
!#####
      end if
!#####
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = .5 + float(i-1)
      end do
      do j = 1 , jlat
        lati(j) = -89.5 + 1.*float(j-1)
      end do
 
!#####
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!#####
!       ****** OISST SST DATA, 1 Deg data, AVAILABLE FROM 12/1981 TO
!       PRESENT ****** GISST SST DATA, 1 Deg data, AVAILABLE FROM
!       12/1947 TO 9/2002
        idate = idateo
        idate0 = idateo
        do while ( idate<=idatef )
          call split_idate(idate, nyear, nmo, nday, nh)
          if ( ssttyp=='GISST' ) then
            nrec = (nyear-1947)*12 + nmo - 11
            read (11,rec=nrec) sst
          else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.           &
                    ssttyp=='OI2ST') then
            write (*,*) idate, idate0
            inpfile = trim(inpglob)//'/SST/sst.mnmean.nc'
            call sst_mn(idate,idate0,ilon,jlat,sst,inpfile)
            if ( ssttyp=='OI2ST' ) then
              inpfile = trim(inpglob)//'/SST/icec.mnmean.nc'
              call ice_mn(idate,idate0,ilon,jlat,ice,inpfile)
            end if
          else
          end if
 
!         ******           PRINT OUT DATA AS A CHECK
          if ( nmo==1 ) call printl(sst,ilon,jlat)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
 
          if ( ssttyp=='OI2ST' ) call bilinx(ice,icemm,xlon,xlat,       &
              & loni,lati,ilon,jlat,iy,jx,1)
          do j = 1 , jx
            do i = 1 , iy
              if ( lsmtyp=='BATS' .and. sstmm(i,j)<-5000. .and.         &
                 & (lu(i,j)>13.5 .and. lu(i,j)<15.5) ) then
                do k = 1 , 20
                  lund(k) = 0.0
                end do
                lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
                lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
                lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
                lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
                lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
                lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
                lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
                lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
                ludom = 18
                lumax = 0
                do k = 1 , 20
                  if ( k<=13 .or. k>=16 ) then
                    if ( lund(k)>lumax ) then
                      ludom = k
                      lumax = lund(k)
                    end if
                  end if
                end do
                lu(i,j) = float(ludom)
                print * , ludom , sstmm(i,j)
              else if ( lsmtyp=='USGS' .and. sstmm(i,j)<-5000. .and.    &
                      & lu(i,j)>24.5 ) then
                do k = 1 , 25
                  lund(k) = 0.0
                end do
                lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
                lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
                lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
                lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
                lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
                lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
                lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
                lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
                ludom = 13
                lumax = 0
                do k = 1 , 24
                  if ( lund(k)>lumax ) then
                    ludom = k
                    lumax = lund(k)
                  end if
                end do
                lu(i,j) = float(ludom)
                print * , ludom , sstmm(i,j)
              else
              end if
              if ( sstmm(i,j)>-100. ) then
                sstmm(i,j) = sstmm(i,j) + 273.15
              else
                sstmm(i,j) = -9999.
              end if
            end do
          end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
          nday = 1
          if ( ssttyp/='OI2ST' ) then
            write (21) nday , nmo , nyear , sstmm
          else
            write (21) nday , nmo , nyear , sstmm , icemm
          end if
          print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
          if ( ssttyp=='OI2ST' )                                        &
           &    write (26,rec=mrec) ((icemm(i,j),j=1,jx),i=1,iy)
          idate = inextmon(idate)
        end do
        write (10,rec=4) ((lu(i,j),j=1,jx),i=1,iy)
!#####
      else
!#####
        idate = idateo
        idate0 = idateo
        if ( idate < 1989123100 ) then
          istep = iwkdiff(idate, 1981102900) - 1
        else
          istep = iwkdiff(idate, 1989123100) - 1
        end if
        do k = 1 , nsteps
          call split_idate(idate, nyear, nmo, nday, nh)
          write (*,*) idate, idate0 , istep+k
          if ( idate<1989123100 ) then
            inpfile = trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc'
          else
            inpfile = trim(inpglob)//'/SST/sst.wkmean.1990-present.nc'
          end if
          call sst_wk(idate,idate0,istep+k,ilon,jlat,sst,inpfile)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
          if ( ssttyp=='OI2WK') then
            if ( idate<1989123100 ) then
              inpfile = trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc'
            else
             inpfile = trim(inpglob)//'/SST/icec.wkmean.1990-present.nc'
            end if
            call ice_wk(idate,idate0,istep+k,ilon,jlat,sst,inpfile)
            call bilinx(ice,icemm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          end if 

          do j = 1 , jx
            do i = 1 , iy
              if ( sstmm(i,j)>-100. ) then
                sstmm(i,j) = sstmm(i,j) + 273.15
              else
                sstmm(i,j) = -9999.
              end if
            end do
          end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
          write (21) nday , nmo , nyear , sstmm
          print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear , idate
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
          idate = inextwk(idate)
        end do
 
!#####
      end if
!#####
      deallocate(lu)
      deallocate(sstmm)
      deallocate(icemm)
      deallocate(xlat)
      deallocate(xlon)
 
      return

 100  continue
      print * , 'ERROR OPENING GISST FILE'
      stop '4810 IN PROGRAM RDSST'
 200  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'

99001 format (a,a,a,a)
      end subroutine sst_1deg
!
!-----------------------------------------------------------------------
!
      subroutine gridml1d(xlon,xlat,lu,iy,jx,idate1,idate2,ibyte,ssttyp)
      implicit none
!
! Dummy arguments
!
      integer :: ibyte , idate1 , idate2 , iy , jx
      real(4) , dimension(iy,jx) :: lu , xlat , xlon
      character(5) :: ssttyp
      intent (in) ibyte , idate1 , idate2 , iy , jx , ssttyp
      intent (out) lu , xlat , xlon
!
! Local variables
!
      real(4) :: truelath , truelatl
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , clat , clon , dsinm , grdfac , plat , plon ,    &
            & ptop , rlatinc , rloninc
      character(3) , dimension(12) :: cmonth
      integer :: i , ibigend , ierr , igrads , iyy , j , jxx , k , kz , &
               & month , nx , ny , period
      character(6) :: iproj
      real(4) , dimension(30) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      read (10,rec=1,iostat=ierr) iyy , jxx , kz , dsinm , clat , clon ,&
                                & plat , plon , grdfac , iproj ,        &
                                & (sigmaf(k),k=1,kz+1) , ptop , igrads ,&
                                & ibigend , truelatl , truelath
      if ( iyy/=iy .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (SST_1DEG.f)'
        print * , '  regcm.in   : ' , iy , jx
        print * , '  DOMAIN.INFO: ' , iyy , jxx
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml1d)'
      end if
      read (10,rec=4,iostat=ierr) ((lu(i,j),j=1,jx),i=1,iy)
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,jx),i=1,iy)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,jx),i=1,iy)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (SST_1DEG.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml1d)'
      end if
!
      if ( igrads==1 ) then
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          do i = 1 , iy
            do j = 1 , jx
              if ( clon>=0.0 ) then
                if ( xlon(i,j)>=0.0 ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)+360.))&
                        & ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else
                  alonmin = amin1(alonmin,xlon(i,j)+360.)
                  alonmax = amax1(alonmax,xlon(i,j)+360.)
                end if
              else if ( xlon(i,j)<0.0 ) then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)-360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else
                alonmin = amin1(alonmin,xlon(i,j)-360.)
                alonmax = amax1(alonmax,xlon(i,j)-360.)
              end if
            end do
          end do
          rlatinc = dsinm*0.001/111./2.
          rloninc = dsinm*0.001/111./2.
          ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = jx/2.
          centeri = iy/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , iy , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , iy , clat , clon , centerj , centeri ,&
                           & truelatl , truelath , clon , dsinm , dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) iy
            write (32,99006) (xlat(i,1),i=1,iy)
          end if 
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , iy , plon , plat , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , iy , plon , plat , dsinm/111000. ,    &
                           & dsinm/111000.*.95238
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        write (31,99008) 1 , 1000.
        month = idate1 - (idate1/100)*100
        period = (idate2/100-idate1/100)*12 + (idate2-(idate2/100)*100) &
               & - (idate1-(idate1/100)*100) + 1
        write (31,99009) period , cmonth(month) , idate1/100
        write (31,99010) 1
        write (31,99011) 'sst ' , 'Sea Surface Temperature    '
        write (31,'(a)') 'endvars'
        close (31)
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          write (32,99008) 1 , 1000.
          write (32,99009) period , cmonth(month) , idate1/100
          write (32,99010) 1
          write (32,99011) 'ice ' , 'Sea Ice fraction           '
          write (32,'(a)') 'endvars'
          close (32)
        end if
      end if
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i1,' levels ',f7.2)
99009 format ('tdef ',i4,' linear 00z16',a3,i4,' 1mo')
99010 format ('vars ',i1)
99011 format (a4,'0 99 ',a26)
!
      end subroutine gridml1d
!
!-----------------------------------------------------------------------
!
      subroutine gridml2(xlon,xlat,lu,iy,jx,idate1,inumber,ibyte,ssttyp)
      implicit none
!
! Dummy arguments
!
      integer :: ibyte , idate1 , iy , jx , inumber
      real(4) , dimension(iy,jx) :: lu , xlat , xlon
      character(5) :: ssttyp
      intent (in) ibyte , idate1 , iy , jx , inumber , ssttyp
      intent (out) lu , xlat , xlon
!
! Local variables
!
      real(4) :: truelath , truelatl
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , clat , clon , dsinm , grdfac , plat , plon ,    &
            & ptop , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: day , i , ibigend , ierr , igrads , iyy , j , jxx , k ,&
               & kz , month , nx , ny
      character(6) :: iproj
      real(4) , dimension(30) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
         & '09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      read (10,rec=1,iostat=ierr) iyy , jxx , kz , dsinm , clat , clon ,&
                                & plat , plon , grdfac , iproj ,        &
                                & (sigmaf(k),k=1,kz+1) , ptop , igrads ,&
                                & ibigend , truelatl , truelath
      if ( iyy/=iy .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (SST_1DEG.f)'
        print * , '  regcm.in   : ' , iy , jx
        print * , '  DOMAIN.INFO: ' , iyy , jxx
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml2)'
      end if
      read (10,rec=4,iostat=ierr) ((lu(i,j),j=1,jx),i=1,iy)
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,jx),i=1,iy)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,jx),i=1,iy)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (SST_1DEG.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml2)'
      end if
!
      if ( igrads==1 ) then
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          do i = 1 , iy
            do j = 1 , jx
              if ( clon>=0.0 ) then
                if ( xlon(i,j)>=0.0 ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)+360.))&
                        & ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else
                  alonmin = amin1(alonmin,xlon(i,j)+360.)
                  alonmax = amax1(alonmax,xlon(i,j)+360.)
                end if
              else if ( xlon(i,j)<0.0 ) then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)-360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else
                alonmin = amin1(alonmin,xlon(i,j)-360.)
                alonmax = amax1(alonmax,xlon(i,j)-360.)
              end if
            end do
          end do
          rlatinc = dsinm*0.001/111./2.
          rloninc = dsinm*0.001/111./2.
          ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = jx/2.
          centeri = iy/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , iy , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , iy , clat , clon , centerj , centeri ,&
                           & truelatl , truelath , clon , dsinm , dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) iy
            write (32,99006) (xlat(i,1),i=1,iy)
          end if 
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , iy , plon , plat , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , iy , plon , plat , dsinm/111000. ,    &
                           & dsinm/111000.*.95238
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        write (31,99008) 1 , 1000.
        month = idate1/100 - (idate1/10000)*100
        day = mod(idate1,100)
        write (31,99009) inumber , cday(day) , cmonth(month) ,          &
                       & idate1/10000
        write (31,99010) 1
        write (31,99011) 'sst ' , 'surface elevation          '
        write (31,'(a)') 'endvars'
        close (31)
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          write (32,99008) 1 , 1000.
          write (32,99009) inumber , cday(day) , cmonth(month) ,        &
                       & idate1/10000
          write (32,99010) 1
          write (32,99011) 'ice ' , 'Sea Ice fraction           '
          write (32,'(a)') 'endvars'
          close (32)
        end if
      end if
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i1,' levels ',f7.2)
99009 format ('tdef ',i4,' linear 00z',a2,a3,i4,' 7dy')
99010 format ('vars ',i1)
99011 format (a4,'0 99 ',a26)
!
      end subroutine gridml2
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      subroutine sst_mn(idate,idate0,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'sst'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = (nyear-1981)*12 + month - 11
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             sst(i,jlat+1-j) = -9999.
          else
             sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do
!
      end subroutine sst_mn
!
!-----------------------------------------------------------------------
!
      subroutine ice_mn(idate,idate0,ilon,jlat,ice,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'ice'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = (nyear-1981)*12 + month - 11
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             ice(i,jlat+1-j) = -9999.
          else
             ice(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do
!
      end subroutine ice_mn
!
!-----------------------------------------------------------------------
!
      subroutine sst_wk(idate,idate0,kkk,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , kkk , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(3) :: varname
      integer :: istatus
      integer(2) , dimension(ilon,jlat) :: work
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'sst'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        istatus = nf90_close(inet)
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = kkk
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             sst(i,jlat+1-j) = -9999.
          else
             sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do

      end subroutine sst_wk
!
!-----------------------------------------------------------------------
!
      subroutine ice_wk(idate,idate0,kkk,ilon,jlat,ice,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , kkk , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(3) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'ice'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        istatus = nf90_close(inet)
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = kkk
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             ice(i,jlat+1-j) = -9999.
          else
             ice(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do

      end subroutine ice_wk
!
      end module mod_sst_1deg
