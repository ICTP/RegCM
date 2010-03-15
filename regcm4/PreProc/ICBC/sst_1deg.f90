!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      program rdsst
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

      use mod_param , only : ix , jx , ssttyp , lsmtyp , ibyte ,        &
            &                idate1 , idate2

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      real(4) , dimension(ilon,jlat) :: sst , ice
      integer :: idate , idate0 , kend , kstart
      integer , dimension(427+1045) :: wkday
      real , dimension(ix,jx) :: lu , sstmm , icemm , xlat , xlon
      integer :: i , idatef , idateo , j , k , ludom , lumax , mrec ,   &
               & nday , nmo , nrec , nyear
      real , dimension(jlat) :: lati
      real , dimension(ilon) :: loni
      integer , dimension(25) :: lund
      real :: truelath , truelatl
      logical :: there
!
      if ( ssttyp=='GISST' ) then
        if ( idate1<1947121512 .or. idate2>2002091512 ) then
          print * , 'GISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , idate1 , idate2
          stop
        end if
        open (11,file='../DATA/SST/GISST_194712_200209',                &
             &form='unformatted',recl=360*180*ibyte,access='direct',    &
            & status='old',err=100)
      else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.               &
            &   ssttyp=='OI2ST' ) then
        if ( idate1<1981121512 .or. idate2<1981121512 ) then
          print * , 'OISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , idate1 , idate2
          stop
        end if
        inquire (file='../DATA/SST/sst.mnmean.nc',exist=there)
        if ( .not.there ) print * , 'sst.mnmean.nc is not available' ,  &
                               &' under ../DATA/SST/'
        if ( ssttyp=='OI2ST' ) then
          inquire (file='../DATA/SST/icec.mnmean.nc',exist=there)
          if ( .not. there )                                            &
            & print *, 'icec.mnmean.nc is not available' ,              &
                   &' under ../DATA/SST/'
        end if
      else if ( ssttyp=='OI_WK' .or. ssttyp=='OI2WK' ) then
        if ( idate1<1981110100 .or. idate2<1981110106 ) then
          print * , 'OI_WK (or OI2WK) data required are not available'
          print * , 'IDATE1, IDATE2 = ' , idate1 , idate2
          stop
        end if
        inquire (file='../DATA/SST/sst.wkmean.1981-1989.nc',exist=there)
        if ( .not.there ) print * ,                                     &
                             &'sst.wkmean.1981-1989.nc is not available'&
                            & , ' under ../DATA/SST/'
        inquire (file='../DATA/SST/sst.wkmean.1990-present.nc',         &
               & exist=there)
        if ( .not.there ) print * ,                                     &
                          &'sst.wkmean.1990-present.nc is not available'&
                         & , ' under ../DATA/SST/'
        call headwk(wkday)
        if ( ssttyp=='OI2WK' ) then
          inquire (file='../DATA/SST/icec.wkmean.1981-1989.nc',         &
                & exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1981-1989.nc is not available',       &
                &    ' under ../DATA/SST/'
          end if
          inquire (file='../DATA/SST/icec.wkmean.1990-present.nc',      &
                &  exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1990-present.nc is not available',    &
                &    ' under ../DATA/SST/'
          end if
        end if
      else
        write (*,*) 'PLEASE SET SSTTYP in domain.param'
        stop
      end if
      open (21,file='SST.RCM',form='unformatted',status='replace')
 
!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      open (10,file='../../Input/DOMAIN.INFO',form='unformatted',       &
          & recl=ix*jx*ibyte,access='direct',status='unknown',err=200)
 
!#####
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!#####
        idate = idate1/10000
        if ( idate-(idate/100)*100==1 ) then
          idate = idate - 89
        else
          idate = idate - 1
        end if
        idateo = idate
        idate0 = idateo*10000 + 100
        idate = idate2/10000
        if ( idate-(idate/100)*100==12 ) then
          idate = idate + 89
        else
          idate = idate + 1
        end if
        idatef = idate
        print * , idate1 , idate2 , idateo , idatef
        call gridml(xlon,xlat,lu,ix,jx,idateo,idatef,ibyte,truelatl,    &
                  & truelath,ssttyp)
!#####
      else
!#####
        idate = idate1/100
        do k = 427 + 1045 , 1 , -1
          if ( wkday(k)<=idate ) then
            kstart = k
            exit
          end if
        end do
        idate = idate2/100
        do k = 1 , 427 + 1045
          if ( wkday(k)>idate ) then
            kend = k
            exit
          end if
        end do
        idateo = wkday(kstart)
        idate0 = wkday(kstart)*100
        idatef = wkday(kend)
        print * , idate1 , idate2 , idateo , idatef , kend - kstart + 1
        call gridml2(xlon,xlat,lu,ix,jx,idateo,kend-kstart+1,ibyte,     &
                   & truelatl,truelath,ssttyp)
!#####
      end if
!#####
      open (25,file='RCM_SST.dat',status='unknown',form='unformatted',  &
          & recl=ix*jx*ibyte,access='direct')
      if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
        open (26,file='RCM_ICE.dat',status='unknown',form='unformatted',&
          & recl=ix*jx*ibyte,access='direct')

      end if
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
        do while ( idate<=idatef )
          nyear = idate/100
          nmo = idate - nyear*100
          if ( ssttyp=='GISST' ) then
            nrec = (nyear-1947)*12 + nmo - 11
            read (11,rec=nrec) sst
          else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.           &
                    ssttyp=='OI2ST') then
            write (*,*) idate*10000 + 100 , idate0
            call sst_mn(idate*10000+100,idate0,ilon,jlat,sst)
            if ( ssttyp=='OI2ST' )                                      &
             &  call ice_mn(idate*10000+100,idate0,ilon,jlat,ice)
          else
          end if
 
!         ******           PRINT OUT DATA AS A CHECK
          if ( nmo==1 ) call printl(sst,ilon,jlat)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,ix,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
 
          if ( ssttyp=='OI2ST' ) call bilinx(ice,icemm,xlon,xlat,       &
              & loni,lati,ilon,jlat,ix,jx,1)
          do j = 1 , jx
            do i = 1 , ix
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
          idate = idate + 1
          if ( nmo==12 ) idate = idate + 88
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,ix)
          if ( ssttyp=='OI2ST' )                                        &
           &    write (26,rec=mrec) ((icemm(i,j),j=1,jx),i=1,ix)
        end do
        write (10,rec=4) ((lu(i,j),j=1,jx),i=1,ix)
!#####
      else
!#####
        do k = kstart , kend
          idate = wkday(k)
          nyear = idate/10000
          nmo = idate/100 - nyear*100
          nday = mod(idate,100)
          write (*,*) idate*100 , idate0 , k
          call sst_wk(idate*100,idate0,k,ilon,jlat,sst)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,ix,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
          if ( ssttyp=='OI2WK') then
            call ice_wk(idate*100,idate0,k)
            call bilinx(ice,icemm,xlon,xlat,loni,lati,ilon,jlat,ix,jx,1)
          end if 

          do j = 1 , jx
            do i = 1 , ix
              if ( sstmm(i,j)>-100. ) then
                sstmm(i,j) = sstmm(i,j) + 273.15
              else
                sstmm(i,j) = -9999.
              end if
            end do
          end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
          write (21) nday , nmo , nyear , sstmm
          idate = wkday(k)
          print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear , idate
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,ix)
        end do
 
!#####
      end if
!#####
 
      stop 99999
 100  continue
      print * , 'ERROR OPENING GISST FILE'
      stop '4810 IN PROGRAM RDSST'
 200  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'
      end program rdsst
!
!-----------------------------------------------------------------------
!
      subroutine gridml(xlon,xlat,lu,ix,jx,idate1,idate2,ibyte,truelatl,&
                      & truelath,ssttyp)
      implicit none
!
! Dummy arguments
!
      integer :: ibyte , idate1 , idate2 , ix , jx
      real :: truelath , truelatl
      real , dimension(ix,jx) :: lu , xlat , xlon
      character(5) :: ssttyp
      intent (in) ibyte , idate1 , idate2 , ix , jx , ssttyp
      intent (out) lu
      intent (inout) truelath , truelatl , xlat , xlon
!
! Local variables
!
      real :: alatmax , alatmin , alonmax , alonmin , centeri ,         &
            & centerj , clat , clon , dsinm , grdfac , plat , plon ,    &
            & ptop , rlatinc , rloninc
      character(3) , dimension(12) :: cmonth
      integer :: i , ibigend , ierr , igrads , ixy , j , jxx , k , kz , &
               & month , nx , ny , period
      character(6) :: iproj
      real , dimension(30) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      read (10,rec=1,iostat=ierr) ixy , jxx , kz , dsinm , clat , clon ,&
                                & plat , plon , grdfac , iproj ,        &
                                & (sigmaf(k),k=1,kz+1) , ptop , igrads ,&
                                & ibigend , truelatl , truelath
      if ( ixy/=ix .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (SST_1DEG.f)'
        print * , '  icbc.param: ' , ix , jx
        print * , '  DOMAIN.INFO: ' , ixy , jxx
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml)'
      end if
      read (10,rec=4,iostat=ierr) ((lu(i,j),j=1,jx),i=1,ix)
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,jx),i=1,ix)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,jx),i=1,ix)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (SST_1DEG.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml)'
      end if
!
      if ( igrads==1 ) then
        open (31,file='RCM_SST.ctl',status='replace')
        write (31,'(a)') 'dset ^RCM_SST.dat'
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           open ( 32,file='RCM_ICE.ctl',status='replace')
           write (32,'(a)') 'dset ^RCM_ICE.dat'
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          alatmin = 999999.
          alatmax = -999999.
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(ix,j)>alatmax ) alatmax = xlat(ix,j)
          end do
          alonmin = 999999.
          alonmax = -999999.
          do i = 1 , ix
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
          centeri = ix/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , ix , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , ix , clat , clon , centerj , centeri ,&
                           & truelatl , truelath , clon , dsinm , dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) ix
          write (31,99006) (xlat(i,1),i=1,ix)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) ix
            write (32,99006) (xlat(i,1),i=1,ix)
          end if 
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , ix , plon , plat , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , ix , plon , plat , dsinm/111000. ,    &
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
      end subroutine gridml
!
!-----------------------------------------------------------------------
!
      subroutine gridml2(xlon,xlat,lu,ix,jx,idate1,inumber,ibyte,       &
                       & truelatl,truelath,ssttyp)
      implicit none
!
! Dummy arguments
!
      integer :: ibyte , idate1 , ix , jx , inumber
      real :: truelath , truelatl
      real , dimension(ix,jx) :: lu , xlat , xlon
      character(5) :: ssttyp
      intent (in) ibyte , idate1 , ix , jx , inumber , ssttyp
      intent (out) lu
      intent (inout) truelath , truelatl , xlat , xlon
!
! Local variables
!
      real :: alatmax , alatmin , alonmax , alonmin , centeri ,         &
            & centerj , clat , clon , dsinm , grdfac , plat , plon ,    &
            & ptop , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: day , i , ibigend , ierr , igrads , ixy , j , jxx , k ,&
               & kz , month , nx , ny
      character(6) :: iproj
      real , dimension(30) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
!
      read (10,rec=1,iostat=ierr) ixy , jxx , kz , dsinm , clat , clon ,&
                                & plat , plon , grdfac , iproj ,        &
                                & (sigmaf(k),k=1,kz+1) , ptop , igrads ,&
                                & ibigend , truelatl , truelath
      if ( ixy/=ix .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (SST_1DEG.f)'
        print * , '  icbc.param: ' , ix , jx
        print * , '  DOMAIN.INFO: ' , ixy , jxx
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml2)'
      end if
      read (10,rec=4,iostat=ierr) ((lu(i,j),j=1,jx),i=1,ix)
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,jx),i=1,ix)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,jx),i=1,ix)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (SST_1DEG.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml2)'
      end if
!
      if ( igrads==1 ) then
        open (31,file='RCM_SST.ctl',status='replace')
        write (31,'(a)') 'dset ^RCM_SST.dat'
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           open ( 32,file='RCM_ICE.ctl',status='replace')
           write (32,'(a)') 'dset ^RCM_ICE.dat'
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          alatmin = 999999.
          alatmax = -999999.
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(ix,j)>alatmax ) alatmax = xlat(ix,j)
          end do
          alonmin = 999999.
          alonmax = -999999.
          do i = 1 , ix
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
          centeri = ix/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , ix , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , ix , clat , clon , centerj , centeri ,&
                           & truelatl , truelath , clon , dsinm , dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) ix
          write (31,99006) (xlat(i,1),i=1,ix)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) ix
            write (32,99006) (xlat(i,1),i=1,ix)
          end if 
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , ix , plon , plat , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , ix , plon , plat , dsinm/111000. ,    &
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
      subroutine headwk(wkday)
      implicit none
!
! Dummy arguments
!
      integer , intent(out) , dimension(427+1045) :: wkday
!
! Local variables
!
      integer :: i , mday , month , myear
!
      wkday(1) = 19811029
      do i = 2 , 427
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      wkday(428) = 19891231
      do i = 429 , 427 + 1045
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      end subroutine headwk
!
!-----------------------------------------------------------------------
!
      subroutine sst_mn(idate,idate0,ilon,jlat,sst)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      intent (in) idate , idate0 , ilon , jlat
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
      character(35) :: pathaddname
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
!
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
!bxq
!bxq_
      data varname/'sst'/
!
      if ( idate==idate0 ) then
        pathaddname = '../DATA/SST/sst.mnmean.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open('../DATA/SST/sst.mnmean.nc',                &
                &           nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
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
      istatus = nf90_get_var(inet,5,work,istart,icount)
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
      istatus = nf90_close(inet)

      end subroutine sst_mn
!
!-----------------------------------------------------------------------
!
      subroutine ice_mn(idate,idate0,ilon,jlat,ice)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      intent (in) idate , idate0 , ilon , jlat
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
      character(35) :: pathaddname
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
!
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
!bxq
!bxq_
      data varname/'ice'/
!
      if ( idate==idate0 ) then
        pathaddname = '../DATA/SST/icec.mnmean.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open('../DATA/SST/icec.mnmean.nc',               &
                &           nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
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
      istatus = nf90_get_var(inet,5,work,istart,icount)
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
      istatus = nf90_close(inet)

      end subroutine ice_mn
!
!-----------------------------------------------------------------------
!
      subroutine sst_wk(idate,idate0,kkk,ilon,jlat,sst)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      intent (in) idate , idate0 , kkk , ilon , jlat
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
      character(38) :: pathaddname
      logical :: there
      character(3) :: varname
      integer(2) , dimension(ilon,jlat) :: work
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
!bxq
!bxq_
      data varname/'sst'/
!
      if ( idate==idate0 ) then
        if ( idate<1989123100 ) then
          pathaddname = '../DATA/SST/sst.wkmean.1981-1989.nc'
        else
          pathaddname = '../DATA/SST/sst.wkmean.1990-present.nc'
        end if
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        pathaddname = '../DATA/SST/sst.wkmean.1990-present.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
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
      if ( kkk>427 ) it = kkk - 427
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,5,work,istart,icount)
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

      istatus = nf90_close(inet)
 
      end subroutine sst_wk
!
!-----------------------------------------------------------------------
!
      subroutine ice_wk(idate,idate0,kkk,ilon,jlat,ice)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      intent (in) idate , idate0 , kkk , ilon , jlat
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
      character(64) :: pathaddname
      logical :: there
      character(3) :: varname
      integer(2) , dimension(ilon,jlat) :: work
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
!bxq
!bxq_
      data varname/'ice'/
!
      if ( idate==idate0 ) then
        if ( idate<1989123100 ) then
          pathaddname = '../DATA/SST/icec.wkmean.1981-1989.nc'
        else
          pathaddname = '../DATA/SST/icec.wkmean.1990-present.nc'
        end if
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        pathaddname = '../DATA/SST/icec.wkmean.1990-present.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        istatus = nf90_get_att(inet,5,'scale_factor',xscale)
        istatus = nf90_get_att(inet,5,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 360
        icount(2) = 180
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
      if ( kkk>427 ) it = kkk - 427
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,5,work,istart,icount)
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

      istatus = nf90_close(inet)
 
      end subroutine ice_wk
