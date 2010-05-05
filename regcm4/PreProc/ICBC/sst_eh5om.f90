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

      program sst_eh5om

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! EH5OM    EH5OM_SST is from the coupled model run by EC-MPI         !
!          6 hourly frequncy, 1.875x1.875                            !
!          for 'RF'       run , from 1941 to 2000,                   !
!          for 'A2'/'B2'/'A1B', from 2001 to 2100.                   !
!                                                                    !
!          ML= 1 is   0.0; ML= 2 is   1.875; => ML=192 is 358.125E   !
!          NL= 1 is  90.0; ML= 2 is  88.75 ; => ML=96 is -89.0625    !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mod_dynparam
      use mod_datenum

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 192 , jlat = 96
!
! Local variables
!
      integer :: it_base
      integer(2) , dimension(192,96) :: ivar
      real(8) :: offset , xscale
      real(4) , dimension(ilon,jlat) :: sst
      integer :: idate
      integer :: nnnend , nstart
      integer :: i , it , j , mrec , nday , nhour , nmo , nyear , ierr
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      logical :: there
      character(256) :: namelistfile, prgname
      real(4) , allocatable , dimension(:,:) :: lu , sstmm , xlat , xlon
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
      allocate(lu(iy,jx))
      allocate(sstmm(iy,jx))
      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))

      if ( ssttyp=='EH5RF' ) then
        there = .false.
        if ( (globidate1>=1941010106 .and. globidate1<=1961123118) .or. &
           & (globidate2>=1941010106 .and. globidate2<=1961123118) ) then
          inquire (file='../DATA/SST/SST_20C_3_1941010106_1961123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1941010106_1961123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=1962010100 .and. globidate1<=1993123118) .or. &
           & (globidate2>=1962010100 .and. globidate2<=1993123118) ) then
          inquire (file='../DATA/SST/SST_20C_3_1962010100_1993123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1962010100_1993123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=1994010100 .and. globidate1<=2001010100) .or. &
           & (globidate2>=1994010100 .and. globidate2<=2001010100) ) then
          inquire (file='../DATA/SST/SST_20C_3_1994010100_2001010100',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1994010100_2001010100 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'EH5RF SST is just available from 1941010106 to' ,  &
               &' 2001010100'
          stop
        end if
      else if ( ssttyp=='EH5A2' ) then
        there = .false.
        if ( (globidate1>=2001010100 .and. globidate1<=2029123118) .or. &
           & (globidate2>=2001010100 .and. globidate2<=2029123118) ) then
          inquire (file='../DATA/SST/SST_A2_1_2001010100_2029123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2001010100_2029123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file='../DATA/SST/SST_A2_1_2030010100_2061123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2030010100_2061123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file='../DATA/SST/SST_A2_1_2062010100_2093123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2062010100_2093123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file='../DATA/SST/SST_A2_1_2094010100_2100123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2094010100_2100123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'EH5A2 SST is just available from 2001010100 to' ,  &
               &' 2100123118'
          stop
        end if
      else if ( ssttyp=='EH5B1' ) then
        there = .false.
        if ( (globidate1>=2001010100 .and. globidate1<=2029123118) .or. &
           & (globidate2>=2001010100 .and. globidate2<=2029123118) ) then
          inquire (file='../DATA/SST/SST_B1_1_2001010100_2029123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2001010100_2029123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file='../DATA/SST/SST_B1_1_2030010100_2061123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2030010100_2061123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file='../DATA/SST/SST_B1_1_2062010100_2093123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2062010100_2093123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file='../DATA/SST/SST_B1_1_2094010100_2100123118',   &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2094010100_2100123118 is not available' &
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'EH5B1 SST is just available from 2001010100 to' ,  &
               &' 2100123118'
          stop
        end if
      else if ( ssttyp=='EHA1B' ) then
        there = .false.
        if ( (globidate1>=2001010100 .and. globidate1<=2029123118) .or. &
           & (globidate2>=2001010100 .and. globidate2<=2029123118) ) then
          inquire (file='../DATA/SST/SST_A1B_3_2001010100_2029123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2001010100_2029123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file='../DATA/SST/SST_A1B_3_2030010100_2061123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2030010100_2061123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file='../DATA/SST/SST_A1B_3_2062010100_2093123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2062010100_2093123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file='../DATA/SST/SST_A1B_3_2094010100_2100123118',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2094010100_2100123118 is not available'&
                & , ' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'EHA1B SST is just available from 2001010100 to' ,  &
               &' 2100123118'
          stop
        end if
      else
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported types are EH5RF EH5A2 EH5B1 EHA1B'
        stop
      end if
      open (21,file='SST.RCM',form='unformatted',status='replace')
 
!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      open (10,file=terfilout,form='unformatted',recl=iy*jx*ibyte,      &
         &  access='direct',status='unknown',err=100)
      call initdate_eh50(ssttyp)
      call finddate_eh50(nstart,globidate1)
      call finddate_eh50(nnnend,globidate2)
      write (*,*) nstart , nnnend
      print * , globidate1 , nnnend - nstart + 1
      call gridml(xlon,xlat,lu,iy,jx,globidate1,nnnend-nstart+1)
      open (25,file='RCM_SST.dat',status='unknown',form='unformatted',  &
          & recl=iy*jx*ibyte,access='direct')
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.875
      end do
      do j = 1 , jlat
        lati(j) = -89.0625 + 1.875*float(j-1)
      end do
 
!     **  REF  SST DATA, 1.875x1.1.25, AVAILABLE FROM 16/1/1959 TO
!     16/1/1991 ** A2&B2 SST DATA, 1.875x1.1.25, AVAILABLE FROM
!     16/1/2069 TO 16/1/2101
      do it = nstart , nnnend
        idate = mdate(it)
        if ( ssttyp=='EH5RF' ) then
          if ( idate>=1941010106 .and. idate<=1961123118 ) then
            open (11,file='../DATA/SST/SST_20C_3_1941010106_1961123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 1
          else if ( idate>=1962010100 .and. idate<=1993123118 ) then
            open (11,file='../DATA/SST/SST_20C_3_1962010100_1993123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 30680
          else if ( idate>=1994010100 .and. idate<=2001010100 ) then
            open (11,file='../DATA/SST/SST_20C_3_1994010100_2001010100',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 30680 + 46752
          else
          end if
        else if ( ssttyp=='EH5A2' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file='../DATA/SST/SST_A2_1_2001010100_2029123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file='../DATA/SST/SST_A2_1_2030010100_2061123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file='../DATA/SST/SST_A2_1_2062010100_2093123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file='../DATA/SST/SST_A2_1_2094010100_2100123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else if ( ssttyp=='EH5B1' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file='../DATA/SST/SST_B1_1_2001010100_2029123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file='../DATA/SST/SST_B1_1_2030010100_2061123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file='../DATA/SST/SST_B1_1_2062010100_2093123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file='../DATA/SST/SST_B1_1_2094010100_2100123118', &
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else if ( ssttyp=='EHA1B' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file='../DATA/SST/SST_A1B_3_2001010100_2029123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file='../DATA/SST/SST_A1B_3_2030010100_2061123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file='../DATA/SST/SST_A1B_3_2062010100_2093123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file='../DATA/SST/SST_A1B_3_2094010100_2100123118',&
                & form='unformatted',recl=(192*96/2+4)*ibyte,           &
                 &access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else
        end if
        nyear = idate/1000000
        nmo = (idate-nyear*1000000)/10000
        nday = (idate-nyear*1000000-nmo*10000)/100
        nhour = idate - nyear*1000000 - nmo*10000 - nday*100
        read (11,rec=it-it_base) offset , xscale , ivar
        write (*,*) offset , xscale
        do j = 1 , jlat
          do i = 1 , ilon
            sst(i,j) = ivar(i,jlat+1-j)*xscale + offset
            if ( sst(i,j)<273.16 ) sst(i,j) = -9999.
          end do
        end do
 
!       ******           PRINT OUT DATA AS A CHECK
        if ( nmo==1 ) call printl(sst,ilon,jlat)
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        write (21) nhour , nday , nmo , nyear , sstmm
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
        mrec = mrec + 1
        write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
      end do
 
      deallocate(lu)
      deallocate(sstmm)
      deallocate(xlat)
      deallocate(xlon)

      stop 99999
!     4810 PRINT *,'ERROR OPENING GISST FILE'
!     STOP '4810 IN PROGRAM RDSST'
 100  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'
      end program sst_eh5om
!
!-----------------------------------------------------------------------
!
      subroutine gridml(xlon,xlat,lu,iy,jx,idate1,numrec)
      implicit none
!
! Dummy arguments
!
      integer :: idate1 , iy , jx , numrec
      real(4) , dimension(iy,jx) :: lu , xlat , xlon
      intent (in) idate1 , iy , jx , numrec
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
      integer :: i , ibigend , igrads , iyy , j , jxx , k , kz ,        &
               & month , nday , nhour , nx , ny , nyear , period
      character(6) :: iproj
      real(4) , dimension(30) :: sigmaf
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      read (10,rec=1) iyy , jxx , kz , dsinm , clat , clon , plat ,     &
                    & plon , grdfac , iproj , (sigmaf(k),k=1,kz+1) ,    &
                    & ptop , igrads , ibigend , truelatl , truelath
      if ( iyy/=iy .or. jxx/=jx ) then
        write (*,*) 'IY,JX,IYY,JXX' , iy , jx , iyy , jxx
        stop
      end if
      read (10,rec=4) ((lu(i,j),j=1,jx),i=1,iy)
      read (10,rec=5) ((xlat(i,j),j=1,jx),i=1,iy)
      read (10,rec=6) ((xlon(i,j),j=1,jx),i=1,iy)
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
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          alatmin = 999999.
          alatmax = -999999.
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          alonmin = 999999.
          alonmax = -999999.
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
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
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
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        write (31,99008) 1 , 1000.
        period = numrec
        nyear = idate1/1000000
        month = (idate1-nyear*1000000)/10000
        nday = (idate1-nyear*1000000-month*10000)/100
        nhour = mod(idate1,100)
        write (31,99009) period , nhour , cday(nday) , cmonth(month) ,  &
                       & nyear
        write (31,99010) 1
        write (31,99011) 'sst ' , 'surface elevation          '
        write (31,'(a)') 'endvars'
        close (31)
      end if
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i1,' levels ',f7.2)
99009 format ('tdef ',i6,' linear ',i2,'z',a2,a3,i4,' 6hr')
99010 format ('vars ',i1)
99011 format (a4,'0 99 ',a26)
!
      end subroutine gridml
