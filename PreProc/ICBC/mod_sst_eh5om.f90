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

      module mod_sst_eh5om

      contains

      subroutine sst_eh5om

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

      use mod_sst_grid
      use mod_date
      use mod_interp , only : bilinx
      use netcdf

      implicit none
!
      integer , parameter :: ilon = 192 , jlat = 96
      integer , parameter :: idtbc = 6
!
      integer :: it_base
      integer(2) , dimension(ilon,jlat) :: ivar
      real(8) :: offset , xscale
      real(4) , dimension(ilon,jlat) :: sst
      integer :: idate , ieh5ostart , ieh5orec, nsteps
      integer :: i , it , j , nday , nhour , nmo , nyear
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      logical :: there
!
      it_base = 0

      if ( ssttyp=='EH5RF' ) then
        there = .false.
        if ( (globidate1>=1941010106 .and. globidate1<=1961123118) .or. &
           & (globidate2>=1941010106 .and. globidate2<=1961123118) ) then
          inquire (file=trim(inpglob)//                                 &
                   '/SST/SST_20C_3_1941010106_1961123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1941010106_1961123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=1962010100 .and. globidate1<=1993123118) .or. &
           & (globidate2>=1962010100 .and. globidate2<=1993123118) ) then
          inquire (file=trim(inpglob)//                                 &
               &   '/SST/SST_20C_3_1962010100_1993123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1962010100_1993123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=1994010100 .and. globidate1<=2001010100) .or. &
           & (globidate2>=1994010100 .and. globidate2<=2001010100) ) then
          inquire (file=trim(inpglob)//                                 &
               &   '/SST/SST_20C_3_1994010100_2001010100',exist=there)
          if ( .not.there ) then
            print * , 'SST_20C_3_1994010100_2001010100 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
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
          inquire (file=trim(inpglob)//                                 &
                 & '/SST/SST_A2_1_2001010100_2029123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2001010100_2029123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file=trim(inpglob)//                                 &
                 & '/SST/SST_A2_1_2030010100_2061123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2030010100_2061123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file=trim(inpglob)//                                 &
           &  '/SST/SST_A2_1_2062010100_2093123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2062010100_2093123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file=trim(inpglob)//                                 &
           &       '/SST/SST_A2_1_2094010100_2100123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A2_1_2094010100_2100123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
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
          inquire (file=trim(inpglob)//                                 &
                &  '/SST/SST_B1_1_2001010100_2029123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2001010100_2029123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file=trim(inpglob)//                                 &
             &     '/SST/SST_B1_1_2030010100_2061123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2030010100_2061123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file=trim(inpglob)//                                 &
              &    '/SST/SST_B1_1_2062010100_2093123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2062010100_2093123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file=trim(inpglob)//                                 &
             &     '/SST/SST_B1_1_2094010100_2100123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_B1_1_2094010100_2100123118 is not available' &
                & , ' under ',trim(inpglob),'/SST/'
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
          inquire (file=trim(inpglob)//                                 &
             &     '/SST/SST_A1B_3_2001010100_2029123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2001010100_2029123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2030010100 .and. globidate1<=2061123118) .or. &
           & (globidate2>=2030010100 .and. globidate2<=2061123118) ) then
          inquire (file=trim(inpglob)//                                 &
                &  '/SST/SST_A1B_3_2030010100_2061123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2030010100_2061123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2062010100 .and. globidate1<=2093123118) .or. &
           & (globidate2>=2062010100 .and. globidate2<=2093123118) ) then
          inquire (file=trim(inpglob)//                                 &
              &    '/SST/SST_A1B_3_2062010100_2093123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2062010100_2093123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( (globidate1>=2094010100 .and. globidate1<=2100123118) .or. &
           & (globidate2>=2094010100 .and. globidate2<=2100123118) ) then
          inquire (file=trim(inpglob)//                                 &
                 & '/SST/SST_A1B_3_2094010100_2100123118',exist=there)
          if ( .not.there ) then
            print * , 'SST_A1B_3_2094010100_2100123118 is not available'&
                & , ' under ',trim(inpglob),'/SST/'
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

      nsteps = idatediff(globidate2,globidate1)/idtbc + 1
      write (*,*) 'GLOBIDATE1 : ' , globidate1
      write (*,*) 'GLOBIDATE2 : ' , globidate2
      write (*,*) 'NSTEPS     : ' , nsteps

      call open_sstfile(globidate1)
 
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
      idate = globidate1
      if ( ssttyp=='EH5RF' ) then
        ieh5ostart = 1989010100
      else
        ieh5ostart = 2001010100
      end if
      do it = 1 , nsteps
        
        if ( ssttyp=='EH5RF' ) then
          if ( idate>=1941010106 .and. idate<=1961123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_20C_3_1941010106_1961123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 1
          else if ( idate>=1962010100 .and. idate<=1993123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_20C_3_1962010100_1993123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 30680
          else if ( idate>=1994010100 .and. idate<=2001010100 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_20C_3_1994010100_2001010100',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 30680 + 46752
          else
          end if
        else if ( ssttyp=='EH5A2' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A2_1_2001010100_2029123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A2_1_2030010100_2061123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A2_1_2062010100_2093123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A2_1_2094010100_2100123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else if ( ssttyp=='EH5B1' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_B1_1_2001010100_2029123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_B1_1_2030010100_2061123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_B1_1_2062010100_2093123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_B1_1_2094010100_2100123118',                &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else if ( ssttyp=='EHA1B' ) then
          if ( idate>=2001010100 .and. idate<=2029123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A1B_3_2001010100_2029123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 0
          else if ( idate>=2030010100 .and. idate<=2061123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A1B_3_2030010100_2061123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368
          else if ( idate>=2062010100 .and. idate<=2093123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A1B_3_2062010100_2093123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752
          else if ( idate>=2094010100 .and. idate<=2100123118 ) then
            open (11,file=trim(inpglob)//                               &
             &    '/SST/SST_A1B_3_2094010100_2100123118',               &
             &    form='unformatted',recl=(ilon*jlat/2+4)*ibyte,        &
             &    access='direct')
            it_base = 42368 + 46752*2
          else
          end if
        else
        end if
        call split_idate(idate, nyear, nmo, nday, nhour)
        ieh5orec = idatediff(idate, ieh5ostart)/dble(ibdyfrq)
        read (11,rec=ieh5orec-it_base) offset , xscale , ivar
        do j = 1 , jlat
          do i = 1 , ilon
            sst(i,j) = real(dble(ivar(i,jlat+1-j))*xscale + offset)
            if ( sst(i,j)<273.16 ) sst(i,j) = -9999.
          end do
        end do
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        call writerec(idate,.false.)
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear

        call addhours(idate, idtbc)

      end do
 
      end subroutine sst_eh5om
!
      end module mod_sst_eh5om
