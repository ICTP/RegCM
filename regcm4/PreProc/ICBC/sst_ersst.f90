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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! ERAIN    ERAIN_SST is provided by ERA-Interim project              !
!          6 hourly frequncy, 1.5x1.5 degree resolution              !
!          from 1989010100 to 2009053118.                            !
!          'ERSST' for using the sea surface temperature;            !
!          'ERSKT' for using the skin temperature.                   !
!                                                                    !
!          ML= 1 is   0.0; ML= 2 is   1.5; => ML=240 is 358.5E       !
!          NL= 1 is  90.0; ML= 2 is  88.5; => ML=121 is -90.         !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use mod_param
      use mod_datenum
      implicit none
!
! PARAMETERS
!
      integer , parameter :: ilon = 240 , jlat = 121
!
! Local variables
!
      integer :: i , it , j , mrec , nday , nhour , nmo , nyear
      real , dimension(jlat) :: lati
      real , dimension(ilon) :: loni
      integer :: idate
      logical :: there
      integer :: nnnend , nstart
      real(4) , dimension(ilon,jlat) :: sst
      real , dimension(iy,jx) :: lu , sstmm , xlat , xlon
!
      if ( ssttyp=='ERSST' ) then
        there = .false.
        if ( (idate1>=1989010100 .and. idate1<=2009053118) .or.         &
           & (idate2>=1989010100 .and. idate2<=2009053118) ) then
          inquire (file='../DATA/SST/sstERAIN.1989-2009.nc',exist=there)
          if ( .not.there ) then
            print * , 'sstERAIN.1989-2009.nc is not available' ,        &
                 &' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'ERSST Sea Surface Temperature is just available' , &
               &' from 1989010100 to 2009053118'
          stop
        end if
      else if ( ssttyp=='ERSKT' ) then
        there = .false.
        if ( (idate1>=1989010100 .and. idate1<=2009053118) .or.         &
           & (idate2>=1989010100 .and. idate2<=2009053118) ) then
          inquire (file='../DATA/SST/tskinERAIN.1989-2009.nc',          &
                 & exist=there)
          if ( .not.there ) then
            print * , 'tskinERAIN.1989-2009.nc is not available' ,      &
                 &' under ../DATA/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'ERSKT Skin Temperature is just available' ,        &
               &' from 1989010100 to 2009053118'
          stop
        end if
      else
        write (*,*) 'PLEASE SET the right SSTTYP in domain.param'
        stop
      end if
      open (21,file='SST.RCM',form='unformatted',status='replace')
 
!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      open (10,file='../../Input/DOMAIN.INFO',form='unformatted',       &
          & recl=iy*jx*ibyte,access='direct',status='unknown',err=100)
      call initdate_era
      call finddate_era(nstart,idate1)
      call finddate_era(nnnend,idate2)
      write (*,*) nstart , nnnend
      print * , idate1 , nnnend - nstart + 1
      call gridml(xlon,xlat,lu,iy,jx,idate1,nnnend-nstart+1,truelatl,   &
                & truelath)
      open (25,file='RCM_SST.dat',status='unknown',form='unformatted',  &
          & recl=iy*jx*ibyte,access='direct')
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.5
      end do
      do j = 1 , jlat
        lati(j) = -90. + 1.5*float(j-1)
      end do
 
      do it = nstart , nnnend
        idate = mdate(it)
        if ( ssttyp=='ERSST' ) then
          call sst_erain(it,nstart,ilon,jlat,sst)
        else if ( ssttyp=='ERSKT' ) then
          call skt_erain(it,nstart,ilon,jlat,sst)
        else
        end if
        nyear = idate/1000000
        nmo = (idate-nyear*1000000)/10000
        nday = (idate-nyear*1000000-nmo*10000)/100
        nhour = idate - nyear*1000000 - nmo*10000 - nday*100
 
!       ******           PRINT OUT DATA AS A CHECK
        if ( nmo==1 ) call printl(sst,ilon,jlat)
 
        call bilinx(sst,loni,lati,ilon,jlat,sstmm,xlon,xlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        write (21) nhour , nday , nmo , nyear , sstmm
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
        mrec = mrec + 1
        write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
      end do
 
      stop 99999
!     4810 PRINT *,'ERROR OPENING GISST FILE'
!     STOP '4810 IN PROGRAM RDSST'
 100  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'
      end program rdsst
!
!-----------------------------------------------------------------------
!
      subroutine gridml(xlon,xlat,lu,iy,jx,idate1,numrec,truelatl,      &
                      & truelath)
      implicit none
!
! Dummy arguments
!
      integer :: idate1 , iy , jx , numrec
      real :: truelath , truelatl
      real , dimension(iy,jx) :: lu , xlat , xlon
      intent (in) idate1 , iy , jx , numrec
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
      integer :: i , ibigend , igrads , iyy , j , jxx , k , kz ,        &
               & month , nday , nhour , nx , ny , nyear , period
      character(6) :: iproj
      real , dimension(30) :: sigmaf
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
!
!-----------------------------------------------------------------------
!
      subroutine skt_erain(it,it0,ilon,jlat,sst)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , it0 , ilon , jlat
      intent (in) it , it0 , ilon , jlat
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
      character(35) :: pathaddname
      logical :: there
      character(3) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
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
      data varname/'skt'/
!
      if ( it==it0 ) then
        pathaddname = '../DATA/SST/tskinERAIN.1989-2009.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open('../DATA/SST/tskinERAIN.1989-2009.nc',       &
               & nf90_nowrite,inet)
        istatus = nf90_get_att(inet,4,'scale_factor',xscale)
        istatus = nf90_get_att(inet,4,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,4,work,istart,icount)
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          sst(i,jlat+1-j) = work(i,j)*xscale + xadd
        end do
      end do
!
      istatus = nf90_close(inet)

      end subroutine skt_erain
!
!-----------------------------------------------------------------------
!
      subroutine sst_erain(it,it0,ilon,jlat,sst)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , it0 , ilon , jlat
      intent (in) it , it0 , ilon , jlat
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
      character(33) :: pathaddname
      logical :: there
      character(3) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer , dimension(10) :: icount , istart
      integer :: inet , istatus
      real(8) :: xadd , xscale
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
      if ( it==it0 ) then
        pathaddname = '../DATA/SST/sstERAIN.1989-2009.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) pathaddname , ' is not available'
          stop
        end if
        istatus = nf90_open('../DATA/SST/sstERAIN.1989-2009.nc',        &
                 &          nf90_nowrite,inet)
        istatus = nf90_get_att(inet,4,'scale_factor',xscale)
        istatus = nf90_get_att(inet,4,'add_offset',xadd)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,4,work,istart,icount)
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          sst(i,jlat+1-j) = work(i,j)*xscale + xadd
        end do
      end do
      istatus = nf90_close(inet)
!
      end subroutine sst_erain
