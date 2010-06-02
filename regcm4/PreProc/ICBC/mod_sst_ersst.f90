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

      module mod_sst_ersst

      contains

      subroutine sst_ersst

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
      use mod_dynparam
      use mod_date
      use mod_interp , only : bilinx
      use mod_printl

      implicit none
!
! PARAMETERS
!
      integer , parameter :: ilon = 240 , jlat = 121
      integer , parameter :: idtbc = 6
!
! Local variables
!
      integer :: i , it , j , mrec , nday , nhour , nmo , nyear
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      integer :: idate , ierastart , ierrec , nsteps
      logical :: there
      real(4) , dimension(ilon,jlat) :: sst
      real(4) , allocatable , dimension(:,:) :: lu , sstmm , xlat , xlon
      character(256) :: terfile , sstfile , inpfile
!
      allocate(lu(iy,jx))
      allocate(sstmm(iy,jx))
      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))
!
      if ( ssttyp=='ERSST' ) then
        there = .false.
        if ( (globidate1>=1989010100 .and. globidate1<=2009053118) .or. &
           & (globidate2>=1989010100 .and. globidate2<=2009053118) ) then
          inquire (file=trim(inpglob)//'/SST/sstERAIN.1989-2009.nc',    &
            &      exist=there)
          if ( .not.there ) then
            print * , 'sstERAIN.1989-2009.nc is not available' ,        &
                 &' under ',trim(inpglob),'/SST/'
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
        if ( (globidate1>=1989010100 .and. globidate1<=2009053118) .or. &
           & (globidate2>=1989010100 .and. globidate2<=2009053118) ) then
          inquire (file=trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc',  &
                 & exist=there)
          if ( .not.there ) then
            print * , 'tskinERAIN.1989-2009.nc is not available' ,      &
                 &' under ',trim(inpglob),'/SST/'
            stop
          end if
        end if
        if ( .not.there ) then
          print * , 'ERSKT Skin Temperature is just available' ,        &
               &' from 1989010100 to 2009053118'
          stop
        end if
      else
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported types are ERSST ERSKT'
        stop
      end if
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,      &
          &  '_SST.RCM'
      open (21,file=sstfile,form='unformatted',status='replace')
 
      nsteps = idatediff(globidate2,globidate1)/idtbc + 1

      write (*,*) 'GLOBIDATE1 : ' , globidate1
      write (*,*) 'GLOBIDATE2 : ' , globidate2
      write (*,*) 'NSTEPS     : ' , nsteps

!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      write (terfile,99001)                                             &
        & trim(dirter), pthsep, trim(domname) , '.INFO'
      open (10,file=terfile,form='unformatted',recl=iy*jx*ibyte,        &
         &  access='direct',status='unknown',err=100)
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,      &
           & '_RCM_SST.dat'
      open (25,file=sstfile,status='unknown',form='unformatted',        &
          & recl=iy*jx*ibyte,access='direct')
      if ( igrads==1 ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname) ,    &
          &  '_RCM_SST.ctl'
        open (31,file=sstfile,status='replace')
        write (31,'(a,a,a)') 'dset ^',trim(domname),'_RCM_SST.dat'
      end if
      call gridmle(xlon,xlat,lu,iy,jx,globidate1,nsteps)
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.5
      end do
      do j = 1 , jlat
        lati(j) = -90. + 1.5*float(j-1)
      end do
 
      idate = globidate1
      ierastart = 1989010100
      do it = 1 , nsteps

        ierrec = idatediff(idate,ierastart)/idtbc+1

        if ( ssttyp=='ERSST' ) then
          inpfile = trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
          call sst_erain(ierrec,ilon,jlat,sst,inpfile)
        else if ( ssttyp=='ERSKT' ) then
          inpfile = trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
          call skt_erain(ierrec,ilon,jlat,sst,inpfile)
        else
        end if

        call split_idate(idate, nyear, nmo, nday, nhour)
 
!       ******           PRINT OUT DATA AS A CHECK
        if ( nmo==1 ) call printl(sst,ilon,jlat)
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        write (21) nhour , nday , nmo , nyear , sstmm
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
        mrec = mrec + 1
        write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)

        call addhours(idate, idtbc)

      end do
 
      deallocate(lu)
      deallocate(sstmm)
      deallocate(xlat)
      deallocate(xlon)

      return

!     4810 PRINT *,'ERROR OPENING GISST FILE'
!     STOP '4810 IN PROGRAM RDSST'
 100  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'

99001 format (a,a,a,a)
      end subroutine sst_ersst
!
!-----------------------------------------------------------------------
!
      subroutine gridmle(xlon,xlat,lu,iy,jx,idate1,numrec)
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
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
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
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
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
      end subroutine gridmle
!
!-----------------------------------------------------------------------
!
      subroutine skt_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
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
      data varname/'skt'/
!
      if ( it==1 ) then
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
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          sst(i,jlat+1-j) = work(i,j)*xscale + xadd
        end do
      end do
!
      end subroutine skt_erain
!
!-----------------------------------------------------------------------
!
      subroutine sst_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
      integer , dimension(10) , save :: icount , istart
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
      if ( it==1 ) then
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
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          sst(i,jlat+1-j) = work(i,j)*xscale + xadd
        end do
      end do
!
      end subroutine sst_erain
!
      end module mod_sst_ersst
