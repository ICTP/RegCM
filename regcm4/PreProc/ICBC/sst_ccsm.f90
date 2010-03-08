      program rdsst
!
!*******************************************************************************
!
! This is a package of subroutines to read CCSM SST 1x1 degree data in
! NETCDF format and interpolate SST at RCM grid
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!
!******************************************************************************
!******************************************************************************
!       DATA PREPARATION
!       Dataset required to use this code can be preapred using NCO utilitiles
!       such as NCKS, NCRCAT etc.
!       We need top level of TEMP for SSTs which can be extarcted as following:
!       ncks -v time,TEMP -d z_t,0 input.nc output.nc
!       Files can be further concatenated using 'ncrcat'
!       Finally, the POP grid can be converted into lat/lon grid at 1x1 degree
!       resolution  using PopLatLon function in NCL
!******************************************************************************
!     NAMING CONVENTION (Global Data File)
!       ccsm_mn.sst.nc  for TEMP
!     PATH /DATA/SST/
!
!******************************************************************************
      use mod_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      real , dimension(jlat) :: glat
      real , dimension(ilon) :: glon
      real , dimension(ilon,jlat) :: sst
      integer :: idate , idate0
      integer :: i , idatef , idateo , j , k , ludom , lumax , mrec ,   &
               & nday , nmo , nyear
      real , dimension(iy,jx) :: lu , sstmm , xlat , xlon
      integer , dimension(20) :: lund
!
      do i = 1 , ilon
        glon(i) = 0.5 + float(i-1)
      end do
      do j = 1 , jlat
        glat(j) = -89.5 + 1.*float(j-1)
      end do
 
      open (21,file='SST.RCM',form='unformatted',status='replace')
 
!     ******    ON WHAT RegCM GRID ARE SST DESIRED?
      open (10,file='../../Input/DOMAIN.INFO',form='unformatted',       &
          & recl=iy*jx*ibyte,access='direct',status='unknown',err=100)
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
 
      call gridml(xlon,xlat,lu,iy,jx,idateo,idatef,truelatl,truelath)
      open (25,file='RCM_SST.dat',status='unknown',form='unformatted',  &
          & recl=iy*jx*ibyte,access='direct')
      mrec = 0
 
      idate = idateo
      do while ( idate<=idatef )
        nyear = idate/100
        nmo = idate - nyear*100
 
        call ccsm_sst(idate*10000+100,idate0,ilon,jlat,sst)
 
!       ******           PRINT OUT DATA AS A CHECK
        if ( nmo==1 ) call printl(sst,jlat,ilon)
        call bilinx(sst,glon,glat,ilon,jlat,sstmm,xlon,xlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1) &
            & + 273.
 
        do j = 1 , jx
          do i = 1 , iy
            if ( sstmm(i,j)<-5000 .and.                                 &
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
            end if
            if ( sstmm(i,j)>-100. ) then
              sstmm(i,j) = sstmm(i,j) + 273.15
            else
              sstmm(i,j) = -9999.
            end if
          end do
        end do
        write (21) nday , nmo , nyear , sstmm
        print * , 'WRITING OUT CCSM SST DATA:' , nmo , nyear
        idate = idate + 1
        if ( nmo==12 ) idate = idate + 88
        mrec = mrec + 1
        write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
!       print*, sstmm
      end do
      write (10,rec=4) ((lu(i,j),j=1,jx),i=1,iy)
 
      stop 99999
 100  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 in PROGRAM RDSST'
      end program rdsst
!
!-----------------------------------------------------------------------
!
      subroutine gridml(xlon,xlat,lu,iy,jx,idate1,idate2,truelatl,      &
                      & truelath)
      implicit none
!
! Dummy arguments
!
      integer :: idate1 , idate2 , iy , jx
      real :: truelath , truelatl
      real , dimension(iy,jx) :: lu , xlat , xlon
      intent (in) idate1 , idate2 , iy , jx , truelath , truelatl
      intent (out) lu
      intent (inout) xlat , xlon
!
! Local variables
!
      real :: alatmax , alatmin , alonmax , alonmin , centeri ,         &
            & centerj , clat , clon , dsinm , grdfac , plat , plon ,    &
            & ptop , rlatinc , rloninc
      character(3) , dimension(12) :: cmonth
      integer :: i , ibigend , igrads , iyy , j , jxx , k , kz ,        &
               & month , nx , ny , period
      character(6) :: iproj
      real , dimension(30) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      read (10,rec=1) iyy , jxx , kz , dsinm , clat , clon , plat ,     &
                    & plon , grdfac , iproj , (sigmaf(k),k=1,kz+1) ,    &
                    & ptop , igrads , ibigend
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
        month = idate1 - (idate1/100)*100
        period = (idate2/100-idate1/100)*12 + (idate2-(idate2/100)*100) &
               & - (idate1-(idate1/100)*100) + 1
        write (31,99009) period , cmonth(month) , idate1/100
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
99009 format ('tdef ',i4,' linear 00z16',a3,i4,' 1mo')
99010 format ('vars ',i1)
99011 format (a4,'0 99 ',a26)
!
      end subroutine gridml
!
!-----------------------------------------------------------------------
!
!     Subroutine to read required records from SST data file
!
      subroutine ccsm_sst(idate,idate0,ilon,jlat,sst)

      use netcdf

      implicit none

      integer , intent (in) :: idate , idate0
      integer , intent (in) :: ilon , jlat
      real, dimension(ilon, jlat) , intent(out) :: sst

      integer, dimension(12) :: ndays
      character(len=4), dimension(2) :: varname
      integer, allocatable ::  work1(:)
      real, dimension (ilon , jlat) :: work2
      real :: imisng

      character(len=26) :: pathaddname
      integer :: nyear , month
      integer :: inet1
      integer, dimension(10) :: istart , icount , istartt , icountt
      integer, dimension(2) :: ivar2
      integer :: it , icode , i , j , npos , nrec
      integer :: latid , lonid , timid
      integer :: latlen , lonlen , timlen
      logical :: there
      integer :: istatus

      data ndays/31,59,90,120,151,181,212,243,273,304,334,365/
      data varname/'time','TEMP'/
      
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      
      if (idate == idate0) then
         pathaddname = '../DATA/SST/ccsm_mn.sst.nc'
         inquire(file=pathaddname,exist=there)
         if (.not.there) then
            print *, pathaddname,' is not available'
            stop
         endif
         istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
         
         write(*,*) inet1 , pathaddname , icode
      endif  
!     GET DIMENSION IDs
      istatus = nf90_inq_dimid(inet1,'lat',latid)
      istatus = nf90_inq_dimid(inet1,'lon',lonid)
      istatus = nf90_inq_dimid(inet1,'time',timid)

!     GET DIMENSION LENGTHS
      istatus = nf90_inquire_dimension(inet1,latid,len=latlen)
      istatus = nf90_inquire_dimension(inet1,lonid,len=lonlen)
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      allocate(work1(timlen))
      
!     MAKE SURE THAT SST DATA IS AT 1X1 DEGREE
      if(latlen /= jlat .or. lonlen /= ilon) then
         print*,'DIMENSIONS DO NOT MATCH'
         print*,'No. of LON in SST file =',lonlen
         print*,'No. of LON in 1x1 degree gloabl grid =',ilon
         print*,'No. of LAT in SST file =',latlen
         print*,'No. of LON in 1x1 degree gloabl grid =',jlat
         STOP   'Check SST data file' 
      endif
!     GET VARIABLE IDs
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
!     GET MISSING DATA VALUE
      istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',imisng)
!     GET TIME VALUES
      istartt(1) = 1
      icountt(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istartt,icountt)
      
!     CHECK FOR THE REQUIRED RECORD IN DATA FILE  
      npos = (nyear - 1000) * 365 + ndays(month)
      i = 0
      print *,  npos
 10   continue
      i = i + 1
      if (npos < work1(i) .or. npos > work1(timlen)) then
         print *, 'Error in finding SST data for',(idate-100)/10000
         print *, 'Required NREC=',npos
         stop    'Check SST data file' 
      else if (work1(i) == npos) then
         nrec=i
         go to 20
      end if
      go to 10
 
 20   it = nrec
      icount(1) = ilon
      icount(2) = jlat
      icount(3) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      do j = 1 , jlat
         do i = 1 , ilon
            if (work2(i,j) > (imisng+10) .and. work2(i,j) < 10000.0) then
               sst(i,j) = work2(i,j)
            else
               sst(i,j) = -9999.
            end if
         end do
      end do

      end subroutine ccsm_sst
