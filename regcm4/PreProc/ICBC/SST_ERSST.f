      PROGRAM RDSST
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
      implicit none
      include 'icbc.param'
      INTEGER MLON
      INTEGER NLAT
      PARAMETER (MLON=240,NLAT=121)
 
      REAL    LONI(MLON), LATI(NLAT)
      REAL*4  SST2
      common /ERAvars/ sst2(mlon,nlat)
      INTEGER IDATE,IDATE0
      COMMON /DATE/IDATE,IDATE0
      REAL    SSTMM(IY,JX)
      REAL    XLON(IY,JX),XLAT(IY,JX), lu(IY,JX)
      COMMON /SAVECM/ SSTMM,XLON,XLAT,lu
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      INTEGER MDATE
      COMMON /DATENUM/MDATE(29824)
      INTEGER NSTART,NNNEND
      COMMON /EH5OM/ NSTART,NNNEND
!
      INTEGER NHOUR,NDAY,NMO,NYEAR,I,J
      INTEGER NREC,MREC,idateo,idatef
      logical there
      integer it
      integer isystm,system
      external system

      IF(SSTTYP.EQ.'ERSST') THEN
        there=.false.
        IF((IDATE1.GE.1989010100.AND.IDATE1.LE.2009053118).OR.
     &     (IDATE2.GE.1989010100.AND.IDATE2.LE.2009053118)) THEN
          inquire(file='../DATA/SST/sstERAIN.1989-2009.nc',exist=there)
          if(.not.there) then
            print *, 'sstERAIN.1989-2009.nc is not available'
     &             , ' under ../DATA/SST/'
            STOP
          endif
        ENDIF
        IF(.not.there) THEN
          print *,'ERSST Sea Surface Temperature is just available'
     &           ,' from 1989010100 to 2009053118'
          STOP
        ENDIF
      ELSE IF(SSTTYP.EQ.'ERSKT') THEN
        there=.false.
        IF((IDATE1.GE.1989010100.AND.IDATE1.LE.2009053118).OR.
     &     (IDATE2.GE.1989010100.AND.IDATE2.LE.2009053118)) THEN
        inquire(file='../DATA/SST/tskinERAIN.1989-2009.nc',exist=there)
          if(.not.there) then
            print *, 'tskinERAIN.1989-2009.nc is not available'
     &             , ' under ../DATA/SST/'
            STOP
          endif
        ENDIF
        IF(.not.there) THEN
          print *,'ERSKT Skin Temperature is just available'
     &           ,' from 1989010100 to 2009053118'
          STOP
        ENDIF
      ELSE
        WRITE(*,*) 'PLEASE SET the right SSTTYP in domain.param'
        STOP
      ENDIF
      inquire(file='SST.RCM',exist=there)
      if(there) isystm=system('/bin/rm SST.RCM')
      OPEN(21,file='SST.RCM',form='unformatted')
 
! ******    ON WHAT RegCM GRID ARE SST DESIRED?
      OPEN(10,file='../../Input/DOMAIN.INFO',form='unformatted'
     &       ,recl=IY*JX*ibyte,access='direct'
     &          ,status='unknown',ERR=4830)
      CALL INITDATE
      CALL FINDDATE(NSTART,IDATE1)
      CALL FINDDATE(NNNEND,IDATE2)
      WRITE(*,*) NSTART, NNNEND
      print*,idate1,NNNEND-NSTART+1
      CALL GRIDML(XLON,XLAT,lu,IY,JX,IDATE1,NNNEND-NSTART+1
     &           ,TRUELATL,TRUELATH)
      OPEN(25,file='RCM_SST.dat',status='unknown',form='unformatted'
     &   ,recl=IY*JX*ibyte,access='direct')
      MREC=0

! ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      DO I=1,MLON
         LONI(I) = FLOAT(I-1)*1.5
      ENDDO
      DO J=1,NLAT
        LATI(J) = -90. + 1.5 * FLOAT(J-1)
      ENDDO
 
      DO IT=NSTART,NNNEND
        IDATE=MDATE(IT)
        IF(SSTTYP.eq.'ERSST') THEN
          call SST_ERAIN(IT,NSTART)
        ELSE IF(SSTTYP.eq.'ERSKT') THEN
          call SKT_ERAIN(IT,NSTART)
        ENDIF
        NYEAR=IDATE/1000000
        NMO  =(IDATE-NYEAR*1000000)/10000
        NDAY =(IDATE-NYEAR*1000000-NMO*10000)/100
        NHOUR= IDATE-NYEAR*1000000-NMO*10000-NDAY*100

! ******           PRINT OUT DATA AS A CHECK
        IF(NMO.EQ.1) CALL PRINTL ( SST2, MLON, NLAT )
 
        CALL BILINX( SST2, LONI, LATI, MLON, NLAT
     &             , SSTMM, XLON, XLAT,   IY,   JX, 1 )
        PRINT *,'XLON,XLAT,SST=', XLON(1,1), XLAT(1,1), SSTMM(1,1)
 
! ******           WRITE OUT SST DATA ON MM4 GRID
        WRITE(21) NHOUR,NDAY,NMO,NYEAR,SSTMM
        PRINT *, 'WRITING OUT MM4 SST DATA:', NMO, NYEAR
        MREC=MREC+1
        WRITE(25,rec=MREC)((SSTMM(I,J),J=1,JX),I=1,IY)
      ENDDO
 
      STOP 99999
 4810 PRINT *,'ERROR OPENING GISST FILE'
      STOP '4810 IN PROGRAM RDSST'
 4830 PRINT *,'ERROR OPENING DOMAIN HEADER FILE'
      STOP '4830 IN PROGRAM RDSST'
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE BILINX(  IN, LONI, LATI, NLONI, NLATI
     A                 , OUT, LONO, LATO, IY, JX, NFLDS )
      implicit none
      INTEGER NLONI,NLATI,NFLDS,IY,JX
!
!  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
!  RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF GRID2.
!  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
!  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT..
!  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
!  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
!  AS WEIGHTING FACTORS IN THE INTERPOLATION.
!  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!  INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
!  THE CROSS POINTS IN THE MM4 MODEL.
!
!  IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!  OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
!  LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!  LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!  P.........EAST-WEST WEIGHTING FACTOR.
!  Q.........NORTH-SOUTH WEIGHTING FACTOR.
!  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
 
      REAL    IN(NLONI,NLATI,NFLDS), LONI(NLONI), LATI(NLATI)
      REAL    OUT(IY,JX,NFLDS), LONO(IY,JX), LATO(IY,JX)
      REAL    LON360
      INTEGER I,J,L,IP,IPP1,JQ,JQP1
      REAL    YIND,Q,XIND,P
      REAL    BAS,SUM
 
      DO 120 J=1,JX
      DO 110 I=1,IY
 
      YIND = (((LATO(I,J)-LATI(1))/(LATI(NLATI)-LATI(1)))
     &     * FLOAT(NLATI-1))+1.
      JQ=INT(YIND)
      JQP1=MIN0(JQ+1,NLATI)
      Q=YIND-JQ
 
      LON360 = LONO(I,J)
      IF(LONO(I,J).LT.0.) LON360 = LONO(I,J) + 360.
      XIND = (((LON360   -LONI(1))/(LONI(NLONI)-LONI(1)))
     &     * FLOAT(NLONI-1))+1.
      IP=INT(XIND)
      IPP1=MIN0(IP+1,NLONI)
      P=XIND-IP
 
      DO 100 L=1,NFLDS
      SUM=0.0
      BAS=0.0
      IF(IN(IP,JQ,L).LT.-9990.0.AND.IN(IPP1,JQ,L).LT.-9990.0.AND.
     &   IN(IPP1,JQP1,L).LT.-9990.0.AND.IN(IP,JQP1,L).LT.-9990.0 ) THEN
         OUT(I,J,L)=-9999.
      ELSE 
         IF(IN(IP,JQ,L).GT.-9990.0) THEN
            SUM=SUM+(1.-Q)*(1.-P)*IN(IP,JQ,L)
            BAS=BAS+(1.-Q)*(1.-P)
         ENDIF
         IF(IN(IPP1,JQ,L).GT.-9990.0) THEN
            SUM=SUM+(1.-Q)*P*IN(IPP1,JQ,L)
            BAS=BAS+(1.-Q)*P
         ENDIF
         IF(IN(IPP1,JQP1,L).GT.-9990.0) THEN
            SUM=SUM+Q*P*IN(IPP1,JQP1,L)
            BAS=BAS+Q*P
         ENDIF
         IF(IN(IP,JQP1,L).GT.-9990.0 ) THEN
            SUM=SUM+Q*(1.-P)*IN(IP,JQP1,L)
            BAS=BAS+Q*(1.-P)
         ENDIF
         OUT(I,J,L) = SUM/BAS
      ENDIF
  100 CONTINUE
  110 CONTINUE
 
  120 CONTINUE
 
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GRIDML(XLON,XLAT,lu,IY,JX,IDATE1,NUMREC
     &                 ,TRUELATL,TRUELATH)
      implicit none
      INTEGER IY,JX,IDATE1,NUMREC,IGRADS,IBIGEND
      REAL    TRUELATL,TRUELATH
      real    XLON(IY,JX),XLAT(IY,JX), lu(IY,JX)
      INTEGER IYY,JXX
      REAL    DSINM,CLAT,CLON,PLAT,PLON
      CHARACTER*6 iproj
      REAL    SIGMAF(30),PTOP
      LOGICAL there
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer nyear,month,nday,nhour,period
      integer i,j,k
      real    alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real    centerj,centeri,GRDFAC
      integer ny,nx,kz
      integer isystm,system
      external system
!
      READ(10,rec=1) IYY,JXX,kz,DSINM,CLAT,CLON,PLAT,PLON,GRDFAC,iproj
     &              ,(SIGMAF(K),K=1,KZ+1),PTOP,igrads,ibigend
     &              ,TRUELATL,TRUELATH
      IF(IYY.NE.IY.OR.JXX.NE.JX) THEN
         WRITE(*,*)'IY,JX,IYY,JXX',IY,JX,IYY,JXX
         STOP
      ENDIF
      READ(10,rec=4) ((lu(I,J),J=1,JX),I=1,IY)
      READ(10,rec=5) ((XLAT(I,J),J=1,JX),I=1,IY)
      READ(10,rec=6) ((XLON(I,J),J=1,JX),I=1,IY)
!
      IF(IGRADS.EQ.1) THEN
         inquire(file='RCM_SST.ctl',exist=there)
         if(there) isystm=system('/bin/rm RCM_SST.ctl')
         OPEN(31,file='RCM_SST.ctl',status='new')
         write(31,'(a)') 'dset ^RCM_SST.dat'
         write(31,'(a)') 'title SST fields for RegCM domain'
         if(ibigend.eq.1) then
            write(31,'(a)') 'options big_endian'
         else
            write(31,'(a)') 'options little_endian'
         endif
         write(31,'(a)') 'undef -9999.'
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx
               if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
               if(xlat(iy,j).gt.alatmax)alatmax=xlat(iy,j)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy
            do j=1,jx
               if(clon.ge.0.0) then
                  if(xlon(i,j).ge.0.0) then
                     alonmin = amin1(alonmin,xlon(i,j))
                     alonmax = amax1(alonmax,xlon(i,j))
                  else
                     if(abs(clon-xlon(i,j)).lt.
     &                  abs(clon-(xlon(i,j)+360.))) then
                        alonmin = amin1(alonmin,xlon(i,j))
                        alonmax = amax1(alonmax,xlon(i,j))
                     else
                        alonmin = amin1(alonmin,xlon(i,j)+360.)
                        alonmax = amax1(alonmax,xlon(i,j)+360.)
                     endif
                  endif
               else
                  if(xlon(i,j).lt.0.0) then
                     alonmin = amin1(alonmin,xlon(i,j))
                     alonmax = amax1(alonmax,xlon(i,j))
                  else
                     if(abs(clon-xlon(i,j)).lt.
     &                  abs(clon-(xlon(i,j)-360.))) then
                        alonmin = amin1(alonmin,xlon(i,j))
                        alonmax = amax1(alonmax,xlon(i,j))
                     else
                        alonmin = amin1(alonmin,xlon(i,j)-360.)
                        alonmax = amax1(alonmax,xlon(i,j)-360.)
                     endif
                  endif
               endif
            enddo
            enddo
            rlatinc=DSINM*0.001/111./2.
            rloninc=DSINM*0.001/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))

            centerj=jx/2.
            centeri=iy/2.
         endif
         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx,iy,clat,clon,centerj,centeri,
     &                    truelatL,truelatH,clon,DSINM,DSINM
 100  format('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx,xlon(1,1),xlon(1,2)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(i,1),i=1,iy)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
            write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &                ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection'
     &                ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            write(31,230) jx,iy,plon,plat,DSINM/111000.
     &                                   ,DSINM/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else            
            write(*,*) 'Are you sure your map projection is correct ?'
            stop
         endif
         write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
         period=NUMREC
         nyear=IDATE1/1000000
         month=(IDATE1-nyear*1000000)/10000
         nday =(IDATE1-nyear*1000000-month*10000)/100
         nhour=mod(IDATE1,100)
         write(31,400) period,nhour,cday(nday),cmonth(month),nyear
 400  format('tdef ',I6,' linear ',I2,'z',A2,A3,I4,' 6hr')
         write(31,500) 1
 500  format('vars ',I1)
 600  format(A4,'0 99 ',A26)
         write(31,600) 'sst ','surface elevation          '
         write(31,'(a)') 'endvars'
         close(31)
      endif
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINTL (A, N1, N2)
      implicit none
      INTEGER N1,N2
      REAL    A(N1,N2)
      INTEGER INC2,INC1,L,K
 
      INC2 = -1
      INC1 =  1
 
      DO 1 L=N2,1,INC2
    1 PRINT 11,L,(A(K,L),K=1,17,INC1)
   11 FORMAT(1X,I3,2X,17F7.2)
 
      RETURN
      END

      SUBROUTINE FEXIST(filnam)
      implicit none

      character filnam*(*), yesno*1
      logical there
      
 1    inquire(file=filnam,exist=there)
      if (there) then
 2      print*,' '
        print*,' '
        print*,'**************************************************'
        print*,'FILE ALREADY EXISTS:  ',filnam
        print*,'Do you want to overwrite the existing file? [y/n/q]'
        read(*,*) yesno
        if (yesno.eq.'y') then
          return
        else if (yesno.eq.'n') then
          print*,'ENTER NEW FILE NAME'
          read(*,*) filnam
          go to 1
        else if (yesno.eq.'q') then
          stop 999
        else
          go to 2
        end if
      end if

      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE INITDATE
      IMPLICIT NONE
!
      INTEGER MDATE
      COMMON /DATENUM/MDATE(29824)
      INTEGER MBASE,NBASE,NREC,NYEAR,MON,NDAY,I,M
!
      NREC=0
      DO NYEAR=1989,2009
        MBASE = NYEAR*1000000
        DO MON=1,12
          MBASE = MBASE+10000
          IF(MON.EQ.1.OR.MON.EQ.3.OR.MON.EQ.5.OR.MON.EQ.7
     &              .OR.MON.EQ.8.OR.MON.EQ.10.OR.MON.EQ.12) THEN
            NDAY=31
          ELSE IF(MON.EQ.4.OR.MON.EQ.6.OR.MON.EQ.9.OR.MON.EQ.11)THEN
            NDAY=30
          ELSE
            IF(MOD(NYEAR,4).EQ.0) THEN
              NDAY=29
              IF(MOD(NYEAR,100).EQ.0) THEN
                NDAY=28
                IF(MOD(NYEAR,400).EQ.0) THEN
                  NDAY=29
                ENDIF
              ENDIF
            ELSE
              NDAY=28
            ENDIF
          ENDIF
          NBASE = MBASE
          DO I=1,NDAY
            NBASE = NBASE+100
            DO M=1,4
              NREC=NREC+1
              if (NREC.gt.29824) goto 333
              IF(M.EQ.1) THEN
                MDATE(NREC)=NBASE
              ELSE IF(M.EQ.2) THEN
                MDATE(NREC)=NBASE+6
              ELSE IF(M.EQ.3) THEN
                MDATE(NREC)=NBASE+12
              ELSE
                MDATE(NREC)=NBASE+18
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
 333  continue
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE FINDDATE(NPOS,IDATE)
      IMPLICIT NONE
      INTEGER MDATE
      COMMON /DATENUM/MDATE(29824)
!
      INTEGER NPOS,IDATE
      INTEGER I
!
      I=0
 10   CONTINUE
      I=I+1
      IF(MDATE(I).EQ.IDATE) THEN
         NPOS=I
         GO TO 200
      ENDIF
      IF(I.GT.29824) GOTO 100
      GO TO 10
 100  WRITE(*,*) 'ERROR IN FINDDATE'
      STOP
 200  RETURN
      END
      subroutine SST_ERAIN(IT,IT0)
      implicit none
      INTEGER IT,IT0
      include '../../Commons/env/include/netcdf.inc'
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat
      parameter (ilon=240)
      parameter (jlat=121)
      integer inet,start,count
      integer status
      COMMON /ECOPEN/ inet,start(10),count(10),status
      real*8  xscale,xadd
      COMMON /EPATCH/ xscale,xadd

!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers. The array 'sst' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer*2 work(ilon,jlat)
!bxq
      real*4  sst
      common /ERAvars/ sst(ilon,jlat)
!bxq_
      character*33 pathaddname
      character*3  varname
      data varname/'sst'/
      logical there
      integer i,j,n
!
      IF(IT.EQ.IT0) then
         pathaddname='../DATA/SST/sstERAIN.1989-2009.nc'
         inquire(file=pathaddname,exist=there)
         if(.not.there) then
            write(*,*) pathaddname,' is not available'
            stop
         endif
         status=nf_open('../DATA/SST/sstERAIN.1989-2009.nc',nf_nowrite
     &                                                     ,inet)
         status=nf_get_att_double(inet,4,'scale_factor',xscale)
         status=nf_get_att_double(inet,4,'add_offset',xadd)
         start(1) = 1
         start(2) = 1
         count(1) = 240
         count(2) = 121
         do n=4,10
            start(n) = 0
            count(n) = 0
         enddo
      ENDIF

!bxq
      start(3) = it
      count(3) = 1
      status=nf_get_vara_int2(inet,4,start,count,work)
!bxq_
!
      do j=1,jlat
      do i=1,ilon
         sst(i,jlat+1-j) = work(i,j)*xscale+xadd
      enddo
      enddo
! 
      return
      end
      subroutine SKT_ERAIN(IT,IT0)
      implicit none
      INTEGER IT,IT0
      include '../../Commons/env/include/netcdf.inc'
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat
      parameter (ilon=240)
      parameter (jlat=121)
      integer inet,start,count
      integer status
      COMMON /ECOPEN/ inet,start(10),count(10),status
      real*8  xscale,xadd
      COMMON /EPATCH/ xscale,xadd

!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers. The array 'sst' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer*2 work(ilon,jlat)
!bxq
      real*4  sst
      common /ERAvars/ sst(ilon,jlat)
!bxq_
      character*35 pathaddname
      character*3  varname
      data varname/'skt'/
      logical there
      integer i,j,n
!
      IF(IT.EQ.IT0) then
         pathaddname='../DATA/SST/tskinERAIN.1989-2009.nc'
         inquire(file=pathaddname,exist=there)
         if(.not.there) then
            write(*,*) pathaddname,' is not available'
            stop
         endif
         status=nf_open('../DATA/SST/tskinERAIN.1989-2009.nc',nf_nowrite
     &                                                       ,inet)
         status=nf_get_att_double(inet,4,'scale_factor',xscale)
         status=nf_get_att_double(inet,4,'add_offset',xadd)
         start(1) = 1
         start(2) = 1
         count(1) = 240
         count(2) = 121
         do n=4,10
            start(n) = 0
            count(n) = 0
         enddo
      ENDIF

!bxq
      start(3) = it
      count(3) = 1
      status=nf_get_vara_int2(inet,4,start,count,work)
!bxq_
!
      do j=1,jlat
      do i=1,ilon
         sst(i,jlat+1-j) = work(i,j)*xscale+xadd
      enddo
      enddo
! 
      return
      end
