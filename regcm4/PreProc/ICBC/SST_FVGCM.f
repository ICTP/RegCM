      PROGRAM RDSST
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Comments on dataset sources and location:                          c
c                                                                    c
c FVGCM    HadAMH_SST in the original netCDF format.                 c
c          for 'RF'          run, from 1959 to 1991, 385 months      c
c          for 'A2' and 'B2' run, from 2069 to 2101, 385 months      c
c                                                                    c
c          ML= 1 is   0.0; ML= 2 is   1.875; => ML=192 is 358.125E   c
c          NL= 1 is  90.0; ML= 2 is  88.75 ; => ML=145 is -90.       c
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'icbc.param'
      INTEGER MLON
      INTEGER NLAT
      PARAMETER ( MLON=192, NLAT=145 )
 
      REAL    SST2(MLON,NLAT), LONI(MLON), LATI(NLAT)
      COMMON /CDCvars/ SST2
      INTEGER IDATE,IDATE0
      COMMON /DATE/IDATE,IDATE0
      REAL    SSTMM(IY,JX)
      REAL    XLON(IY,JX),XLAT(IY,JX), lu(IY,JX)
      COMMON /SAVECM/ SSTMM,XLON,XLAT,lu
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      INTEGER NDAY,NMO,NYEAR,I,J
      INTEGER MREC,idateo,idatef
      logical there
      integer ludom, lund(20), lumax, k
      real(kind=4)  temp
      COMMON /ALLSST/ temp(192,145)
      integer it
      integer isystm,system
!      external system

      IF(SSTTYP.EQ.'FV_RF') THEN
         inquire(file='../DATA/SST/Sst_1959_1991ref.dat',exist=there)
         if(.not.there) then
            print *, 'Sst_1959_1991ref.dat is not available'
     &             , ' under ../DATA/SST/'
         endif
         open(11,file='../DATA/SST/Sst_1959_1991ref.dat'
     &          ,form='unformatted'
     &          ,recl=192*145*ibyte,access='direct')
      ELSE IF(SSTTYP.EQ.'FV_A2') THEN
         inquire(file='../DATA/SST/Sst_2069_2101_A2.dat',exist=there)
         if(.not.there) then
            print *, 'Sst_2069_2101_A2.dat is not available'
     &             , ' under ../DATA/SST/'
         endif
         open(11,file='../DATA/SST/Sst_2069_2101_A2.dat'
     &          ,form='unformatted'
     &          ,recl=192*145*ibyte,access='direct')
      ELSE IF(SSTTYP.EQ.'FV_B2') THEN
         inquire(file='../DATA/SST/Sst_2069_2101_B2.dat',exist=there)
         if(.not.there) then
            print *, 'Sst_2069_2101_B2.dat is not available'
     &             , ' under ../DATA/SST/'
         endif
         open(11,file='../DATA/SST/Sst_2069_2101_B2.dat'
     &          ,form='unformatted'
     &          ,recl=192*145*ibyte,access='direct')
      ELSE
         WRITE(*,*) 'PLEASE SET SSTTYP in domain.param'
         STOP
      ENDIF
      inquire(file='SST.RCM',exist=there)
      if(there) isystm=system('/bin/rm SST.RCM')
      OPEN(21,file='SST.RCM',form='unformatted')
 
C ******    ON WHAT RegCM GRID ARE SST DESIRED?
      OPEN(10,file='../../Input/DOMAIN.INFO',form='unformatted'
     &       ,recl=IY*JX*ibyte,access='direct'
     &          ,status='unknown',ERR=4830)
      idate=IDATE1/10000
      if(idate-(idate/100)*100.eq.1) then
         idate=idate-89
      else
         idate=idate-1
      endif
      idateo=idate
      IDATE0=idateo*10000+100
      idate=IDATE2/10000
      if(idate-(idate/100)*100.eq.12) then
         idate=idate+89
      else
         idate=idate+1
      endif
      idatef=idate
      print*,idate1,idate2,idateo,idatef
      CALL GRIDML(XLON,XLAT,lu,IY,JX,IDATEO,IDATEF
     &           ,TRUELATL,TRUELATH)
      OPEN(25,file='RCM_SST.dat',status='unknown',form='unformatted'
     &   ,recl=IY*JX*ibyte,access='direct')
      MREC=0

C ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      DO I=1,MLON
         LONI(I) = FLOAT(I-1)*1.875
      ENDDO
      DO J=1,NLAT
         LATI(J) = -90. + 1.25 * FLOAT(J-1)
      ENDDO
 
C **  REF  SST DATA, 1.875x1.1.25, AVAILABLE FROM 16/1/1959 TO 16/1/1991
C ** A2&B2 SST DATA, 1.875x1.1.25, AVAILABLE FROM 16/1/2069 TO 16/1/2101
      IDATE=IDATEO
      DO while (IDATE.le.IDATEF)
         NYEAR=IDATE/100
         NMO=IDATE-NYEAR*100
         NDAY=16
         WRITE(*,*) IDATE*10000+100,IDATE0
c        IF(IDATE*10000+100.EQ.IDATE0) CALL SST_MN(SSTTYP)
         IF(SSTTYP.EQ.'FV_RF') THEN
            it=(NYEAR-1959)*12+NMO
         ELSE 
            it=(NYEAR-2069)*12+NMO
         ENDIF
         read(11,rec=it) temp
         DO J=1,NLAT
         DO I=1,MLON
c           if(temp(I,NLAT+1-J).gt.-9000.0.and.
c    &         temp(I,NLAT+1-J).lt.10000.0) then
c              SST2(I,J)=temp(I,NLAT+1-J)
            if(temp(I,J).gt.-9000.0.and.
     &         temp(I,J).lt.10000.0) then
               SST2(I,J)=temp(I,J) 
            else
               SST2(I,J)=-9999.
            endif
         ENDDO
         ENDDO

C ******           PRINT OUT DATA AS A CHECK
         IF(NMO.EQ.1) CALL PRINTL ( SST2, MLON, NLAT )
 
         CALL BILINX( SST2, LONI, LATI, MLON, NLAT
     &              , SSTMM, XLON, XLAT,   IY,   JX, 1 )
         PRINT *,'XLON,XLAT,SST=', XLON(1,1), XLAT(1,1), SSTMM(1,1)
 
         DO J=1,JX
         DO I=1,IY
            if (sstmm(I,J).lt.-5000. .and.
     &          (lu(i,j).gt.13.5.and.lu(i,j).lt.15.5)) then
              do k=1,20
                lund(k) = 0.0
              end do
              lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
              lund(nint(lu(i-1,j)))   = lund(nint(lu(i-1,j)))   + 3
              lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
              lund(nint(lu(i,j-1)))   = lund(nint(lu(i,j-1)))   + 3
              lund(nint(lu(i,j+1)))   = lund(nint(lu(i,j+1)))   + 3
              lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
              lund(nint(lu(i+1,j)))   = lund(nint(lu(i+1,j)))   + 3
              lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
              ludom = 18
              lumax = 0
              do k=1,20
                if (k.le.13 .or. k.ge.16)  then
                if (lund(k).gt.lumax) then
                  ludom = k
                  lumax = lund(k)
                end if
                end if
              end do
              lu(i,j) = float(ludom)
              print*,ludom,sstmm(i,j)
            end if
            IF(SSTMM(I,J).GT.-100.) then
               SSTMM(I,J) = SSTMM(I,J)
            ELSE
               SSTMM(I,J) = -9999.
            ENDIF
         ENDDO
         ENDDO
 
C ******           WRITE OUT SST DATA ON MM4 GRID
         WRITE(21) NDAY,NMO,NYEAR,SSTMM
         PRINT *, 'WRITING OUT MM4 SST DATA:', NMO, NYEAR
         IDATE=IDATE+1
         IF(NMO.eq.12) IDATE=IDATE+88
         MREC=MREC+1
         WRITE(25,rec=MREC)((SSTMM(I,J),J=1,JX),I=1,IY)
      ENDDO
      write(10,rec=4) ((lu(I,J),J=1,JX),I=1,IY)
 
      STOP 99999
! 4810 PRINT *,'ERROR OPENING GISST FILE'
!      STOP '4810 IN PROGRAM RDSST'
 4830 PRINT *,'ERROR OPENING DOMAIN HEADER FILE'
      STOP '4830 IN PROGRAM RDSST'
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BILINX(  IN, LONI, LATI, NLONI, NLATI
     A                 , OUT, LONO, LATO,    IY,    JX, NFLDS )
      implicit none
      INTEGER NLONI,NLATI,NFLDS,IY,JX
C
C  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
C  RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF GRID2.
C  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
C  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT..
C  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
C  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
C  AS WEIGHTING FACTORS IN THE INTERPOLATION.
C  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
C  INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
C  THE CROSS POINTS IN THE MM4 MODEL.
C
C  IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
C  OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
C  LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
C  LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
C  P.........EAST-WEST WEIGHTING FACTOR.
C  Q.........NORTH-SOUTH WEIGHTING FACTOR.
C  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
C  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
 
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRIDML(XLON,XLAT,lu,IY,JX,IDATE1,IDATE2
     &                 ,TRUELATL,TRUELATH)
      implicit none
      INTEGER IY,JX,IDATE1,IDATE2,IGRADS,IBIGEND
      REAL    TRUELATL,TRUELATH
      real    XLON(IY,JX),XLAT(IY,JX), lu(IY,JX)
      INTEGER IYY,JXX
      REAL    DSINM,CLAT,CLON,PLAT,PLON
      CHARACTER*6 iproj
      REAL    SIGMAF(30),PTOP
      LOGICAL there
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer month,period
      integer i,j,k
      real    alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real    centerj,centeri,GRDFAC
      integer ny,nx,kz
      integer isystm,system
!      external system
C
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
C
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
         month=IDATE1-(IDATE1/100)*100
         period=(IDATE2/100-IDATE1/100)*12+
     &          (IDATE2-(IDATE2/100)*100)-(IDATE1-(IDATE1/100)*100)+1
         write(31,400) period,cmonth(month),IDATE1/100
 400  format('tdef ',I4,' linear 00z16',A3,I4,' 1mo')
         write(31,500) 1
 500  format('vars ',I1)
 600  format(A4,'0 99 ',A26)
         write(31,600) 'sst ','surface elevation          '
         write(31,'(a)') 'endvars'
         close(31)
      endif
c
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
