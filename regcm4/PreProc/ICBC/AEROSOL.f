      PROGRAM AEROSOL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Comments on dataset sources and location:                          c
c                                                                    c
c EDGAR                                                              c
c                                                                    c
c LIOUSSE96                                                          c
c                                                                    c
c BOND                                                               c
c                                                                    c
c GEIA                                                               c
c                                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'icbc.param'
      INTEGER MLON
      INTEGER NLAT
      PARAMETER ( MLON=360, NLAT=180 )
 
      REAL    AER2(MLON,NLAT), LONI(MLON), LATI(NLAT)
      COMMON /AERvars/ AER2
      REAL    AERMM(IY,JX)
      REAL    XLON(IY,JX),XLAT(IY,JX)
      COMMON /SAVECM/ AERMM,XLON,XLAT
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      INTEGER I,J
      INTEGER NREC
      logical there
      integer isystm,system
!      external system

      inquire(file='../DATA/AERGLOB/AEROSOL.dat',exist=there)
      if(.not.there) then
         print *, 'AEROSOL.dat is not available'
     &          , ' under ../DATA/AERGLOB/'
      endif
      open(11,file='../DATA/AERGLOB/AEROSOL.dat'
     &       ,form='unformatted',recl=360*180*ibyte
     &       ,access='direct',status='old',ERR=4810)
      inquire(file='../../Input/AERO.dat',exist=there)
      if(there) isystm=system('/bin/rm ../../Input/AERO.dat')
      OPEN(25,file='../../Input/AERO.dat',form='unformatted'
     &   ,recl=IY*JX*ibyte,access='direct')
 
C ******    ON WHAT RegCM GRID ARE AEROSOL DESIRED?
      OPEN(10,file='../../Input/DOMAIN.INFO',form='unformatted'
     &       ,recl=IY*JX*ibyte,access='direct'
     &          ,status='unknown',ERR=4830)

C
      CALL GRIDML(XLON,XLAT,IY,JX,ibyte,truelatL,truelatH)
C

C ******    SET UP LONGITUDES AND LATITUDES FOR AEROSOL DATA
      DO I=1,MLON
         LONI(I) = -179.5 + FLOAT(I-1)
      ENDDO
      DO J=1,NLAT
         LATI(J) = -89.5 + 1. * FLOAT(J-1)
      ENDDO
 
C ****** ALL AEROSOL DATA, 1 Deg data, Climate value
      DO NREC=1,39
         READ(11,rec=NREC) AER2
 
         CALL BILINX( AER2, LONI, LATI, MLON, NLAT
     &              , AERMM, XLON, XLAT,   IY,   JX, 1 )
 
C ******           WRITE OUT AEROSOL DATA ON RegCM GRID
         WRITE(25,rec=NREC)((AERMM(I,J),J=1,JX),I=1,IY)
      ENDDO
 
      STOP 99999
 4810 PRINT *,'ERROR OPENING AEROSOL FILE'
      STOP '4810 IN PROGRAM AEROSOL'
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
      SUBROUTINE GRIDML(XLON,XLAT,IY,JX,ibyte,TRUELATL,TRUELATH)
      implicit none
      INTEGER IY,JX,IGRADS,IBIGEND,ibyte,ierr
      REAL    TRUELATL,TRUELATH
      real    XLON(IY,JX),XLAT(IY,JX)
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
      READ(10,rec=1,iostat=ierr) IYY,JXX,kz,DSINM,CLAT,CLON,PLAT,PLON
     &                          ,GRDFAC,iproj,(SIGMAF(K),K=1,KZ+1),PTOP
     &                          ,igrads,ibigend,TRUELATL,TRUELATH
      IF(IYY.NE.IY.OR.JXX.NE.JX) THEN
        print*,'IMPROPER DIMENSION SPECIFICATION (AEROSOL.f)'
        print*,'  icbc.param: ',IY,JX
        print*,'  DOMAIN.INFO: ',IYY,JXX
        print*,'  Also check ibyte in icbc.param: ibyte= ',ibyte
        STOP 'Dimensions (subroutine gridml)'
      ENDIF
      READ(10,rec=5,iostat=ierr) ((XLAT(I,J),J=1,JX),I=1,IY)
      READ(10,rec=6,iostat=ierr) ((XLON(I,J),J=1,JX),I=1,IY)
      if (ierr.ne.0) then
        print*,'END OF FILE REACHED (AEROSOL.f)'
        print*,'  Check ibyte in icbc.param: ibyte= ',ibyte
        stop 'EOF (subroutine gridml)'
      endif
C
      IF(IGRADS.EQ.1) THEN
         inquire(file='../../Input/AERO.ctl',exist=there)
         if(there) isystm=system('/bin/rm ../../Input/AERO.ctl')
         OPEN(31,file='../../Input/AERO.ctl',status='new')
         write(31,'(a)') 'dset ^AERO.dat'
      write(31,'(a)') 'title AEROSOL fields for RegCM domain, kg/m2/sec'
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
         month= 1
         period=1
         write(31,400) period,cmonth(month),2001
 400  format('tdef ',I4,' linear 00z16',A3,I4,' 1mo')
         write(31,500) 39
 500  format('vars ',I2)
 600  format(A6,'0 99 ',A40)
      write(31,600) 'so2   ','Anthropogenic SO2 emission, EDGAR       '
      write(31,600) 'bc    ','Anthropogenic Black Carbon (BC), EDGAR  '
      write(31,600) 'oc    ','Anthropogenic Organic Carbon (OC), EDGAR'
      write(31,600) 'so201 ','Biomass SO2 emission, EDGAR, January    '
      write(31,600) 'so202 ','Biomass SO2 emission, EDGAR, February   '
      write(31,600) 'so203 ','Biomass SO2 emission, EDGAR, March      '
      write(31,600) 'so204 ','Biomass SO2 emission, EDGAR, April      '
      write(31,600) 'so205 ','Biomass SO2 emission, EDGAR, May        '
      write(31,600) 'so206 ','Biomass SO2 emission, EDGAR, June       '
      write(31,600) 'so207 ','Biomass SO2 emission, EDGAR, July       '
      write(31,600) 'so208 ','Biomass SO2 emission, EDGAR, August     '
      write(31,600) 'so209 ','Biomass SO2 emission, EDGAR, September  '
      write(31,600) 'so210 ','Biomass SO2 emission, EDGAR, October    '
      write(31,600) 'so211 ','Biomass SO2 emission, EDGAR, November   '
      write(31,600) 'so212 ','Biomass SO2 emission, EDGAR, December   '
      write(31,600) 'bc_01 ','Biomass BC emission, LIOUSSE, January   '
      write(31,600) 'bc_02 ','Biomass BC emission, LIOUSSE, February  '
      write(31,600) 'bc_03 ','Biomass BC emission, LIOUSSE, March     '
      write(31,600) 'bc_04 ','Biomass BC emission, LIOUSSE, April     '
      write(31,600) 'bc_05 ','Biomass BC emission, LIOUSSE, May       '
      write(31,600) 'bc_06 ','Biomass BC emission, LIOUSSE, June      '
      write(31,600) 'bc_07 ','Biomass BC emission, LIOUSSE, July      '
      write(31,600) 'bc_08 ','Biomass BC emission, LIOUSSE, August    '
      write(31,600) 'bc_09 ','Biomass BC emission, LIOUSSE, September '
      write(31,600) 'bc_10 ','Biomass BC emission, LIOUSSE, October   '
      write(31,600) 'bc_11 ','Biomass BC emission, LIOUSSE, November  '
      write(31,600) 'bc_12 ','Biomass BC emission, LIOUSSE, December  '
      write(31,600) 'oc_01 ','Biomass OC emission, LIOUSSE, January   '
      write(31,600) 'oc_02 ','Biomass OC emission, LIOUSSE, February  '
      write(31,600) 'oc_03 ','Biomass OC emission, LIOUSSE, March     '
      write(31,600) 'oc_04 ','Biomass OC emission, LIOUSSE, April     '
      write(31,600) 'oc_05 ','Biomass OC emission, LIOUSSE, May       '
      write(31,600) 'oc_06 ','Biomass OC emission, LIOUSSE, June      '
      write(31,600) 'oc_07 ','Biomass OC emission, LIOUSSE, July      '
      write(31,600) 'oc_08 ','Biomass OC emission, LIOUSSE, August    '
      write(31,600) 'oc_09 ','Biomass OC emission, LIOUSSE, September '
      write(31,600) 'oc_10 ','Biomass OC emission, LIOUSSE, October   '
      write(31,600) 'oc_11 ','Biomass OC emission, LIOUSSE, November  '
      write(31,600) 'oc_12 ','Biomass OC emission, LIOUSSE, December  '
         
         write(31,'(a)') 'endvars'
         close(31)
      endif
c
      return
      end
