      program ICBC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   ICBC is the third component of the REGional Climate Modeling       !
!   (RegCM) system version 3.0 and used to access archived global      !
!   analysed datasets at regular latitude-longititude (NNRP1, NNRP2,   !
!   ERA40, ERAIN,EIN75,EIN15, EIN25, GFS11)                            !
!   or original T159 (N80) datasets(ERAHI),                            !
!   or T42 datasets at Gaussian grids (ECWCRP, simply ECMWF), as well  !
!   as NEST run from previous FVGCM run (FVGCM), ECHAM5 run  (EH5OM)   !
!   and RegCM run (FNEST).                                             !
!                                                                      !
!   The present ICBC code could treat NNRP1, NNRP2, ECWCRP, ERA40,     !
!   ERAIN, EIN75, EIN15, EIN25, GFS11, ERAHI, FVGCM, EH5OM, and RegCM  !
!   datasets,  4 times daily.                                          !
!                                                                      !
!                        Xunqiang Bi, ESP group, Abdus Salam ICTP      !
!                                                October 07, 2009      !
!                                                                      !
!   NNRP1: NCEP/NCAR Reanalysis datasets are available at:             !
!          ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/            !
!          Current holdings: 1948 - present, 2.5x2.5L13, netCDF.       !
!   NNRP2: NCEP/DOE AMIP-II Reanalysis (Reanalysis-2) are at:          !
!          ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2/           !
!          Current holdings: 1979 - 2007, 2.5x2.5L13, netCDF.          !
!   NRP2W: Small Window (instead of global) of NNRP1/2 to save disk    !
!          space. (For example, African window: 40W,80E;60S,70N)       !
!   ECMWF: ECMWF TOGA/WCRP Uninitialized Data - (ECWCRP)               !
!          NCAR MSS:/TRENBERT/CTEC/ , ET42yymmdd, where yy=year,       !
!          mm=month, dd=day=01,04,07,10,13,16,19,22,25,28, or 31       !
!          Current holdings: January, 1993 - December, 1997            !
!          Reformatted by PWC/ICTP to direct-access binary,            !
!          T42L15, Gaussian Grid.                                      !
!   EH5OM: EH5OM run by the MPI at Hamburg, T63, Gaussian grid.        !
!          For present day  run: 1941 - 2000;                           !
!          For A1B scenario run: 2001 - 2100.                           !
!          17 pressure levels, 4 times daily, direct-access binary.    !
!   ERA40: ECMWF 40 year reanalysis datasets are available at:         !
!          http://data.ecmwf.int/data/d/era40_daily/                   !
!          Current holdings: 01/09/1957 - 31/08/2002,                  !
!          Pressure levels, 2.5x2.5L23, 4 times daily.                 !
!   ERAIN/EIN15: ECMWF INTERIM 10 year reanalysis datasets             !
!          Current holdings: 01/01/1989 - 31/05/2009,                  !
!          Pressure levels, 1.5x1.5L37, 4 times daily.                 !
!   EIN75: ECMWF INTERIM 10 year reanalysis datasets                   !
!          Current holdings: 01/01/1989 - 31/12/2007,                  !
!          Pressure levels, 0.75x0.75L37, 4 times daily.               !
!   EIN25: ECMWF INTERIM 10 year reanalysis datasets                   !
!          Current holdings: 01/01/1989 - 31/12/1998,                  !
!          Pressure levels, 2.5x2.5L37, 4 times daily.                 !
!   GFS11: NCEP Global Forecast System (GFS) product FNL are           !
!                                                available at:         !
!          http://dss.ucar.edu/datasets/ds083.2/data/fnl-yyyymm/       !
!          Current holdings: 01/01/2000 - present,                     !
!          Pressure levels, 1.0x1.0L27, 4 times daily.                 !
!   ERAHI: ECMWF 40 year reanalysis datasets, origigal model level     !
!          fields: T, U, V and log(Ps) are in spectral coefficients    !
!          Oro and Q are at the reduced Gaussian grids.                !
!          T159L60 (N80L60), 01/09/1957 - 31/08/2002.                  !
!   FVGCM: FVGCM run by the PWC group of Abdus Salam ICTP.             !
!          For present day run: 1960 - 1990;                           !
!          For A2          run: 2070 - 2100.                           !
!          1x1.25L18, 4 times daily, direct-access binary.             !
!   FNEST: do Further oneway NESTing from previous RegCM run.          !
!                                                                      !
!   The code need NetCDF library to access ERA40, ERAIN (EIN75/15/25), !
!   NNRP1 and NNRP2 data.                                              !
!   And we have already provided the NetCDF libraries for many         !
!   platforms, if your platform is unfortunately out of our support,   !
!   you need install the netcdf library by yourself.                   !
!   The code need EMOSLIB library to access ERAHI (T159L60) data, we   !
!   have just provided EMOSLIB library for LINUX PGI5 and IBM AIX.     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      include 'icbc.param'
      INTEGER IDATE
      INTEGER NSTART,NNNEND
      COMMON /ICONTR/ NSTART,NNNEND
      INTEGER MDATE
      COMMON /DATENUM/MDATE(299500)
      INTEGER NOUTREC
      COMMON /OUTCOUNT/ NOUTREC

      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
      INTEGER NNN,ISIZE,NUMFILE,NUMBER
      CHARACTER*26 FINAME

      integer iyr,imon,iday,imonold,imonnew,ifile,idatef
      character a50*50, inrcm*25, inrcm2*25, a120*120, a8*8, a7*7
      logical there
      integer isystm,system
!      external system
!
      CALL INITDATE
      CALL FINDDATE(NSTART,IDATE1)
      CALL FINDDATE(NNNEND,IDATE2)
      WRITE(*,*) 'NSTART,NNNEND: ',NSTART,NNNEND
      WRITE(*,*) 'IDATE1,IDATE2: ',IDATE1,IDATE2

      ISIZE=JX*IY*4*(KZ*4+3)
      NUMFILE=2100000000/ISIZE
      NUMFILE=(NUMFILE/20)*20

      open(60,file='SST.RCM',form='unformatted'
     &       ,status='old',ERR=4810)

      CALL COMMHEAD
      IF(DATTYP.EQ.'NNRP1'.or.DATTYP.EQ.'NNRP2'
     &                    .or.DATTYP.EQ.'NRP2W') THEN
         CALL HEADERNC
      ELSE IF(DATTYP.EQ.'ECMWF') THEN
         CALL HEADEREC
      ELSE IF(DATTYP.EQ.'ERA40') THEN
         CALL HEADERERA
      ELSE IF(DATTYP.EQ.'ERAIN'.or.DATTYP.EQ.'EIN15') THEN
         CALL HEADEREIN15
      ELSE IF(DATTYP.EQ.'EIN75') THEN
         CALL HEADEREIN75
      ELSE IF(DATTYP.EQ.'EIN25') THEN
         CALL HEADEREIN25
      ELSE IF(DATTYP.EQ.'GFS11') THEN
         CALL HEADERGFS
      ELSE IF(DATTYP.EQ.'ERAHI') THEN
         CALL HEADEREHI
      ELSE IF(DATTYP.EQ.'EH5OM') THEN
         CALL HEADERMPI(EHSO4)
      ELSE IF(DATTYP.EQ.'FVGCM') THEN
         CALL HEADERFV
      ELSE IF(DATTYP.EQ.'FNEST') THEN
         CALL HEADNEST
      ENDIF
      IF(SSTTYP.EQ.'OI_WK') CALL HEADWK

      inrcm = '../../Commons/regcm.in'
      inrcm2 = '../../Commons/regcm0.in'
      inquire(file=inrcm,exist=there)
      if (there) then
        inquire(file='tmp.in',exist=there)
        if (there) isystm=system('/bin/rm -f tmp.in')
        a120 = 'cat '//inrcm//' | grep -v rest | grep -v idat > tmp.in'
        isystm=system (a120)
        a120 = '/bin/mv -f '//inrcm//' '//inrcm2
        isystm=system(a120)
        OPEN(99,file=inrcm,status='unknown')
        write(99,*) '&restartparam'
        a7 = '.false.'
        a8 = 'ifrest'
        write(99,101) a8,a7
        a8 = 'idate0'
        write(99,102)a8,idate1
        a8 = 'idate1'
        write(99,102)a8,idate1
        a8 = 'idate2'
        write(99,102)a8,idate2
 101    format(x,a8,x,'=',4x,a7,',')
 102    format(x,a8,x,'=',x,i10,',')
        close (99)
        a120 = 'cat tmp.in >> '//inrcm
        isystm=system (a120)
        isystm=system('/bin/rm -f tmp.in')
      end if


      inrcm = '../../Commons/regcm.x'
      inrcm2 = '../../Commons/regcm0.x'
      a120 = '/bin/mv -f '//inrcm//' '//inrcm2
      inquire(file=inrcm,exist=there)
      if(there) isystm=system(a120)
      OPEN(99,file=inrcm,status='new')
      a50 = '#!/bin/csh -f'
      write(99,199) a50
      a50 = 'set mydir=$PWD'
      write(99,199) a50
      a50 = 'cd ../Main'
      write(99,199) a50
      a50 = 'make clean'
      write(99,199) a50
      a50 = './MAKECODE'
      write(99,199) a50
      a50 = 'make'
      write(99,199) a50
      a50 = 'cd $mydir'
      write(99,199) a50
      a50 = 'mv ../Main/regcm .'
      write(99,199) a50
      a50 = '/bin/ln -sf ../Input/DOMAIN.INFO fort.10'
      write(99,199) a50
 199  format(a50)
      if(nsg.gt.1.and.nsg.lt.10) then
      write(99,299)'/bin/ln -sf ../Input/DOMAIN',nsg,'.INFO fort.11'
      else if(nsg.ge.10) then
      write(99,399)'/bin/ln -sf ../Input/DOMAIN',nsg,'.INFO fort.11'
      endif
      if(AERTYP(4:4).eq.'1'.or.AERTYP(5:5).eq.'1') then
         write(99,499)'/bin/ln -sf ../Input/AERO.dat AERO.dat'
      endif
 299  format(a27,I1,a13)
 399  format(a27,I2,a13)
 499  format(a38)
      imonold = 0
      ifile = 101
      DO NNN=NSTART,NNNEND
         IDATE=MDATE(NNN)
         iyr = idate/1000000
         imon = idate/10000 - iyr*100
!        IF(MOD(NNN-NSTART,NUMFILE).EQ.0 .or.
!    &     (imon.ne.imonold.and.nnn.lt.nnnend.and.nnn.gt.nstart)) THEN
         IF(NNN.eq.NSTART .or.
     &     (imon.ne.imonold.and.nnn.lt.nnnend.and.nnn.gt.nstart)) THEN
            iday=idate/100-iyr*10000-imon*100
            WRITE(FINAME,100) IDATE
 100        FORMAT('../../Input/ICBC',I10)
            IF(NNN.GT.NSTART) THEN
               IF(DATTYP.EQ.'NNRP1'.or.DATTYP.EQ.'NNRP2') THEN
                  CALL GETNCEP(IDATE)
               ELSE IF(DATTYP.EQ.'NRP2W') THEN
                  CALL GETNCEPW(IDATE)
               ELSE IF(DATTYP.EQ.'ECMWF') THEN
                  CALL GETECWCP(IDATE)
               ELSE IF(DATTYP.EQ.'ERA40') THEN
                  CALL GETERA40(IDATE)
               ELSE IF(DATTYP.EQ.'ERAIN'.or.DATTYP.EQ.'EIN15') THEN
                  CALL GETEIN15(IDATE)
               ELSE IF(DATTYP.EQ.'EIN75') THEN
                  CALL GETEIN75(IDATE)
               ELSE IF(DATTYP.EQ.'EIN25') THEN
                  CALL GETEIN25(IDATE)
               ELSE IF(DATTYP.EQ.'GFS11') THEN
                  CALL GETGFS11(IDATE)
               ELSE IF(DATTYP.EQ.'ERAHI') THEN
                  CALL GETERAHI(IDATE)
               ELSE IF(DATTYP.EQ.'EH5OM') THEN
                  CALL GETEH5OM(IDATE)
               ELSE IF(DATTYP.EQ.'FVGCM') THEN
                  CALL GETFVGCM(IDATE)
               ELSE IF(DATTYP.EQ.'FNEST') THEN
                  CALL GET_NEST(IDATE,0)
               ENDIF
            ENDIF
            imonnew = imon + 1
            if (imon.ge.12) then
              imonnew = 1
              iyr = iyr + 1
            end if
            idatef = iyr*1000000 + imonnew*10000 + 100
            IF(imon.eq.1.or.imon.eq.3.or.imon.eq.5.or.imon.eq.7.or.
     &         imon.eq.8.or.imon.eq.10.or.imon.eq.12) THEN
               NUMBER=(32-iday)*4+1
            ELSE IF(imon.eq.4.or.imon.eq.6.or.
     &              imon.eq.9.or.imon.eq.11) THEN
               NUMBER=(31-iday)*4+1
            ELSE
               IF(mod(iyr,4).eq.0) THEN
                  NUMBER=(30-iday)*4+1
               ELSE
                  NUMBER=(29-iday)*4+1
               ENDIF
               IF(mod(iyr,100).eq.0) NUMBER=(29-iday)*4+1
               IF(mod(iyr,400).eq.0) NUMBER=(30-iday)*4+1
            ENDIF
            IF(IGRADS.EQ.1) CALL GRADSCTL(FINAME,IDATE,NUMBER)
            call fexist(FINAME)
            OPEN(64,file=FINAME,form='unformatted',status='unknown'
     &             ,recl=JX*IY*ibyte,access='direct')
            write(a50,198) finame(4:),ifile
 198        format('/bin/ln -sf ',a26,' fort.',i3)
            write(99,199) a50
            imonold = imon
            NOUTREC=0
            ifile = ifile + 1
         ENDIF
         IF(DATTYP.EQ.'NNRP1'.or.DATTYP.EQ.'NNRP2') THEN
            CALL GETNCEP(IDATE)
         ELSE IF(DATTYP.EQ.'NRP2W') THEN
            CALL GETNCEPW(IDATE)
         ELSE IF(DATTYP.EQ.'ECMWF') THEN
            CALL GETECWCP(IDATE)
         ELSE IF(DATTYP.EQ.'ERA40') THEN
            CALL GETERA40(IDATE)
         ELSE IF(DATTYP.EQ.'ERAIN'.or.DATTYP.EQ.'EIN15') THEN
            CALL GETEIN15(IDATE)
         ELSE IF(DATTYP.EQ.'EIN75') THEN
            CALL GETEIN75(IDATE)
         ELSE IF(DATTYP.EQ.'EIN25') THEN
            CALL GETEIN25(IDATE)
         ELSE IF(DATTYP.EQ.'GFS11') THEN
            CALL GETGFS11(IDATE)
         ELSE IF(DATTYP.EQ.'ERAHI') THEN
            CALL GETERAHI(IDATE)
         ELSE IF(DATTYP.EQ.'EH5OM') THEN
            CALL GETEH5OM(IDATE)
         ELSE IF(DATTYP.EQ.'FVGCM') THEN
            CALL GETFVGCM(IDATE)
         ELSE IF(DATTYP.EQ.'FNEST') THEN
            CALL GET_NEST(IDATE,1)
         ENDIF
      ENDDO
      a50 = './regcm<./regcm.in'
      write(99,199) a50
      a50 = 'chmod ugo+x '//inrcm
      isystm=system(a50)
      close(99)
      
      STOP
 4810 PRINT *,'ERROR OPENING SST.RCM FILE'
      STOP '4810 IN PROGRAM ICBC'
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE BILINX(B3,B2,ALON,ALAT,HLON,HLAT,NLON,NLAT,JX,IY,LLEV)
      implicit none
!
!  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
!  RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF GRID2.
!  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
!  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT.
!  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
!  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
!  AS WEIGHTING FACTORS IN THE INTERPOLATION.
!  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!  INTERPOLATED BECAUSE ALAT AND ALON ARE NOT DEFINED FOR
!  THE CROSS POINTS IN THE RegCM MODEL.
!
!  B2(JX,IY,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!  B3(JX,IY,NLEV) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
!  HLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!  HLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!  P.........EAST-WEST WEIGHTING FACTOR.
!  Q.........NORTH-SOUTH WEIGHTING FACTOR.
!  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
!
      INTEGER NLON,NLAT,JX,IY,LLEV
      REAL    B3(JX,IY,LLEV),B2(NLON,NLAT,LLEV)
      REAL    ALON(JX,IY),ALAT(JX,IY),HLON(NLON),HLAT(NLAT)
      INTEGER I,J,I1,II,I2,J1,JJ,J2,L
      REAL    P1,P2,Q1,Q2,AVE
!
      DO 120 J=1,IY
      DO 110 I=1,JX

      I1 = 1000
      DO II=1,NLON-1
         IF(ALON(I,J).GE.HLON(II).AND.ALON(I,J).LT.HLON(II+1))THEN
            P1=ALON(I,J)-HLON(II)
            P2=HLON(II+1)-ALON(I,J)
            I1=II
            I2=II+1
            GOTO 1000
         ELSE IF(ALON(I,J).GE.HLON(II)-360.AND.
     &           ALON(I,J).LT.HLON(II+1)-360.) THEN
            P1=ALON(I,J)-(HLON(II)-360.)
            P2=(HLON(II+1)-360.)-ALON(I,J)
            I1=II
            I2=II+1
            GOTO 1000
         ELSE IF(ALON(I,J).GE.HLON(II)+360.AND.
     &           ALON(I,J).LT.HLON(II+1)+360.) THEN
            P1=ALON(I,J)-(HLON(II)+360.)
            P2=(HLON(II+1)+360.)-ALON(I,J)
            I1=II
            I2=II+1
            GOTO 1000
         ENDIF
      ENDDO
      IF(ALON(I,J).GE.HLON(NLON).AND.
     &   ALON(I,J).LT.HLON(1)+360.) THEN
            P1=ALON(I,J)-HLON(NLON)
            P2=(HLON(1)+360.)-ALON(I,J)
            I1=NLON
            I2=1
            GOTO 1000
      ELSE IF(ALON(I,J).GE.HLON(NLON)+360.AND.
     &        ALON(I,J).LT.HLON(1)+720.) THEN
            P1=ALON(I,J)-(HLON(NLON)+360.)
            P2=(HLON(1)+720.)-ALON(I,J)
            I1=NLON
            I2=1
            GOTO 1000
      ELSE IF(ALON(I,J).GE.HLON(NLON)-360.AND.
     &        ALON(I,J).LT.HLON(1)) THEN
            P1=ALON(I,J)-(HLON(NLON)-360.)
            P2=HLON(1)-ALON(I,J)
            I1=NLON
            I2=1
            GOTO 1000
      ELSE IF(ALON(I,J).GE.HLON(NLON)-720.AND.
     &        ALON(I,J).LT.HLON(1)-360.) THEN
            P1=ALON(I,J)-(HLON(NLON)-720.)
            P2=(HLON(1)-360.)-ALON(I,J)
            I1=NLON
            I2=1
            GOTO 1000
      ENDIF
1000  CONTINUE
      IF(I1.EQ.1000) STOP 'Could not find the right longitute'
      J1 = 1000
      DO JJ=1,NLAT-1
         IF(ALAT(I,J).GE.HLAT(JJ).AND.ALAT(I,J).LT.HLAT(JJ+1))THEN
            Q1=ALAT(I,J)-HLAT(JJ)
            Q2=HLAT(JJ+1)-ALAT(I,J)
            J1=JJ
            J2=JJ+1
         ELSE IF(ALAT(I,J).LT.HLAT(1)) THEN
            IF(ALAT(I,J).GE.-90.) THEN
               Q1=ALAT(I,J)+90.
               Q2=HLAT(1)-ALAT(I,J)
               J1=0
               J2=1
            ENDIF
         ELSE IF(ALAT(I,J).GT.HLAT(NLAT)) THEN
            IF(ALAT(I,J).LE.90.) THEN
               Q1=ALAT(I,J)-HLAT(NLAT)
               Q2=90.-ALAT(I,J)
               J1=NLAT
               J2=NLAT+1
            ENDIF
         ENDIF
      ENDDO
      IF(J1.EQ.1000) STOP 'Could not find the right latitude'
      IF(J1.GT.0.AND.J1.LT.NLAT) THEN
         DO L=1,LLEV
            B3(I,J,L)=( (B2(I1,J1,L)*P2+B2(I2,J1,L)*P1)*Q2
     &                 +(B2(I1,J2,L)*P2+B2(I2,J2,L)*P1)*Q1 )
     &               /(P1+P2)/(Q1+Q2)
         ENDDO
      ELSE IF(J1.EQ.0) THEN
         DO L=1,LLEV
            AVE=0.0
            DO II=1,NLON
               AVE=AVE+B2(II,1,L)
            ENDDO
            AVE=AVE/FLOAT(NLON)
            B3(I,J,L)=( (AVE*(P1+P2))*Q2
     &                 +(B2(I1,J2,L)*P2+B2(I2,J2,L)*P1)*Q1 )
     &               /(P1+P2)/(Q1+Q2)
         ENDDO
      ELSE IF(J1.EQ.NLAT) THEN
         DO L=1,LLEV
            AVE=0.0
            DO II=1,NLON
               AVE=AVE+B2(II,NLAT,L)
            ENDDO
            AVE=AVE/FLOAT(NLON)
            B3(I,J,L)=( (B2(I1,J1,L)*P2+B2(I2,J1,L)*P1)*Q2
     &                + (AVE*(P1+P2))*Q1 )
     &               /(P1+P2)/(Q1+Q2)
         ENDDO
      ENDIF
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE CRESSMCR(B3,B2,ALON,ALAT,GLON,GLAT,JX,IY
     &                   ,I1UR,I1UL,I1DR,I1DL,J1UR,J1UL,J1DR,J1DL
     &                   ,D1XT,D1Xa,D1Xb,D1Xc,D1Xd,NLON,NLAT,NLEV)
      implicit none
!
!  FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
!  THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
!  THE CLOSEST ONE HAS BIGGEST WEIGHT.
!
!  B2(JX,IY,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or irregular GRID.
!  B3(JX,IY,NLEV) IS THE OUTPUT FIELD ON new (regular or irregular) GRID.
!  GLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!  GLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!
      INTEGER NLON,NLAT,JX,IY,NLEV
      REAL    B3(JX,IY,NLEV,5),B2(NLON,NLAT,NLEV,5)
      REAL    ALON(JX,IY),ALAT(JX,IY),GLON(NLON,NLAT),GLAT(NLON,NLAT)
      INTEGER I1UR(JX,IY),I1UL(JX,IY),I1DR(JX,IY),I1DL(JX,IY)
      INTEGER J1UR(JX,IY),J1UL(JX,IY),J1DR(JX,IY),J1DL(JX,IY)
      REAL    D1XT(JX,IY)
      REAL    D1Xa(JX,IY),D1Xb(JX,IY),D1Xc(JX,IY),D1Xd(JX,IY)
      INTEGER I,J,K,L,M,N,MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      REAL    AAA
!
      REAL    GLONMX,GLONMN,GLATMX,GLATMN
      REAL    ALONMX,ALONMN,ALATMX,ALATMN,PI
      INTEGER IMXMN,LCROSS,LDOT
      COMMON /MXNCOM/ GLONMX,GLONMN,GLATMX,GLATMN
     &               ,ALONMX,ALONMN,ALATMX,ALATMN
     &               ,PI,IMXMN,LCROSS,LDOT
      REAL    DIST,DISTa,DISTb,DISTc,DISTd
!
!     Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
!
      IF(IMXMN.EQ.0) THEN
         PI=ATAN(1.)*4.
         GLONMX=-361.
         GLONMN= 361.
         DO N=1,NLAT
         DO M=1,NLON
            IF(GLONMX.LT.GLON(M,N)) GLONMX=GLON(M,N)
            IF(GLONMN.GT.GLON(M,N)) GLONMN=GLON(M,N)
         ENDDO
         ENDDO
         ALONMX=-361.
         ALONMN= 361.
         DO J=1,IY
         DO I=1,JX
            IF(ALONMX.LT.ALON(I,J)) ALONMX=ALON(I,J)
            IF(ALONMN.GT.ALON(I,J)) ALONMN=ALON(I,J)
         ENDDO
         ENDDO
         GLATMX=-91.
         GLATMN= 91.
         DO N=1,NLAT
         DO M=1,NLON
            IF(GLATMX.LT.GLAT(M,N)) GLATMX=GLAT(M,N)
            IF(GLATMN.GT.GLAT(M,N)) GLATMN=GLAT(M,N)
         ENDDO
         ENDDO
         ALATMX=-91.
         ALATMN= 91.
         DO J=1,IY
         DO I=1,JX
            IF(ALATMX.LT.ALAT(I,J)) ALATMX=ALAT(I,J)
            IF(ALATMN.GT.ALAT(I,J)) ALATMN=ALAT(I,J)
         ENDDO
         ENDDO
         WRITE(*,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
         WRITE(*,*)  GLONMN,ALONMN,ALONMX,GLONMX
         WRITE(*,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
         WRITE(*,*)  GLATMN,ALATMN,ALATMX,GLATMX
         IMXMN=1
      ENDIF
      IF(LCROSS.EQ.0) THEN
         DO 120 J=1,IY
         DO 110 I=1,JX

         MUR=1000
         NUR=1000
         MUL=1000
         NUL=1000
         MDR=1000
         NDR=1000
         MDL=1000
         NDL=1000

         DISTa=1.E8
         DISTb=1.E8
         DISTc=1.E8
         DISTd=1.E8
         DO N=2,NLAT
         DO M=2,NLON
        IF((GLON(M,N).GE.ALON(I,J).AND.GLON(M,N)-ALON(I,J).lt.10.).AND.
     &     (GLAT(M,N).GE.ALAT(I,J).AND.GLAT(M,N)-ALAT(I,J).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTa.GT.AAA) THEN
               DISTa=AAA
               MUR=M
               NUR=N
            ENDIF
         ENDIF
        IF((GLON(M,N).LT.ALON(I,J).AND.ALON(I,J)-GLON(M,N).lt.10.).AND.
     &     (GLAT(M,N).GE.ALAT(I,J).AND.GLAT(M,N)-ALAT(I,J).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTb.GT.AAA) THEN
               DISTb=AAA
               MUL=M
               NUL=N
            ENDIF
         ENDIF
        IF((GLON(M,N).GE.ALON(I,J).AND.GLON(M,N)-ALON(I,J).lt.10.).AND.
     &     (GLAT(M,N).LT.ALAT(I,J).AND.ALAT(I,J)-GLAT(M,N).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTc.GT.AAA) THEN
               DISTc=AAA
               MDR=M
               NDR=N
            ENDIF
         ENDIF
        IF((GLON(M,N).LT.ALON(I,J).AND.ALON(I,J)-GLON(M,N).lt.10.).AND.
     &     (GLAT(M,N).LT.ALAT(I,J).AND.ALAT(I,J)-GLAT(M,N).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTd.GT.AAA) THEN
               DISTd=AAA
               MDL=M
               NDL=N
            ENDIF
         ENDIF
         ENDDO
         ENDDO
         DIST=amin1(DISTa,DISTb,DISTc,DISTd)
      
         I1UR(I,J)=MUR
         J1UR(I,J)=NUR
         I1UL(I,J)=MUL
         J1UL(I,J)=NUL
         I1DR(I,J)=MDR
         J1DR(I,J)=NDR
         I1DL(I,J)=MDL
         J1DL(I,J)=NDL
         D1XT(I,J)=DIST
         D1Xa(I,J)=DISTa
         D1Xb(I,J)=DISTb
         D1Xc(I,J)=DISTc
         D1Xd(I,J)=DISTd
         DO L=1,4
         DO K=1,NLEV
         IF(DIST.GT.0.000001) THEN
            B3(I,J,K,L)=(B2(MUR,NUR,K,L)/DISTa+B2(MUL,NUL,K,L)/DISTb
     &                  +B2(MDR,NDR,K,L)/DISTc+B2(MDL,NDL,K,L)/DISTd)
     &               /(1./DISTa+1./DISTb+1./DISTc+1./DISTd)
         ELSE
            IF(DIST.EQ.DISTa) THEN
               B3(I,J,K,L)=B2(MUR,NUR,K,L)
            ELSE IF(DIST.EQ.DISTb) THEN
               B3(I,J,K,L)=B2(MUL,NUL,K,L)
            ELSE IF(DIST.EQ.DISTc) THEN
               B3(I,J,K,L)=B2(MDR,NDR,K,L)
            ELSE IF(DIST.EQ.DISTd) THEN
               B3(I,J,K,L)=B2(MDL,NDL,K,L)
            ENDIF
         ENDIF
         ENDDO
         ENDDO
  110    CONTINUE
  120    CONTINUE
         LCROSS=1
      ELSE
         DO 220 J=1,IY
         DO 210 I=1,JX

         MUR=I1UR(I,J)
         NUR=J1UR(I,J)
         MUL=I1UL(I,J)
         NUL=J1UL(I,J)
         MDR=I1DR(I,J)
         NDR=J1DR(I,J)
         MDL=I1DL(I,J)
         NDL=J1DL(I,J)
         DIST=D1XT(I,J)
         DISTa=D1Xa(I,J)
         DISTb=D1Xb(I,J)
         DISTc=D1Xc(I,J)
         DISTd=D1Xd(I,J)

         DO L=1,4
         DO K=1,NLEV
         IF(DIST.GT.0.000001) THEN
            B3(I,J,K,L)=(B2(MUR,NUR,K,L)/DISTa+B2(MUL,NUL,K,L)/DISTb
     &                  +B2(MDR,NDR,K,L)/DISTc+B2(MDL,NDL,K,L)/DISTd)
     &               /(1./DISTa+1./DISTb+1./DISTc+1./DISTd)
         ELSE
            IF(DIST.EQ.DISTa) THEN
               B3(I,J,K,L)=B2(MUR,NUR,K,L)
            ELSE IF(DIST.EQ.DISTb) THEN
               B3(I,J,K,L)=B2(MUL,NUL,K,L)
            ELSE IF(DIST.EQ.DISTc) THEN
               B3(I,J,K,L)=B2(MDR,NDR,K,L)
            ELSE IF(DIST.EQ.DISTd) THEN
               B3(I,J,K,L)=B2(MDL,NDL,K,L)
            ENDIF
         ENDIF
         ENDDO
         ENDDO
  210    CONTINUE
  220    CONTINUE
      ENDIF

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE CRESSMDT(B3,B2,ALON,ALAT,GLON,GLAT,JX,IY
     &                   ,I1UR,I1UL,I1DR,I1DL,J1UR,J1UL,J1DR,J1DL
     &                   ,D1XT,D1Xa,D1Xb,D1Xc,D1Xd,NLON,NLAT,NLEV)
      implicit none
!
!  FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
!  THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
!  THE CLOSEST ONE HAS BIGGEST WEIGHT.
!
!  B2(JX,IY,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or irregular GRID.
!  B3(JX,IY,NLEV) IS THE OUTPUT FIELD ON new (regular or irregular) GRID.
!  GLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!  GLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
!
      INTEGER NLON,NLAT,JX,IY,NLEV
      REAL    B3(JX,IY,NLEV,5),B2(NLON,NLAT,NLEV,5)
      REAL    ALON(JX,IY),ALAT(JX,IY),GLON(NLON,NLAT),GLAT(NLON,NLAT)
      INTEGER I1UR(JX,IY),I1UL(JX,IY),I1DR(JX,IY),I1DL(JX,IY)
      INTEGER J1UR(JX,IY),J1UL(JX,IY),J1DR(JX,IY),J1DL(JX,IY)
      REAL    D1XT(JX,IY)
      REAL    D1Xa(JX,IY),D1Xb(JX,IY),D1Xc(JX,IY),D1Xd(JX,IY)
      INTEGER I,J,K,L,M,N,MUR,NUR,MUL,NUL,MDR,NDR,MDL,NDL
      REAL    AAA
!
      REAL    GLONMX,GLONMN,GLATMX,GLATMN
      REAL    ALONMX,ALONMN,ALATMX,ALATMN,PI
      INTEGER IMXMN,LCROSS,LDOT
      COMMON /MXNCOM/ GLONMX,GLONMN,GLATMX,GLATMN
     &               ,ALONMX,ALONMN,ALATMX,ALATMN
     &               ,PI,IMXMN,LCROSS,LDOT
      REAL    DIST,DISTa,DISTb,DISTc,DISTd
!
!     Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
!
      IF(IMXMN.EQ.0) THEN
         PI=ATAN(1.)*4.
         GLONMX=-361.
         GLONMN= 361.
         DO N=1,NLAT
         DO M=1,NLON
            IF(GLONMX.LT.GLON(M,N)) GLONMX=GLON(M,N)
            IF(GLONMN.GT.GLON(M,N)) GLONMN=GLON(M,N)
         ENDDO
         ENDDO
         ALONMX=-361.
         ALONMN= 361.
         DO J=1,IY
         DO I=1,JX
            IF(ALONMX.LT.ALON(I,J)) ALONMX=ALON(I,J)
            IF(ALONMN.GT.ALON(I,J)) ALONMN=ALON(I,J)
         ENDDO
         ENDDO
         GLATMX=-91.
         GLATMN= 91.
         DO N=1,NLAT
         DO M=1,NLON
            IF(GLATMX.LT.GLAT(M,N)) GLATMX=GLAT(M,N)
            IF(GLATMN.GT.GLAT(M,N)) GLATMN=GLAT(M,N)
         ENDDO
         ENDDO
         ALATMX=-91.
         ALATMN= 91.
         DO J=1,IY
         DO I=1,JX
            IF(ALATMX.LT.ALAT(I,J)) ALATMX=ALAT(I,J)
            IF(ALATMN.GT.ALAT(I,J)) ALATMN=ALAT(I,J)
         ENDDO
         ENDDO
         WRITE(*,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
         WRITE(*,*)  GLONMN,ALONMN,ALONMX,GLONMX
         WRITE(*,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
         WRITE(*,*)  GLATMN,ALATMN,ALATMX,GLATMX
         IMXMN=1
      ENDIF
      IF(LDOT.EQ.0) THEN
         DO 120 J=1,IY
         DO 110 I=1,JX

         MUR=1000
         NUR=1000
         MUL=1000
         NUL=1000
         MDR=1000
         NDR=1000
         MDL=1000
         NDL=1000

         DISTa=1.E8
         DISTb=1.E8
         DISTc=1.E8
         DISTd=1.E8
         DO N=2,NLAT
         DO M=2,NLON
        IF((GLON(M,N).GE.ALON(I,J).AND.GLON(M,N)-ALON(I,J).lt.10.).AND.
     &     (GLAT(M,N).GE.ALAT(I,J).AND.GLAT(M,N)-ALAT(I,J).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTa.GT.AAA) THEN
               DISTa=AAA
               MUR=M
               NUR=N
            ENDIF
         ENDIF
        IF((GLON(M,N).LT.ALON(I,J).AND.ALON(I,J)-GLON(M,N).lt.10.).AND.
     &     (GLAT(M,N).GE.ALAT(I,J).AND.GLAT(M,N)-ALAT(I,J).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTb.GT.AAA) THEN
               DISTb=AAA
               MUL=M
               NUL=N
            ENDIF
         ENDIF
        IF((GLON(M,N).GE.ALON(I,J).AND.GLON(M,N)-ALON(I,J).lt.10.).AND.
     &     (GLAT(M,N).LT.ALAT(I,J).AND.ALAT(I,J)-GLAT(M,N).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTc.GT.AAA) THEN
               DISTc=AAA
               MDR=M
               NDR=N
            ENDIF
         ENDIF
        IF((GLON(M,N).LT.ALON(I,J).AND.ALON(I,J)-GLON(M,N).lt.10.).AND.
     &     (GLAT(M,N).LT.ALAT(I,J).AND.ALAT(I,J)-GLAT(M,N).lt.10.))THEN
            AAA=((GLON(M,N)-ALON(I,J))
     &          *COS((GLAT(M,N)+ALAT(I,J))/360.*PI))**2
     &         +(GLAT(M,N)-ALAT(I,J))**2
            IF(DISTd.GT.AAA) THEN
               DISTd=AAA
               MDL=M
               NDL=N
            ENDIF
         ENDIF
         ENDDO
         ENDDO
         DIST=amin1(DISTa,DISTb,DISTc,DISTd)
      
         I1UR(I,J)=MUR
         J1UR(I,J)=NUR
         I1UL(I,J)=MUL
         J1UL(I,J)=NUL
         I1DR(I,J)=MDR
         J1DR(I,J)=NDR
         I1DL(I,J)=MDL
         J1DL(I,J)=NDL
         D1XT(I,J)=DIST
         D1Xa(I,J)=DISTa
         D1Xb(I,J)=DISTb
         D1Xc(I,J)=DISTc
         D1Xd(I,J)=DISTd
         DO L=1,3
         DO K=1,NLEV
         IF(DIST.GT.0.000001) THEN
            B3(I,J,K,L)=(B2(MUR,NUR,K,L)/DISTa+B2(MUL,NUL,K,L)/DISTb
     &                  +B2(MDR,NDR,K,L)/DISTc+B2(MDL,NDL,K,L)/DISTd)
     &               /(1./DISTa+1./DISTb+1./DISTc+1./DISTd)
         ELSE
            IF(DIST.EQ.DISTa) THEN
               B3(I,J,K,L)=B2(MUR,NUR,K,L)
            ELSE IF(DIST.EQ.DISTb) THEN
               B3(I,J,K,L)=B2(MUL,NUL,K,L)
            ELSE IF(DIST.EQ.DISTc) THEN
               B3(I,J,K,L)=B2(MDR,NDR,K,L)
            ELSE IF(DIST.EQ.DISTd) THEN
               B3(I,J,K,L)=B2(MDL,NDL,K,L)
            ENDIF
         ENDIF
         ENDDO
         ENDDO
  110    CONTINUE
  120    CONTINUE
         LDOT=1
      ELSE
         DO 220 J=1,IY
         DO 210 I=1,JX

         MUR=I1UR(I,J)
         NUR=J1UR(I,J)
         MUL=I1UL(I,J)
         NUL=J1UL(I,J)
         MDR=I1DR(I,J)
         NDR=J1DR(I,J)
         MDL=I1DL(I,J)
         NDL=J1DL(I,J)
         DIST=D1XT(I,J)
         DISTa=D1Xa(I,J)
         DISTb=D1Xb(I,J)
         DISTc=D1Xc(I,J)
         DISTd=D1Xd(I,J)

         DO L=1,3
         DO K=1,NLEV
         IF(DIST.GT.0.000001) THEN
            B3(I,J,K,L)=(B2(MUR,NUR,K,L)/DISTa+B2(MUL,NUL,K,L)/DISTb
     &                  +B2(MDR,NDR,K,L)/DISTc+B2(MDL,NDL,K,L)/DISTd)
     &               /(1./DISTa+1./DISTb+1./DISTc+1./DISTd)
         ELSE
            IF(DIST.EQ.DISTa) THEN
               B3(I,J,K,L)=B2(MUR,NUR,K,L)
            ELSE IF(DIST.EQ.DISTb) THEN
               B3(I,J,K,L)=B2(MUL,NUL,K,L)
            ELSE IF(DIST.EQ.DISTc) THEN
               B3(I,J,K,L)=B2(MDR,NDR,K,L)
            ELSE IF(DIST.EQ.DISTd) THEN
               B3(I,J,K,L)=B2(MDL,NDL,K,L)
            ENDIF
         ENDIF
         ENDDO
         ENDDO
  210    CONTINUE
  220    CONTINUE
      ENDIF

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE COMMHEAD
      IMPLICIT NONE
!
!  X X X X X   SET 2 : PARAMETERS FOR RegCM DATASET   X X X X
!
!HH   change the vertical levels and the model domain size.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RegCM DATASET.
      include 'icbc.param'
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER K
!
      CALL GRIDML(XLON,XLAT,DLON,DLAT,TOPOGM,TOPOSDGM,XLANDU
     &           ,MSFX,PTOP,SIGMAF,CLON,CLAT,PLON,PLAT,DELX,GRDFAC
     &           ,JX,IY,KZ,DATTYP,LGTYPE,igrads,ibigend,ibyte)
!
      DO K=1,KZ
         SIGMA2(K) = 0.5*(SIGMAF(K+1)+SIGMAF(K))
         DSIGMA(K) = SIGMAF(K+1)-SIGMAF(K)
      ENDDO
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE FINDDATE(NPOS,IDATE)
      IMPLICIT NONE
      INTEGER MDATE
      COMMON /DATENUM/MDATE(299500)
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
      IF(I.GT.299500) GOTO 100
      GO TO 10
 100  WRITE(*,*) 'ERROR IN FINDDATE'
      STOP
 200  RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETECWCP(IDATE)
      implicit none
      INTEGER IDATE
      INTEGER NYEAR,MONTH,NDAY,NHOUR,NREC
      CHARACTER*12 FINM(12,5)
      DATA    FINM/
     & 'ECT421993JAN','ECT421993FEB','ECT421993MAR','ECT421993APR',
     & 'ECT421993MAY','ECT421993JUN','ECT421993JUL','ECT421993AUG',
     & 'ECT421993SEP','ECT421993OCT','ECT421993NOV','ECT421993DEC',
     & 'ECT421994JAN','ECT421994FEB','ECT421994MAR','ECT421994APR',
     & 'ECT421994MAY','ECT421994JUN','ECT421994JUL','ECT421994AUG',
     & 'ECT421994SEP','ECT421994OCT','ECT421994NOV','ECT421994DEC',
     & 'ECT421995JAN','ECT421995FEB','ECT421995MAR','ECT421995APR',
     & 'ECT421995MAY','ECT421995JUN','ECT421995JUL','ECT421995AUG',
     & 'ECT421995SEP','ECT421995OCT','ECT421995NOV','ECT421995DEC',
     & 'ECT421996JAN','ECT421996FEB','ECT421996MAR','ECT421996APR',
     & 'ECT421996MAY','ECT421996JUN','ECT421996JUL','ECT421996AUG',
     & 'ECT421996SEP','ECT421996OCT','ECT421996NOV','ECT421996DEC',
     & 'ECT421997JAN','ECT421997FEB','ECT421997MAR','ECT421997APR',
     & 'ECT421997MAY','ECT421997JUN','ECT421997JUL','ECT421997AUG',
     & 'ECT421997SEP','ECT421997OCT','ECT421997NOV','ECT421997DEC'/
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 : PARAMETERS FOR ECMWF DATASET   X X X X
! A1
      INTEGER NLEV1,NLAT1,NLON1
!as   PARAMETER (NLEV1=14,NLAT1=64,NLON1=128)
      PARAMETER (NLEV1=15,NLAT1=64,NLON1=128)
!
!  NLEV1  = NUMBER OF PRESSURE LEVELS IN ECMWF DATASET.
!  NLAT1  = NUMBER OF LATITUDES ON ECMWF GRID.
!  NLON1  = NUMBER OF LONGITUDES ON ECMWF GRID.
!
!  X X X X X   SET 2 : PARAMETERS FOR RegCM DATASET   X X X X
!
!HH   change the vertical levels and the model domain size.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RegCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    HLON,HLAT,SIGMA1,SIGMAR
      COMMON /SPECFS/HLON(NLON1),HLAT(NLAT1),SIGMA1(NLEV1),SIGMAR(NLEV1)
!
!     DOMAIN VARIABLES FOR RegCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)

! DIMENSION SURFACE TEMPERATURE ON ECMWF SURFACE; NOT GIVEN BY ECMWF DATA
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!B1
      REAL    T1,Q1,Z1,U1,V1,W1
      COMMON /ECVARS/ T1(NLON1,NLAT1,NLEV1),Q1(NLON1,NLAT1,NLEV1),
     &                Z1(NLON1,NLAT1,NLEV1),
     &                U1(NLON1,NLAT1,NLEV1),V1(NLON1,NLAT1,NLEV1),
     &                W1(NLON1,NLAT1,NLEV1)
      REAL    B2(NLON1,NLAT1,NLEV1*3)
      EQUIVALENCE (B2(1,1,1),T1(1,1,1))
      REAL    D2(NLON1,NLAT1,NLEV1*2)
      EQUIVALENCE (D2(1,1,1),U1(1,1,1))
!
!B3
      REAL    T3,Q3,H3,U3,V3,W3,B3PD
      COMMON /VARSB3/ T3(JX,IY,NLEV1),Q3(JX,IY,NLEV1),
     &                H3(JX,IY,NLEV1),
     &                U3(JX,IY,NLEV1),V3(JX,IY,NLEV1),
     &                W3(JX,IY,NLEV1),B3PD(JX,IY)
      REAL    B3(JX,IY,NLEV1*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,NLEV1*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!B4
      REAL    PS4,T4,Q4,H4,TS4,U4,V4
      COMMON /VARSB4/ PS4(JX,IY),
     &                T4(JX,IY,KZ),Q4(JX,IY,KZ),
     &                H4(JX,IY,KZ),TS4(JX,IY),
     &                U4(JX,IY,KZ),V4(JX,IY,KZ)
!
      INTEGER I,J,K
      INTEGER NYRP,NMOP
      REAL    WT
      logical there
!
!  B1 IS FOR DATA RECORDS FROM ECMWF GAUSSIAN GRID
!  B2 IS FOR LAT-LON GRID WITH ECMWF VERTICAL STRUCTURE
!  B3 IS FOR RegCM HORIZONTAL GRID, BUT ECMWF VERTICAL GRID
!  B4 IS FOR RegCM 3-DIMENSIONAL GRID
!
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-(IDATE/10000)*100
      NHOUR=MOD(IDATE,100)
      
      inquire(file='../DATA/ECWCRP/'//FINM(MONTH,NYEAR-1992)
     &       ,exist=there)
      if(.not.there) then
         write(*,*) '../DATA/ECWCRP/'//FINM(MONTH,NYEAR-1992)
     &             ,' is not available'
         stop
      endif
      OPEN(63,file='../DATA/ECWCRP/'//FINM(MONTH,NYEAR-1992)
     &       ,form='unformatted',recl=NLON1*NLAT1*ibyte
     &       ,access='direct')
      nrec=((NDAY-1)*4+NHOUR/6)*(NLEV1*6+1)
      nrec=nrec+1
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((Z1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((T1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((U1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((V1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((W1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      do k=1,NLEV1
         nrec=nrec+1
         read(63,rec=nrec) ((Q1(i,j,k),i=1,NLON1),j=1,NLAT1)
      enddo
      WRITE(*,*) 'READ IN fields at DATE:',IDATE,' from '
     &          ,FINM(MONTH,NYEAR-1992)
      close(21)
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,HLON,HLAT,NLON1,NLAT1,JX,IY,NLEV1*3)
      CALL BILINX(D3,D2,DLON,DLAT,HLON,HLAT,NLON1,NLAT1,JX,IY,NLEV1*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,NLEV1
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,NLEV1)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!       INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,NLEV1)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
         CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
         CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
         CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2     DETERMINE P* AND HEIGHT.
!
! F3     INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV1)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV1)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV1)
 
      CALL HUMID1(T3,Q3,100.,0.0,SIGMA1,JX,IY,NLEV1)
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV1)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4     DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G      WRITE AN INITIAL FILE FOR THE RegCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADEREC
      implicit none
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 : PARAMETERS FOR ECMWF DATASET   X X X X
! A1
      INTEGER NLEV1,NLAT1,NLON1
!as   PARAMETER (NLEV1=14,NLAT1=64,NLON1=128)
      PARAMETER (NLEV1=15,NLAT1=64,NLON1=128)
!
!  NLEV1  = NUMBER OF PRESSURE LEVELS IN ECMWF DATASET.
!  NLAT1  = NUMBER OF LATITUDES  ON ECMWF GRID.
!  NLON1  = NUMBER OF LONGITUDES ON ECMWF GRID.
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    HLON,HLAT,SIGMA1,SIGMAR
      COMMON /SPECFS/HLON(NLON1),HLAT(NLAT1),SIGMA1(NLEV1),SIGMAR(NLEV1)
!
      INTEGER I,K,KR
!
      HLAT( 1)= -87.8638
      HLAT( 2)= -85.0965
      HLAT( 3)= -82.3129
      HLAT( 4)= -79.5256
      HLAT( 5)= -76.7369
      HLAT( 6)= -73.9475
      HLAT( 7)= -71.1578
      HLAT( 8)= -68.3678
      HLAT( 9)= -65.5776
      HLAT(10)= -62.7874
      HLAT(11)= -59.9970
      HLAT(12)= -57.2066
      HLAT(13)= -54.4162
      HLAT(14)= -51.6257
      HLAT(15)= -48.8352
      HLAT(16)= -46.0447
      HLAT(17)= -43.2542
      HLAT(18)= -40.4636
      HLAT(19)= -37.6731
      HLAT(20)= -34.8825
      HLAT(21)= -32.0919
      HLAT(22)= -29.3014
      HLAT(23)= -26.5108
      HLAT(24)= -23.7202
      HLAT(25)= -20.9296
      HLAT(26)= -18.1390
      HLAT(27)= -15.3484
      HLAT(28)= -12.5578
      HLAT(29)=  -9.76715
      HLAT(30)=  -6.97653
      HLAT(31)=  -4.18592
      HLAT(32)=  -1.39531
      HLAT(33)=   1.39531
      HLAT(34)=   4.18592
      HLAT(35)=   6.97653
      HLAT(36)=   9.76715
      HLAT(37)=  12.5578
      HLAT(38)=  15.3484
      HLAT(39)=  18.1390
      HLAT(40)=  20.9296
      HLAT(41)=  23.7202
      HLAT(42)=  26.5108
      HLAT(43)=  29.3014
      HLAT(44)=  32.0919
      HLAT(45)=  34.8825
      HLAT(46)=  37.6731
      HLAT(47)=  40.4636
      HLAT(48)=  43.2542
      HLAT(49)=  46.0447
      HLAT(50)=  48.8352
      HLAT(51)=  51.6257
      HLAT(52)=  54.4162
      HLAT(53)=  57.2066
      HLAT(54)=  59.9970
      HLAT(55)=  62.7874
      HLAT(56)=  65.5776
      HLAT(57)=  68.3678
      HLAT(58)=  71.1578
      HLAT(59)=  73.9475
      HLAT(60)=  76.7369
      HLAT(61)=  79.5256
      HLAT(62)=  82.3129
      HLAT(63)=  85.0965
      HLAT(64)=  87.8638
      DO I=1,NLON1
         HLON(I)=FLOAT(I-1)*2.8125
      ENDDO
!
      SIGMAR(1)  = .01
      SIGMAR(2)  = .03
      SIGMAR(3)  = .05
      SIGMAR(4)  = .07
      SIGMAR(5)  = .1
      SIGMAR(6)  = .15
      SIGMAR(7)  = .2
      SIGMAR(8)  = .25
      SIGMAR(9)  = .3
      SIGMAR(10) = .4
      SIGMAR(11) = .5
      SIGMAR(12) = .7
      SIGMAR(13) = .85
      SIGMAR(14) = .925
      SIGMAR(15) =1.0
!
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
!HH:ECMWF-T42
!:
      DO K=1,NLEV1
         KR=NLEV1-K+1
         SIGMA1(K)=SIGMAR(KR)
      ENDDO
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETEH5OM(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
      INTEGER NYEAR,MONTH,NDAY,NHOUR,NREC,KREC
      CHARACTER*4 YR_RF(61),YR_A2(100)
      CHARACTER*3 CHmon(12)
      DATA    YR_RF/ '1941','1942','1943','1944','1945',
     &               '1946','1947','1948','1949','1950',
     &               '1951','1952','1953','1954','1955',
     &               '1956','1957','1958','1959','1960',
     &               '1961','1962','1963','1964','1965',
     &               '1966','1967','1968','1969','1970',
     &               '1971','1972','1973','1974','1975',
     &               '1976','1977','1978','1979','1980',
     &               '1981','1982','1983','1984','1985',
     &               '1986','1987','1988','1989','1990',
     &               '1991','1992','1993','1994','1995',
     &               '1996','1997','1998','1999','2000','2001'/
      DATA    YR_A2/ '2001','2002','2003','2004','2005',
     &               '2006','2007','2008','2009','2010',
     &               '2011','2012','2013','2014','2015',
     &               '2016','2017','2018','2019','2020',
     &               '2021','2022','2023','2024','2025',
     &               '2026','2027','2028','2029','2030',
     &               '2031','2032','2033','2034','2035',
     &               '2036','2037','2038','2039','2040',
     &               '2041','2042','2043','2044','2045',
     &               '2046','2047','2048','2049','2050',
     &               '2051','2052','2053','2054','2055',
     &               '2056','2057','2058','2059','2060',
     &               '2061','2062','2063','2064','2065',
     &               '2066','2067','2068','2069','2070',
     &               '2071','2072','2073','2074','2075',
     &               '2076','2077','2078','2079','2080',
     &               '2081','2082','2083','2084','2085',
     &               '2086','2087','2088','2089','2090',
     &               '2091','2092','2093','2094','2095',
     &               '2096','2097','2098','2099','2100'/
      DATA    CHmon/ 'JAN','FEB','MAR','APR','MAY','JUN',
     &               'JUL','AUG','SEP','OCT','NOV','DEC'/
      CHARACTER*21 FINM
      CHARACTER*21 PSNM
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR MPI ECHAM5 CLIMATE CHANGE RUN X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=17,jlat=96,ilon=192)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar!,Wvar
      common /EH5vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEH/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /TMPVAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,         T3(JX,IY,klev)
     &       ,         H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY)!, SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
!      INTEGER NYRP,NMOP
!      REAL    WT

      real    lon0,lon1,lat0,lat1
      INTEGER i0,i1,j0
      COMMON /SZwindow/lon0,lon1,lat0,lat1,i0,i1,j0
      integer numx,numy,ii,i2,j2,I,J,K
!      REAL    temp(192,96)
      real(kind=8)  scale,offset
      integer(kind=2) itmp(192,96)
      logical there

      integer mlev
      parameter (mlev=31)
      REAL    hyai,hybi,hyam,hybm
      integer(kind=4) start,count
      real(kind=4)  PSO4_0
      COMMON /GLOBALEH2/hyai(mlev+1),hybi(mlev+1),hyam(mlev),hybm(mlev),
     &                  PSO4_0(ilon,jlat),start(10),count(10)

      real(kind=4)  sulfate
      COMMON /INPUTSO4/sulfate(ilon,jlat,mlev,12)
      REAL    sulfate1,sulfate2,sulfate3,sulfate4
      REAL    PSO4_2(ilon,jlat),PSO4_3(JX,IY)
      COMMON /WRK_SO4/sulfate1(ilon,jlat,mlev,2),
     &                sulfate2(ilon,jlat,mlev),
     &                sulfate3(JX,IY,mlev),
     &                sulfate4(JX,IY,KZ),PSO4_2,PSO4_3
      character*39 fnso4_RF
      character*44 fnso4_A1B
      character*42 fnso4_A2
      character*42 fnso4_B1
      integer(kind=4) ncid, status
      COMMON /WRK_SO4s/ncid, status
      REAL    PRCM,PMPI,PMPJ
      INTEGER L,K0
      include 'netcdf.inc'
      INTEGER NOUTREC
      COMMON /OUTCOUNT/NOUTREC
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER EH5OM (ECHAM5_MPIOM) HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-(IDATE/10000)*100
      NHOUR=MOD(IDATE,100)

      IF(SSTTYP.EQ.'EH5RF') THEN
        IF(IDATE.LT.1941010106) THEN
          WRITE(*,*) 'EH5RF dataset is just available from 1941010106'
          STOP
        ENDIF
        IF(EHSO4) THEN
           IF(IDATE.LT.1950010100) THEN
              WRITE(*,*) 'the monthly EH5RF sulfate data are just',
     &                   ' available from 1950010100'
              STOP
           ENDIF
        ENDIF
        IF(IDATE.GT.2001010100) THEN
          WRITE(*,*) 'EH5RF dataset is just available  to  2001010100'
          STOP
        ENDIF
      ENDIF
      IF(SSTTYP.EQ.'EH5A2') THEN
        IF(IDATE.LT.2001010100) THEN
          WRITE(*,*) 'EH5A2 dataset is just avaiable from 2001010100'
          STOP
        ENDIF
        IF(IDATE.GE.2101010100) THEN
          WRITE(*,*) 'EH5A2 dataset is just avaiable  to  2100123118'
          STOP
        ENDIF
      ENDIF
      IF(SSTTYP.EQ.'EH5B1') THEN
        IF(IDATE.LT.2001010100) THEN
          WRITE(*,*) 'EH5B1 dataset is just avaiable from 2001010100'
          STOP
        ENDIF
        IF(IDATE.GE.2101010100) THEN
          WRITE(*,*) 'EH5B1 dataset is just avaiable  to  2100123118'
          STOP
        ENDIF
      ENDIF
      IF(SSTTYP.EQ.'EHA1B') THEN
        IF(IDATE.LT.2001010100) THEN
          WRITE(*,*) 'EHA1B dataset is just avaiable from 2001010100'
          STOP
        ENDIF
        IF(IDATE.GE.2101010100) THEN
          WRITE(*,*) 'EHA1B dataset is just avaiable  to  2100123118'
          STOP
        ENDIF
      ENDIF
      numx=nint((lon1-lon0)/1.875)+1
      numy=nint((lat1-lat0)/1.875)+1
      IF(.not.(numx.eq.192.and.numy.eq.96)) THEN
        IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
          IF(SSTTYP.EQ.'EH5RF') THEN
            FINM='RF/'//YR_RF(NYEAR-1940)//'/'//'EH_RF'//
     &                  YR_RF(NYEAR-1940)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='RF/'//YR_RF(NYEAR-1940)//'/'//'EH_PS'//
     &                  YR_RF(NYEAR-1940)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EH5A2') THEN
            FINM='A2/'//YR_A2(NYEAR-2000)//'/'//'EH_A2'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='A2/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EH5B1') THEN
            FINM='B1/'//YR_A2(NYEAR-2000)//'/'//'EH_B1'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='B1/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EHA1B') THEN
            FINM='A1B/'//YR_A2(NYEAR-2000)//'/'//'E_A1B'//
     &                   YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                   YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE
            WRITE(*,*)'ERROR IN SSTTYP'
            STOP
          ENDIF
        ELSE
          IF(.NOT.(MONTH.EQ.1)) THEN
            IF(SSTTYP.EQ.'EH5RF') THEN
               FINM='RF/'//YR_RF(NYEAR-1940)//'/'//'EH_RF'//
     &                     YR_RF(NYEAR-1940)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='RF/'//YR_RF(NYEAR-1940)//'/'//'EH_PS'//
     &                     YR_RF(NYEAR-1940)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EH5A2') THEN
               FINM='A2/'//YR_A2(NYEAR-2000)//'/'//'EH_A2'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='A2/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EH5B1') THEN
               FINM='B1/'//YR_A2(NYEAR-2000)//'/'//'EH_B1'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='B1/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EHA1B') THEN
               FINM='A1B/'//YR_A2(NYEAR-2000)//'/'//'E_A1B'//
     &                      YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EH_PS'//
     &                      YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
          ELSE
            IF(SSTTYP.EQ.'EH5RF') THEN
               FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_RF'//
     &                     YR_RF(NYEAR-1941)//CHmon(12)
               IF(EHSO4)
     &         PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_PS'//
     &                     YR_RF(NYEAR-1941)//CHmon(12)
            ELSE IF(SSTTYP.EQ.'EH5A2') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_RF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_PS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='A2/'//YR_A2(NYEAR-2001)//'/'//'EH_A2'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='A2/'//YR_A2(NYEAR-2001)//'/'//'EH_PS'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE IF(SSTTYP.EQ.'EH5B1') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_RF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_PS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='B1/'//YR_A2(NYEAR-2001)//'/'//'EH_B1'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='B1/'//YR_A2(NYEAR-2001)//'/'//'EH_PS'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE IF(SSTTYP.EQ.'EHA1B') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_RF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EH_PS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='A1B/'//YR_A2(NYEAR-2001)//'/'//'E_A1B'//
     &                         YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='A1B/'//YR_A2(NYEAR-2001)//'/'//'EH_PS'//
     &                         YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
          IF(SSTTYP.EQ.'EH5RF') THEN
            FINM='RF/'//YR_RF(NYEAR-1940)//'/'//'EHgRF'//
     &                  YR_RF(NYEAR-1940)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='RF/'//YR_RF(NYEAR-1940)//'/'//'EHgPS'//
     &                  YR_RF(NYEAR-1940)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EH5A2') THEN
            FINM='A2/'//YR_A2(NYEAR-2000)//'/'//'EHgA2'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='A2/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EH5B1') THEN
            FINM='B1/'//YR_A2(NYEAR-2000)//'/'//'EHgB1'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='B1/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                  YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE IF(SSTTYP.EQ.'EHA1B') THEN
            FINM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EgA1B'//
     &                   YR_A2(NYEAR-2000)//CHmon(MONTH)
            IF(EHSO4)
     &      PSNM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                   YR_A2(NYEAR-2000)//CHmon(MONTH)
          ELSE
            WRITE(*,*)'ERROR IN SSTTYP'
            STOP
          ENDIF
        ELSE
          IF(.NOT.(MONTH.EQ.1)) THEN
            IF(SSTTYP.EQ.'EH5RF') THEN
               FINM='RF/'//YR_RF(NYEAR-1940)//'/'//'EHgRF'//
     &                     YR_RF(NYEAR-1940)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='RF/'//YR_RF(NYEAR-1940)//'/'//'EHgPS'//
     &                     YR_RF(NYEAR-1940)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EH5A2') THEN
               FINM='A2/'//YR_A2(NYEAR-2000)//'/'//'EHgA2'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='A2/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EH5B1') THEN
               FINM='B1/'//YR_A2(NYEAR-2000)//'/'//'EHgB1'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='B1/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                     YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'EHA1B') THEN
               FINM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EgA1B'//
     &                      YR_A2(NYEAR-2000)//CHmon(MONTH-1)
               IF(EHSO4)
     &         PSNM='A1B/'//YR_A2(NYEAR-2000)//'/'//'EHgPS'//
     &                      YR_A2(NYEAR-2000)//CHmon(MONTH-1)
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
          ELSE
            IF(SSTTYP.EQ.'EH5RF') THEN
               FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgRF'//
     &                     YR_RF(NYEAR-1941)//CHmon(12)
               IF(EHSO4)
     &         PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgPS'//
     &                     YR_RF(NYEAR-1941)//CHmon(12)
            ELSE IF(SSTTYP.EQ.'EH5A2') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgRF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgPS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='A2/'//YR_A2(NYEAR-2001)//'/'//'EHgA2'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='A2/'//YR_A2(NYEAR-2001)//'/'//'EHgPS'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE IF(SSTTYP.EQ.'EH5B1') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgRF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgPS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='B1/'//YR_A2(NYEAR-2001)//'/'//'EHgB1'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='B1/'//YR_A2(NYEAR-2001)//'/'//'EHgPS'//
     &                        YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE IF(SSTTYP.EQ.'EHA1B') THEN
               IF(NYEAR.eq.2001) THEN
                  FINM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgRF'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='RF/'//YR_RF(NYEAR-1941)//'/'//'EHgPS'//
     &                        YR_RF(NYEAR-1941)//CHmon(12)
               ELSE
                  FINM='A1B/'//YR_A2(NYEAR-2001)//'/'//'EgA1B'//
     &                         YR_A2(NYEAR-2001)//CHmon(12)
                  IF(EHSO4)
     &            PSNM='A1B/'//YR_A2(NYEAR-2001)//'/'//'EHgPS'//
     &                         YR_A2(NYEAR-2001)//CHmon(12)
               ENDIF
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      DO K=1,klev*3
      DO J=1,jlat
      DO I=1,ilon
         B2(I,J,K) = -9999.
      ENDDO
      ENDDO
      ENDDO
      DO K=1,klev*2
      DO J=1,jlat
      DO I=1,ilon
         D2(I,J,K) = -9999.
      ENDDO
      ENDDO
      ENDDO
      IF(EHSO4) THEN
         DO K=1,mlev
         DO J=1,jlat
         DO I=1,ilon
            sulfate2(I,J,K) = -9999.
         ENDDO
         ENDDO
         ENDDO
         DO J=1,jlat
         DO I=1,ilon
            PSO4_2(I,J) = -9999.
         ENDDO
         ENDDO
      ENDIF

      inquire(file='../DATA/EH5OM/'//FINM,exist=there)
      if(.not.there) then
         write(*,*) '../DATA/EH5OM/'//FINM,' is not available'
         write(*,*) 'please copy EH5OM output under ../DATA/EH5OM/'
         stop
      endif
      OPEN(63,file='../DATA/EH5OM/'//FINM,form='unformatted'
     &       ,recl=(numx*numy*2+16)/4*ibyte,access='direct')
      IF(EHSO4) THEN
         inquire(file='../DATA/EH5OM/'//PSNM,exist=there)
         if(.not.there) then
            write(*,*) '../DATA/EH5OM/'//PSNM,' is not available'
            write(*,*) 'please copy EH5OM output under ../DATA/EH5OM/'
            stop
         endif
         OPEN(62,file='../DATA/EH5OM/'//PSNM,form='unformatted'
     &          ,recl=(numx*numy*2+16)/4*ibyte,access='direct')
      ENDIF
      IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
         nrec=((NDAY-1)*4+NHOUR/6-1)*(klev*5)
         IF(EHSO4) krec= (NDAY-1)*4+NHOUR/6
      ELSE
         IF(MONTH.EQ.1.OR.MONTH.EQ.2.OR.MONTH.EQ.4.OR.
     &      MONTH.EQ.6.OR.MONTH.EQ.8.OR.MONTH.EQ.9.OR.
     &      MONTH.EQ.11) THEN
            nrec=(31*4-1)*(klev*5)
            IF(EHSO4) krec= 31*4
         ELSE IF(MONTH.EQ.5.OR.MONTH.EQ.7.OR.MONTH.EQ.10
     &                     .OR.MONTH.EQ.12) THEN
            nrec=(30*4-1)*(klev*5)
            IF(EHSO4) krec= 30*4
         ELSE
            IF(MOD(NYEAR,4).EQ.0.AND.NYEAR.NE.2100) THEN
               nrec=(29*4-1)*(klev*5)
               IF(EHSO4) krec= 29*4
            ELSE
               nrec=(28*4-1)*(klev*5)
               IF(EHSO4) krec= 28*4
            ENDIF
         ENDIF
      ENDIF
      IF(EHSO4) THEN
         read(62,rec=krec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              PSO4_2(ii,49-j)=itmp(i2,j2)*scale+offset+PSO4_0(ii,49-j)
            else
           PSO4_2(ii,j+96/2)=itmp(i2,j2)*scale+offset+PSO4_0(ii,j+96/2)
            endif
         enddo
         enddo
         close(62)
      ENDIF
         
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              Hvar(ii,49-j,k)=itmp(i2,j2)*scale+offset
            else
              Hvar(ii,j+96/2,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              RHvar(ii,49-j,k)=
     &        dmin1(dmax1(itmp(i2,j2)*scale+offset,0.d0),1.d0)
            else
              RHvar(ii,j+96/2,k)=
     &        dmin1(dmax1(itmp(i2,j2)*scale+offset,0.d0),1.d0)
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              Tvar(ii,49-j,k)=itmp(i2,j2)*scale+offset
            else
              Tvar(ii,j+96/2,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              Uvar(ii,49-j,k)=itmp(i2,j2)*scale+offset
            else
              Uvar(ii,j+96/2,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint((lat0+.9375)/1.875),nint((lat1+.9375)/1.875)
         do i=nint(lon0/1.875),nint(lon1/1.875)
            ii=i+1
            if(ii.le.0) ii=ii+192
            if(ii.gt.192) ii=ii-192
            i2=i-nint(lon0/1.875)+1
            j2=j-nint((lat0+.9375)/1.875)+1
            if(numx.eq.192.and.numy.eq.96) then
              Vvar(ii,49-j,k)=itmp(i2,j2)*scale+offset
            else
              Vvar(ii,j+96/2,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      close(63)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL MKSST3(TS4,SST1,TOPOGM,XLANDU,JX,IY,IDATE)
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
      IF(EHSO4) THEN
         IF(SSTTYP.EQ.'EH5RF') THEN
            fnso4_RF='../DATA/EH5OM/SO4/RF/T63L31_skg_'//
     &                YR_RF(NYEAR-1940)//'.nc'
            inquire(file=fnso4_RF,exist=there)
            if(.not.there) then
               write(*,*) fnso4_RF,' is not available'
               stop
            endif
            status=nf_open(fnso4_RF,nf_nowrite,ncid)
            status=nf_get_vara_real(ncid,10,start,count,sulfate)
            status=nf_close(ncid)
            IF(NYEAR.eq.1950.and.MONTH.eq.1.and.NDAY.lt.16) THEN
               do k=1,mlev
               do j=1,jlat
               do i=1,ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,1)
               enddo
               enddo
               enddo
            ELSE IF(NYEAR.eq.2000.and.MONTH.eq.12.and.NDAY.ge.16) THEN
               do k=1,mlev
               do j=1,jlat
               do i=1,ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,12)
               enddo
               enddo
               enddo
            ELSE
               IF((MONTH.eq.1.or.MONTH.eq.3..or.MONTH.eq.5.or.MONTH.eq.7
     &         .or.MONTH.eq.8.or.MONTH.eq.10).and.NDAY.GE.16) THEN
                  IF(MONTH.eq.1) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/30.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/30.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.3.or.MONTH.eq.5.or.
     &                    MONTH.eq.8.or.MONTH.eq.10) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/31.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/31.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.7) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/32.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/32.)
                     enddo
                     enddo
                     enddo
                  ENDIF
               ELSE IF(MONTH.eq.2.and.NDAY.GE.15) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/30.)
     &             + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/30.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.
     &                  MONTH.eq.11).and.NDAY.GE.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/31.)
     &             + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/31.)
                  enddo
                  enddo
                  enddo
               ELSE IF(MONTH.eq.12..and.NDAY.GE.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                  enddo
                  enddo
                  enddo
                  fnso4_RF='../DATA/EH5OM/SO4/RF/T63L31_skg_'//
     &                      YR_RF(NYEAR-1939)//'.nc'
                  inquire(file=fnso4_RF,exist=there)
                  if(.not.there) then
                     write(*,*) fnso4_RF,' is not available'
                     stop
                  endif
                  status=nf_open(fnso4_RF,nf_nowrite,ncid)
                  status=nf_get_vara_real(ncid,10,start,count,sulfate)
                  status=nf_close(ncid)
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                  enddo
                  enddo
                  enddo
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate1(i,j,k,1)*(1.-float(NDAY-16)/32.)
     &             + sulfate1(i,j,k,2)*(float(NDAY-16)/32.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.3.or.MONTH.eq.5.or.MONTH.eq.7.or.
     &                  MONTH.eq.8.or.MONTH.eq.10.or.MONTH.eq.12)
     &                 .and.NDAY.LT.16) THEN
                  IF(MONTH.eq.3) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/30.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/30.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.5.or.MONTH.eq.7.or.MONTH.eq.10.or.
     &                    MONTH.eq.12) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/31.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/31.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.8) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/32.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/32.)
                     enddo
                     enddo
                     enddo
                  ENDIF
               ELSE IF(MONTH.eq.2.and.NDAY.LT.15) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH-1)*(float(15-NDAY)/30.)
     &             + sulfate(i,j,k,MONTH)*(1.-float(15-NDAY)/30.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.
     &                  MONTH.eq.11).and.NDAY.LT.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/31.)
     &             + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/31.)
                  enddo
                  enddo
                  enddo
               ELSE IF(MONTH.eq.1.and.NDAY.LT.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                  enddo
                  enddo
                  enddo
                  fnso4_RF='../DATA/EH5OM/SO4/RF/T63L31_skg_'//
     &                      YR_RF(NYEAR-1941)//'.nc'
                  inquire(file=fnso4_RF,exist=there)
                  if(.not.there) then
                     write(*,*) fnso4_RF,' is not available'
                     stop
                  endif
                  status=nf_open(fnso4_RF,nf_nowrite,ncid)
                  status=nf_get_vara_real(ncid,10,start,count,sulfate)
                  status=nf_close(ncid)
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                  enddo
                  enddo
                  enddo
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate1(i,j,k,1)*(float(16-NDAY)/32.)
     &             + sulfate1(i,j,k,2)*(1.-float(16-NDAY)/32.)
                  enddo
                  enddo
                  enddo
               ENDIF
            ENDIF
         ELSE
            IF(SSTTYP.EQ.'EH5A2') THEN
               fnso4_A2='../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//
     &                   YR_A2(NYEAR-2000)//'.nc'
               inquire(file=fnso4_A2,exist=there)
               if(.not.there) then
                  write(*,*) fnso4_A2,' is not available'
                  stop
               endif
               status=nf_open(fnso4_A2,nf_nowrite,ncid)
            ELSE IF(SSTTYP.EQ.'EHA1B') THEN
               fnso4_A1B='../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//
     &                   YR_A2(NYEAR-2000)//'.nc'
               inquire(file=fnso4_A1B,exist=there)
               if(.not.there) then
                  write(*,*) fnso4_A1B,' is not available'
                  stop
               endif
               status=nf_open(fnso4_A1B,nf_nowrite,ncid)
            ELSE IF(SSTTYP.EQ.'EH5B1') THEN
               fnso4_B1='../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//
     &                   YR_A2(NYEAR-2000)//'.nc'
               inquire(file=fnso4_B1,exist=there)
               if(.not.there) then
                  write(*,*) fnso4_B1,' is not available'
                  stop
               endif
               status=nf_open(fnso4_B1,nf_nowrite,ncid)
            ENDIF
            status=nf_get_vara_real(ncid,10,start,count,sulfate)
            status=nf_close(ncid)
            IF(NYEAR.eq.2001.and.MONTH.eq.1.and.NDAY.lt.16) THEN
               do k=1,mlev
               do j=1,jlat
               do i=1,ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,1)
               enddo
               enddo
               enddo
            ELSE IF(NYEAR.eq.2100.and.MONTH.eq.12.and.NDAY.ge.16) THEN
               do k=1,mlev
               do j=1,jlat
               do i=1,ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,12)
               enddo
               enddo
               enddo
            ELSE
               IF((MONTH.eq.1.or.MONTH.eq.3..or.MONTH.eq.5.or.MONTH.eq.7
     &         .or.MONTH.eq.8.or.MONTH.eq.10).and.NDAY.GE.16) THEN
                  IF(MONTH.eq.1) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/30.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/30.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.3.or.MONTH.eq.5.or.
     &                    MONTH.eq.8.or.MONTH.eq.10) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/31.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/31.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.7) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/32.)
     &                + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/32.)
                     enddo
                     enddo
                     enddo
                  ENDIF
               ELSE IF(MONTH.eq.2.and.NDAY.GE.15) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/30.)
     &             + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/30.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.
     &                  MONTH.eq.11).and.NDAY.GE.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH)*(1.-float(NDAY-16)/31.)
     &             + sulfate(i,j,k,MONTH+1)*(float(NDAY-16)/31.)
                  enddo
                  enddo
                  enddo
               ELSE IF(MONTH.eq.12..and.NDAY.GE.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                  enddo
                  enddo
                  enddo
                  IF(SSTTYP.EQ.'EH5A2') THEN
                     fnso4_A2='../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//
     &                         YR_A2(NYEAR-1999)//'.nc'
                     inquire(file=fnso4_A2,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_A2,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_A2,nf_nowrite,ncid)
                  ELSE IF(SSTTYP.EQ.'EHA1B') THEN
                     fnso4_A1B='../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//
     &                         YR_A2(NYEAR-1999)//'.nc'
                     inquire(file=fnso4_A1B,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_A1B,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_A1B,nf_nowrite,ncid)
                  ELSE IF(SSTTYP.EQ.'EH5B1') THEN
                     fnso4_B1='../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//
     &                         YR_A2(NYEAR-1999)//'.nc'
                     inquire(file=fnso4_B1,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_B1,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_B1,nf_nowrite,ncid)
                  ENDIF
                  status=nf_get_vara_real(ncid,10,start,count,sulfate)
                  status=nf_close(ncid)
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                  enddo
                  enddo
                  enddo
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate1(i,j,k,1)*(1.-float(NDAY-16)/32.)
     &             + sulfate1(i,j,k,2)*(float(NDAY-16)/32.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.3.or.MONTH.eq.5.or.MONTH.eq.7.or.
     &                  MONTH.eq.8.or.MONTH.eq.10.or.MONTH.eq.12)
     &                 .and.NDAY.LT.16) THEN
                  IF(MONTH.eq.3) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/30.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/30.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.5.or.MONTH.eq.7.or.MONTH.eq.10.or.
     &                    MONTH.eq.12) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/31.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/31.)
                     enddo
                     enddo
                     enddo
                  ELSE IF(MONTH.eq.8) THEN
                     do k=1,mlev
                     do j=1,jlat
                     do i=1,ilon
                        sulfate2(i,j,k) =
     &                  sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/32.)
     &                + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/32.)
                     enddo
                     enddo
                     enddo
                  ENDIF
               ELSE IF(MONTH.eq.2.and.NDAY.LT.15) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH-1)*(float(15-NDAY)/30.)
     &             + sulfate(i,j,k,MONTH)*(1.-float(15-NDAY)/30.)
                  enddo
                  enddo
                  enddo
               ELSE IF((MONTH.eq.4..or.MONTH.eq.6.or.MONTH.eq.9.or.
     &                  MONTH.eq.11).and.NDAY.LT.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate(i,j,k,MONTH-1)*(float(16-NDAY)/31.)
     &             + sulfate(i,j,k,MONTH)*(1.-float(16-NDAY)/31.)
                  enddo
                  enddo
                  enddo
               ELSE IF(MONTH.eq.1.and.NDAY.LT.16) THEN
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                  enddo
                  enddo
                  enddo
                  IF(SSTTYP.EQ.'EH5A2') THEN
                     fnso4_A2='../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//
     &                         YR_A2(NYEAR-2001)//'.nc'
                     inquire(file=fnso4_A2,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_A2,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_A2,nf_nowrite,ncid)
                  ELSE IF(SSTTYP.EQ.'EHA1B') THEN
                     fnso4_A1B='../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//
     &                         YR_A2(NYEAR-2001)//'.nc'
                     inquire(file=fnso4_A1B,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_A1B,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_A1B,nf_nowrite,ncid)
                  ELSE IF(SSTTYP.EQ.'EH5B1') THEN
                     fnso4_B1='../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//
     &                         YR_A2(NYEAR-2001)//'.nc'
                     inquire(file=fnso4_B1,exist=there)
                     if(.not.there) then
                        write(*,*) fnso4_B1,' is not available'
                        stop
                     endif
                     status=nf_open(fnso4_B1,nf_nowrite,ncid)
                  ENDIF
                  status=nf_get_vara_real(ncid,10,start,count,sulfate)
                  status=nf_close(ncid)
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                  enddo
                  enddo
                  enddo
                  do k=1,mlev
                  do j=1,jlat
                  do i=1,ilon
                     sulfate2(i,j,k) =
     &               sulfate1(i,j,k,1)*(float(16-NDAY)/32.)
     &             + sulfate1(i,j,k,2)*(1.-float(16-NDAY)/32.)
                  enddo
                  enddo
                  enddo
               ENDIF
            ENDIF
         ENDIF
         CALL BILINX(sulfate3,sulfate2,XLON,XLAT,GLON,GLAT,
     &               ilon,jlat,JX,IY,mlev)
         CALL BILINX(PSO4_3,PSO4_2,XLON,XLAT,GLON,GLAT,
     &               ilon,jlat,JX,IY,1)
         do i=1,iy
         do j=1,jx
            do l=1,KZ
               prcm=((PS4(j,i)-PTOP)*SIGMA2(l)+PTOP)*10.
               k0 = -1
               do k=mlev,1,-1
                  pmpi= (PSO4_3(j,i)*hybm(k)+hyam(k))*0.01
                  k0=k
                  if(prcm.gt.pmpi) goto 1001
               enddo
 1001          continue
               if(k0.eq.mlev) then
                  pmpj=(PSO4_3(j,i)*hybm(mlev-1)+hyam(mlev-1))*0.01
                  pmpi=(PSO4_3(j,i)*hybm(mlev)+hyam(mlev))*0.01
                  sulfate4(j,i,l)=sulfate3(j,i,mlev)
     &          +(sulfate3(j,i,mlev)-sulfate3(j,i,mlev-1))
     &                          *(prcm-pmpi)/(pmpi-pmpj)
               else if(k0.ge.1) then
                  pmpj=(PSO4_3(j,i)*hybm(k0)+hyam(k0))*0.01
                  pmpi=(PSO4_3(j,i)*hybm(k0+1)+hyam(k0+1))*0.01
                  sulfate4(j,i,l)=(sulfate3(j,i,k0+1)*(prcm-pmpj)
     &                            +sulfate3(j,i,k0)*(pmpi-prcm))
     &                           /(pmpi-pmpj)
               endif
            enddo
         enddo
         enddo
!
         DO K=KZ,1,-1
            NOUTREC=NOUTREC+1
            WRITE(64,rec=NOUTREC) ((sulfate4(j,i,K),j=1,jx),i=1,IY)
         ENDDO
      ENDIF
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADERMPI(EHSO4)
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      LOGICAL EHSO4
      INTEGER klev,jlat,ilon
      PARAMETER(klev=17,jlat=96,ilon=192)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEH/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)

      integer mlev
      parameter (mlev=31)
      REAL    hyai,hybi,hyam,hybm
      integer(kind=4) start,count
      real(kind=4)  PSO4_0
      COMMON /GLOBALEH2/hyai(mlev+1),hybi(mlev+1),hyam(mlev),hybm(mlev),
     &                  PSO4_0(ilon,jlat),start(10),count(10)
      logical there
      INTEGER I,K,KR
!
      GLAT( 1)= -88.5719985961914
      GLAT( 2)= -86.7229995727539
      GLAT( 3)= -84.8619995117188
      GLAT( 4)= -82.9990005493164
      GLAT( 5)= -81.1350021362305
      GLAT( 6)= -79.2710037231445
      GLAT( 7)= -77.4059982299805
      GLAT( 8)= -75.5410003662109
      GLAT( 9)= -73.6760025024414
      GLAT(10)= -71.8109970092773
      GLAT(11)= -69.9459991455078
      GLAT(12)= -68.0810012817383
      GLAT(13)= -66.2160034179688
      GLAT(14)= -64.3509979248047
      GLAT(15)= -62.4860000610352
      GLAT(16)= -60.6199989318848
      GLAT(17)= -58.7550010681152
      GLAT(18)= -56.8899993896484
      GLAT(19)= -55.0250015258789
      GLAT(20)= -53.1599998474121
      GLAT(21)= -51.2939987182617
      GLAT(22)= -49.4290008544922
      GLAT(23)= -47.5639991760254
      GLAT(24)= -45.6990013122559
      GLAT(25)= -43.8330001831055
      GLAT(26)= -41.9679985046387
      GLAT(27)= -40.1030006408691
      GLAT(28)= -38.2379989624023
      GLAT(29)= -36.3720016479492
      GLAT(30)= -34.5069999694824
      GLAT(31)= -32.6419982910156
      GLAT(32)= -30.7770004272461
      GLAT(33)= -28.9109992980957
      GLAT(34)= -27.0459995269775
      GLAT(35)= -25.1809997558594
      GLAT(36)= -23.3159999847412
      GLAT(37)= -21.4500007629395
      GLAT(38)= -19.5849990844727
      GLAT(39)= -17.7199993133545
      GLAT(40)= -15.8549995422363
      GLAT(41)= -13.9890003204346
      GLAT(42)= -12.1239995956421
      GLAT(43)= -10.2589998245239
      GLAT(44)= -8.39400005340576
      GLAT(45)= -6.52799987792969
      GLAT(46)= -4.66300010681152
      GLAT(47)= -2.79800009727478
      GLAT(48)= -0.933000028133392
      GLAT(49)=  0.933000028133392
      GLAT(50)=  2.79800009727478
      GLAT(51)=  4.66300010681152
      GLAT(52)=  6.52799987792969
      GLAT(53)=  8.39400005340576
      GLAT(54)= 10.2589998245239
      GLAT(55)= 12.1239995956421
      GLAT(56)= 13.9890003204346
      GLAT(57)= 15.8549995422363
      GLAT(58)= 17.7199993133545
      GLAT(59)= 19.5849990844727
      GLAT(60)= 21.4500007629395
      GLAT(61)= 23.3159999847412
      GLAT(62)= 25.1809997558594
      GLAT(63)= 27.0459995269775
      GLAT(64)= 28.9109992980957
      GLAT(65)= 30.7770004272461
      GLAT(66)= 32.6419982910156
      GLAT(67)= 34.5069999694824
      GLAT(68)= 36.3720016479492
      GLAT(69)= 38.2379989624023
      GLAT(70)= 40.1030006408691
      GLAT(71)= 41.9679985046387
      GLAT(72)= 43.8330001831055
      GLAT(73)= 45.6990013122559
      GLAT(74)= 47.5639991760254
      GLAT(75)= 49.4290008544922
      GLAT(76)= 51.2939987182617
      GLAT(77)= 53.1599998474121
      GLAT(78)= 55.0250015258789
      GLAT(79)= 56.8899993896484
      GLAT(80)= 58.7550010681152
      GLAT(81)= 60.6199989318848
      GLAT(82)= 62.4860000610352
      GLAT(83)= 64.3509979248047
      GLAT(84)= 66.2160034179688
      GLAT(85)= 68.0810012817383
      GLAT(86)= 69.9459991455078
      GLAT(87)= 71.8109970092773
      GLAT(88)= 73.6760025024414
      GLAT(89)= 75.5410003662109
      GLAT(90)= 77.4059982299805
      GLAT(91)= 79.2710037231445
      GLAT(92)= 81.1350021362305
      GLAT(93)= 82.9990005493164
      GLAT(94)= 84.8619995117188
      GLAT(95)= 86.7229995727539
      GLAT(96)= 88.5719985961914
!
      SIGMAR(1) = .01
      SIGMAR(2) = .03
      SIGMAR(3) = .05
      SIGMAR(4) = .07
      SIGMAR(5) = .1
      SIGMAR(6) = .15
      SIGMAR(7) = .2
      SIGMAR(8) = .25
      SIGMAR(9) = .3
      SIGMAR(10)= .4
      SIGMAR(11)= .5
      SIGMAR(12)= .6
      SIGMAR(13)= .7
      SIGMAR(14)= .775
      SIGMAR(15)= .85
      SIGMAR(16)= .925
      SIGMAR(17)=1.0
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*1.875
      ENDDO
!     DO J = 1,jlat
!        GLAT(J) = -89.0625+FLOAT(J-1)*1.875
!     ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      IF(EHSO4) THEN
         hyai( 1) =     0.
         hyai( 2) =  2000.
         hyai( 3) =  4000.
         hyai( 4) =  6000.
         hyai( 5) =  8000.
         hyai( 6) =  9976.13671875
         hyai( 7) = 11820.5400390625
         hyai( 8) = 13431.3896484375
         hyai( 9) = 14736.3603515625
         hyai(10) = 15689.2099609375
         hyai(11) = 16266.6103515625
         hyai(12) = 16465.
         hyai(13) = 16297.6201171875
         hyai(14) = 15791.599609375
         hyai(15) = 14985.26953125
         hyai(16) = 13925.51953125
         hyai(17) = 12665.2900390625
         hyai(18) = 11261.23046875
         hyai(19) =  9771.40625
         hyai(20) =  8253.2109375
         hyai(21) =  6761.33984375
         hyai(22) =  5345.9140625
         hyai(23) =  4050.71801757812
         hyai(24) =  2911.56909179688
         hyai(25) =  1954.80505371094
         hyai(26) =  1195.89001464844
         hyai(27) =   638.14892578125
         hyai(28) =   271.626495361328
         hyai(29) =    72.0635833740234
         hyai(30) =     0.
         hyai(31) =     0.
         hyai(32) =     0.
   
         hybi( 1) =  0.
         hybi( 2) =  0.
         hybi( 3) =  0.
         hybi( 4) =  0.
         hybi( 5) =  0.
         hybi( 6) =  0.000390858185710385
         hybi( 7) =  0.0029197009280324
         hybi( 8) =  0.00919413194060326
         hybi( 9) =  0.0203191600739956
         hybi(10) =  0.0369748584926128
         hybi(11) =  0.0594876408576965
         hybi(12) =  0.0878949835896492
         hybi(13) =  0.122003600001335
         hybi(14) =  0.161441504955292
         hybi(15) =  0.205703303217888
         hybi(16) =  0.254188597202301
         hybi(17) =  0.306235402822495
         hybi(18) =  0.361144989728928
         hybi(19) =  0.418202310800552
         hybi(20) =  0.476688086986542
         hybi(21) =  0.535886585712433
         hybi(22) =  0.595084190368652
         hybi(23) =  0.65356457233429
         hybi(24) =  0.710594415664673
         hybi(25) =  0.765405178070068
         hybi(26) =  0.817166984081268
         hybi(27) =  0.86495578289032
         hybi(28) =  0.907715916633606
         hybi(29) =  0.944213211536407
         hybi(30) =  0.972985208034515
         hybi(31) =  0.992281496524811
         hybi(32) =  1.0
   
         hyam( 1) =  1000.
         hyam( 2) =  3000.
         hyam( 3) =  5000.
         hyam( 4) =  7000.
         hyam( 5) =  8988.068359375
         hyam( 6) = 10898.3383789062
         hyam( 7) = 12625.96484375
         hyam( 8) = 14083.875
         hyam( 9) = 15212.78515625
         hyam(10) = 15977.91015625
         hyam(11) = 16365.8051757812
         hyam(12) = 16381.3100585938
         hyam(13) = 16044.6098632812
         hyam(14) = 15388.4345703125
         hyam(15) = 14455.39453125
         hyam(16) = 13295.4047851562
         hyam(17) = 11963.2602539062
         hyam(18) = 10516.318359375
         hyam(19) =  9012.30859375
         hyam(20) =  7507.275390625
         hyam(21) =  6053.626953125
         hyam(22) =  4698.31604003906
         hyam(23) =  3481.1435546875
         hyam(24) =  2433.18707275391
         hyam(25) =  1575.34753417969
         hyam(26) =   917.019470214844
         hyam(27) =   454.887710571289
         hyam(28) =   171.845039367676
         hyam(29) =    36.0317916870117
         hyam(30) =     0.
         hyam(31) =     0.
   
         hybm( 1) = 0.
         hybm( 2) = 0.
         hybm( 3) = 0.
         hybm( 4) = 0.
         hybm( 5) = 0.000195429092855193
         hybm( 6) = 0.00165527955687139
         hybm( 7) = 0.00605691643431783
         hybm( 8) = 0.0147566460072994
         hybm( 9) = 0.0286470092833042
         hybm(10) = 0.0482312496751547
         hybm(11) = 0.0736913122236729
         hybm(12) = 0.104949291795492
         hybm(13) = 0.141722552478313
         hybm(14) = 0.18357240408659
         hybm(15) = 0.229945950210094
         hybm(16) = 0.280212000012398
         hybm(17) = 0.333690196275711
         hybm(18) = 0.38967365026474
         hybm(19) = 0.447445198893547
         hybm(20) = 0.506287336349487
         hybm(21) = 0.565485388040543
         hybm(22) = 0.624324381351471
         hybm(23) = 0.682079493999481
         hybm(24) = 0.737999796867371
         hybm(25) = 0.791286081075668
         hybm(26) = 0.841061383485794
         hybm(27) = 0.886335849761963
         hybm(28) = 0.925964564085007
         hybm(29) = 0.958599209785461
         hybm(30) = 0.982633352279663
         hybm(31) = 0.996140748262405
   
         start(1) = 1
         start(2) = 1
         start(3) = 1
         start(4) = 1
         start(5) = 0
         start(6) = 0
         start(7) = 0
         start(8) = 0
         start(9) = 0
         start(10)= 0
   
         count(1) = 192
         count(2) = 96
         count(3) = 31
         count(4) = 12
         count(5) = 0
         count(6) = 0
         count(7) = 0
         count(8) = 0
         count(9) = 0
         count(10)= 0
         
         inquire(file='../DATA/EH5OM/EHgPS.dat',exist=there)
         if(.not.there) then
           write(*,*) '../DATA/EH5OM/EHgPS.dat is not available'
           stop
         endif
         open(30,file='../DATA/EH5OM/EHgPS.dat',form='unformatted'
     &          ,recl=ilon*jlat*4,access='direct')
         read(30,rec=1) PSO4_0
         close(30)
      ENDIF

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETERA40(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,Qsoil,TSIce,TSoil,SNOWh
      common /ERAvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        QSoil(ilon,jlat,4)   ,TSIce(ilon,jlat,4)
     &       ,        TSoil(ilon,jlat,4)   ,SNOWh(ilon,jlat)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    S2(ilon,jlat,4*3+1)
      EQUIVALENCE (S2(1,1,1),QSoil(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALERA/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3,QS3,TI3,TS3,SNOW
      COMMON /ERAVAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,         T3(JX,IY,klev),H3(JX,IY,klev),Q3(JX,IY,klev)
     &       ,  QS3(JX,IY,4),TI3(JX,IY,4),TS3(JX,IY,4),SNOW(JX,IY)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
      REAL    S3(JX,IY,4*3+1)
      EQUIVALENCE (S3(1,1,1),QS3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER ERA40 REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL ERA6HOUR(DATTYP,LSMTYP,IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
      IF(LSMTYP.EQ.'USGS') THEN
       CALL BILINX(S3,S2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,4*3+1)
      ENDIF
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEFS(U4,V4,T4,Q4,PS4,TS4,QS3,TI3,TS3,SNOW
     &            ,PTOP,JX,IY,KZ,IDATE,LSMTYP)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADERERA
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALERA/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .001
      SIGMAR(2) = .002
      SIGMAR(3) = .003
      SIGMAR(4) = .005
      SIGMAR(5) = .007
      SIGMAR(6) = .01
      SIGMAR(7) = .02
      SIGMAR(8) = .03
      SIGMAR(9) = .05
      SIGMAR(10)= .07
      SIGMAR(11)= .1
      SIGMAR(12)= .15
      SIGMAR(13)= .2
      SIGMAR(14)= .25
      SIGMAR(15)= .3
      SIGMAR(16)= .4
      SIGMAR(17)= .5
      SIGMAR(18)= .6
      SIGMAR(19)= .7
      SIGMAR(20)= .775
      SIGMAR(21)= .85
      SIGMAR(22)= .925
      SIGMAR(23)=1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*2.5
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*2.5
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETEIN75(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=241,ilon=480)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN75vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,          RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,           Vvar(ilon,jlat,klev)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN75/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /EIN75VAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,           T3(JX,IY,klev),H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER EIN75 REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL EIN756HOUR(DATTYP,IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!     CALL HUMID1(T3,Q3,100.,0.0,SIGMA1,JX,IY,klev)
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADEREIN75
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=241,ilon=480)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN75/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .001
      SIGMAR(2) = .002
      SIGMAR(3) = .003
      SIGMAR(4) = .005
      SIGMAR(5) = .007
      SIGMAR(6) = .01
      SIGMAR(7) = .02
      SIGMAR(8) = .03
      SIGMAR(9) = .05
      SIGMAR(10)= .07
      SIGMAR(11)= .1
      SIGMAR(12)= .15
      SIGMAR(13)= .2
      SIGMAR(14)= .25
      SIGMAR(15)= .3
      SIGMAR(16)= .4
      SIGMAR(17)= .5
      SIGMAR(18)= .6
      SIGMAR(19)= .7
      SIGMAR(20)= .775
      SIGMAR(21)= .85
      SIGMAR(22)= .925
      SIGMAR(23)=1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*0.75
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*0.75
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETEIN15(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=121,ilon=240)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN15vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,          RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,           Vvar(ilon,jlat,klev)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN15/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /EIN15VAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,           T3(JX,IY,klev),H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER EIN15 REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL EIN156HOUR(DATTYP,IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!     CALL HUMID1(T3,Q3,100.,0.0,SIGMA1,JX,IY,klev)
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.EQ.'OI_WK') THEN
         CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ELSE IF(SSTTYP.EQ.'ERSST'.or.SSTTYP.EQ.'ERSKT') THEN
         CALL MKSST3(TS4,SST1,TOPOGM,XLANDU,JX,IY,IDATE)
      ELSE
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
         PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
         CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
         CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADEREIN15
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=121,ilon=240)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN15/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .001
      SIGMAR(2) = .002
      SIGMAR(3) = .003
      SIGMAR(4) = .005
      SIGMAR(5) = .007
      SIGMAR(6) = .01
      SIGMAR(7) = .02
      SIGMAR(8) = .03
      SIGMAR(9) = .05
      SIGMAR(10)= .07
      SIGMAR(11)= .1
      SIGMAR(12)= .15
      SIGMAR(13)= .2
      SIGMAR(14)= .25
      SIGMAR(15)= .3
      SIGMAR(16)= .4
      SIGMAR(17)= .5
      SIGMAR(18)= .6
      SIGMAR(19)= .7
      SIGMAR(20)= .775
      SIGMAR(21)= .85
      SIGMAR(22)= .925
      SIGMAR(23)=1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*1.50
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*1.50
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETEIN25(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN25vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,          RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,           Vvar(ilon,jlat,klev)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN25/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3!,QS3,TI3,TS3,SNOW
      COMMON /EIN25VAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,           T3(JX,IY,klev),H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER EIN25 REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL EIN256HOUR(DATTYP,IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!     CALL HUMID1(T3,Q3,100.,0.0,SIGMA1,JX,IY,klev)
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADEREIN25
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=23,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALEIN25/ GLAT(jlat),GLON(ilon),SIGMA1(klev)
     &                    ,SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .001
      SIGMAR(2) = .002
      SIGMAR(3) = .003
      SIGMAR(4) = .005
      SIGMAR(5) = .007
      SIGMAR(6) = .01
      SIGMAR(7) = .02
      SIGMAR(8) = .03
      SIGMAR(9) = .05
      SIGMAR(10)= .07
      SIGMAR(11)= .1
      SIGMAR(12)= .15
      SIGMAR(13)= .2
      SIGMAR(14)= .25
      SIGMAR(15)= .3
      SIGMAR(16)= .4
      SIGMAR(17)= .5
      SIGMAR(18)= .6
      SIGMAR(19)= .7
      SIGMAR(20)= .775
      SIGMAR(21)= .85
      SIGMAR(22)= .925
      SIGMAR(23)=1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*2.5
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*2.5
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETGFS11(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
      INTEGER NYEAR,MONTH,NDAY,NHOUR,NREC
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP GFS FORECAST DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=26,jlat=181,ilon=360)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar
      common /GFSvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALGFS/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /GFSVAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,         T3(JX,IY,klev),H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
      CHARACTER*17 FINM
      CHARACTER*4 YRGFS(9)
      DATA    YRGFS/'2000','2001','2002','2003','2004','2005'
     &             ,'2006','2007','2008'/
      CHARACTER*3 CHmon(12)
      DATA    CHmon/ 'JAN','FEB','MAR','APR','MAY','JUN',
     &               'JUL','AUG','SEP','OCT','NOV','DEC'/

      real    lon0,lon1,lat0,lat1
      INTEGER i0,i1,j0
      COMMON /SZwindow/lon0,lon1,lat0,lat1,i0,i1,j0
      integer numx,numy,ii,i2,j2,I,J,K
!      REAL    temp(360,181)
      real(kind=8)  scale,offset
      integer(kind=2) itmp(360,181)
      logical there
 
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER GFS11 FORECAST HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-(IDATE/10000)*100
      NHOUR=MOD(IDATE,100)
      IF(IDATE.LT.2000010106) THEN
         WRITE(*,*) 'GFS 1x1 datasets is just avaiable from 2000010106'
         STOP
      ENDIF
!     IF(IDATE.GT.2008060100) THEN
!        WRITE(*,*) 'GFS 1x1 datasets is just avaiable to 2008060100'
!        STOP
!     ENDIF
      numx=nint((lon1-lon0)/1.0)+1
      numy=nint((lat1-lat0)/1.0)+1
      IF(.not.(numx.eq.360.and.numy.eq.181)) THEN
        IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
          FINM=YRGFS(NYEAR-1999)//'/'//'GFS__'//
     &         YRGFS(NYEAR-1999)//CHmon(MONTH)
        ELSE
          IF(.NOT.(MONTH.EQ.1)) THEN
            FINM=YRGFS(NYEAR-1999)//'/'//'GFS__'//
     &           YRGFS(NYEAR-1999)//CHmon(MONTH-1)
          ELSE
            FINM=YRGFS(NYEAR-2000)//'/'//'GFS__'//
     &           YRGFS(NYEAR-2000)//CHmon(12)
          ENDIF
        ENDIF
      ELSE
        IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
          FINM=YRGFS(NYEAR-1999)//'/'//'GFSg_'//
     &         YRGFS(NYEAR-1999)//CHmon(MONTH)
        ELSE
          IF(.NOT.(MONTH.EQ.1)) THEN
            FINM=YRGFS(NYEAR-1999)//'/'//'GFSg_'//
     &           YRGFS(NYEAR-1999)//CHmon(MONTH-1)
          ELSE
            FINM=YRGFS(NYEAR-2000)//'/'//'GFSg_'//
     &           YRGFS(NYEAR-2000)//CHmon(12)
          ENDIF
        ENDIF
      ENDIF
      DO K=1,klev*3
      DO J=1,jlat
      DO I=1,ilon
         B2(I,J,K) = -9999.
      ENDDO
      ENDDO
      ENDDO
      DO K=1,klev*2
      DO J=1,jlat
      DO I=1,ilon
         D2(I,J,K) = -9999.
      ENDDO
      ENDDO
      ENDDO
      inquire(file='../DATA/GFS11/'//FINM,exist=there)
      if(.not.there) then
         write(*,*) '../DATA/GFS11/'//FINM,' is not available'
         write(*,*) 'please copy GFS11 datasets under ../DATA/GFS11/'
         stop
      endif
      OPEN(63,file='../DATA/GFS11/'//FINM,form='unformatted'
     &       ,recl=(numx*numy*2+16)/4*ibyte,access='direct')
      IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
         nrec=((NDAY-1)*4+NHOUR/6-1)*127
      ELSE
         IF(MONTH.EQ.1.OR.MONTH.EQ.2.OR.MONTH.EQ.4.OR.
     &      MONTH.EQ.6.OR.MONTH.EQ.8.OR.MONTH.EQ.9.OR.
     &      MONTH.EQ.11) THEN
            nrec=(31*4-1)*127
         ELSE IF(MONTH.EQ.5.OR.MONTH.EQ.7.OR.MONTH.EQ.10.OR.
     &           MONTH.EQ.12) THEN
            nrec=(30*4-1)*127
         ELSE
            IF(MOD(NYEAR,4).EQ.0) THEN
               nrec=(29*4-1)*127
            ELSE
               nrec=(28*4-1)*127
            ENDIF
         ENDIF
      ENDIf
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0/1.0),nint(lat1/1.0)
         do i=nint(lon0/1.0),nint(lon1/1.0)
            ii=i+1
            if(ii.le.0) ii=ii+360
            if(ii.gt.360) ii=ii-360
            i2=i-nint(lon0/1.0)+1
            j2=j-nint(lat0/1.0)+1
            if(numx.eq.360.and.numy.eq.181) then
               Hvar(ii,91-j,k)=itmp(i2,j2)*scale+offset
            else
               Hvar(ii,j+90,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0/1.0),nint(lat1/1.0)
         do i=nint(lon0/1.0),nint(lon1/1.0)
            ii=i+1
            if(ii.le.0) ii=ii+360
            if(ii.gt.360) ii=ii-360
            i2=i-nint(lon0/1.0)+1
            j2=j-nint(lat0/1.0)+1
            if(numx.eq.360.and.numy.eq.181) then
               Tvar(ii,91-j,k)=itmp(i2,j2)*scale+offset
            else
               Tvar(ii,j+90,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,6,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0/1.0),nint(lat1/1.0)
         do i=nint(lon0/1.0),nint(lon1/1.0)
            ii=i+1
            if(ii.le.0) ii=ii+360
            if(ii.gt.360) ii=ii-360
            i2=i-nint(lon0/1.0)+1
            j2=j-nint(lat0/1.0)+1
            if(numx.eq.360.and.numy.eq.181) then
               RHvar(ii,91-j,k)=itmp(i2,j2)*scale+offset
            else
               RHvar(ii,j+90,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=5,1,-1
         do j=1,jlat
         do i=1,ilon
            RHvar(i,j,k) = RHvar(i,j,k+1)
         enddo
         enddo
      enddo
      do k=1,klev
         do j=1,jlat
         do i=1,ilon
            RHvar(i,j,k) = amax1(RHvar(i,j,k)*0.01,1.0)
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0/1.0),nint(lat1/1.0)
         do i=nint(lon0/1.0),nint(lon1/1.0)
            ii=i+1
            if(ii.le.0) ii=ii+360
            if(ii.gt.360) ii=ii-360
            i2=i-nint(lon0/1.0)+1
            j2=j-nint(lat0/1.0)+1
            if(numx.eq.360.and.numy.eq.181) then
               Uvar(ii,91-j,k)=itmp(i2,j2)*scale+offset
            else
               Uvar(ii,j+90,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      do k=klev,1,-1
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0/1.0),nint(lat1/1.0)
         do i=nint(lon0/1.0),nint(lon1/1.0)
            ii=i+1
            if(ii.le.0) ii=ii+360
            if(ii.gt.360) ii=ii-360
            i2=i-nint(lon0/1.0)+1
            j2=j-nint(lat0/1.0)+1
            if(numx.eq.360.and.numy.eq.181) then
               Vvar(ii,91-j,k)=itmp(i2,j2)*scale+offset
            else
               Vvar(ii,j+90,k)=itmp(i2,j2)*scale+offset
            endif
         enddo
         enddo
      enddo
      close(63)

      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADERGFS
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP GFS FORECAST DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=26,jlat=181,ilon=360)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALGFS/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .01
      SIGMAR(2) = .02
      SIGMAR(3) = .03
      SIGMAR(4) = .05
      SIGMAR(5) = .07
      SIGMAR(6) = .1
      SIGMAR(7) = .15
      SIGMAR(8) = .2
      SIGMAR(9) = .25
      SIGMAR(10)= .3
      SIGMAR(11)= .35
      SIGMAR(12)= .4
      SIGMAR(13)= .45
      SIGMAR(14)= .5
      SIGMAR(15)= .55
      SIGMAR(16)= .6
      SIGMAR(17)= .65
      SIGMAR(18)= .7
      SIGMAR(19)= .75
      SIGMAR(20)= .8
      SIGMAR(21)= .85
      SIGMAR(22)= .9
      SIGMAR(23)= .925
      SIGMAR(24)= .95
      SIGMAR(25)= .975
      SIGMAR(26)=1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*1.0
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*1.0
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETERAHI(IDATE)
      implicit none
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 : PARAMETERS FOR ERA40 original REANALYSIS DATASET   X X X
! A1
      INTEGER NLEVS,NLATS,NLONS,NLEV2
      PARAMETER (NLEVS=60,NLATS=160,NLONS=320,NLEV2=18)
!
!  NLEVS  = NUMBER OF PRESSURE LEVELS IN ERA40 original DATASET.
!  NLATS  = NUMBER OF LATITUDES  ON ERA40 original GRID.
!  NLONS  = NUMBER OF LONGITUDES ON ERA40 original GRID.
!
!  X X X X X   SET 2 : PARAMETERS FOR RegCM DATASET   X X X X
!
!HH   change the vertical levels and the model domain size.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RegCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    SLON,SLAT,SIGMA1,SIGMAR,PPLEV,AK,BK
      COMMON /EHIGRID/SLON(NLONS),SLAT(NLATS),SIGMA1(NLEV2)
     &               ,SIGMAR(NLEV2),PPLEV(NLEV2),AK(NLEVS+1),BK(NLEVS+1)
!
!     DOMAIN VARIABLES FOR RegCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)

! DIMENSION SURFACE TEMPERATURE ON FVGCM SURFACE, GIVEN BY HADAMH DATA
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!B1
      REAL    ZS2,LSM,PS2,T2,Q2,U2,V2
      COMMON /ERVARS/ ZS2(NLONS,NLATS),LSM(NLONS,NLATS),PS2(NLONS,NLATS)
     &               ,T2(NLONS,NLATS,NLEVS),Q2(NLONS,NLATS,NLEVS)
     &               ,U2(NLONS,NLATS,NLEVS),V2(NLONS,NLATS,NLEVS)
!
      REAL    Z1,HP,TP,QP,UP,VP,PP3D
      COMMON /ERVART/ Z1(NLONS,NLATS,NLEVS),TP(NLONS,NLATS,NLEV2),
     &                QP(NLONS,NLATS,NLEV2),HP(NLONS,NLATS,NLEV2),
     &                UP(NLONS,NLATS,NLEV2),VP(NLONS,NLATS,NLEV2),
     &              PP3D(NLONS,NLATS,NLEVS)
      REAL    B2(NLONS,NLATS,NLEV2*3)
      EQUIVALENCE (B2(1,1,1),TP(1,1,1))
      REAL    D2(NLONS,NLATS,NLEV2*2)
      EQUIVALENCE (D2(1,1,1),UP(1,1,1))
!
!B3
      REAL    T3,Q3,H3,U3,V3,W3,B3PD
      COMMON /ER_SB3/ T3(JX,IY,NLEV2),Q3(JX,IY,NLEV2),
     &                H3(JX,IY,NLEV2),
     &                U3(JX,IY,NLEV2),V3(JX,IY,NLEV2),
     &                W3(JX,IY,NLEV2),B3PD(JX,IY)
      REAL    B3(JX,IY,NLEV2*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,NLEV2*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!B4
      REAL    PS4,T4,Q4,H4,TS4,U4,V4
      COMMON /ER_SB4/ PS4(JX,IY),
     &                T4(JX,IY,KZ),Q4(JX,IY,KZ),
     &                H4(JX,IY,KZ),TS4(JX,IY),
     &                U4(JX,IY,KZ),V4(JX,IY,KZ)
!
      INTEGER I,J,K
      INTEGER NYRP,NMOP
      REAL    WT
!      logical there
      character*14 finame
      integer nrec!,lrec
      REAL    XLONMIN,XLONMAX,SLONMIN,SLONMAX
!
!  B1 IS FOR DATA RECORDS FROM FVGCM 1x1.25L18 GRID
!  B2 IS FOR LAT-LON GRID WITH FVGCM VERTICAL STRUCTURE
!  B3 IS FOR RegCM HORIZONTAL GRID, BUT FVGCM VERTICAL GRID
!  B4 IS FOR RegCM 3-DIMENSIONAL GRID
!
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      IF(IDATE.EQ.IDATE1) THEN
         XLONMIN= 400.
         XLONMAX=-400.
         DO J=1,IY
         DO I=1,JX
            IF(XLON(I,J).LT.XLONMIN) XLONMIN=XLON(I,J)
            IF(XLON(I,J).GT.XLONMAX) XLONMAX=XLON(I,J)
         ENDDO
         ENDDO
         WRITE(*,*) 'XLONMIN,XLONMAX= ',XLONMIN,XLONMAX
         SLONMIN= 400.
         SLONMAX=-400.
         DO I=1,NLONS
            IF(SLON(I).LT.SLONMIN) SLONMIN=SLON(I)
            IF(SLON(I).GT.SLONMAX) SLONMAX=SLON(I)
         ENDDO
         WRITE(*,*) 'SLONMIN,SLONMAX= ',SLONMIN,SLONMAX
      ENDIF
      write(finame,19) IDATE
 19   format('EHI_',I10)
      open(61,file=finame,form='unformatted'
     &       ,recl=320*160*IBYTE,access='direct')
      nrec=0
      nrec=nrec+1
      read(61,rec=nrec) ZS2
      nrec=nrec+1
      read(61,rec=nrec) LSM
      nrec=nrec+1
      read(61,rec=nrec) PS2
      do j=1,NLATS
      do i=1,NLONS
         if(LSM(i,j).lt.0.5) THEN
            ZS2(i,j)=0.000
         endif
      enddo
      enddo
      do k=1,NLEVS
         nrec=nrec+1
         read(61,rec=nrec) ((T2(i,j,k),i=1,NLONS),j=1,NLATS)
      enddo
      do k=1,NLEVS
         nrec=nrec+1
         read(61,rec=nrec) ((Q2(i,j,k),i=1,NLONS),j=1,NLATS)
      enddo
      do k=1,NLEVS
         nrec=nrec+1
         read(61,rec=nrec) ((U2(i,j,k),i=1,NLONS),j=1,NLATS)
      enddo
      do k=1,NLEVS
         nrec=nrec+1
         read(61,rec=nrec) ((V2(i,j,k),i=1,NLONS),j=1,NLATS)
      enddo

      WRITE(*,*) 'READ IN fields at DATE:',IDATE
      DO K = 1,NLEVS
      DO J = 1,NLATS
      DO I = 1,NLONS
         IF(PS2(I,J).GT.-9995.) THEN
            PP3d(I,J,K) = PS2(I,J)*0.5*(BK(K)+BK(K+1))
     &                  +          0.5*(AK(K)+AK(K+1))
         ELSE
            PP3d(I,J,K) = -9999.0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!
!     to calculate Heights on sigma surfaces.
      CALL HTSIG(T2,Z1,PP3D,PS2,ZS2,NLONS,NLATS,NLEVS)
!
!     to interpolate H,U,V,T,Q and QC
!        1. For Heights
      CALL HEIGHT(HP,Z1,T2,PS2,PP3D,ZS2,NLONS,NLATS,NLEVS,PPLEV,NLEV2)
!        2. For Zonal and Meridional Winds
      CALL INTLIN(UP,U2,PS2,PP3D,NLONS,NLATS,NLEVS,PPLEV,NLEV2)
      CALL INTLIN(VP,V2,PS2,PP3D,NLONS,NLATS,NLEVS,PPLEV,NLEV2)
!        3. For Temperatures
      CALL INTLOG(TP,T2,PS2,PP3D,NLONS,NLATS,NLEVS,PPLEV,NLEV2)
!        4. For Moisture qva & qca
      CALL HUMID1FV(T2,Q2,PP3D,NLONS,NLATS,NLEVS)
      CALL INTLIN(QP,Q2,PS2,PP3D,NLONS,NLATS,NLEVS,PPLEV,NLEV2)
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,SLON,SLAT,NLONS,NLATS,JX,IY,NLEV2*3)
      CALL BILINX(D3,D2,DLON,DLAT,SLON,SLAT,NLONS,NLATS,JX,IY,NLEV2*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,NLEV2
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,NLEV2)
      CALL TOP2BTM(Q3,JX,IY,NLEV2)
      CALL TOP2BTM(H3,JX,IY,NLEV2)
      CALL TOP2BTM(U3,JX,IY,NLEV2)
      CALL TOP2BTM(V3,JX,IY,NLEV2)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,NLEV2)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
!
! F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!       INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,NLEV2)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2     DETERMINE P* AND HEIGHT.
!
! F3     INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
      CALL HUMID2FV(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4     DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G      WRITE AN INITIAL FILE FOR THE RegCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADEREHI
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR ERA40 original REANALYSIS DATASET X X X
! A1
      INTEGER NLEVS,NLATS,NLONS,NLEV2
      PARAMETER (NLEVS=60,NLATS=160,NLONS=320,NLEV2=18)
!
!  NLEVS  = NUMBER OF PRESSURE LEVELS IN ERA40 original DATASET.
!  NLATS  = NUMBER OF LATITUDES  ON ERA40 original GRID.
!  NLONS  = NUMBER OF LONGITUDES ON ERA40 original GRID.
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    SLON,SLAT,SIGMA1,SIGMAR,PPLEV,AK,BK
      COMMON /EHIGRID/SLON(NLONS),SLAT(NLATS),SIGMA1(NLEV2)
     &               ,SIGMAR(NLEV2),PPLEV(NLEV2),AK(NLEVS+1),BK(NLEVS+1)
!
      INTEGER I,K,KR
!
      SLAT(  1)= -89.142
      SLAT(  2)= -88.029
      SLAT(  3)= -86.911
      SLAT(  4)= -85.791
      SLAT(  5)= -84.670
      SLAT(  6)= -83.549
      SLAT(  7)= -82.428
      SLAT(  8)= -81.307
      SLAT(  9)= -80.185
      SLAT( 10)= -79.064
      SLAT( 11)= -77.943
      SLAT( 12)= -76.821
      SLAT( 13)= -75.700
      SLAT( 14)= -74.578
      SLAT( 15)= -73.457
      SLAT( 16)= -72.336
      SLAT( 17)= -71.214
      SLAT( 18)= -70.093
      SLAT( 19)= -68.971
      SLAT( 20)= -67.850
      SLAT( 21)= -66.728
      SLAT( 22)= -65.607
      SLAT( 23)= -64.485
      SLAT( 24)= -63.364
      SLAT( 25)= -62.242
      SLAT( 26)= -61.121
      SLAT( 27)= -60.000
      SLAT( 28)= -58.878
      SLAT( 29)= -57.757
      SLAT( 30)= -56.635
      SLAT( 31)= -55.514
      SLAT( 32)= -54.392
      SLAT( 33)= -53.271
      SLAT( 34)= -52.149
      SLAT( 35)= -51.028
      SLAT( 36)= -49.906
      SLAT( 37)= -48.785
      SLAT( 38)= -47.663
      SLAT( 39)= -46.542
      SLAT( 40)= -45.420
      SLAT( 41)= -44.299
      SLAT( 42)= -43.177
      SLAT( 43)= -42.056
      SLAT( 44)= -40.934
      SLAT( 45)= -39.813
      SLAT( 46)= -38.691
      SLAT( 47)= -37.570
      SLAT( 48)= -36.448
      SLAT( 49)= -35.327
      SLAT( 50)= -34.205
      SLAT( 51)= -33.084
      SLAT( 52)= -31.962
      SLAT( 53)= -30.841
      SLAT( 54)= -29.719
      SLAT( 55)= -28.598
      SLAT( 56)= -27.476
      SLAT( 57)= -26.355
      SLAT( 58)= -25.234
      SLAT( 59)= -24.112
      SLAT( 60)= -22.991
      SLAT( 61)= -21.869
      SLAT( 62)= -20.748
      SLAT( 63)= -19.626
      SLAT( 64)= -18.505
      SLAT( 65)= -17.383
      SLAT( 66)= -16.262
      SLAT( 67)= -15.140
      SLAT( 68)= -14.019
      SLAT( 69)= -12.897
      SLAT( 70)= -11.776
      SLAT( 71)= -10.654
      SLAT( 72)=  -9.533
      SLAT( 73)=  -8.411
      SLAT( 74)=  -7.290
      SLAT( 75)=  -6.168
      SLAT( 76)=  -5.047
      SLAT( 77)=  -3.925
      SLAT( 78)=  -2.804
      SLAT( 79)=  -1.682
      SLAT( 80)=  -0.561
      SLAT( 81)=   0.561
      SLAT( 82)=   1.682
      SLAT( 83)=   2.804
      SLAT( 84)=   3.925
      SLAT( 85)=   5.047
      SLAT( 86)=   6.168
      SLAT( 87)=   7.290
      SLAT( 88)=   8.411
      SLAT( 89)=   9.533
      SLAT( 90)=  10.654
      SLAT( 91)=  11.776
      SLAT( 92)=  12.897
      SLAT( 93)=  14.019
      SLAT( 94)=  15.140
      SLAT( 95)=  16.262
      SLAT( 96)=  17.383
      SLAT( 97)=  18.505
      SLAT( 98)=  19.626
      SLAT( 99)=  20.748
      SLAT(100)=  21.869
      SLAT(101)=  22.991
      SLAT(102)=  24.112
      SLAT(103)=  25.234
      SLAT(104)=  26.355
      SLAT(105)=  27.476
      SLAT(106)=  28.598
      SLAT(107)=  29.719
      SLAT(108)=  30.841
      SLAT(109)=  31.962
      SLAT(110)=  33.084
      SLAT(111)=  34.205
      SLAT(112)=  35.327
      SLAT(113)=  36.448
      SLAT(114)=  37.570
      SLAT(115)=  38.691
      SLAT(116)=  39.813
      SLAT(117)=  40.934
      SLAT(118)=  42.056
      SLAT(119)=  43.177
      SLAT(120)=  44.299
      SLAT(121)=  45.420
      SLAT(122)=  46.542
      SLAT(123)=  47.663
      SLAT(124)=  48.785
      SLAT(125)=  49.906
      SLAT(126)=  51.028
      SLAT(127)=  52.149
      SLAT(128)=  53.271
      SLAT(129)=  54.392
      SLAT(130)=  55.514
      SLAT(131)=  56.635
      SLAT(132)=  57.757
      SLAT(133)=  58.878
      SLAT(134)=  60.000
      SLAT(135)=  61.121
      SLAT(136)=  62.242
      SLAT(137)=  63.364
      SLAT(138)=  64.485
      SLAT(139)=  65.607
      SLAT(140)=  66.728
      SLAT(141)=  67.850
      SLAT(142)=  68.971
      SLAT(143)=  70.093
      SLAT(144)=  71.214
      SLAT(145)=  72.336
      SLAT(146)=  73.457
      SLAT(147)=  74.578
      SLAT(148)=  75.700
      SLAT(149)=  76.821
      SLAT(150)=  77.943
      SLAT(151)=  79.064
      SLAT(152)=  80.185
      SLAT(153)=  81.307
      SLAT(154)=  82.428
      SLAT(155)=  83.549
      SLAT(156)=  84.670
      SLAT(157)=  85.791
      SLAT(158)=  86.911
      SLAT(159)=  88.029
      SLAT(160)=  89.142

      DO I=1,NLONS
         SLON(I)=FLOAT(I-1)*1.125
      ENDDO

      PPLEV(1)  =  30.
      PPLEV(2)  =  50.
      PPLEV(3)  =  70.
      PPLEV(4)  = 100.
      PPLEV(5)  = 150.
      PPLEV(6)  = 200.
      PPLEV(7)  = 250.
      PPLEV(8)  = 300.
      PPLEV(9)  = 350.
      PPLEV(10) = 420.
      PPLEV(11) = 500.
      PPLEV(12) = 600.
      PPLEV(13) = 700.
      PPLEV(14) = 780.
      PPLEV(15) = 850.
      PPLEV(16) = 920.
      PPLEV(17) = 960.
      PPLEV(18) =1000.

      DO K=1,NLEV2
         SIGMAR(K)= PPLEV(K)*0.001
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,NLEV2
         KR=NLEV2-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      AK( 1)=    0.00000000
      AK( 2)=    0.20000000
      AK( 3)=    0.38425343
      AK( 4)=    0.63647804
      AK( 5)=    0.95636963
      AK( 6)=    1.34483307
      AK( 7)=    1.80584351
      AK( 8)=    2.34779053
      AK( 9)=    2.98495789
      AK(10)=    3.73971924
      AK(11)=    4.64618134
      AK(12)=    5.75651001
      AK(13)=    7.13218079
      AK(14)=    8.83660522
      AK(15)=   10.94834717
      AK(16)=   13.56474609
      AK(17)=   16.80640259
      AK(18)=   20.82273926
      AK(19)=   25.79888672
      AK(20)=   31.96421631
      AK(21)=   39.60291504
      AK(22)=   49.06708496
      AK(23)=   60.18019531
      AK(24)=   73.06631348
      AK(25)=   87.65053711
      AK(26)=  103.76126953
      AK(27)=  120.77446289
      AK(28)=  137.75325195
      AK(29)=  153.79805664
      AK(30)=  168.19474609
      AK(31)=  180.45183594
      AK(32)=  190.27695313
      AK(33)=  197.55109375
      AK(34)=  202.22205078
      AK(35)=  204.29863281
      AK(36)=  203.84480469
      AK(37)=  200.97402344
      AK(38)=  195.84330078
      AK(39)=  188.64750000
      AK(40)=  179.61357422
      AK(41)=  168.99468750
      AK(42)=  157.06447266
      AK(43)=  144.11124023
      AK(44)=  130.43218750
      AK(45)=  116.32758789
      AK(46)=  102.09500977
      AK(47)=   88.02356445
      AK(48)=   74.38803223
      AK(49)=   61.44314941
      AK(50)=   49.41778320
      AK(51)=   38.50913330
      AK(52)=   28.87696533
      AK(53)=   20.63779785
      AK(54)=   13.85912598
      AK(55)=    8.55361755
      AK(56)=    4.67333588
      AK(57)=    2.10393890
      AK(58)=    0.65889244
      AK(59)=    0.07367743
      AK(60)=    0.00000000
      AK(61)=    0.00000000

      BK( 1)=   0.00000000
      BK( 2)=   0.00000000
      BK( 3)=   0.00000000
      BK( 4)=   0.00000000
      BK( 5)=   0.00000000
      BK( 6)=   0.00000000
      BK( 7)=   0.00000000
      BK( 8)=   0.00000000
      BK( 9)=   0.00000000
      BK(10)=   0.00000000
      BK(11)=   0.00000000
      BK(12)=   0.00000000
      BK(13)=   0.00000000
      BK(14)=   0.00000000
      BK(15)=   0.00000000
      BK(16)=   0.00000000
      BK(17)=   0.00000000
      BK(18)=   0.00000000
      BK(19)=   0.00000000
      BK(20)=   0.00000000
      BK(21)=   0.00000000
      BK(22)=   0.00000000
      BK(23)=   0.00000000
      BK(24)=   0.00000000
      BK(25)=   0.00007582
      BK(26)=   0.00046139
      BK(27)=   0.00181516
      BK(28)=   0.00508112
      BK(29)=   0.01114291
      BK(30)=   0.02067788
      BK(31)=   0.03412116
      BK(32)=   0.05169041
      BK(33)=   0.07353383
      BK(34)=   0.09967469
      BK(35)=   0.13002251
      BK(36)=   0.16438432
      BK(37)=   0.20247594
      BK(38)=   0.24393314
      BK(39)=   0.28832296
      BK(40)=   0.33515489
      BK(41)=   0.38389215
      BK(42)=   0.43396294
      BK(43)=   0.48477158
      BK(44)=   0.53570992
      BK(45)=   0.58616841
      BK(46)=   0.63554746
      BK(47)=   0.68326861
      BK(48)=   0.72878581
      BK(49)=   0.77159661
      BK(50)=   0.81125343
      BK(51)=   0.84737492
      BK(52)=   0.87965691
      BK(53)=   0.90788388
      BK(54)=   0.93194032
      BK(55)=   0.95182151
      BK(56)=   0.96764523
      BK(57)=   0.97966272
      BK(58)=   0.98827010
      BK(59)=   0.99401945
      BK(60)=   0.99763012
      BK(61)=   1.00000000

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETFVGCM(IDATE)
      implicit none
      INTEGER IDATE
      INTEGER NYEAR,MONTH,NDAY,NHOUR,NREC,MREC
      CHARACTER*5 FN_RF,FN_A2
      CHARACTER*5 PN_RF,PN_A2
      CHARACTER*4 YR_RF(30),YR_A2(30)
      CHARACTER*3 CHmon(12)
      DATA    FN_RF/'FV_RF'/  ,FN_A2/'FV_A2'/
      DATA    PN_RF/'PS_RF'/  ,PN_A2/'PS_A2'/
      DATA    YR_RF/
     & '1961','1962','1963','1964','1965','1966','1967','1968',
     & '1969','1970','1971','1972','1973','1974','1975','1976',
     & '1977','1978','1979','1980','1981','1982','1983','1984',
     & '1985','1986','1987','1988','1989','1990'/
      DATA    YR_A2/
     & '2071','2072','2073','2074','2075','2076','2077','2078',
     & '2079','2080','2081','2082','2083','2084','2085','2086',
     & '2087','2088','2089','2090','2091','2092','2093','2094',
     & '2095','2096','2097','2098','2099','2100'/
      DATA    CHmon/ 'JAN','FEB','MAR','APR','MAY','JUN',
     &               'JUL','AUG','SEP','OCT','NOV','DEC'/
      CHARACTER*20 FINM
      CHARACTER*20 FIPS
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 : PARAMETERS FOR FVGCM DATASET   X X X X
! A1
      INTEGER NLEV2,NLAT2,NLON2
      PARAMETER (NLEV2=18,NLAT2=181,NLON2=288)
!
!  NLEV1  = NUMBER OF PRESSURE LEVELS IN ECMWF DATASET.
!  NLAT1  = NUMBER OF LATITUDES ON ECMWF GRID.
!  NLON1  = NUMBER OF LONGITUDES ON ECMWF GRID.
!
!  X X X X X   SET 2 : PARAMETERS FOR RegCM DATASET   X X X X
!
!HH   change the vertical levels and the model domain size.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RegCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    VLON,VLAT,SIGMA1,SIGMAR,PPLEV,AK,BK
      COMMON /FVGRID/VLON(NLON2),VLAT(NLAT2),SIGMA1(NLEV2),SIGMAR(NLEV2)
     &              ,PPLEV(NLEV2),AK(NLEV2+1),BK(NLEV2+1)
!
!     DOMAIN VARIABLES FOR RegCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)

! DIMENSION SURFACE TEMPERATURE ON FVGCM SURFACE, GIVEN BY HADAMH DATA
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!B1
      REAL    ZS2,PS2,T2,Q2,U2,V2
      COMMON /FVVARS/ ZS2(NLON2,NLAT2),PS2(NLON2,NLAT2),
     &                T2(NLON2,NLAT2,NLEV2),Q2(NLON2,NLAT2,NLEV2),
     &                U2(NLON2,NLAT2,NLEV2),V2(NLON2,NLAT2,NLEV2)
      REAL    BB(NLON2,NLAT2,NLEV2*4+1)
      EQUIVALENCE(BB(1,1,1),PS2(1,1))
!
      REAL    Z1,HP,TP,QP,UP,VP,PP3D
      COMMON /FVVART/ Z1(NLON2,NLAT2,NLEV2),TP(NLON2,NLAT2,NLEV2),
     &                QP(NLON2,NLAT2,NLEV2),HP(NLON2,NLAT2,NLEV2),
     &                UP(NLON2,NLAT2,NLEV2),VP(NLON2,NLAT2,NLEV2),
     &              PP3D(NLON2,NLAT2,NLEV2)
      REAL    B2(NLON2,NLAT2,NLEV2*3)
      EQUIVALENCE (B2(1,1,1),TP(1,1,1))
      REAL    D2(NLON2,NLAT2,NLEV2*2)
      EQUIVALENCE (D2(1,1,1),UP(1,1,1))
!
!B3
      REAL    T3,Q3,H3,U3,V3,W3,B3PD
      COMMON /FV_SB3/ T3(JX,IY,NLEV2),Q3(JX,IY,NLEV2),
     &                H3(JX,IY,NLEV2),
     &                U3(JX,IY,NLEV2),V3(JX,IY,NLEV2),
     &                W3(JX,IY,NLEV2),B3PD(JX,IY)
      REAL    B3(JX,IY,NLEV2*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,NLEV2*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!B4
      REAL    PS4,T4,Q4,H4,TS4,U4,V4
      COMMON /FV_SB4/ PS4(JX,IY),
     &                T4(JX,IY,KZ),Q4(JX,IY,KZ),
     &                H4(JX,IY,KZ),TS4(JX,IY),
     &                U4(JX,IY,KZ),V4(JX,IY,KZ)
!
      real    lon0,lon1,lat0,lat1
      INTEGER i0,i1,j0
      COMMON /SZwindow/lon0,lon1,lat0,lat1,i0,i1,j0
      integer numx,numy,ii,i2,j2
      REAL    temp(288,181)
      real(kind=8)  scale,offset
      integer(kind=2) itmp(288,181)
!
      INTEGER I,J,K
      INTEGER NYRP,NMOP
      REAL    WT
      logical there
!
!  B1 IS FOR DATA RECORDS FROM FVGCM 1x1.25L18 GRID
!  B2 IS FOR LAT-LON GRID WITH FVGCM VERTICAL STRUCTURE
!  B3 IS FOR RegCM HORIZONTAL GRID, BUT FVGCM VERTICAL GRID
!  B4 IS FOR RegCM 3-DIMENSIONAL GRID
!
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      IF(IDATE.EQ.IDATE1) THEN
      numx=nint((lon1-lon0)/1.25)+1
      numy=nint(lat1-lat0)+1
      inquire(file='../DATA/FVGCM/HT_SRF',exist=there)
      if(.not.there) then
        write(*,*) '../DATA/FVGCM/HT_SRF is not available'
        stop
      endif
      OPEN(61,file='../DATA/FVGCM/HT_SRF',form='unformatted'
     &       ,recl=numx*numy*ibyte,access='direct')
      read(61,rec=1) ((temp(i,j),i=1,numx),j=1,numy)
      do j=nint(lat0),nint(lat1)
      do i=nint(lon0/1.25),nint(lon1/1.25)
         ii=i+1
         if(ii.le.0) ii=ii+288
         if(ii.gt.288) ii=ii-288
         i2=i-nint(lon0/1.25)+1
         j2=j-nint(lat0)+1
         ZS2(ii,j+91)=temp(i2,j2)/9.80616
      enddo
      enddo
      CLOSE(61)
      ENDIF

      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-(IDATE/10000)*100
      NHOUR=MOD(IDATE,100)

      IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
         IF(SSTTYP.EQ.'FV_RF') THEN
            FINM='RF/'//YR_RF(NYEAR-1960)//'/'//FN_RF//
     &                  YR_RF(NYEAR-1960)//CHmon(MONTH)
            FIPS='RF/'//YR_RF(NYEAR-1960)//'/'//PN_RF//
     &                  YR_RF(NYEAR-1960)//CHmon(MONTH)
         ELSE IF(SSTTYP.EQ.'FV_A2') THEN
            FINM='A2/'//YR_A2(NYEAR-2070)//'/'//FN_A2//
     &                  YR_A2(NYEAR-2070)//CHmon(MONTH)
            FIPS='A2/'//YR_A2(NYEAR-2070)//'/'//PN_A2//
     &                  YR_A2(NYEAR-2070)//CHmon(MONTH)
         ELSE
            WRITE(*,*)'ERROR IN SSTTYP'
            STOP
         ENDIF
      ELSE
         IF(.NOT.(MONTH.EQ.1)) THEN
            IF(SSTTYP.EQ.'FV_RF') THEN
               FINM='RF/'//YR_RF(NYEAR-1960)//'/'//FN_RF//
     &                     YR_RF(NYEAR-1960)//CHmon(MONTH-1)
               FIPS='RF/'//YR_RF(NYEAR-1960)//'/'//PN_RF//
     &                     YR_RF(NYEAR-1960)//CHmon(MONTH-1)
            ELSE IF(SSTTYP.EQ.'FV_A2') THEN
               FINM='A2/'//YR_A2(NYEAR-2070)//'/'//FN_A2//
     &                     YR_A2(NYEAR-2070)//CHmon(MONTH-1)
               FIPS='A2/'//YR_A2(NYEAR-2070)//'/'//PN_A2//
     &                     YR_A2(NYEAR-2070)//CHmon(MONTH-1)
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
         ELSE
            IF(SSTTYP.EQ.'FV_RF') THEN
               IF(NYEAR.EQ.1961) THEN
                  WRITE(*,*) 'Fields on 00z01jan1961 is not saved'
                  WRITE(*,*) 'Please run from 00z02jan1961'
                  STOP
               ENDIF
               FINM='RF/'//YR_RF(NYEAR-1961)//'/'//FN_RF//
     &                     YR_RF(NYEAR-1961)//CHmon(12)
               FIPS='RF/'//YR_RF(NYEAR-1961)//'/'//PN_RF//
     &                     YR_RF(NYEAR-1961)//CHmon(12)
            ELSE IF(SSTTYP.EQ.'FV_A2') THEN
               IF(NYEAR.EQ.2071) THEN
                  WRITE(*,*) 'Fields on 00z01jan2071 is not saved'
                  WRITE(*,*) 'Please run from 00z02jan2071'
                  STOP
               ENDIF
               FINM='A2/'//YR_A2(NYEAR-2071)//'/'//FN_A2//
     &                     YR_A2(NYEAR-2071)//CHmon(12)
               FIPS='A2/'//YR_A2(NYEAR-2071)//'/'//PN_A2//
     &                     YR_A2(NYEAR-2071)//CHmon(12)
            ELSE
               WRITE(*,*)'ERROR IN SSTTYP'
               STOP
            ENDIF
         ENDIF
      ENDIF
      numx=nint((lon1-lon0)/1.25)+1
      numy=nint(lat1-lat0)+1
      DO K=1,NLEV2*4+1
      DO J=1,NLAT2
      DO I=1,NLON2
         BB(I,J,K) = -9999.
      ENDDO
      ENDDO
      ENDDO
      inquire(file='../DATA/FVGCM/'//FINM,exist=there)
      if(.not.there) then
         write(*,*) '../DATA/FVGCM/'//FINM,' is not available'
         write(*,*) 'please copy FVGCM output under ../DATA/FVGCM/'
         stop
      endif
      OPEN(63,file='../DATA/FVGCM/'//FINM,form='unformatted'
     &       ,recl=(numx*numy*2+16)/4*ibyte,access='direct')
      OPEN(62,file='../DATA/FVGCM/'//FIPS,form='unformatted'
     &       ,recl=numx*numy*ibyte,access='direct')
      IF(.NOT.(NDAY.EQ.1.AND.NHOUR.EQ.0)) THEN
         nrec=((NDAY-1)*4+NHOUR/6-1)*(NLEV2*4)
         mrec= (NDAY-1)*4+NHOUR/6-1
      ELSE
         IF(MONTH.EQ.1.OR.MONTH.EQ.2.OR.MONTH.EQ.4.OR.
     &      MONTH.EQ.6.OR.MONTH.EQ.8.OR.MONTH.EQ.9.OR.
     &      MONTH.EQ.11) THEN
            nrec=(31*4-1)*(NLEV2*4)
            mrec= 31*4-1
         ELSE IF(MONTH.EQ.5.OR.MONTH.EQ.7.OR.MONTH.EQ.10
     &                     .OR.MONTH.EQ.12) THEN
            nrec=(30*4-1)*(NLEV2*4)
            mrec= 30*4-1
         ELSE
            IF(MOD(NYEAR,4).EQ.0.AND.NYEAR.NE.2100) THEN
               nrec=(29*4-1)*(NLEV2*4)
               mrec= 29*4-1
            ELSE
               nrec=(28*4-1)*(NLEV2*4)
               mrec= 28*4-1
            ENDIF
         ENDIF
      ENDIF
      mrec=mrec+1
      read(62,rec=mrec) ((temp(i,j),i=1,numx),j=1,numy)
      do j=nint(lat0),nint(lat1)
      do i=nint(lon0/1.25),nint(lon1/1.25)
         ii=i+1
         if(ii.le.0) ii=ii+288
         if(ii.gt.288) ii=ii-288
         i2=i-nint(lon0/1.25)+1
         j2=j-nint(lat0)+1
         PS2(ii,j+91)=temp(i2,j2)*0.01
      enddo
      enddo
      do k=1,NLEV2
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0),nint(lat1)
         do i=nint(lon0/1.25),nint(lon1/1.25)
            ii=i+1
            if(ii.le.0) ii=ii+288
            if(ii.gt.288) ii=ii-288
            i2=i-nint(lon0/1.25)+1
            j2=j-nint(lat0)+1
            U2(ii,j+91,k)=itmp(i2,j2)*scale+offset
         enddo
         enddo
      enddo
      do k=1,NLEV2
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0),nint(lat1)
         do i=nint(lon0/1.25),nint(lon1/1.25)
            ii=i+1
            if(ii.le.0) ii=ii+288
            if(ii.gt.288) ii=ii-288
            i2=i-nint(lon0/1.25)+1
            j2=j-nint(lat0)+1
            V2(ii,j+91,k)=itmp(i2,j2)*scale+offset
         enddo
         enddo
      enddo
      do k=1,NLEV2
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0),nint(lat1)
         do i=nint(lon0/1.25),nint(lon1/1.25)
            ii=i+1
            if(ii.le.0) ii=ii+288
            if(ii.gt.288) ii=ii-288
            i2=i-nint(lon0/1.25)+1
            j2=j-nint(lat0)+1
            T2(ii,j+91,k)=itmp(i2,j2)*scale+offset
         enddo
         enddo
      enddo
      do k=1,NLEV2
         nrec=nrec+1
         read(63,rec=nrec) offset,scale,((itmp(i,j),i=1,numx),j=1,numy)
         do j=nint(lat0),nint(lat1)
         do i=nint(lon0/1.25),nint(lon1/1.25)
            ii=i+1
            if(ii.le.0) ii=ii+288
            if(ii.gt.288) ii=ii-288
            i2=i-nint(lon0/1.25)+1
            j2=j-nint(lat0)+1
            Q2(ii,j+91,k)=itmp(i2,j2)*scale+offset
         enddo
         enddo
      enddo
      close(63)
      close(62)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE, ' from ',FINM
      DO K = 1,NLEV2
      DO J = 1,NLAT2
      DO I = 1,NLON2
         IF(PS2(I,J).GT.-9995.) THEN
            PP3d(I,J,K) = PS2(I,J)*0.5*(BK(K)+BK(K+1))
     &                  +          0.5*(AK(K)+AK(K+1))
         ELSE
            PP3d(I,J,K) = -9999.0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!
!     to calculate Heights on sigma surfaces.
      CALL HTSIG(T2,Z1,PP3D,PS2,ZS2,NLON2,NLAT2,NLEV2)
!
!     to interpolate H,U,V,T,Q and QC
!        1. For Heights
      CALL HEIGHT(HP,Z1,T2,PS2,PP3D,ZS2,NLON2,NLAT2,NLEV2,PPLEV,NLEV2)
!        2. For Zonal and Meridional Winds
      CALL INTLIN(UP,U2,PS2,PP3D,NLON2,NLAT2,NLEV2,PPLEV,NLEV2)
      CALL INTLIN(VP,V2,PS2,PP3D,NLON2,NLAT2,NLEV2,PPLEV,NLEV2)
!        3. For Temperatures
      CALL INTLOG(TP,T2,PS2,PP3D,NLON2,NLAT2,NLEV2,PPLEV,NLEV2)
!        4. For Moisture qva & qca
      CALL HUMID1FV(T2,Q2,PP3D,NLON2,NLAT2,NLEV2)
      CALL INTLIN(QP,Q2,PS2,PP3D,NLON2,NLAT2,NLEV2,PPLEV,NLEV2)
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,VLON,VLAT,NLON2,NLAT2,JX,IY,NLEV2*3)
      CALL BILINX(D3,D2,DLON,DLAT,VLON,VLAT,NLON2,NLAT2,JX,IY,NLEV2*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,NLEV2
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,NLEV2)
      CALL TOP2BTM(Q3,JX,IY,NLEV2)
      CALL TOP2BTM(H3,JX,IY,NLEV2)
      CALL TOP2BTM(U3,JX,IY,NLEV2)
      CALL TOP2BTM(V3,JX,IY,NLEV2)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,NLEV2)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
!
! F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!       INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,NLEV2)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2     DETERMINE P* AND HEIGHT.
!
! F3     INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,NLEV2)
      CALL HUMID2FV(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4     DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G      WRITE AN INITIAL FILE FOR THE RegCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADERFV
      implicit none
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 : PARAMETERS FOR FVGCM DATASET   X X X X
! A1
      INTEGER NLEV2,NLAT2,NLON2
      PARAMETER (NLEV2=18,NLAT2=181,NLON2=288)
!
!  NLEV2  = NUMBER OF PRESSURE LEVELS IN FVGCM DATASET.
!  NLAT2  = NUMBER OF LATITUDES  ON FVGCM GRID.
!  NLON2  = NUMBER OF LONGITUDES ON FVGCM GRID.
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    VLON,VLAT,SIGMA1,SIGMAR,PPLEV,AK,BK
      COMMON /FVGRID/VLON(NLON2),VLAT(NLAT2),SIGMA1(NLEV2),SIGMAR(NLEV2)
     &              ,PPLEV(NLEV2),AK(NLEV2+1),BK(NLEV2+1)
!
      INTEGER I,J,K,KR
!
      DO J=1,NLAT2
         VLAT(J)=FLOAT(J-1)-90.0
      ENDDO
!
      DO I=1,NLON2
         VLON(I)=FLOAT(I-1)*1.25
      ENDDO

      PPLEV(1)=  30.
      PPLEV(2)=  50.
      PPLEV(3)=  70.
      PPLEV(4)= 100.
      PPLEV(5)= 150.
      PPLEV(6)= 200.
      PPLEV(7)= 250.
      PPLEV(8)= 300.
      PPLEV(9)= 350.
      PPLEV(10) = 420.
      PPLEV(11) = 500.
      PPLEV(12) = 600.
      PPLEV(13) = 700.
      PPLEV(14) = 780.
      PPLEV(15) = 850.
      PPLEV(16) = 920.
      PPLEV(17) = 960.
      PPLEV(18) =1000.

      DO K=1,NLEV2
         SIGMAR(K)= PPLEV(K)*0.001
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,NLEV2
         KR=NLEV2-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE
!
      AK(1) =   2.9170
      AK(2) =   7.9292
      AK(3) =  21.5539
      AK(4) =  49.1834
      AK(5) =  83.1425
      AK(6) =  79.9308
      AK(7) =  75.7738
      AK(8) =  70.5752
      AK(9) =  64.2963
      AK(10)=  56.9838
      AK(11)=  48.7913
      AK(12)=  39.9895
      AK(13)=  30.9631
      AK(14)=  22.1902
      AK(15)=  14.2039
      AK(16)=   7.5413
      AK(17)=   2.6838
      AK(18)=   0.0
      AK(19)=   0.0
      
      BK(1) = 0.0
      BK(2) = 0.0
      BK(3) = 0.0
      BK(4) = 0.0
      BK(5) = 0.0
      BK(6) = 0.0380541
      BK(7) = 0.0873088
      BK(8) = 0.1489307
      BK(9) = 0.2232996
      BK(10)= 0.3099406
      BK(11)= 0.4070096
      BK(12)= 0.5112977
      BK(13)= 0.6182465
      BK(14)= 0.7221927
      BK(15)= 0.8168173
      BK(16)= 0.8957590
      BK(17)= 0.9533137
      BK(18)= 0.9851222
      BK(19)= 1.0
      
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETNCEP(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=13,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,PSvar
      common /CDCvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        PSvar(ilon,jlat)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALNC/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /NCPVAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,         T3(JX,IY,klev)
     &       ,         H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER NCEP/NCAR REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL CDC6HOUR(DATTYP,IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GETNCEPW(IDATE)
      IMPLICIT NONE
      INTEGER IDATE
!
! A      SET PARAMETERS
!
!  X X X X X X X                                      X X X X X X X X X
!  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=13,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
!  X X X X X   SET 2 : MODEL DOMAIN PARAMETERS FOR RCM DATASET   X X X X
!
!HH   CHANGE THE VERTICAL LEVELS AND THE MODEL DOMAIN SIZE.
!
!  JX  = NUMBER OF GRID POINTS ALONG LONGITUDES ON OUTPUT GRID.
!  IY  = NUMBER OF GRID POINTS ALONG LATITUDES ON OUTPUT GRID.
!  KZ  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN RCM DATASET.
      include 'icbc.param'
!
!  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
!  X X X X X X X                                      X X X X X X X X X
!
      REAL    Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,PSvar
      common /CDCvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        PSvar(ilon,jlat)
      REAL    B2(ilon,jlat,klev*3)
      EQUIVALENCE (B2(1,1,1),Tvar(1,1,1))
      REAL    D2(ilon,jlat,klev*2)
      EQUIVALENCE (D2(1,1,1),Uvar(1,1,1))
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALNC/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
!
! A7     DIMENSION VARIABLES FOR RCM HORIZONTAL GRID (P-LEVELS)
      REAL    U3,V3,H3,Q3,T3
      COMMON /NCPVAR3/ U3(JX,IY,klev),V3(JX,IY,klev)
     &       ,         T3(JX,IY,klev)
     &       ,         H3(JX,IY,klev),Q3(JX,IY,klev)
      REAL    B3(JX,IY,klev*3)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      REAL    D3(JX,IY,klev*2)
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
! A8     DIMENSION VARIABLES FOR RCM INPUT FILE

      REAL    U4,V4,T4,Q4,H4,PS4,TS4
      COMMON /RCMVAR4/ U4(JX,IY,KZ),V4(JX,IY,KZ),T4(JX,IY,KZ)
     &       ,         Q4(JX,IY,KZ),H4(JX,IY,KZ),PS4(JX,IY)
     &       ,         TS4(JX,IY)
      REAL    B3PD(JX,IY)
!
!----------------------------------------------------------------------
!
! DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF DATA
! READ FROM THE OUTPUT OF SST STEP
      REAL    SST1(JX,IY), SST2(JX,IY)
 
! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)
 
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
!
      INTEGER NYRP,NMOP
      REAL    WT
 
!
!  B2 IS FOR LAT-LON GRID WITH PRESSURE LEVEL STRUCTURE
!  B3 IS FOR RCM HORIZONTAL GRID, BUT WITH P-LEVEL STRUCTURE
!  B4 IS FOR RCM 3-DIMENSIONAL GRID
!                            S T A R T
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
! C  BEGIN LOOP OVER NCEP/NCAR REALALYSIS HISTORY-FILE VOLUMES
!
!HH
!HH:OVER
!
! D      BEGIN LOOP OVER NTIMES
!
      CALL CDC6HOUR2(IDATE,IDATE1)
      WRITE(*,*) 'READ IN fields at DATE:',IDATE
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL BILINX(B3,B2,XLON,XLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*3)
      CALL BILINX(D3,D2,DLON,DLAT,GLON,GLAT,ilon,jlat,JX,IY,klev*2)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,klev
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,klev)
      CALL TOP2BTM(Q3,JX,IY,klev)
      CALL TOP2BTM(H3,JX,IY,klev)
      CALL TOP2BTM(U3,JX,IY,klev)
      CALL TOP2BTM(V3,JX,IY,klev)
!HH:OVER
!
! ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR,JX,IY,klev)
 
      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
 
!
! F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR,PTOP,JX,IY,klev)
 
      IF(SSTTYP.NE.'OI_WK') THEN
! F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
      PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
      CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
      CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF
 
! F2  DETERMINE P* AND HEIGHT.
!
! F3  INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
 
      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR,PTOP,JX,IY,KZ,klev)
      CALL HUMID2(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
!
! F4  DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G   WRITE AN INITIAL FILE FOR THE RCM
      CALL WRITEF(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HEADERNC
      IMPLICIT NONE
!
!  X X X X X   SET 1 :PARAMETERS FOR NCEP/NCAR REALALYSIS DATASET X X X
! A1
      INTEGER klev,jlat,ilon
      PARAMETER(klev=13,jlat=73,ilon=144)
!
!  ilon  = NUMBER OF LONGITUDES ON NCEP GRID.
!  jlat  = NUMBER OF LATITUDES ON NCEP GRID.
!  klev  = NUMBER OF PRESSURE LEVELS IN NCEP DATASET.
!
      REAL    GLAT,GLON,SIGMA1,SIGMAR
      COMMON /GLOBALNC/ GLAT(jlat),GLON(ilon),SIGMA1(klev),SIGMAR(klev)
      INTEGER I,J,K,KR
!
      SIGMAR(1) = .07
      SIGMAR(2) = .1
      SIGMAR(3) = .15
      SIGMAR(4) = .2
      SIGMAR(5) = .25
      SIGMAR(6) = .3
      SIGMAR(7) = .4
      SIGMAR(8) = .5
      SIGMAR(9) = .6
      SIGMAR(10)= .7
      SIGMAR(11)= .85
      SIGMAR(12)= .925
      SIGMAR(13)=1.0
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      DO I = 1,ilon
         GLON(I) = FLOAT(I-1)*2.5
      ENDDO
      DO J = 1,jlat
         GLAT(J) = -90.0+FLOAT(J-1)*2.5
      ENDDO
!HH:OVER
! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      DO 116 K=1,klev
         KR=klev-K+1
         SIGMA1(K)=SIGMAR(KR)
  116 CONTINUE

      RETURN
      END
!
      SUBROUTINE GET_NEST(IDATE,NCR)
      implicit none
      INTEGER IDATE,NCR
      include 'icbc.param'
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)

! DIMENSION SURFACE TEMPERATURE ON FVGCM SURFACE, GIVEN BY HADAMH DATA
      REAL    SST1(JX,IY), SST2(JX,IY)

! ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL    PA(JX,IY), ZA(JX,IY)
      REAL    TLAYER(JX,IY)

      REAL    XLON_O,XLAT_O,HT_O,SIGF_O,SIG_O,PTOP_O
      REAL    PLEV_O,SIGMAR_O
      REAL    CLON_O,CLAT_O,GRDFAC_O,PLON_O,PLAT_O
      CHARACTER*6 LGTYPE_O
      INTEGER IOTYP
      COMMON /DOMAIN_O/ XLON_O(JX_O,IY_O),XLAT_O(JX_O,IY_O)
     &                 ,HT_O(JX_O,IY_O),SIGF_O(KZ_O+1),SIG_O(KZ_O)
     &                 ,PTOP_O,PLEV_O(NP),SIGMAR_O(NP),IOTYP
     &                 ,CLON_O,CLAT_O,GRDFAC_O,PLON_O,PLAT_O
     &                 ,LGTYPE_O
      INTEGER IDATE0
      COMMON /DATE00/ IDATE0

      INTEGER I1UR(JX,IY),I1UL(JX,IY),I1DR(JX,IY),I1DL(JX,IY)
      INTEGER J1UR(JX,IY),J1UL(JX,IY),J1DR(JX,IY),J1DL(JX,IY)
      REAL    D1XT(JX,IY)
      REAL    D1Xa(JX,IY),D1Xb(JX,IY),D1Xc(JX,IY),D1Xd(JX,IY)
      COMMON /CRESM1/I1UR,I1UL,I1DR,I1DL,J1UR,J1UL,J1DR,J1DL
     &              ,D1XT,D1Xa,D1Xb,D1Xc,D1Xd

      INTEGER I2UR(JX,IY),I2UL(JX,IY),I2DR(JX,IY),I2DL(JX,IY)
      INTEGER J2UR(JX,IY),J2UL(JX,IY),J2DR(JX,IY),J2DL(JX,IY)
      REAL    D2XT(JX,IY)
      REAL    D2Xa(JX,IY),D2Xb(JX,IY),D2Xc(JX,IY),D2Xd(JX,IY)
      COMMON /CRESM2/I2UR,I2UL,I2DR,I2DL,J2UR,J2UL,J2DR,J2DL
     &              ,D2XT,D2Xa,D2Xb,D2Xc,D2Xd
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
      REAL    U_O,V_O,T_O,Q_O,C_O,PS_O
      COMMON /FNEST1/ U_O(JX_O,IY_O,KZ_O),V_O(JX_O,IY_O,KZ_O)
     &               ,T_O(JX_O,IY_O,KZ_O),Q_O(JX_O,IY_O,KZ_O)
     &               ,C_O(JX_O,IY_O,KZ_O),PS_O(JX_O,IY_O)
!
!B1
      REAL    Z1,HP,TP,QP,CP,UP,VP,PP3D
      COMMON /RCMVAR/ Z1(JX_O,IY_O,KZ_O),TP(JX_O,IY_O,NP),
     &                QP(JX_O,IY_O,NP),CP(JX_O,IY_O,NP),
     &                HP(JX_O,IY_O,NP),
     &                UP(JX_O,IY_O,NP),VP(JX_O,IY_O,NP),
     &              PP3D(JX_O,IY_O,NP)
      REAL    B2(JX_O,IY_O,NP*4),D2(JX_O,IY_O,NP*2)
      EQUIVALENCE (B2(1,1,1),TP(1,1,1))
      EQUIVALENCE (D2(1,1,1),UP(1,1,1))
!
!B3 
      REAL    T3,Q3,C3,H3,U3,V3,W3,B3PD 
      COMMON /FNEST3/ T3(JX,IY,NP),Q3(JX,IY,NP),
     &                C3(JX,IY,NP),H3(JX,IY,NP),
     &                U3(JX,IY,NP),V3(JX,IY,NP),
     &                W3(JX,IY,NP),B3PD(JX,IY)
      REAL    B3(JX,IY,NP*4),D3(JX,IY,NP*2)
      EQUIVALENCE (B3(1,1,1),T3(1,1,1))
      EQUIVALENCE (D3(1,1,1),U3(1,1,1))
!
!B4 
      REAL    PS4,T4,Q4,C4,H4,TS4,U4,V4
      COMMON /FNEST4/ PS4(JX,IY),T4(JX,IY,KZ),
     &                Q4(JX,IY,KZ),C4(JX,IY,KZ),
     &                H4(JX,IY,KZ),TS4(JX,IY),
     &                U4(JX,IY,KZ),V4(JX,IY,KZ)
!
      INTEGER NREC
      COMMON /COUNT/ NREC
!
      INTEGER NY0,MN0,ND0,NH0,NY1,MN1,ND1,NH1
      CHARACTER*14 fillin
      logical there
      integer I,J,K,idateK
      INTEGER NYRP,NMOP
      REAL    WT
!

 100  FORMAT('ATM.',I10)
      IF(IDATE.EQ.IDATE0) THEN
        WRITE(fillin,100) IDATE
        inquire(file='../DATA/RegCM/'//fillin,exist=there)
        if(.not.there) then
          write(*,*) '../DATA/RegCM/'//fillin,' is not available'
          write(*,*) 'please copy (or link)',fillin
          stop
        endif
        if(iotyp.eq.1) then
          open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &           ,recl=IY_O*JX_O*ibyte,access='direct')
          nrec=0
        else if(iotyp.eq.2) then
          open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
          rewind(55)
        endif
      ELSE IF(IDATE.EQ.IDATE1) THEN
        NY0=IDATE0/1000000
        MN0=MOD(IDATE0/10000,100)
        ND0=MOD(IDATE0/100,100)
        NH0=MOD(IDATE0,100)

        NY1=IDATE1/1000000
        MN1=MOD(IDATE1/10000,100)
        ND1=MOD(IDATE1/100,100)
        NH1=MOD(IDATE1,100)

        IF(NY0.eq.NY1.and.MN0.eq.MN1) THEN
          WRITE(fillin,100) IDATE0
          inquire(file='../DATA/RegCM/'//fillin,exist=there)
          if(.not.there) then
            write(*,*) '../DATA/RegCM/'//fillin,' is not available'
            write(*,*) 'please copy (or link)',fillin
            stop
          endif
          if(iotyp.eq.1) then
            open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &             ,recl=IY_O*JX_O*ibyte,access='direct')
            nrec=((ND1-ND0)*4+(NH1-NH0)/6)*(KZ_O*6+5)
          else if(iotyp.eq.2) then
            open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
            rewind(55)
          endif
        ELSE IF(ND1.eq.1.and.NH1.eq.0) then
          IF((NY1-NY0)*12+(MN1-MN0).eq.1) then
            WRITE(fillin,100) IDATE0
            inquire(file='../DATA/RegCM/'//fillin,exist=there)
            if(.not.there) then
              write(*,*) '../DATA/RegCM/'//fillin,' is not available'
              write(*,*) 'please copy (or link)',fillin
              stop
            endif
            if(iotyp.eq.1) then
              open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &               ,recl=IY_O*JX_O*ibyte,access='direct')
              if(MN0.eq.1.or.MN0.eq.3.or.MN0.eq.5.or.
     &           MN0.eq.7.or.MN0.eq.8.or.MN0.eq.10.or.
     &           MN0.eq.12) THEN
                nrec=(124-(ND0-1)*4+NH0/6)*(KZ_O*6+5)
              else if(MN0.eq.4.or.MN0.eq.6.or.MN0.eq.9.or.
     &                MN0.eq.11) THEN
                nrec=(120-(ND0-1)*4+NH0/6)*(KZ_O*6+5)
              else
                nrec=112-(ND0-1)*4+NH0/6
                if(MOD(NY0,4).eq.0) nrec=nrec+4
                if(MOD(NY0,100).eq.0) nrec=nrec-4
                if(MOD(NY0,400).eq.0) nrec=nrec+4
                nrec=nrec*(KZ_O*6+5)
              endif
            else if(iotyp.eq.2) then
              open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
              rewind(55)
            endif
          ELSE
            IF(MN1.gt.1) THEN
              WRITE(fillin,100) NY1*1000000+(MN1-1)*10000+100
            ELSE
              WRITE(fillin,100) (NY1-1)*1000000+120100
            ENDIF
            inquire(file='../DATA/RegCM/'//fillin,exist=there)
            if(.not.there) then
              write(*,*) '../DATA/RegCM/'//fillin,' is not available'
              write(*,*) 'please copy (or link)',fillin
              stop
            endif
            if(iotyp.eq.1) then
              open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &               ,recl=IY_O*JX_O*ibyte,access='direct')
              if(MN0.eq.1.or.MN0.eq.3.or.MN0.eq.5.or.
     &           MN0.eq.7.or.MN0.eq.8.or.MN0.eq.10.or.
     &           MN0.eq.12) THEN
                nrec=123*(KZ_O*6+5)
              else if(MN0.eq.4.or.MN0.eq.6.or.MN0.eq.9.or.
     &                MN0.eq.11) THEN
                nrec=119*(KZ_O*6+5)
              else
                nrec=111
                if(MOD(NY0,4).eq.0) nrec=nrec+4
                if(MOD(NY0,100).eq.0) nrec=nrec-4
                if(MOD(NY0,400).eq.0) nrec=nrec+4
                nrec=nrec*(KZ_O*6+5)
              endif
            else if(iotyp.eq.2) then
              open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
              rewind(55)
            endif
          ENDIF
        ELSE
          WRITE(fillin,100) NY1*1000000+MN1*10000+100
          inquire(file='../DATA/RegCM/'//fillin,exist=there)
          if(.not.there) then
            write(*,*) '../DATA/RegCM/'//fillin,' is not available'
            write(*,*) 'please copy (or link)',fillin
            stop
          endif
          if(iotyp.eq.1) then
            open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &             ,recl=IY_O*JX_O*ibyte,access='direct')
            nrec=((ND1-1)*4+NH1/6-1)*(KZ_O*6+5)
          else if(iotyp.eq.2) then
            open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
            rewind(55)
          endif
        ENDIF
      ENDIF
!     WRITE(*,*) 'Open ATM file:', '../DATA/RegCM/'//fillin

      if(iotyp.eq.1) then
         IF(IDATE.NE.IDATE1.AND.MOD(IDATE,10000).EQ.100
     &                     .AND.NCR.EQ.1) nrec=nrec-(KZ_O*6+5)
         idateK=IDATE
         do k=KZ_O,1,-1
            nrec=nrec+1
            read(55,rec=nrec) ((U_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            nrec=nrec+1
            read(55,rec=nrec) ((V_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         nrec=nrec+KZ_O            ! skip omega
         do k=KZ_O,1,-1
            nrec=nrec+1
            read(55,rec=nrec) ((T_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            nrec=nrec+1
            read(55,rec=nrec) ((Q_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            nrec=nrec+1
            read(55,rec=nrec) ((C_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         nrec=nrec+1
         read(55,rec=nrec) ((PS_O(i,j),i=1,JX_O),j=1,IY_O)
         nrec=nrec+4
      else if(iotyp.eq.2) then
         IF(IDATE.NE.IDATE1.AND.MOD(IDATE,10000).EQ.100
     &                     .AND.NCR.EQ.1) rewind(55)
1000     read(55) idateK
         if(idateK.ne.idate) then
         do k=1,KZ_O*6+5
            read(55)
         enddo
!        WRITE(*,*) 'READ IN fields at DATE:',idateK
         goto 1000
         endif
!        idate=idateK
         
!        print*,' IDATE = ',idate
         do k=KZ_O,1,-1
            read(55) ((U_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            read(55) ((V_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            read(55) 
         enddo
         do k=KZ_O,1,-1
            read(55) ((T_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=KZ_O,1,-1
            read(55) ((Q_O(i,j,k),i=1,JX_O),j=1,IY_O)
         enddo
         do k=1,KZ_O
            read(55) ((C_O(i,j,k),i=1,JX_O),j=1,IY_O) 
         enddo
         read(55) ((PS_O(i,j),i=1,JX_O),j=1,IY_O)
         do k=1,4
            read(55)
         enddo
      endif
      WRITE(*,*) 'READ IN fields at DATE:',idateK,' from ',fillin
      IF(IDATE.NE.IDATE1.AND.MOD(IDATE,10000).EQ.100
     &                  .AND.NCR.EQ.1) THEN
        WRITE(fillin,100) IDATE
        inquire(file='../DATA/RegCM/'//fillin,exist=there)
        if(.not.there) then
          write(*,*) '../DATA/RegCM/'//fillin,' is not available'
          write(*,*) 'please copy (or link)',fillin
          stop
        endif
        if(iotyp.eq.1) then
          open(55,file='../DATA/RegCM/'//fillin,form='unformatted'
     &           ,recl=IY_O*JX_O*ibyte,access='direct')
          nrec=0
        else if(iotyp.eq.2) then
          open(55,file='../DATA/RegCM/'//fillin,form='unformatted')
          rewind(55)
        endif
!       WRITE(*,*) 'Open ATM file:', fillin
      ENDIF
!
!     to calculate Heights on sigma surfaces.
      CALL HTSIG_O(T_O,Z1,PS_O,HT_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O)
!
!     to interpolate H,U,V,T,Q and QC
!        1. For Heights
      CALL HEIGHT_O(HP,Z1,T_O,PS_O,HT_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O
     &             ,PLEV_O,NP)
!        2. For Zonal and Meridional Winds
      CALL INTLIN_O(UP,U_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O,PLEV_O,NP)
      CALL INTLIN_O(VP,V_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O,PLEV_O,NP)
!        3. For Temperatures
      CALL INTLOG_O(TP,T_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O,PLEV_O,NP)
!        4. For Moisture qva & qca
      CALL HUMID1_O(T_O,Q_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O)
      CALL INTLIN_O(QP,Q_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O,PLEV_O,NP)
      CALL INTLOG_O(CP,C_O,PS_O,SIG_O,PTOP_O,JX_O,IY_O,KZ_O,PLEV_O,NP)
      CALL UVROT4NX(UP,VP,XLON_O,XLAT_O,CLON_O,CLAT_O,GRDFAC_O
     &             ,JX_O,IY_O,NP,PLON_O,PLAT_O,LGTYPE_O)
!
! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      CALL CRESSMCR(B3,B2,XLON,XLAT,XLON_O,XLAT_O,JX,IY
     &             ,I1UR,I1UL,I1DR,I1DL,J1UR,J1UL,J1DR,J1DL
     &             ,D1XT,D1Xa,D1Xb,D1Xc,D1Xd,JX_O,IY_O,NP)
      CALL CRESSMDT(D3,D2,DLON,DLAT,XLON_O,XLAT_O,JX,IY
     &             ,I2UR,I2UL,I2DR,I2DL,J2UR,J2UL,J2DR,J2DL
     &             ,D2XT,D2Xa,D2Xb,D2Xc,D2Xd,JX_O,IY_O,NP)
!
! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      CALL UVROT4(U3,V3,DLON,DLAT,CLON,CLAT,GRDFAC,JX,IY,NP
     &           ,PLON,PLAT,LGTYPE)
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!                  V E R T I C A L   I N T E R P O L A T I O N
! 
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
!HH: CHANGE THE VERTICAL ORDER.
      CALL TOP2BTM(T3,JX,IY,NP)
      CALL TOP2BTM(Q3,JX,IY,NP)
      CALL TOP2BTM(C3,JX,IY,NP)
      CALL TOP2BTM(H3,JX,IY,NP)
      CALL TOP2BTM(U3,JX,IY,NP)
      CALL TOP2BTM(V3,JX,IY,NP)
!HH:OVER
! 
! ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      CALL INTGTB(PA,ZA,TLAYER,TOPOGM,T3,H3,SIGMAR_O,JX,IY,NP)

      CALL INTPSN(PS4,TOPOGM,PA,ZA,TLAYER,PTOP,JX,IY)
      CALL P1P2(B3PD,PS4,JX,IY)
! 
! F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!       INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      CALL INTV3(TS4,T3,PS4,SIGMAR_O,PTOP,JX,IY,NP)

      IF(SSTTYP.EQ.'EH5RF'.or.SSTTYP.EQ.'EH5A2'.or.
     &   SSTTYP.EQ.'EH5B1'.or.SSTTYP.EQ.'EHA1B') THEN
         CALL MKSST3(TS4,SST1,TOPOGM,XLANDU,JX,IY,IDATE)
      ELSE IF(SSTTYP.NE.'OI_WK') THEN
! F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
         CALL JULIAN( IDATE, NYRP, NMOP, WT )
!
         CALL MKSST(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      ELSE
         CALL MKSST2(TS4,SST1,SST2,TOPOGM,XLANDU,JX,IY,IDATE/100)
      ENDIF

! F2     DETERMINE P* AND HEIGHT.
!
! F3     INTERPOLATE U, V, T, AND Q.
      CALL INTV1(U4,U3,B3PD,SIGMA2,SIGMAR_O,PTOP,JX,IY,KZ,NP)
      CALL INTV1(V4,V3,B3PD,SIGMA2,SIGMAR_O,PTOP,JX,IY,KZ,NP)
!
      CALL INTV2(T4,T3,PS4,SIGMA2,SIGMAR_O,PTOP,JX,IY,KZ,NP)

      CALL INTV1(Q4,Q3,PS4,SIGMA2,SIGMAR_O,PTOP,JX,IY,KZ,NP)
      CALL HUMID2FV(T4,Q4,PS4,PTOP,SIGMA2,JX,IY,KZ)
      CALL INTV1(C4,C3,PS4,SIGMA2,SIGMAR_O,PTOP,JX,IY,KZ,NP)
!
! F4     DETERMINE H
      CALL HYDROST(H4,T4,TOPOGM,PS4,PTOP,SIGMAF,SIGMA2,DSIGMA,JX,IY,KZ)
!
! G      WRITE AN INITIAL FILE FOR THE RegCM
      CALL WRITEF2(U4,V4,T4,Q4,PS4,TS4,PTOP,JX,IY,KZ,IDATE)
!
      RETURN
      END
!
      SUBROUTINE HEADNEST
      implicit none
      include 'icbc.param'
      REAL    XLON_O,XLAT_O,HT_O,SIGF_O,SIG_O,PTOP_O
      REAL    PLEV_O,SIGMAR_O
      REAL    CLON_O,CLAT_O,GRDFAC_O,PLON_O,PLAT_O
      CHARACTER*6 LGTYPE_O
      INTEGER IOTYP
      COMMON /DOMAIN_O/ XLON_O(JX_O,IY_O),XLAT_O(JX_O,IY_O)
     &                 ,HT_O(JX_O,IY_O),SIGF_O(KZ_O+1),SIG_O(KZ_O)
     &                 ,PTOP_O,PLEV_O(NP),SIGMAR_O(NP),IOTYP
     &                 ,CLON_O,CLAT_O,GRDFAC_O,PLON_O,PLAT_O
     &                 ,LGTYPE_O
      INTEGER IDATE0
      COMMON /DATE00/ IDATE0
!
      REAL    GLONMX,GLONMN,GLATMX,GLATMN
      REAL    ALONMX,ALONMN,ALATMX,ALATMN,PI
      INTEGER IMXMN,LCROSS,LDOT
      COMMON /MXNCOM/ GLONMX,GLONMN,GLATMX,GLATMN
     &               ,ALONMX,ALONMN,ALATMX,ALATMN
     &               ,PI,IMXMN,LCROSS,LDOT
!
      integer mdate0,ibltyp,icup,ipptls,iboudy,il,jl,kl
      real    dxsp,ptsp,clat,clon,plat,plon
      character*6 proj
      real    dto,dtb,dtr,dtc
      real    TRUELAT1,TRUELAT2
      real    d2r,SIGN
      integer k
      logical there
!
      PLEV_O(1) =  50.
      PLEV_O(2) =  70.
      PLEV_O(3) = 100.
      PLEV_O(4) = 150.
      PLEV_O(5) = 200.
      PLEV_O(6) = 250.
      PLEV_O(7) = 300.
      PLEV_O(8) = 400.
      PLEV_O(9) = 500.
      PLEV_O(10)= 600.
      PLEV_O(11)= 700.
      PLEV_O(12)= 775.
      PLEV_O(13)= 850.
      PLEV_O(14)= 925.
      PLEV_O(15)=1000.

      DO K=1,NP
         SIGMAR_O(K) = PLEV_O(K)*0.001
      ENDDO

      inquire(file='../DATA/RegCM/OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) '../DATA/RegCM/OUT_HEAD is not available'
         write(*,*) 'please copy (or link) the previous output OUT_HEAD'
         stop
      endif
      open(49,file='../DATA/RegCM/OUT_HEAD',form='unformatted'
     &         ,access='direct',recl=IY_O*JX_O*ibyte)
      read(49,rec=1) mdate0,ibltyp,icup,ipptls,iboudy
     &    , il,jl,kl,(SIGF_O(k),k=KZ_O+1,1,-1),dxsp,ptsp,clat,clon
     &    , plat,plon,proj,dto,dtb,dtr,dtc,iotyp,TRUELAT1,TRUELAT2
      if(il-2.ne.IY_O.or.jl-2.ne.JX_O.or.kl.ne.KZ_O) then
         write(*,*) 'domain parameters are not consistent'
         write(*,*) 'il,jl,kl,iy,jx,kx',il,jl,kl,IY_O,JX_O,KZ_O
         stop
      endif
      IDATE0=mdate0
      PTOP_O=ptsp*10.
      CLAT_O=clat
      CLON_O=clon
      PLAT_O=plat
      PLON_O=plon
      LGTYPE_O=proj
      IF(LGTYPE_O.EQ.'LAMCON') then
        d2r=ATAN(1.)*4./180.
        IF(clat.LT.0.) THEN
          SIGN=-1.         ! SOUTH HEMESPHERE
        ELSE
          SIGN= 1.         ! NORTH HEMESPHERE
        ENDIF
        IF (ABS(TRUELAT1-TRUELAT2) .GT. 1.E-1) THEN
         GRDFAC_O=(ALOG10(COS(TRUELAT1*d2r))-ALOG10(COS(TRUELAT2*d2r)))
     &           /(ALOG10(TAN((45.0-SIGN*TRUELAT1/2.0)*d2r))-
     &             ALOG10(TAN((45.0-SIGN*TRUELAT2/2.0)*d2r)))
        ELSE
         GRDFAC_O=SIGN*SIN(TRUELAT1*d2r)
        ENDIF
      ELSE IF(LGTYPE_O.EQ.'POLSTR') then
         GRDFAC_O=1.0
      ELSE IF(LGTYPE_O.EQ.'NORMER') then
         GRDFAC_O=0.0
      ELSE
         GRDFAC_O=0.0
      ENDIF
      read(49,rec=2) HT_O
      read(49,rec=6) XLAT_O
      read(49,rec=7) XLON_O
      close(49)

      do k=1,KZ_O
         SIG_O(k) = 0.5*(SIGF_O(k)+SIGF_O(k+1))
      enddo

      IMXMN=0
      LCROSS=0
      LDOT=0

      return
      end
!
      SUBROUTINE HEIGHT_O(HP,H,T,PSTAR,HT,SIG,PTOP,IM,JM,KM,P,KP)

!  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
!     ON INPUT:
!        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
!        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
!        SIG = SIGMA LEVELS.
!        P = PRESSURE LEVELS DESIRED.
!     ON OUTPUT:
!        ALL FIELDS EXCEPT H ARE UNCHANGED.
!        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
!
!  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
!  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
!     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
!     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
!     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
!     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
!      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
!
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    T(IM,JM,KM),H(IM,JM,KM),HP(IM,JM,KP)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      REAL    PSIG(100)
      INTEGER I,J,K,KBC,N,KT,KB
      REAL    PSFC,TEMP,WT,WB
!
      RGAS   = 287.04
      GRAV   =   9.80616
      BLTOP  =    .96
      TLAPSE =  -6.5E-3
!
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
!     PRINT *,'FIRST SIGMA LEVEL ABOVE BNDY LAYER:', SIG(KBC)
!
      DO J=1,JM
        DO I=1,IM
          DO K=1,KM
            PSIG(K) = SIG(K) * (PSTAR(I,J)-PTOP) + PTOP
          ENDDO
          PSFC = PSTAR(I,J)
          DO N = 1,KP
            KT = 1
            DO K=1,KM
              IF (PSIG(K).LT.P(N)) KT=K
            ENDDO
            KB = KT + 1
            IF(P(N).LE.PSIG(1)) THEN
              TEMP = T(I,J,1)
              HP(I,J,N) =H(I,J,1)+RGAS*TEMP*LOG(PSIG(1)/P(N))/GRAV
            ELSE IF((P(N).GT.PSIG(1)).AND.(P(N).LT.PSIG(KM))) THEN
              WT = LOG(PSIG(KB)/P(N)) / LOG(PSIG(KB)/PSIG(KT))
              WB = LOG(P(N)/PSIG(KT)) / LOG(PSIG(KB)/PSIG(KT))
              TEMP = WT * T(I,J,KT) + WB * T(I,J,KB)
              TEMP = ( TEMP + T(I,J,KB) ) / 2.
              HP(I,J,N) =H(I,J,KB)+RGAS*TEMP*LOG(PSIG(KB)/P(N))/GRAV
            ELSE IF((P(N).GE.PSIG(KM)).AND.(P(N).LE.PSFC)) THEN
              TEMP = T(I,J,KM)
              HP(I,J,N) =HT(I,J)+RGAS*TEMP*LOG(PSFC/P(N))/GRAV
            ELSE IF(P(N).GT.PSFC) THEN
              TEMP = T(I,J,KBC) - TLAPSE * (H(I,J,KBC)-HT(I,J))
              HP(I,J,N) =HT(I,J)-(TEMP/TLAPSE)
     &              * ( 1.-EXP(-RGAS*TLAPSE*LOG(P(N)/PSFC)/GRAV))
!
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE HTSIG_O(T,H,PSTAR,HT,SIG,PTOP,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM)
      REAL    PSTAR(IM,JM),HT(IM,JM)
      REAL    SIG(KM)
      REAL    PTOP,RGAS,GRAV
      INTEGER I,J,K
      REAL    TBAR
!
      RGAS   = 287.04
      GRAV   =   9.80616
!
      DO J=1,JM
      DO I=1,IM
         H(I,J,KM) = HT(I,J) + RGAS/GRAV*T(I,J,KM)
     &             * LOG(PSTAR(I,J)/((PSTAR(I,J)-PTOP)*SIG(KM)+PTOP))
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM
      DO I=1,IM
         TBAR = 0.5*( T(I,J,K)+T(I,J,K+1) )
         H(I,J,K) = H(I,J,K+1) +RGAS/GRAV*TBAR
     &            * LOG(((PSTAR(I,J)-PTOP)*SIG(K+1)+PTOP)
     &                 /((PSTAR(I,J)-PTOP)*SIG(K)+PTOP))
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE HUMID1_O(T,Q,PS,SIGMA,PTOP,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    TR,QMIN
      PARAMETER (TR=1./273.16)
      PARAMETER (QMIN=0.0)   ! MINIMUM VALUE OF SPECIFIC HUMIDITY
      REAL    T(IM,JM,KM),Q(IM,JM,KM)
      REAL    PS(IM,JM)
      REAL    SIGMA(KM)
      REAL    PTOP
      INTEGER I,J,K
      REAL    HL,SATVP,QS,P
!
!  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!  DATA ON SIGMA LEVELS
!
      DO K=1,KM
        DO J=1,JM
          DO I=1,IM
            P=SIGMA(K)*(PS(I,J)-PTOP)+PTOP
            HL=597.3-.566*(T(I,J,K)-273.16)           ! LATENT HEAT OF EVAP.
            SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K))) ! SATURATION VAP PRESS.
            QS=.622*SATVP/(P-SATVP)                   ! SAT. MIXING RATIO
            Q(I,J,K)=amax1(Q(I,J,K)/QS,0.0)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE INTLIN_O(FP,F,PSTAR,SIG,PTOP,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP
      INTEGER I,J,K,N
      INTEGER K1,K1P
      REAL    SIGP,WP,W1
!
!  INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
!        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
!        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.

      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = (SIGP-SIG(K1))/(SIG(K1P)-SIG(K1))
              W1  = 1.-WP
              FP(I,J,N)  = W1*F(I,J,K1)+WP*F(I,J,K1P)
            ELSE IF(SIGP.GE.SIG(KM)) THEN
              FP(I,J,N)  = F(I,J,KM)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE INTLOG_O(FP,F,PSTAR,SIG,PTOP,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PSTAR(IM,JM)
      REAL    SIG(KM),P(KP)
      REAL    PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P,KBC
      REAL    SIGP,WP,W1
!
      RGAS   = 287.04
      GRAV   =   9.80616
      BLTOP  =    .96
      TLAPSE =  -6.5E-3
!
!  INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATURES IN THE LAYER.

!
!** FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      DO K=1,KM
        IF(SIG(K).LT.BLTOP) KBC = K
      ENDDO
      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / (PSTAR(I,J)-PTOP)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = LOG(SIGP/SIG(K1)) / LOG(SIG(K1P)/SIG(K1))
              W1  = 1. - WP
              FP(I,J,N)= W1*F(I,J,K1) + WP*F(I,J,K1P)
            ELSE IF((SIGP.GE.SIG(KM)).AND.(SIGP.LE.1.))THEN
              FP(I,J,N)= F(I,J,KM)
            ELSE IF(SIGP.GT.1.) THEN
              FP(I,J,N) = F(I,J,KBC) 
     &        * EXP(-RGAS*TLAPSE*LOG(SIGP/SIG(KBC))/GRAV)
!                  ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GRADSCTL(FINAME,IDATE,NUMBER)
      IMPLICIT NONE
      CHARACTER FINAME*(*)
      INTEGER IDATE,NUMBER
      include 'icbc.param'
      LOGICAL there
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer nyear,month,nday,nhour
      integer i,j
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer ny,nx
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      INTEGER K
      integer isystm,system
!      external system
!
      inquire(file=FINAME//'.CTL',exist=there)
      if(there) isystm=system('/bin/rm '//FINAME//'.CTL')
      OPEN(71,file=FINAME//'.CTL',status='new')
      write(71,'(a)') 'dset ^'//FINAME(13:26)
      write(71,'(a)') 'title ICBC fields for RegCM domain'
      if(ibigend.eq.1) then
         write(71,'(a)') 'options big_endian'
      else
         write(71,'(a)') 'options little_endian'
      endif
      write(71,'(a)') 'undef -9999.'
      if(LGTYPE.eq.'LAMCON'.or.LGTYPE.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx
            if(xlat(j,1 ).lt.alatmin)alatmin=xlat(j,1 )
            if(xlat(j,iy).gt.alatmax)alatmax=xlat(j,iy)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy
         do j=1,jx
            if(clon.ge.0.0) then
               if(xlon(j,i).ge.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt.
     &               abs(clon-(xlon(j,i)+360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)+360.)
                     alonmax = amax1(alonmax,xlon(j,i)+360.)
                  endif
               endif
            else
               if(xlon(j,i).lt.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt.
     &               abs(clon-(xlon(j,i)-360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)-360.)
                     alonmax = amax1(alonmax,xlon(j,i)-360.)
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=DELX*0.001/111./2.
         rloninc=DELX*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

         centerj=jx/2.
         centeri=iy/2.
      endif
      if(LGTYPE.eq.'LAMCON') then        ! Lambert projection
         write(71,100) jx,iy,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,DELX,DELX
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(71,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(71,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(LGTYPE.eq.'POLSTR') then    !
      elseif(LGTYPE.eq.'NORMER') then
         write(71,200)  jx,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(71,210) iy
 210  format('ydef ',I3,' levels')
         write(71,220) (xlat(1,i),i=1,iy)
 220  format(10f7.2)
      elseif(LGTYPE.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(71,230) jx,iy,plon,plat,DELX/111000.
     &                                ,DELX/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(71,110) nx+2,alonmin-rloninc,rloninc
         write(71,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is correct ?'
         stop
      endif
      write(71,300)KZ,((1013.25-PTOP*10.)*SIGMA2(K)+PTOP*10.,K=KZ,1,-1)
 300  format('zdef ',I2,' levels ',30f7.2)
      nyear=IDATE/1000000
      month=(IDATE-nyear*1000000)/10000
      nday =(IDATE-nyear*1000000-month*10000)/100
      nhour=MOD(IDATE,100)
      write(71,400) NUMBER,nhour,cday(nday),cmonth(month),nyear
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' 6hr')
      if(DATTYP.eq.'EH5OM') THEN
         if(EHSO4) THEN
            if(LSMTYP.eq.'USGS') THEN
               write(71,501) 21
            else
               write(71,500) 8
            endif
         else
            if(LSMTYP.eq.'USGS') THEN
               write(71,501) 20
            else
               write(71,500) 7
            endif
         endif
      else
         if(LSMTYP.eq.'USGS') THEN
            write(71,501) 20
         else
            write(71,500) 7
         endif
      endif
 500  format('vars ',I1)
 501  format('vars ',I2)
      write(71,'(a)') 'date 0 99 header information'
 600  format(A4,I2,' 0 ',A17)
 611  format(A4,I2,' 33,100 ',A17)
 612  format(A4,I2,' 34,100 ',A17)
      if(LGTYPE.eq.'LAMCON') then        ! Lambert projection
         write(71,611) 'u   ',KZ,'westerly wind    '
         write(71,612) 'v   ',KZ,'southerly wind   '
      else
         write(71,600) 'u   ',KZ,'westerly wind    '
         write(71,600) 'v   ',KZ,'southerly wind   '
      endif
      write(71,600) 't   ',KZ,'air temperature  '
      write(71,600) 'q   ',KZ,'specific moisture'
 700  format(A4,'0 99 ',A26)
      write(71,700) 'px  ','surface pressure           '
      write(71,700) 'ts  ','surface air temperature    '
      if(DATTYP.eq.'EH5OM'.and.EHSO4) THEN
         write(71,600) 'so4 ',KZ,'sulfate amount   '
      endif
      if(LSMTYP.eq.'USGS') THEN
         write(71,700) 'qs1 ','soil moisture level 1      '
         write(71,700) 'qs2 ','soil moisture level 2      '
         write(71,700) 'qs3 ','soil moisture level 3      '
         write(71,700) 'qs4 ','soil moisture level 4      '
         write(71,700) 'ti1 ','ice  temperature level 1   '
         write(71,700) 'ti2 ','ice  temperature level 2   '
         write(71,700) 'ti3 ','ice  temperature level 3   '
         write(71,700) 'ti4 ','ice  temperature level 4   '
         write(71,700) 'ts1 ','soil temperature level 1   '
         write(71,700) 'ts2 ','soil temperature level 2   '
         write(71,700) 'ts3 ','soil temperature level 3   '
         write(71,700) 'ts4 ','soil temperature level 4   '
         write(71,700) 'snd ','snow depth (in metre)      '
      endif
      write(71,'(a)') 'endvars'
      close(71)
!
      return
      end
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GRADSCTLb(FINAME,IDATE,NUMBER)
      IMPLICIT NONE
      CHARACTER FINAME*(*)
      INTEGER IDATE,NUMBER
      include 'icbc.param'
      LOGICAL there
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer nyear,month,nday,nhour
      integer i,j
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer ny,nx
!
      CHARACTER*6 LGTYPE
      REAL    PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
      INTEGER IGRADS,IBIGEND
      COMMON /LGRID2/ PTOP,CLAT,CLON,PLAT,PLON,DELX,GRDFAC
     &              , IGRADS,IBIGEND,LGTYPE
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
      REAL    XLON,XLAT,DLON,DLAT,CORIOL,XLANDU,SNOWCV,TOPOGM,TOPOSDGM
      REAL    MSFX,SIGMA2,SIGMAF,DSIGMA
      COMMON /DOMAIN/ XLON(JX,IY),XLAT(JX,IY),DLON(JX,IY),DLAT(JX,IY)
     &       ,CORIOL(JX,IY),XLANDU(JX,IY),SNOWCV(JX,IY),TOPOGM(JX,IY)
     &       ,TOPOSDGM(JX,IY)
     &       ,MSFX(JX,IY),SIGMA2(KZ),SIGMAF(KZ+1),DSIGMA(KZ)
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      INTEGER K
      integer isystm,system
!      external system
!
      inquire(file=FINAME//'.CTL',exist=there)
      if(there) isystm=system('/bin/rm '//FINAME//'.CTL')
      OPEN(71,file=FINAME//'.CTL',status='new')
      write(71,'(a)') 'dset ^'//FINAME(13:26)
      write(71,'(a)') 'title ICBC fields for RegCM domain'
      if(ibigend.eq.1) then
         write(71,'(a)') 'options big_endian'
      else
         write(71,'(a)') 'options little_endian'
      endif
      write(71,'(a)') 'undef -9999.'
      if(LGTYPE.eq.'LAMCON'.or.LGTYPE.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx
            if(xlat(j,1 ).lt.alatmin)alatmin=xlat(j,1 )
            if(xlat(j,iy).gt.alatmax)alatmax=xlat(j,iy)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy
         do j=1,jx
            if(clon.ge.0.0) then
               if(xlon(j,i).ge.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt.
     &               abs(clon-(xlon(j,i)+360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)+360.)
                     alonmax = amax1(alonmax,xlon(j,i)+360.)
                  endif
               endif
            else
               if(xlon(j,i).lt.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt.
     &               abs(clon-(xlon(j,i)-360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)-360.)
                     alonmax = amax1(alonmax,xlon(j,i)-360.)
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=DELX*0.001/111./2.
         rloninc=DELX*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

         centerj=jx/2.
         centeri=iy/2.
      endif
      if(LGTYPE.eq.'LAMCON') then        ! Lambert projection
         write(71,100) jx,iy,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,DELX,DELX
 100  format('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
         write(71,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(71,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(LGTYPE.eq.'POLSTR') then    !
      elseif(LGTYPE.eq.'NORMER') then
         write(71,200)  jx,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(71,210) iy
 210  format('ydef ',I3,' levels')
         write(71,220) (xlat(1,i),i=1,iy)
 220  format(10f7.2)
      elseif(LGTYPE.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(71,230) jx,iy,plon,plat,DELX/111000.
     &                                ,DELX/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(71,110) nx+2,alonmin-rloninc,rloninc
         write(71,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is correct ?'
         stop
      endif
      write(71,300)KZ,((1013.25-PTOP*10.)*SIGMA2(K)+PTOP*10.,K=KZ,1,-1)
 300  format('zdef ',I2,' levels ',30f7.2)
      nyear=IDATE/1000000
      month=(IDATE-nyear*1000000)/10000
      nday =(IDATE-nyear*1000000-month*10000)/100
      nhour=MOD(IDATE,100)
      write(71,400) NUMBER,nhour,cday(nday),cmonth(month),nyear
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' 6hr')
      if(LSMTYP.eq.'USGS') THEN
         write(71,501) 20
      else
         write(71,500) 7
      endif
 500  format('vars ',I1)
 501  format('vars ',I2)
      write(71,'(a)') 'date 0 99 header information'
 600  format(A4,I2,' 0 ',A17)
      write(71,600) 'u   ',KZ,'westerly wind    '
      write(71,600) 'v   ',KZ,'southerly wind   '
      write(71,600) 't   ',KZ,'air temperature  '
      write(71,600) 'q   ',KZ,'specific moisture'
 700  format(A4,'0 99 ',A26)
      write(71,700) 'px  ','surface pressure           '
      write(71,700) 'ts  ','surface air temperature    '
      if(LSMTYP.eq.'USGS') THEN
         write(71,700) 'qs1 ','soil moisture level 1      '
         write(71,700) 'qs2 ','soil moisture level 2      '
         write(71,700) 'qs3 ','soil moisture level 3      '
         write(71,700) 'qs4 ','soil moisture level 4      '
         write(71,700) 'ti1 ','ice  temperature level 1   '
         write(71,700) 'ti2 ','ice  temperature level 2   '
         write(71,700) 'ti3 ','ice  temperature level 3   '
         write(71,700) 'ti4 ','ice  temperature level 4   '
         write(71,700) 'ts1 ','soil temperature level 1   '
         write(71,700) 'ts2 ','soil temperature level 2   '
         write(71,700) 'ts3 ','soil temperature level 3   '
         write(71,700) 'ts4 ','soil temperature level 4   '
         write(71,700) 'snd ','snow depth (in metre)      '
      endif
      write(71,'(a)') 'endvars'
      close(71)
!
      return
      end
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE GRIDML(XLON,XLAT,DLON,DLAT,ZS,ZSSD,XLANDU,MSFX,
     &                  PTOP,SIGMAF,CLON,CLAT,PLON,PLAT,DELX,GRDFAC,
     &                  JX,IY,KZ,DATTYP,CGTYPE,igrads,ibigend,ibyte)
      IMPLICIT NONE
      CHARACTER*6 CGTYPE
      CHARACTER*5 DATTYP
      INTEGER JX,IY,KZ
      REAL    XLON(JX,IY), XLAT(JX,IY),DLON(JX,IY), DLAT(JX,IY)
     B,       ZS(JX,IY), ZSSD(JX,IY), XLANDU(JX,IY),
     &        MSFX(JX,IY),SIGMAF(KZ+1)
      INTEGER igrads,ibigend,ibyte,ierr
      REAL    PTOP,CLON,CLAT,PLON,PLAT,DELX,GRDFAC
      INTEGER IYY,JXX,KZZ
      INTEGER K
      real    lon0,lon1,lat0,lat1
      INTEGER i0,i1,j0
      COMMON /SZwindow/lon0,lon1,lat0,lat1,i0,i1,j0
      REAL    TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
!
!  THIS SUBROUTINE CALLS ROUTINES TO PRODUCE THE MAP FACTORS
!  IT ALSO READS A FILE OF TOPOGRAPHY AND LANDUSE APPROPRIATE FOR GRID
!  FOR EXPLANATION OF VARIABLES SEE SUBROUTINE MAPRON.
!
!  READ APPROPRIATE FILE OF TERRAIN AND LANDUSE FOR THIS GRID
!
      open(10,file='../../Input/DOMAIN.INFO',form='unformatted'
     &       ,recl=JX*IY*ibyte,access='direct')
      IF(DATTYP.EQ.'FVGCM'.OR.DATTYP.EQ.'NRP2W'.OR.DATTYP.EQ.'GFS11'
     &                    .OR.DATTYP.EQ.'EH5OM') THEN
      READ(10,rec=1,iostat=ierr) IYY,JXX,KZZ,DELX,CLAT,CLON,PLAT,PLON
     &                          ,GRDFAC,CGTYPE,(SIGMAF(K),K=1,KZ+1),PTOP
     &                          ,igrads,ibigend
     &                          ,TRUELATL,TRUELATH,lon0,lon1,lat0,lat1
      ELSE
      READ(10,rec=1,iostat=ierr) IYY,JXX,KZZ,DELX,CLAT,CLON,PLAT,PLON
     &                          ,GRDFAC,CGTYPE,(SIGMAF(K),K=1,KZ+1),PTOP
     &                          ,igrads,ibigend,TRUELATL,TRUELATH
      ENDIF
      IF(IYY.NE.IY.OR.JXX.NE.JX.OR.KZZ.NE.KZ) THEN
         print*,'IMPROPER DIMENSION SPECIFICATION (ICBC.f)'
         print*,'  icbc.param: ',IY,JX,KZ
         print*,'  DOMAIN.INFO: ',IYY,JXX,KZZ
         print*,'  Also check ibyte in icbc.param: ibyte= ',ibyte
         STOP 'Dimensions (subroutine gridml)'
      ENDIF
      READ(10,rec=2,iostat=ierr) ZS
      READ(10,rec=3,iostat=ierr) ZSSD
      READ(10,rec=4,iostat=ierr) XLANDU
      READ(10,rec=5,iostat=ierr) XLAT
      READ(10,rec=6,iostat=ierr) XLON
      READ(10,rec=7,iostat=ierr) DLAT
      READ(10,rec=8,iostat=ierr) DLON
      READ(10,rec=9,iostat=ierr) MSFX
      CLOSE(10)
      if (ierr.ne.0) then
        print*,'END OF FILE REACHED (ICBC.f)'
        print*,'  Check ibyte in icbc.param: ibyte= ',ibyte
        stop 'EOF (subroutine gridml)'
      endif
!
      RETURN
      END
      SUBROUTINE HEADWK
      implicit none
      integer WKDAY
      COMMON /DATEWK/ WKDAY(427+1045)
      integer I,MYEAR,MONTH,MDAY
!
      WKDAY(1)=19811029
      DO I=2,427
         WKDAY(I)=WKDAY(I-1)+7
         MYEAR=WKDAY(I)/10000
         MONTH=WKDAY(I)/100-MYEAR*100
         MDAY =MOD(WKDAY(I),10000)-MONTH*100
         IF(MONTH.EQ.1.OR.MONTH.EQ.3.OR.MONTH.EQ.5.OR.
     &      MONTH.EQ.7.OR.MONTH.EQ.8.OR.MONTH.EQ.10) THEN
            IF(MDAY.GT.31) THEN
               MDAY =MDAY-31
               MONTH=MONTH+1
            ENDIF
         ELSE IF(MONTH.EQ.12) THEN
            IF(MDAY.GT.31) THEN
               MDAY =MDAY-31
               MONTH=1
               MYEAR=MYEAR+1
            ENDIF
         ELSE IF(MONTH.EQ.4.OR.MONTH.EQ.6.OR.MONTH.EQ.9.OR.
     &           MONTH.EQ.11)THEN
            IF(MDAY.GT.30) THEN
               MDAY =MDAY-30
               MONTH=MONTH+1
            ENDIF
         ELSE
            IF(MOD(MYEAR,4).NE.0) THEN
               IF(MDAY.GT.28) THEN
                  MDAY =MDAY-28
                  MONTH=MONTH+1
               ENDIF
            ELSE IF(MOD(MYEAR,400).EQ.0) THEN
               IF(MDAY.GT.29) THEN
                  MDAY =MDAY-29
                  MONTH=MONTH+1
               ENDIF
            ELSE IF(MOD(MYEAR,100).EQ.0) THEN
               IF(MDAY.GT.28) THEN
                  MDAY =MDAY-28
                  MONTH=MONTH+1
               ENDIF
            ELSE
               IF(MDAY.GT.29) THEN
                  MDAY =MDAY-29
                  MONTH=MONTH+1
               ENDIF
            ENDIF
         ENDIF
         WKDAY(I)=MYEAR*10000+MONTH*100+MDAY
      ENDDO
!
      WKDAY(428)=19891231
      DO I=429,427+1045
         WKDAY(I)=WKDAY(I-1)+7
         MYEAR=WKDAY(I)/10000
         MONTH=WKDAY(I)/100-MYEAR*100
         MDAY =MOD(WKDAY(I),10000)-MONTH*100
         IF(MONTH.EQ.1.OR.MONTH.EQ.3.OR.MONTH.EQ.5.OR.
     &      MONTH.EQ.7.OR.MONTH.EQ.8.OR.MONTH.EQ.10) THEN
            IF(MDAY.GT.31) THEN
               MDAY =MDAY-31
               MONTH=MONTH+1
            ENDIF
         ELSE IF(MONTH.EQ.12) THEN
            IF(MDAY.GT.31) THEN
               MDAY =MDAY-31
               MONTH=1
               MYEAR=MYEAR+1
            ENDIF
         ELSE IF(MONTH.EQ.4.OR.MONTH.EQ.6.OR.MONTH.EQ.9.OR.
     &           MONTH.EQ.11)THEN
            IF(MDAY.GT.30) THEN
               MDAY =MDAY-30
               MONTH=MONTH+1
            ENDIF
         ELSE
            IF(MOD(MYEAR,4).NE.0) THEN
               IF(MDAY.GT.28) THEN
                  MDAY =MDAY-28
                  MONTH=MONTH+1
               ENDIF
            ELSE IF(MOD(MYEAR,400).EQ.0) THEN
               IF(MDAY.GT.29) THEN
                  MDAY =MDAY-29
                  MONTH=MONTH+1
               ENDIF
            ELSE IF(MOD(MYEAR,100).EQ.0) THEN
               IF(MDAY.GT.28) THEN
                  MDAY =MDAY-28
                  MONTH=MONTH+1
               ENDIF
            ELSE
               IF(MDAY.GT.29) THEN
                  MDAY =MDAY-29
                  MONTH=MONTH+1
               ENDIF
            ENDIF
         ENDIF
         WKDAY(I)=MYEAR*10000+MONTH*100+MDAY
      ENDDO
!
      RETURN
      END
!
      SUBROUTINE HEIGHT(HP,H,T,PS,P3D,HT,IM,JM,KM,P,KP)

!  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
!     ON INPUT:
!        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
!        PS = SURFACE PRESSURE
!        PTOP = MODEL TOP PRESSURE.
!        SIG = SIGMA LEVELS.
!        P = PRESSURE LEVELS DESIRED.
!     ON OUTPUT:
!        ALL FIELDS EXCEPT H ARE UNCHANGED.
!        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
!
!  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
!  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
!     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
!     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
!     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
!     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
!      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
!
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    T(IM,JM,KM),H(IM,JM,KM),HP(IM,JM,KP)
      REAL    PS(IM,JM),P3D(IM,JM,KM),HT(IM,JM)
      REAL    SIG(60),P(KP)
      REAL    RGAS,GRAV,BLTOP,TLAPSE
      REAL    PSIG(61)
      INTEGER I,J,K,KBC,N,KT,KB
      REAL    PSFC,TEMP,WT,WB
!
      RGAS   =287.04
      GRAV   =  9.80616
      BLTOP  =  0.96
      TLAPSE = -6.5E-3
!
      DO J=1,JM
      DO I=1,IM
         PSFC = PS(I,J)
         IF(PSFC.GT.-9995.0) THEN
         DO K=1,KM
            SIG(K) = P3D(I,J,K) / PS(I,J)
            IF(SIG(K).LT.BLTOP) KBC = K
            PSIG(K) = P3D(I,J,K)
         ENDDO
         DO N = 1,KP
            KT = 1
            DO K=1,KM
              IF (PSIG(K).LT.P(N)) KT=K
            ENDDO
            KB = KT + 1
            IF(P(N).LE.PSIG(1)) THEN
              TEMP = T(I,J,1)
              HP(I,J,N) =H(I,J,1)+RGAS*TEMP*LOG(PSIG(1)/P(N))/GRAV
            ELSE IF((P(N).GT.PSIG(1)).AND.(P(N).LT.PSIG(KM))) THEN
              WT = LOG(PSIG(KB)/P(N)) / LOG(PSIG(KB)/PSIG(KT))
              WB = LOG(P(N)/PSIG(KT)) / LOG(PSIG(KB)/PSIG(KT))
              TEMP = WT * T(I,J,KT) + WB * T(I,J,KB)
              TEMP = ( TEMP + T(I,J,KB) ) / 2.
              HP(I,J,N) =H(I,J,KB)+RGAS*TEMP*LOG(PSIG(KB)/P(N))/GRAV
            ELSE IF((P(N).GE.PSIG(KM)).AND.(P(N).LE.PSFC)) THEN
              TEMP = T(I,J,KM)
              HP(I,J,N) =HT(I,J)+RGAS*TEMP*LOG(PSFC/P(N))/GRAV
            ELSE IF(P(N).GT.PSFC) THEN
              TEMP = T(I,J,KBC) - TLAPSE * (H(I,J,KBC)-HT(I,J))
              HP(I,J,N) =HT(I,J)-(TEMP/TLAPSE)
     &              * ( 1.-EXP(-RGAS*TLAPSE*LOG(P(N)/PSFC)/GRAV))
!
            ENDIF
         ENDDO
         ELSE
         DO N = 1,KP
            HP(I,J,N) = -9999.0
         ENDDO
         ENDIF
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE HTSIG(T,H,P3D,PS,HT,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL    T(IM,JM,KM),H(IM,JM,KM),P3D(IM,JM,KM)
      REAL    PS(IM,JM),HT(IM,JM)
      REAL    RGAS,GRAV
      INTEGER I,J,K
      REAL    TBAR
!
      RGAS   =287.04
      GRAV   =  9.80616
!
      DO J=1,JM
      DO I=1,IM
         IF(PS(I,J).GT.-9995.0) THEN
         H(I,J,KM)=HT(I,J)+RGAS/GRAV*T(I,J,KM)*LOG(PS(I,J)/P3D(I,J,KM))
         ELSE
         H(I,J,KM)=-9999.0
         ENDIF
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM
      DO I=1,IM
         IF(H(I,J,K+1).GT.-9995.0) THEN
         TBAR = 0.5*( T(I,J,K)+T(I,J,K+1) )
         H(I,J,K)=H(I,J,K+1)+RGAS/GRAV*TBAR*LOG(P3D(I,J,K+1)/P3D(I,J,K))
         ELSE
         H(I,J,K)=-9999.0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HUMID1(T,Q,PS,PT,SIGMA,NI,NJ,NK)
      IMPLICIT NONE
      REAL    TR
      PARAMETER(TR=1./273.16)
      INTEGER NI,NJ,NK
      REAL    T(NI,NJ,NK),Q(NI,NJ,NK),PS,SIGMA(NK),PT
      INTEGER I,J,K
      REAL    P,HL,SATVP,QS
!
!  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
      DO 1 I=1,NI
      DO 1 J=1,NJ
      DO 1 K=1,NK
      P=(PT+SIGMA(K)*PS)*10.                    ! PRESSURE AT LEVEL K
      HL=597.3-.566*(T(I,J,K)-273.16)           ! LATENT HEAT OF EVAP.
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K))) ! SATURATION VAP PRESS.
      QS=.622*SATVP/(P-SATVP)                   ! SAT. MIXING RATIO
!AS MIXING RATIO MOD
!AS   Q(I,J,K)=Q(I,J,K)/QS/(1.-Q(I,J,K))        !CONVERTS SP.HUM TO MIXING RATIO
      Q(I,J,K)=amax1(Q(I,J,K)/QS,0.0)           !ALREADY MIXING RATIO
    1 CONTINUE
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HUMID1FV(T,Q,P3D,NI,NJ,NK)
      IMPLICIT NONE
      REAL    TR
      PARAMETER(TR=1./273.16)
      INTEGER NI,NJ,NK
      REAL    T(NI,NJ,NK),Q(NI,NJ,NK),P3D(NI,NJ,NK)
      INTEGER I,J,K
      REAL    HL,SATVP,QS
!
!  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
      DO 1 I=1,NI
      DO 1 J=1,NJ
      DO 1 K=1,NK
      IF(P3D(I,J,K).GT.-9990.) THEN
      HL=597.3-.566*(T(I,J,K)-273.16)              ! LATENT HEAT OF EVAP.
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))    ! SATURATION VAP PRESS.
      QS=.622*SATVP/(P3D(I,J,K)-SATVP)             ! SAT. MIXING RATIO
!AS MIXING RATIO MOD
!AS   Q(I,J,K)=amax1(Q(I,J,K)/QS/(1.-Q(I,J,K)),0.) ! CONVERTS SP.HUM TO
!                                                  ! MIXING RATIO
      Q(I,J,K)=amax1(Q(I,J,K)/QS,0.0)              !ALREADY MIXING RATIO
      ELSE
         Q(I,J,K)=-9999.
      ENDIF
    1 CONTINUE
      RETURN
      END
!
      SUBROUTINE HUMID2(T,Q,PS,PT,SIGMA,NI,NJ,NK)
      IMPLICIT NONE
      REAL    TR
      PARAMETER(TR=1./273.16)
      INTEGER NI,NJ,NK
      REAL    T(NI,NJ,NK),Q(NI,NJ,NK),PS(NI,NJ),SIGMA(NK),PT
      INTEGER I,J,K
      REAL    P,HL,SATVP,QS
!
!  THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
!
      DO 2 I=1,NI
      DO 2 J=1,NJ
      DO 2 K=1,NK
      P=(PT+SIGMA(K)*PS(I,J))*10.
      HL=597.3-.566*(T(I,J,K)-273.16)
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))
!     IF(P.LT.300.) P = 300.           ! GB MOD: KEEP Q SMALL FOR P<300
      QS=.622*SATVP/(P-SATVP)
!     IF (Q(I,J,K).LT.0.1) Q(I,J,K)=0.1
      Q(I,J,K)=amax1(Q(I,J,K)*QS,0.0)
    2 CONTINUE
!
      RETURN
      END
!
      SUBROUTINE HUMID2FV(T,Q,PS,PT,SIGMA,NI,NJ,NK)
      IMPLICIT NONE
      REAL    TR
      PARAMETER(TR=1./273.16)
      INTEGER NI,NJ,NK
      REAL    T(NI,NJ,NK),Q(NI,NJ,NK),PS(NI,NJ),SIGMA(NK),PT
      INTEGER I,J,K
      REAL    P,HL,SATVP,QS
!
!  THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
!
      DO 2 I=1,NI
      DO 2 J=1,NJ
      DO 2 K=1,NK
      P=(PT+SIGMA(K)*PS(I,J))*10.
      HL=597.3-.566*(T(I,J,K)-273.16)
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))
!     IF(P.LT.300.) P = 300.           ! GB MOD: KEEP Q SMALL FOR P<300
      QS=.622*SATVP/(P-SATVP)
!     IF (Q(I,J,K).LT.0.1) Q(I,J,K)=0.1
      Q(I,J,K)=amax1(Q(I,J,K)*QS,0.0)
    2 CONTINUE
!
      RETURN
      END
!
! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE HYDROST(H,T,PHIS,PS,PT,SIGMAF,SIGMAH,DSIGMA,NI,NJ,NK)
      IMPLICIT NONE
!
! ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
! THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
! HEIGHT IS COMPUTED IN THE RCM MODEL.
!
      REAL    RGAS,GRAV
      PARAMETER(RGAS=287.04, GRAV=9.8)
      INTEGER NI,NJ,NK
      REAL    H(NI,NJ,NK),T(NI,NJ,NK),PHIS(NI,NJ),PS(NI,NJ)
      REAL    SIGMAF(NK+1),SIGMAH(NK),DSIGMA(NK),PT
      REAL    RG,PF,TBAR
      INTEGER I,J,K,K1,K2
!
      RG=RGAS/GRAV
      DO 1 K=1,NK
    1 DSIGMA(K)=SIGMAF(K+1)-SIGMAF(K)
!
! SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
! DEFINED AT CROSS POINTS AND AT HALF LEVELS.
!
      DO 2 I=1,NI
      DO 2 J=1,NJ
      PF=PT/PS(I,J)
      H(I,J,NK)=PHIS(I,J)+RG*T(I,J,NK)*LOG((1.+PF)/(SIGMAH(NK)+PF))
      DO 2 K2=1,NK-1
      K=NK-K2
      K1=K+1
      TBAR=(T(I,J,K)*DSIGMA(K)+T(I,J,K1)*DSIGMA(K1))/(DSIGMA(K)
     A     +DSIGMA(K1))
      H(I,J,K)=H(I,J,K1)+RG*TBAR*LOG((SIGMAH(K1)+PF)/
     A         (SIGMAH(K)+PF))
    2 CONTINUE
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INITDATE
      IMPLICIT NONE
!
      INTEGER MDATE
      COMMON /DATENUM/MDATE(299500)
      INTEGER MBASE,NBASE,NREC,NYEAR,MON,NDAY,I,M
!
      NREC=0
      DO NYEAR=1941,2145
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
               IF(MOD(NYEAR,100).EQ.0) NDAY=NDAY-1
               IF(MOD(NYEAR,400).EQ.0) NDAY=NDAY+1
            ELSE
               NDAY=28
            ENDIF
         ENDIF
         NBASE = MBASE
         DO I=1,NDAY
            NBASE = NBASE+100
            DO M=1,4
               NREC=NREC+1
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
      WRITE(*,*) 'NREC = ',NREC
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INTGTB( PA, ZA, TLAYER, ZRCM
     A                 , TP, ZP, SCCM, NI, NJ, NLEV1 )
      IMPLICIT NONE
 
      REAL    RGAS,RLAPSE,GRAV,RGAS2,B1
      PARAMETER(RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER(RGAS2=RGAS/2., B1=GRAV/RLAPSE )
 
      INTEGER NI,NJ,NLEV1
      REAL    PA(NI,NJ), ZA(NI,NJ), TLAYER(NI,NJ), ZRCM(NI,NJ)
      REAL    TP(NI,NJ,NLEV1), ZP(NI,NJ,NLEV1), SCCM(NLEV1)
 
      INTEGER I,J,K
      INTEGER KT,KB
 
!  INTGTB CALCULATES ALL VARIABLES NEEDED TO COMPUTE P* ON THE RCM
!        TOPOGRAPHY.  THE MEAN TEMPERATURE IN THE LAYER BETWEEN
!        THE TOPOGRAPHY AND THE PRESSURE LEVEL ABOVE IS CALULATED
!        BY LINEARLY INTERPOLATING WITH HEIGHT THE TEMPS ON
!        PRESSURE LEVELS.
!        INPUT:    TP        TEMPS ON ECMWF PRESSURE LEVELS
!                  ZP        HEIGHTS OF ECMWF PRESSURE LEVELS
!                  ZRCM      RCM TOPOGRAPHY
!                  SCCM      ECMWF PRESSURE LEVELS (DIVIDED BY 1000.)
!        OUTPUT:   TLAYER    MEAN LAYER TEMP ABOVE RCM SURFACE
!                  PA        PRESSURE AT TOP OF LAYER
!                  ZA        HEIGHT AT PRESSURE PA
!
      DO 40 I=1,NI
      DO 40 J=1,NJ
 
      KT = 0
      DO 30 K=1,NLEV1-1
      IF(ZRCM(I,J).LE.ZP(I,J,NLEV1+1-K).AND.
     &   ZRCM(I,J).GT.ZP(I,J,NLEV1-K)) KT=K
   30 CONTINUE
      KB = KT + 1
 
      IF(KT.NE.0) THEN
         TLAYER(I,J) =
     &    ( TP(I,J,NLEV1+1-KT) * (ZRCM(I,J)-ZP(I,J,NLEV1+1-KB))
     &    + TP(I,J,NLEV1+1-KB) * (ZP(I,J,NLEV1+1-KT)-ZRCM(I,J)) )
     &                / (ZP(I,J,NLEV1+1-KT)-ZP(I,J,NLEV1+1-KB))
         TLAYER(I,J) = ( TP(I,J,NLEV1+1-KT)+TLAYER(I,J) ) / 2.
         ZA(I,J) = ZP(I,J,NLEV1+1-KT)
         PA(I,J) = 100. * SCCM(KT)
      ELSE
         TLAYER(I,J) = TP(I,J,1)
         ZA(I,J) = ZP(I,J,1)
         PA(I,J) = 100.
      ENDIF
 
   40 CONTINUE
 
!     PRINT *, 'ZRCM, ZP(6)   =', ZRCM(5,5), ZP(5,5,NLEV1+1-6)
!     PRINT *, '      TP(6)   =',            TP(5,5,NLEV1+1-6)
!     PRINT *, 'TLAYER, ZA, PA =', TLAYER(5,5), ZA(5,5), PA(5,5)
 
      RETURN
      END
!
      SUBROUTINE INTLIN(FP,F,PS,P3D,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PS(IM,JM),P3D(IM,JM,KM)
      REAL    SIG(61),P(KP)
      INTEGER I,J,K,N
      INTEGER K1,K1P
      REAL    SIGP,WP,W1
!
!  INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
!        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
!        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.

      DO J=1,JM
      DO I=1,IM
         IF(PS(I,J).GT.-9995.0) THEN
         DO K=1,KM
            SIG(K) = P3D(I,J,K)/PS(I,J)
         ENDDO
         DO N=1,KP
            SIGP = P(N) / PS(I,J)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = (SIGP-SIG(K1))/(SIG(K1P)-SIG(K1))
              W1  = 1.-WP
              FP(I,J,N)  = W1*F(I,J,K1)+WP*F(I,J,K1P)
            ELSE IF(SIGP.GE.SIG(KM)) THEN
              FP(I,J,N)  = F(I,J,KM)
            ENDIF
         ENDDO
         ELSE
         DO N=1,KP
            FP(I,J,N)  = -9999.0
         ENDDO
         ENDIF
      ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE INTLOG(FP,F,PS,P3D,IM,JM,KM,P,KP)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL    FP(IM,JM,KP),F(IM,JM,KM)
      REAL    PS(IM,JM),P3D(IM,JM,KM)
      REAL    SIG(61),P(KP)
      INTEGER I,J,K,N
      INTEGER K1,K1P,KBC
      REAL    SIGP,WP,W1
      REAL    BLTOP,RGAS,TLAPSE,GRAV
!
      RGAS   =287.04
      GRAV   =  9.80616
      BLTOP  =  0.96
      TLAPSE = -6.5E-3
!
!  INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATURES IN THE LAYER.

!
!** FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      DO J=1,JM
      DO I=1,IM
         IF(PS(I,J).GT.-9995.0) THEN
         DO K=1,KM
            SIG(K) = P3D(I,J,K)/PS(I,J)
            IF(SIG(K).LT.BLTOP) KBC = K
         ENDDO
         DO N=1,KP
            SIGP = P(N) / PS(I,J)
            K1=0
            DO K=1,KM
              IF (SIGP.GT.SIG(K)) K1=K
            ENDDO
            IF(SIGP.LE.SIG(1)) THEN
              FP(I,J,N) = F(I,J,1)
            ELSE IF((SIGP.GT.SIG(1)).AND.(SIGP.LT.SIG(KM))) THEN
              K1P = K1 + 1
              WP  = LOG(SIGP/SIG(K1)) / LOG(SIG(K1P)/SIG(K1))
              W1  = 1. - WP
              FP(I,J,N)= W1*F(I,J,K1) + WP*F(I,J,K1P)
            ELSE IF((SIGP.GE.SIG(KM)).AND.(SIGP.LE.1.))THEN
              FP(I,J,N)= F(I,J,KM)
            ELSE IF(SIGP.GT.1.) THEN
              FP(I,J,N) = F(I,J,KBC)
     &        * EXP(-RGAS*TLAPSE*LOG(SIGP/SIG(KBC))/GRAV)
!                  ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
            ENDIF
         ENDDO
         ELSE
         DO N=1,KP
            FP(I,J,N) = -9999.0
         ENDDO
         ENDIF
      ENDDO
      ENDDO

      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INTPSN( PSRCM, ZRCM, PA, ZA, TLAYER,PT,NI,NJ)
      IMPLICIT NONE
 
      REAL    RGAS,RLAPSE,GRAV
      PARAMETER(RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8)
      INTEGER NI,NJ
      REAL    PSRCM(NI,NJ), ZRCM(NI,NJ), PT
      REAL    PA(NI,NJ), ZA(NI,NJ), TLAYER(NI,NJ)
      REAL    RL2,GDRM,TB
      INTEGER I,J
 
!  EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
!        USE TLAYER CALCULATED IN INTGTB.
!        PSRCM = SURFACE PRESSURE - PTOP
!
      RL2=RLAPSE/2.
      GDRM=-GRAV/RGAS
      DO 1 I=1,NI
      DO 1 J=1,NJ
!!!   TB=T(I,J)+RL2*(ZRCM(I,J)-ZCCM(I,J))
      TB = TLAYER(I,J)
      PSRCM(I,J) = PA(I,J) * EXP(GDRM*(ZRCM(I,J)-ZA(I,J))/TB)-PT
    1 CONTINUE
 
!     PRINT *, 'ZRCM, ZA, PA, PT =', ZRCM(5,5), ZA(5,5), PA(5,5), PT
!     PRINT *, 'TLAYER(5,5), PSRCM(5,5) = ', TLAYER(5,5), PSRCM(5,5)
 
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INTV1(FRCM,FCCM,PSRCM,SRCM,SCCM,PT,NI,NJ,KRCM,KCCM)
      IMPLICIT NONE
      REAL    RGAS,RLAPSE,GRAV,RGAS2,B1,PSCCM
      PARAMETER(RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER(RGAS2=RGAS/2., B1=GRAV/RLAPSE,PSCCM=100. )

      INTEGER NI,NJ,KRCM,KCCM
      REAL    PSRCM(NI,NJ),SRCM(KRCM),SCCM(KCCM)
     A,       FRCM(NI,NJ,KRCM), FCCM(NI,NJ,KCCM)
      REAL    PT
 
      INTEGER I,J,K,N,K1,K1P
      REAL    DP1,PT1,SC,RC,RC1

!  INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
!        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
!        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!  INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATUES IN THE LAYER.
 
      DO 20 I=1,NI
      DO 20 J=1,NJ
      DP1=PSRCM(I,J)/PSCCM
      PT1=PT/PSCCM
      DO 20 N=1,KRCM
      SC=SRCM(N)*DP1+PT1
      K1=0
      DO 10 K=1,KCCM
      IF (SC.GT.SCCM(K)) K1=K
   10 CONTINUE
! 
!  CONDITION FOR SC .LT. SCCM(1) FOLLOWS
      IF (K1.EQ.0) THEN
         FRCM(I,J,N)=FCCM(I,J,KCCM)
!
!  CONDITION FOR SCCM(1) .LT. SC .LT. SCCM(KCCM) FOLLOWS
      ELSE IF (K1.NE.KCCM) THEN
         K1P=K1+1
         RC=(SC-SCCM(K1))/(SCCM(K1)-SCCM(K1P))
         RC1=RC+1.
         FRCM(I,J,N)=RC1*FCCM(I,J,KCCM+1-K1)-RC*FCCM(I,J,KCCM+1-K1P)
!
!  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
      ELSE
         FRCM(I,J,N)=FCCM(I,J,1)
!
      ENDIF
   20 CONTINUE
 
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INTV2(FRCM,FCCM,PSRCM,SRCM,SCCM,PT,NI,NJ,KRCM,KCCM)
      IMPLICIT NONE
      REAL    RGAS,RLAPSE,GRAV,RGAS2,B1,PSCCM
      PARAMETER(RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER(RGAS2=RGAS/2., B1=GRAV/RLAPSE,PSCCM=100. )

      INTEGER NI,NJ,KRCM,KCCM
      REAL    PSRCM(NI,NJ),SRCM(KRCM),SCCM(KCCM)
     &       ,FRCM(NI,NJ,KRCM), FCCM(NI,NJ,KCCM)
      REAL    PT
 
      INTEGER I,J,K,N,K1,K1P
      REAL    DP1,PT1,SC,RC,RC1,A1

!  INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
!        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
!        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!  INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATUES IN THE LAYER.
! 
      DO 40 I=1,NI
      DO 40 J=1,NJ
      DP1=PSRCM(I,J)/PSCCM
      PT1=PT/PSCCM
      DO 40 N=1,KRCM
      SC=SRCM(N)*DP1+PT1
      K1=0
      DO 30 K=1,KCCM
      IF (SC.GT.SCCM(K)) K1=K
   30 CONTINUE
! 
!  CONDITION FOR SC .LT. SCCM(1) FOLLOWS
      IF (K1.EQ.0) THEN
         FRCM(I,J,N)=FCCM(I,J,KCCM)
!
!  CONDITION FOR SCCM(1) .LT. SC .LT. SCCM(KCCM) FOLLOWS
      ELSE IF (K1.NE.KCCM) THEN
         K1P=K1+1
         RC=LOG(SC/SCCM(K1))/LOG(SCCM(K1)/SCCM(K1P))
         RC1=RC+1.
         FRCM(I,J,N)=RC1*FCCM(I,J,KCCM+1-K1)-RC*FCCM(I,J,KCCM+1-K1P)
!
!  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
      ELSE
         A1=RGAS2*LOG(SC/SCCM(KCCM))
         FRCM(I,J,N)=FCCM(I,J,1)*(B1-A1)/(B1+A1)
!
      ENDIF
   40 CONTINUE
 
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE INTV3(FSCCM,FCCM,PSRCCM,SCCM,PTOP,NI,NJ,KCCM)
      IMPLICIT NONE
 
      REAL    PTOP
      REAL    RGAS,RLAPSE,GRAV,RGAS2,B1
      PARAMETER(RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER(RGAS2=RGAS/2., B1=GRAV/RLAPSE )

      INTEGER NI,NJ,KCCM
      REAL    FSCCM(NI,NJ),FCCM(NI,NJ,KCCM),PSRCCM(NI,NJ)
      REAL    SCCM(KCCM)
 
      INTEGER I,J,K,K1,K1P
      REAL    SC,A1,RC,RC1
 
!** INTV3 IS FOR VERTICAL INTERPOLATION OF TSCCM.  THE INTERPOLATION IS
!        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!        CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!        TWO EXTREME TEMPERATUES IN THE LAYER.
!
      DO 40 I=1,NI
      DO 40 J=1,NJ
!bug  SC=PSRCCM(I,J)/100.
      SC=(PSRCCM(I,J)+PTOP)/100.
      DO 30 K=1,KCCM-1
      IF (SC.LE.SCCM(K+1) .AND. SC.GE.SCCM(K)) K1=K
   30 CONTINUE
 
!  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
         IF(SC.GT.SCCM(KCCM)) THEN
            A1=RGAS2*LOG(SC/SCCM(KCCM))
            FSCCM(I,J)=FCCM(I,J,KCCM+1-KCCM)*(B1-A1)/(B1+A1)
            GO TO 38
         END IF
!
!  CONDITION FOR SC .LT. SCCM(KCCM) FOLLOWS
         K1P=K1+1
         RC=LOG(SC/SCCM(K1))/LOG(SCCM(K1)/SCCM(K1P))
         RC1=RC+1.
         FSCCM(I,J)=RC1*FCCM(I,J,KCCM+1-K1)-RC*FCCM(I,J,KCCM+1-K1P)
!
   38 CONTINUE
   40 CONTINUE
 
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE JULIAN( MDATE, NYRP, NMOP, WT )
      IMPLICIT NONE
 
      INTEGER MDATE,NYRP,NMOP
      REAL    WT
      INTEGER LENMON(12), MIDMON(12), JULMID(12), JPREV(12)
 
      DATA    LENMON / 31,28,31,30,31,30,31,31,30,31,30,31 /
      DATA    MIDMON / 16,14,16,15,16,15,16,16,15,16,15,16 /
      INTEGER IDATE,IYR,IMO,IDAY,ILEAP,J,JULDAY,NYR,NMO
      REAL    FNUMER,FDENOM
 
! ******           INITIALIZE NMOP, NYRP
      NMOP = 1
      NYRP = 0
 
      IDATE = MDATE / 100
      IYR = IDATE / 10000
      IMO = ( IDATE - IYR*10000 ) / 100
      IDAY = MOD( IDATE, 100 )

      ILEAP = MOD( IYR, 4 )
      LENMON(2) = 28
      IF(ILEAP.EQ.0) LENMON(2) = 29
 
      JPREV(1) = 0
      DO 10 J=2,12
      JPREV(J)  = JPREV(J-1) + LENMON(J-1)
   10 CONTINUE
      DO 15 J=1,12
      JULMID(J) = JPREV(J) + MIDMON(J)
   15 CONTINUE
      JULDAY = IDAY + JPREV(IMO)
 
!     PRINT *, 'MDATE, IYR, IMO, IDAY, JULDAY = '
!    A       ,  MDATE, IYR, IMO, IDAY, JULDAY
 
      DO 50  NYR=1948, 2145 !94
      DO 50  NMO=1,12
 
      IF( (NYR.EQ.IYR) .AND. (JULMID(NMO).GT.JULDAY) ) GOTO 150
      IF  (NYR.GT.IYR) GOTO 150
 
      NMOP = NMO
      NYRP = NYR
 
   50 CONTINUE
 
  150 CONTINUE
      FNUMER = FLOAT( JULDAY      - JULMID(NMOP) )
      IF(FNUMER.LT.0.) FNUMER = FNUMER + 365.
      FDENOM = FLOAT( JULMID(NMO) - JULMID(NMOP) )
      IF(FDENOM.LE.0.) FDENOM = FDENOM + 365.
      WT = FNUMER / FDENOM
 
!     PRINT *, 'JULMID(NMOP), JULDAY, JULMID(NMO), WT ='
!    A       ,  JULMID(NMOP), JULDAY, JULMID(NMO), WT
 
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE MKSST(TSCCM,SST1,SST2,TOPOGM,XLANDU,JX,IY,NYRP,NMOP,WT)
      IMPLICIT NONE
 
      INTEGER NYRP,NMOP
      INTEGER JX,IY
      REAL    XLANDU(JX,IY)
      REAL    TSCCM(JX,IY), TOPOGM(JX,IY)
      REAL    SST1(JX,IY), SST2(JX,IY)
      REAL    WT
      INTEGER LON,LAT,NDAY,NMO,NYEAR
      INTEGER I,J
 
! ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
      DO 5 LON=1,JX
      DO 5 LAT=1,IY
      SST1(LON,LAT) = 0.
      SST2(LON,LAT) = 0.
    5 CONTINUE
 
      IF(NYRP.EQ.0) THEN
         WT = 1.
         GOTO 15
      ENDIF
 
! ******           READ IN RCM MONTHLY SST DATASET
   10 READ(60,END=998) NDAY,NMO,NYEAR,((SST1(I,J),J=1,IY),I=1,JX)
      IF(NYEAR.LT.100) NYEAR=NYEAR+1900
      IF( (NYEAR.NE.NYRP) .OR. (NMO.NE.NMOP) ) GOTO 10
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
 
! ******           READ IN RCM MONTHLY SST DATASET
   15 READ(60,END=998) NDAY,NMO,NYEAR,((SST2(I,J),J=1,IY),I=1,JX)
      IF(NYEAR.LT.100) NYEAR=NYEAR+1900
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
      REWIND(60)
 
      DO 20 I=1,JX
      DO 20 J=1,IY
      IF( (TOPOGM(I,J).LE.1.) .AND.
     &    (XLANDU(I,J).GT.13.9.AND.XLANDU(I,J).LT.15.1) .AND.
     &    (SST1(I,J).GT.-900.0.AND.SST2(I,J).GT.-900.0) ) THEN
         TSCCM(I,J) = (1.-WT) * SST1(I,J) + WT * SST2(I,J)
      ENDIF
   20 CONTINUE
 
      REWIND(60)
 
      RETURN
  998 PRINT *, 'SST file is not the right one'
      STOP 12
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE MKSST2(TSCCM,SST1,SST2,TOPOGM,XLANDU,JX,IY,KDATE)
      IMPLICIT NONE
 
      INTEGER JX,IY,KDATE
      REAL    XLANDU(JX,IY)
      REAL    TSCCM(JX,IY), TOPOGM(JX,IY)
      REAL    SST1(JX,IY), SST2(JX,IY)
      integer WKDAY
      COMMON /DATEWK/ WKDAY(427+1045)
      INTEGER MDATE
      COMMON /DATENUM/MDATE(299500)
      REAL    WT
      INTEGER LON,LAT,NDAY,NMO,NYEAR
      INTEGER KDATE1,KDATE2
      INTEGER I,J,K,KS,KS1,KS2
!
      DO K=427+1045,1,-1
         IF(WKDAY(K).LE.KDATE) THEN
            KS=K
            GOTO 100
         ENDIF
      ENDDO
 100  CONTINUE
      KDATE1=WKDAY(KS)
!
      DO K=1,427+1045
         IF(WKDAY(K).GT.KDATE) THEN
            KS=K
            GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE
      KDATE2=WKDAY(KS)
      CALL FINDDATE(KS1,KDATE1*100)
      CALL FINDDATE(KS ,KDATE *100)
      CALL FINDDATE(KS2,KDATE2*100)
      WT = FLOAT(KS-KS1)/FLOAT(KS2-KS1)
 
! ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
      DO 5 LON=1,JX
      DO 5 LAT=1,IY
      SST1(LON,LAT) = 0.
      SST2(LON,LAT) = 0.
    5 CONTINUE
 
! ******           READ IN RCM MONTHLY SST DATASET
   10 READ(60,END=998) NDAY,NMO,NYEAR,((SST1(I,J),J=1,IY),I=1,JX)
      IF(NYEAR*10000+NMO*100+NDAY.NE.KDATE1) GOTO 10
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
 
! ******           READ IN RCM MONTHLY SST DATASET
!   15 READ(60,END=998) NDAY,NMO,NYEAR,((SST2(I,J),J=1,IY),I=1,JX)
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
 
      DO 20 I=1,JX
      DO 20 J=1,IY
      IF( (TOPOGM(I,J).LE.1.) .AND.
     &    (XLANDU(I,J).GT.13.9.AND.XLANDU(I,J).LT.15.1) .AND.
     &    (SST1(I,J).GT.-900.0.AND.SST2(I,J).GT.-900.0) ) THEN
         TSCCM(I,J) = (1.-WT) * SST1(I,J) + WT * SST2(I,J)
      ENDIF
   20 CONTINUE
 
      REWIND(60)
 
      RETURN
  998 PRINT *, 'SST file is not the right one'
      STOP 12
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE MKSST3(TSCCM,SST1,TOPOGM,XLANDU,JX,IY,KDATE)
      IMPLICIT NONE
 
      INTEGER JX,IY,KDATE
      REAL    XLANDU(JX,IY)
      REAL    TSCCM(JX,IY), TOPOGM(JX,IY)
      REAL    SST1(JX,IY)
      INTEGER NHOUR,NDAY,NMO,NYEAR
      INTEGER         MHOUR,MDAY,MMO,MYEAR
      INTEGER I,J
!
      NYEAR= KDATE/1000000
      NMO  =(KDATE-NYEAR*1000000)/10000
      NDAY =(KDATE-NYEAR*1000000-NMO*10000)/100
      NHOUR= KDATE-NYEAR*1000000-NMO*10000-NDAY*100
!
! ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
 
! ******           READ IN RCM 6 HOUR SST DATASET
   10 READ(60,END=998) MHOUR,MDAY,MMO,MYEAR,
     &                 ((SST1(I,J),J=1,IY),I=1,JX)
      IF(.NOT.(NYEAR.EQ.MYEAR.AND.NMO.EQ.MMO.AND.
     &         NDAY.EQ.MDAY.AND.NHOUR.EQ.MHOUR)) GOTO 10
 
      DO 20 I=1,JX
      DO 20 J=1,IY
      IF( (TOPOGM(I,J).LE.1.) .AND.
     &    (XLANDU(I,J).GT.13.9.AND.XLANDU(I,J).LT.15.1) .AND.
     &    (SST1(I,J).GT.-900.0) ) THEN
         TSCCM(I,J) = SST1(I,J)
      ENDIF
   20 CONTINUE
 
      REWIND(60)
 
      RETURN
  998 PRINT *, 'SST file is not the right one'
      STOP 12
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE P1P2(PD,PX,NI,NJ)
      IMPLICIT NONE
      INTEGER NI,NJ
      REAL    PD(NI,NJ),PX(NI,NJ)
      INTEGER NI1,NJ1
      INTEGER I,J
!
!  THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
!  ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
!  SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE I'S.
!
      NI1=NI-1
      NJ1=NJ-1
!
      DO 1 J=2,NJ1
      DO 1 I=2,NI1
    1 PD(I,J)=0.25*(PX(I,J)+PX(I-1,J)+PX(I,J-1)+PX(I-1,J-1))
!
      DO 2 I=2,NI1
      PD(I,1)=0.5*(PX(I,1)+PX(I-1,1))
    2 PD(I,NJ)=0.5*(PX(I,NJ1)+PX(I-1,NJ1))
!
      DO 3 J=2,NJ1
      PD(1,J)=0.5*(PX(1,J)+PX(1,J-1))
    3 PD(NI,J)=0.5*(PX(NI1,J)+PX(NI1,J-1))
!
      PD(1,1)=PX(1,1)
      PD(1,NJ)=PX(1,NJ1)
      PD(NI,1)=PX(NI1,1)
      PD(NI,NJ)=PX(NI1,NJ1)
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE TOP2BTM( X,NLON1,NLAT1,NLEV1 )
      IMPLICIT NONE
      INTEGER NLON1,NLAT1,NLEV1
      REAL    X(NLON1,NLAT1,NLEV1),  WORK(30)
      INTEGER I,J,K,KR
      DO I=1,NLON1
      DO J=1,NLAT1
         DO K=1,NLEV1
            WORK(K)=X(I,J,K)
         ENDDO
         DO K=1,NLEV1
            KR=NLEV1-K+1
            X(I,J,K)=WORK(KR)
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE UVROT4(U,V,DLON,DLAT,CLON,CLAT,GRIDFC,JX,IY,LL
     &                 ,POLLON,POLLAT,LGTYPE)
      IMPLICIT NONE
      INTEGER JX,IY,LL
      REAL    U(JX,IY,LL),V(JX,IY,LL),DLON(JX,IY),DLAT(JX,IY)
      REAL    CLON,CLAT,GRIDFC,POLLON,POLLAT
      CHARACTER*6 LGTYPE
      REAL    PIR180,X,XS,XC,D
      REAL    POLLAM,POLPHI,POLCPHI,POLSPHI,ZPHI,ZRLA,ZRLAP,ZARG1,ZARG2
      REAL    ZNORM,SINDEL,COSDEL,US,VS
      INTEGER I,J,L
!
! CHANGE U AND V FROM TRUE (N,E) TO MAP VALUES (X,Y)
!
!      FOR ROTATED MERCATOR PROJECTION
!UVUSVS   -   ROTATES THE TWO WINDCOMPONENTS U AND V AT POINT
!             DLON,DLAT TO THE WINDCOMPONENTS US AND VS IN A
!             ROTATED POLE GRID WHERE THE ORIGIN IS LOCATED
!             AT POLLON,POLLAT
!**   CALL  :   CALL UVUSVS(U,V,US,VS,DLON,DLAT,POLLON,POLLAT)
!**   AUTHOR:   D.MAJEWSKI
!
      PIR180=ATAN(1.)/45.
      IF(LGTYPE.EQ.'ROTMER') THEN
         IF(POLLAT.GT.0.) then
            POLLAM = POLLON + 180.
            POLPHI = 90. - POLLAT
         ELSE
            POLPHI = 90. + POLLAT
            POLLAM = POLLON
         ENDIF
         IF (POLLAM.GT.180.) POLLAM = POLLAM - 360.

         POLCPHI = COS ( PIR180*POLPHI )
         POLSPHI = SIN ( PIR180*POLPHI )

         DO J=1,IY
         DO I=1,JX
            ZPHI   = DLAT(I,J)*PIR180
            ZRLA   = DLON(I,J)*PIR180
            IF ( DLAT(I,J).GT.89.999999 ) ZRLA = 0.0
            ZRLAP  = POLLAM*PIR180 - ZRLA
            ZARG1  = POLCPHI*SIN(ZRLAP)
            ZARG2  = POLSPHI*COS(ZPHI) - POLCPHI*SIN(ZPHI)*COS(ZRLAP)
            ZNORM  = 1.0/SQRT ( ZARG1**2 + ZARG2**2 )
            SINDEL = ZARG1*ZNORM
            COSDEL = ZARG2*ZNORM
            DO L=1,LL
               US = U(I,J,L)*COSDEL - V(I,J,L)*SINDEL
               VS = U(I,J,L)*SINDEL + V(I,J,L)*COSDEL
               U(I,J,L) = US
               V(I,J,L) = VS
            ENDDO
         ENDDO
         ENDDO
      ELSE
         DO J=1,IY
         DO I=1,JX
            IF((CLON.GE.0.0.AND.DLON(I,J).GE.0.).OR.
     &         (CLON.LT.0.0.AND.DLON(I,J).LT.0.)) THEN
               X=(CLON-DLON(I,J))*PIR180*GRIDFC
            ELSE
               IF(CLON.GE.0.0) THEN
                  IF(ABS(CLON-(DLON(I,J)+360.)).LT.
     &               ABS(CLON- DLON(I,J)) ) THEN
                     X=(CLON-(DLON(I,J)+360.))*PIR180*GRIDFC
                  ELSE
                     X=(CLON-DLON(I,J))*PIR180*GRIDFC
                  ENDIF
               ELSE
                  IF(ABS(CLON-(DLON(I,J)-360.)).LT.
     &               ABS(CLON- DLON(I,J)) ) THEN
                     X=(CLON-(DLON(I,J)-360.))*PIR180*GRIDFC
                  ELSE
                     X=(CLON-DLON(I,J))*PIR180*GRIDFC
                  ENDIF
               ENDIF
            ENDIF
            XS=SIN(X)
            XC=COS(X)
            IF(CLAT.GE.0.) THEN
               DO L=1,LL
                  D=V(I,J,L)*XS+U(I,J,L)*XC
                  V(I,J,L)=V(I,J,L)*XC-U(I,J,L)*XS
                  U(I,J,L)=D
               ENDDO
            ELSE
               DO L=1,LL
                  D=-V(I,J,L)*XS+U(I,J,L)*XC
                  V(I,J,L)=V(I,J,L)*XC+U(I,J,L)*XS
                  U(I,J,L)=D
               ENDDO
            ENDIF
         ENDDO
         ENDDO
      ENDIF
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE UVROT4NX(U,V,DLON,DLAT,CLON,CLAT,GRIDFC,JX,IY,LL
     &                 ,POLLON,POLLAT,LGTYPE)
      IMPLICIT NONE
      INTEGER JX,IY,LL
      REAL    U(JX,IY,LL),V(JX,IY,LL),DLON(JX,IY),DLAT(JX,IY)
      REAL    CLON,CLAT,GRIDFC,POLLON,POLLAT
      CHARACTER*6 LGTYPE
      REAL    PIR180,X,XS,XC,D
      REAL    POLLAM,POLPHI,POLCPHI,POLSPHI,ZPHI,ZRLA,ZRLAP,ZARG1,ZARG2
      REAL    ZNORM,SINDEL,COSDEL,US,VS
      INTEGER I,J,L
!
! CHANGE U AND V FROM TRUE (N,E) TO MAP VALUES (X,Y)
!
!      FOR ROTATED MERCATOR PROJECTION
!UVUSVS   -   ROTATES THE TWO WINDCOMPONENTS U AND V AT POINT
!             DLON,DLAT TO THE WINDCOMPONENTS US AND VS IN A
!             ROTATED POLE GRID WHERE THE ORIGIN IS LOCATED
!             AT POLLON,POLLAT
!**   CALL  :   CALL UVUSVS(U,V,US,VS,DLON,DLAT,POLLON,POLLAT)
!**   AUTHOR:   D.MAJEWSKI
!
      PIR180=ATAN(1.)/45.
      IF(LGTYPE.EQ.'ROTMER') THEN
         IF(POLLAT.GT.0.) then
            POLLAM = POLLON + 180.
            POLPHI = 90. - POLLAT
         ELSE
            POLPHI = 90. + POLLAT
            POLLAM = POLLON
         ENDIF
         IF (POLLAM.GT.180.) POLLAM = POLLAM - 360.

         POLCPHI = COS ( PIR180*POLPHI )
         POLSPHI = SIN ( PIR180*POLPHI )

         DO J=1,IY
         DO I=1,JX
            ZPHI   = DLAT(I,J)*PIR180
            ZRLA   = DLON(I,J)*PIR180
            IF ( DLAT(I,J).GT.89.999999 ) ZRLA = 0.0
            ZRLAP  = POLLAM*PIR180 - ZRLA
            ZARG1  = POLCPHI*SIN(ZRLAP)
            ZARG2  = POLSPHI*COS(ZPHI) - POLCPHI*SIN(ZPHI)*COS(ZRLAP)
            ZNORM  = 1.0/SQRT ( ZARG1**2 + ZARG2**2 )
            SINDEL = ZARG1*ZNORM
            COSDEL = ZARG2*ZNORM
            DO L=1,LL
               US = U(I,J,L)*COSDEL + V(I,J,L)*SINDEL
               VS =-U(I,J,L)*SINDEL + V(I,J,L)*COSDEL
               U(I,J,L) = US
               V(I,J,L) = VS
            ENDDO
         ENDDO
         ENDDO
      ELSE
         DO J=1,IY
         DO I=1,JX
            IF((CLON.GE.0.0.AND.DLON(I,J).GE.0.).OR.
     &         (CLON.LT.0.0.AND.DLON(I,J).LT.0.)) THEN
               X=(CLON-DLON(I,J))*PIR180*GRIDFC
            ELSE
               IF(CLON.GE.0.0) THEN
                  IF(ABS(CLON-(DLON(I,J)+360.)).LT.
     &               ABS(CLON- DLON(I,J)) ) THEN
                     X=(CLON-(DLON(I,J)+360.))*PIR180*GRIDFC
                  ELSE
                     X=(CLON-DLON(I,J))*PIR180*GRIDFC
                  ENDIF
               ELSE
                  IF(ABS(CLON-(DLON(I,J)-360.)).LT.
     &               ABS(CLON- DLON(I,J)) ) THEN
                     X=(CLON-(DLON(I,J)-360.))*PIR180*GRIDFC
                  ELSE
                     X=(CLON-DLON(I,J))*PIR180*GRIDFC
                  ENDIF
               ENDIF
            ENDIF
            XS=SIN(X)
            XC=COS(X)
            IF(CLAT.GE.0.) THEN
               DO L=1,LL
                  D=U(I,J,L)*XC-V(I,J,L)*XS
                  V(I,J,L)=U(I,J,L)*XS+V(I,J,L)*XC
                  U(I,J,L)=D
               ENDDO
            ELSE
               DO L=1,LL
                  D=U(I,J,L)*XC+V(I,J,L)*XS
                  V(I,J,L)=V(I,J,L)*XC-U(I,J,L)*XS
                  U(I,J,L)=D
               ENDDO
            ENDIF
         ENDDO
         ENDDO
      ENDIF
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE WRITEF(U,V,T,Q,PX,TS,PTOP,NI,NJ,NK,IDATE)
      IMPLICIT NONE
      INTEGER IDATE
      INTEGER NI,NJ,NK
      REAL    U(NI,NJ,NK),V(NI,NJ,NK),T(NI,NJ,NK),Q(NI,NJ,NK)
     A,       PX(NI,NJ), TS(NI,NJ)
      REAL    PTOP
      INTEGER NOUTREC
      COMMON /OUTCOUNT/ NOUTREC
      INTEGER I,J,K
!
!  THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) IDATE,NI,NJ,NK
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((U(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((V(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((T(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((Q(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) ((PX(I,J)+PTOP,I=1,NI),J=1,NJ)
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) TS
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE WRITEF2(U,V,T,Q,PX,TS,PTOP,NI,NJ,NK,IDATE)
      IMPLICIT NONE
      INTEGER IDATE
      INTEGER NI,NJ,NK
      REAL    U(NI,NJ,NK),V(NI,NJ,NK),T(NI,NJ,NK),Q(NI,NJ,NK)
     &       ,PX(NI,NJ), TS(NI,NJ)
      REAL    PTOP
      INTEGER NOUTREC
      COMMON /OUTCOUNT/ NOUTREC
      INTEGER I,J,K
!
!  THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) IDATE,NI,NJ,NK
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((U(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((V(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((T(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((Q(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
!     DO K=NK,1,-1
!        NOUTREC=NOUTREC+1
!        WRITE(64,rec=NOUTREC)((C(I,J,K),I=1,NI),J=1,NJ)
!     ENDDO
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) ((PX(I,J)+PTOP,I=1,NI),J=1,NJ)
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) TS
!
      RETURN
      END
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      SUBROUTINE WRITEFS(U,V,T,Q,PX,TS,QS3,TI3,TS3,SNOW
     &                  ,PTOP,NI,NJ,NK,IDATE,LSMTYP)
      IMPLICIT NONE
      INTEGER IDATE
      CHARACTER*4 LSMTYP
      INTEGER NI,NJ,NK
      REAL    U(NI,NJ,NK),V(NI,NJ,NK),T(NI,NJ,NK),Q(NI,NJ,NK)
     A,       PX(NI,NJ), TS(NI,NJ)
      REAL    QS3(NI,NJ,4),TI3(NI,NJ,4),TS3(NI,NJ,4),SNOW(NI,NJ)
      REAL    PTOP
      INTEGER NOUTREC
      COMMON /OUTCOUNT/ NOUTREC
      INTEGER I,J,K
!
!  THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) IDATE,NI,NJ,NK
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((U(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((V(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((T(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      DO K=NK,1,-1
         NOUTREC=NOUTREC+1
         WRITE(64,rec=NOUTREC)((Q(I,J,K),I=1,NI),J=1,NJ)
      ENDDO
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) ((PX(I,J)+PTOP,I=1,NI),J=1,NJ)
      NOUTREC=NOUTREC+1
      WRITE(64,rec=NOUTREC) TS

      IF(LSMTYP.EQ.'USGS') THEN
        DO K=1,4
           NOUTREC=NOUTREC+1
           WRITE(64,rec=NOUTREC) ((QS3(I,J,K),I=1,NI),J=1,NJ)
        ENDDO
        DO K=1,4
           NOUTREC=NOUTREC+1
           WRITE(64,rec=NOUTREC) ((TI3(I,J,K),I=1,NI),J=1,NJ)
        ENDDO
        DO K=1,4
           NOUTREC=NOUTREC+1
           WRITE(64,rec=NOUTREC) ((TS3(I,J,K),I=1,NI),J=1,NJ)
        ENDDO
        NOUTREC=NOUTREC+1
        WRITE(64,rec=NOUTREC) ((SNOW(I,J),I=1,NI),J=1,NJ)
      ENDIF
!
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
      SUBROUTINE MXMN3D(VAR,CVAR,JX,IY,NP)
      implicit none
      integer JX,IY,NP
      REAL    VAR(JX,IY,NP)
      CHARACTER*2 CVAR
      integer i,j,k
      REAL    SMAX,SMIN
!
      DO K=1,NP
         SMAX=-1.E8
         SMIN= 1.E8
         DO J=1,IY
         DO I=1,JX
            IF(SMAX.LT.VAR(I,J,K)) SMAX=VAR(I,J,K)
            IF(SMIN.GT.VAR(I,J,K)) SMIN=VAR(I,J,K)
         ENDDO
         ENDDO
         WRITE(*,*) CVAR,K,SMAX,SMIN
      ENDDO
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cdc6hour(DATTYP,IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
      CHARACTER*5 DATTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=144)
      parameter (jlat=73)
      parameter (klev=13)

!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.
!
!    DATA ARRAY AND WORK ARRAY
      integer(kind=2) work(ilon,jlat,klev)
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,PSvar
      common /CDCvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        PSvar(ilon,jlat)
!bxq_
      character*21 inname
      character*35 pathaddname
      character*5  varname(7)
      data varname/'air','hgt','rhum','uwnd','vwnd','omega','pres'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet7,start,count
      integer status
      COMMON /NNOPEN/ inet7(7),start(10),count(10)
      real(kind=8)  xscl,xoff
      COMMON /EPATCH/ xscl(7),xoff(7)
      real(kind=8)  xscale,xadd
      integer kkrec,nlev,it,m,inet,ilev,i,j,k
      logical there
!
! Below in the ncopen call is the file name of the netCDF file.
! You may want to add code to read in the file name and the variable name.
!
!    OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
!fix  do kkrec=1,7
      do kkrec=1,5
         if(DATTYP.EQ.'NNRP1') then
            if(kkrec.eq.1.or.kkrec.eq.2.or.
     &         kkrec.eq.4.or.kkrec.eq.5) nlev=klev
            if(kkrec.eq.6) nlev=12
            if(kkrec.eq.3) nlev=8
            if(kkrec.eq.7) nlev=0
         else if(DATTYP.EQ.'NNRP2') then
            if(kkrec.le.6) nlev=klev
            if(kkrec.eq.7) nlev=0
         endif
         IF(IDATE.EQ.IDATE0.OR.
     &     (MOD(IDATE,100000).EQ.10100.AND.
     &      MOD(IDATE,1000000).NE.110100)) THEN
            if(kkrec.eq.1) then
               write(inname,100) NYEAR,  'air.',NYEAR 
            else if(kkrec.eq.2) then
               write(inname,100) NYEAR,  'hgt.',NYEAR
            else if(kkrec.eq.3) then
               write(inname,200) NYEAR, 'rhum.',NYEAR
            else if(kkrec.eq.4) then
               write(inname,200) NYEAR, 'uwnd.',NYEAR
            else if(kkrec.eq.5) then
               write(inname,200) NYEAR, 'vwnd.',NYEAR
            else if(kkrec.eq.6) then
               write(inname,300) NYEAR,'omega.',NYEAR
            else if(kkrec.eq.7) then
               write(inname,400) NYEAR,'pres.sfc.',NYEAR
            endif
 100        format(I4,'/',A4,I4,'.nc')
 200        format(I4,'/',A5,I4,'.nc')
 300        format(I4,'/',A6,I4,'.nc')
 400        format(I4,'/',A9,I4,'.nc')
    
            if(DATTYP.EQ.'NNRP1') then
               pathaddname = '../DATA/NNRP1/'//inname
            else if(DATTYP.EQ.'NNRP2') then
               pathaddname = '../DATA/NNRP2/'//inname
            endif
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet7(kkrec))
            status=nf_get_att_double(inet7(kkrec),5,
     &                               'scale_factor',xscl(kkrec))
            status=nf_get_att_double(inet7(kkrec),5,
     &                               'add_offset',xoff(kkrec))
            write(*,*) inet7(kkrec),pathaddname,xscl(kkrec),xoff(kkrec)
         ENDIF

         it=(NDAY-1)*4+NHOUR/6+1
         if(MONTH.EQ.2) it=it+ 31*4
         if(MONTH.EQ.3) it=it+ 59*4
         if(MONTH.EQ.4) it=it+ 90*4
         if(MONTH.EQ.5) it=it+120*4
         if(MONTH.EQ.6) it=it+151*4
         if(MONTH.EQ.7) it=it+181*4
         if(MONTH.EQ.8) it=it+212*4
         if(MONTH.EQ.9) it=it+243*4
         if(MONTH.EQ.10)it=it+273*4
         if(MONTH.EQ.11)it=it+304*4
         if(MONTH.EQ.12)it=it+334*4
         if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+4
         if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-4
         if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+4
!bxq_
         do m=1,4
            start(m)=1
         enddo
         do m=5,10
            start(m)=0
            count(m)=0
         enddo
         count(1)=ilon
         count(2)=jlat
         count(4)=1460
         if(MOD(NYEAR,4)  .EQ.0) count(4)=1464
         if(MOD(NYEAR,100).EQ.0) count(4)=1460
         if(MOD(NYEAR,400).EQ.0) count(4)=1464
         start(4)=it
         count(4)=1
         inet=inet7(kkrec)
         if(nlev.gt.0) then
            count(3)=nlev
            status=nf_get_vara_int2(inet,5,start,count,work)
            xscale = xscl(kkrec)
            xadd   = xoff(kkrec)
            do ilev = 1, nlev
               if(kkrec.eq.1) then
                  do j=1,jlat
                  do i=1,ilon
                  Tvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)*xscale+xadd
                  enddo
                  enddo
               else if(kkrec.eq.2) then
                  do j=1,jlat
                  do i=1,ilon
                  Hvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)*xscale+xadd
                  enddo
                  enddo
               else if(kkrec.eq.3) then
                  do j=1,jlat
                  do i=1,ilon
                  RHvar(i,jlat+1-j,14-ilev)=
     &                    dmin1((work(i,j,ilev)*xscale+xadd)*0.01,1.d0)
                  enddo
                  enddo
               else if(kkrec.eq.4) then
                  do j=1,jlat
                  do i=1,ilon
                  Uvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)*xscale+xadd
                  enddo
                  enddo
               else if(kkrec.eq.5) then
                  do j=1,jlat
                  do i=1,ilon
                  Vvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)*xscale+xadd
                  enddo
                  enddo
               else if(kkrec.eq.6) then
                  do j=1,jlat
                  do i=1,ilon
                  Wvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)*xscale+xadd
                  enddo
                  enddo
               endif
            enddo
         else if(nlev.eq.0)then
            count(3)=1
            status=nf_get_vara_int2(inet,5,start,count,work)
            if(kkrec.eq.7) then
               do j=1,jlat
               do i=1,ilon
                  PSvar(i,jlat+1-j) = work(i,j,1)*xscale+xadd
               enddo
               enddo
            endif
         endif
         IF(DATTYP.EQ.'NNRP1') then
! It's a pity that we have to nudge the values by the following way
            DO k=5,1,-1
               do j=1,jlat
               do i=1,ilon
                  RHvar(i,j,k) = RHvar(i,j,k+1)
               enddo
               enddo
            ENDDO

            do j=1,jlat
            do i=1,ilon
               Wvar(i,j,1) = 0.0
            enddo
            enddo
         ENDIF
      enddo
! 
      return
      end
      subroutine cdc6hour2(IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
!      CHARACTER*5 DATTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=144)
      parameter (jlat=73)
      parameter (klev=13)

!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.
!
!    DATA ARRAY AND WORK ARRAY
      include 'icbcWIN.param'
      integer(kind=2) work(III,JJJ,klev+1)
      real    lon0,lon1,lat0,lat1
      INTEGER i0,i1,j0
      COMMON /SZwindow/lon0,lon1,lat0,lat1,i0,i1,j0
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,PSvar
      common /CDCvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        PSvar(ilon,jlat)
!bxq_
      character*24 inname
      character*38 pathaddname
      character*5  varname(7)
      data varname/'air','hgt','rhum','uwnd','vwnd','omega','pres'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet7,start,count
      integer status
      COMMON /NNOPEN/ inet7(7),start(10),count(10)
      real(kind=8)  xscl,xoff
      COMMON /EPATCH/ xscl(7),xoff(7)
      real(kind=8)  xscale,xadd
      integer kkrec,nlev,it,m,inet,ilev,i,j
      logical there
      integer ii,jj
!
      IF(IDATE.EQ.IDATE0) THEN
         i0=lon0/2.5+1
         if(i0.le.0)    i0=i0+ilon
         if(i0.gt.ilon) i0=i0-ilon
         i1=lon1/2.5+1
         if(i1.le.0)   i1=i1+ilon
         if(i1.gt.ilon) i1=i1-ilon
         j0=lat0/2.5+36
      ENDIF

!
!bxq
      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
!
!fix  do kkrec=1,7
      do kkrec=1,5
         if(kkrec.le.6) nlev=klev
         if(kkrec.eq.7) nlev=0
         IF(IDATE.EQ.IDATE0.OR.
     &     (MOD(IDATE,100000).EQ.10100.AND.
     &      MOD(IDATE,1000000).NE.110100)) THEN
            if(kkrec.eq.1) then
               write(inname,100) NYEAR,  'air.WIN.',NYEAR 
            else if(kkrec.eq.2) then
               write(inname,100) NYEAR,  'hgt.WIN.',NYEAR
            else if(kkrec.eq.3) then
               write(inname,200) NYEAR, 'rhum.WIN.',NYEAR
            else if(kkrec.eq.4) then
               write(inname,200) NYEAR, 'uwnd.WIN.',NYEAR
            else if(kkrec.eq.5) then
               write(inname,200) NYEAR, 'vwnd.WIN.',NYEAR
            else if(kkrec.eq.6) then
               write(inname,300) NYEAR,'omega.WIN.',NYEAR
            else if(kkrec.eq.7) then
               write(inname,400) NYEAR,'pres.sfc.WIN.',NYEAR
            endif
 100        format(I4,'/',A8 ,I4,'.nc')
 200        format(I4,'/',A9 ,I4,'.nc')
 300        format(I4,'/',A10,I4,'.nc')
 400        format(I4,'/',A13,I4,'.nc')
    
            pathaddname = '../DATA/NNRP2/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet7(kkrec))
            status=nf_get_att_double(inet7(kkrec),5,
     &                               'scale_factor',xscl(kkrec))
            status=nf_get_att_double(inet7(kkrec),5,
     &                               'add_offset',xoff(kkrec))
            write(*,*) inet7(kkrec),pathaddname,xscl(kkrec),xoff(kkrec)
         ENDIF

         it=(NDAY-1)*4+NHOUR/6+1
         if(MONTH.EQ.2) it=it+ 31*4
         if(MONTH.EQ.3) it=it+ 59*4
         if(MONTH.EQ.4) it=it+ 90*4
         if(MONTH.EQ.5) it=it+120*4
         if(MONTH.EQ.6) it=it+151*4
         if(MONTH.EQ.7) it=it+181*4
         if(MONTH.EQ.8) it=it+212*4
         if(MONTH.EQ.9) it=it+243*4
         if(MONTH.EQ.10)it=it+273*4
         if(MONTH.EQ.11)it=it+304*4
         if(MONTH.EQ.12)it=it+334*4
         if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+4
         if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-4
         if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+4
!bxq_
         do m=1,4
            start(m)=1
         enddo
         do m=5,10
            start(m)=0
            count(m)=0
         enddo
         count(1)=III
         count(2)=JJJ
         count(4)=1460
         if(MOD(NYEAR,4)  .EQ.0) count(4)=1464
         if(MOD(NYEAR,100).EQ.0) count(4)=1460
         if(MOD(NYEAR,400).EQ.0) count(4)=1464
         start(4)=it
         count(4)=1
         inet=inet7(kkrec)
         if(nlev.gt.0) then
            count(3)=nlev+1
            status=nf_get_vara_int2(inet,5,start,count,work)
            xscale = xscl(kkrec)
            xadd   = xoff(kkrec)
            do ilev = 1, nlev
               if(kkrec.eq.1) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     Tvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     Tvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     Tvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  endif
                  enddo
                
               else if(kkrec.eq.2) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     Hvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     Hvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     Hvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  endif
                  enddo
               else if(kkrec.eq.3) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     RHvar(ii,jj,ilev) =
     &                 dmin1((work(i,j,ilev+1)*xscale+xadd)*0.01,1.d0)
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     RHvar(ii,jj,ilev) =
     &                 dmin1((work(i,j,ilev+1)*xscale+xadd)*0.01,1.d0)
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     RHvar(ii,jj,ilev) =
     &                 dmin1((work(i,j,ilev+1)*xscale+xadd)*0.01,1.d0)
                     enddo
                  endif
                  enddo
               else if(kkrec.eq.4) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     Uvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     Uvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     Uvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  endif
                  enddo
               else if(kkrec.eq.5) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     Vvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     Vvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     Vvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  endif
                  enddo
               else if(kkrec.eq.6) then
                  do j=1,JJJ
                  jj=j0+j
                  if(i0.gt.i1) then
                     do ii=i0,ilon
                        i=ii-i0+1
                     Wvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                     do ii=1,i1
                        i=ii+(ilon-i0)+1
                     Wvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  else
                     do ii=i0,i1
                        i=ii-i0+1
                     Wvar(ii,jj,ilev) = work(i,j,ilev+1)*xscale+xadd
                     enddo
                  endif
                  enddo
               endif
            enddo
         else if(nlev.eq.0)then
            count(3)=nlev
            status=nf_get_vara_int2(inet,5,start,count,work)
            if(kkrec.eq.7) then
               do j=1,JJJ
               jj=j0+j
               if(i0.gt.i1) then
                  do ii=i0,ilon
                     i=ii-i0+1
                     PSvar(ii,jj) = work(i,j,1)*xscale+xadd
                  enddo
                  do ii=1,i1
                     i=ii+(ilon-i0)+1
                     PSvar(ii,jj) = work(i,j,1)*xscale+xadd
                  enddo
               else
                  do ii=i0,i1
                     i=ii-i0+1
                     PSvar(ii,jj) = work(i,j,1)*xscale+xadd
                  enddo
               endif
               enddo
            endif
         endif
      enddo
! 
      return
      end
!-----------------------------------------------------------------------
!    MAIN CODE for reading ERA40 data
!  readnet.f
!  This file is a fortran template file designed to read the given
!  netCDF file into memory.
!-----------------------------------------------------------------------
      subroutine era6hour(DATTYP,LSMTYP,IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
      CHARACTER*5 DATTYP
      CHARACTER*4 LSMTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=144)
      parameter (jlat=73)
      parameter (klev=23)
!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.  The array 'x' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer(kind=2) work(ilon,jlat,klev),work2D(ilon,jlat)
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar,Wvar,Qsoil,TSIce,TSoil,
     & SNOWh
      common /ERAvars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev), Wvar(ilon,jlat,klev)
     &       ,        QSoil(ilon,jlat,4),TSIce(ilon,jlat,4)
     &       ,        TSoil(ilon,jlat,4),SNOWh(ilon,jlat)
!bxq_
      character*24 inname
      character*38 pathaddname
      character*5  varname(6),sarname(3,4)
      character*2  snownm
      data varname/'air','hgt','rhum','uwnd','vwnd','omega'/
      data sarname/'swvl1','istl1','stl1','swvl2','istl2','stl2',
     &             'swvl3','istl3','stl3','swvl4','istl4','stl4'/
      data snownm/'sd'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet6,isnet3,isnow,start,count
      integer status
      COMMON /ECOPEN/ inet6(6,4),isnet3(3,4,4),isnow(4)
     &               ,start(10),count(10)
      real(kind=8)  xscl,xoff,xscl_s,xoff_s,xscl_sn,xoff_sn
      COMMON /EPATCH/ xscl(6,4),xoff(6,4),xscl_s(3,4,4),xoff_s(3,4,4)
     &               ,xscl_sn(4),xoff_sn(4)
      real(kind=8)  xscale,xadd
      integer k4,l4,kkrec,it,i,j,k,inet
      logical there
!
! Below in the ncopen call is the file name of the netCDF file.
! You may want to add code to read in the file name and the variable name.
!
!    OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      IF(IDATE.lt.1957090100.or.IDATE.gt.2002083118) THEN
         write(*,*) 'ERA40 datasets is just available from'
     &             ,' 1957090100 to 2002083118'
         STOP
      ENDIF

      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      IF(IDATE.EQ.IDATE0.OR.
     &   (MOD(IDATE,100000).EQ.10100.AND.
     &    MOD(IDATE,1000000).NE.110100)) THEN
         do k4=1,4
         do kkrec=1,5
            if(kkrec.eq.1) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'air.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'air.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'air.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'air.',NYEAR 
               endif
            else if(kkrec.eq.2) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'hgt.',NYEAR 
               endif
            else if(kkrec.eq.3) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'rhum.',NYEAR 
               endif
            else if(kkrec.eq.4) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'uwnd.',NYEAR 
               endif
            else if(kkrec.eq.5) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'vwnd.',NYEAR 
               endif
            else if(kkrec.eq.6) then
               if(k4.eq.1) then
                  write(inname,300) NYEAR,'omega.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,301) NYEAR,'omega.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,302) NYEAR,'omega.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,303) NYEAR,'omega.',NYEAR 
               endif
            endif
    
            pathaddname = '../DATA/'//DATTYP//'/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'scale_factor',xscl(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'add_offset',xoff(kkrec,k4))
            write(*,*) inet6(kkrec,k4),pathaddname
     &                 ,xscl(kkrec,k4),xoff(kkrec,k4)
         enddo
         enddo

         if(LSMTYP.eq.'USGS') THEN
         do k4=1,4
         do l4=1,4
         do kkrec=1,3
            if(kkrec.eq.1) then
               if(k4.eq.1) then
                  write(inname,400) 'Qsoil_',l4
               else if(k4.eq.2) then
                  write(inname,401) 'Qsoil_',l4
               else if(k4.eq.3) then
                  write(inname,402) 'Qsoil_',l4
               else if(k4.eq.4) then
                  write(inname,403) 'Qsoil_',l4
               endif
            else if(kkrec.eq.2) then
               if(k4.eq.1) then
                  write(inname,500) 'Tice_',l4
               else if(k4.eq.2) then
                  write(inname,501) 'Tice_',l4
               else if(k4.eq.3) then
                  write(inname,502) 'Tice_',l4
               else if(k4.eq.4) then
                  write(inname,503) 'Tice_',l4
               endif
            else if(kkrec.eq.3) then
               if(k4.eq.1) then
                  write(inname,400) 'Tsoil_',l4
               else if(k4.eq.2) then
                  write(inname,401) 'Tsoil_',l4
               else if(k4.eq.3) then
                  write(inname,402) 'Tsoil_',l4
               else if(k4.eq.4) then
                  write(inname,403) 'Tsoil_',l4
               endif
            endif
    
            pathaddname = '../DATA/'//DATTYP//'/0surface/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,isnet3(kkrec,l4,k4))
            status=nf_get_att_double(isnet3(kkrec,l4,k4),4,
     &                              'scale_factor',xscl_s(kkrec,l4,k4))
            status=nf_get_att_double(isnet3(kkrec,l4,k4),4,
     &                              'add_offset',xoff_s(kkrec,l4,k4))
            write(*,*) isnet3(kkrec,l4,k4),pathaddname
     &                ,xscl_s(kkrec,l4,k4),xoff_s(kkrec,l4,k4)
         enddo
         enddo

         if(k4.eq.1) then
            write(inname,600) 'snowdpth'
         else if(k4.eq.2) then
            write(inname,601) 'snowdpth'
         else if(k4.eq.3) then
            write(inname,602) 'snowdpth'
         else if(k4.eq.4) then
            write(inname,603) 'snowdpth'
         endif
         pathaddname = '../DATA/'//DATTYP//'/0surface/'//inname
         inquire(file=pathaddname,exist=there)
         if(.not.there) then
            print*,pathaddname,' is not available'
            stop
         endif
        status=nf_open(pathaddname,nf_nowrite,isnow(k4))
        status=nf_get_att_double(isnow(k4),4,'scale_factor',xscl_sn(k4))
        status=nf_get_att_double(isnow(k4),4,'add_offset',xoff_sn(k4))
        write(*,*) isnow(k4),pathaddname,xscl_sn(k4),xoff_sn(k4)

         enddo
         ENDIF
      ENDIF
 100  format(I4,'/',A4,I4,'.00.nc')
 101  format(I4,'/',A4,I4,'.06.nc')
 102  format(I4,'/',A4,I4,'.12.nc')
 103  format(I4,'/',A4,I4,'.18.nc')
 200  format(I4,'/',A5,I4,'.00.nc')
 201  format(I4,'/',A5,I4,'.06.nc')
 202  format(I4,'/',A5,I4,'.12.nc')
 203  format(I4,'/',A5,I4,'.18.nc')
 300  format(I4,'/',A6,I4,'.00.nc')
 301  format(I4,'/',A6,I4,'.06.nc')
 302  format(I4,'/',A6,I4,'.12.nc')
 303  format(I4,'/',A6,I4,'.18.nc')
 400  format(A6,I1,'L.00.nc')
 401  format(A6,I1,'L.06.nc')
 402  format(A6,I1,'L.12.nc')
 403  format(A6,I1,'L.18.nc')
 500  format(A5,I1,'L.00.nc')
 501  format(A5,I1,'L.06.nc')
 502  format(A5,I1,'L.12.nc')
 503  format(A5,I1,'L.18.nc')
 600  format(A8,'.00.nc')
 601  format(A8,'.06.nc')
 602  format(A8,'.12.nc')
 603  format(A8,'.18.nc')

      k4=NHOUR/6+1
      it=NDAY
      if(MONTH.EQ.2) it=it+ 31
      if(MONTH.EQ.3) it=it+ 59
      if(MONTH.EQ.4) it=it+ 90
      if(MONTH.EQ.5) it=it+120
      if(MONTH.EQ.6) it=it+151
      if(MONTH.EQ.7) it=it+181
      if(MONTH.EQ.8) it=it+212
      if(MONTH.EQ.9) it=it+243
      if(MONTH.EQ.10)it=it+273
      if(MONTH.EQ.11)it=it+304
      if(MONTH.EQ.12)it=it+334
      if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+1
      if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-1
      if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+1
      do k=1,4
         start(k)=1
      enddo
      do k=5,10
         start(k)=0
         count(k)=0
      enddo
      count(1)=ilon
      count(2)=jlat
      count(3)=klev
      count(4)=365
      if(MOD(NYEAR,4)  .EQ.0) count(4)=366
      if(MOD(NYEAR,100).EQ.0) count(4)=365
      if(MOD(NYEAR,400).EQ.0) count(4)=366
      if(NYEAR.eq.2002) count(4)=243
      if(NYEAR.eq.1957) count(4)=122
      if(NYEAR.eq.1957) it=it-243
      start(4)=it
      count(4)=1
!bxq_
      do kkrec=1,5
         inet=inet6(kkrec,k4)
         status=nf_get_vara_int2(inet,5,start,count,work)
         xscale = xscl(kkrec,k4)
         xadd   = xoff(kkrec,k4)
         if(kkrec.eq.1) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Tvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
            enddo
            enddo
            enddo
         else if(kkrec.eq.2) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Hvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Hvar(i,jlat+1-j,k)=Hvar(i,jlat+1-j,k)/9.80616
            enddo
            enddo
            enddo
         else if(kkrec.eq.3) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               RHvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               RHvar(i,jlat+1-j,k)=RHvar(i,jlat+1-j,k)*0.01
!              RHvar(i,jlat+1-j,k)=amax1(RHvar(i,jlat+1-j,k),1.05)
            enddo
            enddo
            enddo
         else if(kkrec.eq.4) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Uvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
            enddo
            enddo
            enddo
         else if(kkrec.eq.5) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Vvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
            enddo
            enddo
            enddo
         else if(kkrec.eq.6) then
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Wvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
            enddo
            enddo
            enddo
         endif
      enddo

      if(LSMTYP.eq.'USGS') THEN
      k4=NHOUR/6+1
      it=NDAY
      if(MONTH.EQ.2) it=it+ 31
      if(MONTH.EQ.3) it=it+ 59
      if(MONTH.EQ.4) it=it+ 90
      if(MONTH.EQ.5) it=it+120
      if(MONTH.EQ.6) it=it+151
      if(MONTH.EQ.7) it=it+181
      if(MONTH.EQ.8) it=it+212
      if(MONTH.EQ.9) it=it+243
      if(MONTH.EQ.10)it=it+273
      if(MONTH.EQ.11)it=it+304
      if(MONTH.EQ.12)it=it+334
      if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+1
      if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-1
      if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+1
      do k=1957,NYEAR-1
         it=it+365
         if(mod(k,4).EQ.0) it=it+1
      enddo
      it=it-243

      do k=1,3
         start(k)=1
      enddo
      do k=4,10
         start(k)=0
         count(k)=0
      enddo
      count(1)=ilon
      count(2)=jlat
      start(3)=it
      count(3)=1
!bxq_
      do l4=1,4
      do kkrec=1,3
         inet=isnet3(kkrec,l4,k4)
         status=nf_get_vara_int2(inet,4,start,count,work2D)
         xscale = xscl_s(kkrec,l4,k4)
         xadd   = xoff_s(kkrec,l4,k4)
         if(kkrec.eq.1) then
            do j=1,jlat
            do i=1,ilon
               QSoil(i,jlat+1-j,l4)=work2D(i,j)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.2) then
            do j=1,jlat
            do i=1,ilon
               TSIce(i,jlat+1-j,l4)=work2D(i,j)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.3) then
            do j=1,jlat
            do i=1,ilon
               TSoil(i,jlat+1-j,l4)=work2D(i,j)*xscale+xadd
            enddo
            enddo
         endif
      enddo
      enddo
      inet=isnow(k4)
      status=nf_get_vara_int2(inet,4,start,count,work2D)
      xscale = xscl_sn(k4)
      xadd   = xoff_sn(k4)
      do j=1,jlat
      do i=1,ilon
         SNOWh(i,jlat+1-j)=work2D(i,j)*xscale+xadd
      enddo
      enddo
      ENDIF
! 
      return
      end
      subroutine EIN756HOUR(DATTYP,IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
      CHARACTER*5 DATTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=480)
      parameter (jlat=241)
      parameter (klev=23)
!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.  The array 'x' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer(kind=2) work(ilon,jlat,37)
      COMMON /AAA/ work
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN75vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev)
!bxq_
      character*24 inname
      character*38 pathaddname
      character*1  varname(5)
      data varname/'t','z','r','u','v'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet6,start,count
      integer status
      COMMON /ECINOPEN/ inet6(5,4),start(10),count(10)
      real(kind=8)  xscl,xoff
      COMMON /EINPATCH/ xscl(5,4),xoff(5,4)
      real(kind=8)  xscale,xadd
      integer k4,kkrec,it,i,j,k,inet
      logical there
!
! Below in the ncopen call is the file name of the netCDF file.
! You may want to add code to read in the file name and the variable name.
!
!    OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      IF(IDATE.lt.1989010100.or.IDATE.gt.2007123118) THEN
         write(*,*) 'EIN75 datasets is just available from'
     &             ,' 1989010100 to 2007123118'
         STOP
      ENDIF

      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      IF(IDATE.EQ.IDATE0.OR.
     &   (MOD(IDATE,100000).EQ.10100.AND.
     &    MOD(IDATE,1000000).NE.110100)) THEN
         do k4=1,4
         do kkrec=1,5
            if(kkrec.eq.1) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'air.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'air.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'air.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'air.',NYEAR 
               endif
            else if(kkrec.eq.2) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'hgt.',NYEAR 
               endif
            else if(kkrec.eq.3) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'rhum.',NYEAR 
               endif
            else if(kkrec.eq.4) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'uwnd.',NYEAR 
               endif
            else if(kkrec.eq.5) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'vwnd.',NYEAR 
               endif
            endif
    
            pathaddname = '../DATA/'//DATTYP//'/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'scale_factor',xscl(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'add_offset',xoff(kkrec,k4))
            write(*,*) inet6(kkrec,k4),pathaddname
     &                 ,xscl(kkrec,k4),xoff(kkrec,k4)
         enddo
         enddo

      ENDIF
 100  format(I4,'/',A4,I4,'.00.nc')
 101  format(I4,'/',A4,I4,'.06.nc')
 102  format(I4,'/',A4,I4,'.12.nc')
 103  format(I4,'/',A4,I4,'.18.nc')
 200  format(I4,'/',A5,I4,'.00.nc')
 201  format(I4,'/',A5,I4,'.06.nc')
 202  format(I4,'/',A5,I4,'.12.nc')
 203  format(I4,'/',A5,I4,'.18.nc')

      k4=NHOUR/6+1
      it=NDAY
      if(MONTH.EQ.2) it=it+ 31
      if(MONTH.EQ.3) it=it+ 59
      if(MONTH.EQ.4) it=it+ 90
      if(MONTH.EQ.5) it=it+120
      if(MONTH.EQ.6) it=it+151
      if(MONTH.EQ.7) it=it+181
      if(MONTH.EQ.8) it=it+212
      if(MONTH.EQ.9) it=it+243
      if(MONTH.EQ.10)it=it+273
      if(MONTH.EQ.11)it=it+304
      if(MONTH.EQ.12)it=it+334
      if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+1
      if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-1
      if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+1
      do k=1,4
         start(k)=1
      enddo
      do k=5,10
         start(k)=0
         count(k)=0
      enddo
      count(1)=ilon
      count(2)=jlat
      count(3)=37
      count(4)=365
      if(MOD(NYEAR,4)  .EQ.0) count(4)=366
      if(MOD(NYEAR,100).EQ.0) count(4)=365
      if(MOD(NYEAR,400).EQ.0) count(4)=366
      start(4)=it
      count(4)=1
!bxq_
      do kkrec=1,5
         inet=inet6(kkrec,k4)
         status=nf_get_vara_int2(inet,5,start,count,work)
         xscale = xscl(kkrec,k4)
         xadd   = xoff(kkrec,k4)
         if(kkrec.eq.1) then
            do j=1,jlat
            do i=1,ilon
!              Tvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Tvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Tvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Tvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Tvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Tvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Tvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Tvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Tvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Tvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Tvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Tvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Tvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Tvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Tvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Tvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Tvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Tvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Tvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Tvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Tvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Tvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Tvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Tvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.2) then
            do j=1,jlat
            do i=1,ilon
!              Hvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Hvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Hvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Hvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Hvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Hvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Hvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Hvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Hvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Hvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Hvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Hvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Hvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Hvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Hvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Hvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Hvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Hvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Hvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Hvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Hvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Hvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Hvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Hvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Hvar(i,j,k)=Hvar(i,j,k)/9.80616
            enddo
            enddo
            enddo
         else if(kkrec.eq.3) then
            do j=1,jlat
            do i=1,ilon
!              RHvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               RHvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               RHvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               RHvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               RHvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               RHvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               RHvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               RHvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               RHvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               RHvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               RHvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               RHvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               RHvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               RHvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               RHvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               RHvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               RHvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               RHvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               RHvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               RHvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               RHvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               RHvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               RHvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               RHvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,23
            do j=1,jlat
            do i=1,ilon
               RHvar(i,j,k)=amax1(RHvar(i,j,k)*0.01,0.00)
            enddo
            enddo
            enddo
         else if(kkrec.eq.4) then
            do j=1,jlat
            do i=1,ilon
!              Uvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Uvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Uvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Uvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Uvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Uvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Uvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Uvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Uvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Uvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Uvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Uvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Uvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Uvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Uvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Uvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Uvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Uvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Uvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Uvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Uvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Uvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Uvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Uvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.5) then
            do j=1,jlat
            do i=1,ilon
!              Vvar(i,jlat+1-j,k) =work(i,j,k)*xscale+xadd
               Vvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Vvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Vvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Vvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Vvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Vvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Vvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Vvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Vvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Vvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Vvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Vvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Vvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Vvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Vvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Vvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Vvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Vvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Vvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Vvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Vvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Vvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Vvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         endif
      enddo
! 
      return
      end
      subroutine EIN156HOUR(DATTYP,IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
      CHARACTER*5 DATTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=240)
      parameter (jlat=121)
      parameter (klev=23)
!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.  The array 'x' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer(kind=2) work(ilon,jlat,37)
      COMMON /AAA/ work
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN15vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev)
!bxq_
      character*24 inname
      character*38 pathaddname
      character*1  varname(5)
      data varname/'t','z','r','u','v'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet6,start,count
      integer status
      COMMON /ECINOPEN/ inet6(5,4),start(10),count(10)
      real(kind=8)  xscl,xoff
      COMMON /EINPATCH/ xscl(5,4),xoff(5,4)
      real(kind=8)  xscale,xadd
      integer k4,kkrec,it,i,j,k,inet
      logical there
!
! Below in the ncopen call is the file name of the netCDF file.
! You may want to add code to read in the file name and the variable name.
!
!    OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      IF(IDATE.lt.1989010100.or.IDATE.gt.2009053118) THEN
         write(*,*) 'EIN15 datasets is just available from'
     &             ,' 1989010100 to 2009053118'
         STOP
      ENDIF

      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      IF(IDATE.EQ.IDATE0.OR.
     &   (MOD(IDATE,100000).EQ.10100.AND.
     &    MOD(IDATE,1000000).NE.110100)) THEN
         do k4=1,4
         do kkrec=1,5
            if(kkrec.eq.1) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'air.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'air.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'air.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'air.',NYEAR 
               endif
            else if(kkrec.eq.2) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'hgt.',NYEAR 
               endif
            else if(kkrec.eq.3) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'rhum.',NYEAR 
               endif
            else if(kkrec.eq.4) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'uwnd.',NYEAR 
               endif
            else if(kkrec.eq.5) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'vwnd.',NYEAR 
               endif
            endif
    
            pathaddname = '../DATA/'//DATTYP//'/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'scale_factor',xscl(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'add_offset',xoff(kkrec,k4))
            write(*,*) inet6(kkrec,k4),pathaddname
     &                 ,xscl(kkrec,k4),xoff(kkrec,k4)
         enddo
         enddo

      ENDIF
 100  format(I4,'/',A4,I4,'.00.nc')
 101  format(I4,'/',A4,I4,'.06.nc')
 102  format(I4,'/',A4,I4,'.12.nc')
 103  format(I4,'/',A4,I4,'.18.nc')
 200  format(I4,'/',A5,I4,'.00.nc')
 201  format(I4,'/',A5,I4,'.06.nc')
 202  format(I4,'/',A5,I4,'.12.nc')
 203  format(I4,'/',A5,I4,'.18.nc')

      k4=NHOUR/6+1
      it=NDAY
      if(MONTH.EQ.2) it=it+ 31
      if(MONTH.EQ.3) it=it+ 59
      if(MONTH.EQ.4) it=it+ 90
      if(MONTH.EQ.5) it=it+120
      if(MONTH.EQ.6) it=it+151
      if(MONTH.EQ.7) it=it+181
      if(MONTH.EQ.8) it=it+212
      if(MONTH.EQ.9) it=it+243
      if(MONTH.EQ.10)it=it+273
      if(MONTH.EQ.11)it=it+304
      if(MONTH.EQ.12)it=it+334
      if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+1
      if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-1
      if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+1
      do k=1,4
         start(k)=1
      enddo
      do k=5,10
         start(k)=0
         count(k)=0
      enddo
      count(1)=ilon
      count(2)=jlat
      count(3)=37
      count(4)=365
      if(MOD(NYEAR,4)  .EQ.0) count(4)=366
      if(MOD(NYEAR,100).EQ.0) count(4)=365
      if(MOD(NYEAR,400).EQ.0) count(4)=366
      start(4)=it
      count(4)=1
!bxq_
      do kkrec=1,5
         inet=inet6(kkrec,k4)
         status=nf_get_vara_int2(inet,5,start,count,work)
         xscale = xscl(kkrec,k4)
         xadd   = xoff(kkrec,k4)
         if(kkrec.eq.1) then
            do j=1,jlat
            do i=1,ilon
!              Tvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Tvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Tvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Tvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Tvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Tvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Tvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Tvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Tvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Tvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Tvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Tvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Tvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Tvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Tvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Tvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Tvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Tvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Tvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Tvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Tvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Tvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Tvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Tvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.2) then
            do j=1,jlat
            do i=1,ilon
!              Hvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Hvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Hvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Hvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Hvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Hvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Hvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Hvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Hvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Hvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Hvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Hvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Hvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Hvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Hvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Hvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Hvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Hvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Hvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Hvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Hvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Hvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Hvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Hvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Hvar(i,j,k)=Hvar(i,j,k)/9.80616
            enddo
            enddo
            enddo
         else if(kkrec.eq.3) then
            do j=1,jlat
            do i=1,ilon
!              RHvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               RHvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               RHvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               RHvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               RHvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               RHvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               RHvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               RHvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               RHvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               RHvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               RHvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               RHvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               RHvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               RHvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               RHvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               RHvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               RHvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               RHvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               RHvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               RHvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               RHvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               RHvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               RHvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               RHvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,23
            do j=1,jlat
            do i=1,ilon
               RHvar(i,j,k)=amax1(RHvar(i,j,k)*0.01,0.00)
            enddo
            enddo
            enddo
         else if(kkrec.eq.4) then
            do j=1,jlat
            do i=1,ilon
!              Uvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Uvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Uvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Uvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Uvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Uvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Uvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Uvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Uvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Uvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Uvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Uvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Uvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Uvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Uvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Uvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Uvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Uvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Uvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Uvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Uvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Uvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Uvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Uvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.5) then
            do j=1,jlat
            do i=1,ilon
!              Vvar(i,jlat+1-j,k) =work(i,j,k)*xscale+xadd
               Vvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Vvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Vvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Vvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Vvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Vvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Vvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Vvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Vvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Vvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Vvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Vvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Vvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Vvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Vvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Vvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Vvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Vvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Vvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Vvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Vvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Vvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Vvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         endif
      enddo
! 
      return
      end
      subroutine EIN256HOUR(DATTYP,IDATE,IDATE0)
      implicit none
      include 'netcdf.inc'
      CHARACTER*5 DATTYP
      INTEGER IDATE,IDATE0
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the netCDF file.
!
      integer ilon,jlat,klev
      parameter (ilon=144)
      parameter (jlat=73)
      parameter (klev=23)
!
! The data are packed into short integers (INTEGER*2).  The array work
! will be used to hold the packed integers.  The array 'x' will contain
! the unpacked data.
!
!    DATA ARRAY AND WORK ARRAY
      integer(kind=2) work(ilon,jlat,37)
      COMMON /AAA/ work
!bxq
      real(kind=4)  Tvar,Hvar,RHvar,Uvar,Vvar
      common /EIN25vars/ Tvar(ilon,jlat,klev), Hvar(ilon,jlat,klev)
     &       ,        RHvar(ilon,jlat,klev), Uvar(ilon,jlat,klev)
     &       ,         Vvar(ilon,jlat,klev)
!bxq_
      character*24 inname
      character*38 pathaddname
      character*1  varname(5)
      data varname/'t','z','r','u','v'/
      INTEGER NYEAR,MONTH,NDAY,NHOUR
      integer inet6,start,count
      integer status
      COMMON /ECINOPEN/ inet6(5,4),start(10),count(10)
      real(kind=8)  xscl,xoff
      COMMON /EINPATCH/ xscl(5,4),xoff(5,4)
      real(kind=8)  xscale,xadd
      integer k4,kkrec,it,i,j,k,inet
      logical there
!
! Below in the ncopen call is the file name of the netCDF file.
! You may want to add code to read in the file name and the variable name.
!
!    OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      IF(IDATE.lt.1989010100.or.IDATE.gt.1998123118) THEN
         write(*,*) 'EIN25 datasets is just available from'
     &             ,' 1989010100 to 1998123118'
         STOP
      ENDIF

      NYEAR=IDATE/1000000
      MONTH=IDATE/10000-NYEAR*100
      NDAY =IDATE/100-NYEAR*10000-MONTH*100
      NHOUR=IDATE-NYEAR*1000000-MONTH*10000-NDAY*100
      IF(IDATE.EQ.IDATE0.OR.
     &   (MOD(IDATE,100000).EQ.10100.AND.
     &    MOD(IDATE,1000000).NE.110100)) THEN
         do k4=1,4
         do kkrec=1,5
            if(kkrec.eq.1) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'air.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'air.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'air.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'air.',NYEAR 
               endif
            else if(kkrec.eq.2) then
               if(k4.eq.1) then
                  write(inname,100) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,101) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,102) NYEAR,  'hgt.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,103) NYEAR,  'hgt.',NYEAR 
               endif
            else if(kkrec.eq.3) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'rhum.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'rhum.',NYEAR 
               endif
            else if(kkrec.eq.4) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'uwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'uwnd.',NYEAR 
               endif
            else if(kkrec.eq.5) then
               if(k4.eq.1) then
                  write(inname,200) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.2) then
                  write(inname,201) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.3) then
                  write(inname,202) NYEAR, 'vwnd.',NYEAR 
               else if(k4.eq.4) then
                  write(inname,203) NYEAR, 'vwnd.',NYEAR 
               endif
            endif
    
            pathaddname = '../DATA/'//DATTYP//'/'//inname
            inquire(file=pathaddname,exist=there)
            if(.not.there) then
               print*,pathaddname,' is not available'
               stop
            endif
            status=nf_open(pathaddname,nf_nowrite,inet6(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'scale_factor',xscl(kkrec,k4))
            status=nf_get_att_double(inet6(kkrec,k4),5,
     &                               'add_offset',xoff(kkrec,k4))
            write(*,*) inet6(kkrec,k4),pathaddname
     &                 ,xscl(kkrec,k4),xoff(kkrec,k4)
         enddo
         enddo

      ENDIF
 100  format(I4,'/',A4,I4,'.00.nc')
 101  format(I4,'/',A4,I4,'.06.nc')
 102  format(I4,'/',A4,I4,'.12.nc')
 103  format(I4,'/',A4,I4,'.18.nc')
 200  format(I4,'/',A5,I4,'.00.nc')
 201  format(I4,'/',A5,I4,'.06.nc')
 202  format(I4,'/',A5,I4,'.12.nc')
 203  format(I4,'/',A5,I4,'.18.nc')

      k4=NHOUR/6+1
      it=NDAY
      if(MONTH.EQ.2) it=it+ 31
      if(MONTH.EQ.3) it=it+ 59
      if(MONTH.EQ.4) it=it+ 90
      if(MONTH.EQ.5) it=it+120
      if(MONTH.EQ.6) it=it+151
      if(MONTH.EQ.7) it=it+181
      if(MONTH.EQ.8) it=it+212
      if(MONTH.EQ.9) it=it+243
      if(MONTH.EQ.10)it=it+273
      if(MONTH.EQ.11)it=it+304
      if(MONTH.EQ.12)it=it+334
      if(MOD(NYEAR,4)  .EQ.0.AND.MONTH.GT.2) it=it+1
      if(MOD(NYEAR,100).EQ.0.AND.MONTH.GT.2) it=it-1
      if(MOD(NYEAR,400).EQ.0.AND.MONTH.GT.2) it=it+1
      do k=1,4
         start(k)=1
      enddo
      do k=5,10
         start(k)=0
         count(k)=0
      enddo
      count(1)=ilon
      count(2)=jlat
      count(3)=37
      count(4)=365
      if(MOD(NYEAR,4)  .EQ.0) count(4)=366
      if(MOD(NYEAR,100).EQ.0) count(4)=365
      if(MOD(NYEAR,400).EQ.0) count(4)=366
      start(4)=it
      count(4)=1
!bxq_
      do kkrec=1,5
         inet=inet6(kkrec,k4)
         status=nf_get_vara_int2(inet,5,start,count,work)
         xscale = xscl(kkrec,k4)
         xadd   = xoff(kkrec,k4)
         if(kkrec.eq.1) then
            do j=1,jlat
            do i=1,ilon
!              Tvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Tvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Tvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Tvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Tvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Tvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Tvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Tvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Tvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Tvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Tvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Tvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Tvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Tvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Tvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Tvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Tvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Tvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Tvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Tvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Tvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Tvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Tvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Tvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.2) then
            do j=1,jlat
            do i=1,ilon
!              Hvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Hvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Hvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Hvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Hvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Hvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Hvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Hvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Hvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Hvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Hvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Hvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Hvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Hvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Hvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Hvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Hvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Hvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Hvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Hvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Hvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Hvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Hvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Hvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,klev
            do j=1,jlat
            do i=1,ilon
               Hvar(i,j,k)=Hvar(i,j,k)/9.80616
            enddo
            enddo
            enddo
         else if(kkrec.eq.3) then
            do j=1,jlat
            do i=1,ilon
!              RHvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               RHvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               RHvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               RHvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               RHvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               RHvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               RHvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               RHvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               RHvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               RHvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               RHvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               RHvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               RHvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               RHvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               RHvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               RHvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               RHvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               RHvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               RHvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               RHvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               RHvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               RHvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               RHvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               RHvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
            do k=1,23
            do j=1,jlat
            do i=1,ilon
               RHvar(i,j,k)=amax1(RHvar(i,j,k)*0.01,0.00)
            enddo
            enddo
            enddo
         else if(kkrec.eq.4) then
            do j=1,jlat
            do i=1,ilon
!              Uvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
               Uvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Uvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Uvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Uvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Uvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Uvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Uvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Uvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Uvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Uvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Uvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Uvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Uvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Uvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Uvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Uvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Uvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Uvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Uvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Uvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Uvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Uvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Uvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         else if(kkrec.eq.5) then
            do j=1,jlat
            do i=1,ilon
!              Vvar(i,jlat+1-j,k) =work(i,j,k)*xscale+xadd
               Vvar(i,jlat+1-j,1) =work(i,j,1)*xscale+xadd
               Vvar(i,jlat+1-j,2) =work(i,j,2)*xscale+xadd
               Vvar(i,jlat+1-j,3) =work(i,j,3)*xscale+xadd
               Vvar(i,jlat+1-j,4) =work(i,j,4)*xscale+xadd
               Vvar(i,jlat+1-j,5) =work(i,j,5)*xscale+xadd
               Vvar(i,jlat+1-j,6) =work(i,j,6)*xscale+xadd
               Vvar(i,jlat+1-j,7) =work(i,j,7)*xscale+xadd
               Vvar(i,jlat+1-j,8) =work(i,j,8)*xscale+xadd
               Vvar(i,jlat+1-j,9) =work(i,j,9)*xscale+xadd
               Vvar(i,jlat+1-j,10)=work(i,j,10)*xscale+xadd
               Vvar(i,jlat+1-j,11)=work(i,j,11)*xscale+xadd
               Vvar(i,jlat+1-j,12)=work(i,j,13)*xscale+xadd
               Vvar(i,jlat+1-j,13)=work(i,j,15)*xscale+xadd
               Vvar(i,jlat+1-j,14)=work(i,j,17)*xscale+xadd
               Vvar(i,jlat+1-j,15)=work(i,j,18)*xscale+xadd
               Vvar(i,jlat+1-j,16)=work(i,j,20)*xscale+xadd
               Vvar(i,jlat+1-j,17)=work(i,j,22)*xscale+xadd
               Vvar(i,jlat+1-j,18)=work(i,j,24)*xscale+xadd
               Vvar(i,jlat+1-j,19)=work(i,j,26)*xscale+xadd
               Vvar(i,jlat+1-j,20)=work(i,j,28)*xscale+xadd
               Vvar(i,jlat+1-j,21)=work(i,j,31)*xscale+xadd
               Vvar(i,jlat+1-j,22)=work(i,j,34)*xscale+xadd
               Vvar(i,jlat+1-j,23)=work(i,j,37)*xscale+xadd
            enddo
            enddo
         endif
      enddo
! 
      return
      end
