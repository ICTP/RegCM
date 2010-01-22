c ** PGI compile: pgf77 CMAP2RCM.f -L../../Commons/env/liblinux -lnetcdf
c ** IFC compile: ifc -tpp7 -O3 -cm -w -w90 -w95 CMAP2RCM.f -L../../Commons/env/liblinux -lnetcdf -o cmap
      PROGRAM RWDATA

      implicit none
      include 'domain.param'
      include 'cmap.param'

      integer ny, nx
      PARAMETER (ny=iy-2,nx=jx-2)
 
      real xlon(ny,nx), xlat(ny,nx), xlonx(iy,jx), xlatx(iy,jx)
     &   , GRIDRCM(NY,NX), GRIDNC(nx,ny)
     &   , clatx, clonx, platx, plonx, dsx

      real*4 gridlalo(NLON,NLAT), LONI(NLON), LATI(NLAT)

C ******    NetCDF Output Stuff
      real xlat1d(ny), xlon1d(nx), vmisdat, xmin, xmax
      INTEGER idim(3), idt0(4), idt1(4), idt2(4), ndim, ncid, julnc
     &      , idcdf, idcdfout, ncopn, ncdid, ncvid, ncnowrit
      integer*2 imisdat
      parameter (vmisdat =-32767.,xmin=-10., xmax=100., ndim=3)
      parameter(ncnowrit = 0)
      character lnam*20, units*13, vnamnew*10, lnamnew*20
      REAL vvarmin(3), vvarmax(3), sigf(2), scale, offset, sig1d(1)
      real*8 xhro

c ****** Other stuff
      real aaa, grdfac, ddeg
      character iprojx*6
      integer i, j, idot, nlo, nla, ierr, idatex0, idatex1, idatex2
     &      , idatex, imo, iyr, idy, ihr, idyr, idmo, itim
     &      , iyy, jxx, kzz

      print*,'INPUT FILE:',infile
      print*,'OUTPUT FILE:',outfile

C ******  SETUP LONGITUDES AND LATITUDES RCM GRID
      open (unit=10,file='fort.10',status='old',form='unformatted'
     &    ,recl=iy*jx*ibyte,access='direct')
      read (10,rec=1) iyy,jxx,kzz,dsx,clatx,clonx,platx,plonx
     &              , grdfac,iprojx
      if (iyy.ne.iy .or. jxx.ne.jx) then
        print*,'DOMAIN.INFO is inconsistent with domain.param'
        print*,'  domain.param:    iy=',iy,   '   jx=',jx
        print*,'  DOMAIN.INFO:     iy=',iyy,  '   jx=',jxx
        stop 780
      end if
      if (abs(ds*1000.-dsx).gt.0.01) then
        print*,'DOMAIN.INFO is inconsistent with domain.param'
        print*,'  domain.param: ds=',ds*1000.
        print*,'  DOMAIN.INFO:  ds=',dsx
        stop 781
      end if
      if (clatx.ne.clat .or. clonx.ne.clon .or.
     &    platx.ne.plat .or. plonx.ne.plon) then
        print*,'DOMAIN.INFO is inconsistent with domain.param'
        print*,'  domain.param:  clat=',clat, ' clon=',clon
        print*,'  DOMAIN.INFO:   clat=',clatx,' clon=',clonx
        print*,'  domain.param:  plat=',plat, ' plon=',plon
        print*,'  DOMAIN.INFO:   plat=',platx,' plon=',plonx
        stop 782
      end if
      if (iprojx.ne.iproj) then
        print*,'DOMAIN.INFO is inconsistent with domain.param'
        print*,'  domain.param: iproj=',iproj
        print*,'  DOMAIN.INFO:  iproj=',iprojx
        stop 783
      end if
      read (10,rec=5) ((xlatx(i,j),j=1,jx),i=1,iy)
      read (10,rec=6) ((xlonx(i,j),j=1,jx),i=1,iy)
      close(10)
      idot = 0
      do j=1,nx
      do i=1,ny
        xlat(i,j) = xlatx(i+1,j+1)
        xlon(i,j) = xlonx(i+1,j+1)
      end do
      end do

C ****** SET UP LONGITUDES AND LATITUDES FOR LAT/LON GRID
      DO NLO=1,NLON
        LONI(NLO) = GLON1 + dlon*FLOAT(NLO-1)
      end do
      DO NLA=1,NLAT
        LATI(NLA) = GLAT1 + dlat*FLOAT(NLA-1)
      end do

C ****** SET UP NetCDF STUFF
      print*,'CALL PARAM'
      CALL PARAM(nx,ny,1,ds,ddeg,clat,clon,plat,plon
     &   , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim)
      print*,vvarmin,vvarmax
      sigf(1) = 1.
      sigf(2) = 1.
      sig1d(1) = 1.
      idatex0 = idatecmap0*10000 + 1500
      idatex1 = idatecmap1*10000 + 1500
      idatex2 = idatecmap2*10000 + 1500
      CALL SPLITDATE (idatex0,idt0)
      CALL SPLITDATE (idatex1,idt1)
      CALL SPLITDATE (idatex2,idt2)

      aaa=2.**16.-1.

C ****** READ/WRITE DATA
      CALL rcrecdf (outfile,ncid,vvarmin,vvarmax,3,ierr)
      do ifld=1,nfld
        print*,'OPENING NetCDF FILE: ',outfile(ifld)
        idcdf = NCOPN(infile(ifld),ncnowrit,ierr) ! Open input file
        idatex = idatex1
        imo = idt1(2)
        iyr = idt1(1)
        idyr = idt1(1)-idt0(1)
        if (idyr.lt.0) idyr = 0
        idmo = idt1(2)-idt0(2)
        itim = idyr*12 + idmo + 1
        print*,'ITIM=',itim
        if (ifld.eq.1) then
          lnamnew ='CMAP Standard Precip'
          vnamnew = 'precip1'
        else if (ifld.eq.2) then
          lnamnew ='CMAP Enhanced Precip'
          vnamnew = 'precip2'
        end if
        do while (idatex.le.idatex2)
          CALL JULIAN(idatex,julnc,iyr,imo,idy,ihr)
          xhro = float(julnc)
          PRINT*,'READ/WRITE: ',vnam(ifld),idatex,itim,xhro
          CALL READCDF3D(idcdf,infile(ifld),itim,vnam(ifld)
     &           , gridlalo,scale,offset,vmisdat,lnam,units)

          CALL BILINX(GRIDRCM,GRIDLALO,XLON,XLAT,LONI,LATI
     &       , NY,NX,NLON,NLAT,1)
c         CALL BILINX_OLD(GRIDLALO,LONI,LATI,NLON,NLAT
c    A       , GRIDRCM,XLON,XLAT,NY,NX,1,xmin,vmisdat)

          do j=1,ny
          do i=1,nx
            gridnc(i,j) = gridrcm(j,i)
          end do
          end do
          scale=(xmax-xmin)/aaa
          offset=(xmax+xmin)/2.
          CALL WRITECDFI2(ncid,vnamnew,gridnc,nx,ny,1,idim
     &       , xhro,lnamnew,units,scale,offset
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sig1d,0,xmin)
          imo = imo + 1
          if (imo.gt.12) then
            imo = 1
            iyr = iyr + 1
          end if
          idatex = iyr*1000000 + imo*10000 + 1500
          itim = itim + 1
        end do
      end do
      call clscdf(ncid,ierr)

      END

      SUBROUTINE JULIAN(idate,julnc,iyr,imo,idy,ihr)

      implicit none
      integer lenmon(12), jprev(12), ileap, j, julnc, julday
      integer iyr,imo,idy,ihr,idate,iyrm1

      data lenmon / 31, 28, 31, 30, 31, 30
     &            , 31, 31, 30, 31, 30, 31 /

      iyr = idate/1000000
      imo = idate/10000-iyr*100
      idy = idate/100-iyr*10000-imo*100
      ihr = idate-idate/100*100
      ileap = mod(iyr,4)
      if(ileap.eq.0) then
        lenmon(2) = 29
      else
        lenmon(2) = 28
      end if

      if (ihr.gt.23) then
        ihr = ihr - 24
        idy = idy + 1
      end if
      if (idy.gt.lenmon(imo)) then
        idy = 1
        imo = imo + 1
      end if
      if (imo.gt.12) then
        imo = 1
        iyr = iyr + 1
      end if
      idate = iyr*1000000 + imo*10000 + idy*100 + ihr

      iyrm1 = iyr - 1


      jprev(1) = 0
      do j=2,12
        jprev(j)  = jprev(j-1) + lenmon(j-1)
      end do

      julday = idy + jprev(imo) - 1

      julnc = ((iyr-1900)*365 + julday + int((iyrm1-1900)/4))*24 + ihr

      return
      end

      SUBROUTINE READCDF3D(idcdf,infil,itim,varnam
     &         , vals,scale,offset,vmisdat,lnam,units)

      implicit none
      include 'cmap.param'
      character infil*(*),varnam*(*),dumnam*31,lnam*(*),units*(*)
	
      real*4 scale,offset,error,vmisdat,vgrid2d(nlon,nlat),misdat
      real vals(nlon,nlat)
      integer*4 istart(3),icount(3)
     &	 ,idcdf,iflag,latid,lonid,timid,invarid
     &	 ,latsiz,lonsiz,timsiz
     &      , idcdfout, ncopn, ncdid, ncvid, ncnowrit
      integer itim,i,j,jj

c *** Check to make sure dimensions match ***
      latid = NCDID(idcdf, 'lat', iflag)
      CALL NCDINQ(idcdf, latid, dumnam, latsiz, iflag)
      lonid = NCDID(idcdf, 'lon', iflag)
      CALL NCDINQ(idcdf, lonid, dumnam, lonsiz, iflag)
      timid = NCDID(idcdf, 'time', iflag)
      CALL NCDINQ(idcdf, timid, dumnam, timsiz, iflag)
      if (latsiz.ne.nlat .or. lonsiz.ne.nlon .or. itim.gt.timsiz) then
        print*,'DECLARED DIMENSIONS AND DATA DEMENSIONS DO NOT MATCH'
        print*,'  LAT: latsiz=',latsiz,' nlat=',nlat
        print*,'  LON: lonsiz=',lonsiz,' nlon=',nlon
        print*,' TIME: timsiz=',timsiz,' itim=',itim
        print*,'latid,lonid,timid'
        print*,lonid,latid,timid
        print*,varnam
        stop 999
      end if
      istart(timid) = itim
      istart(lonid) = 1
      istart(latid) = 1
      icount(timid) = 1
      icount(lonid) = nlon
      icount(latid) = nlat

c  /*get variable and attributes*/
      invarid = NCVID(idcdf, varnam, iflag)
      CALL NCVGT(idcdf, invarid, istart, icount, vgrid2d, iflag)
      CALL NCAGT (idcdf, invarid, 'add_offset', offset, iflag)
      CALL NCAGT (idcdf, invarid, 'scale_factor', scale, iflag)
      CALL NCAGT (idcdf, invarid, 'missing_value', misdat, iflag)
      CALL NCAGTC (idcdf, invarid, 'long_name', lnam, 50, iflag)
      CALL NCAGTC (idcdf, invarid, 'units', units, 50, iflag)

      do j=1,nlat
      do i=1,nlon
c       if (vgrid2d(i,j).lt.0.0) print*,vgrid2d(i,j)
        jj = nlat - j + 1
        if (vgrid2d(i,j).gt.misdat) then
          vals(i,jj) = vgrid2d(i,j)
        else
          vals(i,jj) = vmisdat
        end if
      end do
      end do

      return
      end

      SUBROUTINE BILINX(B3,B2,ALON,ALAT,GLON,GLAT,JX,IY,NLON,NLAT,NLEV)
      implicit none
C
C  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
C  RECTANGULAR GRID TO A GRID DESCRIBED BY ALON AND ALAT OF GRID2.
C  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
C  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT.
C  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
C  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
C  AS WEIGHTING FACTORS IN THE INTERPOLATION.
C  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
C  INTERPOLATED BECAUSE ALAT AND ALON ARE NOT DEFINED FOR
C  THE CROSS POINTS IN THE RegCM MODEL.
C
C  B2(JX,IY,NLEV) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
C  B3(JX,IY,NLEV) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
C  GLON......LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
C  GLAT......LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
C  P.........EAST-WEST WEIGHTING FACTOR.
C  Q.........NORTH-SOUTH WEIGHTING FACTOR.
C  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
C  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
C
      INTEGER NLON,NLAT,JX,IY,NLEV
      REAL    B3(JX,IY,NLEV),B2(NLON,NLAT,NLEV)
      REAL    ALON(JX,IY),ALAT(JX,IY),GLON(NLON),GLAT(NLAT)
      INTEGER I,J,I1,II,I2,J1,JJ,J2,K
      REAL    P1,P2,Q1,Q2,DELT
C
      DELT = GLON(2)-GLON(1)
      DO 120 J=1,IY
      DO 110 I=1,JX

      I1 = 1000
      DO II=1,NLON-1
         IF(ALON(I,J).GE.GLON(II).AND.ALON(I,J).LT.GLON(II+1))THEN
            P1=ALON(I,J)-GLON(II)
            P2=GLON(II+1)-ALON(I,J)
            I1=II
            I2=II+1
         ELSE IF(ALON(I,J).GE.GLON(II)-360.AND.
     &           ALON(I,J).LT.GLON(II+1)-360.) THEN
            P1=ALON(I,J)-(GLON(II)-360.)
            P2=(GLON(II+1)-360.)-ALON(I,J)
            I1=II
            I2=II+1
         ELSE IF(ALON(I,J).GE.360.-DELT.AND.ALON(I,J).LT.360.) THEN
            P1=ALON(I,J)-(360.-DELT)
            P2=360.-ALON(I,J)
            I1=NLON
            I2=1
         ELSE IF(ALON(I,J).GE.-DELT.AND.ALON(I,J).LT.0.) THEN
            P1=ALON(I,J)+DELT
            P2=-ALON(I,J)
            I1=NLON
            I2=1
         ENDIF
      ENDDO
      IF(I1.EQ.1000) STOP 'Could not find the right longitute'
      J1 = 1000
      DO JJ=1,NLAT-1
         IF(ALAT(I,J).GE.GLAT(JJ).AND.ALAT(I,J).LT.GLAT(JJ+1))THEN
            Q1=ALAT(I,J)-GLAT(JJ)
            Q2=GLAT(JJ+1)-ALAT(I,J)
            J1=JJ
            J2=JJ+1
         ELSE IF(ALAT(I,J).EQ.GLAT(NLAT)) THEN
            J1=NLAT
         ENDIF
      ENDDO
      IF(J1.EQ.1000) STOP 'Could not find the right latitude'
      IF(J1.GT.0.AND.J1.LT.NLAT) THEN
         DO K=1,NLEV
            B3(I,J,K)=( (B2(I1,J1,K)*P2+B2(I2,J1,K)*P1)*Q2
     &                 +(B2(I1,J2,K)*P2+B2(I2,J2,K)*P1)*Q1 )
     &               /(P1+P2)/(Q1+Q2)
         ENDDO
      ELSE IF(J1.EQ.NLAT) THEN
         DO K=1,NLEV
            B3(I,J,K)= B2(2,NLAT,K)
         ENDDO
      ELSE
         STOP 'STOP IN BILINX'
      ENDIF
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BILINX_OLD(  IN, LONI, LATI, NLONI, NLATI
     A                 , OUT, LONO, LATO,    IY,    JX, NFLDS
     &                 , xming, vmisdat)
CCC   SUBROUTINE BILINX (B3,B2,XLONS,XLATS,G4LON,G4LAT,NLON2,NLAT2
CCC  A,                 NLON4,NLAT4,NVECTS,NIN,NJN,INDX)
C
C  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
C  RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF GRID2.
C  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
C  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT.
C  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
C  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
C  AS WEIGHTING FACTORS IN THE INTERPOLATION.
C  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
C  INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
C  THE CROSS POINTS IN THE RCM MODEL.
C
C  IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
C  OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
C  LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
C  LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
C  P.........EAST-WEST WEIGHTING FACTOR.
C  Q.........NORTH-SOUTH WEIGHTING FACTOR.
C  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
C  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
 
      REAL  IN(NLONI,NLATI,NFLDS), LONI(NLONI), LATI(NLATI)
      REAL OUT(IY,JX,NFLDS), LONO(IY,JX), LATO(IY,JX)
      REAL LON360
 
      DO 120 J=1,JX
      DO 110 I=1,IY
 
      YIND = (((LATO(I,J)-LATI(1))/(LATI(NLATI)-LATI(1)))
     A     * FLOAT(NLATI-1))+1.
      JQ=INT(YIND)
      JQP1=MIN0(JQ+1,NLATI)
      Q=YIND-JQ
 
      LON360 = LONO(I,J)
cjsp  IF(LONO(I,J).LT.0.) LON360 = LONO(I,J) + 360.
      XIND = (((LON360   -LONI(1))/(LONI(NLONI)-LONI(1)))
     A     * FLOAT(NLONI-1))+1.
      IP=INT(XIND)
      IPP1=MIN0(IP+1,NLONI)
      P=XIND-IP
 
        DO 100 L=1,NFLDS

          if (in(ip,jq,l).le.xming.and.in(ipp1,jq,l).gt.xming) then
            temp1 = in(ipp1,jq,l) 
          elseif(in(ip,jq,l).le.xming.and.in(ipp1,jq,l).le.xming) then
            temp1 = vmisdat
          elseif(in(ip,jq,l).gt.xming.and.in(ipp1,jq,l).le.xming) then
            temp1 = in(ip,jq,l)
          else
            temp1 = (1.0-p)*in(ip,jq,l) + p*in(ipp1,jq,l)
          endif
    
          if (in(ipp1,jqp1,l).le.xming.and.in(ip,jqp1,l).gt.xming) then
            temp2 = in(ip,jqp1,l) 
          elseif(in(ipp1,jqp1,l).le.xming
     &           .and.in(ip,jqp1,l).le.xming) then
            temp2 = vmisdat
          elseif(in(ipp1,jqp1,l).gt.xming
     &           .and.in(ip,jqp1,l).le.xming) then
            temp2 = in(ipp1,jqp1,l)
          else
            temp2 = p*in(ipp1,jqp1,l) + (1.0-p)*in(ip,jqp1,l)
          endif
  
          if(temp1.eq.xming.and.temp2.ne.xming) then
            out(i,j,l) = temp2
          elseif(temp1.eq.xming.and.temp2.eq.xming) then
            out(i,j,l) = xming
          elseif(temp1.ne.xming.and.temp2.eq.xming) then
            out(i,j,l) = temp1
          else
            out(i,j,l) = (1. - q)*temp1 + q*temp2
          endif

  100   CONTINUE
  110 CONTINUE
 
  120 CONTINUE
 
      RETURN
      END

      SUBROUTINE MAPLAM( XLON, XLAT, XMAP, CORIOL, IY, JX
     A                 , CLON, CLAT, DELX, IDOT )
 
      REAL XLON(IY,JX), XLAT(IY,JX), XMAP(IY,JX), CORIOL(IY,JX)
      REAL CLON, CLAT, DELX
 
C  COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
C  LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND 60.N.
 
C  IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
C  CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
C  DELX IS GRID SPACING IN METERS.
C  IDOT=1 FOR DOT GRID, 0 FOR CROSS GRID.
C  CORIOLIS FACTOR DEFINED ON DOT GRID ONLY.

      XOMEGA=7.2722E-5              ! ANG. ROT OF EARTH IN S**-1
      GRDFAC=.716                   ! FACTOR FOR LAMBERT CONFORMAL
      EARTHR=6.37E6                 ! METERS
      PSI1=30.                      ! LAT AT WHICH PROJ. IS TRUE
C                                   !   TRUE AT COLAT OF PSI1 TOO
      PIFAC=ATAN(1.)/45.            ! CONVERT DEGREES TO RADIANS
      PSI0=90.-CLAT                 ! COLAT OF CENTER
C
      PSI1=PIFAC*PSI1               ! CONVERT TO RADIANS
      PSI0=PIFAC*PSI0               ! CONVERT TO RADIANS
      RGRDF=1./GRDFAC               ! RECIPROCAL OF GRID FACTOR
      ANS=EARTHR*SIN(PSI1)*RGRDF    ! A*SIN(PS1)/N
      C1=-CLON-90.*RGRDF
      C2=ANS*(TAN(PSI0*0.5)/TAN(PSI1*0.5))**GRDFAC
      PF=TAN(PSI1*0.5)/(ANS**RGRDF)
      XMF=SIN(PSI1)/(TAN(0.5*PSI1)**GRDFAC)
      PI90=2.*ATAN(1.)              ! PI OVER TWO (90 DEGREES)
C
      YC=0.5*FLOAT(IY+IDOT)
      XC=0.5*FLOAT(JX+IDOT)
      DO 10 J=1,JX
      X=(J-XC)*DELX
      DO 10 I=1,IY
      Y=(I-YC)*DELX
      IF (X.NE.0.) THEN
         XLP=ATAN2(Y-C2,X)
         R=ABS(X/COS(XLP))
      ELSE
         XLP=-PI90
         R=ABS(Y-C2)
      ENDIF
      PSI=ABS(2.*ATAN(PF*(R**RGRDF)))        ! CO LATITUDE (RADIANS)
      PFI=PI90-PSI
      XLON(I,J)=XLP*RGRDF/PIFAC-C1           ! LONGITUDE DEGREES EAST
      XLAT(I,J)=PFI/PIFAC                    ! LATITUDE DEGREES NORTH
      XMAP(I,J)=(XMF/SIN(PSI))*(TAN(0.5*PSI)**GRDFAC)
   10 CONTINUE
C
      IF (IDOT.EQ.1) THEN
         XOMEG2=2.*XOMEGA
         DO 20 I=1,IY
         DO 20 J=1,JX
         CORIOL(I,J)=XOMEG2*SIN(XLAT(I,J)*PIFAC)
   20    CONTINUE
      ENDIF
 
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MAPMER( XLON, XLAT, XMAP, CORIOL, IY, JX 
     A                 , CLON,CLAT, POLLON, POLLAT, DS, IDOT)    
      implicit none
C---------------------------------------------------------------------  
C  COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
C  ROTATED POLE MAP CENTERED AT CLON,CLAT. ORIGIN OF ROTATED GRID IS
C  GIVEN BY POLLON AND POLLAT.
C  IMX,JMX,KSV,KSH AND LSIG2 MUST CORRESPOND TO IY,JMAX,KSTRPV,KSTRPH 
C  AND NVERT2 IN THE MASTER INPUT FILE; OTHERWISE THE PROGRAM WILL      
C  ABORT:LMX MUST BE THE MAXIMUM NUMBER OF LEVELS (PRESSURE--OR--SIGMA) 
C  IMAXN AND IMXC ARE NESTED AND NON-EXPANDED GRID DIMENSIONS. IMXC IS  
C  EQUAL TO IMX IF YOU ARE NOT USING THE EXPANDED GRID. SAME FOR J.     

      INTEGER IY,JX,IDOT
      REAL*4  XLAT(IY,JX),XLON(IY,JX), CORIOL(IY,JX),XMAP(IY,JX)
      REAL*4  CLON,CLAT,POLLON,POLLAT,DS
C
      REAL*4  XOMEGA,PIFAC,PFINV,A,CNTRJ,CNTRI,DDEG,XOFF,YOFF
      REAL*4  XR,YR,X,Y,FAI,XOMEGA2
      INTEGER I,J
      LOGICAL SOUTHH
      REAL*4  ZLAT,ZPLAT,AAA
C
      IF(CLAT.LT.0.) THEN
         SOUTHH=.true.
         ZLAT  =-CLAT
         ZPLAT =-POLLAT
      ELSE
         SOUTHH=.false.
         ZLAT  = CLAT
         ZPLAT = POLLAT
      ENDIF
C
      XOMEGA = 7.2722E-5                       ! ANG. ROT OF EARTH IN S**-1
      PIFAC = ATAN(1.)/45.                     ! CONVERT DEGREES TO RADIANS
      PFINV = 1./PIFAC                         ! CONVERT RADIANS TO DEGREES
      A = 6371229.                             ! RADIUS OF EARTH IN METERS
C-----CENTER OF GRID
      CNTRJ = (JX+IDOT)/2.
      CNTRI = (IY+IDOT)/2.

      DDEG=DS*PFINV/A                          ! GRID SPACING IN DEGREES
      XOFF=CLON-POLLON
      YOFF=ZLAT-ZPLAT
C-----CALCULATE X AND Y POSITIONS OF GRID                               
      DO 41 I=1,IY
      DO 41 J=1,JX                                                    
        XR = XOFF + (J-CNTRJ)*DDEG                                            
        YR = YOFF + (I-CNTRI)*DDEG
        YR = 2*PFINV*ATAN(EXP(PIFAC*YR))-90.    
C-----NOW CALCULATE LAT AND LON OF THIS POINT                           
C-----  ROTATE COORDINATES BACK TO NONRATED POLE
        CALL ROT2NROT(XR,YR,POLLON,ZPLAT,X,Y)
        XLON(I,J) = X
        XLAT(I,J) = Y
        FAI = PIFAC*YR
        XMAP(I,J) = 1.0/COS(FAI)
   41 CONTINUE

      IF(IDOT.EQ.1) THEN
         XOMEGA2=2.*XOMEGA
         DO 45 I=1,IY
         DO 45 J=1,JX
            CORIOL(I,J)=XOMEGA2*SIN(XLAT(I,J)*PIFAC)
   45    CONTINUE
      END IF
      IF(SOUTHH) THEN
         DO J=1,JX
         DO I=1,IY/2
            AAA = XLON(IY+1-I,J)
            XLON(IY+1-I,J) = XLON(I,J)
            XLON(I,J) = AAA

            AAA =-XLAT(IY+1-I,J)
            XLAT(IY+1-I,J) =-XLAT(I,J)
            XLAT(I,J) = AAA

            AAA = XMAP(IY+1-I,J)
            XMAP(IY+1-I,J) = XMAP(I,J)
            XMAP(I,J) = AAA

            IF(IDOT.EQ.1) THEN
            AAA =-CORIOL(IY+1-I,J)
            CORIOL(IY+1-I,J) =-CORIOL(I,J)
            CORIOL(I,J) = AAA
            ENDIF
         ENDDO
         ENDDO
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE rot2nrot(LAMS,PHIS,POLLON,POLLAT,LAM,PHI)
      implicit none
c----------------------------------------------------------------------------
c Purpose:
c     Adaption of the DWD-Functions to convert rotated pole coordinates
c     (PHIS,LAMS) into geogrphic coordinates (PHI,LAM). The location of 
c     the rotated pole is passed trough POLLON and POLLAT. POLLON and
c     POLLAT give the origin of the rotated grid. The 
c     first four arguments are input, the last two are output. All angles 
c     are in degrees (north>0, east>0)
c History:
c     05/90   D.MAJEWSKI (DWD)
c     03/93   D.BRESCH (ETHZ)
c     09/96   D.LUETHI (ETHZ)
 
c declaration of arguments:
      REAL*4  LAMS,PHIS,PPHI,PLAM,PHI,LAM,POLLON,POLLAT
 
c declaration of internal vars:
      real*4  zarg1,zarg2
      real*4  zphis,zlams,arg
 
      real*4  ZRPI18,ZPIR18
      DATA    ZRPI18,ZPIR18  / 57.2957795 , 0.0174532925 /
      real*4  ZSINPOL,ZCOSPOL,ZLAMPOL
 
      PLAM=POLLON + 180.
cjsp  IF (POLLAT.GT.0.) THEN
      IF (POLLAT.GE.0.) THEN
        PPHI = 90. - POLLAT
      ELSE
        PPHI = 90. + POLLAT
        PLAM = PLAM - 180.
      ENDIF

      IF (PLAM.GT.180.) PLAM = PLAM-360.
      ZSINPOL = SIN(ZPIR18*PPHI)
      ZCOSPOL = COS(ZPIR18*PPHI)
      ZLAMPOL =     ZPIR18*PLAM

c do case without rotated pole
c      if (abs(pphi-90.).lt.1.e-3) then
c        phi=phis
c        lam=lams
c        return
c      endif
 
C first, the conversion of PHIS to PHI:
      ZPHIS  = ZPIR18*PHIS
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = ZPIR18*ZLAMS
      ARG     = ZCOSPOL*COS(ZPHIS)*COS(ZLAMS) + ZSINPOL*SIN(ZPHIS)
      PHI = ZRPI18*ASIN(ARG)
 
c follows conversion of LAMS to LAM:
      ZPHIS   = ZPIR18*PHIS
      ZLAMS   = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS   = ZPIR18*ZLAMS
      ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     &                          ZCOSPOL*           SIN(ZPHIS)) -
     &          COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     &                          ZCOSPOL*           SIN(ZPHIS)) +
     &          SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LAM =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
          LAM =  90.0
        ELSE
          LAM = -90.0
        ENDIF
      ELSE
        LAM = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF

      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MAPPOL( XLON, XLAT, XMAP, CORIOL, IY, JX
     A                 , CLON, CLAT, DELX, idot )
      implicit none
      INTEGER IY,JX
      REAL*4  CLON,CLAT,DELX

      REAL*4  XLON(IY,JX), XLAT(IY,JX)
      REAL*4  XMAP(IY,JX), CORIOL(IY,JX)
      REAL*4  XN,PSI1,PI
C
C         XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
C         PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
C         PI IS PI.
C
C---------------------------------------------------------------------
      REAL*4  AA,POLE,DEGRAN
      REAL*4  CNTRJ,CNTRI,PSX,CELL,CELL2,R,XCNTR,YCNTR
      REAL*4  X,Y,FLP,FLPP,CELL1,PSIX,xomega,xomeg2,pifac
      INTEGER II1,JJ1,I,J,idot

C  COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
C  LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND 60.N.

C  IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
C  CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
C  DELX IS GRID SPACING IN METERS.
C  ALWAYS FOR CROSS GRID.

C
      AA=6.371229E6
      XN=1.0
      PI=ATAN(1.)*4.
      DEGRAN=PI/180.
C
      POLE=90.
      PSI1=30.
      PSI1=PSI1*DEGRAN
      IF (CLAT.LT.0.0) THEN
         POLE = -90.0
         PSI1 = -30.
         PSI1=PSI1*DEGRAN
      ENDIF
      CNTRJ = FLOAT(JX+idot)/2.
      CNTRI = FLOAT(IY+idot)/2.
C
      PSX = (POLE-CLAT)*DEGRAN
      CELL=AA*SIN(PSX)/XN
      CELL2 = (1. + COS(PSI1))/(1. + COS(PSX))
      R=CELL*(CELL2)**XN
      XCNTR=0.
      YCNTR=-R
C
      II1 = IY
      JJ1 = JX
      DO 70 I=1,II1
         Y=YCNTR+(I-CNTRI)*DELX
         DO 70 J=1,JJ1
            X=XCNTR+(J-CNTRJ)*DELX
            R = SQRT(X*X + Y*Y)
            IF(Y.EQ.0.) THEN
               IF(X.GE.0.) THEN
                  FLP=90.*DEGRAN
               ELSE
                  FLP=-90.*DEGRAN
               ENDIF
            ELSE
               IF (CLAT.LT.0.0) THEN
                  FLP = ATAN2(X,Y)
               ELSE
                  FLP = ATAN2(X,-Y)
               END IF
            END IF
            FLPP = FLP/XN/DEGRAN+CLON
            IF (FLPP.GT.180.0) FLPP = FLPP-360.0
            IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
            XLON(I,J) = FLPP
            IF (CLAT.LT.0.0) R = -R
            CELL = R/AA
            CELL1 = CELL/(1.0 + COS(PSI1))
            CELL2=ATAN(CELL1)
            PSX=2.*CELL2/DEGRAN
            XLAT(I,J) = POLE - PSX
            PSIX=PSX*DEGRAN
            XMAP(I,J) = ((1.0 + COS(PSI1))/(1.0 + COS(PSIX)))**XN
 70   CONTINUE

      IF (IDOT.EQ.1) THEN
         XOMEGA = 7.2722E-5           ! ANG. ROT OF EARTH IN S**-1
         PIFAC = ATAN(1.)/45.         ! CONVERT DEGREES TO RADIANS
         XOMEG2=2.*XOMEGA
         DO I=1,IY
         DO J=1,JX
           CORIOL(I,J)=XOMEG2*SIN(XLAT(I,J)*PIFAC)
         end do
         end do
      ENDIF
      RETURN
      END



      SUBROUTINE MAPROT( XLON, XLAT, XMAP, CORIOL, IMAX, JMAX 
     A                 , XLONC, PHIC, POLLON, POLLAT, DS, IDOT)    
C---------------------------------------------------------------------  
C  COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
C  ROTATED POLE MAP CENTERED AT CLON,CLAT. ORIGIN OF ROTATED GRID IS
C  GIVEN BY POLLON AND POLLAT.
C  IMX,JMX,KSV,KSH AND LSIG2 MUST CORRESPOND TO IMAX,JMAX,KSTRPV,KSTRPH 
C  AND NVERT2 IN THE MASTER INPUT FILE; OTHERWISE THE PROGRAM WILL      
C  ABORT:LMX MUST BE THE MAXIMUM NUMBER OF LEVELS (PRESSURE--OR--SIGMA) 
C  IMAXN AND IMXC ARE NESTED AND NON-EXPANDED GRID DIMENSIONS. IMXC IS  
C  EQUAL TO IMX IF YOU ARE NOT USING THE EXPANDED GRID. SAME FOR J.     

      real XLAT(IMAX,JMAX),XLON(IMAX,JMAX),CORIOL(IMAX,JMAX)
     &    ,XMAP(IMAX,JMAX), XLONC, PHIC, POLLON, POLLAT, DS
C                                                                       
      XOMEGA = 7.2722E-5                       ! ANG. ROT OF EARTH IN S**-1
      PIFAC = ATAN(1.)/45.                     ! CONVERT DEGREES TO RADIANS
      PFINV = 1./PIFAC                         ! CONVERT RADIANS TO DEGREES
      A = 6370000.                             ! RADIUS OF EARTH IN METERS
C-----CENTER OF GRID                                                    
      CNTRJ = JMAX/2.                                                    
      CNTRI = IMAX/2.                                                    
      PRINT 10,PHIC,XLONC                                               
   10 FORMAT(1X,'LATITUDE AND LONGITUDE OF CENTER OF GRID = ',F8.3,'N', 
     1 F8.3,'W')                                                        

      DDEG=DS*PFINV/A                    ! GRID SPACING IN DEGREES
      XOFF=XLONC-POLLON
      YOFF=PHIC-POLLAT
C-----CALCULATE X AND Y POSITIONS OF GRID                               
      DO 41 I=1,IMAX                                                    
      DO 41 J=1,JMAX                                                    
        XR = XOFF + (J-CNTRJ)*DDEG                                            
        YR = YOFF + (I-CNTRI)*DDEG                                            
C-----NOW CALCULATE LAT AND LON OF THIS POINT                           
C-----  ROTATE COORDINATES BACK TO NONRATED POLE
        CALL ROT2NROT(XR,YR,POLLON,POLLAT,X,Y)
        XLON(I,J) = X
        XLAT(I,J) = Y
        FAI = PIFAC*YR
        XMAP(I,J) = 1.0/COS(FAI)
   41 CONTINUE                                                          

      IF(IDOT.EQ.1) THEN
         XOMEGA2=2.*XOMEGA
         DO 45 I=1,IMAX
         DO 45 J=1,JMAX
            CORIOL(I,J)=XOMEGA2*SIN(XLAT(I,J)*PIFAC)
   45    CONTINUE
      END IF

      RETURN                                                            
      END

      SUBROUTINE SPLITDATE (idate,idt)

      integer idt(4)

      idt(1) = idate/1000000
      idt(2) = idate/10000-idt(1)*100
      idt(3) = idate/100-idt(1)*10000-idt(2)*100
      idt(4) = idate-idt(1)*1000000-idt(2)*10000-idt(3)*100

      return
      end


      SUBROUTINE PARAM(nx,ny,kz,ds,ddeg,clat,clon,xplat,xplon
     &         , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim)

      real   vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real xlat(ny,nx), xlon(ny,nx)
      integer idim(ndim)

      ddeg=ds*180./3.14159265359/6370000.
      print*,'DS=',ds,'DDEG=',ddeg
      vvarmin(1) = xlon(ny/2,1)
      vvarmin(2) = xlat(1,nx/2)
      vvarmin(3) = 1050.
      vvarmax(1) = xlon(ny/2,nx)
      vvarmax(2) = xlat(ny,nx/2)
      vvarmax(3) = 1050.
      idim(1) = nx
      idim(2) = ny
      idim(3) = kz
      do i=1,nx
        xlon1d(i) = xlon(ny/2,i)
      end do
      do j=1,ny
        xlat1d(j) = xlat(j,nx/2)
      end do

      return
      end


      subroutine clscdfCOARDS (cdfid, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine closes an open netCDF file.
c     Aguments :
c        cdfid  int  input   the id of the file to be closed.
c        error  int  output  indicates possible errors found in this
c                            routine.
c                            error = 0   no errors detected.
c                            error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer      cdfid, error

c     Local variable declarations.
      integer      ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     Make sure netCDF errors do not abort execution.
      call ncpopt (NCVERBOS)

c     Close requested file.
      call ncclos (cdfid, error)

c     Reset error options.
      call ncpopt (ncopts)

      end


      subroutine crecdfCOARDS (filnam, cdfid, phymin, phymax, ndim, cfn, 
     &                   error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to create a netCDF file for use with 
c        the UWGAP plotting package.
c           Any netCDF file written to must be closed with the call
c        'call clscdf(cdfid,error)', where cdfid and error are
c        as in the argumentlist below. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  input   the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  input   the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   input   the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        cfn     char  input   constants file name 
c                              ('0' = no constants file).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created crecdf.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, error
      character *(*) filnam,cfn
      real           phymin(*), phymax(*)

c     Local variable declarations.
      character *(20) attnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, k, ibeg, iend, lenfil, ncopts
      data            chrid/'x','y','z','a'/

c     External function declarations.
      integer        strbeg, strend

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     create the netCDF file
      cdfid = nccre (filnam(1:strend(filnam)), NCCLOB, error)
      if (error.ne.0) go to 920

c     define global attributes
      do k=1,ndim
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymin(k),error)
        if (error.gt.0) goto 920

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymax(k),error)
        if (error.gt.0) goto 920
      enddo

c     define constants file name
      if (cfn.ne.'0') then
        ibeg = strbeg(cfn)
        iend = strend(cfn)
        lenfil = iend - ibeg + 1
        call ncaptc (cdfid, NCGLOBAL, 'constants_file_name', 
     &             NCCHAR, lenfil+1, cfn // char(0) , error)
        if (error.gt.0) goto 920
      endif

c     End variable definitions.
      call ncendf (cdfid, error)
      if (error.gt.0) goto 920

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'create the data file in subroutine crecdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end



      function strbeg (string)
c-----------------------------------------------------------------------
c     Purpose:
c        This function returns the position of the first nonblank
c           character in the input string.  Returns 0 if the entire
c           string is blank.
c     Arguments:
c        string  char  input  string to be examined.
c     History:
c-----------------------------------------------------------------------

c     Function declaration
      integer strbeg

c     Argument declarations
      character*(*) string

c     Local variable declarations.
      integer i

c     External function declarations.
      integer len

      do 10 i = 1, len(string)
         strbeg = i
         if (string(i:i).ne.' ' .and. string(i:i).ne.char(0)) then
            return
         endif
 10   continue
      strbeg = 0
      end


      function strend (string)
c-----------------------------------------------------------------------
c     Purpose:
c        This function returns the position of the last nonblank
c           character in the input string.  Returns 0 if the entire
c           string is blank.
c     Arguments:
c        string  char  input  string to be examined.
c     History:
c-----------------------------------------------------------------------

c     Function declaration
      integer strend

c     Argument declarations
      character*(*) string

c     Local variable declarations.
      integer i

c     External function declarations.
      integer len

      do 10 i = len(string), 1, -1
         strend = i
         if (string(i:i) .ne. ' ' .and. string(i:i) .ne. char(0)) then
            return
         endif
 10   continue
      strend = 0
      end


      subroutine getdefCOARDS (cdfid, varnam, ndim, misdat, 
     &                              vardim, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions and attributes of 
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              Is dimensioned at least (ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),vnam
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
cjsp  if (vartyp.ne.NCFLOAT) error=1
      if (vartyp.ne.NCSHORT) error=1
      if (error.gt.0) goto 920

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_value',misdat,error)
      if (error.gt.0) goto 920     

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) '*ERROR*: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end

      subroutine getdefCOARDSFLT (cdfid, varnam, ndim, misdat, 
     &                              vardim, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions and attributes of 
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              Is dimensioned at least (ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(ndim), error, cdfid
      real           misdat

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),vnam
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) then
	print*,'cdfid,ndims,nvars,ngatts,recdim,error'
	print*,cdfid,ndims,nvars,ngatts,recdim,error
        goto 920
      end if

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) then
          print*,'cdfid,i,dimnam(i),dimsiz(i),error'
          print*,cdfid,i,dimnam(i),dimsiz(i),error
          goto 920
	end if
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCFLOAT) error=1
cjsp  if (vartyp.ne.NCSHORT) error=1
      if (error.gt.0) then
        print*,'cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error'
        print*,cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error
        goto 920
      end if

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_value',misdat,error)
      if (error.gt.0) then
	print*,'cdfid,id,misdat,error'
	print*,cdfid,id,misdat,error
        goto 920     
      end if

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) '*ERROR*: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end


      subroutine putdatCOARDS(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to write the data of a variable
c        to an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation), the 
c        variable must be defined with a call to putdef.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be written to the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). If 'time' is not
c                              yet known to the file, a knew time-level is
c                              allocated and appended to the times-array.
c        level   int input     the horizontal level(s) to be written 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call writes the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call writes the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int output    indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is new, but
c                                          appending it would yield a non
c                                          ascending times-array.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

C     Declaration of local variables

      character*(*) varnam
      character*(20) chars
      integer cdfid


      integer*2 	dat(*)
      real	misdat
      real*8    time, timeval
      real*8    dvrange(2)

      integer	corner(4),edgeln(4),did(4),vardim(4),ndims
      integer	error, ierr
      integer	level,ntime
      integer	idtime,idvar,iflag
      integer	i,ik

      integer   strend
 
      call ncpopt(NCVERBOS)

c     get definitions of data
      call getdefCOARDS (cdfid, varnam(1:strend(varnam)), ndims, misdat, 
     &                           vardim, ierr)
      if (ierr.ne.0)  print *,'*ERROR* in getdef in putdat'

c     get id of variable
      idvar=ncvid(cdfid,varnam(1:strend(varnam)),ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'

c     get times-array 
      did(4)=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* did(4) in putdat'
      call ncdinq(cdfid,did(4),chars,ntime,ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncdinq in putdat'
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'
C     Check if a new time step is starting
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in putdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then		! new time step
        ntime=ntime+1
        iflag=ntime
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvid in putdat'
        call ncvpt1(cdfid,idtime,ntime,time,ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvpt1 in putdat'
        call ncagt (cdfid,idtime,'actual_range',dvrange,ierr)
        if ((dvrange(1).gt.time).or.(dvrange(1).eq.0.))
     &         dvrange(1) = time 
        if ((dvrange(2).lt.time).or.(dvrange(2).eq.0.))
     &         dvrange(2) = time
        call ncapt (cdfid,idtime,'actual_range',NCDOUBLE,
     &                2,dvrange,ierr)
      endif

C     Define data volume to write on the NetCDF file in index space
      corner(1)=1               ! starting corner of data volume
      corner(2)=1
      edgeln(1)=vardim(1)       ! edge lengthes of data volume
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        ik = 3
      else
        ik = 4
        corner(3)=level
        edgeln(3)=1
      endif
      corner(ik)=iflag
      edgeln(ik)=1
      
C     Put data on NetCDF file

c      print *,'vor Aufruf ncvpt d.h. Daten schreiben in putdat '
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln
     
c      print *,'dat(1)=', dat(1)

      call ncvpt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvpt in putdat - Put data on NetCDF file'
      endif

C     Synchronize output to disk and close the files

      call ncsnc(cdfid,ierr)
      if (ierr.ne.0) print *, '*ERROR* in ncsnc in putdat'
      end

      subroutine putdatCOARDSFLT(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to write the data of a variable
c        to an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation), the 
c        variable must be defined with a call to putdef.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be written to the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). If 'time' is not
c                              yet known to the file, a knew time-level is
c                              allocated and appended to the times-array.
c        level   int input     the horizontal level(s) to be written 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call writes the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call writes the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int output    indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is new, but
c                                          appending it would yield a non
c                                          ascending times-array.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

C     Declaration of local variables

      character*(*) varnam
      character*(20) chars
      integer cdfid


      real      	dat(*)
cjsp  integer*2 	dat(*)
      real	misdat
      real*8    time, timeval
      real*8    dvrange(2)

      integer	corner(4),edgeln(4),did(4),vardim(4),ndims
      integer	error, ierr
      integer	level,ntime
      integer	idtime,idvar,iflag
      integer	i,ik

      integer   strend
 
      call ncpopt(NCVERBOS)

c     get definitions of data
      call getdefCOARDSFLT(cdfid, varnam(1:strend(varnam)), ndims
     &   , misdat, vardim, ierr)
      if (ierr.ne.0)  print *,'*ERROR* in getdef in putdat'

c     get id of variable
      idvar=ncvid(cdfid,varnam(1:strend(varnam)),ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'

c     get times-array 
      did(4)=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* did(4) in putdat'
      call ncdinq(cdfid,did(4),chars,ntime,ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncdinq in putdat'
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'
C     Check if a new time step is starting
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) then
	  print *,'*ERROR* in ncvgt1 in putdat'
	  print*,ntime,cdfid,idtime,i,timeval,ierr
	end if
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then		! new time step
        ntime=ntime+1
        iflag=ntime
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvid in putdat'
        call ncvpt1(cdfid,idtime,ntime,time,ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvpt1 in putdat'
        call ncagt (cdfid,idtime,'actual_range',dvrange,ierr)
        if ((dvrange(1).gt.time).or.(dvrange(1).eq.0.))
     &         dvrange(1) = time 
        if ((dvrange(2).lt.time).or.(dvrange(2).eq.0.))
     &         dvrange(2) = time
        call ncapt (cdfid,idtime,'actual_range',NCDOUBLE,
     &                2,dvrange,ierr)
      endif

C     Define data volume to write on the NetCDF file in index space
      corner(1)=1               ! starting corner of data volume
      corner(2)=1
      edgeln(1)=vardim(1)       ! edge lengthes of data volume
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        ik = 3
      else
        ik = 4
        corner(3)=level
        edgeln(3)=1
      endif
      corner(ik)=iflag
      edgeln(ik)=1
      
C     Put data on NetCDF file

c      print *,'vor Aufruf ncvpt d.h. Daten schreiben in putdat '
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln
     
c      print *,'dat(1)=', dat(1)

      call ncvpt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvpt in putdat - Put data on NetCDF file'
      endif

C     Synchronize output to disk and close the files

      call ncsnc(cdfid,ierr)
      if (ierr.ne.0) print *, '*ERROR* in ncsnc in putdat'
      end


      subroutine putdefCOARDS(cdfid,varnam,ndim,misdat,clname,clunits
     &           , offset, scale, vardim, varmin, varmax, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to define the dimensions and the
c        attributes of a variable on an IVE-NetCDF file for use with the
c        IVE plotting package. Prior to calling this routine, the file must
c        be opened with a call to opncdf (extend an existing file) or
c        crecdf (create a new file).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c        ndim    int   input   the number of dimensions (ndim<=4). 
c                              Upon ndim=4, the fourth dimension of the
c                              variable is specified as 'unlimited'
c                              on the file (time-dimension). It can 
c                              later be extended to arbitrary length.
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              Is dimensioned at least Min(ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error =10   other errors detected.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam,clname,clunits
      integer        vardim(*), ndim, error, cdfid
      real           varmin(*), varmax(*)
      real           scale,offset
cjsp  real           misdat
      integer*2      misdat

c     Local variable declarations.
      character *(20) dimnam,dimchk
      character *(20) dimnams(MAXNCDIM)
      character *(5)  rdim(MAXDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),idtime,i,k,ik,ierr
      integer         ncopts

      integer         ibeg,iend
      integer         idcoor
      real            vrange(2)
      real*8          dvrange(2)
      character*(20)  long_name(MAXDIM)
      character*(30)  units(MAXDIM)
      data            rdim /'lon  ','lat  ','level','time '/
      data            long_name /'Longitude', 'Latitude', 
     &                           'Height_Index', 'Time' /

      data            units / 'degrees_east', 'degrees_north',
     &                  'level' , 'hours since 1900-1-1 00:00:0.0' /


c     External function declarations.
      integer         strbeg, strend

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 10
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c         print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      print*,cdfid
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define spatial dimensions
      ik = 0 
      do k=1,min0(3,ndim)
       if (vardim(k).gt.1) then
        ik = ik + 1
c       define the default dimension-name
        ibeg=strbeg(rdim(k))
        iend=strend(rdim(k))
        dimnam = '                    '
        dimnam(1:1+iend-ibeg)=rdim(k)(ibeg:iend)
        did(ik)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a new dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((dimnam(1:3).eq.dimchk(1:3))) then 
              did(ik)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if ((did(ik).lt.0)) then
c         define the dimension and an array with name of dimension
          did(ik)=ncddef(cdfid,dimnam,vardim(k),error)
          idcoor = ncvdef(cdfid,dimnam,NCFLOAT,1,did(ik),ierr)
          if ((error.ne.0).or.(ierr.ne.0)) goto 921
          vrange(1) = varmin(k)
          vrange(2) = varmax(k)
          iend=strend(long_name(k))
          call ncaptc(cdfid,idcoor,'long_name',
     &                 NCCHAR,iend,long_name(k)(1:iend),ierr)
          iend=strend(units(k))
          call ncaptc(cdfid,idcoor,'units',
     &                 NCCHAR,iend,units(k)(1:iend),ierr)
          call ncapt (cdfid,idcoor,'actual_range',
     &                 NCFLOAT,2,vrange,ierr)
        endif
       endif
      enddo

c     define the times-array
      if (ndim.eq.4) then
        ik = ik +1
c       define dimension 'time'
        did(ik)=ncdid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         this dimension must first be defined
          did(ik) = ncddef (cdfid,'time',NCUNLIM,ierr)
        endif
c       define array 'time'
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         define the times-array
          idtime = ncvdef (cdfid,'time',NCDOUBLE,1,did(ik),ierr)
          call ncaptc(cdfid,idtime,'long_name',NCCHAR,4,'Time',ierr)
          iend = strend(units(4))
          call ncaptc(cdfid,idtime,'units',
     &                 NCCHAR,iend,units(4)(1:iend),ierr)
          dvrange(1)=0.
          dvrange(2)=0.
          call ncapt (cdfid,idtime,'actual_range',
     &                 NCDOUBLE,2,dvrange,ierr)
        endif
      endif

c     define variable
      id=ncvdef(cdfid,varnam,NCSHORT,ik,did,error)
      if (error.ne.0) goto 922

c     define long_name
      iend=strend(clname)
      call ncaptc(cdfid,id,'long_name',NCCHAR,iend,clname,error)
      if (error.gt.0) goto 923     

c     define units
      iend = strend(clunits)
      call ncaptc(cdfid,id,'units',NCCHAR,iend,clunits,error)
      if (error.gt.0) goto 924     

c     define missing data value
cjsp  call ncapt(cdfid,id,'missing_value',NCFLOAT,1,misdat,error)
      call ncapt(cdfid,id,'missing_value',NCSHORT,1,misdat,error)
      if (error.gt.0) goto 925     

c     define offset_value
      call ncapt(cdfid,id,'add_offset',NCFLOAT,1,offset,error)
      if (error.gt.0) goto 926     

c     define scale_factor
      call ncapt(cdfid,id,'scale_factor',NCFLOAT,1,scale,error)
      if (error.gt.0) goto 927     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 928     

c     synchronise output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write(6,*)'ERROR 920 (NCREDF): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 921  write(6,*)'ERROR 921 (NCDDEF/NCVDEF): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 922  write(6,*)'ERROR 922 (NCVDEF): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 923  write(6,*)'ERROR 923 (LONG NAME): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 924  write(6,*)'ERROR 924 (UNITS): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 925  write(6,*)'ERROR 925 (MISSING DATA): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 926  write(6,*)'ERROR 926 (OFFSET): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 927  write(6,*)'ERROR 927 (SCALE): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
 928  write(6,*)'ERROR 928 (NCENDF): An error occurred while'
     &,' attempting to write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end

      subroutine putdefCOARDSFLT(cdfid,varnam,ndim,misdat,clname
     &         , clunits,offset,scale,vardim,varmin,varmax,error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to define the dimensions and the
c        attributes of a variable on an IVE-NetCDF file for use with the
c        IVE plotting package. Prior to calling this routine, the file must
c        be opened with a call to opncdf (extend an existing file) or
c        crecdf (create a new file).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c        ndim    int   input   the number of dimensions (ndim<=4). 
c                              Upon ndim=4, the fourth dimension of the
c                              variable is specified as 'unlimited'
c                              on the file (time-dimension). It can 
c                              later be extended to arbitrary length.
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              Is dimensioned at least Min(ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error =10   other errors detected.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam,clname,clunits
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  varmin(*), varmax(*)
      real           scale,offset

c     Local variable declarations.
      character *(20) dimnam,dimchk
      character *(20) dimnams(MAXNCDIM)
      character *(5)  rdim(MAXDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),idtime,i,k,ik,ierr
      integer         ncopts

      integer         ibeg,iend
      integer         idcoor
      real            vrange(2)
      real*8          dvrange(2)
      character*(20)  long_name(MAXDIM)
      character*(30)  units(MAXDIM)
      data            rdim /'lon  ','lat  ','level','time '/
      data            long_name /'Longitude', 'Latitude', 
     &                           'Height_Index', 'Time' /

      data            units / 'degrees_east', 'degrees_north',
     &                 'level' , 'hours since 1900-1-1 00:00:0.0' /


c     External function declarations.
      integer         strbeg, strend

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 10
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c         print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define spatial dimensions
      ik = 0 
      do k=1,min0(3,ndim)
       if (vardim(k).gt.1) then
        ik = ik + 1
c       define the default dimension-name
        ibeg=strbeg(rdim(k))
        iend=strend(rdim(k))
        dimnam(1:1+iend-ibeg)=rdim(k)(ibeg:iend)
        did(ik)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a new dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((dimnam(1:3).eq.dimchk(1:3))) then 
              did(ik)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if ((did(ik).lt.0)) then
c         define the dimension and an array with name of dimension
          did(ik)=ncddef(cdfid,dimnam,vardim(k),error)
          idcoor = ncvdef(cdfid,dimnam,NCFLOAT,1,did(ik),ierr)
          if ((error.ne.0).or.(ierr.ne.0)) goto 920
          vrange(1) = varmin(k)
          vrange(2) = varmax(k)
          iend=strend(long_name(k))
          call ncaptc(cdfid,idcoor,'long_name',
     &                 NCCHAR,iend,long_name(k)(1:iend),ierr)
          iend=strend(units(k))
          call ncaptc(cdfid,idcoor,'units',
     &                 NCCHAR,iend,units(k)(1:iend),ierr)
          call ncapt (cdfid,idcoor,'actual_range',
     &                 NCFLOAT,2,vrange,ierr)
        endif
       endif
      enddo

c     define the times-array
      if (ndim.eq.4) then
        ik = ik +1
c       define dimension 'time'
        did(ik)=ncdid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         this dimension must first be defined
          did(ik) = ncddef (cdfid,'time',NCUNLIM,ierr)
        endif
c       define array 'time'
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         define the times-array
          idtime = ncvdef (cdfid,'time',NCDOUBLE,1,did(ik),ierr)
          call ncaptc(cdfid,idtime,'long_name',NCCHAR,4,'Time',ierr)
          iend = strend(units(4))
          call ncaptc(cdfid,idtime,'units',
     &                 NCCHAR,iend,units(4)(1:iend),ierr)
          dvrange(1)=0.
          dvrange(2)=0.
          call ncapt (cdfid,idtime,'actual_range',
     &                 NCDOUBLE,2,dvrange,ierr)
        endif
      endif

c     define variable
cjsp  id=ncvdef(cdfid,varnam,NCSHORT,ik,did,error)
      id=ncvdef(cdfid,varnam,NCFLOAT,ik,did,error)
      if (error.ne.0) goto 920

c     define long_name
      iend=strend(clname)
      call ncaptc(cdfid,id,'long_name',NCCHAR,iend,clname,error)
      if (error.gt.0) goto 920     

c     define units
      iend = strend(clunits)
      call ncaptc(cdfid,id,'units',NCCHAR,iend,clunits,error)
      if (error.gt.0) goto 920     

c     define missing data value
      call ncapt(cdfid,id,'missing_value',NCFLOAT,1,misdat,error)
      if (error.gt.0) goto 920     

c     define offset_value
cjsp  call ncapt(cdfid,id,'add_offset',NCFLOAT,1,offset,error)
cjsp  if (error.gt.0) goto 920     

c     define scale_factor
cjsp  call ncapt(cdfid,id,'scale_factor',NCFLOAT,1,scale,error)
cjsp  if (error.gt.0) goto 920     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 920     

c     synchronise output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end


      subroutine gettimesCOARDS(cdfid,times,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real output  array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  output  number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C------------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

      integer	ierr,i
      real*8    times(*)
      integer   didtim,ntimes

      integer	cdfid,idtime
      integer	ncopts
      character*(20) dimnam

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

      didtim=ncdid(cdfid,'time',ierr)	! inquire id for time dimension
      if (ierr.ne.0) goto 900
      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) goto 900
      call ncdinq(cdfid,didtim,dimnam,ntimes,ierr)      ! inquire # of times
      if (ierr.ne.0) goto 900
  
      do 10 i=1,ntimes
        call ncvgt1(cdfid,idtime,i,times(i),ierr) ! get times
        if (ierr.ne.0) goto 900
   10 continue
  
c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 900  ntimes=1
      times(1)=0.
      call ncpopt (ncopts)
      end


      subroutine puttimesCOARDS(cdfid,times,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real input   array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  input   number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C        Christoph Schaer, ETHZ
C     Note:
C        This preliminary version does not define the times-array, but only
C        overwrites or extends an existing times-array.
C------------------------------------------------------------------------

      integer	ierr,i
      real*8    times(*)
      integer   didtim,ntimes

      integer	cdfid,idtime,nfiltim
      integer	ncdid,ncvid
      character*(20) tnam

      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) return
      didtim=ncdid(cdfid,'time',ierr)	! inquire id for time dimension
      if (ierr.ne.0) return

      call ncdinq(cdfid,didtim,tnam,nfiltim,ierr)   ! inquire # of times
      if (ierr.ne.0) return
      if (nfiltim.lt.ntimes) then
        print *,'Warning: puttimes is extending times-array'
      else if (nfiltim.gt.ntimes) then
        print *,'Warning: puttimes does not cover range of times-array'
      endif

      do 10 i=1,ntimes
        call ncvpt1(cdfid,idtime,i,times(i),ierr)
        if (ierr.ne.0) return
   10 continue
      end


      subroutine putcoords(cdfid,ndim,vardim,xlat1d,xlon1d,sigh,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real input   array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  input   number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C        Christoph Schaer, ETHZ
C     Note:
C        This preliminary version does not define the times-array, but only
C        overwrites or extends an existing times-array.
C------------------------------------------------------------------------

      integer   vardim(*)
      integer	ierr,i,k,ndim,ii
      real      xlat1d(*),xlon1d(*),sigh(*)

      integer	cdfid,idcoor
      integer	ncvid
      real      actval
      character*(5)  rdim(3)

      data            rdim /'lon  ','lat  ','level'/


      do k=1,max0(ndim,3)
        if (vardim(k).gt.1) then
          idcoor = ncvid(cdfid,rdim(k),ierr)
	  if (k.eq.1) then
            if (ierr.eq.0) then
              do i=1,vardim(k)
                actval=xlon1d(i)
                CALL NCVPT1(cdfid,idcoor,i,actval,ierr)
              enddo
            endif
	  else if (k.eq.2) then
            if (ierr.eq.0) then
              do i=1,vardim(k)
                actval=xlat1d(i)
                CALL NCVPT1(cdfid,idcoor,i,actval,ierr)
              enddo
            endif
	  else
            if (ierr.eq.0) then
              do i=1,vardim(k)
		ii = vardim(k)-i+1
		actval=sigh(ii)
                CALL NCVPT1(cdfid,idcoor,i,actval,ierr)
              enddo
            endif
	  end if
        endif
      enddo


      return
      end


c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RCRECDF (filnam,cdfid,varmin,varmax,ndim,ierr)
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      character*(*) filnam
      integer       cdfid,ierr
      real*4        varmin(ndim),varmax(ndim)

      call crecdf (filnam,cdfid,varmin,varmax,3,ierr)

      end

      FUNCTION RMAX(ar,n)
c     ===================
      dimension ar(n)
      rmax=ar(1)
      do i=2,n
        if (ar(i).gt.rmax) rmax=ar(i)
      enddo
      end


      FUNCTION RMIN(ar,n)
c     ===================
      dimension ar(n)
      rmin=ar(1)
      do i=2,n
        if (ar(i).lt.rmin) rmin=ar(i)
      enddo
      end

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRITECDF(cdfid,varnam,arr,ie,je,ke,idim,vtstep,
     &                     vvarmin,vvarmax,xlat1d,xlon1d,izstag,misdat)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c write data in array arr to netcdf-file for timestep vtstep
c 
c arguments an its meaning
c    cdfid     integer            ID returned by crecdf      
c    varnam    character*10       name of the variable to be written
c    arr       real*8(ie,je,ke)   array containing data to be written
c    ie,je,ke  integer            dimensions of arr as used in declar
c    idim      integer(3)         array with used dimensions
c    vtstep    real*8             timestep to be written
c    vvarmin   real*8(3)          minimum values of dim in phys. space
c    vvarmax   real*8(3)          maximum values of dim in phys. space
c    izstag    integer            1 for full levels   0 otherwise
c    vmisdat   real*8             value for missing data


      integer  cdfid,idim(3),izstag,ie,je,ke
      real     arr(ie,je,ke)
      real     xlat1d(je), xlon1d(ie)
      real     vvarmin(3),vvarmax(3)   ! ,vmisdat
      real*8 vtstep
      character*(10)  varnam

      parameter  (ijmax=100000)

      real*4   ufeld(ijmax)

      include '../../Commons/env/include/netcdf.inc'

c     declarations needed for netcdf-stuff
      real*4    varmin(3),varmax(3),varstg(3),misdat,tstep
      integer   strbeg,strend
      integer   idtest,it
      integer   vardim(3)

c *********** declare some auxiliary variables **********

      integer i,j,k,iputlev
      integer ievar,jevar,kevar

c     convert real*8 input variables to real*4 variables
      do i=1,3
        varmin(i)=vvarmin(i)
        varmax(i)=vvarmax(i)
      enddo

cjsp  misdat = vmisdat
      tstep  = vtstep
      ievar  = idim(1)
      jevar  = idim(2)
      kevar  = idim(3)

c     set dimensions used in call to putdef and putdat
      ndims=4
      vardim(1)=ievar
      vardim(2)=jevar
      vardim(3)=kevar
      varstg(1)=0.
      varstg(2)=0.
      varstg(3)=-0.5
      if (izstag.ne.0) varstg(3)=0.

      call ncpopt(NCVERBOS)

c     test if variable varnam already exists
c     if not then define it

      idtest = ncvid(cdfid,varnam(strbeg(varnam):strend(varnam)),it)
      if (it.ne.0) then
        call putdef(cdfid, varnam, ndims, misdat,
     &  vardim, varmin, varmax, varstg, ierr)
        if (ierr.ne.0) stop '*ERROR* in putdef in writecdf'
      endif

c     loop over all vertical levels and write the data
      do k=1,kevar

c     write data from inputarray to outputarray
        ij=0
        do j=1,jevar
          do i=1,ievar
            ij=ij+1
            ufeld(ij)=arr(i,j,k)
          enddo
        enddo
      
        if (kevar.eq.1) then
          iputlev=0
        else
          iputlev=k
        endif

        call putdat(cdfid,varnam(strbeg(varnam):strend(varnam)),
     &         tstep, iputlev, ufeld, ierr) 
        if (ierr.ne.0) then
          print *,'*ERROR* in putdat: variable ',
     &          varnam(strbeg(varnam):strend(varnam))
        endif
      
      enddo   ! loop over vertical levels

cjsp  if (ierr.eq.0) print*,varnam,':  Data written for timestep ',tstep

      end

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WRITECDFI2(cdfid,varnam,arr,ie,je,ke,idim,vtstep
     &         , lname,unit,factor,offset,vvarmin,vvarmax
     &         , xlat1d,xlon1d,sigh,izstag,misdat)

c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c write data in array arr to netcdf-file for timestep vtstep
c 
c arguments an its meaning
c    cdfid     integer            ID returned by crecdf      
c    varnam    character*10       name of the variable to be written
c    arr       real*8(ie,je,ke)   array containing data to be written
c    ie,je,ke  integer            dimensions of arr as used in declar
c    idim      integer(3)         array with used dimensions
c    vtstep    real*8             timestep to be written
c    lname     character*20       long name of variable
c    unit      character*13       units of variable
c    factor    real*8             scaling factor for output
c    offset    real*8             offset for output
c    vvarmin   real*8(3)          minimum values of dim in phys. space
c    vvarmax   real*8(3)          maximum values of dim in phys. space
c    izstag    integer            1 for full levels   0 otherwise
c    vmisdat   real*8             value for missing data


      integer    cdfid,idim(3),izstag,ie,je,ke
      real     arr(ie,je,ke)
      real     xlat1d(je), xlon1d(ie)
      real     vvarmin(3),vvarmax(3)
      real*8     vtstep
      character*(*)  varnam
      character*(*)  lname
      character*(*)  unit

      parameter  (ijmax=100000000)

      integer*2       ufeld(ijmax)

CJSP ADD THIS LINE
      integer*2       mdat
CJSP

      include '../../Commons/env/include/netcdf.inc'

c     declarations needed for netcdf-stuff
      real*4        varmin(3),varmax(3),varstg(3),misdat,tstep,sigh(ke)
      integer       strbeg,strend
      integer       idtest,it
      integer       vardim(3)
      real*4        offset,factor
      real*8        rfac

c *********** declare some auxiliary variables **********

      integer i,j,k,iputlev
      integer ievar,jevar,kevar

c     convert real*8 input variables to real*4 variables
      do i=1,3
        varmin(i)=vvarmin(i)
        varmax(i)=vvarmax(i)
      enddo

      tstep  = vtstep
      ievar  = idim(1)
      jevar  = idim(2)
      kevar  = idim(3)

c     set dimensions used in call to putdef and putdat
      ndims=4
      vardim(1)=ievar
      vardim(2)=jevar
      vardim(3)=kevar
      varstg(1)=0.
      varstg(2)=0.
      varstg(3)=-0.5
      if (izstag.ne.0) varstg(3)=0.

CJSP ADD THIS LINE
c     mdat=nint((misdat-offset)/factor)
      mdat = -32767
CJSP

c     test if variable varnam already exists
c     if not then define it
      call ncpopt(NCVERBOS)


      idtest = ncvid(cdfid,varnam(strbeg(varnam):strend(varnam)),it)
      if (it.ne.0) then

CJSP CHANGE THIS CALL
cjsp    call putdefCOARDS(cdfid, varnam, ndims, misdat, lname, unit
cjsp &            ,  offset, factor, vardim, varmin, varmax, ierr)
        call putdefCOARDS(cdfid, varnam, ndims, mdat, lname, unit
     &            ,  offset, factor, vardim, varmin, varmax, ierr)
CJSP

        if (ierr.ne.0) stop '*ERROR* in putdef in writecdf'
        call putcoords(cdfid,ndims,vardim,xlat1d,xlon1d,sigh,ierr)
        if (ierr.ne.0) stop '*ERROR* in putdef in grbtst'
      endif

c     loop over all vertical levels and write the data
      do k=1,kevar

c     write data from inputarray to outputarray
        rfac = 1./factor
        ij = 0
        do j=1,jevar
          do i=1,ievar
            ij = ij + 1
            if (arr(i,j,k).lt.misdat) then 
               ufeld(ij) = mdat  ! the lowest integer value, should match min
            else 
              ufeld(ij)=nint((arr(i,j,k)-offset)*rfac)
            endif
          enddo
        enddo
      
        if (kevar.eq.1) then
          iputlev=0
        else
          iputlev=k
        endif

        call putdatCOARDS(cdfid,varnam(strbeg(varnam):strend(varnam)),
     &         vtstep, iputlev, ufeld, ierr) 
        if (ierr.ne.0) then
          print *,'*ERROR* in putdat: variable ',
     &          varnam(strbeg(varnam):strend(varnam))
        endif
      
      enddo   ! loop over vertical levels

      end


      subroutine putdat(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to write the data of a variable
c        to an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation), the 
c        variable must be defined with a call to putdef.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be written to the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). If 'time' is not
c                              yet known to the file, a knew time-level is
c                              allocated and appended to the times-array.
c        level   int input     the horizontal level(s) to be written 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call writes the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call writes the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int output    indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is new, but
c                                          appending it would yield a non
c                                          ascending times-array.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c-----------------------------------------------------------------------

      include '../../Commons/env/include/netcdf.inc'

C     Declaration of local variables

      character*(*) varnam
      character*(20) chars
      integer cdfid


      real 	dat(*)
      real	misdat,varmin(3),varmax(3),stag(3)
      real	time, timeval
      data	stag/0.,0.,0./

      integer	corner(4),edgeln(4),did(4),vardim(4),ndims
      integer	error, ierr
      integer	level,ntime
      integer	idtime,idvar,iflag
      integer	i

      integer strend
 
      call ncpopt(NCVERBOS)

c     get definitions of data
      call getdef (cdfid, varnam(1:strend(varnam)), ndims, misdat, 
     &                           vardim, varmin, varmax, stag, ierr)
      if (ierr.ne.0)  print *,'*ERROR* in getdef in putdat'

c     get id of variable
      idvar=ncvid(cdfid,varnam(1:strend(varnam)),ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'

c     get times-array 
      did(4)=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* did(4) in putdat'
      call ncdinq(cdfid,did(4),chars,ntime,ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncdinq in putdat'
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'
C     Check if a new time step is starting
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in putdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then		! new time step
        ntime=ntime+1
        iflag=ntime
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvid in putdat'
        call ncvpt1(cdfid,idtime,ntime,time,ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvpt1 in putdat'
      endif

C     Define data volume to write on the NetCDF file in index space
      corner(1)=1               ! starting corner of data volume
      corner(2)=1
      edgeln(1)=vardim(1)       ! edge lengthes of data volume
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        corner(3)=1
        edgeln(3)=vardim(3)
      else
        corner(3)=level
        edgeln(3)=1
      endif
      corner(4)=iflag
      edgeln(4)=1
      
C     Put data on NetCDF file

c      print *,'vor Aufruf ncvpt d.h. Daten schreiben in putdat '
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln

      call ncvpt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvpt in putdat - Put data on NetCDF file'
      endif

C     Synchronize output to disk and close the files

      call ncsnc(cdfid,ierr)
      if (ierr.ne.0) print *, '*ERROR* in ncsnc in putdat'
      end



      subroutine putdef (cdfid, varnam, ndim, misdat, 
     &                            vardim, varmin, varmax, stag, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to define the dimensions and the
c        attributes of a variable on an IVE-NetCDF file for use with the
c        IVE plotting package. Prior to calling this routine, the file must
c        be opened with a call to opncdf (extend an existing file) or
c        crecdf (create a new file).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c        ndim    int   input   the number of dimensions (ndim<=4). 
c                              Upon ndim=4, the fourth dimension of the
c                              variable is specified as 'unlimited'
c                              on the file (time-dimension). It can 
c                              later be extended to arbitrary length.
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              Is dimensioned at least Min(ndim). 
c        stag    real  input   the grid staggering for each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error =10   other errors detected.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*)

c     Local variable declarations.
      character *(20) dimnam,attnam,dimchk
      character *(1)  chrid(MAXDIM)
      character *(20) dimnams(MAXNCDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),idtime,i,k,ierr
      integer         ncopts
      integer         ibeg,iend
      data            chrid/'x','y','z','t'/

c     External function declarations.
      integer         strbeg, strend

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 10
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c         print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define spatial dimensions
      ibeg=strbeg(varnam)
      iend=strend(varnam)
      do k=1,min0(3,ndim)
c       define the default dimension-name
        dimnam(1:3)='dim'
        dimnam(4:4)=chrid(k)
        dimnam(5:5)='_'
        dimnam(6:6+iend-ibeg)=varnam(ibeg:iend)
        dimnam=dimnam(1:6+iend-ibeg)
        did(k)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a new dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((vardim(k).eq.dimvals(i)).and.
     &        (dimnam(1:4).eq.dimchk(1:4))) then 
              did(k)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if (did(k).lt.0) then
c         define the dimension
          did(k)=ncddef(cdfid,dimnam,vardim(k),error)
          if (error.ne.0) goto 920
        endif
      enddo

c     define the times-array
      if (ndim.eq.4) then
c       define dimension 'time'
        did(4)=ncdid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         this dimension must first be defined
          did(4) = ncddef (cdfid,'time',NCUNLIM,ierr)
        endif
c       define array 'time'
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) then 
c         define the times-array
          idtime = ncvdef (cdfid,'time',NCFLOAT,1,did(4),ierr)
        endif
      endif

c     define variable
      id=ncvdef(cdfid,varnam,NCFLOAT,ndim,did,error)
      if (error.ne.0) goto 920

c     define attributes
      do k=1,min0(ndim,3)
c       min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmin(k),error)
        if (error.gt.0) goto 920
c       max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmax(k),error)
        if (error.gt.0) goto 920     
c       staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,stag(k),error)
        if (error.gt.0) goto 920
      enddo

c     define missing data value
      call ncapt(cdfid,id,'missing_data',NCFLOAT,1,misdat,error)
      if (error.gt.0) goto 920     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 920     

c     synchronyse output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end



      subroutine clscdf (cdfid, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine closes an open netCDF file.
c     Aguments :
c        cdfid  int  input   the id of the file to be closed.
c        error  int  output  indicates possible errors found in this
c                            routine.
c                            error = 0   no errors detected.
c                            error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer      cdfid, error

c     Local variable declarations.
      integer      ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     Make sure netCDF errors do not abort execution.
      call ncpopt (NCVERBOS)

c     Close requested file.
      call ncclos (cdfid, error)

c     Reset error options.
      call ncpopt (ncopts)

      end


      subroutine crecdf (filnam, cdfid, phymin, phymax, ndim, error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to create a netCDF file for use with 
c        the UWGAP plotting package.
c           Any netCDF file written to must be closed with the call
c        'call clscdf(cdfid,error)', where cdfid and error are
c        as in the argumentlist below. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  input   the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  input   the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   input   the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        cfn     char  input   constants file name 
c                              ('0' = no constants file).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created crecdf.
c-----------------------------------------------------------------------

      include  '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, error
      character *(*) filnam
      real           phymin(*), phymax(*)

c     Local variable declarations.
      character *(20) attnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, k, ncopts
      data            chrid/'x','y','z','a'/

c     External function declarations.
      integer        strend

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     create the netCDF file
      cdfid = nccre (filnam(1:strend(filnam)), NCCLOB, error)
      if (error.ne.0) go to 920

c     define global attributes
      do k=1,ndim
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymin(k),error)
        if (error.gt.0) goto 920

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymax(k),error)
        if (error.gt.0) goto 920
      enddo

c     End variable definitions.
      call ncendf (cdfid, error)
      if (error.gt.0) goto 920

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'create the data file in subroutine crecdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end




      subroutine getdef (cdfid, varnam, ndim, misdat, 
     &                              vardim, varmin, varmax, stag, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions and attributes of 
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              Is dimensioned at least (ndim). 
c        varmin  real  output  the location in physical space of the
c                              origin of each variable. 
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  output  the extend of each variable in physical
c                              space.
c                              Is dimensioned at least Min(3,ndim). 
c        stag    real  output  the grid staggering for each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include '../../Commons/env/include/netcdf.inc'

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*)

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),attnam,vnam
      character *(1)  chrid(MAXDIM)
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts
      data            chrid/'x','y','z','t'/

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCFLOAT) error=1
      if (error.gt.0) goto 920

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get attributes
      do k=1,min0(3,ndim)
c       get min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmin(k),error)
        if (error.gt.0) goto 920
c       get max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmax(k),error)
        if (error.gt.0) goto 920     
c       get staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncagt(cdfid,id,attnam,stag(k),error)
        if (error.gt.0) goto 920
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_data',misdat,error)
      if (error.gt.0) goto 920     

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) '*ERROR*: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end

