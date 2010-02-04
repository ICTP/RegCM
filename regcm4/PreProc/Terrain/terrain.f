       PROGRAM TERRAIN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   TERRAIN is the first component of the REGional Climate Modeling    !
!   (RegCM) system version 3.0 and used to access archived terrain     !
!   height and landuse charactistics data at regular latitude-         !
!   longititude intervals and interpolate to the mesoscale grid for    !
!   a specified map projection.                                        !
!                                                                      !
!                                     PWC group, Abdus Salam ICTP      !
!                                                   May. 27, 2006      !
!                                                                      !
!   The authors wish to acknowledge the authors of NCAR MM4/5          !
!   terrain codes, their works is our code base.                       !
!                                                                      !
!       MM4 terrain code: A. Mcnab, T. Tarbell and N. Seaman           !
!                         C. Larkin and N. Seaman  (before 1986)       !
!       MM5 terrain code: Yong-Run Guo and Sue Chen  10/21/1993        !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  In the present version of RegCM, all the variables, which directly
!  related to map factors, are just calculated and stored in TERRAIN
!  for later use.
!
!  This program reads terrain height data from the NCAR air force
!  terrain tapes (1 deg., 30', or 5') or terrain and landuse data
!  from the PSU/NCAR combined landuse tapes (1 deg., 30', 10'), or 10'
!  GLCC landuse data (Loveland et al 1999) and analyzes heights and/or
!  landuse values to a given grid.
!     
!  terrain will be read from unit 12, landuse from unit 10 or 11, 
!  output variables unit 9.
!---------------------------------------------------------------------
      implicit none
      include 'domain.param'
!---------------------------------------------------------------------
! note: internally set params!!!!!
!
!   iblk = dimension of arrays xobs,yobs,ht,htsd
!     estimate iblk = ihmax*jhmax where ihmax = (xmaxlat-xminlat)/xnc
!     and jhmax = (xmaxlon-xminlon)/xnc.  xnc = 1,0.5,1./6.or 5./60. for
!     1 deg, 30 min, 10 min or 5 min data, respectively.  add 1-2000 to
!     estimate of iblk for safety.  if you underestimate the
!     dimensions the program will abort with a diagnostic message
!     indicating the correct dimensions required .
!   iter,jter = dimensions of array lnd8.
!     within the search region. at present, iter and jter must be
!     equal and >= max(ihmax,jhmax) where ihmax and jhmax are
!     as calculated above.
!
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay
      real*8  lnd8(iter,jter)
      common /block8/lnd8
!--------------------------------------------------------------------
      real*4  xlat(iy,jx),xlon(iy,jx),xmap(iy,jx),coriol(iy,jx)
      real*4  dlat(iy,jx),dlon(iy,jx),dmap(iy,jx)
      real*4  htsdgrid(iy,jx)
      real*4  htgrid(iy,jx),lndout(iy,jx),texout(iy,jx)
      real*4  sanda(iy,jx),sandb(iy,jx),claya(iy,jx),clayb(iy,jx)
      real*4  frac_lnd(iy,jx,nveg),frac_tex(iy,jx,ntex)
      character*1 ch(iy,jx)
      real*4  snowam(iy,jx),mask(iy,jx)
      integer lnduse(iy,jx),intext(iy,jx)
      common /maps/ xlat,xlon,xmap,coriol,dlat,dlon,dmap
     &     , htgrid,htsdgrid,sanda,sandb,claya,clayb
     &     , frac_lnd,frac_tex,snowam,mask,lndout,texout
     &     , lnduse,intext
      real*4  sigma(kz+1)
      common /verts/ sigma
!--------------------------------------------------------------------
      real*4  xlat_s(iy*nsg,jx*nsg),xlon_s(iy*nsg,jx*nsg)
      real*4  dlat_s(iy*nsg,jx*nsg),dlon_s(iy*nsg,jx*nsg)
      real*4  xmap_s(iy*nsg,jx*nsg),dmap_s(iy*nsg,jx*nsg)
      real*4  coriol_s(iy*nsg,jx*nsg)
      real*4  htsdgrid_s(iy*nsg,jx*nsg),htgrid_s(iy*nsg,jx*nsg)
      real*4  lndout_s(iy*nsg,jx*nsg),texout_s(iy*nsg,jx*nsg)
      real*4  sanda_s(iy*nsg,jx*nsg),sandb_s(iy*nsg,jx*nsg)
      real*4  claya_s(iy*nsg,jx*nsg),clayb_s(iy*nsg,jx*nsg)
      real*4  frac_lnd_s(iy*nsg,jx*nsg,nveg)
      real*4  frac_tex_s(iy*nsg,jx*nsg,ntex)
      character*1 ch_s(iy*nsg,jx*nsg)
      real*4  snowam_s(iy*nsg,jx*nsg),mask_s(iy*nsg,jx*nsg)
      integer lnduse_s(iy*nsg,jx*nsg),intext_s(iy*nsg,jx*nsg)
      common /maps_s/ xlat_s,xlon_s,xmap_s,coriol_s,dlat_s,dlon_s,dmap_s
     &     , htgrid_s,htsdgrid_s,sanda_s,sandb_s,claya_s,clayb_s
     &     , frac_lnd_s,frac_tex_s,snowam_s,mask_s,lndout_s,texout_s
     &     , lnduse_s,intext_s
!--------------------------------------------------------------------
      integer nunitc,nunitc_s
      common /icontr/nunitc,nunitc_s
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc

      real*4  hscr1(iy,jx),land(iy,jx,2),itex(iy,jx,2)
      real*4  hscr1_s(iy*nsg,jx*nsg)
      real*4  land_s(iy*nsg,jx*nsg,2),itex_s(iy*nsg,jx*nsg,2)
      logical ibndry
      real*4  dxcen,dycen
      real*4  corc(iy,jx),sumc(iy,jx)
      real*4  wtmaxc(iy,jx),htsavc(iy,jx)
      integer nsc(iy,jx)
      real*4  corc_s(iy*nsg,jx*nsg),sumc_s(iy*nsg,jx*nsg)
      real*4  wtmaxc_s(iy*nsg,jx*nsg),htsavc_s(iy*nsg,jx*nsg)
      integer nsc_s(iy*nsg,jx*nsg)
      character*50 filout_s ! Terrain Output filename subgrid
      character*50 filctl_s ! Grads control filename for subgrid output
      character*10 CHAR_LND,CHAR_TEX
      real*4  clong
      common /addstack/dxcen,dycen
      common /addstack2/hscr1,land,itex,hscr1_s,land_s,itex_s
     &       , corc,sumc,wtmaxc,htsavc,corc_s,sumc_s,wtmaxc_s,htsavc_s
     &       , clong,nsc,nsc_s
      integer i,j,k

!
      dxcen = 0.0
      dycen = 0.0
!
      clong = clon
      if(clong.gt. 180.) clong=clong-360.
      if(clong.le.-180.) clong=clong+360.
!
 1011 format(A18,I1,A5)
 1012 format(A18,I2,A5)
 1013 format(A18,I1,A4)
 1014 format(A18,I2,A4)
      nunitc=9
      if(nsg.gt.1) then
        nunitc_s=19
        if(nsg.lt.10) then
          write(filout_s,1011) filout(1:18),nsg,filout(19:23)
          write(filctl_s,1013) filctl(1:18),nsg,filctl(19:22)
        else
          write(filout_s,1012) filout(1:18),nsg,filout(19:23)
          write(filctl_s,1014) filctl(1:18),nsg,filctl(19:22)
        endif
        call setup(nunitc_s,iy*nsg,jx*nsg,ntypec_s,nveg,iproj,ds/nsg
     &     , clat,clong,igrads,ibyte,filout_s,filctl_s)
        if (iproj.eq.'LAMCON') then
          CALL LAMBRT(XLON_S,XLAT_S,XMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,0,XN,TRUELATL,TRUELATH)
          CALL LAMBRT(DLON_S,DLAT_S,DMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,1,XN,TRUELATL,TRUELATH)
          write(*,*)'XN,TRUELATL,TRUELATH = ',XN,TRUELATL,TRUELATH
        elseif (iproj.eq.'POLSTR') then
          CALL MAPPOL(XLON_S,XLAT_S,XMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,0)
          CALL MAPPOL(DLON_S,DLAT_S,DMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,1)
          xn=1.
        elseif (iproj.eq.'NORMER') then
          CALL NORMER(XLON_S,XLAT_S,XMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,0)
          CALL NORMER(DLON_S,DLAT_S,DMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLONG,CLAT,dsinm,1)
          xn=0.
        elseif (iproj.eq.'ROTMER') then
          CALL ROTMER(XLON_S,XLAT_S,XMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLON,CLAT,PLON,PLAT,dsinm,0)
          CALL ROTMER(DLON_S,DLAT_S,DMAP_S,CORIOL_S,iy*nsg,jx*nsg
     &        ,CLON,CLAT,PLON,PLAT,dsinm,1)
          xn=0.
        else
          print*,'iproj MAP PROJECTION IS NOT AN OPTION'
          stop 999
        endif
        print*, 'after calling MAP PROJECTION, for subgrid'
!
!       reduce the search area for the domain
!              [minlat:maxlat,minlon:maxlon]
        CALL MXMNLL(iy*nsg,jx*nsg,CLONG,XLON_S,XLAT_S,ntypec_s)
        print*, 'after calling MXMNLL, for subgrid'
!
!       read in the terrain & landuse data
        CALL RDLDTR_s
        print*, 'after calling RDLDTR_s, for subgrid'
        if (ifanal) then
!       convert xobs and yobs from LON and LAT to x and y in mesh
          CALL XYOBSLL(iy*nsg,jx*nsg,iproj,clat,clong,plat,plon
     &                ,truelatL,truelatH)
          print*, 'after calling XYOBSLL, for subgrid'
!
!       create the terrain height fields
          call anal2(htsdgrid_s,ht2,nobs,iy*nsg,jx*nsg
     &        ,ntypec_s,corc_s,sumc_s,nsc_s,wtmaxc_s,htsavc_s)
          call anal2(htgrid_s,ht,nobs,iy*nsg,jx*nsg
     &        ,ntypec_s,corc_s,sumc_s,nsc_s,wtmaxc_s,htsavc_s)
          print*, 'after calling ANAL2, for subgrid'
          do j=1,jx*nsg
          do i=1,iy*nsg
            htgrid_s(i,j) = amax1(htgrid_s(i,j)*100.,0.0)
            htsdgrid_s(i,j)=amax1(htsdgrid_s(i,j)*100000.,0.0)
            htsdgrid_s(i,j)=
     &                 sqrt(amax1(htsdgrid_s(i,j)-htgrid_s(i,j)**2,0.0))
          enddo
          enddo
        else
          call interp_s
          print*, 'after calling INTERP, for subgrid'
!         print*, '  Note that the terrain standard deviation is'
!         print*, '  underestimated using INTERP. (I dont know why?)'
        endif
!       create surface landuse types
        call surf(xlat_s,xlon_s,lnduse_s,iy*nsg,jx*nsg,iter,jter,nnc,xnc
     &      ,lndout_s,land_s,lnd8,nobs,grdltmn,grdlnmn,h2opct
     &      ,LSMTYP,sanda_s,sandb_s,claya_s,clayb_s,frac_lnd_s,nveg
     &      ,AERTYP,intext_s,texout_s,frac_tex_s,ntex)
        print*, 'after calling SURF, for subgrid'
! **** Adjust the Great Lake Heights to their actual values.
        if (lakadj) then
        print*,'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
        CALL LAKEADJ(LSMTYP,lnduse_s,htgrid_s,xlat_s,xlon_s
     &               ,iy*nsg,jx*nsg)
        print*, 'after calling LAKEADJ, for subgrid'
        end if
        CALL SMTH121(htgrid_s,iy*nsg,jx*nsg,hscr1_s)
        CALL SMTH121(htsdgrid_s,iy*nsg,jx*nsg,hscr1_s)
! **** Readjust the Great Lake Heights to their actual values again.
        if (lakadj) then
        print*,'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
        CALL LAKEADJ(LSMTYP,lnduse_s,htgrid_s,xlat_s,xlon_s
     &               ,iy*nsg,jx*nsg)
        print*, 'after calling LAKEADJ, for subgrid'
        end if
        ibndry=.TRUE.
        if (ibndry) then
          do j=2,jx*nsg-1
            htgrid_s(1,j)      = htgrid_s(2,j)
            htgrid_s(iy*nsg,j) = htgrid_s(iy*nsg-1,j)
            lnduse_s(1,j)      = lnduse_s(2,j)
            lnduse_s(iy*nsg,j) = lnduse_s(iy*nsg-1,j)
            lndout_s(1,j)      = lndout_s(2,j)
            lndout_s(iy*nsg,j) = lndout_s(iy*nsg-1,j)

            if(LSMTYP.eq.'USGS') then
              sanda_s(1, j)     = sanda_s(2,   j)
              sanda_s(iy*nsg,j) = sanda_s(iy*nsg-1,j)
              sandb_s(1, j)     = sandb_s(2,   j)
              sandb_s(iy*nsg,j) = sandb_s(iy*nsg-1,j)
              claya_s(1, j)     = claya_s(2,   j)
              claya_s(iy*nsg,j) = claya_s(iy*nsg-1,j)
              clayb_s(1, j)     = clayb_s(2,   j)
              clayb_s(iy*nsg,j) = clayb_s(iy*nsg-1,j)
              do k = 1, nveg
                frac_lnd_s( 1,j,k)    = frac_lnd_s(2,   j,k)
                frac_lnd_s(iy*nsg,j,k) = frac_lnd_s(iy*nsg-1,j,k)
              enddo
            endif
            if(AERTYP(7:7).eq.'1') then
              intext_s(1,j)     = intext_s(2,j)
              intext_s(iy*nsg,j)= intext_s(iy*nsg-1,j)
              texout_s(1,j)     = texout_s(2,j)
              texout_s(iy*nsg,j)= texout_s(iy*nsg-1,j)
              do k = 1, ntex
                frac_tex_s( 1,j,k)     = frac_tex_s(2,   j,k)
                frac_tex_s(iy*nsg,j,k) = frac_tex_s(iy*nsg-1,j,k)
              enddo
            endif
          enddo
          do i=1,iy*nsg
            htgrid_s(i,1)     = htgrid_s(i,2)
            htgrid_s(i,jx*nsg)= htgrid_s(i,jx*nsg-1)
            lnduse_s(i,1)     = lnduse_s(i,2)
            lnduse_s(i,jx*nsg)= lnduse_s(i,jx*nsg-1)
            lndout_s(i,1)     = lndout_s(i,2)
            lndout_s(i,jx*nsg)= lndout_s(i,jx*nsg-1)

            if(LSMTYP.eq.'USGS') then
              sanda_s(i, 1)     = sanda_s(i,2)
              sanda_s(i,jx*nsg) = sanda_s(i,jx*nsg-1)
              sandb_s(i, 1)     = sandb_s(i,2)
              sandb_s(i,jx*nsg) = sandb_s(i,jx*nsg-1)
              claya_s(i, 1)     = claya_s(i,2)
              claya_s(i,jx*nsg) = claya_s(i,jx*nsg-1)
              clayb_s(i, 1)     = clayb_s(i,2)
              clayb_s(i,jx*nsg) = clayb_s(i,jx*nsg-1)
              do k = 1, nveg
                frac_lnd_s(i, 1,k)     = frac_lnd_s(i,2,k)
                frac_lnd_s(i,jx*nsg,k) = frac_lnd_s(i,jx*nsg-1,k)
              enddo
            endif
            if(AERTYP(7:7).eq.'1') then
              intext_s(i,1)     = intext_s(i,2)
              intext_s(i,jx*nsg)= intext_s(i,jx*nsg-1)
              texout_s(i,1)     = texout_s(i,2)
              texout_s(i,jx*nsg)= texout_s(i,jx*nsg-1)
              do k = 1, ntex
                frac_tex_s(i, 1,k)     = frac_tex_s(i,2,k)
                frac_tex_s(i,jx*nsg,k) = frac_tex_s(i,jx*nsg-1,k)
              enddo
            endif
          enddo
        endif
        do i=1,iy*nsg
        do j=1,jx*nsg
          snowam_s(i,j) = 0.0
        enddo
        enddo
!       land/sea mask fudging
        if(nsg.ge.10) then
          write(CHAR_LND,1122) 'LANDUSE_',NSG
          write(CHAR_TEX,1122) 'TEXTURE_',NSG
        else
          write(CHAR_LND,1121) 'LANDUSE_',NSG
          write(CHAR_TEX,1121) 'TEXTURE_',NSG
        endif
 1122   format(A8,I2)
 1121   format(A8,I1)
        CALL LNDFUDGE(FUDGE_LND_s,ch_s,lndout_s,htgrid_s,iy*nsg,jx*nsg
     &               ,LSMTYP,CHAR_LND)
        IF(AERTYP(7:7).eq.'1') 
     &  CALL TEXFUDGE(FUDGE_TEX_s,ch_s,texout_s,htgrid_s,iy*nsg,jx*nsg
     &               ,CHAR_TEX)
        print*, 'after calling FUDGE, for subgrid'
!       output terrestrial fields
!all    CALL OUTPUT(nunitc_s,iy*nsg,jx*nsg,1,dsinm,clat,clong,plat,plon
!all &    ,iproj,htgrid_s,htsdgrid_s,lndout_s,xlat_s,xlon_s
!all &    ,dlat_s,dlon_s,xmap_s,DATTYP,dmap_s,coriol_s,snowam_s
!all &    ,igrads,ibigend,kz,sigma,mask_s,ptop
!all &    ,htgrid_s,lndout_s,ibyte,nsg,truelatL,truelatH,XN,filout_s
!all &    ,LSMTYP,sanda_s,sandb_s,claya_s,clayb_s,frac_lnd_s,nveg
!all &    ,AERTYP,texout_s,frac_tex_s,ntex)
!all    print*, 'after calling OUTPUT, for subgrid'
!
!       OUTPUT2 is used to output also the fraction of each
!                LANDUSE legend and TEXTURE type
        CALL OUTPUT2(nunitc_s,iy*nsg,jx*nsg,1,dsinm,clat,clong,plat,plon
     &    ,iproj,htgrid_s,htsdgrid_s,lndout_s,xlat_s,xlon_s
     &    ,dlat_s,dlon_s,xmap_s,DATTYP,dmap_s,coriol_s,snowam_s
     &    ,igrads,ibigend,kz,sigma,mask_s,ptop
     &    ,htgrid_s,lndout_s,ibyte,nsg,truelatL,truelatH,XN,filout_s
     &    ,LSMTYP,sanda_s,sandb_s,claya_s,clayb_s,frac_lnd_s,nveg
     &    ,AERTYP,texout_s,frac_tex_s,ntex)
        print*, 'after calling OUTPUT, for subgrid'
      endif
!
      dxcen = 0.0
      dycen = 0.0
!
!     set up the parameters and constants
      call setup(nunitc,iy,jx,ntypec,nveg,iproj,ds,clat,clong,igrads
     &   , ibyte,filout,filctl)
      print*, 'after calling SETUP'
!
!-----calling the map projection subroutine
      if (iproj.eq.'LAMCON') then
        CALL LAMBRT(XLON,XLAT,XMAP,CORIOL,iy,jx,CLONG,CLAT,dsinm,0
     &             ,XN,TRUELATL,TRUELATH)
        CALL LAMBRT(DLON,DLAT,DMAP,CORIOL,iy,jx,CLONG,CLAT,dsinm,1
     &             ,XN,TRUELATL,TRUELATH)
        write(*,*)'XN,TRUELATL,TRUELATH = ',XN,TRUELATL,TRUELATH
      elseif (iproj.eq.'POLSTR') then
        CALL MAPPOL(XLON,XLAT,XMAP,CORIOL,iy,jx,CLONG,CLAT,dsinm,0)
        CALL MAPPOL(DLON,DLAT,DMAP,CORIOL,iy,jx,CLONG,CLAT,dsinm,1)
        xn=1.
      elseif (iproj.eq.'NORMER') then
        CALL NORMER(XLON,XLAT,XMAP,CORIOL,IY,JX,CLONG,CLAT,dsinm,0)
        CALL NORMER(DLON,DLAT,DMAP,CORIOL,IY,JX,CLONG,CLAT,dsinm,1)
        xn=0.
      elseif (iproj.eq.'ROTMER') then
        CALL ROTMER(XLON,XLAT,XMAP,CORIOL,iy,jx,CLONG,CLAT
     &             ,PLON,PLAT,dsinm,0)
        CALL ROTMER(DLON,DLAT,DMAP,CORIOL,iy,jx,CLONG,CLAT
     &             ,PLON,PLAT,dsinm,1)
        xn=0.
      else
        print*,'iproj MAP PROJECTION IS NOT AN OPTION'
        stop 999
      endif
      print*, 'after calling MAP PROJECTION'
!
!     reduce the search area for the domain 
!            [minlat:maxlat,minlon:maxlon]
      CALL MXMNLL(IY,JX,CLONG,XLON,XLAT,ntypec)
      print*, 'after calling MXMNLL'
!
!     read in the terrain & landuse data
      CALL RDLDTR
      print*, 'after calling RDLDTR'

!     compute the scaled standard deviation of terrain height.
!     (must be called before XYOBSLL because xobs/yobs are modified)
!     CALL SCALESD
!     print*, 'after calling SCALESD'
 
      if (ifanal) then
!     convert xobs and yobs from LON and LAT to x and y in mesh
        CALL XYOBSLL(iy,jx,iproj,clat,clong,plat,plon
     &              ,truelatL,truelatH)
        print*, 'after calling XYOBSLL'
 
!     create the terrain height fields
        call anal2(htsdgrid,ht2,nobs,iy,jx,
     &             ntypec,corc,sumc,nsc,wtmaxc,htsavc)
        call anal2(htgrid,ht,nobs,iy,jx,
     &             ntypec,corc,sumc,nsc,wtmaxc,htsavc)
        print*, 'after calling ANAL2'
        do j=1,jx
        do i=1,iy
          htgrid(i,j) = amax1(htgrid(i,j)*100.,0.0)
          htsdgrid(i,j)=amax1(htsdgrid(i,j)*100000.,0.0)
          htsdgrid(i,j)=sqrt(amax1(htsdgrid(i,j)-htgrid(i,j)**2,0.0))
        enddo
        enddo
      else
        call interp
!       print*, 'after calling INTERP'
!       print*, '  Note that the terrain standard deviation is'
!       print*, '  underestimated using INTERP. (I dont know why!)'
      endif


!     create surface landuse types
      call surf(xlat,xlon,lnduse,iy,jx,iter,jter,nnc,xnc
     &         ,lndout,land,lnd8,nobs,grdltmn,grdlnmn,h2opct
     &         ,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
     &         ,AERTYP,intext,texout,frac_tex,ntex)
      print*, 'after calling SURF'

! **** Adjust the Great Lake Heights to their actual values.
      if (lakadj) then
        print*,'CALLING LAKEADJ FOR THE FIRST TIME (before 2dx pass)'
        CALL LAKEADJ(LSMTYP,lnduse,htgrid,xlat,xlon,iy,jx)
        print*, 'after calling LAKEADJ'
      end if

! ******           preliminary heavy smoothing of boundaries
      if (smthbdy) call smthtr( htgrid, iy, jx )

! ******           grell smoothing to eliminate 2 delx wave (6/90):
      CALL SMTH121(htgrid,iy,jx,hscr1)
      CALL SMTH121(htsdgrid,iy,jx,hscr1)

! **** Readjust the Great Lake Heights to their actual values again.
      if (lakadj) then
        print*,'CALLING LAKEADJ FOR THE SECOND TIME (after 2dx pass)'
        CALL LAKEADJ(LSMTYP,lnduse,htgrid,xlat,xlon,iy,jx)
      end if

      ibndry=.TRUE.
      if (ibndry) then
         do j=2,jx-1
            htgrid(1,j) = htgrid(2,j)
            htgrid(iy,j)= htgrid(iy-1,j)
            lnduse(1,j) = lnduse(2,j)
            lnduse(iy,j)= lnduse(iy-1,j)
            lndout(1,j) = lndout(2,j)
            lndout(iy,j)= lndout(iy-1,j)

            if(LSMTYP.eq.'USGS') then
               sanda(1, j) = sanda(2,   j)
               sanda(iy,j) = sanda(iy-1,j)
               sandb(1, j) = sandb(2,   j)
               sandb(iy,j) = sandb(iy-1,j)
               claya(1, j) = claya(2,   j)
               claya(iy,j) = claya(iy-1,j)
               clayb(1, j) = clayb(2,   j)
               clayb(iy,j) = clayb(iy-1,j)
               do k = 1, nveg
                  frac_lnd( 1,j,k) = frac_lnd(2,   j,k)
                  frac_lnd(iy,j,k) = frac_lnd(iy-1,j,k)
               enddo
            endif
            if(AERTYP(7:7).eq.'1') then
               intext(1,j) = intext(2,j)
               intext(iy,j)= intext(iy-1,j)
               texout(1,j) = texout(2,j)
               texout(iy,j)= texout(iy-1,j)
               do k = 1, ntex
                  frac_tex( 1,j,k) = frac_tex(2,   j,k)
                  frac_tex(iy,j,k) = frac_tex(iy-1,j,k)
               enddo
            endif
         enddo
         do i=1,iy
            htgrid(i,1) = htgrid(i,2)
            htgrid(i,jx)= htgrid(i,jx-1)
            lnduse(i,1) = lnduse(i,2)
            lnduse(i,jx)= lnduse(i,jx-1)
            lndout(i,1) = lndout(i,2)
            lndout(i,jx)= lndout(i,jx-1)

            if(LSMTYP.eq.'USGS') then
               sanda(i, 1) = sanda(i,2)
               sanda(i,jx) = sanda(i,jx-1)
               sandb(i, 1) = sandb(i,2)
               sandb(i,jx) = sandb(i,jx-1)
               claya(i, 1) = claya(i,2)
               claya(i,jx) = claya(i,jx-1)
               clayb(i, 1) = clayb(i,2)
               clayb(i,jx) = clayb(i,jx-1)
               do k = 1, nveg
                  frac_lnd(i, 1,k) = frac_lnd(i,2,k)
                  frac_lnd(i,jx,k) = frac_lnd(i,jx-1,k)
               enddo
            endif
            if(AERTYP(7:7).eq.'1') then
               intext(i,1) = intext(i,2)
               intext(i,jx)= intext(i,jx-1)
               texout(i,1) = texout(i,2)
               texout(i,jx)= texout(i,jx-1)
               do k = 1, ntex
                  frac_tex(i, 1,k) = frac_tex(i,2,k)
                  frac_tex(i,jx,k) = frac_tex(i,jx-1,k)
               enddo
            endif
         enddo
      endif

      do i=1,iy
      do j=1,jx
         snowam(i,j) = 0.0
      enddo
      enddo
!     land/sea mask fudging
      CHAR_LND='LANDUSE'
      CHAR_TEX='TEXTURE'
      CALL LNDFUDGE(FUDGE_LND,ch,lndout,htgrid,iy,jx,LSMTYP,CHAR_LND)
      IF(AERTYP(7:7).eq.'1') 
     &CALL TEXFUDGE(FUDGE_TEX,ch,texout,htgrid,iy,jx,CHAR_TEX)
      print*, 'after calling FUDGE'
!     output terrestrial fields
!all  CALL OUTPUT(nunitc,iy,jx,nsg,dsinm,clat,clong,plat,plon,iproj
!all &    ,htgrid,htsdgrid,lndout,xlat,xlon,dlat,dlon,xmap,DATTYP
!all &    ,dmap,coriol,snowam,igrads,ibigend,kz,sigma,mask,ptop
!all &    ,htgrid_s,lndout_s,ibyte,nsg,truelatL,truelatH,XN,filout
!all &    ,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
!all &    ,AERTYP,texout,frac_tex,ntex)
!all  print*, 'after calling OUTPUT'
!
!     OUTPUT2 is used to output also the fraction of each
!                LANDUSE legend and TEXTURE type
      CALL OUTPUT2(nunitc,iy,jx,nsg,dsinm,clat,clong,plat,plon,iproj
     &    ,htgrid,htsdgrid,lndout,xlat,xlon,dlat,dlon,xmap,DATTYP
     &    ,dmap,coriol,snowam,igrads,ibigend,kz,sigma,mask,ptop
     &    ,htgrid_s,lndout_s,ibyte,nsg,truelatL,truelatH,XN,filout
     &    ,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
     &    ,AERTYP,texout,frac_tex,ntex)
      print*, 'after calling OUTPUT'

!     prepare domain and time parameters for ICBC step
      CALL FORICBC(iy,jx,kz,nsg,IDATE1,IDATE2,igrads,ibigend,ibyte
     &            ,DATTYP,SSTTYP,EHSO4,LSMTYP,AERTYP,filout)

!     prepare domain parameters for RegCM step
      CALL FORMAIN(iy,jx,kz,nsg,ibyte,DATTYP,EHSO4,LSMTYP,AERTYP,NPROC)
!
      stop 9999
      end
!
      subroutine anal2(a2,asta,nsta,iy,jx,ntype,cor,sum,ns,wtmax,htsav)
      implicit none
!
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      integer nsta
      real*4  asta(nsta)
      integer iy,jx,ntype
      real*4  dxcen,dycen
      common /addstack/ dxcen,dycen
      real*4  a2(iy,jx),htsav(iy,jx),cor(iy,jx),
     &        sum(iy,jx),wtmax(iy,jx)
      integer ns(iy,jx)
!
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
!
      real*4  deltas,xcntr,ycntr,dy,dx,ris,ris2
      integer i,j,kk,ie,je,nscan,nskip
      real*4  dxobs,dyobs,rjobs,riobs,ymaxi,ymini,xmaxj,xminj
      real*4  x,y,rx,ry,rsq,wt
      integer maxi,mini,maxj,minj
!
!     objective analysis to fill a grid based on observations
!     xobs and yobs are x and y positions on observations, not
!     necessarily grid points.
!
      ie     = iy-1
      je     = jx-1
      nscan  = 1
      deltas = dsinm
      print 10,rin,deltas
   10 format(1x,' rin,ds(m) =',2e12.3)
!
!-----grid lengths in x and y directions are unity.
!
      xcntr = jx/2.
      ycntr = iy/2.
      dy    = deltas
      dx    = deltas
      ris   = rin**2
!-----rin is radius of influence in grid units
      ris2  = 2.*ris
!
      nskip = 1
      do 30 j = 1,jx
      do 30 i = 1,iy
      cor(i,j) = 0.0
      sum(i,j) = 0.0
      ns(i,j)  = 0
      wtmax(i,j) = 0.0
      htsav(i,j) = 0.0
   30 continue
!
!-----begin to process the nobs observations with loop 22
!
      do 80 kk = 1,nobs,nskip
      if (asta(kk) .gt. 400.) go to 80
!
!-----define max and min i and j values to limit the number of points
!-----must be considered.
!
!-----dind obs. location in terms of dx and dy
!
      dxobs = xobs(kk)/dx
      dyobs = yobs(kk)/dy
!-----convert obs. location to grid increment values of i and j
      rjobs = dxobs + xcntr - dxcen
      riobs = dyobs + ycntr - dycen
!
      ymaxi = riobs + rin
      maxi  = int(ymaxi + 0.99)
      maxi  = min0(maxi,ie)
!
      ymini = riobs - rin
      mini  = int(ymini)
      mini  = max0(mini,1)
!
      xmaxj = rjobs + rin
      maxj  = int(xmaxj + 0.99)
      maxj  = min0(maxj,je)
!
      xminj = rjobs - rin
      minj  = int(xminj)
      minj  = max0(minj,1)
!
      do 70 i=mini,maxi
      do 70 j=minj,maxj
         x = j-xcntr + dxcen
         y = i-ycntr + dycen
!-----compute distance of k_th station from i,jth grid point
         rx  = x-xobs(kk)/dx
         ry  = y-yobs(kk)/dy
         rsq = rx**2+ry**2
         if (rsq.ge.ris2) goto 70
         wt   = (ris-rsq)/(ris+rsq)
   60    continue
!
!-----save max. weighting factor and terrain height to check if grid
!-----point should be treated as a land or sea point.
!
         if (wt.gt.0.0) then
            wtmax(i,j) = amax1(wt,wtmax(i,j))
            if ((wt - wtmax(i,j)).eq.0.0) htsav(i,j) = asta(kk)
            cor(i,j)   = cor(i,j) + wt*asta(kk)
            sum(i,j)   = sum(i,j) + wt
            ns(i,j)    = ns(i,j) + 1
         endif
  70  continue
  80  continue
!
!-----now apply summed weights and weighted observations to determine
!-----terrain value at i,j points
!
      do 90 i = 1,ie
      do 90 j = 1,je
         if (ns(i,j) .ne. 0) then
            cor(i,j) = cor(i,j)/sum(i,j)
            a2(i,j)  = cor(i,j)
!--------if closest observation to i,j is ocean, override a2(i,j) to
!--------preserve the coastline.
!Sara
!           if (htsav(i,j) .le. 0.001) a2(i,j) = htsav(i,j)
!Sara_
         endif
   90 continue
   26 format(' no observations are within rin=',f7.2,
     & ' grid lengths of i=',i3,' j=',i3)

!-----may want to smooth final field a2 here
      return
      end

      subroutine interp

      implicit none
      include 'domain.param'
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay
      real*4  xlat(iy,jx),xlon(iy,jx),xmap(iy,jx),coriol(iy,jx)
      real*4  dlat(iy,jx),dlon(iy,jx),dmap(iy,jx)
      real*4  htsdgrid(iy,jx)
      real*4  htgrid(iy,jx),lndout(iy,jx),texout(iy,jx)
      real*4  sanda(iy,jx),sandb(iy,jx),claya(iy,jx),clayb(iy,jx)
      real*4  frac_lnd(iy,jx,nveg),frac_tex(iy,jx,ntex)
      character*1 ch(iy,jx)
      real*4  snowam(iy,jx),mask(iy,jx)
      integer lnduse(iy,jx),intext(iy,jx)
      common /maps/ xlat,xlon,xmap,coriol,dlat,dlon,dmap
     &     , htgrid,htsdgrid,sanda,sandb,claya,clayb
     &     , frac_lnd,frac_tex,snowam,mask,lndout,texout
     &     , lnduse,intext
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
!
      logical flag
      integer i,j,ii,iindex,jindex
      real*4  dsgrid
      real*8  xx,yy
      real*8  bint
      external bint
      real*8  xin1(iter,jter), xin2(iter,jter)
      real*8  h1, h2, v21
!
      do i = 1,iter
      do j = 1,jter
        xin1(i,j) = 0.0
        xin2(i,j) = 0.0
      end do
      end do
      do ii=1,nobs
        jindex = (xobs(ii)-grdlnmn)*nnc + 1.1
        iindex = (yobs(ii)-grdltmn)*nnc + 1.1
        if (iindex.gt.iter .or. jindex.gt.jter) then
          print 410,ii,xobs(ii),nnc,yobs(ii),iindex,jindex
          stop 400
        endif
        h1 = max(ht(ii),0.0)
        xin1(iindex,jindex) = h1/100.
        h2 = max(ht2(ii),0.0)
        xin2(iindex,jindex) = h2/100000.
      enddo

      flag  = .false.
      dsgrid = float(ntypec)/60.
 
      do i=1,iy-1
      do j=1,jx-1
 
         yy = -(grdltmn-xlat(i,j))/dsgrid + 1.0
         if (grdlnmn.le.-180.0.and.xlon(i,j).gt.0.0)
     1       xlon(i,j)=xlon(i,j)-360.
         xx = -(grdlnmn-xlon(i,j))/dsgrid + 1.0
 
!        yy and xx are the exact index values of a point i,j of the
!        mesoscale mesh when projected onto an earth-grid of lat_s
!        and lon_s for which terrain observations are available.  it
!        is assumed that the earth grid has equal spacing in both
!        latitude and longitude.
 
         h1 = max(bint(yy,xx,xin1,iter,jter,flag),0.d0)*100.
         h2 = max(bint(yy,xx,xin2,iter,jter,flag),0.d0)*100000.
         htgrid(i,j) = h1
         v21 = h2 - h1**2
         htsdgrid(i,j) = sqrt(max(v21,0.d0))

      end do
      end do

  410 format(1x,'ii = ',i6,' xobs(ii) = ',f10.4,' incr = ',i3,
     1  'yobs(ii) = ',f10.4,' iindex = ',i10,'  jindex = ',i10)
 
      return
      end
      subroutine interp_s

      implicit none
      include 'domain.param'
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay
      real*4  xlat(iy*nsg,jx*nsg),xlon(iy*nsg,jx*nsg)
      real*4  dlat(iy*nsg,jx*nsg),dlon(iy*nsg,jx*nsg)
      real*4  xmap(iy*nsg,jx*nsg),dmap(iy*nsg,jx*nsg)
      real*4  coriol(iy*nsg,jx*nsg)
      real*4  htsdgrid(iy*nsg,jx*nsg),htgrid(iy*nsg,jx*nsg)
      real*4  lndout(iy*nsg,jx*nsg),texout(iy*nsg,jx*nsg)
      real*4  sanda(iy*nsg,jx*nsg),sandb(iy*nsg,jx*nsg)
      real*4  claya(iy*nsg,jx*nsg),clayb(iy*nsg,jx*nsg)
      real*4  frac_lnd(iy*nsg,jx*nsg,nveg),frac_tex(iy*nsg,jx*nsg,ntex)
      character*1 ch(iy*nsg,jx*nsg)
      real*4  snowam(iy*nsg,jx*nsg),mask(iy*nsg,jx*nsg)
      integer lnduse(iy*nsg,jx*nsg),intext(iy*nsg,jx*nsg)
      common /maps_s/ xlat,xlon,xmap,coriol,dlat,dlon,dmap
     &     , htgrid,htsdgrid,sanda,sandb,claya,clayb
     &     , frac_lnd,frac_tex,snowam,mask,lndout,texout
     &     , lnduse,intext
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
!
      logical flag
      integer i,j,ii,iindex,jindex
      real*4  dsgrid
      real*8  xx,yy
      real*8  bint
      external bint
      real*8  xin1(iter,jter), xin2(iter,jter)
      real*8  h1, h2, v21
!
      do i = 1,iter
      do j = 1,jter
        xin1(i,j) = 0.0
        xin2(i,j) = 0.0
      end do
      end do
      do ii=1,nobs
        jindex = (xobs(ii)-grdlnmn)*nnc + 1.1
        iindex = (yobs(ii)-grdltmn)*nnc + 1.1
        if (iindex.gt.iter .or. jindex.gt.jter) then
          print 410,ii,xobs(ii),nnc,yobs(ii),iindex,jindex
          stop 400
        endif
        h1 = max(ht(ii),0.0)
        xin1(iindex,jindex) = h1/100.
        h2 = max(ht2(ii),0.0)
        xin2(iindex,jindex) = h2/100000.
      enddo

      flag  = .false.
      dsgrid = float(ntypec)/60.
 
      do i=1,iy*nsg-1
      do j=1,jx*nsg-1
 
         yy = -(grdltmn-xlat(i,j))/dsgrid + 1.0
         if (grdlnmn.le.-180.0.and.xlon(i,j).gt.0.0)
     1       xlon(i,j)=xlon(i,j)-360.
         xx = -(grdlnmn-xlon(i,j))/dsgrid + 1.0
 
!        yy and xx are the exact index values of a point i,j of the
!        mesoscale mesh when projected onto an earth-grid of lat_s
!        and lon_s for which terrain observations are available.  it
!        is assumed that the earth grid has equal spacing in both
!        latitude and longitude.
 
         h1 = max(bint(yy,xx,xin1,iter,jter,flag),0.d0)*100.
         h2 = max(bint(yy,xx,xin2,iter,jter,flag),0.d0)*100000.
         htgrid(i,j) = h1
         v21 = h2 - h1**2
         htsdgrid(i,j) = sqrt(max(v21,0.d0))

      end do
      end do

  410 format(1x,'ii = ',i6,' xobs(ii) = ',f10.4,' incr = ',i3,
     1  'yobs(ii) = ',f10.4,' iindex = ',i10,'  jindex = ',i10)
 
      return
      end
 
      function bint(xx,yy,list,iii,jjj,flag)
      implicit none
!
!-----bilinear interpolation among four grid values
!
      integer iii,jjj
      real*8  xx,yy
      real*8  list(iii,jjj),stl(4,4)
      logical flag
      integer i,j,k,kk,l,n,ll,knear,lnear
      real*8  x,y,a,b,c,d,e,f,g,h
      real*8  bint
      real*8  oned
      external oned
!
      bint = 0.0
      n    = 0
      i    = int(xx+0.00001)
      j    = int(yy+0.00001)
      x    = xx-i
      y    = yy-j
      if (abs(x).gt.0.0001.or.abs(y).gt.0.0001) goto 10
      bint = list(i,j)
      return
   10 continue
      do 20 k=1,4
         kk = i+k-2
         do 20 l=1,4
            stl(k,l) = 0.
            if (flag.and.(l.eq.1)) goto 20
            if (flag.and.(l.eq.4)) goto 20
            if (flag.and.(k.eq.1)) goto 20
            if (flag.and.(k.eq.4)) goto 20
            ll = j+l-2
            if (kk.lt.1.or.kk.gt.iii) goto 20
            if (ll.gt.jjj.or.ll.lt.1) goto 20
            stl(k,l) = list(kk,ll)
            n = n+1
            if (stl(k,l).le.0.0) stl(k,l) = -1.e-20
   20 continue
!
!-----find index of closest point to xx,yy.
!
      knear = float(2) + x + 0.5
      lnear = float(2) + y + 0.5
      a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
      b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
      c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
      d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
      bint = oned(y,a,b,c,d)
!
!--------if closest point is ocean, automatically reset terrain to
!--------preserve coastline.
!
      if (.not.flag.and.stl(knear,lnear).le.0.001) bint = -0.00001
      if (n.eq.16) return
      if (flag.and.n.eq.4) return
      e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
      f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
      g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
      h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
      bint = (bint+oned(x,e,f,g,h))/2.
      if (.not.flag.and.stl(knear,lnear).le.0.001) bint = -0.00001
      return
      end

      subroutine foricbc(iy,jx,kz,nsg,IDATE1,IDATE2,igrads,ibigend
     &            ,ibyte,DATTYP,SSTTYP,EHSO4,LSMTYP,AERTYP,filout)
      implicit none
      integer iy,jx,kz,nsg,IDATE1,IDATE2,igrads,ibigend,ibyte
      CHARACTER*5 DATTYP, SSTTYP
      LOGICAL     EHSO4
      CHARACTER*4 LSMTYP
      CHARACTER*7 AERTYP
      integer JX_O,IY_O,KZ_O
      character a80*80, filout*(*)
      integer isystm,system
!      external system
!
      open(23,file='../ICBC/icbc.param')
      write(23,'(a)') '      INTEGER JX'
      write(23,'(a)') '      INTEGER IY'
      write(23,'(a)') '      INTEGER KZ'
      write(23,'(a)') '      INTEGER NSG'
      write(23,'(a)') '      INTEGER JX_O'
      write(23,'(a)') '      INTEGER IY_O'
      write(23,'(a)') '      INTEGER KZ_O'
      write(23,'(a)') '      INTEGER NP'
      write(23,'(a)') '      INTEGER IDATE1'
      write(23,'(a)') '      INTEGER IDATE2'
      write(23,'(a)') '      INTEGER IBYTE'
      write(23,'(a)') '      CHARACTER*5 DATTYP'
      write(23,'(a)') '      CHARACTER*5 SSTTYP'
      write(23,'(a)') '      LOGICAL     EHSO4 '
      write(23,'(a)') '      CHARACTER*4 LSMTYP'
      write(23,'(a)') '      CHARACTER*7 AERTYP'
      write(23,100) 'JX     =',JX
      write(23,100) 'IY     =',IY
      write(23,100) 'KZ     =',kz
      write(23,100) 'NSG    =',nsg
      IF(DATTYP.EQ.'FNEST') THEN
         WRITE(*,*) 'You are preparing the NEST run from previous run'
         WRITE(*,*) 'Please input the original size of ARRAYs'
         WRITE(*,*) 'JX_O, IY_O, KZ_O'
         WRITE(*,*) 'You need cut 2 rows for both JX_O and IY_O'
         READ(*,*)   JX_O, IY_O, KZ_O
         write(23,100) 'JX_O   =',JX_O
         write(23,100) 'IY_O   =',IY_O
         write(23,100) 'KZ_O   =',KZ_O
         write(23,100) 'NP     =',15
      ELSE
         write(23,100) 'JX_O   =',1
         write(23,100) 'IY_O   =',1
         write(23,100) 'KZ_O   =',14
         write(23,100) 'NP     =',15
      ENDIF
      write(23,200) 'IDATE1 =',IDATE1
      write(23,200) 'IDATE2 =',IDATE2
      write(23,200) 'IBYTE  =',IBYTE
 100  format('      parameter(',A8,I4,')')
 200  format('      parameter(',A8,I10,')')
      if(DATTYP.EQ.'ECMWF') then
         write(23,'(a)') "      parameter(DATTYP='ECMWF')"
      else if(DATTYP.EQ.'ERA40') then
         write(23,'(a)') "      parameter(DATTYP='ERA40')"
      else if(DATTYP.EQ.'ERAIN'.or.DATTYP.EQ.'EIN15') then
         write(23,'(a)') "      parameter(DATTYP='EIN15')"
      else if(DATTYP.EQ.'EIN75') then
         write(23,'(a)') "      parameter(DATTYP='EIN75')"
      else if(DATTYP.EQ.'EIN25') then
         write(23,'(a)') "      parameter(DATTYP='EIN25')"
      else if(DATTYP.EQ.'ERAHI') then
         write(23,'(a)') "      parameter(DATTYP='ERAHI')"
      else if(DATTYP.EQ.'NNRP1') then
         write(23,'(a)') "      parameter(DATTYP='NNRP1')"
      else if(DATTYP.EQ.'NNRP2') then
         write(23,'(a)') "      parameter(DATTYP='NNRP2')"
      else if(DATTYP.EQ.'NRP2W') then
         write(23,'(a)') "      parameter(DATTYP='NRP2W')"
      else if(DATTYP.EQ.'GFS11') then
         write(23,'(a)') "      parameter(DATTYP='GFS11')"
      else if(DATTYP.EQ.'EH5OM') then
         write(23,'(a)') "      parameter(DATTYP='EH5OM')"
      else if(DATTYP.EQ.'FVGCM') then
         write(23,'(a)') "      parameter(DATTYP='FVGCM')"
      else if(DATTYP.EQ.'FNEST') then
         write(23,'(a)') "      parameter(DATTYP='FNEST')"
      else
         print*,'BOUNDARY CONDITION DATASET DOES NOT EXIST'
         stop 'subroutine foricbc'
      endif
      if(SSTTYP.EQ.'OISST'.or.SSTTYP.EQ.'OI_NC') then
         write(23,'(a)') "      parameter(SSTTYP='OISST')"
      else if(SSTTYP.EQ.'OI_WK') then
         write(23,'(a)') "      parameter(SSTTYP='OI_WK')"
      else if(SSTTYP.EQ.'GISST') then
         write(23,'(a)') "      parameter(SSTTYP='GISST')"
      else if(SSTTYP.EQ.'FV_RF') then
         write(23,'(a)') "      parameter(SSTTYP='FV_RF')"
      else if(SSTTYP.EQ.'FV_A2') then
         write(23,'(a)') "      parameter(SSTTYP='FV_A2')"
      else if(SSTTYP.EQ.'FV_B2') then
         write(23,'(a)') "      parameter(SSTTYP='FV_B2')"
      else if(SSTTYP.EQ.'EH5RF') then
         write(23,'(a)') "      parameter(SSTTYP='EH5RF')"
      else if(SSTTYP.EQ.'EH5A2') then
         write(23,'(a)') "      parameter(SSTTYP='EH5A2')"
      else if(SSTTYP.EQ.'EH5B1') then
         write(23,'(a)') "      parameter(SSTTYP='EH5B1')"
      else if(SSTTYP.EQ.'EHA1B') then
         write(23,'(a)') "      parameter(SSTTYP='EHA1B')"
      else if(SSTTYP.EQ.'ERSST') then
         write(23,'(a)') "      parameter(SSTTYP='ERSST')"
      else if(SSTTYP.EQ.'ERSKT') then
         write(23,'(a)') "      parameter(SSTTYP='ERSKT')"
      else
         print*,'SST DATASET DOES NOT EXIST'
         stop 'subroutine foricbc'
      endif
      if(DATTYP.EQ.'EH5OM'.and.EHSO4) then
         write(23,'(a)') "      parameter(EHSO4 =.true. )"
      else
         write(23,'(a)') "      parameter(EHSO4 =.false.)"
      endif
      if(LSMTYP.EQ.'BATS') then
         write(23,'(a)') "      parameter(LSMTYP='BATS')"
      else if(LSMTYP.EQ.'USGS') then
         write(23,'(a)') "      parameter(LSMTYP='USGS')"
      else
         print*,'LANDUSE LEGEND DOES NOT EXIST'
         stop 'subroutine foricbc'
      endif
      if(AERTYP.EQ.'AER00D0') then
         write(23,'(a)') "      parameter(AERTYP='AER00D0')"
      else if(AERTYP.EQ.'AER01D0') then
         write(23,'(a)') "      parameter(AERTYP='AER01D0')"
      else if(AERTYP.EQ.'AER10D0') then
         write(23,'(a)') "      parameter(AERTYP='AER10D0')"
      else if(AERTYP.EQ.'AER11D0') then
         write(23,'(a)') "      parameter(AERTYP='AER11D0')"
      else if(AERTYP.EQ.'AER00D1') then
         write(23,'(a)') "      parameter(AERTYP='AER00D1')"
      else if(AERTYP.EQ.'AER01D1') then
         write(23,'(a)') "      parameter(AERTYP='AER01D1')"
      else if(AERTYP.EQ.'AER10D1') then
         write(23,'(a)') "      parameter(AERTYP='AER10D1')"
      else if(AERTYP.EQ.'AER11D1') then
         write(23,'(a)') "      parameter(AERTYP='AER11D1')"
      else
         print*,'AEROSOL TYPE DOES NOT EXIST'
         stop 'subroutine foricbc'
      endif
      close(23)
      open(23,file='../ICBC/icbc.x')
      write(23,'(a)') '#!/bin/csh -f'
      write(23,'(a)') 'make clean'
      write(23,'(a)') 'foreach FILE (RCM_SST.ctl RCM_SST.dat SST.RCM)'
      write(23,'(a)') 'if ( -f $FILE ) /bin/rm $FILE'
      write(23,'(a)') 'end'
      IF(DATTYP.EQ.'FVGCM'.or.
     &        (DATTYP.EQ.'FNEST'.and.
     &         (SSTTYP.EQ.'FV_RF'.or.SSTTYP.EQ.'FV_A2'))) THEN
         write(23,'(a)') 'make SST_FVGCM'
         write(23,'(a)') './SST_FVGCM'
         write(23,'(a)') '/bin/rm -f SST_FVGCM*.o SST_FVGCM'
      ELSE IF( DATTYP.EQ.'EH5OM'.or.
     &        (DATTYP.EQ.'FNEST'.and.
     &         (SSTTYP.EQ.'EH5RF'.or.SSTTYP.EQ.'EH5A2'.or.
     &          SSTTYP.EQ.'EH5B1'.or.SSTTYP.EQ.'EHA1B'))) THEN
         write(23,'(a)') 'make SST_EH5OM'
         write(23,'(a)') './SST_EH5OM'
         write(23,'(a)') '/bin/rm -f SST_EH5OM*.o SST_EH5OM'
      ELSE IF((DATTYP.EQ.'EIN15'.or.DATTYP.EQ.'ERAIN').and.
     &        (SSTTYP.EQ.'ERSST'.or.SSTTYP.EQ.'ERSKT')) THEN
         write(23,'(a)') 'make SST_ERSST'
         write(23,'(a)') './SST_ERSST'
         write(23,'(a)') '/bin/rm -f SST_ERSST*.o SST_ERSST'
      ELSE
         write(23,'(a)') 'make SST'
         write(23,'(a)') './SST'
         write(23,'(a)') '/bin/rm -f SST_1DEG*.o SST'
      ENDIF
      IF(.not.(AERTYP(4:5).eq.'00')) THEN
         write(23,'(a)') 'make AEROSOL'
         write(23,'(a)') './AEROSOL'
         write(23,'(a)') '/bin/rm -f AEROSOL.o AEROSOL'
      ENDIF
      IF(DATTYP.EQ.'ERAHI') THEN
         write(23,'(a)') 'cp ../../Commons/tools/srcERAHI/*.f .'
         write(23,'(a)') 'make ERAHI_HT'
         write(23,'(a)') './ERAHI_HT'
         write(23,'(a)') 'make ERAHI_PS'
         write(23,'(a)') './ERAHI_PS'
         write(23,'(a)') 'make ERAHI_T'
         write(23,'(a)') './ERAHI_T'
         write(23,'(a)') 'make ERAHI_Q'
         write(23,'(a)') './ERAHI_Q'
         write(23,'(a)') 'make ERAHI_U'
         write(23,'(a)') './ERAHI_U'
         write(23,'(a)') 'make ERAHI_V'
         write(23,'(a)') './ERAHI_V'
      ENDIF
      write(23,'(a)') 'make ICBC'
      write(23,'(a)') './ICBC'
      write(23,'(a)') '/bin/rm -f ICBC.o ICBC SST.RCM'
      close(23)
      IF(DATTYP.EQ.'ERAHI') THEN
         open(23,file='../ICBC/ERAHI.param')
         write(23,'(a)') '      INTEGER IDATE1'
         write(23,'(a)') '      INTEGER IDATE2'
         write(23,'(a)') '      INTEGER IBYTE'
         write(23,200) 'IDATE1 =',IDATE1
         write(23,200) 'IDATE2 =',IDATE2
         write(23,200) 'IBYTE  =',IBYTE
         close(23)
      ENDIF

      return
      end
      subroutine formain(iy,jx,kz,nsg,ibyte,DATTYP,EHSO4,LSMTYP,AERTYP
     &                  ,NPROC)
      implicit none
      integer iy,jx,kz,nsg,ibyte,NPROC
      CHARACTER*5 DATTYP
      LOGICAL     EHSO4
      CHARACTER*4 LSMTYP
      CHARACTER*7 AERTYP
!
      open(23,file='../../Main/regcm.param')
      write(23,'(a)') '      INTEGER IX'
      write(23,'(a)') '      INTEGER JX'
      write(23,'(a)') '      INTEGER KX'
      write(23,'(a)') '      INTEGER NSG'
      write(23,'(a)') '      INTEGER NNSG'
      write(23,'(a)') '      INTEGER IBYTE'
      write(23,'(a)') '      CHARACTER*5 DATTYP'
      write(23,'(a)') '      LOGICAL     EHSO4 '
      write(23,'(a)') '      CHARACTER*4 LSMTYP'
      write(23,'(a)') '      CHARACTER*7 AERTYP'
      write(23,'(a)') '      integer jlx,jlxm'
      write(23,100) 'IX     =',IY
      write(23,100) 'JX     =',JX
      write(23,100) 'KX     =',KZ
      write(23,100) 'NSG    =',NSG
      write(23,100) 'NNSG   =',NSG*NSG
      write(23,100) 'IBYTE  =',IBYTE
      write(23,200) "DATTYP='",DATTYP
      if(DATTYP.EQ.'EH5OM'.and.EHSO4) then
         write(23,'(a)') "      parameter(EHSO4 =.true. )"
      else
         write(23,'(a)') "      parameter(EHSO4 =.false.)"
      endif
      write(23,300) "LSMTYP='",LSMTYP
      write(23,400) "AERTYP='",AERTYP
      write(23,'(a)') '      parameter(jlx=jx-1,jlxm=jx-2)'
      close(23)
 100  format('      parameter(',A8,I6,')')
 200  format('      parameter(',A8,A5,"')")
 300  format('      parameter(',A8,A4,"')")
 400  format('      parameter(',A8,A7,"')")
 500  format('      parameter(',A18,")")
!
      if(NPROC.gt.0) then
        open(23,file='../../Main/regcm.param2')
        write(23,'(a)') '      INTEGER IX'
        write(23,'(a)') '      INTEGER NPROC'
        write(23,'(a)') '      INTEGER MJX'
        write(23,'(a)') '      INTEGER KX'
        write(23,'(a)') '      INTEGER NSG'
        write(23,'(a)') '      INTEGER NNSG'
        write(23,'(a)') '      INTEGER IBYTE'
        write(23,'(a)') '      INTEGER JXP'
        write(23,'(a)') '      CHARACTER*5 DATTYP'
        write(23,'(a)') '      LOGICAL     EHSO4 '
        write(23,'(a)') '      CHARACTER*4 LSMTYP'
        write(23,'(a)') '      CHARACTER*7 AERTYP'
        write(23,'(a)') '      integer jxbb'
        write(23,100) 'IX     =',IY
        write(23,100) 'NPROC  =',NPROC
        write(23,100) 'MJX    =',JX
        if(mod(JX,NPROC).eq.0) then
          write(23,500) 'JXP    = MJX/NPROC'
        else
          write(*,*) 'The present parallel code requirea that'
          write(*,*) 'MJX can be divided by NPROC'
          stop 'subroutine formain' 
        endif
        write(23,100) 'KX     =',KZ
        write(23,100) 'NSG    =',NSG
        write(23,100) 'NNSG   =',NSG*NSG
        write(23,100) 'IBYTE  =',IBYTE
        write(23,200) "DATTYP='",DATTYP
        if(DATTYP.EQ.'EH5OM'.and.EHSO4) then
           write(23,'(a)') "      parameter(EHSO4 =.true. )"
        else
           write(23,'(a)') "      parameter(EHSO4 =.false.)"
        endif
        write(23,300) "LSMTYP='",LSMTYP
        write(23,400) "AERTYP='",AERTYP
        write(23,'(a)') '      parameter(jxbb=mjx-1)'
        close(23)
      endif
!
      return
      end

      SUBROUTINE LNDFUDGE(FUDGE,ch,lndout,htgrid,iy,jx,LSMTYP,CHAR_LND)
      implicit none
      logical FUDGE
      integer iy,jx
      real*4  htgrid(iy,jx)
      real*4  lndout(iy,jx)
      character*1 ch(iy,jx)
      CHARACTER*4 LSMTYP
      CHARACTER*10 CHAR_LND
      integer i,j
!
      if(LSMTYP.eq.'BATS') then
         if(FUDGE) then
            open(13,file=CHAR_LND,form='formatted')
            do i=iy,1,-1
               read(13,100)(ch(i,j),j=1,jx)
            enddo
            close(13)
            do j=1,jx
            do i=1,iy
               if(ch(i,j).eq.' ') then
                  lndout(i,j) = 15.
               else if(ch(i,j).eq.'1') then
                  lndout(i,j) = 1.
               else if(ch(i,j).eq.'2') then
                  lndout(i,j) = 2.
               else if(ch(i,j).eq.'3') then
                  lndout(i,j) = 3.
               else if(ch(i,j).eq.'4') then
                  lndout(i,j) = 4.
               else if(ch(i,j).eq.'5') then
                  lndout(i,j) = 5.
               else if(ch(i,j).eq.'6') then
                  lndout(i,j) = 6.
               else if(ch(i,j).eq.'7') then
                  lndout(i,j) = 7.
               else if(ch(i,j).eq.'8') then
                  lndout(i,j) = 8.
               else if(ch(i,j).eq.'9') then
                  lndout(i,j) = 9.
               else if(ch(i,j).eq.'A') then
                  lndout(i,j) = 10.
               else if(ch(i,j).eq.'B') then
                  lndout(i,j) = 11.
               else if(ch(i,j).eq.'C') then
                  lndout(i,j) = 12.
               else if(ch(i,j).eq.'D') then
                  lndout(i,j) = 13.
               else if(ch(i,j).eq.'E') then
                  lndout(i,j) = 14.
               else if(ch(i,j).eq.'F') then
                  lndout(i,j) = 15.
               else if(ch(i,j).eq.'G') then
                  lndout(i,j) = 16.
               else if(ch(i,j).eq.'H') then
                  lndout(i,j) = 17.
               else if(ch(i,j).eq.'I') then
                  lndout(i,j) = 18.
               else if(ch(i,j).eq.'J') then
                  lndout(i,j) = 19.
               else if(ch(i,j).eq.'K') then
                  lndout(i,j) = 20.
               else if(nint(lndout(i,j)).eq.0) then
!                 ch(i,j) = 'X'
                  ch(i,j) = ' '
               else
                  write(*,*) 'LANDUSE MASK exceed the limit'
                  STOP
               endif
!_fix          if(nint(lndout(i,j)).eq.15) htgrid(i,j) = 0.0
               if(htgrid(i,j).lt.0.1.and.nint(lndout(i,j)).eq.15)
     &            htgrid(i,j) = 0.0
            enddo
            enddo
         else
            do j=1,jx
            do i=1,iy
               if(nint(lndout(i,j)).eq.15) then
                  ch(i,j) = ' '
               else if(nint(lndout(i,j)).eq.1) then
                  ch(i,j) = '1'
               else if(nint(lndout(i,j)).eq.2) then
                  ch(i,j) = '2'
               else if(nint(lndout(i,j)).eq.3) then
                  ch(i,j) = '3'
               else if(nint(lndout(i,j)).eq.4) then
                  ch(i,j) = '4'
               else if(nint(lndout(i,j)).eq.5) then
                  ch(i,j) = '5'
               else if(nint(lndout(i,j)).eq.6) then
                  ch(i,j) = '6'
               else if(nint(lndout(i,j)).eq.7) then
                  ch(i,j) = '7'
               else if(nint(lndout(i,j)).eq.8) then
                  ch(i,j) = '8'
               else if(nint(lndout(i,j)).eq.9) then
                  ch(i,j) = '9'
               else if(nint(lndout(i,j)).eq.10) then
                  ch(i,j) = 'A'
               else if(nint(lndout(i,j)).eq.11) then
                  ch(i,j) = 'B'
               else if(nint(lndout(i,j)).eq.12) then
                  ch(i,j) = 'C'
               else if(nint(lndout(i,j)).eq.13) then
                  ch(i,j) = 'D'
               else if(nint(lndout(i,j)).eq.14) then
                  ch(i,j) = 'E'
               else if(nint(lndout(i,j)).eq.16) then
                  ch(i,j) = 'G'
               else if(nint(lndout(i,j)).eq.17) then
                  ch(i,j) = 'H'
               else if(nint(lndout(i,j)).eq.18) then
                  ch(i,j) = 'I'
               else if(nint(lndout(i,j)).eq.19) then
                  ch(i,j) = 'J'
               else if(nint(lndout(i,j)).eq.20) then
                  ch(i,j) = 'K'
               else
                  write(*,*) 'LANDUSE MASK',nint(lndout(i,j)),
     &                       'exceed the limit'
                  STOP
               endif
            enddo
            enddo
            open(13,file=CHAR_LND,form='formatted')
            do i=iy,1,-1
               write(13,100)(ch(i,j),j=1,jx)
            enddo
            close(13)
         endif
      else if(LSMTYP.eq.'USGS') then
         if(FUDGE) then
            open(13,file=CHAR_LND,form='formatted')
            do i=iy,1,-1
               read(13,100)(ch(i,j),j=1,jx)
            enddo
            close(13)
            do j=1,jx
            do i=1,iy
               if(ch(i,j).eq.' ') then
                  lndout(i,j) = 25.
               else if(ch(i,j).eq.'1') then
                  lndout(i,j) = 1.
               else if(ch(i,j).eq.'2') then
                  lndout(i,j) = 2.
               else if(ch(i,j).eq.'3') then
                  lndout(i,j) = 3.
               else if(ch(i,j).eq.'4') then
                  lndout(i,j) = 4.
               else if(ch(i,j).eq.'5') then
                  lndout(i,j) = 5.
               else if(ch(i,j).eq.'6') then
                  lndout(i,j) = 6.
               else if(ch(i,j).eq.'7') then
                  lndout(i,j) = 7.
               else if(ch(i,j).eq.'8') then
                  lndout(i,j) = 8.
               else if(ch(i,j).eq.'9') then
                  lndout(i,j) = 9.
               else if(ch(i,j).eq.'A') then
                  lndout(i,j) = 10.
               else if(ch(i,j).eq.'B') then
                  lndout(i,j) = 11.
               else if(ch(i,j).eq.'C') then
                  lndout(i,j) = 12.
               else if(ch(i,j).eq.'D') then
                  lndout(i,j) = 13.
               else if(ch(i,j).eq.'E') then
                  lndout(i,j) = 14.
               else if(ch(i,j).eq.'F') then
                  lndout(i,j) = 15.
               else if(ch(i,j).eq.'G') then
                  lndout(i,j) = 16.
               else if(ch(i,j).eq.'H') then
                  lndout(i,j) = 17.
               else if(ch(i,j).eq.'I') then
                  lndout(i,j) = 18.
               else if(ch(i,j).eq.'J') then
                  lndout(i,j) = 19.
               else if(ch(i,j).eq.'K') then
                  lndout(i,j) = 20.
               else if(ch(i,j).eq.'L') then
                  lndout(i,j) = 21.
               else if(ch(i,j).eq.'M') then
                  lndout(i,j) = 22.
               else if(ch(i,j).eq.'N') then
                  lndout(i,j) = 23.
               else if(ch(i,j).eq.'O') then
                  lndout(i,j) = 24.
               else
                  write(*,*) 'LANDUSE MASK exceed the limit'
                  STOP
               endif
               if(htgrid(i,j).lt.0.1.and.
     &            (nint(lndout(i,j)).eq.25.or.nint(lndout(i,j)).eq.0))
     &            htgrid(i,j) = 0.0
            enddo
            enddo
         else
            do j=1,jx
            do i=1,iy
               if(nint(lndout(i,j)).eq.25.or.
     &            nint(lndout(i,j)).eq.0) then
                  ch(i,j) = ' '
               else if(nint(lndout(i,j)).eq.1) then
                  ch(i,j) = '1'
               else if(nint(lndout(i,j)).eq.2) then
                  ch(i,j) = '2'
               else if(nint(lndout(i,j)).eq.3) then
                  ch(i,j) = '3'
               else if(nint(lndout(i,j)).eq.4) then
                  ch(i,j) = '4'
               else if(nint(lndout(i,j)).eq.5) then
                  ch(i,j) = '5'
               else if(nint(lndout(i,j)).eq.6) then
                  ch(i,j) = '6'
               else if(nint(lndout(i,j)).eq.7) then
                  ch(i,j) = '7'
               else if(nint(lndout(i,j)).eq.8) then
                  ch(i,j) = '8'
               else if(nint(lndout(i,j)).eq.9) then
                  ch(i,j) = '9'
               else if(nint(lndout(i,j)).eq.10) then
                  ch(i,j) = 'A'
               else if(nint(lndout(i,j)).eq.11) then
                  ch(i,j) = 'B'
               else if(nint(lndout(i,j)).eq.12) then
                  ch(i,j) = 'C'
               else if(nint(lndout(i,j)).eq.13) then
                  ch(i,j) = 'D'
               else if(nint(lndout(i,j)).eq.14) then
                  ch(i,j) = 'E'
               else if(nint(lndout(i,j)).eq.15) then
                  ch(i,j) = 'F'
               else if(nint(lndout(i,j)).eq.16) then
                  ch(i,j) = 'G'
               else if(nint(lndout(i,j)).eq.17) then
                  ch(i,j) = 'H'
               else if(nint(lndout(i,j)).eq.18) then
                  ch(i,j) = 'I'
               else if(nint(lndout(i,j)).eq.19) then
                  ch(i,j) = 'J'
               else if(nint(lndout(i,j)).eq.20) then
                  ch(i,j) = 'K'
               else if(nint(lndout(i,j)).eq.21) then
                  ch(i,j) = 'L'
               else if(nint(lndout(i,j)).eq.22) then
                  ch(i,j) = 'M'
               else if(nint(lndout(i,j)).eq.23) then
                  ch(i,j) = 'N'
               else if(nint(lndout(i,j)).eq.24) then
                  ch(i,j) = 'O'
               else
                  write(*,*) 'LANDUSE MASK',nint(lndout(i,j)),
     &                       'exceed the limit'
                  STOP
               endif
            enddo
            enddo
            open(13,file=CHAR_LND,form='formatted')
            do i=iy,1,-1
               write(13,100)(ch(i,j),j=1,jx)
            enddo
            close(13)
         endif
      else
         print*,'LANDUSE LEGEND DOES NOT EXIST'
         stop 'subroutine LNDFUDGE'
      endif
 100  format(132A1)
      return
      end

      SUBROUTINE TEXFUDGE(FUDGE,ch,texout,htgrid,iy,jx,CHAR_TEX)
      implicit none
      logical FUDGE
      integer iy,jx
      real*4  htgrid(iy,jx)
      real*4  texout(iy,jx)
      character*1 ch(iy,jx)
      character*10 CHAR_TEX
      integer i,j
!
      if(FUDGE) then
         open(13,file=CHAR_TEX,form='formatted')
         do i=iy,1,-1
            read(13,100)(ch(i,j),j=1,jx)
         enddo
         close(13)
         do j=1,jx
         do i=1,iy
            if(ch(i,j).eq.' ') then
               texout(i,j) = 14.
            else if(ch(i,j).eq.'1') then
               texout(i,j) = 1.
            else if(ch(i,j).eq.'2') then
               texout(i,j) = 2.
            else if(ch(i,j).eq.'3') then
               texout(i,j) = 3.
            else if(ch(i,j).eq.'4') then
               texout(i,j) = 4.
            else if(ch(i,j).eq.'5') then
               texout(i,j) = 5.
            else if(ch(i,j).eq.'6') then
               texout(i,j) = 6.
            else if(ch(i,j).eq.'7') then
               texout(i,j) = 7.
            else if(ch(i,j).eq.'8') then
               texout(i,j) = 8.
            else if(ch(i,j).eq.'9') then
               texout(i,j) = 9.
            else if(ch(i,j).eq.'A') then
               texout(i,j) = 10.
            else if(ch(i,j).eq.'B') then
               texout(i,j) = 11.
            else if(ch(i,j).eq.'C') then
               texout(i,j) = 12.
            else if(ch(i,j).eq.'D') then
               texout(i,j) = 13.
            else if(ch(i,j).eq.'E') then
               texout(i,j) = 14.
            else if(ch(i,j).eq.'F') then
               texout(i,j) = 15.
            else if(ch(i,j).eq.'G') then
               texout(i,j) = 16.
            else if(ch(i,j).eq.'H') then
               texout(i,j) = 17.
            else if(nint(texout(i,j)).eq.0) then
!              ch(i,j) = 'X'
               ch(i,j) = ' '
            else
               write(*,*) 'TEXTURE TYPE exceed the limit'
               STOP
            endif
            if(nint(texout(i,j)).eq.14) htgrid(i,j) = 0.0
         enddo
         enddo
      else
         do j=1,jx
         do i=1,iy
            if(nint(texout(i,j)).eq.14) then
               ch(i,j) = ' '
            else if(nint(texout(i,j)).eq.1) then
               ch(i,j) = '1'
            else if(nint(texout(i,j)).eq.2) then
               ch(i,j) = '2'
            else if(nint(texout(i,j)).eq.3) then
               ch(i,j) = '3'
            else if(nint(texout(i,j)).eq.4) then
               ch(i,j) = '4'
            else if(nint(texout(i,j)).eq.5) then
               ch(i,j) = '5'
            else if(nint(texout(i,j)).eq.6) then
               ch(i,j) = '6'
            else if(nint(texout(i,j)).eq.7) then
               ch(i,j) = '7'
            else if(nint(texout(i,j)).eq.8) then
               ch(i,j) = '8'
            else if(nint(texout(i,j)).eq.9) then
               ch(i,j) = '9'
            else if(nint(texout(i,j)).eq.10) then
               ch(i,j) = 'A'
            else if(nint(texout(i,j)).eq.11) then
               ch(i,j) = 'B'
            else if(nint(texout(i,j)).eq.12) then
               ch(i,j) = 'C'
            else if(nint(texout(i,j)).eq.13) then
               ch(i,j) = 'D'
            else if(nint(texout(i,j)).eq.15) then
               ch(i,j) = 'F'
            else if(nint(texout(i,j)).eq.16) then
               ch(i,j) = 'G'
            else if(nint(texout(i,j)).eq.17) then
               ch(i,j) = 'H'
            else
               write(*,*) 'TEXTURE TYPE',nint(texout(i,j)),
     &                    'exceed the limit'
               STOP
            endif
         enddo
         enddo
         open(13,file=CHAR_TEX,form='formatted')
         do i=iy,1,-1
            write(13,100)(ch(i,j),j=1,jx)
         enddo
         close(13)
      endif
 100  format(132A1)
      return
      end

      SUBROUTINE LAKEADJ(LSMTYP,lnduse,htgrid,xlat,xlon,imx,jmx)

      implicit none
      CHARACTER*4 LSMTYP
      integer imx, jmx, i, j
      integer lnduse(imx,jmx)
      real*4  htgrid(imx,jmx)
      real*4  xlat(imx,jmx),xlon(imx,jmx)

      real*4  zerie,zhuron,zontar,zsup,zmich,xx,yy
      parameter (zerie=174.,zhuron=177.,zontar=75.,zsup=183.,zmich=177.)

! ****  ADJUST GREAT LAKE ELEVATION **** C

      do i=1,imx
      do j=1,jmx
        if ((LSMTYP.eq.'BATS'.and.lnduse(i,j).eq.14).or.
     &      (LSMTYP.eq.'USGS'.and.lnduse(i,j).eq.16)) then
          xx=xlon(i,j)
          yy=xlat(i,j)
          IF (yy.LE.43.2 .AND. yy.GE.41.0 .AND.         ! LAKE ERIE
     &        xx.LE.-78.0 .AND. xx.GE.-84.0) THEN
            print*,'**** ADUJUSTING LAKE ERIE LEVEL ****'
            print*,'     NEW:',zerie,'    OLD:',htgrid(i,j),i,j
            htgrid(i,j)=zerie
          elseif (yy.LE.46.4 .AND. yy.GE.43.0 .AND.     ! LAKE HURON
     &            xx.LE.-79.9 .AND. yy.GE.-85.0) THEN
            print*,'**** ADUJUSTING LAKE HURON LEVEL ****'
            print*,'     NEW:',zhuron,'    OLD:',htgrid(i,j),i,j
            htgrid(i,j)=zhuron
          elseif (yy.LE.44.5 .AND. yy.GE.43.2 .AND.     ! LAKE ONTARIO
     &            xx.LE.-75.0 .AND. yy.GE.-79.9) THEN
            print*,'**** ADUJUSTING LAKE ONTARIO LEVEL ****'
            print*,'     NEW:',zontar,'    OLD:',htgrid(i,j),i,j
            htgrid(i,j)=zontar
          elseif (yy.LE.49.4 .AND. yy.GE.46.2 .AND.     ! LAKE SUPERIOR
     &            xx.LE.-84.2 .AND. xx.GE.-93.0) THEN
            print*,'**** ADUJUSTING LAKE SUPERIOR LEVEL ****'
            print*,'     NEW:',zsup,'    OLD:',htgrid(i,j),i,j
            htgrid(i,j)=zsup
          elseif (yy.LE.46.2 .AND. yy.GE.41.0 .AND.     ! LAKE MICHIGAN
     &            xx.LE.-84.8 .AND. xx.GE.-89.0) THEN
            print*,'**** ADUJUSTING LAKE MICHIGAN LEVEL ****'
            print*,'     NEW:',zmich,'    OLD:',htgrid(i,j),i,j
            htgrid(i,j)=zmich
          endif
        endif
      end do
      end do

      return
      end
      SUBROUTINE LAMBRT(XLON,XLAT,SMAP,CORIOL,IY,JX,CLON,CLAT,DS,IDOT
     &                 ,XN,TRUELATL,TRUELATH)
      implicit none
      REAL*4 CLAT,CLON,DS
!           CLAT IS THE CENTRAL LATITUDE OF THE COARSE DOMAIN.
!           CLON IS THE CENTRAL LONGITUDE OF THE COARSE DOMAIN.
      INTEGER IY,JX
!
!  THIS ROUTINE CALCULATES MESO MAP(LAT,LONG,CORIOLIS,MAP SCALE)
!  FOR LAMBERT CONFORMAL PROJECTION
!
!     IY IS THE I DIMENSION FOR THIS DOMAIN.
!     JX IS THE J DIMENSION FOR THIS DOMAIN.
!     IDOT IS ICROSS ( = 1)  OR IDOT ( = 0).
!
!---------------------------------------------------------------------
!
      REAL*4 AA,XN,PSI1,PI,TRUELATL,TRUELATH,TRUELAT1,TRUELAT2
!
!           XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!           PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!           PI IS PI.
!
!---------------------------------------------------------------------
!
      REAL*4  CORIOL(IY,JX),SMAP(IY,JX)
      REAL*4  XLAT(IY,JX),XLON(IY,JX)
      INTEGER IDOT
      REAL*4  POLE,d2r,OMEGA2
      REAL*4  CNTRJ,CNTRI,PSX,CELL,CELL2,R,XCNTR,YCNTR
      REAL*4  X,Y,FLP,FLPP,CELL1,PSIX,SIGN
      INTEGER I,J

      AA=6.371229E6
      IF(CLAT.LT.0.) THEN
         SIGN=-1.        ! SOUTH HEMESPHERE
      ELSE
         SIGN= 1.        ! NORTH HEMESPHERE
      ENDIF
      POLE=SIGN*90.0
      PI=ATAN(1.)*4.
      d2r=PI/180.

      TRUELAT1=TRUELATH
      TRUELAT2=TRUELATL
      IF (ABS(TRUELAT1-TRUELAT2) .GT. 1.E-1) THEN
        XN=(ALOG10(COS(TRUELAT1*d2r))-ALOG10(COS(TRUELAT2*d2r)))
     &    /(ALOG10(TAN((45.0-SIGN*TRUELAT1/2.0)*d2r))-
     &      ALOG10(TAN((45.0-SIGN*TRUELAT2/2.0)*d2r)))
      ELSE
        XN=SIGN*SIN(TRUELAT1*d2r)
      ENDIF
!     XN=0.716

      PSI1=90.-SIGN*TRUELAT1
      IF(CLAT.LT.0.) THEN
         PSI1=-PSI1
      ENDIF
!
      PSI1=PSI1*d2r
      OMEGA2=1.454441E-4
      cntrj = (JX+IDOT)/2.
      cntri = (IY+IDOT)/2.
!
      PSX = (POLE-CLAT)*d2r
      CELL=AA*SIN(PSI1)/XN
      write(*,*) 'PSX,PSI1 = ',PSX,PSI1
      CELL2=(TAN(PSX/2.))/(TAN(PSI1/2.))
      R=CELL*(CELL2)**XN
      XCNTR=0.
      YCNTR=-R
!
      DO 70 J=1,JX
         X=XCNTR+(J-CNTRJ)*DS
         DO 70 I=1,IY
            Y=YCNTR+(I-CNTRI)*DS
            R = SQRT(X*X + Y*Y)
            IF(Y.EQ.0.) THEN
               IF(X.GE.0.) THEN
                  FLP=90.*d2r
               ELSE
                  FLP=-90.*d2r
               ENDIF
            ELSE
               IF (CLAT.LT.0.0) THEN
                  FLP = ATAN2(X,Y)
               ELSE
                  FLP = ATAN2(X,-Y)
               END IF
            END IF
            FLPP = FLP/XN/d2r+CLON
!           IF (FLPP.GT.180.0) FLPP = FLPP-360.0
!           IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
            XLON(I,J) = FLPP
            IF (CLAT.LT.0.0) R = -R
            CELL=R*XN/(AA*SIN(PSI1))
            CELL1=TAN(PSI1/2.)*CELL**(1./XN)
            CELL2=ATAN(CELL1)
            PSX=2.*CELL2/d2r
            XLAT(I,J) = POLE - PSX
            PSIX=PSX*d2r
            SMAP(I,J)=(SIN(PSI1)/SIN(PSIX))*
     &           ((TAN(PSIX/2.)/TAN(PSI1/2.))**XN)
 70   CONTINUE
      IF (IDOT.EQ.1) THEN
         DO 20 I=1,IY
         DO 20 J=1,JX
         CORIOL(I,J)=OMEGA2*SIN(XLAT(I,J)*d2r)
   20    CONTINUE
      ENDIF
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ROTMER( XLON, XLAT, XMAP, CORIOL, IY, JX 
     A                 , CLON,CLAT, POLLON, POLLAT, DS, IDOT)    
      implicit none
!---------------------------------------------------------------------  
!  COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
!  ROTATED POLE MAP CENTERED AT CLON,CLAT. ORIGIN OF ROTATED GRID IS
!  GIVEN BY POLLON AND POLLAT.
!  IMX,JMX,KSV,KSH AND LSIG2 MUST CORRESPOND TO IY,JMAX,KSTRPV,KSTRPH 
!  AND NVERT2 IN THE MASTER INPUT FILE; OTHERWISE THE PROGRAM WILL      
!  ABORT:LMX MUST BE THE MAXIMUM NUMBER OF LEVELS (PRESSURE--OR--SIGMA) 
!  IMAXN AND IMXC ARE NESTED AND NON-EXPANDED GRID DIMENSIONS. IMXC IS  
!  EQUAL TO IMX IF YOU ARE NOT USING THE EXPANDED GRID. SAME FOR J.     

      INTEGER IY,JX,IDOT
      REAL*4  XLAT(IY,JX),XLON(IY,JX), CORIOL(IY,JX),XMAP(IY,JX)
      REAL*4  CLON,CLAT,POLLON,POLLAT,DS
!
      REAL*4  XOMEGA,d2r,r2d,A,CNTRJ,CNTRI,DDEG,XOFF,YOFF
      REAL*4  XR,YR,X,Y,FAI,XOMEGA2
      INTEGER I,J
!
      XOMEGA = 7.2722E-5                       ! ANG. ROT OF EARTH IN S**-1
      d2r = ATAN(1.)/45.                     ! CONVERT DEGREES TO RADIANS
      r2d = 1./d2r                         ! CONVERT RADIANS TO DEGREES
      A = 6371229.                             ! RADIUS OF EARTH IN METERS
!-----CENTER OF GRID
      CNTRJ = (JX+IDOT)/2.
      CNTRI = (IY+IDOT)/2.

      DDEG=DS*r2d/A                          ! GRID SPACING IN DEGREES
      XOFF=CLON-POLLON
      YOFF=clat-pollat
!-----CALCULATE X AND Y POSITIONS OF GRID                               
      DO 41 I=1,IY
      DO 41 J=1,JX                                                    
        XR = XOFF + (J-CNTRJ)*DDEG
        YR = YOFF + (I-CNTRI)*DDEG
!-----NOW CALCULATE LAT AND LON OF THIS POINT                           
!-----  ROTATE COORDINATES BACK TO NONRATED POLE
        CALL ROT2NROT(XR,YR,POLLON,pollat,X,Y)
        XLON(I,J) = X
        XLAT(I,J) = Y
        FAI = d2r*YR
        XMAP(I,J) = 1.0/COS(FAI)
   41 CONTINUE

      IF(IDOT.EQ.1) THEN
         XOMEGA2=2.*XOMEGA
         DO 45 I=1,IY
         DO 45 J=1,JX
            CORIOL(I,J)=XOMEGA2*SIN(XLAT(I,J)*d2r)
   45    CONTINUE
      END IF

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MAPPOL( XLON, XLAT, XMAP, CORIOL, IY, JX
     A                 , CLON, CLAT, DELX, idot )
      implicit none
      INTEGER IY,JX
      REAL*4  CLON,CLAT,DELX

      REAL*4  XLON(IY,JX), XLAT(IY,JX)
      REAL*4  XMAP(IY,JX), CORIOL(IY,JX)
      REAL*4  XN,PSI1,PI
!
!         XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!         PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!         PI IS PI.
!
!---------------------------------------------------------------------
      REAL*4  AA,POLE,d2r
      REAL*4  CNTRJ,CNTRI,PSX,CELL,CELL2,R,XCNTR,YCNTR
      REAL*4  X,Y,FLP,FLPP,CELL1,PSIX,xomega,xomeg2
      INTEGER II1,JJ1,I,J,idot

!  COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
!  LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND 60.N.

!  IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
!  CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
!  DELX IS GRID SPACING IN METERS.
!  ALWAYS FOR CROSS GRID.

!
      AA=6.371229E6
      XN=1.0
      PI=ATAN(1.)*4.
      d2r=PI/180.
!
      POLE=90.
      PSI1=30.
      PSI1=PSI1*d2r
      IF (CLAT.LT.0.0) THEN
         POLE = -90.0
         PSI1 = -30.
         PSI1=PSI1*d2r
      ENDIF
      CNTRJ = FLOAT(JX+idot)/2.
      CNTRI = FLOAT(IY+idot)/2.
!
      PSX = (POLE-CLAT)*d2r
      CELL=AA*SIN(PSX)/XN
      CELL2 = (1. + COS(PSI1))/(1. + COS(PSX))
      R=CELL*(CELL2)**XN
      XCNTR=0.
      YCNTR=-R
!
      II1 = IY
      JJ1 = JX
      DO 70 I=1,II1
         Y=YCNTR+(I-CNTRI)*DELX
         DO 70 J=1,JJ1
            X=XCNTR+(J-CNTRJ)*DELX
            R = SQRT(X*X + Y*Y)
            IF(Y.EQ.0.) THEN
               IF(X.GE.0.) THEN
                  FLP=90.*d2r
               ELSE
                  FLP=-90.*d2r
               ENDIF
            ELSE
               IF (CLAT.LT.0.0) THEN
                  FLP = ATAN2(X,Y)
               ELSE
                  FLP = ATAN2(X,-Y)
               END IF
            END IF
            FLPP = FLP/XN/d2r+CLON
!           IF (FLPP.GT.180.0) FLPP = FLPP-360.0
!           IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
            XLON(I,J) = FLPP
            IF (CLAT.LT.0.0) R = -R
            CELL = R/AA
            CELL1 = CELL/(1.0 + COS(PSI1))
            CELL2=ATAN(CELL1)
            PSX=2.*CELL2/d2r
            XLAT(I,J) = POLE - PSX
            PSIX=PSX*d2r
            XMAP(I,J) = ((1.0 + COS(PSI1))/(1.0 + COS(PSIX)))**XN
 70   CONTINUE

      IF (IDOT.EQ.1) THEN
         XOMEGA = 7.2722E-5           ! ANG. ROT OF EARTH IN S**-1
         d2r = ATAN(1.)/45.         ! CONVERT DEGREES TO RADIANS
         XOMEG2=2.*XOMEGA
         DO I=1,IY
         DO J=1,JX
           CORIOL(I,J)=XOMEG2*SIN(XLAT(I,J)*d2r)
         end do
         end do
      ENDIF
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE NORMER( XLON, XLAT, XMAP, CORIOL, IY, JX
     A                 , CLON, CLAT, DELX, IDOT )
      implicit none
      INTEGER IY,JX,IDOT
      REAL*4  CLON,CLAT,DELX

      REAL*4  XLON(IY,JX), XLAT(IY,JX), CORIOL(IY,JX)
      REAL*4  XMAP(IY,JX)
      REAL*4  PI
      REAL*4  AA,POLE,d2r
      REAL*4  CNTRJ,CNTRI,CELL,XCNTR,YCNTR,C2,PHICTR
      REAL*4  PHI1,X,Y,DEGLAT,XOMEGA,XOMEGA2
      INTEGER II1,JJ1,I,J

!  COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
!  LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND 60.N.

!  IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
!  CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
!  DELX IS GRID SPACING IN METERS.
!  ALWAYS FOR CROSS GRID.

      XOMEGA = 7.2722E-5                       ! ANG. ROT OF EARTH IN S**-1
      AA=6.371229E6
      PI=ATAN(1.)*4.
      d2r = ATAN(1.)/45.                     ! CONVERT DEGREES TO RADIANS
      POLE=90.
      CNTRJ = (JX+IDOT)/2.
      CNTRI = (IY+IDOT)/2.
      IF (CLAT.LT.0.0) POLE = -90.0
!
!     FOR MERCATOR PROJECTION TRUE AT PHI1
!
      PHI1=0.
      PHI1=PHI1*d2r
      C2=AA*COS(PHI1)
      XCNTR=0.
      PHICTR=CLAT*d2r
      CELL=COS(PHICTR)/(1.+SIN(PHICTR))
      YCNTR=-C2*LOG(CELL)
!
      II1 = IY
      JJ1 = JX
      DO 70 I=1,II1
         Y=YCNTR+(I-CNTRI)*DELX
         DO 70 J=1,JJ1
            X=XCNTR+(J-CNTRJ)*DELX
!
!           CALCULATIONS FOR MERCATOR
!
            XLON(I,J)=CLON+((X-XCNTR)/C2)/d2r
            CELL=EXP(Y/C2)
            XLAT(I,J)=2.*(ATAN(CELL)/d2r)-90.
            DEGLAT=XLAT(I,J)*d2r
            XMAP(I,J)=COS(PHI1)/COS(DEGLAT)
 70   CONTINUE

      IF(IDOT.EQ.1) THEN
        XOMEGA2=2.*XOMEGA
        DO I=1,IY
        DO J=1,JX
          CORIOL(I,J)=XOMEGA2*SIN(XLAT(I,J)*d2r)
        END DO
        END DO
      ENDIF

      RETURN
      END
!
      SUBROUTINE MXMNLL(IY,JX,CLON,XLON,XLAT,ntypec)
      implicit none
!
!     PURPOSE : FINDS THE MAXIMUM AND MINIMUM LATITUDE AND LONGITUDE
!
      INTEGER IY,JX,ntypec
      REAL*4  CLON,XLON(IY,JX),XLAT(IY,JX)
!
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      integer i,j
!
      xmaxlat =  -90
      xminlat =   90
      xminlon =  999999.
      xmaxlon = -999999.
!
      do i=1,iy
      do j=1,jx
         xminlat = amin1(xminlat,xlat(i,j))
         xmaxlat = amax1(xmaxlat,xlat(i,j))
      enddo
      enddo
      do i=1,iy
      do j=1,jx
         if(clon.ge.0.0) then
            if(xlon(i,j).ge.0.0) then
               xminlon = amin1(xminlon,xlon(i,j))
               xmaxlon = amax1(xmaxlon,xlon(i,j))
            else
               if(abs(clon-xlon(i,j)).lt.
     &            abs(clon-(xlon(i,j)+360.))) then
                  xminlon = amin1(xminlon,xlon(i,j))
                  xmaxlon = amax1(xmaxlon,xlon(i,j))
               else
                  xminlon = amin1(xminlon,xlon(i,j)+360.)
                  xmaxlon = amax1(xmaxlon,xlon(i,j)+360.)
               endif
            endif
         else
            if(xlon(i,j).lt.0.0) then
               xminlon = amin1(xminlon,xlon(i,j))
               xmaxlon = amax1(xmaxlon,xlon(i,j))
            else
               if(abs(clon-xlon(i,j)).lt.
     &            abs(clon-(xlon(i,j)-360.))) then
                  xminlon = amin1(xminlon,xlon(i,j))
                  xmaxlon = amax1(xmaxlon,xlon(i,j))
               else
                  xminlon = amin1(xminlon,xlon(i,j)-360.)
                  xmaxlon = amax1(xmaxlon,xlon(i,j)-360.)
               endif
            endif
         endif
      enddo
      enddo

      print 100,xminlat,xmaxlat,xminlon,xmaxlon,ntypec
 100  format(1x,'xminlat,xmaxlat,xminlon,xmaxlon,ntypec= ',4f10.2,i10)
!--------initialize minimum lat and lon of data from tape
      grdltmn=xminlat + 5.
      grdlnmn=xminlon + 5.
!
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE nrot2rot(LAM,PHI,POLLON,POLLAT,LAMS,PHIS)
      implicit none
      REAL*4  POLLON,POLLAT
! -----------------------------------------------------------------------------
! Purpose:
!     Adaption of the DWD-Functions to convert real geographical coordinates
!     (PHI,LAM) into coordinates in the roteted system (PHIS,LAMS).
!     The rotated pole is passed trough POLLON and POLLAT. POLLON and
!     POLLAT give the origin of the new rotated grid. The first four
!     arguments are input, the second two are output. All angles are in 
!     degrees (north>0, east>0)
! History:
!     05/90   D.MAJEWSKI (DWD), G. DE MORSIER (SMA) 
!     03/93   D.BRESCH (ETHZ)
!     11/97   D.LUETHI (ETHZ)


! decalration of input/output vars:
      REAL*4  LAM,PHI,LAMS,PHIS,PLAM,PPHI
 
! decalration of internal vars:
      real*4  zlam,zphi
      real*4  zarg,zarg1,zarg2

      real*4  r2d , d2r
      real*4  ZSINPOL,ZCOSPOL,ZLAMPOL

      r2d = 45./atan(1.)
      d2r = atan(1.)/45.

      PLAM=POLLON + 180.
      PPHI = 90. - POLLAT

      IF (PLAM.GT.180.) PLAM = PLAM-360.

      ZSINPOL = SIN(d2r*PPHI)
      ZCOSPOL = COS(d2r*PPHI)
      ZLAMPOL =     d2r*PLAM

! first, the conversion of PHI to PHIS:
      ZPHI    = d2r*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = d2r*ZLAM
      ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
      PHIS = ASIN(ZARG)
      phis = log(tan(phis/2.+atan(1.)))*r2d
 
! now, the conversion for LAMS follws: 
      ZPHI    =     d2r*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = d2r*ZLAM
      ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
      ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LAMS =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
          LAMS =  90.0
        ELSE
          LAMS = -90.0
        ENDIF
      ELSE
        LAMS = r2d*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      END

      function oned(x,a,b,c,d)
      implicit none
      real*8  x,a,b,c,d
      real*8  oned
!
      oned = 0.
      if (x.eq.0.) oned = b
      if (x.eq.1.) oned = c
      if (b*c.eq.0.) return
      if (a*d.eq.0.) goto 20
      oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))
     &     +x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      return
   20 oned = b*(1.0-x)+c*x
      if (a.ne.0.0) oned = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))
      if (d.ne.0.0) oned = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
      return
      end
!
      SUBROUTINE OUTPUT(nunitc,iy,jx,nsg,dsinm,clat,clon,plat,plon,iproj
     &          ,htgrid,htsdgrid,lndout,xlat,xlon,dlat,dlon,xmap,DATTYP
     &          ,dmap,coriol,snowam,igrads,ibigend,kz,sigma,mask,ptop
     &          ,htgrid_s,lndout_s,ibyte,ngrid,truelatL,truelatH,GRDFAC
     &          ,filout,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
     &          ,AERTYP,texout,frac_tex,ntex)
      implicit none
      integer nunitc,iy,jx,nsg,kz,igrads,ibigend,ibyte,ngrid,nveg,ntex
      real*4  dsinm,clat,clon,plat,plon,GRDFAC,ptop
      real*4  truelatL,truelatH
      character*6 iproj
      character*50 filout
      character*5 DATTYP
      real*4  htgrid(iy,jx),htsdgrid(iy,jx)
      real*4  lndout(iy,jx),texout(iy,jx)
      real*4  htgrid_s(iy*nsg,jx*nsg),htgrid_a
      real*4  lndout_s(iy*nsg,jx*nsg)
      character*4 LSMTYP
      character*7 AERTYP
      real*4  sanda(iy,jx),sandb(iy,jx),claya(iy,jx),clayb(iy,jx)
      real*4  frac_lnd(iy,jx,nveg),frac_tex(iy,jx,ntex)
      real*4  xlat(iy,jx),xlon(iy,jx),xmap(iy,jx)
      real*4  dlat(iy,jx),dlon(iy,jx),dmap(iy,jx)
      real*4  coriol(iy,jx),snowam(iy,jx),sigma(kz+1),mask(iy,jx)
      integer i,j,i0,j0,m,n,k0,kk,lnd(25)
      real*4  htave
      character*25 fsubname
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx,k
      real*4  lon0,lon1,lat0,lat1
      integer III,JJJ
      logical there
!
      if(kz.eq.14) then                       ! RegCM2
         sigma(1) = 0.0
         sigma(2) = 0.04
         sigma(3) = 0.10
         sigma(4) = 0.17
         sigma(5) = 0.25
         sigma(6) = 0.35
         sigma(7) = 0.46
         sigma(8) = 0.56
         sigma(9) = 0.67
         sigma(10)= 0.77
         sigma(11)= 0.86
         sigma(12)= 0.93
         sigma(13)= 0.97
         sigma(14)= 0.99
         sigma(15)= 1.0
      else if(kz.eq.18) then                  ! RegCM3, default
         sigma(1) = 0.0
         sigma(2) = 0.05
         sigma(3) = 0.10
         sigma(4) = 0.16
         sigma(5) = 0.23
         sigma(6) = 0.31
         sigma(7) = 0.39
         sigma(8) = 0.47
         sigma(9) = 0.55
         sigma(10)= 0.63
         sigma(11)= 0.71
         sigma(12)= 0.78
         sigma(13)= 0.84
         sigma(14)= 0.89
         sigma(15)= 0.93
         sigma(16)= 0.96
         sigma(17)= 0.98
         sigma(18)= 0.99
         sigma(19)= 1.0
      else if(kz.eq.23) then                  ! MM5V3
         sigma(1) = 0.0
         sigma(2) = 0.05
         sigma(3) = 0.1
         sigma(4) = 0.15
         sigma(5) = 0.2
         sigma(6) = 0.25
         sigma(7) = 0.3
         sigma(8) = 0.35
         sigma(9) = 0.4
         sigma(10)= 0.45
         sigma(11)= 0.5
         sigma(12)= 0.55
         sigma(13)= 0.6
         sigma(14)= 0.65
         sigma(15)= 0.7
         sigma(16)= 0.75
         sigma(17)= 0.8
         sigma(18)= 0.85
         sigma(19)= 0.89
         sigma(20)= 0.93
         sigma(21)= 0.96
         sigma(22)= 0.98
         sigma(23)= 0.99
         sigma(24)= 1.0
      else
         write(*,*) 'You vertical level number is not 14, 18, or 23'
         write(*,*) 'Please set your sigma parameters in OUTPUT'
         stop
      endif

      if(NSG.gt.1) then
         if(NSG.lt.10) then
            write(fsubname,9) NSG
         else
            write(fsubname,11) NSG
         endif
 9    format('../../Input/DOMAIN',I1,'.INFO')
11    format('../../Input/DOMAIN',I2,'.INFO')
         write(*,*) fsubname
         inquire(file=fsubname,exist=there)
         if(.not.there) then
            write(*,*) 'Subgrid Terrain and Landuse must be available'
            stop
         endif
         open(19,file=fsubname,form='unformatted'
     &          ,recl=iy*jx*nsg*nsg*ibyte,access='direct')
         read(19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         read(19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         do i=1,iy*nsg
         do j=1,jx*nsg
          if((LSMTYP.eq.'BATS'.and.(lndout_s(i,j).gt.20..or.
     &                              lndout_s(i,j).lt.0.))
     &   .or.(LSMTYP.eq.'USGS'.and.(lndout_s(i,j).gt.25..or.
     &                              lndout_s(i,j).lt.0.))) then
           print*,i,j,lndout_s(i,j)
           stop 999
          end if
         end do
         end do
         do i=1,iy
         do j=1,jx
            i0=(i-1)*nsg
            j0=(j-1)*nsg
            htave=0.0
            do m=1,nsg
            do n=1,nsg
               htave=htave+htgrid_s(i0+m,j0+n)
            enddo
            enddo
            htgrid_a=htave/float(nsg*nsg)
            do m=1,nsg
            do n=1,nsg
               htgrid_s(i0+m,j0+n)=htgrid_s(i0+m,j0+n)-htgrid_a
     &                                                +htgrid(i,j)
            enddo
            enddo
!           htgrid(i,j)=htave/float(nsg*nsg)
         enddo
         enddo
         do i=1,iy
         do j=1,jx
            i0=(i-1)*nsg
            j0=(j-1)*nsg
            do k=1,nveg
               lnd(k)=0
            enddo
            do m=1,nsg
            do n=1,nsg
               k=nint(lndout_s(i0+m,j0+n))
               lnd(k)=lnd(k)+1
            enddo
            enddo
            k0=0
            kk=0
            do k=1,nveg
               if(LSMTYP.eq.'BATS') then
                  if(k.ne.15) then
                     if(lnd(k).gt.k0) then
                        k0=lnd(k)
                        kk=k
                     endif
                  endif
               else if(LSMTYP.eq.'USGS') then
                  if(k.ne.25) then
                     if(lnd(k).gt.k0) then
                        k0=lnd(k)
                        kk=k
                     endif
                  endif
               endif
            enddo
            if(LSMTYP.eq.'BATS'.and.
     &         float(lnd(15))/float(nsg*nsg).gt.0.75) then
               lndout(i,j)=15.
            else if(LSMTYP.eq.'USGS'.and.
     &         float(lnd(25))/float(nsg*nsg).gt.0.75) then
               lndout(i,j)=25.
            endif
            if(LSMTYP.eq.'BATS') then
               if(htgrid(i,j).gt.1.0.and.(lndout(i,j).gt.14.5.and.
     &                                    lndout(i,j).lt.15.5)) then
                  lndout(i,j)=float(kk)
               else if(htgrid(i,j).lt.0.1.and.(lndout(i,j).gt.14.5.and.
     &                                         lndout(i,j).lt.15.5))then
                  do m=1,nsg
                  do n=1,nsg
                     lndout_s(i0+m,j0+n)=15.0
                  enddo
                  enddo
               endif
            else if(LSMTYP.eq.'USGS') then
               if(htgrid(i,j).gt.1.0.and.lndout(i,j).gt.24.5) then
                  lndout(i,j)=float(kk)
               else if(htgrid(i,j).lt.0.1.and.lndout(i,j).gt.24.5) then
                  do m=1,nsg
                  do n=1,nsg
                     lndout_s(i0+m,j0+n)=25.0
                  enddo
                  enddo
               endif
            endif
         enddo
         enddo
         write(19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         write(19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         close(19)
      endif
         
      do i=1,iy
      do j=1,jx
         if ((LSMTYP.eq.'BATS'.and.
     &         (lndout(i,j).gt.20. .or. lndout(i,j).lt.0.))
     &   .or.(LSMTYP.eq.'USGS'.and.
     &        (lndout(i,j).gt.25. .or. lndout(i,j).lt.0.))) then
           print*,i,j,lndout(i,j)
           stop 999
         end if
      end do
      end do
      if(ngrid.ne.nsg.or.(DATTYP.NE.'FVGCM'.and.DATTYP.NE.'NRP2W'
     &               .and.DATTYP.NE.'GFS11'.and.DATTYP.NE.'EH5OM')) then
      write(nunitc,rec=1) iy,jx,kz,dsinm,clat,clon,plat,plon,GRDFAC
     &                  , iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend
     &                  , truelatL,truelatH
      else
         write(*,*) 'please input lon0,lon1,lat0,lat1'
         write(*,*) 'Note: lon0 < lon1, and lat0 < lat1'
         read(*,*) lon0,lon1,lat0,lat1
      write(nunitc,rec=1) iy,jx,kz,dsinm,clat,clon,plat,plon,GRDFAC
     &                  , iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend
     &                  , truelatL,truelatH,lon0,lon1,lat0,lat1
      endif
      if(ngrid.eq.nsg) then
        open(23,file='../ICBC/icbcWIN.param')
        write(23,'(a)') '      INTEGER III'
        write(23,'(a)') '      INTEGER JJJ'
        IF(DATTYP.EQ.'NRP2W') THEN
           III=NINT((lon1-lon0)/2.5)+1
           JJJ=NINT((lat1-lat0)/2.5)+1
           write(23,101) 'III   =',III
           write(23,101) 'JJJ   =',JJJ
        ELSE
           write(23,101) 'III   =',1
           write(23,101) 'JJJ   =',1
        ENDIF
 101    format('      parameter(',A8,I4,')')
        close(23)
      endif
      write(nunitc,rec=2) ((htgrid(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=3) ((htsdgrid(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=4) ((lndout(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=5) ((xlat(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=6) ((xlon(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=7) ((dlat(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=8) ((dlon(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=9) ((xmap(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=10) ((dmap(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=11) ((coriol(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=12) ((snowam(i,j),j=1,jx),i=1,iy)
      do i=1,iy
      do j=1,jx
         if(LSMTYP.eq.'BATS') then
            if(lndout(i,j).gt.13.5.and.lndout(i,j).lt.15.5) then
               mask(i,j)=0.0
            else
               mask(i,j)=2.0
            endif
         else if(LSMTYP.eq.'USGS') then
            if(lndout(i,j).gt.24.5) then
               mask(i,j)=0.0
            else
               mask(i,j)=2.0
            endif
         endif
      enddo
      enddo
      write(nunitc,rec=13) ((mask(i,j),j=1,jx),i=1,iy)
      if(LSMTYP.eq.'USGS') then
        write(nunitc,rec=14) ((sanda(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=15) ((sandb(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=16) ((claya(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=17) ((clayb(i,j),j=1,jx),i=1,iy)
      endif
      if(AERTYP(7:7).eq.'1') then
        if(LSMTYP.eq.'BATS') then
          write(nunitc,rec=14) ((texout(i,j),j=1,jx),i=1,iy)
        else if(LSMTYP.eq.'USGS') then
          write(nunitc,rec=18) ((texout(i,j),j=1,jx),i=1,iy)
        endif
      endif
      close(nunitc)

      if(igrads.eq.1) then
         if(ngrid.ne.nsg) then
           if(ngrid.lt.10) then
             write(31,12) filout(13:24)
           else
             write(31,13) filout(13:25)
           endif
         else
           write(31,10) 
         endif
  10  format('dset ^DOMAIN.INFO')
  12  format('dset ^',A12)
  13  format('dset ^',A13)
         write(31,20)
  20  format('title RegCM domain information')
         if(ibigend.eq.1) then
            write(31,30)
  30  format('options big_endian')
         else
            write(31,40)
  40  format('options little_endian')
         endif
         write(31,50)
  50  format('undef -9999.')
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
            rlatinc=dsinm*0.001/111./2.
            rloninc=dsinm*0.001/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))

            centerj=jx/2.
            centeri=iy/2.
         endif
         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx,iy,clat,clon,centerj,centeri,
     &                    truelatL,truelatH,clon,dsinm,dsinm
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
            write(31,230) jx,iy,plon,plat,dsinm/111000.
     &                                   ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
         write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
         write(31,400) 1
 400  format('tdef ',I1,' linear 00z01Jan2001 1mo')
         if(LSMTYP.eq.'BATS') then
            if(AERTYP(7:7).eq.'1') then
               write(31,500) 13+1
            else
               write(31,500) 13
            endif
         else if(LSMTYP.eq.'USGS') then
            if(AERTYP(7:7).eq.'1') then
               write(31,500) 17+1
            else
               write(31,500) 17
            endif
         endif
 500  format('vars ',I2)
         write(31,600) 'head    ','header information         '
 600  format(A8,'0 99 ',A26)
         write(31,600) 'ht      ','surface elevation          '
         write(31,600) 'htsd    ','surface elevation std. dev.'
         write(31,600) 'landuse ','surface landuse type       '
         write(31,600) 'xlat    ','latitude  of cross points  '
         write(31,600) 'xlon    ','longitude of cross points  '
         write(31,600) 'dlat    ','latitude  of dot points    '
         write(31,600) 'dlon    ','longitude of dot points    '
         write(31,600) 'xmap    ','map factors of cross points'
         write(31,600) 'dmap    ','map factors of dot points  '
         write(31,600) 'coriol  ','coriol force               '
         write(31,600) 'snowam  ','initial snow amount        '
         write(31,600) 'mask    ','land/sea mask              '
         if(LSMTYP.eq.'USGS') then
           write(31,600) 'sanda   ','sand percentage (0-30cm)   '
           write(31,600) 'sandb   ','sand percentage (30-100cm) '
           write(31,600) 'claya   ','clay percentage (0-30cm)   '
           write(31,600) 'clayb   ','clay percentage (30-100cm) '
         endif
         if(AERTYP(7:7).eq.'1') then
           write(31,600) 'texture ','soil texture               '
         endif
         write(31,700)
 700  format('endvars')
      close(31)
      endif
!
      return
      end
      SUBROUTINE OUTPUT2(nunitc,iy,jx,nsg,dsinm,clat,clon,plat,plon
     &    ,iproj,htgrid,htsdgrid,lndout,xlat,xlon,dlat,dlon,xmap,DATTYP
     &          ,dmap,coriol,snowam,igrads,ibigend,kz,sigma,mask,ptop
     &          ,htgrid_s,lndout_s,ibyte,ngrid,truelatL,truelatH,GRDFAC
     &          ,filout,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
     &          ,AERTYP,texout,frac_tex,ntex)
      implicit none
      integer nunitc,iy,jx,nsg,kz,igrads,ibigend,ibyte,ngrid,nveg,ntex
      real*4  dsinm,clat,clon,plat,plon,GRDFAC,ptop
      real*4  truelatL,truelatH
      character*6 iproj
      character*50 filout
      character*5 DATTYP
      real*4  htgrid(iy,jx),htsdgrid(iy,jx)
      real*4  lndout(iy,jx),texout(iy,jx)
      real*4  htgrid_s(iy*nsg,jx*nsg),htgrid_a
      real*4  lndout_s(iy*nsg,jx*nsg)
      character*4 LSMTYP
      character*7 AERTYP
      real*4  sanda(iy,jx),sandb(iy,jx),claya(iy,jx),clayb(iy,jx)
      real*4  frac_lnd(iy,jx,nveg),frac_tex(iy,jx,ntex)
      real*4  xlat(iy,jx),xlon(iy,jx),xmap(iy,jx)
      real*4  dlat(iy,jx),dlon(iy,jx),dmap(iy,jx)
      real*4  coriol(iy,jx),snowam(iy,jx),sigma(kz+1),mask(iy,jx)
      integer i,j,i0,j0,m,n,k0,kk,lnd(25)
      real*4  htave
      character*25 fsubname
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx,k
      real*4  lon0,lon1,lat0,lat1
      integer III,JJJ
      logical there
!
      if(kz.eq.14) then                       ! RegCM2
         sigma(1) = 0.0
         sigma(2) = 0.04
         sigma(3) = 0.10
         sigma(4) = 0.17
         sigma(5) = 0.25
         sigma(6) = 0.35
         sigma(7) = 0.46
         sigma(8) = 0.56
         sigma(9) = 0.67
         sigma(10)= 0.77
         sigma(11)= 0.86
         sigma(12)= 0.93
         sigma(13)= 0.97
         sigma(14)= 0.99
         sigma(15)= 1.0
      else if(kz.eq.18) then                  ! RegCM3, default
         sigma(1) = 0.0
         sigma(2) = 0.05
         sigma(3) = 0.10
         sigma(4) = 0.16
         sigma(5) = 0.23
         sigma(6) = 0.31
         sigma(7) = 0.39
         sigma(8) = 0.47
         sigma(9) = 0.55
         sigma(10)= 0.63
         sigma(11)= 0.71
         sigma(12)= 0.78
         sigma(13)= 0.84
         sigma(14)= 0.89
         sigma(15)= 0.93
         sigma(16)= 0.96
         sigma(17)= 0.98
         sigma(18)= 0.99
         sigma(19)= 1.0
      else if(kz.eq.23) then                  ! MM5V3
         sigma(1) = 0.0
         sigma(2) = 0.05
         sigma(3) = 0.1
         sigma(4) = 0.15
         sigma(5) = 0.2
         sigma(6) = 0.25
         sigma(7) = 0.3
         sigma(8) = 0.35
         sigma(9) = 0.4
         sigma(10)= 0.45
         sigma(11)= 0.5
         sigma(12)= 0.55
         sigma(13)= 0.6
         sigma(14)= 0.65
         sigma(15)= 0.7
         sigma(16)= 0.75
         sigma(17)= 0.8
         sigma(18)= 0.85
         sigma(19)= 0.89
         sigma(20)= 0.93
         sigma(21)= 0.96
         sigma(22)= 0.98
         sigma(23)= 0.99
         sigma(24)= 1.0
      else
         write(*,*) 'You vertical level number is not 14, 18, or 23'
         write(*,*) 'Please set your sigma parameters in OUTPUT'
         stop
      endif

      if(NSG.gt.1) then
         if(NSG.lt.10) then
            write(fsubname,9) NSG
         else
            write(fsubname,11) NSG
         endif
 9    format('../../Input/DOMAIN',I1,'.INFO')
11    format('../../Input/DOMAIN',I2,'.INFO')
         write(*,*) fsubname
         inquire(file=fsubname,exist=there)
         if(.not.there) then
            write(*,*) 'Subgrid Terrain and Landuse must be available'
            stop
         endif
         open(19,file=fsubname,form='unformatted'
     &          ,recl=iy*jx*nsg*nsg*ibyte,access='direct')
         read(19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         read(19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         do i=1,iy*nsg
         do j=1,jx*nsg
          if((LSMTYP.eq.'BATS'.and.(lndout_s(i,j).gt.20..or.
     &                              lndout_s(i,j).lt.0.))
     &   .or.(LSMTYP.eq.'USGS'.and.(lndout_s(i,j).gt.25..or.
     &                              lndout_s(i,j).lt.0.))) then
           print*,i,j,lndout_s(i,j)
           stop 999
          end if
         end do
         end do
         do i=1,iy
         do j=1,jx
            i0=(i-1)*nsg
            j0=(j-1)*nsg
            htave=0.0
            do m=1,nsg
            do n=1,nsg
               if(LSMTYP.eq.'BATS') then
                  if(htgrid(i,j).lt.0.1.and.(lndout(i,j).gt.14.5.and.
     &                                       lndout(i,j).lt.15.5)) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 15.
                  endif
               else if(LSMTYP.eq.'USGS') then
                  if(htgrid(i,j).lt.0.1.and.lndout(i,j).gt.24.5) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 25.
                  endif
               endif
               htave=htave+htgrid_s(i0+m,j0+n)
            enddo
            enddo
            htgrid_a=htave/float(nsg*nsg)
            do m=1,nsg
            do n=1,nsg
               htgrid_s(i0+m,j0+n)=htgrid_s(i0+m,j0+n)-htgrid_a
     &                                                +htgrid(i,j)
            enddo
            enddo
         enddo
         enddo
         write(19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         write(19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
         close(19)
      endif
         
      do i=1,iy
      do j=1,jx
         if ((LSMTYP.eq.'BATS'.and.
     &         (lndout(i,j).gt.20. .or. lndout(i,j).lt.0.))
     &   .or.(LSMTYP.eq.'USGS'.and.
     &        (lndout(i,j).gt.25. .or. lndout(i,j).lt.0.))) then
           print*,i,j,lndout(i,j)
           stop 999
         end if
      end do
      end do
      if(ngrid.ne.nsg.or.(DATTYP.NE.'FVGCM'.and.DATTYP.NE.'NRP2W'
     &               .and.DATTYP.NE.'GFS11'.and.DATTYP.NE.'EH5OM')) then
      write(nunitc,rec=1) iy,jx,kz,dsinm,clat,clon,plat,plon,GRDFAC
     &                  , iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend
     &                  , truelatL,truelatH
      else
         write(*,*) 'please input lon0,lon1,lat0,lat1'
         write(*,*) 'Note: lon0 < lon1, and lat0 < lat1'
         read(*,*) lon0,lon1,lat0,lat1
      write(nunitc,rec=1) iy,jx,kz,dsinm,clat,clon,plat,plon,GRDFAC
     &                  , iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend
     &                  , truelatL,truelatH,lon0,lon1,lat0,lat1
      endif
      if(ngrid.eq.nsg) then
        open(23,file='../ICBC/icbcWIN.param')
        write(23,'(a)') '      INTEGER III'
        write(23,'(a)') '      INTEGER JJJ'
        IF(DATTYP.EQ.'NRP2W') THEN
           III=NINT((lon1-lon0)/2.5)+1
           JJJ=NINT((lat1-lat0)/2.5)+1
           write(23,101) 'III   =',III
           write(23,101) 'JJJ   =',JJJ
        ELSE
           write(23,101) 'III   =',1
           write(23,101) 'JJJ   =',1
        ENDIF
 101    format('      parameter(',A8,I4,')')
        close(23)
      endif
      write(nunitc,rec=2) ((htgrid(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=3) ((htsdgrid(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=4) ((lndout(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=5) ((xlat(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=6) ((xlon(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=7) ((dlat(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=8) ((dlon(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=9) ((xmap(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=10) ((dmap(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=11) ((coriol(i,j),j=1,jx),i=1,iy)
      write(nunitc,rec=12) ((snowam(i,j),j=1,jx),i=1,iy)
      do i=1,iy
      do j=1,jx
         if(LSMTYP.eq.'BATS') then
            if(lndout(i,j).gt.13.5.and.lndout(i,j).lt.15.5) then
               mask(i,j)=0.0
            else
               mask(i,j)=2.0
            endif
         else if(LSMTYP.eq.'USGS') then
            if(lndout(i,j).gt.24.5) then
               mask(i,j)=0.0
            else
               mask(i,j)=2.0
            endif
         endif
      enddo
      enddo
      write (nunitc,rec=13) ((mask(i,j),j=1,jx),i=1,iy)
      if(LSMTYP.eq.'USGS') then
        write(nunitc,rec=14) ((sanda(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=15) ((sandb(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=16) ((claya(i,j),j=1,jx),i=1,iy)
        write(nunitc,rec=17) ((clayb(i,j),j=1,jx),i=1,iy)
        do k=1,nveg
           write(nunitc,rec=17+k)((frac_lnd(i,j,k),j=1,jx),i=1,iy)
        enddo
      endif
      if(AERTYP(7:7).eq.'1') then
        if(LSMTYP.eq.'BATS') then
          write(nunitc,rec=14) ((texout(i,j),j=1,jx),i=1,iy)
          do k=1,ntex
             write(nunitc,rec=14+k)((frac_tex(i,j,k),j=1,jx),i=1,iy)
          enddo
        else if(LSMTYP.eq.'USGS') then
          write(nunitc,rec=18+nveg) ((texout(i,j),j=1,jx),i=1,iy)
          do k=1,ntex
          write(nunitc,rec=18+nveg+k)((frac_tex(i,j,k),j=1,jx),i=1,iy)
          enddo
        endif
      endif
      close(nunitc)

      if(igrads.eq.1) then
         if(ngrid.ne.nsg) then
           if(ngrid.lt.10) then
             write(31,12) filout(13:24)
           else
             write(31,13) filout(13:25)
           endif
         else
           write(31,10) 
         endif
  10  format('dset ^DOMAIN.INFO')
  12  format('dset ^',A12)
  13  format('dset ^',A13)
         write(31,20)
  20  format('title RegCM domain information')
         if(ibigend.eq.1) then
            write(31,30)
  30  format('options big_endian')
         else
            write(31,40)
  40  format('options little_endian')
         endif
         write(31,50)
  50  format('undef -9999.')
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
            rlatinc=dsinm*0.001/111./2.
            rloninc=dsinm*0.001/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))

            centerj=jx/2.
            centeri=iy/2.
         endif
         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx,iy,clat,clon,centerj,centeri,
     &                    truelatL,truelatH,clon,dsinm,dsinm
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
            write(31,230) jx,iy,plon,plat,dsinm/111000.
     &                                   ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
         write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
         write(31,400) 1
 400  format('tdef ',I1,' linear 00z01Jan2001 1mo')
         if(LSMTYP.eq.'BATS') then
            if(AERTYP(7:7).eq.'1') then
               write(31,500) 13+ntex+1
            else
               write(31,500) 13
            endif
         else if(LSMTYP.eq.'USGS') then
            if(AERTYP(7:7).eq.'1') then
               write(31,500) 42+ntex+1
            else
               write(31,500) 42
            endif
         endif
 500  format('vars ',I2)
         write(31,600) 'head    ','header information         '
 600  format(A8,'0 99 ',A26)
         write(31,600) 'ht      ','surface elevation          '
         write(31,600) 'htsd    ','surface elevation std. dev.'
         write(31,600) 'landuse ','surface landuse type       '
         write(31,600) 'xlat    ','latitude  of cross points  '
         write(31,600) 'xlon    ','longitude of cross points  '
         write(31,600) 'dlat    ','latitude  of dot points    '
         write(31,600) 'dlon    ','longitude of dot points    '
         write(31,600) 'xmap    ','map factors of cross points'
         write(31,600) 'dmap    ','map factors of dot points  '
         write(31,600) 'coriol  ','coriol force               '
         write(31,600) 'snowam  ','initial snow amount        '
         write(31,600) 'mask    ','land/sea mask              '
         if(LSMTYP.eq.'USGS') then
           write(31,600) 'sanda   ','sand percentage (0-30cm)   '
           write(31,600) 'sandb   ','sand percentage (30-100cm) '
           write(31,600) 'claya   ','clay percentage (0-30cm)   '
           write(31,600) 'clayb   ','clay percentage (30-100cm) '
           write(31,600) 'per1    ','percentage landuse type 1  '
           write(31,600) 'per2    ','percentage landuse type 2  '
           write(31,600) 'per3    ','percentage landuse type 3  '
           write(31,600) 'per4    ','percentage landuse type 4  '
           write(31,600) 'per5    ','percentage landuse type 5  '
           write(31,600) 'per6    ','percentage landuse type 6  '
           write(31,600) 'per7    ','percentage landuse type 7  '
           write(31,600) 'per8    ','percentage landuse type 8  '
           write(31,600) 'per9    ','percentage landuse type 9  '
           write(31,600) 'per10   ','percentage landuse type 10 '
           write(31,600) 'per11   ','percentage landuse type 11 '
           write(31,600) 'per12   ','percentage landuse type 12 '
           write(31,600) 'per13   ','percentage landuse type 13 '
           write(31,600) 'per14   ','percentage landuse type 14 '
           write(31,600) 'per15   ','percentage landuse type 15 '
           write(31,600) 'per16   ','percentage landuse type 16 '
           write(31,600) 'per17   ','percentage landuse type 17 '
           write(31,600) 'per18   ','percentage landuse type 18 '
           write(31,600) 'per19   ','percentage landuse type 19 '
           write(31,600) 'per20   ','percentage landuse type 20 '
           write(31,600) 'per21   ','percentage landuse type 21 '
           write(31,600) 'per22   ','percentage landuse type 22 '
           write(31,600) 'per23   ','percentage landuse type 23 '
           write(31,600) 'per24   ','percentage landuse type 24 '
           write(31,600) 'per25   ','percentage landuse type 24 '
         endif
         if(AERTYP(7:7).eq.'1') then
           write(31,600) 'texture ','soil texture               '
           write(31,600) 'text01  ','Sand              frac.    '
           write(31,600) 'text02  ','Loamy Sand        frac.    '
           write(31,600) 'text03  ','Sandy Loam        frac.    '
           write(31,600) 'text04  ','Silt Loam         frac.    '
           write(31,600) 'text05  ','Silt              frac.    '
           write(31,600) 'text06  ','Loam              frac.    '
           write(31,600) 'text07  ','Sandy Clay Loam   frac.    '
           write(31,600) 'text08  ','Silty Clay Loam   frac.    '
           write(31,600) 'text09  ','Clay Loam         frac.    '
           write(31,600) 'text10  ','Sandy Clay        frac.    '
           write(31,600) 'text11  ','Silty Clay        frac.    '
           write(31,600) 'text12  ','Clay              frac.    '
           write(31,600) 'text13  ','OM                frac.    '
           write(31,600) 'text14  ','Water             frac.    '
           write(31,600) 'text15  ','Bedrock           frac.    '
           write(31,600) 'text16  ','Other             frac.    '
           write(31,600) 'text17  ','No data           frac.    '
         endif
         write(31,700)
 700  format('endvars')
      close(31)
      endif
!
      return
      end
!
      SUBROUTINE RDLDTR
      implicit none
!     include '../../Commons/env/include/netcdf.inc'
      include 'domain.param'
!
      integer nxmax, nymax
      parameter (nxmax=21600/ntypec,nymax=nxmax/2)
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay
!
      integer ltype,lrec
      real*4  stores
      common/a/ stores(50),ltype
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
!
      integer ilat1, ilat2
      integer iflag, ierr
      character filelev*36, filcat*36, filusgs*11, filsand*7, filclay*7
      character filtext*10

      integer*2 isand(nxmax,2)
      integer*2 iclay(nxmax,2)
      integer*2 iusgs(nxmax,nveg)
      character*1 ch_tex(nxmax,ntex)
      character*1 ch_cat(nxmax,nveg)
      character*2 ch_topo(nxmax),char_2,ch_htsd(nxmax)

!
      real*4  center,rlat,rlon
      integer irec,irect,j,k
      integer ihmax,jhmax
      logical there
!
      if (ntypec.lt.10) then
        write (filelev,101) ntypec
 101    format('../DATA/SURFACE/GTOPO30_',i1,'MIN.dat')
      else if (ntypec.lt.100) then
        write (filelev,102) ntypec
 102    format('../DATA/SURFACE/GTOPO30_',i2,'MIN.dat')
      else
        write(*,*) 'For terrain, ntypec is not set correctly ',ntypec
        stop 'subroutine RDLDTR'
      end if
      inquire(file=filelev,exist=there)
      if(.not.there) then
        print*,'ERROR OPENING ',filelev,' FILE:  FILE DOES NOT EXIST'
        STOP '4810 IN SUBROUTINE RDLDTR'
      endif
      open(46,file=filelev,form='unformatted'
     &       ,recl=nxmax*ibyte/2,access='direct')
      
      if(LSMTYP.eq.'BATS') then
         if (ntypec.lt.10) then
           write (filcat,201) ntypec
 201       format('../DATA/SURFACE/GLCC',i1,'MIN_BATS.dat')
         else if (ntypec.lt.100) then
           write (filcat,202) ntypec
 202       format('../DATA/SURFACE/GLCC',i2,'MIN_BATS.dat')
         else
           write(*,*) 'For landuse, ntypec is not set correctly ',ntypec
           stop 'subroutine RDLDTR'
         end if
         inquire(file=filcat,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ',filcat, ' FILE:  FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR'
         endif
         open(47,file=filcat,form='unformatted'
     &          ,recl=nxmax*nveg*ibyte/4,access='direct')
      else if(LSMTYP.eq.'USGS') then
         if (ntypec.lt.10) then
           write (filusgs,301) ntypec
 301       format('VEG-USGS.0',i1)
           write (filsand,401) ntypec
 401       format('SAND.0',i1)
           write (filclay,501) ntypec
 501       format('CLAY.0',i1)
         else if (ntypec.lt.100) then
           write (filusgs,302) ntypec
 302       format('VEG-USGS.',i2)
           write (filsand,402) ntypec
 402       format('SAND.',i2)
           write (filclay,502) ntypec
 502       format('CLAY.',i2)
         else
           write(*,*) 'Error input of ntypec: ',ntypec
         endif
         if(ntypec.eq.2) then
            inquire(file='../DATA/SURFACE/'//filusgs//'A',exist=there)
            if(.not.there) then
               print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &               ,' FILE: FILE DOES NOT EXIST'
               STOP '4820 IN SUBROUTINE RDLDTR'
            endif
            inquire(file='../DATA/SURFACE/'//filusgs//'B',exist=there)
            if(.not.there) then
               print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &               ,' FILE: FILE DOES NOT EXIST'
               STOP '4820 IN SUBROUTINE RDLDTR'
            endif
         else
         inquire(file='../DATA/SURFACE/'//filusgs,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR'
         endif
         endif
         inquire(file='../DATA/SURFACE/'//filsand,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filsand
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR'
         endif
         inquire(file='../DATA/SURFACE/'//filclay,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filclay
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR'
         endif
         if(ntypec.eq.2) then
            open(41,file='../DATA/SURFACE/'//filusgs//'A'
     &      ,form='unformatted',recl=nxmax*ibyte/2,access='direct')
            open(44,file='../DATA/SURFACE/'//filusgs//'B'
     &      ,form='unformatted',recl=nxmax*ibyte/2,access='direct')
         else
            open(41,file='../DATA/SURFACE/'//filusgs,form='unformatted'
     &             ,recl=nxmax*ibyte/2,access='direct')
         endif
         open(42,file='../DATA/SURFACE/'//filsand,form='unformatted'
     &         ,recl=nxmax*ibyte/2,access='direct')
         open(43,file='../DATA/SURFACE/'//filclay,form='unformatted'
     &         ,recl=nxmax*ibyte/2,access='direct')
      endif
      if(AERTYP(7:7).eq.'1') then
         if (ntypec.lt.10) then
           write (filtext,303) ntypec
 303       format('SOILCAT.0',i1)
         else if (ntypec.lt.100) then
           write (filtext,304) ntypec
 304       format('SOILCAT.',i2)
         else
           write(*,*) 'Error input of ntypec: ',ntypec
         endif
         inquire(file='../DATA/SURFACE/'//filtext,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filtext
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4830 IN SUBROUTINE RDLDTR'
         endif
         open(45,file='../DATA/SURFACE/'//filtext,form='unformatted'
     &         ,recl=nxmax*ntex*ibyte/4,access='direct')
      endif

      lrec=0
      rewind(48)
      center = xnc
      ilat1 = nint((90.-xmaxlat)/xnc)
      ilat1 = max(min(nymax,ilat1),1)
      ilat2 = nint((90.-xminlat)/xnc)+1
      ilat2 = max(min(nymax,ilat2),1)
      do 30 irec=ilat1,ilat2
        rlat = 90.-center*irec+center/2.
        irect = nymax - irec + 1
        read(46,rec=irect)(ch_topo(j),j=1,nxmax)
!       print*,'   ELEVATION READ IN', ilat1,ilat2,irec
        read(46,rec=nymax+irect)(ch_htsd(j),j=1,nxmax)
!       print*,'   ELEVATION STD DEV READ IN', ilat1,ilat2,irec
        if(LSMTYP.eq.'BATS') then
          read(47,rec=irect)((ch_cat(j,k),k=1,nveg),j=1,nxmax)
        else if(LSMTYP.eq.'USGS') then
          if(ntypec.eq.2) then
             do k=1,13
                read(41,rec=nymax*(k-1)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
             do k=14,nveg
                read(44,rec=nymax*(k-14)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
          else
             do k=1,nveg
                read(41,rec=nymax*(k-1)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
          endif
          read(42,rec=      irec)(isand(j,1),j=1,nxmax)
          read(42,rec=nymax+irec)(isand(j,2),j=1,nxmax)
          read(43,rec=      irec)(iclay(j,1),j=1,nxmax)
          read(43,rec=nymax+irec)(iclay(j,2),j=1,nxmax)
        endif
        if(AERTYP(7:7).eq.'1')
     &    read(45,rec=irec) ((ch_tex(j,k),k=1,ntex),j=1,nxmax)
!       print*,'   LANDUSE READ IN', ilat1,ilat2,irec
!........process slice and store in array stores
        do 20 j=1,nxmax
          rlon = j*center - center/2. - 180.
          if (xminlon.lt.-180..and.rlon.gt.0.) rlon=rlon-360.
          if (xmaxlon.gt. 180..and.rlon.lt.0.) rlon=rlon+360.
          if (rlon.ge.xminlon .and. rlon.le.xmaxlon) then
            lrec = lrec + 1
            stores(1) = rlat
            stores(2) = rlon
            char_2    = ch_topo(j)
            if (ichar(char_2(1:1))*256+ichar(char_2(2:2))-1000.lt.-200)
     &                then ! OCEAN/UNDEFINED
              stores(3) = 0.0
            else
              stores(3) = ichar(char_2(1:1))*256+ichar(char_2(2:2))-1000
            endif
            char_2    = ch_htsd(j)
            stores(4) = ichar(char_2(1:1))*256+ichar(char_2(2:2))
            if(LSMTYP.eq.'BATS') then
              do k=1,nveg
                 stores(k+4) = float(ichar(ch_cat(j,k)))
              enddo
            else if(LSMTYP.eq.'USGS') then
              do k=1,nveg
                 stores(k+4) = float(iusgs(j,k))*50./32767.+50.
              enddo
              stores(nveg+5) = float(isand(j,1))*50./32767.+50.
              stores(nveg+6) = float(isand(j,2))*50./32767.+50.
              stores(nveg+7) = float(iclay(j,1))*50./32767.+50.
              stores(nveg+8) = float(iclay(j,2))*50./32767.+50.
            endif
            if(AERTYP(7:7).eq.'1') then
              if(LSMTYP.eq.'BATS') then
                do k=1,ntex
                   stores(nveg+4+k)=float(ichar(ch_tex(j,k)))
                enddo
              else if(LSMTYP.eq.'USGS') then
                do k=1,ntex
                   stores(nveg+8+k)=float(ichar(ch_tex(j,k)))
                enddo
              endif
            endif
            write(48) stores
!add
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec)   = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec)  = stores(4)**2+stores(3)**2
            if (grdlnmn.le.-180.0.and.xobs(lrec).gt.0.0)
     &                                xobs(lrec)=xobs(lrec)-360.
            if (xobs(lrec).lt.grdlnmn) grdlnmn=xobs(lrec)
            if (yobs(lrec).lt.grdltmn) grdltmn=yobs(lrec)
!add_
          end if
   20   continue
   30 continue

      close(46)
      if(LSMTYP.eq.'BATS') then
        close(47)
      else if(LSMTYP.eq.'USGS') then
        close(41)
        close(42)
        close(43)
        if(ntypec.eq.2) close(44)
      endif
!
      print 300,lrec
  300 format(1x,i10,' terrain heights read from land use volume')
      ihmax = (xmaxlat-xminlat)/xnc
      jhmax = (xmaxlon-xminlon)/xnc
      if (ihmax.gt.iter .or. jhmax.gt.jter)
     &      print 270,iter,ihmax,jter,jhmax
  270 format(1x,'***array dimension error***',/,'     iter = ',i5,
     &   ' must be greater than ',i5,10x,'jter = ',i5,
     &   ' must be greater than ',i5)
      if (ihmax*jhmax.gt.iblk)
     &      print 271,ihmax*jhmax,iblk
  271 format(1x,'***array dimension error***',/,'  ihmax*jhmax = ',i5,
     &   ' must be greater than ',i5,10x,'iblk = ',i5)

      nobs = lrec

      rewind(48)
      return
      end
!
      SUBROUTINE RDLDTR_s
      implicit none
!     include '../../Commons/env/include/netcdf.inc'
      include 'domain.param'
!
      integer nxmax, nymax
      parameter (nxmax=21600/ntypec_s,nymax=nxmax/2)
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay
!
      integer ltype,lrec
      real*4  stores
      common/a/ stores(50),ltype
      real*4  xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      common/aa/xminlat,xminlon,xmaxlat,xmaxlon,grdltmn,grdlnmn
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
!
      integer ilat1, ilat2
      integer iflag, ierr
      character filelev*36, filcat*36, filusgs*11, filsand*7, filclay*7
      character filtext*10

      integer*2 isand(nxmax,2)
      integer*2 iclay(nxmax,2)
      integer*2 iusgs(nxmax,nveg)
      character*1 ch_tex(nxmax,ntex)
      character*1 ch_cat(nxmax,nveg)
      character*2 ch_topo(nxmax),char_2,ch_htsd(nxmax)
!
      real*4  center,rlat,rlon
      integer irec,irect,j,k
      integer ihmax,jhmax
      logical there
!
      if (ntypec_s.lt.10) then
        write (filelev,101) ntypec_s
 101    format('../DATA/SURFACE/GTOPO30_',i1,'MIN.dat')
      else if (ntypec_s.lt.100) then
        write (filelev,102) ntypec_s
 102    format('../DATA/SURFACE/GTOPO30_',i2,'MIN.dat')
      else
      write(*,*) 'For terrain, ntypec_s is not set correctly ',ntypec_s
        stop 'subroutine RDLDTR_s'
      end if
      inquire(file=filelev,exist=there)
      if(.not.there) then
        print*,'ERROR OPENING ',filelev,' FILE:  FILE DOES NOT EXIST'
        STOP '4810 IN SUBROUTINE RDLDTR_s'
      endif
      open(46,file=filelev,form='unformatted'
     &       ,recl=nxmax*ibyte/2,access='direct')
      
      if(LSMTYP.eq.'BATS') then
         if (ntypec_s.lt.10) then
           write (filcat,201) ntypec_s
 201       format('../DATA/SURFACE/GLCC',i1,'MIN_BATS.dat')
         else if (ntypec_s.lt.100) then
           write (filcat,202) ntypec_s
 202       format('../DATA/SURFACE/GLCC',i2,'MIN_BATS.dat')
         else
       write(*,*) 'For landuse, ntypec_s is not set correctly ',ntypec_s
           stop 'subroutine RDLDTR_s'
         end if
         inquire(file=filcat,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ',filcat, ' FILE:  FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR_s'
         endif
         open(47,file=filcat,form='unformatted'
     &          ,recl=nxmax*nveg*ibyte/4,access='direct')
      else if(LSMTYP.eq.'USGS') then
         if (ntypec_s.lt.10) then
           write (filusgs,301) ntypec_s
 301       format('VEG-USGS.0',i1)
           write (filsand,401) ntypec_s
 401       format('SAND.0',i1)
           write (filclay,501) ntypec_s
 501       format('CLAY.0',i1)
         else if (ntypec_s.lt.100) then
           write (filusgs,302) ntypec_s
 302       format('VEG-USGS.',i2)
           write (filsand,402) ntypec_s
 402       format('SAND.',i2)
           write (filclay,502) ntypec_s
 502       format('CLAY.',i2)
         else
           write(*,*) 'Error input of ntypec_s: ',ntypec_s
         endif
         if(ntypec_s.eq.2) then
            inquire(file='../DATA/SURFACE/'//filusgs//'A',exist=there)
            if(.not.there) then
               print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &               ,' FILE: FILE DOES NOT EXIST'
               STOP '4820 IN SUBROUTINE RDLDTR_s'
            endif
            inquire(file='../DATA/SURFACE/'//filusgs//'B',exist=there)
            if(.not.there) then
               print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &               ,' FILE: FILE DOES NOT EXIST'
               STOP '4820 IN SUBROUTINE RDLDTR_s'
            endif
         else
         inquire(file='../DATA/SURFACE/'//filusgs,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filusgs
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR_s'
         endif
         endif
         inquire(file='../DATA/SURFACE/'//filsand,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filsand
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR_s'
         endif
         inquire(file='../DATA/SURFACE/'//filclay,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filclay
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4820 IN SUBROUTINE RDLDTR_s'
         endif
         if(ntypec_s.eq.2) then
            open(41,file='../DATA/SURFACE/'//filusgs//'A'
     &      ,form='unformatted',recl=nxmax*ibyte/2,access='direct')
            open(44,file='../DATA/SURFACE/'//filusgs//'B'
     &      ,form='unformatted',recl=nxmax*ibyte/2,access='direct')
         else
            open(41,file='../DATA/SURFACE/'//filusgs,form='unformatted'
     &             ,recl=nxmax*ibyte/2,access='direct')
         endif
         open(42,file='../DATA/SURFACE/'//filsand,form='unformatted'
     &         ,recl=nxmax*ibyte/2,access='direct')
         open(43,file='../DATA/SURFACE/'//filclay,form='unformatted'
     &         ,recl=nxmax*ibyte/2,access='direct')
      endif
      if(AERTYP(7:7).eq.'1') then
         if (ntypec_s.lt.10) then
           write (filtext,303) ntypec_s
 303       format('SOILCAT.0',i1)
         else if (ntypec_s.lt.100) then
           write (filtext,304) ntypec_s
 304       format('SOILCAT.',i2)
         else
           write(*,*) 'Error input of ntypec_s: ',ntypec_s
         endif
         inquire(file='../DATA/SURFACE/'//filtext,exist=there)
         if(.not.there) then
           print*,'ERROR OPENING ','../DATA/SURFACE/'//filtext
     &           ,' FILE: FILE DOES NOT EXIST'
           STOP '4830 IN SUBROUTINE RDLDTR_s'
         endif
         open(45,file='../DATA/SURFACE/'//filtext,form='unformatted'
     &         ,recl=nxmax*ntex*ibyte/4,access='direct')
      endif

      lrec=0
      rewind(48)
      center = xnc
      ilat1 = nint((90.-xmaxlat)/xnc)
      ilat1 = max(min(nymax,ilat1),1)
      ilat2 = nint((90.-xminlat)/xnc)+1
      ilat2 = max(min(nymax,ilat2),1)
      do 30 irec=ilat1,ilat2
        rlat = 90.-center*irec+center/2.
        irect = nymax - irec + 1
        read(46,rec=irect)(ch_topo(j),j=1,nxmax)
!       print*,'   ELEVATION READ IN', ilat1,ilat2,irec
        read(46,rec=irect)(ch_htsd(j),j=1,nxmax)
!       print*,'   ELEVATION STD DEV READ IN', ilat1,ilat2,irec
        if(LSMTYP.eq.'BATS') then
          read(47,rec=irect)((ch_cat(j,k),k=1,nveg),j=1,nxmax)
        else if(LSMTYP.eq.'USGS') then
          if(ntypec_s.eq.2) then
             do k=1,13
                read(41,rec=nymax*(k-1)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
             do k=14,nveg
                read(44,rec=nymax*(k-14)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
          else
             do k=1,nveg
                read(41,rec=nymax*(k-1)+irec)(iusgs(j,k),j=1,nxmax)
             enddo
          endif
          read(42,rec=      irec)(isand(j,1),j=1,nxmax)
          read(42,rec=nymax+irec)(isand(j,2),j=1,nxmax)
          read(43,rec=      irec)(iclay(j,1),j=1,nxmax)
          read(43,rec=nymax+irec)(iclay(j,2),j=1,nxmax)
        endif
        if(AERTYP(7:7).eq.'1')
     &    read(45,rec=irec) ((ch_tex(j,k),k=1,ntex),j=1,nxmax)
!       print*,'   LANDUSE READ IN', ilat1,ilat2,irec
!........process slice and store in array stores
        do 20 j=1,nxmax
          rlon = j*center - center/2. - 180.
          if (xminlon.lt.-180..and.rlon.gt.0.) rlon=rlon-360.
          if (xmaxlon.gt. 180..and.rlon.lt.0.) rlon=rlon+360.
          if (rlon.ge.xminlon .and. rlon.le.xmaxlon) then
            lrec = lrec + 1
            stores(1) = rlat
            stores(2) = rlon
            char_2    = ch_topo(j)
            if (ichar(char_2(1:1))*256+ichar(char_2(2:2))-1000.lt.-200)
     &              then ! OCEAN/UNDEFINED
              stores(3) = 0.0
            else
              stores(3) = ichar(char_2(1:1))*256+ichar(char_2(2:2))-1000
            endif
            char_2    = ch_htsd(j)
            stores(4) = ichar(char_2(1:1))*256+ichar(char_2(2:2))
            if(LSMTYP.eq.'BATS') then
              do k=1,nveg
                 stores(k+4) = float(ichar(ch_cat(j,k)))
              enddo
            else if(LSMTYP.eq.'USGS') then
              do k=1,nveg
                 stores(k+4) = float(iusgs(j,k))*50./32767.+50.
              enddo
              stores(nveg+5) = float(isand(j,1))*50./32767.+50.
              stores(nveg+6) = float(isand(j,2))*50./32767.+50.
              stores(nveg+7) = float(iclay(j,1))*50./32767.+50.
              stores(nveg+8) = float(iclay(j,2))*50./32767.+50.
            endif
            if(AERTYP(7:7).eq.'1') then
              if(LSMTYP.eq.'BATS') then
                do k=1,ntex
                   stores(nveg+4+k)=float(ichar(ch_tex(j,k)))
                enddo
              else if(LSMTYP.eq.'USGS') then
                do k=1,ntex
                   stores(nveg+8+k)=float(ichar(ch_tex(j,k)))
                enddo
              endif
            endif
            write(48) stores
!add
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec)   = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec)  = stores(4)**2+stores(3)**2
            if (grdlnmn.le.-180.0.and.xobs(lrec).gt.0.0)
     &                                xobs(lrec)=xobs(lrec)-360.
            if (xobs(lrec).lt.grdlnmn) grdlnmn=xobs(lrec)
            if (yobs(lrec).lt.grdltmn) grdltmn=yobs(lrec)
!add_
          end if
   20   continue
   30 continue

      close(46)
      if(LSMTYP.eq.'BATS') then
        close(47)
      else if(LSMTYP.eq.'USGS') then
        close(41)
        close(42)
        close(43)
        if(ntypec.eq.2) close(44)
      endif
!
      print 300,lrec
  300 format(1x,i10,' terrain heights read from land use volume')
      ihmax = (xmaxlat-xminlat)/xnc
      jhmax = (xmaxlon-xminlon)/xnc
      if (ihmax.gt.iter .or. jhmax.gt.jter)
     &      print 270,iter,ihmax,jter,jhmax
  270 format(1x,'***array dimension error***',/,'     iter = ',i5,
     &   ' must be greater than ',i5,10x,'jter = ',i5,
     &   ' must be greater than ',i5)
      if (ihmax*jhmax.gt.iblk)
     &      print 271,ihmax*jhmax,iblk
  271 format(1x,'***array dimension error***',/,'  ihmax*jhmax = ',i5,
     &   ' must be greater than ',i5,10x,'iblk = ',i5)

      nobs = lrec

      rewind(48)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE rot2nrot(LAMS,PHIS,POLLON,POLLAT,LAM,PHI)
      implicit none
!----------------------------------------------------------------------------
! Purpose:
!     Adaption of the DWD-Functions to convert rotated pole coordinates
!     (PHIS,LAMS) into geogrphic coordinates (PHI,LAM). The location of 
!     the rotated pole is passed trough POLLON and POLLAT. POLLON and
!     POLLAT give the origin of the rotated grid. The 
!     first four arguments are input, the last two are output. All angles 
!     are in degrees (north>0, east>0)
! History:
!     05/90   D.MAJEWSKI (DWD)
!     03/93   D.BRESCH (ETHZ)
!     09/96   D.LUETHI (ETHZ)
 
! declaration of arguments:
      REAL*4  LAMS,PHIS,PPHI,PLAM,PHI,LAM,POLLON,POLLAT
 
! declaration of internal vars:
      real*4  zarg1,zarg2
      real*4  zphis,zlams,arg
 
      real*4  r2d,d2r
      real*4  ZSINPOL,ZCOSPOL,ZLAMPOL

      r2d = 45./atan(1.)
      d2r = atan(1.)/45.
 
      PLAM=POLLON + 180.
      PPHI = 90. - POLLAT

      IF (PLAM.GT.180.) PLAM = PLAM-360.
      ZSINPOL = SIN(d2r*PPHI)
      ZCOSPOL = COS(d2r*PPHI)
      ZLAMPOL =     d2r*PLAM

      zphis = 2*ATAN(EXP(d2r*phis))-atan(1.)*2.
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = d2r*ZLAMS

! first, the conversion of PHIS to PHI:
      ARG     = ZCOSPOL*COS(ZPHIS)*COS(ZLAMS) + ZSINPOL*SIN(ZPHIS)
      PHI = r2d*ASIN(ARG)
 
! follows conversion of LAMS to LAM:
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
        LAM = r2d*ATAN2(ZARG1,ZARG2)
      ENDIF

      END
!
      subroutine setup(nunit,iy,jx,ntypec,nveg,iproj,ds,clat,clon
     &         , igrads,ibyte,filout,filctl)
      implicit none
      integer nunit,iy,jx,ntypec,nveg,igrads,ibyte
      character iproj*6, filout*50, filctl*50
      real*4  ds,clat,clon
!
      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
!
      rin =  1.5         ! 1.5 rad of influence-coarse mesh

      write(6,*) 'ntypec=',ntypec
      write(6,*) 'iy=',iy
      write(6,*) 'jx=',jx
      write(6,*) 'ds=',ds
      write(6,*) 'clat=',clat
      write(6,*) 'clon=',clon
      write(6,*) 'rin=',rin
      write(6,*) 'iproj=',iproj
!
      call fexist(filout)
      open(nunit,file=filout,status='unknown'
     &    ,form='unformatted',access='direct',recl=iy*jx*ibyte)     
      if(igrads.eq.1) then
        call fexist(filctl)
        open(31,file=filctl,status='unknown')
      endif

!
      dsinm = ds*1000.
!
      nnc=nint(60./float(ntypec))
      xnc=float(ntypec)/60.
      print*,'***** Terrain resolution (min): ',xnc*60.
!
      return
      end
!
      SUBROUTINE SMTH121(htgrid,iy,jx,hscr1)
      implicit none
!
!   PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
!              2DX WAVES FROM THE FIELDS htgrid
!
      integer iy,jx
      real*4  htgrid(iy,jx),hscr1(iy,jx)
      integer i,j
!
      do j=1,jx
      do i=1,iy
         hscr1(i,j) = htgrid(i,j)
      enddo
      enddo
      do i=1,iy
      do j=2,jx-1
         if((htgrid(i,j).le.-.1).or.(htgrid(i,j).gt.0.)) then
            hscr1(i,j) = .25 * ( 2.*htgrid(i,j) + htgrid(i,j+1)
     &                                          + htgrid(i,j-1) )
         endif
      enddo
      enddo
      do j=1,jx
      do i=2,iy-1
         if((hscr1(i,j).le.-.1).or.(hscr1(i,j).gt.0.)) then
            htgrid(i,j) = .25 * (2.*hscr1(i,j)+hscr1(i+1,j)
     &                                        +hscr1(i-1,j))
         endif
      enddo
      enddo
      return
      end
      subroutine smther(slab, is1, is2, npass, point, iflg)
      implicit none
!
!        purpose: spatially smooth data in slab to dampen short
!                 wavelength components
!
      integer is1,is2,npass,iflg
      real*4  slab(is1,is2),xnu(2)
      character point*5
!
      integer icross,ie,je,iem,jem,i,j,k,kp
      real*4  asv,aplus,cell
!
      icross=0
      if(point.eq.'cross') icross=1
      ie = is1
      je = is2
      iem = ie - 5
      jem = je - 5
      xnu(1)=0.5
      xnu(2) = -.50
      do 50 k=1,npass
      do 50 kp=1,2
!   first smooth in the is1 direction
      do 1 i=1,ie
      asv=slab(i,1)
      do 1 j=2,je-1
      aplus=slab(i,j+1)
      cell=slab(i,j)
      slab(i,j)=slab(i,j)+xnu(kp)*((asv+aplus)/2.0-slab(i,j))
      if(iflg.eq.0)then
         if((i.gt.6).and.(i.lt.iem).and.(j.gt.6).and.(j.lt.jem))
     a   slab(i,j) = cell
      elseif(iflg.eq.1)then
         if(i.gt.20) slab(i,j) = cell
      endif
      asv=cell
    1 continue
!   smooth in the is2 direction
      do 2 j=1,je
      asv=slab(1,j)
      do 2 i=2,ie-1
      aplus=slab(i+1,j)
      cell=slab(i,j)
      slab(i,j)=slab(i,j)+xnu(kp)*((asv+aplus)/2.0-slab(i,j))
      if(iflg.eq.0)then
         if((i.gt.6).and.(i.lt.iem).and.(j.gt.6).and.(j.lt.jem))
     a   slab(i,j) = cell
      elseif(iflg.eq.1)then
         if(i.gt.20) slab(i,j) = cell
      endif
      asv=cell
    2 continue
   40 continue
   50 continue
      return
      end
!
      subroutine smthtr(slab1,is1,is2)
      implicit none
!
!  smooth terrain arrays
!
      integer is1,is2
      integer nocean
      parameter (nocean=20000)
      integer ii(nocean), jj(nocean)
      real*4  slab1(is1,is2)
      character point*5
      integer i,j,n,k
      integer n1,npass,iflg
!
      n=1
      do 10 i = 1,is1
      do 10 j = 1,is2
         if(slab1(i,j).lt.0.0) then
            ii(n)=i
            jj(n)=j
            slab1(i,j)=0.0
            n=n+1
         end if
   10 continue
      n1 = n - 1
      print 15, n1
   15 format(5x,'  there are a total of ',i5,' points of ocean')
      if(n.gt.nocean) print 16
   16 format(1h0,2x,'dimension exceeded in subr smthtr')
      point='cross'
      npass = 10
      iflg = 0          ! 0 = smoothing only at boundary
      call smther( slab1, is1, is2, npass, point, iflg )
!     npass = 2000
!     iflg = 0          ! 1 = extensive smoothing over south boundary
!     call smther( slab1, is1, is2, npass, point, iflg )
      do 20 i = 1, is1
      do 20 j = 1, is2
         slab1(i,j)=slab1(i,j)
         if(slab1(i,j).lt.0.0) slab1(i,j)=0.0
   20 continue
      do 25 k=1,n-1
         i=ii(k)
         j=jj(k)
         slab1(i,j)=-0.001
   25 continue
      return
      end
!
      subroutine surf(xlat,xlon,lnduse,iy,jx,iter,jter,incr,dsgrid
     &              ,lndout,land,lnd8,nrec,grdltmn,grdlnmn,h2opct
     &              ,LSMTYP,sanda,sandb,claya,clayb,frac_lnd,nveg
     &              ,AERTYP,intext,texout,frac_tex,ntex)
      implicit none
!---------------------------------------------------------------------
!
!  imx,jmx must correspond to iy,jx in the master input file;
!  otherwise the program will abort.
!
      integer iy,jx,iter,jter,nveg,ntex,nrec
      real*4  grdltmn,grdlnmn,h2opct
      CHARACTER*4 LSMTYP
      CHARACTER*7 AERTYP
      real*4  sanda(iy,jx),sandb(iy,jx),claya(iy,jx),clayb(iy,jx)
      real*4  frac_lnd(iy,jx,nveg),frac_tex(iy,jx,ntex)
!-----------------------------------------------------------------------
      real*4  stores
      integer ltype
      common/a/ stores(50),ltype
      real*4  xlat(iy,jx),xlon(iy,jx)
      integer lnduse(iy,jx),intext(iy,jx)
      real*4  lndout(iy,jx),land(iy,jx,2),itex(iy,jx,2),texout(iy,jx)
      logical flag
      integer i,j,k
      integer ii,jj,ilev,iindex,jindex,lrec,incr,lengdo,nbase
      real*4  dsgrid
      real*8  xx,yy
      real*8  lnd8(iter,jter)
      real*8  bint
      integer isystm,system
      external bint
!
      flag    = .true.
      lrec    = 0
      do 5 k=1,2
      do 5 j=1,jx
      do 5 i=1,iy
      land(i,j,k) = 0
    5 itex(i,j,k) = 0
!
!
!-----grid the data.  grdltmn=minimum latitude of incoming data.
!-----grdlnmn = minimum longitude of incoming data.  point(1,1)
!-----is value at (grdltmn,grdlnmn)
!
   60 format(1x,'*** iindex = ',i3,'   jindex = ',i3,'   lrec = ',i5,
     1 '   lat = ',f10.3,3x,'lon = ',f10.3,3x,'grdltmn = ',f10.3,5x,
     2 'grdlnmn = ',f10.3)
!
      if(AERTYP(7:7).eq.'1') then
         if(LSMTYP.eq.'BATS') then
            lengdo = nveg+ntex
         else if(LSMTYP.eq.'USGS') then
            lengdo = nveg+4+ntex
         endif
      else
         if(LSMTYP.eq.'BATS') then
            lengdo = nveg
         else if(LSMTYP.eq.'USGS') then
            lengdo = nveg+4
         endif
      endif
      do 110 ilev=1,lengdo
         rewind(48)
         do 90 lrec=1,nrec
            read(48) stores
            jindex = (stores(2) - grdlnmn)*incr + 1.1
            iindex = (stores(1) - grdltmn)*incr + 1.1
            if (iindex.gt.iter .or. jindex.gt.jter) then
               print 60,iindex,jindex,lrec,stores(1),stores(2),
     &            grdltmn,grdlnmn
               stop 60
            endif
            lnd8(iindex,jindex) = stores(ilev+4)
   90    continue
!
         if(ilev.le.nveg) then
         do 100 ii=1,iy
         do 100 jj=1,jx
            yy = - (grdltmn-xlat(ii,jj))/dsgrid + 1.0
            if (grdlnmn.le.-180.0.and.xlon(ii,jj).gt.0.0)
     &          xlon(ii,jj)=xlon(ii,jj)-360.
            xx = - (grdlnmn-xlon(ii,jj))/dsgrid + 1.0
            lndout(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
            frac_lnd(ii,jj,ilev)=lndout(ii,jj)
!
!     note: it is desirable to force grid boxes with less
!           than 75 percent water to be a land category,
!           even if water is the largest single category.
!
            if ( LSMTYP.eq.'BATS'.and.(ilev.eq.14.or.ilev.eq.15)
     &          .and. lndout(ii,jj).lt.h2opct) goto 100
            if ( LSMTYP.eq.'USGS'.and.ilev.eq.25
     &          .and. lndout(ii,jj).lt.h2opct) goto 100

            if (lndout(ii,jj).gt.land(ii,jj,1)) then
               land(ii,jj,1) = lndout(ii,jj)
               land(ii,jj,2) = ilev
            endif
  100    continue
         else if(((LSMTYP.eq.'USGS'.and.ilev.gt.nveg+4).or.
     &            (LSMTYP.eq.'BATS'.and.ilev.gt.nveg)).and.
     &           AERTYP(7:7).eq.'1') then
         if(LSMTYP.eq.'BATS') then
            nbase = nveg
         else if(LSMTYP.eq.'USGS') then
            nbase = nveg+4
         endif
         do 117 ii=1,iy
         do 117 jj=1,jx
            yy = - (grdltmn-xlat(ii,jj))/dsgrid + 1.0
            if (grdlnmn.le.-180.0.and.xlon(ii,jj).gt.0.0)
     &          xlon(ii,jj)=xlon(ii,jj)-360.
            xx = - (grdlnmn-xlon(ii,jj))/dsgrid + 1.0
            texout(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
            frac_tex(ii,jj,ilev-nbase)=texout(ii,jj)
!
!     note: it is desirable to force grid boxes with less
!           than 75 percent water to be a land category,
!           even if water is the largest single category.
!
            if ( LSMTYP.eq.'BATS'.and.ilev.eq.nveg+14
     &          .and. texout(ii,jj).lt.h2opct) goto 117
            if ( LSMTYP.eq.'USGS'.and.ilev.eq.nveg+18
     &          .and. texout(ii,jj).lt.h2opct) goto 117
            if (texout(ii,jj).gt.itex(ii,jj,1)) then
               itex(ii,jj,1) = texout(ii,jj)
               if(LSMTYP.eq.'BATS') then
                  itex(ii,jj,2) = ilev-nbase
               else if(LSMTYP.eq.'USGS') then
                  itex(ii,jj,2) = ilev-nbase
               endif
            endif
  117    continue
         else 
           if(LSMTYP.eq.'USGS') then
              do ii=1,iy
              do jj=1,jx
                 yy = - (grdltmn-xlat(ii,jj))/dsgrid + 1.0
                 if (grdlnmn.le.-180.0.and.xlon(ii,jj).gt.0.0)
     &               xlon(ii,jj)=xlon(ii,jj)-360.
                 xx = - (grdlnmn-xlon(ii,jj))/dsgrid + 1.0
                 if(ilev.eq.nveg+1) then
                    sanda(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
                 else if(ilev.eq.nveg+2) then
                    sandb(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
                 else if(ilev.eq.nveg+3) then
                    claya(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
                 else if(ilev.eq.nveg+4) then
                    clayb(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
                 endif
              enddo
              enddo
            endif
         endif
  110 continue
      close(48)
      isystm=system('/bin/rm fort.48')
!
      do i=1,iy
      do j=1,jx
         lndout(i,j) = land(i,j,2)
         lnduse(i,j) = int(land(i,j,2))
         if(AERTYP(7:7).eq.'1') then
            texout(i,j) = itex(i,j,2)
            intext(i,j) = int(itex(i,j,2))
         endif
      enddo
      enddo
!
      return
      end

      SUBROUTINE XYOBSLL(iy,jx,iproj,clat,clon,plat,plon
     &                  ,truelatL,truelatH)
      implicit none
      integer iy,jx
      character*6 iproj
      real*4  clat,clon,plat,plon,truelatL,truelatH
!
      integer iter,jter,iblk
      parameter(iter=2400,jter=2400,iblk=2880000)
      integer nobs
      real*4  xobs(iblk),yobs(iblk)
      common /block0/ xobs,yobs,nobs
      real*4  ht(iblk),htsd(iblk),ht2(iblk),sand(iblk,2),clay(iblk,2)
      common /block/ ht,htsd,ht2,sand,clay

      real*4  dsinm,rin,xn,xnc
      common /const/ dsinm,rin,xn,xnc
      integer nnc
      common /const_int/ nnc
!
      integer ilen,ie,je
      real*4  r2d,d2r,c1,psi1,pole,a,psx,cell,cell2,r,xcntr,ycntr
      real*4  cntrj,cntri,pi,xrot,phix,xlonx,flpp
      real*4  flp,xnr,ynr,xr,yr,phi1,phir,c2,phic
      integer ii,im
!
!
      ilen = iy*jx+2
      ie = iy - 1
      je = jx - 1
      r2d = 57.29578
      c1   = 1.454441e-4
      psi1 = 1.0e36
      pole = 90.
      if (iproj .eq. 'LAMCON') psi1 = 90.-TRUELATH
      if (iproj .eq. 'POLSTR') psi1 = 30.0
      psi1 = psi1/r2d
!-----psi1 is colatitude of lat where cone or plane intersects earth
      a = 6371.229
      if (clat.lt.0.)then
         if(truelatH.gt.0.) then
            psi1 = -(90.-TRUELATH)
         else
            psi1 = -(90.+TRUELATH)
         endif
         pole = -90.
         psi1 = psi1/r2d
      endif
      if (iproj.eq.'LAMCON'.or.iproj.eq.'POLSTR') then
         psx = (pole-clat)/r2d
         if (iproj.eq.'LAMCON') then
            cell  = a*sin(psi1)/xn
            cell2 = (tan(psx/2.))/(tan(psi1/2.))
         else if (iproj.eq.'POLSTR') then
            cell  = a*sin(psx)/xn
            cell2 = (1. + cos(psi1))/(1. + cos(psx))
         endif
         r = cell*(cell2)**xn
         xcntr = 0.0
         ycntr = -r
      endif
      cntrj = float(je)/2.
      cntri = float(ie)/2.
!
!-----grid incoming data.  grdltmn=minimum latitude of incoming data.
!-----grdlnmn=minimum longitude of incoming data.
!
      pi = 4.0*atan(1.0)
      r2d = 45./atan(1.0)
      d2r = atan(1.0)/45.
      a = 6371.229
      do 460 ii=1,nobs
         im = ii-1
         if (iproj.eq.'LAMCON'.or.iproj.eq.'POLSTR') then
            xrot = clon + 90./xn
            phix  = yobs(ii)
            xlonx = xobs(ii)
            flpp  = (xlonx-xrot)/r2d
            flp   = xn*flpp
            psx   = (pole-phix)/r2d
            if (iproj.eq.'LAMCON') then
               cell  = a*sin(psi1)/xn
               cell2 = (tan(psx/2.))/(tan(psi1/2.))
            else if (iproj.eq.'POLSTR') then
               cell  = a*sin(psx)/xn
               cell2 = (1.+cos(psi1))/(1.+cos(psx))
            endif
            r = cell*(cell2)**xn
            xobs(ii) = (r*cos(flp)-xcntr)*1000.
            yobs(ii) = (r*sin(flp)-ycntr)*1000.
            if (clat.lt.0.0) xobs(ii) = -xobs(ii)
         endif
         if (iproj.eq.'NORMER') then
            phi1 = 0.0 ! plat/r2d
            phir = yobs(ii)/r2d
            phic = clat/r2d
            c2 = a*cos(phi1)
            cell = cos(phir)/(1.0+sin(phir))
            cell2 = cos(phic)/(1.0+sin(phic))
            ycntr = -c2*log(cell2)
            xobs(ii) = (c2*(xobs(ii)-clon)/r2d)*1000.
            yobs(ii) = (-c2*log(cell)-ycntr)*1000.
         endif
         if (iproj.eq.'ROTMER') then
            xcntr  = plon-clon
            ycntr  = plat-clat
            xnr=xobs(ii)
            ynr=yobs(ii)
            call NROT2ROT(xnr,ynr,plon,plat,xr,yr)
            xobs(ii) = a*d2r*(xcntr+xr)*1000.
            yobs(ii) = a*d2r*(ycntr+yr)*1000.
         endif
         ht(ii)  = ht(ii)/100.
         ht2(ii) = ht2(ii)/100000.
 460  continue
      return
      end

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
