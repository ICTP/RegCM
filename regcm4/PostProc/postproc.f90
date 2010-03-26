      program postproc
      use mod_regcm_param , only : ix , jx , ixm2 , jxm2 , kx , nsg ,   &
                             &     ibyte
      use mod_batflds
      use mod_bcflds
      use mod_cheflds
      use mod_outflds
      use mod_radflds
      use mod_subflds
      use mod_point
      use mod_postproc_param
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: vmisdat = -1.E30
      integer , parameter :: iin = 10 , ndim = 3
!
! Local variables
!
      character(256) :: adate , filinfo , icbcdir , icbcheadnam ,        &
                     & inheadnam , rcmdir
      integer :: add , batproc , bcproc , brec , cheproc , icount ,     &
               & crec , i , idate , idate0 , idate1 , idate2 ,          &
               & idatenew , idateold , idatex , idday , idirect ,       &
               & idout , idy , idy0 , idy1 , idy2 , idyx , ierr , ifil ,&
               & ihr , ihr0 , ihr1 , ihr2 , ihrx , imo , imo0 , imo1 ,  &
               & imo2 , imox , iotyp , irec , iyr , iyr0 , iyr1 , iyr2 ,&
               & iyrx , j , julmid , julnc , julnc0 , julnc1 , julnc2 , &
               & julncx , l , mdate0 , nb , nday , nfiles , nr1b ,      &
               & nr1c , nr1i , nr1o , nr1r , nr1s , nr2b , nr2c , nr2i ,&
               & nr2o , nr2r , nr2s , nr3b , nr3c , nr3i , nr3o , nr3r ,&
               & nr3s , nr4b , nr4c , nr4i , nr4o , nr4r , nr4s , ns ,  &
               & ntime
      logical :: batavg , batday , batdiur , bats , bc , bcavg , bcday ,&
               & bcdiur , che , cheavg , cheday , chediur , ioall ,     &
               & ioavg , ioday , iodiur , lout , outavg , outday ,      &
               & outdiur , outhead , plv , rad , radavg , radday ,      &
               & raddiur , sub , subavg , subday , subdiur , there ,    &
               & usgs
      real(4) :: clat , clon , ds , dssb , pt , xplat , xplon
      real(4) , dimension(jxm2,ixm2) :: dlat , dlon , dmap , f , ls ,   &
                                  & xlat , xlon , xmap , zs , zssd
      real(4) , dimension(jxsg,ixsg) :: dmapsb , fsb ,                  &
                                   & lssb , xlatsb , xlonsb , xmapsb ,  &
                                   & zssb , zssdsb
      real(4) , dimension(nbat2) :: factbat , offsetbat , xmaxbat ,     &
                                  & xminbat
      real(4) , dimension(nitot) :: factbc , offsetbc , xmaxbc , xminbc
      real(4) , dimension(nctot) :: factche , offsetche , xmaxche ,     &
                                  & xminche
      real(4) , dimension(notot) :: factout , offsetout , xmaxout ,     &
                                  & xminout
      real(4) , dimension(nrtot) :: factrad , offsetrad , xmaxrad ,     &
                                  & xminrad
      real(4) , dimension(nsub2) :: factsub , offsetsub , xmaxsub ,     &
                                  & xminsub
      character(256) :: fhdbat , fhdout , fhdrad , fhdsub , filavgbat , &
                     & filavgbc , filavgche , filavgout , filavgrad ,   &
                     & filavgsub , filbat , filbc , filche , fildaybat ,&
                     & fildaybc , fildayche , fildayout , fildayrad ,   &
                     & fildaysub , fildiurbat , fildiurbc , fildiurche ,&
                     & fildiurout , fildiurrad , fildiursub , filout ,  &
                     & filrad , filsub , icbchead , icbcheadsb , inhead
      character(4) :: filext
      character(256) , dimension(nfmax) :: filrcm
      integer , dimension(ndim) :: iadm
      integer , dimension(nfmax) :: ircmext
      character(64) , dimension(nbat2) :: lnambat
      character(64) , dimension(nitot) :: lnambc
      character(64) , dimension(nctot) :: lnamche
      character(64) , dimension(notot) :: lnamout
      character(64) , dimension(nrtot) :: lnamrad
      character(64) , dimension(nsub2) :: lnamsub
      integer , dimension(nhrbat) :: nbattime
      integer , dimension(nhrbc) :: nbctime
      integer , dimension(nhrche) :: nchetime
      integer , dimension(nhrout) :: nouttime
      integer , dimension(nhrrad) :: nradtime
      integer , dimension(nhrsub) :: nsubtime
      integer :: orec , outproc , radproc , rrec , srec , subproc ,     &
               & ubc3d , ubctot , ubtot , uc3d , uctot , un1b ,         &
               & un1c , un1i , un1o , un1r , un1s , un2b , un2c , un2i ,&
               & un2o , un2r , un2s , un3b , un3c , un3i , un3o , un3r ,&
               & un3s , un4b , un4c , un4i , un4o , un4r , un4s , uo3d ,&
               & uotot , ur3d , urtot , userin , ustot
      character(256) , dimension(nfmax) :: rcmext
      real(4) , dimension(2) :: sigb
      real(4) , dimension(kx+1) :: sigf
      real(4) , dimension(kx) :: sigh , sighrev
      character(64) , dimension(nbat2) :: ubat
      character(64) , dimension(nitot) :: ubc
      character(64) , dimension(nctot) :: uche
      character(64) , dimension(notot) :: uout
      character(64) , dimension(nrtot) :: urad
      character(64) , dimension(nsub2) :: usub
      integer , dimension(nbat2) :: u_bat
      integer , dimension(nitot) :: u_bc
      integer , dimension(nctot) :: u_che
      integer , dimension(notot) :: u_out
      integer , dimension(nrtot) :: u_rad
      integer , dimension(nsub2) :: u_sub
      character(64) , dimension(nbat2) :: vnambat
      character(64) , dimension(nitot) :: vnambc
      character(64) , dimension(nctot) :: vnamche
      character(64) , dimension(notot) :: vnamout
      character(64) , dimension(nrtot) :: vnamrad
      character(64) , dimension(nsub2) :: vnamsub
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(8) :: xhr , xhr0 , xhr1 , xhr2 , xhrdy , xhrm
      real(4) , dimension(ixm2) :: xlat1d
      real(4) , dimension(ixsg) :: xlatsb1d
      real(4) , dimension(jxm2) :: xlon1d
      real(4) , dimension(jxsg) :: xlonsb1d
!
      data un1i , un2i , un3i , un4i/80 , 70 , 60 , 50/
      data un1o , un2o , un3o , un4o/80 , 70 , 60 , 50/
      data un1b , un2b , un3b , un4b/80 , 70 , 60 , 50/
      data un1s , un2s , un3s , un4s/80 , 70 , 60 , 50/
      data un1r , un2r , un3r , un4r/80 , 70 , 60 , 50/
      data un1c , un2c , un3c , un4c/80 , 70 , 60 , 50/
!
      print * , 'ENTER THE TYPE OF REGCM OUTPUT TO BE PROCESSED:'
      print * , '  ICBC (0=no; 1=yes)'
      read (*,*) bcproc
      print * , '  ATM (0=no; 1=yes)'
      read (*,*) outproc
      print * , '  SRF (0=no; 1=yes)'
      read (*,*) batproc
      print * , '  RAD (0=no; 1=yes)'
      read (*,*) radproc
      print * , '  CHE (0=no; 1=yes)'
      read (*,*) cheproc
      print * , '  SUB (0=no; 1=yes)'
      read (*,*) subproc
 
      do i = 1 , nfmax
        filrcm(i) = 'FILE NOT NAMED'
      end do
 
      open (11,file=plist,status='old')
      read (11,*) idate0
      read (11,*) idate1
      read (11,*) idate2
      read (11,*) iotyp
      read (11,*) plv
      read (11,*) usgs
      read (11,*) outhead
      read (11,*) ioall
      read (11,*) ioavg
      read (11,*) iodiur
      read (11,*) ioday
      read (11,*) nday
      read (11,*) filinfo
      read (11,*) icbcdir
      read (11,*) rcmdir
      read (11,*) icbcheadnam
      read (11,*) inheadnam
      ierr = 0
      ifil = 0
      u_bc(:) = 0
      open (13,file=ulist,status='old')
      icount = 0
      do i = 1 , 7
        read (13,*)
      end do
      do i = 1 , 10
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_bc(icount) = 1
          add = add + 1
        end if
      end do
      ubc3d = add
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 3
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_bc(icount) = 1
          add = add + 1
        end if
      end do
      ubctot = add
      do i = 1 , 13
        if ( usgs ) then
          read (13,*) userin
          icount = icount + 1
          if ( userin==1 ) then
            u_bc(icount) = 1
            add = add + 1
          end if
        else
          read (13,*)
        end if
      end do
      ubctot = add
      icount = 0
      add = 0
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 12
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_out(icount) = 1
          add = add + 1
          print * , u_out(icount) , icount , '000'
!
        end if
      end do
      uo3d = add
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 6
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_out(icount) = 1
          add = add + 1
          print * , u_out(icount) , icount
        end if
      end do
!
      uotot = add
      icount = 0
      add = 0
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 28
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_bat(icount) = 1
          add = add + 1
        end if
      end do
      ustot = add
      icount = 0
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 21
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_sub(icount) = 1
          add = add + 1
        end if
      end do
      ubtot = add
      icount = 0
      add = 0
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 4
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_rad(icount) = 1
          add = add + 1
        end if
      end do
      ur3d = add
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 9
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_rad(icount) = 1
          add = add + 1
        end if
      end do
      urtot = add
      icount = 0
      add = 0
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 13
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_che(icount) = 1
          add = add + 1
        end if
      end do
      uc3d = add
      do i = 1 , 4
        read (13,*)
      end do
      do i = 1 , 72
        read (13,*) userin
        icount = icount + 1
        if ( userin==1 ) then
          u_che(icount) = 1
          add = add + 1
        end if
      end do
      uctot = add
 
      do while ( ifil<=nfmax .and. ierr==0 )
        ifil = ifil + 1
        read (11,*,iostat=ierr) rcmext(ifil)
      end do
      if ( ierr==0 ) then
        print * , 'You have exceeded the set maximum of allowed'//      &
             &' files that can be processed (nfmax)'
        print * , '   nfmax= ' , nfmax , 'ifil= ' , ifil
        print * , '   To increase the maximum number, change'//         &
             &' nfmax in postproc2.param'
        stop 'NFMAX EXCEEDED'
      end if
      nfiles = ifil - 1
      do ifil = 1 , nfiles
        adate = rcmext(ifil)
        ircmext(ifil) = (iachar(adate(1:1))-48)*1000000000 +            &
                      & (iachar(adate(2:2))-48)*100000000 +             &
                      & (iachar(adate(3:3))-48)*10000000 +              &
                      & (iachar(adate(4:4))-48)*1000000 +               &
                      & (iachar(adate(5:5))-48)                         &
                      & *100000 + (iachar(adate(6:6))-48)               &
                      & *10000 + (iachar(adate(7:7))-48)                &
                      & *1000 + (iachar(adate(8:8))-48)                 &
                      & *100 + (iachar(adate(9:9))-48)                  &
                      & *10 + (iachar(adate(10:10))-48)
      end do
 
      inhead = trim(rcmdir)//'/'//trim(inheadnam)
      icbchead = trim(icbcdir)//'/'//trim(icbcheadnam)
      icbcheadsb = 'fort.11'
      if ( iotyp==1 .or. iotyp==2 ) then
        filext = '.nc'
      else if ( iotyp==3 ) then
        filext = '.DAT'
      else if ( iotyp==4 ) then
        filext = '.V5D'
      else
      end if
      bc = .false.
      bcavg = .false.
      bcdiur = .false.
      bcday = .false.
      lout = .false.
      outavg = .false.
      outdiur = .false.
      outday = .false.
      bats = .false.
      batavg = .false.
      batdiur = .false.
      batday = .false.
      sub = .false.
      subavg = .false.
      subdiur = .false.
      subday = .false.
      rad = .false.
      radavg = .false.
      raddiur = .false.
      radday = .false.
      che = .false.
      cheavg = .false.
      chediur = .false.
      cheday = .false.
      if ( bcproc==1 ) then
        bc = ioall
        bcavg = ioavg
        bcdiur = iodiur
        bcday = ioday
      end if
      if ( outproc==1 ) then
        lout = ioall
        outavg = ioavg
        outdiur = iodiur
        outday = ioday
      end if
      if ( batproc==1 ) then
        bats = ioall
        batavg = ioavg
        batdiur = iodiur
        batday = ioday
      end if
      if ( subproc==1 ) then
        sub = ioall
        subavg = ioavg
        subdiur = iodiur
        subday = ioday
      end if
      if ( radproc==1 ) then
        rad = ioall
        radavg = ioavg
        raddiur = iodiur
        radday = ioday
      end if
      if ( cheproc==1 ) then
        che = ioall
        cheavg = ioavg
        chediur = iodiur
        cheday = ioday
      end if
 
      if ( idate0>idate1 .or. idate0>idate2 .or. idate1>idate2 ) then
        print * , 'ERROR IN DATE SPECIFICATION'
        print * , 'IDATE0 =' , idate0 , 'IDATE1 =' , idate1 ,           &
            & 'IDATE2 =' , idate2
        stop 965
      end if
 
      print * , ' '
      print * , '**************************************************'
      print * , 'IDATE0 =' , idate0 , 'IDATE1 =' , idate1 , 'IDATE2 =' ,&
          & idate2
      print * , '**************************************************'
      print * , ' '
      if ( bcproc==1 ) then
        print * , ' '
        print * , '  PROCESS ICBC DATA'
        if ( bc ) print * , '  WRITING ALL ICBC DATA'
        if ( bcavg ) print * , '  AVERAGING ALL ICBC DATA'
        if ( bcdiur ) print * , '  AVERAGING ALL DIURNAL ICBC DATA'
        if ( bcday ) print * , '  WRITING DAILY ICBC DATA'
      end if
      if ( outproc==1 ) then
        print * , ' '
        print * , '  PROCESS ATMOS OUTPUT FROM MODEL'
        if ( lout ) print * , '  WRITING ALL OUTPUT DATA'
        if ( outavg ) print * , '  AVERAGING ALL OUTPUT DATA'
        if ( outdiur ) print * , '  AVERAGING ALL DIURNAL OUTPUT DATA'
        if ( outday ) print * , '  WRITING DAILY OUTPUT DATA'
      end if
      if ( batproc==1 ) then
        print * , ' '
        print * , '  PROCESS SURFACE OUTPUT FROM MODEL'
        if ( bats ) print * , '  WRITING ALL BATS DATA'
        if ( batavg ) print * , '  AVERAGING ALL BATS DATA'
        if ( batdiur ) print * , '  AVERAGING ALL DIURNAL BATS DATA'
        if ( batday ) print * , '  WRITING DAILY BATS DATA'
      end if
      if ( radproc==1 ) then
        print * , ' '
        print * , 'PROCESS RADATION OUTPUT FROM MODEL'
        if ( rad ) print * , '  WRITING ALL RAD DATA'
        if ( radavg ) print * , '  AVERAGING ALL RAD DATA'
        if ( raddiur ) print * , '  AVERAGING ALL DIURNAL RAD DATA'
        if ( radday ) print * , '  WRITING DAILY RAD DATA'
      end if
      if ( cheproc==1 ) then
        print * , ' '
        print * , 'PROCESS CHEMTRAC OUTPUT FROM MODEL'
        if ( che ) print * , '  WRITING ALL CHE DATA'
        if ( cheavg ) print * , '  AVERAGING ALL CHE DATA'
        if ( chediur ) print * , '  AVERAGING ALL DIURNAL CHE DATA'
        if ( cheday ) print * , '  WRITING DAILY CHE DATA'
      end if
      if ( subproc==1 ) then
        print * , ' '
        print * , '  PROCESS SURFACE OUTPUT FROM MODEL'
        if ( sub ) print * , '  WRITING ALL SUB (SUB-BATS) DATA'
        if ( subavg ) print * , '  AVERAGING ALL SUB DATA'
        if ( subdiur ) print * , '  AVERAGING ALL DIURNAL SUB DATA'
        if ( subday ) print * , '  WRITING DAILY SUB DATA'
      end if
      print * , ' '
      print * , '**************************************************'
      print * , ' '
 
      print * , idate0 , idate1 , idate2
      call julian(idate0,julnc0,iyr0,imo0,idy0,ihr0)
      xhr0 = float(julnc0)
      call julian(idate1,julnc1,iyr1,imo1,idy1,ihr1)
      xhr1 = float(julnc1)
      call julian(idate2,julnc2,iyr2,imo2,idy2,ihr2)
      xhr2 = float(julnc2)
      julmid = (julnc1+julnc2)/2
      xhrm = float(julmid-mod(julmid,24))
 
!     **** READ HEADER FILE **** C
      if ( outhead ) then
        fhdout = 'HEAD_OUT.NC'
        fhdbat = 'HEAD_BAT.NC'
        fhdsub = 'HEAD_SUB.NC'
        fhdrad = 'HEAD_RAD.NC'
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
        print * , 'HEADER READ IN'
        call param(jxm2,ixm2,kx,kx,xlat,xlon,vvarmin,vvarmax,           &
                 & xlat1d,xlon1d,iadm,ndim,plv)
        call rcrecdf(fhdout,idout,vvarmin,vvarmax,ndim,ierr)
        call writehead(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,      &
                     & vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,idout,   &
                     & xhr0,iotyp)
        call clscdf(idout,ierr)
        call param(jxm2,ixm2,1,1,xlat,xlon,vvarmin,vvarmax,             &
                &  xlat1d,xlon1d,iadm,ndim,plv)
        call rcrecdf(fhdbat,idout,vvarmin,vvarmax,ndim,ierr)
        call writebathead(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls,vvarmin,&
                        & vvarmax,xlat1d,xlon1d,iadm,ndim,idout,xhr0,   &
                        & iotyp)
        call clscdf(idout,ierr)
!SUB    CALL RCRECDF(fhdsub,idout,vvarmin,vvarmax,ndim,ierr)
!SUB    CALL WRITESUBHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls
!SUB    &     , vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,idout,xhr0
!SUB    &     , iotyp)
!SUB    CALL CLSCDF(idout,ierr)
        call param(jxm2,ixm2,kx,kx,xlat,xlon,vvarmin,vvarmax,           &
             &     xlat1d,xlon1d,iadm,ndim,plv)
        call rcrecdf(fhdrad,idout,vvarmin,vvarmax,ndim,ierr)
        call writehead(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,      &
                     & vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,idout,   &
                     & xhr0,iotyp)
        call clscdf(idout,ierr)
      end if
!     **** COMPUTE VVARMIN AND VVARMAX **** C
      close (iin)
 
 
!     *********************** C
!     ****   ICBC FILE   **** C
!     *********************** C
      if ( bc .or. bcavg .or. bcdiur .or. bcday ) then
        if ( bcday .and. (bcavg .or. bcdiur) ) then
          print * , 'MUST RUN BCDAY SEPARATE FROM BCAVG AND BCDIUR'
          stop 999
        end if
        if ( bc ) then
          filbc = 'ICBC'//trim(filinfo)//filext
          call fexist(filbc)
          print * , 'ICBC FILE: ' , filbc
        end if
        if ( bcavg ) then
          filavgbc = 'ICBC'//trim(filinfo)//'AVG'//filext
          call fexist(filavgbc)
          print * , 'ICBC AVERAGE FILE: ' , filavgbc
        end if
        if ( bcdiur ) then
          fildiurbc = 'ICBC'//trim(filinfo)//'DIUR'//filext
          call fexist(fildiurbc)
          print * , 'ICBC AVERAGE FILE: ' , fildiurbc
        end if
        if ( bcday ) then
          fildaybc = 'ICBC'//trim(filinfo)//'AVG'//filext
          call fexist(fildaybc)
          print * , 'ICBC DAILY AVERAGE FILE: ' , fildaybc
        end if
        call rdheadicbc(ix,jx,jxm2,ixm2,kx,clat,clon,ds,pt,sigf,sigh,   &
                      & sighrev,xplat,xplon,f,xmap,dmap,xlat,xlon,      &
                      & zs,zssd,ls,iin,icbchead,ibyte)
        call param(jxm2,ixm2,kx,npl,xlat,xlon,vvarmin,vvarmax,          &
                 & xlat1d,xlon1d,iadm,ndim,plv)
        ifil = 1
        irec = 0
        filrcm(ifil) = trim(icbcdir)//'/ICBC'//trim(rcmext(ifil))
        open (iin,file=filrcm(ifil),status='old',form='unformatted',    &
            & access='direct',recl=ix*jx*ibyte)
        print * , 'INPUT (ICBC) FILE: ' , filrcm(ifil)
        if ( bc ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filbc,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un1i,file=filbc,status='unknown',form='unformatted',  &
                & recl=jxm2*ixm2*ibyte,access='direct')
            nr1i = 0
          else
          end if
        end if
        if ( bcday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildaybc,idday,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un2i,file=fildaybc,status='unknown',                  &
                & form='unformatted',recl=jxm2*ixm2*ibyte,              &
                & access='direct')
            nr2i = 0
          else
          end if
        end if
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvlubc(vnambc,lnambc,ubc,xminbc,xmaxbc,factbc,offsetbc)
!       **** ZERO OUT AVERAGE ARRAYS **** C
        if ( bcavg .or. bcdiur .or. bcday ) then
          if ( .not.plv ) then
            call setconst(o3davg,vmisdat,jxm2,ixm2,kx,nbc2d,nhrbc,      &
                      &   1,jxm2,1,ixm2)
          else
            call setconst(o3davg_p,vmisdat,jxm2,ixm2,npl,nbc2d,nhrbc,   &
                      &   1,jxm2,1,ixm2)
          end if
          call setconst(o2davg,vmisdat,jxm2,ixm2,nbc2d,nhrbc,1,         &
                     &  1,jxm2,1,ixm2)
          if ( .not.plv ) then
            call setconst(o3davg,0.0,jxm2,ixm2,kx,nbc2d,nhrbc,          &
                     &    1,jxm2,1,ixm2)
          else
            call setconst(o3davg_p,0.0,jxm2,ixm2,npl,nbc2d,nhrbc,       &
                     &    1,jxm2,1,ixm2)
          end if
          call setconst(o2davg,0.0,jxm2,ixm2,nbc2d,nhrbc,1,             &
                     &  1,jxm2,1,ixm2)
          do l = 1 , nhrbc
            nbctime(l) = 0
          end do
        end if
!       **** READ IN DATA **** C
        idate = idate0
        idateold = idate
        do while ( idate<=idate2 )
          idatenew = idateold + nint(dtbc)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdicbc(idate,iin,irec,ierr)
!         **** FOR END OF FILES **** C
          if ( ierr/=0 .or. idate>idatenew .or.                         &
             & (idate>ircmext(ifil+1) .and. ifil+1<=nfiles) ) then
            print * , 'END OF FILE REACHED: ifil=' , ifil , 'ierr=' ,   &
                & ierr
 10         continue
            ifil = ifil + 1
            filrcm(ifil) = trim(icbcdir)//'/ICBC'//trim(rcmext(ifil))
            print * , '  OPENING NEW FILE: ' , filrcm(ifil)
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) go to 100
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 10
            end if
            close (iin)
            irec = 0
            open (iin,file=filrcm(ifil),status='old',form='unformatted',&
                & access='direct',recl=ix*jx*ibyte)
            idate = 0
            do while ( idate<idatenew )
              print * , 'SEARCHING FOR PROPER DAY:' , idate , idatenew
              call rdicbc(idate,iin,irec,ierr)
              if ( ierr/=0 ) stop 'READ ICBC ERROR'
              if ( idate>idatenew ) then
                print * , 'FILE ERROR (ICBC): DATE EXCEEDED'
                print * , idate , idatenew
                stop 999
              end if
            end do
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtbc)
          if ( ihr==0 ) ihr = 24/nint(dtbc)
          call calcrh(ifld2d,ifld3d,jxm2,ixm2,kx,nbc2d,nbc3d,sigh,pt,   &
                    & nti,nqvi,npsi,nrhi,ntdi,nthi,jxm2,ixm2)
          call htsig(ifld3d,ifld2d,nti,nhgti,npsi,zs,sigh,pt,jxm2,ixm2, &
                   & kx,jxm2,ixm2,nbc3d,nbc2d)
!         CALL CALCHGT(ifld2d,ifld3d,jxm2,ixm2,kx,nbc2d,nbc3d
!         &       , zs,sigf,sigh,pt,nti,nqvi,npsi,nhgti,jxm2,ixm2)
          call calcslp(ifld3d,ifld2d,nhgti,nti,npsi,zs,nslpi,sigh,      &
                     & jxm2,ixm2,kx,nbc3d,nbc2d,jxm2,ixm2)
          call calcvd(ifld3d,jxm2,ixm2,kx,nbc3d,ds,dmap,xmap,nui,nvi,   &
                    & nvori,ndivi,jxm2,ixm2)
          if ( plv ) then
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nui,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nvi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nqvi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nrhi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nvori,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlin(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,ndivi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlog(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nti,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlog(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,jxm2,ixm2,  &
                     &  kx,nthi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call intlog(ifld3d_p,ifld3d,ifld2d,npsa,pt,sigh,jxm2,ixm2,  &
                     &  kx,ntdi,plev,npl,nbc3d,nbc2d,jxm2,ixm2)
            call height(ifld3d_p,ifld3d,ifld2d,nti,npsi,zs,sigh,jxm2,   &
                      & ixm2,kx,npl,nhgti,plev,nbc3d,nbc2d,pt,jxm2,ixm2)
          end if
!         **** AVERAGE DATA **** c
          if ( bcavg .or. bcdiur .or. bcday ) then
            if ( idate>idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
              call avgdata2d(i2davg,ifld2d,jxm2,ixm2,nbc2d,nhrbc,ihr,   &
                           & vmisdat)
              if ( .not.plv ) then
                call avgdata3d(i3davg,ifld3d,jxm2,ixm2,kx,nbc3d,nhrbc,  &
                             & ihr,vmisdat)
              else
                call avgdata3d(i3davg_p,ifld3d_p,jxm2,ixm2,npl,nbc3d,   &
                             & nhrbc,ihr,vmisdat)
              end if
              nbctime(ihr) = nbctime(ihr) + 1
              if ( bcday ) then
                ntime = 0
                do l = 1 , nhrbc
                  ntime = ntime + nbctime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtbc) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtbc - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgbc(sighrev,vnambc,lnambc,ubc,xminbc,     &
                                & xmaxbc,factbc,offsetbc,vvarmin,       &
                                & vvarmax,xlat1d,xlon1d,iadm,ndim,xhrdy,&
                                & nbctime,idday,vmisdat,iotyp,un2i,nr2i,&
                                & plv,u_bc)
                  do l = 1 , nhrbc
                    nbctime(l) = 0
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE ICBC DATA IN NETCDF FORMAT AT EACH TIME STEP **** c
          if ( bc ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writebc(vvarmin,vvarmax,vnambc,lnambc,ubc,xminbc,    &
                         & xmaxbc,factbc,offsetbc,iadm,ndim,xlat1d,     &
                         & xlon1d,sighrev,vmisdat,idout,xhr,iotyp,un1i, &
                         & nr1i,plv,u_bc)
              print * , 'DATA WRITTEN: ' , xhr , idate
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtbc)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        if ( bc ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1i)
          else
          end if
          print * , 'DONE WRITING ICBC DATA!!!'
        end if
        if ( bcday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2i)
          else
          end if
          print * , 'DONE WRITING ICBC DATA!!!'
        end if
      end if
 
!     **** CLOSE INPUT (ICBC) FILE **** c
 100  continue
      close (iin)
 
!     **** WRITE ICBC AVERAGED OUTPUT FIELDS **** C
      if ( bcavg .or. bcdiur ) then
        do ihr = 1 , nhrbc
          if ( nbctime(ihr)<=0 ) then
            print * , 'Not enough data for average.'
            print * , 'nbctime must have values great than zero.'
            print * , nbctime
            stop 'BCAVG or BCDIUR'
          end if
        end do
      end if
      if ( bcavg ) then
        if ( iotyp==1 .or. iotyp==2 ) then
          call rcrecdf(filavgbc,idout,vvarmin,vvarmax,ndim,ierr)
        else if ( iotyp==3 ) then
          open (un3i,file=filavgbc,status='unknown',form='unformatted', &
              & recl=jxm2*ixm2*ibyte,access='direct')
          nr3i = 0
        else
        end if
        call writeavgbc(sighrev,vnambc,lnambc,ubc,xminbc,xmaxbc,factbc, &
                      & offsetbc,vvarmin,vvarmax,xlat1d,xlon1d,iadm,    &
                      & ndim,xhrm,nbctime,idout,vmisdat,iotyp,un3i,nr3i,&
                      & plv,u_bc)
        if ( iotyp==1 .or. iotyp==2 ) then
          call clscdf(idout,ierr)
        else if ( iotyp==3 ) then
          close (un3i)
        else
        end if
        print * , 'DONE WRITING AVERAGED DATA IN COARDS CONVENTIONS!!!'
      end if
      if ( bcdiur ) then
        if ( iotyp==1 .or. iotyp==2 ) then
          call rcrecdf(fildiurbc,idout,vvarmin,vvarmax,ndim,ierr)
        else if ( iotyp==3 ) then
          open (un4i,file=fildiurbc,status='unknown',form='unformatted',&
              & recl=jxm2*ixm2*ibyte,access='direct')
          nr4i = 0
        else
        end if
        call writediurbc(sighrev,vnambc,lnambc,ubc,xminbc,xmaxbc,factbc,&
                       & offsetbc,vvarmin,vvarmax,xlat1d,xlon1d,iadm,   &
                       & ndim,xhrm,nbctime,idout,vmisdat,iotyp,un4i,    &
                       & nr4i,plv,u_bc)
        if ( iotyp==1 .or. iotyp==2 ) then
          call clscdf(idout,ierr)
        else if ( iotyp==3 ) then
          close (un4i)
        else
        end if
        print * , 'DONE WRITING AVERAGED DIURNAL DATA!!!'
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
 
 
!     ********************** C
!     ****   OUT FILE   **** C
!     ********************** C
      if ( lout .or. outavg .or. outdiur .or. outday ) then
        if ( outday .and. (outavg .or. outdiur) ) then
          print * , 'MUST RUN OUTDAY SEPARATE FROM OUTAVG AND OUTDIUR'
          stop 999
        end if
        if ( lout ) then
          filout = 'ATM'//trim(filinfo)//filext
          call fexist(filout)
          print * , 'OUTPUT FILE: ' , filout
        end if
        if ( outavg ) then
          filavgout = 'ATM'//trim(filinfo)//'AVG'//filext
          call fexist(filavgout)
          print * , 'OUTPUT AVERAGE FILE: ' , filavgout
        end if
        if ( outdiur ) then
          fildiurout = 'ATM'//trim(filinfo)//'DIUR'//filext
          call fexist(fildiurout)
          print * , 'OUTPUT AVERAGE FILE: ' , fildiurout
        end if
        if ( outday ) then
          fildayout = 'ATM'//trim(filinfo)//'AVG'//filext
          call fexist(fildayout)
          print * , 'OUTPUT DAILY AVERAGE FILE: ' , fildayout
        end if
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
        call param(jxm2,ixm2,kx,npl,xlat,xlon,vvarmin,vvarmax,          &
                &  xlat1d,xlon1d,iadm,ndim,plv)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/ATM.'//trim(rcmext(ifil))
        print * , 'INPUT (OUT) FILE: ' , filrcm(ifil)
        if ( idirect==1 ) then
          orec = 0
          open (iin,file=filrcm(ifil),status='old',form='unformatted',  &
              & recl=ixm2*jxm2*ibyte,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
        if ( lout ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filout,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un1o,file=filout,status='unknown',form='unformatted', &
                & recl=jxm2*ixm2*ibyte,access='direct')
            nr1o = 0
          else
          end if
        end if
        if ( outday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildayout,idday,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un2o,file=fildayout,status='unknown',                 &
                 & form='unformatted',recl=jxm2*ixm2*ibyte,             &
                 & access='direct')
            nr2o = 0
          else
          end if
        end if
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvluout(vnamout,lnamout,uout,xminout,xmaxout,factout,     &
                    & offsetout)
!       **** ZERO OUT AVERAGE ARRAYS **** C
        if ( outavg .or. outdiur .or. outday ) then
          if ( .not.plv ) then
            call setconst(o3davg,vmisdat,jxm2,ixm2,kx,nout2d,nhrout,    &
                        & 1,jxm2,1,ixm2)
          else
            call setconst(o3davg_p,vmisdat,jxm2,ixm2,npl,nout2d,nhrout, &
                        & 1,jxm2,1,ixm2)
          end if
          call setconst(o2davg,vmisdat,jxm2,ixm2,nout2d,nhrout,1,       &
                     &  1,jxm2,1,ixm2)
          if ( .not.plv ) then
            call setconst(o3davg,0.0,jxm2,ixm2,kx,nout2d,nhrout,        &
                     &    1,jxm2,1,ixm2)
          else
            call setconst(o3davg,0.0,jxm2,ixm2,npl,nout2d,nhrout,       &
                     &    1,jxm2,1,ixm2)
          end if
          call setconst(o2davg,0.0,jxm2,ixm2,nout2d,nhrout,1,           &
                     &  1,jxm2,1,ixm2)
          do l = 1 , nhrout
            nouttime(l) = 0
          end do
        end if
!       **** READ IN DATA **** C
        idate = idate0
        if ( mdate0==idate0 ) then
          idate = idate0
        else
          idate = idate0 + nint(dtout)
        end if
        idateold = idate
        do while ( idate<=idate2 )
          idatenew = idateold + nint(dtout)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdatm(idate,iin,orec,idirect,ierr)
!         **** FOR END OF FILES **** C
          if ( ierr/=0 .or. idate>idatenew .or.                         &
             & (idate>ircmext(ifil+1) .and. ifil+1<=nfiles) ) then
            print * , 'END OF FILE REACHED: ifil=' , ifil , 'ierr=' ,   &
                & ierr
 110        continue
            ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/ATM.'//trim(rcmext(ifil))
            print * , '  OPENING NEW FILE: ' , filrcm(ifil)
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) go to 200
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 110
            end if
            close (iin)
            idate = idatenew
            if ( idirect==1 ) then
              orec = 0
              open (iin,file=filrcm(ifil),status='old',                 &
                   & form='unformatted',recl=ixm2*jxm2*ibyte,           &
                   & access='direct')
            else
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted')
            end if
            call rdatm(idate,iin,orec,idirect,ierr)
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtout)
          if ( ihr==0 ) ihr = 24/nint(dtout)
          call calcrh(ofld2d,ofld3d,jxm2,ixm2,kx,nout2d,nout3d,sigh,pt, &
                    & nta,nqva,npsa,nrh,ntda,ntha,jxm2,ixm2)
          call htsig(ofld3d,ofld2d,nta,nhgt,npsa,zs,sigh,pt,jxm2,ixm2,  &
                   & kx,jxm2,ixm2,nout3d,nout2d)
!         CALL CALCHGT(ifld2d,ifld3d,jxm2,ixm2,kx,nbc2d,nbc3d
!         &       , zs,sigf,sigh,pt,nti,nqvi,npsi,nhgti,jxm2,ixm2)
          call calcslp(ofld3d,ofld2d,nhgt,nta,npsa,zs,nslp,sigh,jxm2,   &
                     & ixm2,kx,nout3d,nout2d,jxm2,ixm2)
          call calcvd(ofld3d,jxm2,ixm2,kx,nout3d,ds,dmap,xmap,nua,nva,  &
                    & nvora,ndiva,jxm2,ixm2)
          if ( plv ) then
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nua,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nva,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nqva,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nrh,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nvora,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlin(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,ndivi,plev,npl,nout3d,nout2d,jxm2,ixm2)
!           print*, ifld3d(:,:,:,nti)
            call intlog(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,nta,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlog(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,ntha,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call intlog(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,jxm2,ixm2,  &
                      & kx,ntda,plev,npl,nout3d,nout2d,jxm2,ixm2)
            call height(ofld3d_p,ofld3d,ofld2d,nti,npsa,zs,sigh,jxm2,   &
                      & ixm2,kx,npl,nhgt,plev,nout3d,nout2d,pt,         &
                      & jxm2,ixm2)
          end if
!         **** AVERAGE DATA **** c
          if ( outavg .or. outdiur .or. outday ) then
            if ( idate>idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
              call avgdata2d(o2davg,ofld2d,jxm2,ixm2,nout2d,nhrout,ihr, &
                           & vmisdat)
              if ( .not.plv ) then
                call avgdata3d(o3davg,ofld3d,jxm2,ixm2,kx,nout3d,nhrout,&
                            &  ihr,vmisdat)
              else
                call avgdata3d(o3davg_p,ofld3d_p,jxm2,ixm2,npl,nout3d,  &
                             & nhrout,ihr,vmisdat)
              end if
              nouttime(ihr) = nouttime(ihr) + 1
              if ( outday ) then
                ntime = 0
                do l = 1 , nhrout
                  ntime = ntime + nouttime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtout) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtout - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgout(sighrev,vnamout,lnamout,uout,xminout,&
                                 & xmaxout,factout,offsetout,vvarmin,   &
                                 & vvarmax,xlat1d,xlon1d,iadm,ndim,     &
                                 & xhrdy,nouttime,idday,vmisdat,iotyp,  &
                                 & un2o,nr2o,plv,u_out)
                  do l = 1 , nhrout
                    nouttime(l) = 0
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE OUT DATA IN NETCDF FORMAT AT EACH TIME STEP **** c
          if ( lout ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writeout(vvarmin,vvarmax,vnamout,lnamout,uout,       &
                          & xminout,xmaxout,factout,offsetout,iadm,ndim,&
                          & xlat1d,xlon1d,sighrev,vmisdat,idout,xhr,    &
                          & iotyp,un1o,nr1o,plv,u_out)
              print * , 'DATA WRITTEN: ' , xhr , idate
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtout)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        if ( lout ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1o)
          else
          end if
          print * , 'DONE WRITING OUT DATA!!!'
        end if
        if ( outday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2o)
          else
          end if
          print * , 'DONE WRITING OUT DATA!!!'
        end if
      end if
 
!     **** CLOSE INPUT (OUT) FILE **** c
 200  continue
      close (iin)
 
!     **** WRITE OUT AVERAGED OUTPUT FIELDS **** C
      if ( outavg .or. outdiur ) then
        do ihr = 1 , nhrout
          if ( nouttime(ihr)<=0 ) then
            print * , 'Not enough data for average.'
            print * , 'nouttime must have values great than zero.'
            print * , nouttime
            stop 'OUTAVG-OUTDIUR'
          end if
        end do
      end if
      if ( outavg ) then
        if ( iotyp==1 .or. iotyp==2 ) then
          call rcrecdf(filavgout,idout,vvarmin,vvarmax,ndim,ierr)
        else if ( iotyp==3 ) then
          open (un3o,file=filavgout,status='unknown',form='unformatted',&
              & recl=jxm2*ixm2*ibyte,access='direct')
          nr3o = 0
        else
        end if
        call writeavgout(sighrev,vnamout,lnamout,uout,xminout,xmaxout,  &
                       & factout,offsetout,vvarmin,vvarmax,xlat1d,      &
                       & xlon1d,iadm,ndim,xhrm,nouttime,idout,vmisdat,  &
                       & iotyp,un3o,nr3o,plv,u_out)
        if ( iotyp==1 .or. iotyp==2 ) then
          call clscdf(idout,ierr)
        else if ( iotyp==3 ) then
          close (un3o)
        else
        end if
        print * , 'DONE WRITING AVERAGED DATA IN COARDS CONVENTIONS!!!'
      end if
      if ( outdiur ) then
        if ( iotyp==1 .or. iotyp==2 ) then
          call rcrecdf(fildiurout,idout,vvarmin,vvarmax,ndim,ierr)
        else if ( iotyp==3 ) then
          open (un4o,file=fildiurout,status='unknown',                  &
              & form='unformatted',recl=jxm2*ixm2*ibyte,access='direct')
          nr4o = 0
        else
        end if
        call writediurout(sighrev,vnamout,lnamout,uout,xminout,xmaxout, &
                        & factout,offsetout,vvarmin,vvarmax,xlat1d,     &
                        & xlon1d,iadm,ndim,xhrm,nouttime,idout,vmisdat, &
                        & iotyp,un4o,nr4o,plv,u_out)
        if ( iotyp==1 .or. iotyp==2 ) then
          call clscdf(idout,ierr)
        else if ( iotyp==3 ) then
          close (un4o)
        else
        end if
        print * , 'DONE WRITING AVERAGED DIURNAL DATA!!!'
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
!     ************************* C
!     ****    BATS FILE    **** C
!     ************************* C
      if ( bats .or. batavg .or. batdiur .or. batday ) then
        if ( batday .and. (batavg .or. batdiur) ) then
          print * , 'MUST RUN BATDAY SEPARATE FROM BATAVG AND BATDIUR'
          stop 999
        end if
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
!       **** OPEN NetCDF FILE **** C
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/SRF.'//trim(rcmext(ifil))
        print * , 'INPUT (BATS) FILE: ' , filrcm(ifil)
        if ( idirect==1 ) then
          brec = 0
          open (iin,file=filrcm(ifil),status='old',form='unformatted',  &
              & recl=ixm2*jxm2*ibyte,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvlubat(vnambat,lnambat,ubat,xminbat,xmaxbat,factbat,     &
                    & offsetbat)
!       **** COMPUTE VVARMIN AND VVARMAX **** C
        call param(jxm2,ixm2,1,1,xlat,xlon,vvarmin,vvarmax,             &
              &    xlat1d,xlon1d,iadm,ndim,plv)
        vvarmax(3) = 1050.
        iadm(3) = 1
        call setconst(sigb,1.0,2,1,1,1,1,2,1,1,1)
        if ( bats ) then
          filbat = 'SRF'//trim(filinfo)//filext
          call fexist(filbat)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filbat,idout,vvarmin,vvarmax,ndim,ierr)
            print * , 'OPENING BATS NetCDF FILE' , idout
          else if ( iotyp==3 ) then
            open (un1b,file=filbat,status='unknown',form='unformatted', &
                & recl=jxm2*ixm2*ibyte,access='direct')
            nr1b = 0
            print * , 'OPENING BATS GrADS FILE' , un1b
          else
          end if
        end if
        if ( batavg .or. batdiur .or. batday ) then
          if ( batavg ) then
            filavgbat = 'SRF'//trim(filinfo)//'AVG'//filext
            call fexist(filavgbat)
          end if
          if ( batdiur ) then
            fildiurbat = 'SRF'//trim(filinfo)//'DIUR'//filext
            call fexist(fildiurbat)
          end if
          if ( batday ) then
            fildaybat = 'SRF'//trim(filinfo)//'AVG'//filext
            call fexist(fildaybat)
            if ( iotyp==1 .or. iotyp==2 ) then
              call rcrecdf(fildaybat,idday,vvarmin,vvarmax,ndim,ierr)
            else if ( iotyp==3 ) then
              open (un2b,file=fildaybat,status='unknown',               &
                   & form='unformatted',recl=jxm2*ixm2*ibyte,           &
                   & access='direct')
              nr2b = 0
            else
            end if
          end if
          call setconst(b2davg,0.0,jxm2,ixm2,nbat2,nhrbat,1,            &
                    &   1,jxm2,1,ixm2)
          do l = 1 , nhrbat
            nbattime(l) = 0
          end do
        end if
        if ( mdate0==idate0 ) then
          idate = idate0
        else
          idate = idate0 + nint(dtbat)
        end if
        idateold = idate
!       **** Initialize Max and Min Temperatures
!       do j=1,ixm2
!       do i=1,jxm2
!       bfld2d(i,j,ntmax) = bfld2d(i,j,ntanm)
!       bfld2d(i,j,ntmin) = bfld2d(i,j,ntanm)
!       end do
!       end do
        print * , idate , idate0 , idate1 , idate2
        do while ( idate<=idate2 )
!         **** READ BATS FILE **** C
          idatenew = idateold + nint(dtbat)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdsrf(idate,iin,brec,idirect,ierr)
          if ( ierr/=0 .or. idate>idatenew .or.                         &
             & (idate>ircmext(ifil+1) .and. ifil+1<=nfiles) ) then
            print * , 'READ ERROR OR END OF FILE REACHED' , ierr
 210        continue
            ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/SRF.'//trim(rcmext(ifil))
            print * , '  OPENING NEW FILE: ' , filrcm(ifil)
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) exit
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 210
            end if
            close (iin)
            idate = idatenew
            if ( idirect==1 ) then
              brec = 0
              open (iin,file=filrcm(ifil),status='old',                 &
                   & form='unformatted',recl=jxm2*ixm2*ibyte,           &
                   & access='direct')
            else
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted')
            end if
            call rdsrf(idate,iin,brec,idirect,ierr)
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtbat) + 1
          if ( ihr==0 ) ihr = 24/nint(dtbat)
!         CALL TMINMAX(bfld2d,b2davg,jxm2,ixm2,nbat2,nhrbat,ihr
!         &       , ntanm,ntmax,ntmin)
!         CALL CALCMSE2D(bfld2d,zs,jxm2,ixm2,nbat2,ntanm,nqanm,nmsea)
          call calcrh2d(bfld2d,jxm2,ixm2,nbat2,ntanm,nqanm,npsrf,nrha,  &
                      & vmisdat)
!         **** AVERAGE DATA **** c
          if ( batavg .or. batdiur .or. batday ) then
            if ( idate>idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
              call avgdatabat(ihr,vmisdat)
              nbattime(ihr) = nbattime(ihr) + 1
              if ( batday ) then
                ntime = 0
                do l = 1 , nhrbat
                  ntime = ntime + nbattime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtbat) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtbat - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgbat(vmisdat,vnambat,lnambat,ubat,xminbat,&
                                 & xmaxbat,factbat,offsetbat,vvarmin,   &
                                 & vvarmax,xlat1d,xlon1d,iadm,ndim,     &
                                 & xhrdy,nbattime,idday,iotyp,un2b,nr2b,&
                                 & u_bat)
                  do l = 1 , nhrbat
                    nbattime(l) = 0
                    do nb = 1 , nbat2
                      do j = 1 , ixm2
                        do i = 1 , jxm2
                          b2davg(i,j,nb,l) = 0.0
                        end do
                      end do
                    end do
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE BATS IN COARDS NetCDF CONVENTIONS **** C
          if ( bats ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writebat(vnambat,lnambat,ubat,xminbat,xmaxbat,       &
                          & factbat,offsetbat,vvarmin,vvarmax,xlat1d,   &
                          & xlon1d,iadm,ndim,vmisdat,xhr,idout,         &
                          & iotyp,un1b,nr1b,u_bat)
              print * , 'BATS DATA WRITTEN:  ' , idate
              print * , ''
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtbat)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        close (iin)
        if ( bats ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1b)
          else
          end if
          print * , 'DONE WRITING BATS DATA!!!'
        end if
        if ( batday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2b)
          else
          end if
          print * , 'DONE WRITING DAILY BATS DATA!!!'
        end if
        if ( batavg .or. batdiur ) then
          do ihr = 1 , nhrbat
            if ( nbattime(ihr)<=0 ) then
              print * , 'Not enough data for average.'
              print * , 'nbattime must have values great than zero.'
              print * , nbattime
              stop 'BATAVG-BATDIUR'
            end if
          end do
        end if
        if ( batavg ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filavgbat,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un3b,file=filavgbat,status='unknown',                 &
                 & form='unformatted',recl=jxm2*ixm2*ibyte,             &
                 & access='direct')
            nr3b = 0
          else
          end if
          call writeavgbat(vmisdat,vnambat,lnambat,ubat,xminbat,xmaxbat,&
                         & factbat,offsetbat,vvarmin,vvarmax,xlat1d,    &
                         & xlon1d,iadm,ndim,xhrm,nbattime,idout,iotyp,  &
                         & un3b,nr3b,u_bat)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un3b)
          else
          end if
          print * , 'DONE WRITING AVERAGED BATS DATA'
        end if
        if ( batdiur ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildiurbat,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un4b,file=fildiurbat,status='unknown',                &
                 & form='unformatted',recl=jxm2*ixm2*ibyte,             &
                 & access='direct')
            nr4b = 0
          else
          end if
          call writediurbat(vmisdat,vnambat,lnambat,ubat,xminbat,       &
                          & xmaxbat,factbat,offsetbat,vvarmin,vvarmax,  &
                          & xlat1d,xlon1d,iadm,ndim,xhrm,nbattime,idout,&
                          & iotyp,un4b,nr4b,u_bat)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un4b)
          else
          end if
          print * , 'DONE WRITING AVERAGED DIURNAL BATS DATA'
        end if
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
 
!     ****************************** C
!     ****    RADIATION FILE    **** C
!     ****************************** C
      if ( rad .or. radavg .or. raddiur .or. radday ) then
        if ( radday .and. (radavg .or. raddiur) ) then
          print * , 'MUST RUN RADDAY SEPARATE FROM RADAVG AND RADDIUR'
          stop 999
        end if
!       **** OPEN NetCDF FILE **** C
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/RAD.'//trim(rcmext(ifil))
        if ( idirect==1 ) then
          rrec = 0
          open (iin,file=filrcm(ifil),status='old',form='unformatted',  &
              & recl=jxm2*ixm2*ibyte,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
        print * , 'INPUT (RAD) FILE: ' , filrcm(ifil)
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvlurad(vnamrad,lnamrad,urad,xminrad,xmaxrad,factrad,     &
                    & offsetrad)
!       **** COMPUTE VVARMIN AND VVARMAX **** C
        call param(jxm2,ixm2,kx,npl,xlat,xlon,vvarmin,vvarmax,          &
               &   xlat1d,xlon1d,iadm,ndim,plv)
        iadm(3) = kx
        if ( rad ) then
          filrad = 'RAD'//trim(filinfo)//filext
          call fexist(filrad)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filrad,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un1r,file=filrad,status='unknown',form='unformatted', &
                & recl=jxm2*ixm2*ibyte,access='direct')
            nr1r = 0
          else
          end if
        end if
        if ( radday ) then
          fildayrad = 'RAD'//trim(filinfo)//'AVG'//filext
          call fexist(fildayrad)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildayrad,idday,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un2r,file=fildayrad,status='unknown',                 &
                 & form='unformatted',recl=jxm2*ixm2*ibyte,             &
                 & access='direct')
            nr2r = 0
          else
          end if
        end if
        if ( radavg .or. raddiur .or. radday ) then
          if ( radavg ) then
            filavgrad = 'RAD'//trim(filinfo)//'AVG'//filext
            call fexist(filavgrad)
          end if
          if ( raddiur ) then
            fildiurrad = 'RAD'//trim(filinfo)//'DIUR'//filext
            call fexist(fildiurrad)
          end if
          if ( .not.plv ) then
            call setconst(r3davg,vmisdat,jxm2,ixm2,kx,nr2d,nhrrad,      &
                       &  1,jxm2,1,ixm2)
          else
            call setconst(r3davg_p,vmisdat,jxm2,ixm2,npl,nr2d,nhrrad,   &
                       &  1,jxm2,1,ixm2)
          end if
          call setconst(r2davg,vmisdat,jxm2,ixm2,nr2d,nhrrad,1,         &
                      & 1,jxm2,1,ixm2)
          if ( .not.plv ) then
            call setconst(r3davg,0.0,jxm2,ixm2,kx,nr2d,nhrrad,          &
                      &   1,jxm2,1,ixm2)
          else
            call setconst(r3davg_p,0.0,jxm2,ixm2,npl,nr2d,nhrrad,       &
                      &   1,jxm2,1,ixm2)
          end if
          call setconst(r2davg,0.0,jxm2,ixm2,nr2d,nhrrad,1,             &
                      & 1,jxm2,1,ixm2)
          do l = 1 , nhrrad
            nradtime(l) = 0
          end do
        end if
        if ( mdate0==idate0 ) then
          idate = idate0
        else
          idate = idate0 + nint(dtrad)
        end if
        idateold = idate
        do while ( idate<=idate2 )
!         **** READ RADIATION FILE **** C
          idatenew = idateold + nint(dtrad)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdrad(iin,idate,rrec,idirect,ierr)
          if ( ierr/=0 .or. idate>idatenew .or.                         &
             & (idate>ircmext(ifil+1) .and. ifil+1<=nfiles) ) then
            print * , 'END OF FILE REACHED: ifil=' , ifil , 'ierr=' ,   &
                & ierr
 220        continue
            ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/RAD.'//trim(rcmext(ifil))
            print * , '  OPENING NEW FILE: ' , filrcm(ifil)
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) exit
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 220
            end if
            close (iin)
            idate = idatenew
            if ( idirect==1 ) then
              rrec = 0
              open (iin,file=filrcm(ifil),status='old',                 &
                   & form='unformatted',recl=jxm2*ixm2*ibyte,           &
                   & access='direct')
            else
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted')
            end if
            call rdrad(iin,idate,rrec,idirect,ierr)
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtrad)
          if ( ihr==0 ) ihr = 24/nint(dtrad)
          if ( plv ) then
            call intlin(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,jxm2,ixm2, &
                      & kx,ncld,plev,npl,nr3d,nr2d,jxm2,ixm2)
            call intlin(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,jxm2,ixm2, &
                      & kx,nclwp,plev,npl,nr3d,nr2d,jxm2,ixm2)
            call intlin(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,jxm2,ixm2, &
                      & kx,nqrs,plev,npl,nr3d,nr2d,jxm2,ixm2)
            call intlin(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,jxm2,ixm2, &
                      & kx,nqrl,plev,npl,nr3d,nr2d,jxm2,ixm2)
          end if
 
!         **** AVERAGE DATA **** c
          if ( radavg .or. raddiur .or. radday ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
              call avgdata2d(r2davg,rfld2d,jxm2,ixm2,nr2d,nhrout,ihr,   &
                           & vmisdat)
              if ( .not.plv ) then
                call avgdata3d(r3davg,rfld3d,jxm2,ixm2,kx,nr3d,nhrout,  &
                             & ihr,vmisdat)
              else
                call avgdata3d(r3davg_p,rfld3d_p,jxm2,ixm2,npl,nr3d,    &
                             & nhrout,ihr,vmisdat)
              end if
              nradtime(ihr) = nradtime(ihr) + 1
              if ( radday ) then
                ntime = 0
                do l = 1 , nhrrad
                  ntime = ntime + nradtime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtrad) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtrad - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgrad(xhrdy,sighrev,vnamrad,lnamrad,urad,  &
                                 & xminrad,xmaxrad,factrad,offsetrad,   &
                                 & vvarmin,vvarmax,xlat1d,xlon1d,iadm,  &
                                 & ndim,vmisdat,nradtime,idday,iotyp,   &
                                 & un2r,nr2r,plv,u_rad)
                  print * , 'DAILY RAD DATA WRITTEN: ' , xhr , idate
                  do l = 1 , nhrrad
                    nradtime(l) = 0
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE RADIATION IN IVE NetCDF CONVENTIONS **** C
          if ( rad ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writerad(vnamrad,lnamrad,urad,xminrad,xmaxrad,       &
                          & factrad,offsetrad,vvarmin,vvarmax,xlat1d,   &
                          & xlon1d,iadm,ndim,sighrev,vmisdat,idout,xhr, &
                          & iotyp,un1r,nr1r)
              print * , 'RADIATION DATA WRITTEN:  ' , idate
              print * , ''
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtrad)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        close (iin)
        if ( rad ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1r)
          else
          end if
          print * , 'DONE RADIATION DATA!!!'
        end if
        if ( radday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2r)
          else
          end if
          print * , 'DONE DAILY RADIATION DATA!!!'
        end if
        if ( radavg .or. raddiur ) then
          do ihr = 1 , nhrrad
            print * , nradtime , 'is nradtime'
            if ( nradtime(ihr)<=0 ) then
              print * , 'Not enough data for average.'
              print * , 'nradtime must have values great than zero.'
              print * , nradtime
              stop 'RADAVG-RADDIUR'
            end if
          end do
        end if
        if ( radavg ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filavgrad,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un3r,file=filavgrad,status='unknown',                 &
                & form='unformatted',recl=jxm2*ixm2*ibyte,              &
                & access='direct')
            nr3r = 0
          else
          end if
          call writeavgrad(xhrm,sighrev,vnamrad,lnamrad,urad,xminrad,   &
                         & xmaxrad,factrad,offsetrad,vvarmin,vvarmax,   &
                         & xlat1d,xlon1d,iadm,ndim,vmisdat,nradtime,    &
                         & idout,iotyp,un3r,nr3r,plv,u_rad)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un3r)
          else
          end if
          print * , 'DONE WRITING AVERAGED RADIATION DATA'
        end if
        if ( raddiur ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildiurrad,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un4r,file=filavgrad,status='unknown',                 &
                & form='unformatted',recl=jxm2*ixm2*ibyte,              &
                & access='direct')
            nr4r = 0
          else
          end if
          call writediurrad(xhrm,sighrev,vnamrad,lnamrad,urad,xminrad,  &
                          & xmaxrad,factrad,offsetrad,vvarmin,vvarmax,  &
                          & xlat1d,xlon1d,iadm,ndim,vmisdat,nradtime,   &
                          & idout,iotyp,un4r,nr4r,plv,u_rad)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un4r)
          else
          end if
          print * , 'DONE WRITING AVERAGED DIURNAL RADIATION DATA'
        end if
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
 
!     ****************************** C
!     ****    CHEM-TRACER FILE  **** C
!     ****************************** C
      if ( che .or. cheavg .or. chediur .or. cheday ) then
        if ( cheday .and. (cheavg .or. chediur) ) then
          print * , 'MUST RUN RADDAY SEPARATE FROM CHEAVG AND CHEDIUR'
          stop 999
        end if
!       **** OPEN NetCDF FILE **** C
        print * , 'toto' , trim(rcmext(ifil))
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/CHE.'//trim(rcmext(ifil))
        print * , 'INPUT CHEM FILE: ' , filrcm
        if ( idirect==1 ) then
          crec = 0
          open (iin,file=filrcm(ifil),status='old',form='unformatted',  &
              & recl=ixm2*jxm2*ibyte,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
        print * , 'INPUT (CHE) FILE: ' , filrcm(ifil)
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvluche(vnamche,lnamche,uche,xminche,xmaxche,factche,     &
                    & offsetche)
        print * , vnamche
 
!       **** COMPUTE VVARMIN AND VVARMAX **** C
        call param(jxm2,ixm2,kx,npl,xlat,xlon,vvarmin,vvarmax,          &
               &   xlat1d,xlon1d,iadm,ndim,plv)
!       iadm(3) = kx
        if ( che ) then
          filche = 'CHE'//trim(filinfo)//filext
          print * , filche
          call fexist(filche)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filche,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un1c,file=filche,status='unknown',form='unformatted', &
                & recl=jxm2*ixm2*ibyte,access='direct')
            nr1c = 0
          else
          end if
        end if
        if ( cheday ) then
          fildayche = 'CHE'//trim(filinfo)//'AVG'//filext
          call fexist(fildayche)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildayche,idday,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un2c,file=fildayche,status='unknown',                 &
                 &form='unformatted',recl=jxm2*ixm2*4,access='direct')
            nr2c = 0
          else
          end if
        end if
        if ( cheavg .or. chediur .or. cheday ) then
          if ( cheavg ) then
            filavgche = 'CHE'//trim(filinfo)//'AVG'//filext
            call fexist(filavgche)
          end if
          if ( chediur ) then
            fildiurche = 'CHE'//trim(filinfo)//'DIUR'//filext
            call fexist(fildiurche)
          end if
          if ( .not.plv ) then
            call setconst(c3davg,vmisdat,jxm2,ixm2,kx,nc3d,nhrche,      &
                       &  1,jxm2,1,ixm2)
          else
            call setconst(c3davg_p,vmisdat,jxm2,ixm2,npl,nc3d,nhrche,   &
                       &  1,jxm2,1,ixm2)
          end if
          call setconst(c2davg,vmisdat,jxm2,ixm2,nc2d,nhrche,1,         &
                      & 1,jxm2,1,ixm2)
          if ( .not.plv ) then
            call setconst(c3davg,0.0,jxm2,ixm2,kx,nc3d,nhrche,          &
                      &   1,jxm2,1,ixm2)
          else
            call setconst(c3davg_p,0.0,jxm2,ixm2,npl,nc3d,nhrche,       &
                      &   1,jxm2,1,ixm2)
          end if
          call setconst(c2davg,0.0,jxm2,ixm2,nc2d,nhrche,1,             &
                      & 1,jxm2,1,ixm2)
          do l = 1 , nhrche
            nchetime(l) = 0
          end do
        end if
        if ( mdate0==idate0 ) then
          idate = idate0
        else
          idate = idate0 + nint(dtche)
        end if
        idateold = idate
        print * , idate , idate0 , idate1 , idate2
        do while ( idate<=idate2 )
!         **** READ CHEM-TRACER FILE **** C
          idatenew = idateold + nint(dtche)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdche(iin,idate,crec,idirect,ierr)
          if ( idate>idatenew ) then
            print * , 'END OF FILE REACHED: ifil=' , ifil , 'ierr=' ,   &
                & ierr
 230        continue
            ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/CHE.'//trim(rcmext(ifil))
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) exit
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 230
            end if
            close (iin)
            idate = idatenew
            if ( idirect==1 ) then
              crec = 0
              open (iin,file=filrcm(ifil),status='old',                 &
                  & form='unformatted',recl=ixm2*jxm2*ibyte,            &
                  & access='direct')
            else
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted')
            end if
            call rdche(iin,idate,crec,idirect,ierr)
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
!sr       maybe add one here?
          ihr = ihr/nint(dtche)
          if ( ihr==0 ) ihr = 24/nint(dtche)
!         **** AVERAGE DATA **** c
          if ( cheavg .or. chediur .or. cheday ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
!sr           nhrcche
              call avgdata2d(c2davg,cfld2d,jxm2,ixm2,nc2d,nhrche,ihr,   &
                           & vmisdat)
              if ( .not.plv ) then
                call avgdata3d(c3davg,cfld3d,jxm2,ixm2,kx,nc3d,nhrche,  &
                             & ihr,vmisdat)
              else
                call avgdata3d(c3davg_p,cfld3d_p,jxm2,ixm2,npl,nc3d,    &
                            &  nhrche,ihr,vmisdat)
              end if
              nchetime(ihr) = nchetime(ihr) + 1
              if ( cheday ) then
                ntime = 0
                do l = 1 , nhrche
                  ntime = ntime + nchetime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtche) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtche - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgche(xhrdy,sighrev,vnamche,lnamche,uche,  &
                                 & xminche,xmaxche,factche,offsetche,   &
                                 & vvarmin,vvarmax,xlat1d,xlon1d,iadm,  &
                                 & ndim,vmisdat,nchetime,idday,iotyp,   &
                                 & un2c,nr2c,plv,u_che)
                  print * , 'DAILY CHE DATA WRITTEN: ' , xhr , idate
                  do l = 1 , nhrche
                    nchetime(l) = 0
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE RADIATION IN IVE NetCDF CONVENTIONS **** C
          if ( che ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writeche(vnamche,lnamche,uche,xminche,xmaxche,       &
                          & factche,offsetche,vvarmin,vvarmax,xlat1d,   &
                          & xlon1d,iadm,ndim,sighrev,vmisdat,idout,xhr, &
                          & iotyp,un1c,nr1c,u_che)
              print * , 'CHEM-TRACER DATA WRITTEN:  ' , idate
              print * , ''
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtche)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        close (iin)
        if ( che ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1c)
          else
          end if
          print * , 'DONE CHEM-TRACER DATA!!!'
        end if
        if ( cheday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2c)
          else
          end if
          print * , 'DONE DAILY CHEM-TRACER DATA!!!'
        end if
        if ( cheavg .or. chediur ) then
          do ihr = 1 , nhrche
            print * , nchetime , ' is nchetime'
            print * , ihr
            if ( nchetime(ihr)<=0 ) then
              print * , 'Not enough data for average.'
              print * , 'nchetime must have values great than zero.'
              print * , cheavg
              print * , chediur
              print * , 'nchetime is' , nchetime
              stop 'CHEAVG-CHEDIUR'
            end if
          end do
        end if
        if ( cheavg ) then
          print * , 'YES' , filavgche
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filavgche,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un3c,file=filavgche,status='unknown',                 &
                & form='unformatted',recl=jxm2*ixm2*ibyte,              &
                & access='direct')
            nr3c = 0
          else
          end if
          call writeavgche(xhrm,sighrev,vnamche,lnamche,uche,xminche,   &
                         & xmaxche,factche,offsetche,vvarmin,vvarmax,   &
                         & xlat1d,xlon1d,iadm,ndim,vmisdat,nchetime,    &
                         & idout,iotyp,un3c,nr3c,plv,u_che)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un3c)
          else
          end if
          print * , 'DONE WRITING AVERAGED CHE-TRACER DATA'
        end if
        if ( chediur ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildiurche,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un4c,file=filavgche,status='unknown',                 &
                 &form='unformatted',recl=jxm2*ixm2*4,access='direct')
            nr4c = 0
          else
          end if
          call writediurche(xhrm,sighrev,vnamche,lnamche,uche,xminche,  &
                          & xmaxche,factche,offsetche,vvarmin,vvarmax,  &
                          & xlat1d,xlon1d,iadm,ndim,vmisdat,nchetime,   &
                          & idout,iotyp,un4c,nr4c,plv,u_che)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un4c)
          else
          end if
          print * , 'DONE WRITING AVERAGED DIURNAL CHE-TRACER DATA'
        end if
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
 
!     ************************* C
!     ****    SUB FILE     **** C
!     ************************* C
      if ( sub .or. subavg .or. subdiur .or. subday ) then
        if ( subday .and. (subavg .or. subdiur) ) then
          print * , 'MUST RUN SUBDAY SEPARATE FROM SUBAVG AND SUBDIUR'
          stop 999
        end if
        call rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,    &
                  & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,iin,          &
                  & inhead,idirect)
        call rdheadicbc(jx*nsg,ix*nsg,jxsg,ixsg,kx,clat,clon,dssb,pt,   &
                      & sigf,sigh,sighrev,xplat,xplon,fsb,xmapsb,dmapsb,&
                      & xlatsb,xlonsb,zssb,zssdsb,lssb,mdate0,iin,      &
                      & icbcheadsb,ibyte)
        call param(jxsg,ixsg,1,1,xlatsb,xlonsb,vvarmin,vvarmax,         &
               &   xlatsb1d,xlonsb1d,iadm,ndim,plv)
        print * , '            '
        print * , '            '
        print * , vvarmin
        print * , vvarmax
        print * , dssb
        print * , '            '
        print * , '            '
!       **** OPEN NetCDF FILE **** C
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/SUB.'//trim(rcmext(ifil))
        print * , 'INPUT (SUB) FILE: ' , filrcm(ifil)
        if ( idirect==1 ) then
          srec = 0
          open (iin,file=filrcm(ifil),status='old',form='unformatted',  &
              & recl=jxsg*ixsg*ibyte,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
!       **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        call mmvlusub(vnamsub,lnamsub,usub,xminsub,xmaxsub,factsub,     &
                    & offsetsub)
!       **** COMPUTE VVARMIN AND VVARMAX **** C
        vvarmax(3) = 1050.
        iadm(3) = 1
        call setconst(sigb,1.0,2,1,1,1,1,2,1,1,1)
        if ( sub ) then
          filsub = 'SUB'//trim(filinfo)//filext
          call fexist(filsub)
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filsub,idout,vvarmin,vvarmax,ndim,ierr)
            print * , 'OPENING SUB NetCDF FILE' , idout
          else if ( iotyp==3 ) then
            open (un1s,file=filsub,status='unknown',form='unformatted', &
                & recl=jxsg*ixsg*ibyte,access='direct')
            nr1s = 0
            print * , 'OPENING SUB GrADS FILE' , un1s
          else
          end if
        end if
        if ( subavg .or. subdiur .or. subday ) then
          if ( subavg ) then
            filavgsub = 'SUB'//trim(filinfo)//'AVG'//filext
            call fexist(filavgsub)
          end if
          if ( subdiur ) then
            fildiursub = 'SUB'//trim(filinfo)//'DIUR'//filext
            call fexist(fildiursub)
          end if
          if ( subday ) then
            fildaysub = 'SUB'//trim(filinfo)//'AVG'//filext
            call fexist(fildaysub)
            if ( iotyp==1 .or. iotyp==2 ) then
              call rcrecdf(fildaysub,idday,vvarmin,vvarmax,ndim,ierr)
            else if ( iotyp==3 ) then
              open (un2s,file=fildaysub,status='unknown',               &
                   &form='unformatted',recl=jxsg*ixsg*ibyte,            &
                   &access='direct')
              nr2s = 0
            else
            end if
          end if
          call setconst(s2davg,0.0,jxsg,ixsg,nsub2,nhrsub,1,1,jxsg,1,   &
                      & ixsg)
          do l = 1 , nhrsub
            nsubtime(l) = 0
          end do
        end if
        if ( mdate0==idate0 ) then
          idate = idate0
        else
          idate = idate0 + nint(dtsub)
        end if
        idateold = idate
!       **** Initialize Max and Min Temperatures
!       do j=1,ixsg
!       do i=1,jxsg
!       sfld2d(i,j,nstmax) = sfld2d(i,j,nstanm)
!       sfld2d(i,j,nstmin) = sfld2d(i,j,nstanm)
!       end do
!       end do
        print * , idate , idate0 , idate1 , idate2
        do while ( idate<=idate2 )
!         **** READ SUB FILE **** C
          idatenew = idateold + nint(dtsub)
          call julian(idatenew,julnc,iyr,imo,idy,ihr)
          call rdsub(idate,iin,srec,idirect,ierr)
          if ( ierr/=0 .or. idate>idatenew .or.                         &
             & (idate>ircmext(ifil+1) .and. ifil+1<=nfiles) ) then
            print * , 'READ ERROR OR END OF FILE REACHED' , ierr
            print * , idate , idatenew , idateold , ircmext(ifil+1)
 240        continue
            ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/SUB.'//trim(rcmext(ifil))
            call fexistnew(filrcm(ifil),there)
            if ( .not.there ) exit
            if ( idateold/=ircmext(ifil) ) then
              print * , ' '
              print * , 'INCONSISTENT DATES: ifil=' , ifil
              print * , '  idateold=' , idateold , 'rcmext=' ,          &
                  & ircmext(ifil)
              go to 240
            end if
            close (iin)
            idate = idatenew
            if ( idirect==1 ) then
              srec = 0
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted',recl=jxsg*ixsg*ibyte,            &
                   &access='direct')
            else
              open (iin,file=filrcm(ifil),status='old',                 &
                   &form='unformatted')
            end if
            call rdsub(idate,iin,srec,idirect,ierr)
          end if
          call julian(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtsub) + 1
          if ( ihr==0 ) ihr = 24/nint(dtsub)
!         CALL TMINMAX(sfld2d,s2davg,jxsg,ixsg,nsub2,nhrsub,ihr
!         &       , nstanm,nstmax,nstmin)
!         CALL CALCMSE2D(sfld2d,zssb,jxsg,ixsg,nsub2
!         &       , nstanm,nsqanm,nsmsea)
          call calcrh2d(sfld2d,jxsg,ixsg,nsub2,nstanm,nsqanm,nspsrf,    &
                      & nsrha,vmisdat)
!         **** AVERAGE DATA **** c
          if ( subavg .or. subdiur .or. subday ) then
            if ( idate>idate1 .and. idate<=idate2 ) then
              print * , 'AVERAGING DATA: ' , idate , xhr , ihr
              call avgdatasub(ihr,vmisdat)
              nsubtime(ihr) = nsubtime(ihr) + 1
              if ( subday ) then
                ntime = 0
                do l = 1 , nhrsub
                  ntime = ntime + nsubtime(l)
                end do
                call julian(idateold,julncx,iyrx,imox,idyx,ihrx)
                if ( (nday==-1 .and. imo/=imox) .or.                    &
                   & (ntime>=nint(nday*24./dtsub) .and. nday>0) ) then
                  idatex = iyrx*1000000 + imox*10000 + 1500
                  if ( nday==-1 ) then
                    call julian(idatex,julncx,iyrx,imox,idyx,ihrx)
                    xhrdy = float(julncx)
                  else
                    xhrdy = float(julncx) + dtsub - 24. - (nday-1)*24.
                  end if
                  print * , 'WRITING CONTINUAL OUTPUT' , xhrdy , idatex
                  call writeavgsub(vmisdat,vnamsub,lnamsub,usub,xminsub,&
                                 & xmaxsub,factsub,offsetsub,vvarmin,   &
                                 & vvarmax,xlatsb1d,xlonsb1d,iadm,ndim, &
                                 & xhrdy,nsubtime,idday,iotyp,un2s,nr2s)
                  do l = 1 , nhrsub
                    nsubtime(l) = 0
                    do ns = 1 , nsub2
                      do j = 1 , ixsg
                        do i = 1 , jxsg
                          s2davg(i,j,ns,l) = 0.0
                        end do
                      end do
                    end do
                  end do
                end if
              end if
            end if
          end if
!         **** WRITE SUB IN COARDS NetCDF CONVENTIONS **** C
          if ( sub ) then
            if ( idate>=idate1 .and. idate<=idate2 ) then
              call writesub(vnamsub,lnamsub,usub,xminsub,xmaxsub,       &
                          & factsub,offsetsub,vvarmin,vvarmax,xlatsb1d, &
                          & xlonsb1d,iadm,ndim,vmisdat,xhr,idout,       &
                          & iotyp,un1s,nr1s)
              print * , 'SUB DATA WRITTEN:  ' , idate
              print * , ''
            end if
          end if
!         **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtsub)
          call julian(idate,julnc,iyr,imo,idy,ihr)
        end do
        close (iin)
        if ( sub ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un1s)
          else
          end if
          print * , 'DONE WRITING SUB DATA!!!'
        end if
        if ( subday ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idday,ierr)
          else if ( iotyp==3 ) then
            close (un2s)
          else
          end if
          print * , 'DONE WRITING DAILY SUB DATA!!!'
        end if
        if ( subavg .or. subdiur ) then
          do ihr = 1 , nhrsub
            if ( nsubtime(ihr)<=0 ) then
              print * , 'Not enough data for average.'
              print * , 'nsubtime must have values great than zero.'
              print * , nsubtime
              stop 'SUBAVG-SUBDIUR'
            end if
          end do
        end if
        if ( subavg ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(filavgsub,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un3s,file=filavgsub,status='unknown',                 &
                 &form='unformatted',recl=jxsg*ixsg*ibyte,              &
                 &access='direct')
            nr3s = 0
          else
          end if
          call writeavgsub(vmisdat,vnamsub,lnamsub,usub,xminsub,xmaxsub,&
                         & factsub,offsetsub,vvarmin,vvarmax,xlatsb1d,  &
                         & xlonsb1d,iadm,ndim,xhrm,nsubtime,idout,iotyp,&
                         & un3s,nr3s)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un3s)
          else
          end if
          print * , 'DONE WRITING AVERAGED SUB DATA'
        end if
        if ( subdiur ) then
          if ( iotyp==1 .or. iotyp==2 ) then
            call rcrecdf(fildiursub,idout,vvarmin,vvarmax,ndim,ierr)
          else if ( iotyp==3 ) then
            open (un4s,file=fildiursub,status='unknown',                &
                 &form='unformatted',recl=jxsg*ixsg*ibyte,              &
                 &access='direct')
            nr4s = 0
          else
          end if
          call writediursub(vmisdat,vnamsub,lnamsub,usub,xminsub,       &
                          & xmaxsub,factsub,offsetsub,vvarmin,vvarmax,  &
                          & xlatsb1d,xlonsb1d,iadm,ndim,xhrm,nsubtime,  &
                          & idout,iotyp,un4s,nr4s)
          if ( iotyp==1 .or. iotyp==2 ) then
            call clscdf(idout,ierr)
          else if ( iotyp==3 ) then
            close (un4s)
          else
          end if
          print * , 'DONE WRITING AVERAGED DIURNAL SUB DATA'
        end if
        print * , '***** if you see <variable not found>, do NOT worry.'
      end if
 
      print * , 'DONE!!!'
 
      end program postproc
