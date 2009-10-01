      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      logical there
      character icbcdir*50, rcmdir*50, rcmext(nfmax)*50, filinfo*50
     &        , filext*4, filrcm(nfmax)*70
     &        , filbc*70, filout*70, filbat*70
     &        , filsub*70, filrad*70, filche*70
     &        , inheadnam*50, inhead*70, icbcheadnam*50, icbchead*70
     &        , icbcheadsb*70
     &        , fhdout*70, fhdbat*70, fhdsub*70, fhdrad*70
     &        , filavgbc*70, filavgout*70, filavgbat*70
     &        , filavgsub*70, filavgrad*70, filavgche*70
     &        , fildiurbc*70, fildiurout*70, fildiurbat*70
     &        , fildiursub*70,fildiurrad*70,fildiurche*70
     &        , fildaybc*70, fildayout*70, fildaybat*70
     &        , fildaysub*70, fildayrad*70,fildayche*70
     &        , vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13
     &        , vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
     &        , vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
     &        , vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
     &        , vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
     &        , vnamche(nctot)*10, lnamche(nctot)*20, uche(nctot)*13
     &        , adate*50
      real xminbc(nitot), xmaxbc(nitot), factbc(nitot)
     &   , offsetbc(nitot)
     &   , xminout(notot), xmaxout(notot), factout(notot)
     &   , offsetout(notot)
     &   , xminbat(nbat2), xmaxbat(nbat2), factbat(nbat2)
     &   , offsetbat(nbat2)
     &   , xminsub(nsub2), xmaxsub(nsub2), factsub(nsub2)
     &   , offsetsub(nsub2)
     &   , xminrad(nrtot), xmaxrad(nrtot), factrad(nrtot)
     &   , offsetrad(nrtot)
     &   , xminche(nctot), xmaxche(nctot), factche(nctot)
     &   , offsetche(nctot)

      real vmisdat,clat,clon,ds,pt,xplat,xplon,dssb
      integer bcproc, outproc, batproc, subproc, radproc, cheproc
     &      , iotyp, idirect, irec, orec, brec, srec, rrec, crec
     &      , nr1i, nr2i, nr3i, nr4i
     &      , nr1o, nr2o, nr3o, nr4o
     &      , nr1b, nr2b, nr3b, nr4b
     &      , nr1s, nr2s, nr3s, nr4s
     &      , nr1r, nr2r, nr3r, nr4r
     &      , nr1c, nr2c, nr3c, nr4c  
      parameter (vmisdat=-1.E30)

      integer iin,ndim,idout,idday,nfiles,ircmext(nfmax)
     &      , un1i,un2i,un3i,un4i,un1o,un2o,un3o,un4o
     &      , un1b,un2b,un3b,un4b,un1s,un2s,un3s,un4s
     &      , un1r,un2r,un3r,un4r,un1c,un2c,un3c,un4c
     &      , iyr,iyr0,iyr1,iyr2,iyrx,imo,imo0,imo1,imo2,imox
     &      , idy,idy0,idy1,idy2,idyx,ihr,ihr0,ihr1,ihr2,ihrx
     &      , mdate0,idate,idate0,idate1,idate2,idatex,idateold,idatenew
     &      , julnc,julnc0,julnc1,julnc2,julncx,julmid
     &      , ierr,i,j,ifil,l,nb,ns
      data un1i, un2i, un3i, un4i
     &    /  80,   70,   60,   50 /
      data un1o, un2o, un3o, un4o
     &    /  80,   70,   60,   50 /
      data un1b, un2b, un3b, un4b
     &    /  80,   70,   60,   50 /
      data un1s, un2s, un3s, un4s
     &    /  80,   70,   60,   50 /
      data un1r, un2r, un3r, un4r
     &    /  80,   70,   60,   50 /
      data un1c, un2c, un3c, un4c
     &    /  80,   70,   60,   50 /
      parameter(iin=10,ndim=3)

      real ifld2d, ifld3d, i2davg, i3davg
      common /bcflds/ ifld2d(nx,ny,nbc2d), ifld3d(nx,ny,nz,nbc3d)
     &   , i2davg(nx,ny,nbc2d,nhrbc), i3davg(nx,ny,nz,nbc3d,nhrbc)
      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)
      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)
      real ifld3d_p, i3davg_p
      common /bcflds_p/ ifld3d_p(nx,ny,npl,nbc3d)
     &   , i3davg_p(nx,ny,npl,nbc3d,nhrbc)
      real ofld3d_p, o3davg_p
      common /outflds_p/ ofld3d_p(nx,ny,npl,nout3d)
     &   , o3davg_p(nx,ny,npl,nout3d,nhrout)
	real rfld3d_p, r3davg_p
      common /radflds_p/ rfld3d_p(nx,ny,npl,nr3d)
     &   , r3davg_p(nx,ny,npl,nr3d,nhrrad)
      real cfld3d_p, c3davg_p
      common /cheflds_p/ cfld3d_p(nx,ny,npl,nc3d)
     &   , c3davg_p(nx,ny,npl,nc3d,nhrche)


      real f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny)
     &    , xlon(nx,ny), dlat(nx,ny), dlon(nx,ny), zs(nx,ny)
     &    , zssd(nx,ny), ls(nx,ny)
      real fsb(nxsb,nysb), xmapsb(nxsb,nysb), dmapsb(nxsb,nysb)
     &   , xlatsb(nxsb,nysb), xlonsb(nxsb,nysb), dlatsb(nxsb,nysb)
     &   , dlonsb(nxsb,nysb), zssb(nxsb,nysb), zssdsb(nxsb,nysb)
     &   , lssb(nxsb,nysb), xlatsb1d(nysb), xlonsb1d(nxsb)
      real sigh(nz), sighrev(nz), sigf(nz+1), sigb(2)
     &   , vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)

      real*8 xhr, xhr0, xhr1, xhr2, xhrm, xhrdy
      integer idim(ndim), nday, ntime
     &      , nbctime(nhrbc), nouttime(nhrout), nbattime(nhrbat)
     &      , nsubtime(nhrsub), nradtime(nhrrad), nchetime(nhrche)

      integer   nui,  nvi,  nqvi, nrhi, nti, ntdi, nthi, nvori
     &      , ndivi,nhgti, npsi, ntgi, nslpi
      common /ipoint3d/  nui,  nvi,  nqvi, nrhi,nti, ntdi, nthi
     &      , nvori, ndivi, nhgti
      data               nui,  nvi,  nti, nqvi, nrhi,  ntdi, nthi
     &            /       1,    2,    3,    4,    5,    6,    7 /
	data               nvori,  ndivi,  nhgti
     &            /       8,    9,    10    /

      common /ipoint2d/ npsi, ntgi, nslpi
      data              npsi, ntgi, nslpi
     &            /       1,    2,  3/

      integer   nua,  nva, nomega,  nta, nqva, nqca,   nrh, nhgt
     &      , ntha,ntda,nvora,ndiva, npsa,  nrt,ntgb, nsmt, nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca,  nrh, nhgt
     &                  , ntha, ntda, nvora, ndiva
      data               nua,  nva,  nomega, nta, nqva, nqca, nrh, nhgt
     &            /       1,    2,    3,    4,    5,    6,    7,    8 /
      data               ntha, ntda, nvora, ndiva
     &            /       9,    10,    11,    12 /
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf, nslp
      data              npsa,  nrt, ntgb, nsmt,  nbf, nslp
     &            /       1,    2,    3,    4,    5,   6 /

      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      data              nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &            /       1,     2,     3,     4,     5,     6,     7 /
      data             nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &            /       8,     9,    10,    11,    12,    13,    14 /
      data             nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &            /      15,    16,    17,    18,    19,    20,    21 /
      data            ntgmax,  ntgmin,  ntamax, ntamin, w10max, psmin
     &            /      22,    23,    24,    25     , 26,      27/
      data	      nrha /28/

      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      data              nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &            /        1,      2,      3,      4,      5,      6 /
      data            nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &            /        7,      8,      9,     10,     11,     12 /
      data            nssnow,   nssh,  nsprc, nspsrf, nsrha
     &            /       13,     14,     15,     16,     17 /

     

      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      data               ncld,  nclwp,   nqrs,   nqrl
     &            /       1,      2,      3,      4 /
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp
      data               nfsw,   nflw, nclrst, nclrss, nclrlt
     &            /       1,      2,      3,      4,      5 /
      data             nclrls, nsolin, nsabtp, nfirtp
     &            /       6,      7,      8,      9 /
      integer ubctot, uotot, ustot, ubtot, urtot, uctot
	integer count,userin, add
	integer ubc3d, uo3d, us3d, ub3d, ur3d, uc3d
	integer u_bc(nitot),u_out(notot),u_rad(nrtot),u_che(nctot)
     &         , u_sub(nsub2), u_bat(nbat2)

      logical ioall, ioavg, iodiur, ioday, outhead
     &      , bc, bcavg, bcdiur, bcday, plv ,usgs
     &      , out, outavg, outdiur, outday
     &      , bats, batavg, batdiur, batday
     &      , sub, subavg, subdiur, subday
     &      , rad, radavg, raddiur, radday
     &      , che, cheavg, chediur, cheday
      print*,'ENTER THE TYPE OF REGCM OUTPUT TO BE PROCESSED:'
      print*,'  ICBC (0=no; 1=yes)'
      read(*,*) bcproc
      print*,'  ATM (0=no; 1=yes)'
      read(*,*) outproc
      print*,'  SRF (0=no; 1=yes)'
      read(*,*) batproc
      print*,'  RAD (0=no; 1=yes)'
      read(*,*) radproc
      print*,'  CHE (0=no; 1=yes)'
      read(*,*) cheproc
      print*,'  SUB (0=no; 1=yes)'
      read(*,*) subproc

      do i=1,nfmax
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
      u_bc(:)=0
	open (13,file=ulist,status='old')
	count = 0
	do i=1, 7
	read (13,*)
	enddo
	do i=1,10
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_bc(count) = 1
	  add = add + 1
	endif
	enddo
	ubc3d = add
	do i=1,4
	read (13,*)
	enddo
	do i=1,3
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_bc(count) = 1
	  add = add + 1
	endif
	enddo
	ubctot = add
	do i =1, 13
        if (USGS) then
        read(13,*) userin
	  count = count + 1
		if(userin.eq.1) then
		u_bc(count)= 1
		add = add + 1
        endif
        else
        read(13,*)
        endif
	enddo
	ubctot=add
	count = 0
	add = 0
	do i=1,4
	read (13,*)
	enddo
	do i=1,12
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_out(count) = 1
	  add = add + 1
	print*, u_out(count), count, '000'
	
	endif
	enddo
	uo3d = add
	do i=1,4
	read (13,*)
	enddo
	do i=1,6
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_out(count) = 1
	  add = add + 1
	print*, u_out(count), count
	endif
	enddo
	
	uotot = add
	count = 0
      	add = 0
	do i=1,4
	read (13,*)
	enddo
	do i=1,28
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_bat(count) = 1
	  add = add + 1
	endif
	enddo
	ustot = add
	count = 0
	do i=1,4
	read (13,*)
	enddo
	do i=1,21
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_sub(count) = 1
	  add = add + 1
	endif
	enddo
	ubtot = add
	count = 0
	add = 0
	do i=1,4
	read (13,*)
	enddo
	do i=1,4
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_rad(count) = 1
	  add = add + 1
	endif
	enddo
	ur3d = add
        do i=1,4
        read (13,*)
        enddo
	do i=1,9
	read (13,*) userin
	count = count + 1
        if (userin.eq.1) then
           u_rad(count)=1
	   add = add + 1
        endif
	enddo
        urtot=add
	count = 0
	add = 0
	do i=1,4
	read (13,*)
	enddo
	do i=1,13
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_che(count) = 1
	  add = add + 1
	endif
	enddo
	uc3d = add
	do i=1,4
	read (13,*)
	enddo
	do i=1,72
	read (13,*) userin
	count = count + 1
	if (userin.eq.1) then
	  u_che(count) = 1
	  add = add + 1
	endif
	enddo
	uctot = add

      do while (ifil.le.nfmax .and. ierr.eq.0)
        ifil = ifil + 1
        read (11,*,iostat=ierr) rcmext(ifil)
      end do
      if (ierr.eq.0) then
        print*,'You have exceeded the set maximum of allowed'
     &       //' files that can be processed (nfmax)'
        print*,'   nfmax= ',nfmax,'ifil= ',ifil
        print*,'   To increase the maximum number, change'
     &       //' nfmax in postproc2.param'
        stop 'NFMAX EXCEEDED'
      end if
      nfiles = ifil - 1
      do ifil=1,nfiles
        adate=rcmext(ifil)
        ircmext(ifil) = (iachar(adate(1:1))-48)*1000000000
     &                + (iachar(adate(2:2))-48)*100000000
     &                + (iachar(adate(3:3))-48)*10000000
     &                + (iachar(adate(4:4))-48)*1000000
     &                + (iachar(adate(5:5))-48)*100000
     &                + (iachar(adate(6:6))-48)*10000
     &                + (iachar(adate(7:7))-48)*1000
     &                + (iachar(adate(8:8))-48)*100
     &                + (iachar(adate(9:9))-48)*10
     &                + (iachar(adate(10:10))-48)
      end do

      inhead = trim(rcmdir)//'/'//trim(inheadnam)
      icbchead = trim(icbcdir)//'/'//trim(icbcheadnam)
      icbcheadsb = 'fort.11'
      if (iotyp.eq.1 .or. iotyp.eq.2) then
        filext = '.nc'
      else if (iotyp.eq.3) then
        filext = '.DAT'
      else if (iotyp.eq.4) then
        filext = '.V5D'
      end if
      bc = .false.
      bcavg = .false.
      bcdiur = .false.
      bcday = .false.
      out = .false.
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
      if (bcproc.eq.1) then
        bc = ioall
        bcavg = ioavg
        bcdiur = iodiur
        bcday = ioday
      end if
      if (outproc.eq.1) then
        out = ioall
        outavg = ioavg
        outdiur = iodiur
        outday = ioday
      end if
      if (batproc.eq.1) then
        bats = ioall
        batavg = ioavg
        batdiur = iodiur
        batday = ioday
      end if
      if (subproc.eq.1) then
        sub = ioall
        subavg = ioavg
        subdiur = iodiur
        subday = ioday
      end if
      if (radproc.eq.1) then
        rad = ioall
        radavg = ioavg
        raddiur = iodiur
        radday = ioday
      end if
      if(cheproc.eq.1) then
        che = ioall
        cheavg = ioavg
        chediur = iodiur
        cheday = ioday
      end if

      if (idate0.gt.idate1 .or. idate0.gt.idate2 .or.
     &    idate1.gt.idate2) then
        print*,'ERROR IN DATE SPECIFICATION'
        print*,'IDATE0 =',idate0,'IDATE1 =',idate1,'IDATE2 =',idate2
        stop 965
      end if

      print*,' '
      print*,'**************************************************'
      print*,'IDATE0 =',idate0,'IDATE1 =',idate1,'IDATE2 =',idate2
      print*,'**************************************************'
      print*,' '
      if (bcproc.eq.1) then
        print*,' '
        print*,'  PROCESS ICBC DATA'
        if (bc) then
          print*,'  WRITING ALL ICBC DATA'
        end if
        if (bcavg) then
          print*,'  AVERAGING ALL ICBC DATA'
        end if
        if (bcdiur) then
          print*,'  AVERAGING ALL DIURNAL ICBC DATA'
        end if
        if (bcday) then
          print*,'  WRITING DAILY ICBC DATA'
        end if
      end if
      if (outproc.eq.1) then
        print*,' '
        print*,'  PROCESS ATMOS OUTPUT FROM MODEL'
        if (out) then
          print*,'  WRITING ALL OUTPUT DATA'
        end if
        if (outavg) then
          print*,'  AVERAGING ALL OUTPUT DATA'
        end if
        if (outdiur) then
          print*,'  AVERAGING ALL DIURNAL OUTPUT DATA'
        end if
        if (outday) then
          print*,'  WRITING DAILY OUTPUT DATA'
        end if
      end if
      if (batproc.eq.1) then
        print*,' '
        print*,'  PROCESS SURFACE OUTPUT FROM MODEL'
        if (bats) then
          print*,'  WRITING ALL BATS DATA'
        end if
        if (batavg) then
          print*,'  AVERAGING ALL BATS DATA'
        end if
        if (batdiur) then
          print*,'  AVERAGING ALL DIURNAL BATS DATA'
        end if
        if (batday) then
          print*,'  WRITING DAILY BATS DATA'
        end if
      end if
      if (radproc.eq.1) then
        print*,' '
        print*,'PROCESS RADATION OUTPUT FROM MODEL'
        if (rad) then
          print*,'  WRITING ALL RAD DATA'
        end if
        if (radavg) then
          print*,'  AVERAGING ALL RAD DATA'
        end if
        if (raddiur) then
          print*,'  AVERAGING ALL DIURNAL RAD DATA'
        end if
        if (radday) then
          print*,'  WRITING DAILY RAD DATA'
        end if
      end if
      if (cheproc.eq.1) then
        print*,' '
        print*,'PROCESS CHEMTRAC OUTPUT FROM MODEL'
        if (che) then
          print*,'  WRITING ALL CHE DATA'
        end if
        if (cheavg) then
          print*,'  AVERAGING ALL CHE DATA'
        end if
        if (chediur) then
          print*,'  AVERAGING ALL DIURNAL CHE DATA'
        end if
        if (cheday) then
          print*,'  WRITING DAILY CHE DATA'
        end if
      end if
      if (subproc.eq.1) then
        print*,' '
        print*,'  PROCESS SURFACE OUTPUT FROM MODEL'
        if (sub) then
          print*,'  WRITING ALL SUB (SUB-BATS) DATA'
        end if
        if (subavg) then
          print*,'  AVERAGING ALL SUB DATA'
        end if
        if (subdiur) then
          print*,'  AVERAGING ALL DIURNAL SUB DATA'
        end if
        if (subday) then
          print*,'  WRITING DAILY SUB DATA'
        end if
      end if
      print*,' '
      print*,'**************************************************'
      print*,' '

      print*,idate0,idate1,idate2
      CALL JULIAN(idate0,julnc0,iyr0,imo0,idy0,ihr0)
      xhr0 = float(julnc0)
      CALL JULIAN(idate1,julnc1,iyr1,imo1,idy1,ihr1)
      xhr1 = float(julnc1)
      CALL JULIAN(idate2,julnc2,iyr2,imo2,idy2,ihr2)
      xhr2 = float(julnc2)
      julmid = (julnc1+julnc2)/2
      xhrm = float(julmid-mod(julmid,24))

C **** READ HEADER FILE **** C
      if (outhead) then
        fhdout = 'HEAD_OUT.NC'
        fhdbat = 'HEAD_BAT.NC'
        fhdsub = 'HEAD_SUB.NC'
        fhdrad = 'HEAD_RAD.NC'
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon
     &     , f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,mdate0
     &     , iin,inhead,idirect)
        print*,'HEADER READ IN'
        CALL PARAM(nx,ny,nz,nz,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,plv)
        CALL RCRECDF(fhdout,idout,vvarmin,vvarmax,ndim,ierr)
        CALL WRITEHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls
     &     , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhr0
     &     , iotyp)
        CALL CLSCDF(idout,ierr)
        CALL PARAM(nx,ny,1,1,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,plv)
        CALL RCRECDF(fhdbat,idout,vvarmin,vvarmax,ndim,ierr)
        CALL WRITEBATHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls
     &     , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhr0
     &     , iotyp)
        CALL CLSCDF(idout,ierr)
cSUB    CALL RCRECDF(fhdsub,idout,vvarmin,vvarmax,ndim,ierr)
cSUB    CALL WRITESUBHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls
cSUB &     , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhr0
cSUB &     , iotyp)
cSUB    CALL CLSCDF(idout,ierr)
        CALL PARAM(nx,ny,nz,nz,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,plv)
        CALL RCRECDF(fhdrad,idout,vvarmin,vvarmax,ndim,ierr)
        CALL WRITEHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls
     &     , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhr0
     &     , iotyp)
        CALL CLSCDF(idout,ierr)
      end if
C **** COMPUTE VVARMIN AND VVARMAX **** C
      close (iin)


C *********************** C
C ****   ICBC FILE   **** C
C *********************** C
      if (bc .or. bcavg .or. bcdiur .or. bcday) then
        if (bcday .and. (bcavg .or. bcdiur)) then
          print*,'MUST RUN BCDAY SEPARATE FROM BCAVG AND BCDIUR'
          stop 999
        end if
        if (bc) then
          filbc = 'ICBC'//trim(filinfo)//filext
          CALL FEXIST(filbc)
          print*,'ICBC FILE: ',filbc
        end if
        if (bcavg) then
          filavgbc = 'ICBC'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(filavgbc)
          print*,'ICBC AVERAGE FILE: ',filavgbc
        end if
        if (bcdiur) then
          fildiurbc = 'ICBC'//trim(filinfo)//'DIUR'//filext
          CALL FEXIST(fildiurbc)
          print*,'ICBC AVERAGE FILE: ',fildiurbc
        end if
        if (bcday) then
          fildaybc = 'ICBC'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(fildaybc)
          print*,'ICBC DAILY AVERAGE FILE: ',fildaybc
        end if
        CALL RDHEADICBC(nxf,nyf,nx,ny,nz,clat,clon,ds
     &     , pt,sigf,sigh,sighrev,xplat,xplon,f,xmap,dmap
     &     , xlat,xlon,dlat,dlon,zs,zssd,ls
     &     , iin,icbchead,ibyte)
        CALL PARAM(nx,ny,nz,npl,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &     , plv)
        ifil = 1
        irec = 0
        filrcm(ifil) = trim(icbcdir)//'/ICBC'//trim(rcmext(ifil))
        open (iin,file=filrcm(ifil),status='old',form='unformatted'
     &       ,access='direct',recl=nxf*nyf*ibyte)
        print*,'INPUT (ICBC) FILE: ',filrcm(ifil)
        if (bc) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filbc,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un1i,file=filbc,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr1i = 0
          end if
        end if
        if (bcday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildaybc,idday,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un2i,file=fildaybc,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr2i = 0
          end if
        end if
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUBC(vnambc,lnambc,ubc,xminbc,xmaxbc
     &     , factbc,offsetbc)
C **** ZERO OUT AVERAGE ARRAYS **** C
        if (bcavg .or. bcdiur .or. bcday) then
	    if (.not.plv) then
          CALL SETCONST(o3davg,vmisdat,nx,ny,nz,nbc2d,nhrbc,1,nx,1,ny)
	    else
	  CALL SETCONST(o3davg_p,vmisdat,nx,ny,npl,nbc2d,nhrbc,1,nx,1,ny)
	    endif
          CALL SETCONST(o2davg,vmisdat,nx,ny,nbc2d,nhrbc,1,1,nx,1,ny)
	    if(.not.plv) then
          CALL SETCONST(o3davg,0.0,nx,ny,nz,nbc2d,nhrbc,1,nx1,1,ny1)
	    else
	  CALL SETCONST(o3davg_p,0.0,nx,ny,npl,nbc2d,nhrbc,1,nx1,1,ny1)
	    endif
          CALL SETCONST(o2davg,0.0,nx,ny,nbc2d,nhrbc,1,1,nx1,1,ny1)
          do l=1,nhrbc
            nbctime(l) = 0
          end do
        end if
C **** READ IN DATA **** C
        idate = idate0
        idateold=idate
        do while(idate.le.idate2)
          idatenew = idateold + nint(dtbc)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)
          CALL RDICBC(vnambc,lnambc,ubc,idate,iin,irec,ierr)
C **** FOR END OF FILES **** C
          if (ierr.ne.0 .or. idate.gt.idatenew .or.
     &        (idate.gt.ircmext(ifil+1).and.ifil+1.le.nfiles)) then
            print*,'END OF FILE REACHED: ifil=',ifil,'ierr=',ierr
 50         ifil = ifil + 1
            filrcm(ifil) = trim(icbcdir)//'/ICBC'//trim(rcmext(ifil))
            print*,'  OPENING NEW FILE: ',filrcm(ifil)
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 95
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 50
            end if
            close (iin)
            irec = 0
            open (iin,file=filrcm(ifil),status='old',form='unformatted'
     &           ,access='direct',recl=nxf*nyf*ibyte)
            idate = 0
            do while(idate.lt.idatenew)
              print*,'SEARCHING FOR PROPER DAY:',idate,idatenew
              CALL RDICBC(vnambc,lnambc,ubc,idate,iin,irec,ierr)
              if (ierr.ne.0) stop 'READ ICBC ERROR'
              if (idate.gt.idatenew) then
                print*,'FILE ERROR (ICBC): DATE EXCEEDED'
                print*,idate,idatenew
                stop 999
              end if
            end do
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtbc)
          if (ihr.eq.0) ihr=24/nint(dtbc)
          CALL CALCRH(ifld2d,ifld3d,nx,ny,nz,nbc2d,nbc3d
     &       , sigh,pt,nti,nqvi,npsi,nrhi,ntdi,nthi,nx1,ny1)
          CALL HTSIG(ifld3d,ifld2d,nti,nhgti,npsi,zs,sigh,pt,nx,ny,nz
     &       ,  nx1,ny1,nbc3d,nbc2d)
c          CALL CALCHGT(ifld2d,ifld3d,nx,ny,nz,nbc2d,nbc3d
c     &       , zs,sigf,sigh,pt,nti,nqvi,npsi,nhgti,nx1,ny1)
	    CALL CALCSLP(ifld3d,ifld2d,nhgti,nti,npsi,pt,zs,nslpi,sigh
     &       , nx,ny,nz,nbc3d,nbc2d,nx1,ny1)
	    CALL CALCVD(ifld3d,nx,ny,nz,nbc3d,ds,dmap,xmap
     &       , nui,nvi,nvori,ndivi,nx1,ny1)
	    if (plv) then
	    CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nui,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nvi,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nqvi,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nrhi,plev,npl,nbc3d,nbc2d,nx1,ny1)
          CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nvori,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLIN(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , ndivi,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLOG(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nti,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLOG(ifld3d_p,ifld3d,ifld2d,npsi,pt,sigh,nx,ny,nz
     &       , nthi,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL INTLOG(ifld3d_p,ifld3d,ifld2d,npsa,pt,sigh,nx,ny,nz
     &       , ntdi,plev,npl,nbc3d,nbc2d,nx1,ny1)
	    CALL HEIGHT(ifld3d_p,ifld3d,ifld2d,nti,npsi,zs,sigh,nx,ny
     &       , nz,npl,nhgti,plev,nbc3d,nbc2d,pt,nx1,ny1)
	    endif
c **** AVERAGE DATA **** c
          if (bcavg .or. bcdiur .or. bcday) then
          if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(i2davg,ifld2d,nx,ny,nbc2d,nhrbc,ihr
     &         , vmisdat)
            if (.not.plv) then
            CALL AVGDATA3D(i3davg,ifld3d,nx,ny,nz,nbc3d,nhrbc,ihr
     &         , vmisdat)
		else
		CALL AVGDATA3D(i3davg_p,ifld3d_p,nx,ny,npl,nbc3d,nhrbc,ihr
     &         , vmisdat)
		endif
            nbctime(ihr) = nbctime(ihr) + 1
            if (bcday) then
              ntime = 0
              do l=1,nhrbc
                ntime = ntime + nbctime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtbc).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
		if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
		else
                  xhrdy = float(julncx)+dtbc-24.-(nday-1)*24.
		end if
		print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGBC(sighrev,vnambc,lnambc,ubc,xminbc
     &             , xmaxbc,factbc,offsetbc,vvarmin,vvarmax,xlat1d
     &             , xlon1d,idim,ndim,xhrdy,nbctime,idday,vmisdat
     &             , iotyp,un2i,nr2i,plv,u_bc)
                do l=1,nhrbc
                  nbctime(l) = 0
                end do
              end if
            end if
          end if
          end if
C **** WRITE ICBC DATA IN NETCDF FORMAT AT EACH TIME STEP **** c
          if (bc) then
          if (idate.ge.idate1.and.idate.le.idate2) then
            CALL WRITEBC(vvarmin,vvarmax,vnambc,lnambc,ubc,xminbc
     &         , xmaxbc,factbc,offsetbc,idim,ndim,xlat1d,xlon1d
     &         , sighrev,vmisdat,idout,xhr,iotyp,un1i,nr1i,plv,u_bc)
            print*,'DATA WRITTEN: ', xhr,idate
          end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtbc)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
        if (bc) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1i)
          end if
          print*,'DONE WRITING ICBC DATA!!!'
        end if
        if (bcday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2i)
          end if
          print*,'DONE WRITING ICBC DATA!!!'
        end if
      end if

c **** CLOSE INPUT (ICBC) FILE **** c
 95   close (iin)

C **** WRITE ICBC AVERAGED OUTPUT FIELDS **** C
      if (bcavg .or. bcdiur) then
        do ihr=1,nhrbc
          if (nbctime(ihr).le.0) then
            print*,'Not enough data for average.'
            print*,'nbctime must have values great than zero.'
            print*,nbctime
            stop 'BCAVG or BCDIUR'
          end if
        end do
      end if
      if (bcavg) then
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL RCRECDF(filavgbc,idout,vvarmin,vvarmax,ndim,ierr)
        elseif (iotyp.eq.3) then
          open(un3i,file=filavgbc,status='unknown'
     &        ,form='unformatted',recl=nx*ny*ibyte,access='direct')
          nr3i = 0
        end if
        CALL WRITEAVGBC(sighrev,vnambc,lnambc,ubc,xminbc,xmaxbc
     &     , factbc,offsetbc,vvarmin,vvarmax,xlat1d,xlon1d
     &     , idim,ndim,xhrm,nbctime,idout,vmisdat,iotyp,un3i,nr3i
     &     , plv,u_bc)
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL CLSCDF(idout,ierr)
        elseif (iotyp.eq.3) then
          close(un3i)
        end if
        print*,'DONE WRITING AVERAGED DATA IN COARDS CONVENTIONS!!!'
      end if
      if (bcdiur) then
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL RCRECDF(fildiurbc,idout,vvarmin,vvarmax,ndim,ierr)
        elseif (iotyp.eq.3) then
          open(un4i,file=fildiurbc,status='unknown'
     &        ,form='unformatted',recl=nx*ny*ibyte,access='direct')
          nr4i = 0
        end if
        CALL WRITEDIURBC(sighrev,vnambc,lnambc,ubc,xminbc,xmaxbc
     &     , factbc,offsetbc,vvarmin,vvarmax,xlat1d,xlon1d
     &     , idim,ndim,xhrm,nbctime,idout,vmisdat,iotyp,un4i,nr4i
     &     , plv,u_bc)
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL CLSCDF(idout,ierr)
        elseif (iotyp.eq.3) then
          close(un4i)
        end if
        print*,'DONE WRITING AVERAGED DIURNAL DATA!!!'
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if



C ********************** C
C ****   OUT FILE   **** C
C ********************** C
      if (out .or. outavg .or. outdiur .or. outday) then
        if (outday .and. (outavg .or. outdiur)) then
          print*,'MUST RUN OUTDAY SEPARATE FROM OUTAVG AND OUTDIUR'
          stop 999
        end if
        if (out) then
          filout = 'ATM'//trim(filinfo)//filext
          CALL FEXIST(filout)
          print*,'OUTPUT FILE: ',filout
        end if
        if (outavg) then
          filavgout = 'ATM'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(filavgout)
          print*,'OUTPUT AVERAGE FILE: ',filavgout
        end if
        if (outdiur) then
          fildiurout = 'ATM'//trim(filinfo)//'DIUR'//filext
          CALL FEXIST(fildiurout)
          print*,'OUTPUT AVERAGE FILE: ',fildiurout
        end if
        if (outday) then
          fildayout = 'ATM'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(fildayout)
          print*,'OUTPUT DAILY AVERAGE FILE: ',fildayout
        end if
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &     , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &     , zs,zssd,ls,mdate0,iin,inhead,idirect)
        CALL PARAM(nx,ny,nz,npl,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &     , plv)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/ATM.'//trim(rcmext(ifil))
        print*,'INPUT (OUT) FILE: ',filrcm(ifil)
        if (idirect.eq.1) then
          orec = 0
          open(iin,file=filrcm(ifil),status='old'
     &        ,form='unformatted',recl=ny*nx*ibyte
     &        ,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        endif
        if (out) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filout,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un1o,file=filout,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr1o = 0
          end if
        end if
        if (outday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildayout,idday,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un2o,file=fildayout,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr2o = 0
          end if
        end if
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUOUT(vnamout,lnamout,uout,xminout,xmaxout
     &     , factout,offsetout)
C **** ZERO OUT AVERAGE ARRAYS **** C
        if (outavg .or. outdiur .or. outday) then
	if (.not.plv) then
          CALL SETCONST(o3davg,vmisdat,nx,ny,nz,nout2d,nhrout,1,nx,1,ny)
	else
	  CALL SETCONST(o3davg_p,vmisdat,nx,ny,npl,nout2d,nhrout,1,nx,1,ny)
	end if
          CALL SETCONST(o2davg,vmisdat,nx,ny,nout2d,nhrout,1,1,nx,1,ny)
	if (.not.plv) then
          CALL SETCONST(o3davg,0.0,nx,ny,nz,nout2d,nhrout,1,nx1,1,ny1)
	else
	  CALL SETCONST(o3davg,0.0,nx,ny,npl,nout2d,nhrout,1,nx1,1,ny1)
	end if
          CALL SETCONST(o2davg,0.0,nx,ny,nout2d,nhrout,1,1,nx1,1,ny1)
          do l=1,nhrout
            nouttime(l) = 0
          end do
        end if
C **** READ IN DATA **** C
        idate = idate0
        if (mdate0.eq.idate0) then
          idate = idate0
        else
          idate = idate0 + nint(dtout)
        end if
        idateold=idate
        do while(idate.le.idate2)
          idatenew = idateold + nint(dtout)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)
          CALL RDATM(vnamout,lnamout,uout,idate,iin,orec,idirect,ierr)
C **** FOR END OF FILES **** C
          if (ierr.ne.0 .or. idate.gt.idatenew .or.
     &        (idate.gt.ircmext(ifil+1).and.ifil+1.le.nfiles)) then
            print*,'END OF FILE REACHED: ifil=',ifil,'ierr=',ierr
 51         ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/ATM.'//trim(rcmext(ifil))
            print*,'  OPENING NEW FILE: ',filrcm(ifil)
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 99
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 51
            end if
            close (iin)
            idate = idatenew
            if (idirect.eq.1) then
              orec = 0
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted',recl=ny*nx*ibyte
     &            ,access='direct')
            else
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted')
            endif
            CALL RDATM(vnamout,lnamout,uout,idate,iin,orec,idirect
     &         , ierr)
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtout)
          if (ihr.eq.0) ihr=24/nint(dtout)
          CALL CALCRH(ofld2d,ofld3d,nx,ny,nz,nout2d,nout3d
     &       , sigh,pt,nta,nqva,npsa,nrh,ntda,ntha,nx1,ny1)
          CALL HTSIG(ofld3d,ofld2d,nta,nhgt,npsa,zs,sigh,pt,nx,ny,nz
     &       ,  nx1,ny1,nout3d,nout2d)
c          CALL CALCHGT(ifld2d,ifld3d,nx,ny,nz,nbc2d,nbc3d
c     &       , zs,sigf,sigh,pt,nti,nqvi,npsi,nhgti,nx1,ny1)
	    CALL CALCSLP(ofld3d,ofld2d,nhgt,nta,npsa,pt,zs,nslp,sigh
     &       , nx,ny,nz,nout3d,nout2d,nx1,ny1)
	    CALL CALCVD(ofld3d,nx,ny,nz,nout3d,ds,dmap,xmap
     &       , nua,nva,nvora,ndiva,nx1,ny1)
	    if (plv) then
	    CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nua,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nva,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nqva,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nrh,plev,npl,nout3d,nout2d,nx1,ny1)
          CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nvora,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLIN(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , ndivi,plev,npl,nout3d,nout2d,nx1,ny1)
c            print*, ifld3d(:,:,:,nti)
	    CALL INTLOG(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , nta,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLOG(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , ntha,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL INTLOG(ofld3d_p,ofld3d,ofld2d,npsa,pt,sigh,nx,ny,nz
     &       , ntda,plev,npl,nout3d,nout2d,nx1,ny1)
	    CALL HEIGHT(ofld3d_p,ofld3d,ofld2d,nti,npsa,zs,sigh,nx,ny
     &       , nz,npl,nhgt,plev,nout3d,nout2d,pt,nx1,ny1)
	    endif
c **** AVERAGE DATA **** c
          if (outavg .or. outdiur .or. outday) then
          if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(o2davg,ofld2d,nx,ny,nout2d,nhrout,ihr
     &         , vmisdat)
     	if (.not.plv) then
            CALL AVGDATA3D(o3davg,ofld3d,nx,ny,nz,nout3d,nhrout,ihr
     &         , vmisdat)
     	else
     	    CALL AVGDATA3D(o3davg_p,ofld3d_p,nx,ny,npl,nout3d,nhrout,ihr
     &         , vmisdat)
     	end if
            nouttime(ihr) = nouttime(ihr) + 1
            if (outday) then
              ntime = 0
              do l=1,nhrout
                ntime = ntime + nouttime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtout).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
		if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
		else
                  xhrdy = float(julncx)+dtout-24.-(nday-1)*24.
		end if
		print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGOUT(sighrev,vnamout,lnamout,uout,xminout
     &             , xmaxout,factout,offsetout,vvarmin,vvarmax,xlat1d
     &             , xlon1d,idim,ndim,xhrdy,nouttime,idday,vmisdat
     &             , iotyp,un2o,nr2o,plv,u_out)
                do l=1,nhrout
                  nouttime(l) = 0
                end do
              end if
            end if
          end if
          end if
C **** WRITE OUT DATA IN NETCDF FORMAT AT EACH TIME STEP **** c
          if (out) then
          if (idate.ge.idate1.and.idate.le.idate2) then
            CALL WRITEOUT(vvarmin,vvarmax,vnamout,lnamout,uout,xminout
     &         , xmaxout,factout,offsetout,idim,ndim,xlat1d,xlon1d
     &         , sighrev,vmisdat,idout,xhr,iotyp,un1o,nr1o,plv,u_out)
            print*,'DATA WRITTEN: ', xhr,idate
          end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtout)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
        if (out) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1o)
          end if
          print*,'DONE WRITING OUT DATA!!!'
        end if
        if (outday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2o)
          end if
          print*,'DONE WRITING OUT DATA!!!'
        end if
      end if

c **** CLOSE INPUT (OUT) FILE **** c
 99   close (iin)

C **** WRITE OUT AVERAGED OUTPUT FIELDS **** C
      if (outavg .or. outdiur) then
        do ihr=1,nhrout
          if (nouttime(ihr).le.0) then
            print*,'Not enough data for average.'
            print*,'nouttime must have values great than zero.'
            print*,nouttime
            stop 'OUTAVG-OUTDIUR'
          end if
        end do
      end if
      if (outavg) then
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL RCRECDF(filavgout,idout,vvarmin,vvarmax,ndim,ierr)
        elseif (iotyp.eq.3) then
          open(un3o,file=filavgout,status='unknown'
     &        ,form='unformatted',recl=nx*ny*ibyte,access='direct')
          nr3o = 0
        end if
        CALL WRITEAVGOUT(sighrev,vnamout,lnamout,uout,xminout,xmaxout
     &     , factout,offsetout,vvarmin,vvarmax,xlat1d,xlon1d
     &     , idim,ndim,xhrm,nouttime,idout,vmisdat,iotyp,un3o,nr3o
     &	   , plv,u_out)
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL CLSCDF(idout,ierr)
        elseif (iotyp.eq.3) then
          close(un3o)
        end if
        print*,'DONE WRITING AVERAGED DATA IN COARDS CONVENTIONS!!!'
      end if
      if (outdiur) then
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL RCRECDF(fildiurout,idout,vvarmin,vvarmax,ndim,ierr)
        elseif (iotyp.eq.3) then
          open(un4o,file=fildiurout,status='unknown'
     &        ,form='unformatted',recl=nx*ny*ibyte,access='direct')
          nr4o = 0
        end if
        CALL WRITEDIUROUT(sighrev,vnamout,lnamout,uout,xminout,xmaxout
     &     , factout,offsetout,vvarmin,vvarmax,xlat1d,xlon1d
     &     , idim,ndim,xhrm,nouttime,idout,vmisdat,iotyp,un4o,nr4o
     &     , plv,u_out)
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL CLSCDF(idout,ierr)
        elseif (iotyp.eq.3) then
          close(un4o)
        end if
        print*,'DONE WRITING AVERAGED DIURNAL DATA!!!'
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if

C ************************* C
C ****    BATS FILE    **** C
C ************************* C
      if (bats .or. batavg .or. batdiur .or. batday) then
        if (batday .and. (batavg .or. batdiur)) then
          print*,'MUST RUN BATDAY SEPARATE FROM BATAVG AND BATDIUR'
          stop 999
        end if
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &     , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &     , zs,zssd,ls,mdate0,iin,inhead,idirect)
C **** OPEN NetCDF FILE **** C
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/SRF.'//trim(rcmext(ifil))
        print*,'INPUT (BATS) FILE: ',filrcm(ifil)
        if (idirect.eq.1) then
          brec = 0
          open(iin,file=filrcm(ifil),status='old'
     &        ,form='unformatted',recl=ny*nx*ibyte
     &        ,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        endif
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUBAT(vnambat,lnambat,ubat,xminbat,xmaxbat
     &     , factbat,offsetbat,nbat2)
C **** COMPUTE VVARMIN AND VVARMAX **** C
        CALL PARAM(nx,ny,1,1,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &     , plv)
        vvarmax(3) = 1050.
        idim(3) = 1
        CALL SETCONST(sigb,1.0,2,1,1,1,1,2,1,1,1)
        if (bats) then
          filbat = 'SRF'//trim(filinfo)//filext
          CALL FEXIST(filbat)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (filbat,idout,vvarmin,vvarmax,ndim,ierr)
            print*,'OPENING BATS NetCDF FILE',idout
          elseif (iotyp.eq.3) then
            open(un1b,file=filbat,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr1b = 0
            print*,'OPENING BATS GrADS FILE',un1b
          end if
        end if
        if (batavg .or. batdiur .or. batday) then
          if (batavg) then
            filavgbat = 'SRF'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(filavgbat)
          end if
          if (batdiur) then
            fildiurbat = 'SRF'//trim(filinfo)//'DIUR'//filext
            CALL FEXIST(fildiurbat)
          end if
          if (batday) then
            fildaybat = 'SRF'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(fildaybat)
            if (iotyp.eq.1 .or. iotyp.eq.2) then
              CALL RCRECDF (fildaybat,idday,vvarmin,vvarmax,ndim,ierr)
            elseif (iotyp.eq.3) then
              open(un2b,file=fildaybat,status='unknown'
     &            ,form='unformatted',recl=nx*ny*ibyte,access='direct')
              nr2b = 0
            end if
          end if
          CALL SETCONST(b2davg,0.0,nx,ny,nbat2,nhrbat,1,1,nx,1,ny)
          do l=1,nhrbat
            nbattime(l) = 0
          end do
        end if
        if (mdate0.eq.idate0) then
          idate = idate0
        else
          idate = idate0 + nint(dtbat)
        end if
        idateold = idate
C **** Initialize Max and Min Temperatures
c        do j=1,ny
c        do i=1,nx
c        bfld2d(i,j,ntmax) = bfld2d(i,j,ntanm)
c          bfld2d(i,j,ntmin) = bfld2d(i,j,ntanm)
c        end do
c        end do
        print*,idate,idate0,idate1,idate2
        do while(idate.le.idate2)
C **** READ BATS FILE **** C
          idatenew = idateold + nint(dtbat)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)
          CALL RDSRF(vnambat,lnambat,ubat,idate,iin,brec,idirect,ierr)
          if (ierr.ne.0 .or. idate.gt.idatenew .or.
     &        (idate.gt.ircmext(ifil+1).and.ifil+1.le.nfiles)) then
            print*,'READ ERROR OR END OF FILE REACHED',ierr
 52         ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/SRF.'//trim(rcmext(ifil))
            print*,'  OPENING NEW FILE: ',filrcm(ifil)
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 98
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 52
            end if
            close (iin)
            idate = idatenew
            if (idirect.eq.1) then
              brec = 0
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted',recl=nx*ny*ibyte
     &            ,access='direct')
            else
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted')
            end if
            CALL RDSRF(vnambat,lnambat,ubat,idate,iin,brec,idirect
     &         , ierr)
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtbat)+1
          if (ihr.eq.0) ihr=24/nint(dtbat)
c          CALL TMINMAX(bfld2d,b2davg,nx,ny,nbat2,nhrbat,ihr
c     &       , ntanm,ntmax,ntmin)
c          CALL CALCMSE2D(bfld2d,zs,nx,ny,nbat2,ntanm,nqanm,nmsea)
          CALL CALCRH2D(bfld2d,nx,ny,nbat2
     &       , ntanm,nqanm,npsrf,nrha,vmisdat)
c **** AVERAGE DATA **** c
          if (batavg .or. batdiur .or. batday) then
          if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
            CALL AVGDATABAT(ihr,vmisdat)
            nbattime(ihr) = nbattime(ihr) + 1
            if (batday) then
              ntime = 0
              do l=1,nhrbat
                ntime = ntime + nbattime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtbat).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
                if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
                else
                  xhrdy = float(julncx)+dtbat-24.-(nday-1)*24.
                end if
                print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGBAT(vmisdat,vnambat,lnambat,ubat,xminbat
     &             , xmaxbat,factbat,offsetbat,vvarmin,vvarmax
     &             , xlat1d,xlon1d,idim,ndim,xhrdy,nbattime,idday
     &             , iotyp,un2b,nr2b,u_bat)
                do l=1,nhrbat
                  nbattime(l) = 0
                  do nb=1,nbat2
                  do j=nb1,nyb
                  do i=nb1,nxb
                    b2davg(i,j,nb,l) = 0.0
                  end do
                  end do
                  end do
                end do
              end if
            end if
          end if
          end if
C **** WRITE BATS IN COARDS NetCDF CONVENTIONS **** C
          if (bats) then
            if (idate.ge.idate1.and.idate.le.idate2) then
             CALL WRITEBAT(vnambat,lnambat,ubat,xminbat,xmaxbat
     &          , factbat,offsetbat,vvarmin,vvarmax,xlat1d,xlon1d
     &          , idim,ndim,vmisdat,xhr,ihr,idout,iotyp,un1b,nr1b
     &          , u_bat)
              print*,'BATS DATA WRITTEN:  ', idate
              print*,''
            end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtbat)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
 98     continue
        close (iin)
        if (bats) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1b)
          end if
          print*,'DONE WRITING BATS DATA!!!'
        end if
        if (batday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2b)
          end if
          print*,'DONE WRITING DAILY BATS DATA!!!'
        end if
        if (batavg .or. batdiur) then
          do ihr=1,nhrbat
            if (nbattime(ihr).le.0) then
              print*,'Not enough data for average.'
              print*,'nbattime must have values great than zero.'
              print*,nbattime
              stop 'BATAVG-BATDIUR'
            end if
          end do
        end if
        if (batavg) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filavgbat,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un3b,file=filavgbat,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr3b = 0
          end if
          CALL WRITEAVGBAT(vmisdat,vnambat,lnambat,ubat,xminbat
     &       , xmaxbat,factbat,offsetbat,vvarmin,vvarmax,xlat1d
     &       , xlon1d,idim,ndim,xhrm,nbattime,idout,iotyp,un3b,nr3b
     &       , u_bat)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un3b)
          end if
          print*,'DONE WRITING AVERAGED BATS DATA'
        end if
        if (batdiur) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildiurbat,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un4b,file=fildiurbat,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr4b = 0
          end if
          CALL WRITEDIURBAT(vmisdat,vnambat,lnambat,ubat,xminbat
     &       , xmaxbat,factbat,offsetbat,vvarmin,vvarmax,xlat1d
     &       , xlon1d,idim,ndim,xhrm,nbattime,idout,iotyp,un4b,nr4b
     &       , u_bat)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un4b)
          end if
          print*,'DONE WRITING AVERAGED DIURNAL BATS DATA'
        end if
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if


C ****************************** C
C ****    RADIATION FILE    **** C
C ****************************** C
      if (rad .or. radavg .or. raddiur .or. radday) then
        if (radday .and. (radavg .or. raddiur)) then
          print*,'MUST RUN RADDAY SEPARATE FROM RADAVG AND RADDIUR'
          stop 999
        end if
C **** OPEN NetCDF FILE **** C
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &     , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &     , zs,zssd,ls,mdate0,iin,inhead,idirect)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/RAD.'//trim(rcmext(ifil))
        if (idirect.eq.1) then
          rrec = 0
          open(iin,file=filrcm(ifil),status='old'
     &        ,form='unformatted',recl=nx*ny*ibyte
     &        ,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        end if
        print*,'INPUT (RAD) FILE: ',filrcm(ifil)
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLURAD(vnamrad,lnamrad,urad,xminrad,xmaxrad
     &     , factrad,offsetrad)
C **** COMPUTE VVARMIN AND VVARMAX **** C
        CALL PARAM(nx,ny,nz,npl,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &     ,plv)
        idim(3) = nz
        if (rad) then
          filrad = 'RAD'//trim(filinfo)//filext
          CALL FEXIST(filrad)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (filrad,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un1r,file=filrad,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr1r = 0
          end if
        end if
        if (radday) then
          fildayrad = 'RAD'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(fildayrad)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (fildayrad,idday,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un2r,file=fildayrad,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr2r = 0
          end if
        end if
        if (radavg .or. raddiur .or. radday) then
          if (radavg) then
            filavgrad = 'RAD'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(filavgrad)
          end if
          if (raddiur) then
            fildiurrad = 'RAD'//trim(filinfo)//'DIUR'//filext
            CALL FEXIST(fildiurrad)
          end if
	  if (.not.plv) then
          CALL SETCONST(r3davg,vmisdat,nx,ny,nz,nr2d,nhrrad,1,nx,1,ny)
	  else
	  CALL SETCONST(r3davg_p,vmisdat,nx,ny,npl,nr2d,nhrrad,1,nx,1,ny)
	  end if
          CALL SETCONST(r2davg,vmisdat,nx,ny,nr2d,nhrrad,1,1,nx,1,ny)
	  if (.not.plv) then
          CALL SETCONST(r3davg,0.0,nx,ny,nz,nr2d,nhrrad,1,nx1,1,ny1)
	  else
	  CALL SETCONST(r3davg_p,0.0,nx,ny,npl,nr2d,nhrrad,1,nx1,1,ny1)
	  end if
          CALL SETCONST(r2davg,0.0,nx,ny,nr2d,nhrrad,1,1,nx1,1,ny1)
          do l=1,nhrrad
            nradtime(l) = 0
          end do
        end if
        if (mdate0.eq.idate0) then
          idate = idate0
        else
          idate = idate0 + nint(dtrad)
        end if
        idateold = idate
        do while(idate.le.idate2)
C **** READ RADIATION FILE **** C
          idatenew = idateold + nint(dtrad)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)
          CALL RDRAD(vnamrad,lnamrad,iin,idate,rrec,idirect,ierr)
          if (ierr.ne.0 .or. idate.gt.idatenew .or.
     &        (idate.gt.ircmext(ifil+1).and.ifil+1.le.nfiles)) then
            print*,'END OF FILE REACHED: ifil=',ifil,'ierr=',ierr
 53         ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/RAD.'//trim(rcmext(ifil))
            print*,'  OPENING NEW FILE: ',filrcm(ifil)
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 97
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 53
            end if
            close (iin)
            idate = idatenew
            if (idirect.eq.1) then
              rrec = 0
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted',recl=nx*ny*ibyte
     &            ,access='direct')
            else
              open(iin,file=filrcm(ifil),status='old'
     &           ,form='unformatted')
            end if
            CALL RDRAD(vnamrad,lnamrad,iin,idate,rrec,idirect,ierr)
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtrad)
          if (ihr.eq.0) ihr=24/nint(dtrad)
	  if (plv) then
	  CALL INTLIN(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,nx,ny,nz
     &       , ncld,plev,npl,nr3d,nr2d,nx1,ny1)
	  CALL INTLIN(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,nx,ny,nz
     &       , nclwp,plev,npl,nr3d,nr2d,nx1,ny1)
	  CALL INTLIN(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,nx,ny,nz
     &       , nqrs,plev,npl,nr3d,nr2d,nx1,ny1)
	  CALL INTLIN(rfld3d_p,rfld3d,rfld2d,npsrf,pt,sigh,nx,ny,nz
     &       , nqrl,plev,npl,nr3d,nr2d,nx1,ny1)
          end if

c **** AVERAGE DATA **** c
          if (radavg .or. raddiur .or. radday) then
          if (idate.ge.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(r2davg,rfld2d,nx,ny,nr2d,nhrout,ihr
     &         , vmisdat)
          if (.not.plv) then
            CALL AVGDATA3D(r3davg,rfld3d,nx,ny,nz,nr3d,nhrout,ihr
     &         , vmisdat)
          else
	    CALL AVGDATA3D(r3davg_p,rfld3d_p,nx,ny,npl,nr3d,nhrout,ihr
     &         , vmisdat)
          end if
            nradtime(ihr) = nradtime(ihr) + 1
            if (radday) then
              ntime = 0
              do l=1,nhrrad
                ntime = ntime + nradtime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtrad).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
		if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
		else
                  xhrdy = float(julncx)+dtrad-24.-(nday-1)*24.
		end if
		print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGRAD(xhrdy,sighrev,vnamrad,lnamrad,urad
     &       ,       xminrad,xmaxrad,factrad,offsetrad,vvarmin,vvarmax
     &             , xlat1d,xlon1d,idim,ndim,vmisdat,nradtime,idday
     &             , iotyp,un2r,nr2r,plv,u_rad)
                print*,'DAILY RAD DATA WRITTEN: ',xhr,idate
                do l=1,nhrrad
                  nradtime(l) = 0
                end do
              end if
            end if
            end if
          end if
C **** WRITE RADIATION IN IVE NetCDF CONVENTIONS **** C
          if (rad) then
            if (idate.ge.idate1.and.idate.le.idate2) then
             CALL WRITERAD(vnamrad,lnamrad,urad,xminrad,xmaxrad,factrad
     &          , offsetrad,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &          , sighrev,vmisdat,idout,xhr,iotyp,un1r,nr1r,plv,u_rad)
              print*,'RADIATION DATA WRITTEN:  ', idate
              print*,''
            end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtrad)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
 97     continue
        close (iin)
        if (rad) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1r)
          end if
          print*,'DONE RADIATION DATA!!!'
        end if
        if (radday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2r)
          end if
          print*,'DONE DAILY RADIATION DATA!!!'
        end if
        if (radavg .or. raddiur) then
          do ihr=1,nhrrad
             print*,nradtime, 'is nradtime'
          if (nradtime(ihr).le.0) then
              print*,'Not enough data for average.'
              print*,'nradtime must have values great than zero.'
              print*,nradtime
              stop 'RADAVG-RADDIUR'
            end if
          end do
        end if
        if (radavg) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filavgrad,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un3r,file=filavgrad,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr3r = 0
          end if
          CALL WRITEAVGRAD(xhrm,sighrev,vnamrad,lnamrad,urad,xminrad
     &       , xmaxrad,factrad,offsetrad,vvarmin,vvarmax,xlat1d,xlon1d
     &       , idim,ndim,vmisdat,nradtime,idout,iotyp,un3r,nr3r
     &       , plv,u_rad)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un3r)
          end if
          print*,'DONE WRITING AVERAGED RADIATION DATA'
        end if
        if (raddiur) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildiurrad,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un4r,file=filavgrad,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr4r = 0
          end if
          CALL WRITEDIURRAD(xhrm,sighrev,vnamrad,lnamrad,urad,xminrad
     &       , xmaxrad,factrad,offsetrad,vvarmin,vvarmax,xlat1d,xlon1d
     &       , idim,ndim,vmisdat,nradtime,idout,iotyp,un4r,nr4r
     &       , plv,u_rad)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un4r)
          end if
          print*,'DONE WRITING AVERAGED DIURNAL RADIATION DATA'
        end if
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if


C ****************************** C
C ****    CHEM-TRACER FILE  **** C
C ****************************** C
      if (che .or. cheavg .or. chediur .or. cheday) then
        if (cheday .and. (cheavg .or. chediur)) then
          print*,'MUST RUN RADDAY SEPARATE FROM CHEAVG AND CHEDIUR'
          stop 999
        end if
C **** OPEN NetCDF FILE **** C
        print*,'toto',trim(rcmext(ifil))
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &     , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &     , zs,zssd,ls,mdate0,iin,inhead,idirect)
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/CHE.'//trim(rcmext(ifil))
        print*,'INPUT CHEM FILE: ',filrcm
        if (idirect.eq.1) then
          crec = 0
          open(iin,file=filrcm(ifil),status='old'
     &        ,form='unformatted',recl=ny*nx*ibyte
     &        ,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        endif
        print*,'INPUT (CHE) FILE: ',filrcm(ifil)
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUCHE(vnamche,lnamche,uche,xminche,xmaxche
     &     , factche,offsetche)
        print*,vnamche

C **** COMPUTE VVARMIN AND VVARMAX **** C
        CALL PARAM(nx,ny,nz,npl,ds,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &     ,plv)
c        idim(3) = nz
        if (che) then
          filche = 'CHE'//trim(filinfo)//filext
          print*,filche
          CALL FEXIST(filche)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (filche,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un1c,file=filche,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr1c = 0
          end if
        end if
        if (cheday) then
          fildayche = 'CHE'//trim(filinfo)//'AVG'//filext
          CALL FEXIST(fildayche)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (fildayche,idday,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un2c,file=fildayche,status='unknown'
     &          ,form='unformatted',recl=nx*ny*4,access='direct')
            nr2c = 0
          end if
        end if
        if (cheavg .or. chediur .or. cheday) then
          if (cheavg) then
            filavgche = 'CHE'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(filavgche)
          end if
          if (chediur) then
             fildiurche = 'CHE'//trim(filinfo)//'DIUR'//filext
            CALL FEXIST(fildiurche)
          end if
	  if (.not.plv) then
          CALL SETCONST(c3davg,vmisdat,nx,ny,nz,nc3d,nhrche,1,nx,1,ny)
	  else
	  CALL SETCONST(c3davg_p,vmisdat,nx,ny,npl,nc3d,nhrche,1,nx,1,ny)
	  end if
          CALL SETCONST(c2davg,vmisdat,nx,ny,nc2d,nhrche,1,1,nx,1,ny)
	  if (.not.plv) then
          CALL SETCONST(c3davg,0.0,nx,ny,nz,nc3d,nhrche,1,nx1,1,ny1)
	  else
	  CALL SETCONST(c3davg_p,0.0,nx,ny,npl,nc3d,nhrche,1,nx1,1,ny1)
	  end if
          CALL SETCONST(c2davg,0.0,nx,ny,nc2d,nhrche,1,1,nx1,1,ny1)
          do l=1,nhrche
            nchetime(l) = 0
          end do
        end if
        if (mdate0.eq.idate0) then
          idate = idate0
        else
          idate = idate0 + nint(dtche)
        end if
        idateold = idate
        print*,idate,idate0,idate1,idate2
         do while(idate.le.idate2)
C **** READ CHEM-TRACER FILE **** C
          idatenew = idateold + nint(dtche)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)    
          CALL RDCHE(iin,idate,crec,idirect,ierr)       
          if (idate.gt.idatenew) then
            print*,'END OF FILE REACHED: ifil=',ifil,'ierr=',ierr
 54         ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/CHE.'//trim(rcmext(ifil))
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 103
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 54
            end if
            close (iin)
            idate = idatenew
            if (idirect.eq.1) then
              crec = 0
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted',recl=ny*nx*ibyte
     &            ,access='direct')
            else
              open(iin,file=filrcm(ifil),status='old'
     &           ,form='unformatted')
            endif
            CALL RDCHE(iin,idate,crec,idirect,ierr)
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
csr maybe add one here? 
          ihr = ihr/nint(dtche) 
          if (ihr.eq.0) ihr=24/nint(dtche)
c **** AVERAGE DATA **** c
          if (cheavg .or. chediur .or. cheday) then
          if (idate.ge.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
Csr nhrcche
            CALL AVGDATA2D(c2davg,cfld2d,nx,ny,nc2d,nhrche,ihr
     &         , vmisdat)
         if (.not.plv) then
         CALL AVGDATA3D(c3davg,cfld3d,nx,ny,nz,nc3d,nhrche,ihr
     &         , vmisdat)
         else
	 CALL AVGDATA3D(c3davg_p,cfld3d_p,nx,ny,npl,nc3d,nhrche,ihr
     &         , vmisdat)
         end if
            nchetime(ihr) = nchetime(ihr) + 1
            if (cheday) then
              ntime = 0
              do l=1,nhrche
                ntime = ntime + nchetime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtche).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
		if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
		else
                  xhrdy = float(julncx)+dtche-24.-(nday-1)*24.
		end if
		print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGCHE(xhrdy,sighrev,vnamche,lnamche,uche
     &       ,       xminche,xmaxche,factche,offsetche,vvarmin,vvarmax
     &             , xlat1d,xlon1d,idim,ndim,vmisdat,nchetime,idday
     &             , iotyp,un2c,nr2c,plv,u_che)
                print*,'DAILY CHE DATA WRITTEN: ',xhr,idate
                do l=1,nhrche
                  nchetime(l) = 0
                end do
              end if
            end if
            end if
          end if
C **** WRITE RADIATION IN IVE NetCDF CONVENTIONS **** C
          if (che) then
            if (idate.ge.idate1.and.idate.le.idate2) then
             CALL WRITECHE(vnamche,lnamche,uche,xminche,xmaxche,factche
     &          , offsetche,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &          , sighrev,vmisdat,idout,xhr,iotyp,un1c,nr1c,plv,u_che)
              print*,'CHEM-TRACER DATA WRITTEN:  ', idate
              print*,''
            end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtche)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
 103     continue
        close (iin)
        if (che) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1c)
          end if
          print*,'DONE CHEM-TRACER DATA!!!'
        end if
        if (cheday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2c)
          end if
          print*,'DONE DAILY CHEM-TRACER DATA!!!'
        end if
        if (cheavg .or. chediur) then
          do ihr=1,nhrche
                print*,nchetime, ' is nchetime'
                print*,ihr
                if (nchetime(ihr).le.0) then
              print*,'Not enough data for average.'
              print*,'nchetime must have values great than zero.'
              print*,cheavg
              print*,chediur
              print*,'nchetime is', nchetime
              stop 'CHEAVG-CHEDIUR'
            end if
          end do
        end if
        if (cheavg) then
          print*, 'YES',filavgche 
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filavgche,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un3c,file=filavgche,status='unknown'
     &          ,form='unformatted',recl=nx*ny*ibyte,access='direct')
            nr3c = 0
          end if
          CALL WRITEAVGCHE(xhrm,sighrev,vnamche,lnamche,uche,xminche
     &       , xmaxche,factche,offsetche,vvarmin,vvarmax,xlat1d,xlon1d
     &       , idim,ndim,vmisdat,nchetime,idout,iotyp,un3c,nr3c
     &       , plv,u_che)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un3c)
          end if
          print*,'DONE WRITING AVERAGED CHE-TRACER DATA'
        end if
        if (chediur) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildiurche,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un4c,file=filavgche,status='unknown'
     &          ,form='unformatted',recl=nx*ny*4,access='direct')
            nr4c = 0
          end if
          CALL WRITEDIURCHE(xhrm,sighrev,vnamche,lnamche,uche,xminche
     &       , xmaxche,factche,offsetche,vvarmin,vvarmax,xlat1d,xlon1d
     &       , idim,ndim,vmisdat,nchetime,idout,iotyp,un4c,nr4c
     &       , plv,u_che)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un4c)
          end if
          print*,'DONE WRITING AVERAGED DIURNAL CHE-TRACER DATA'
        end if
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if


C ************************* C
C ****    SUB FILE     **** C
C ************************* C
      if (sub .or. subavg .or. subdiur .or. subday) then
        if (subday .and. (subavg .or. subdiur)) then
          print*,'MUST RUN SUBDAY SEPARATE FROM SUBAVG AND SUBDIUR'
          stop 999
        end if
        CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &     , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &     , zs,zssd,ls,mdate0,iin,inhead,idirect)
        CALL RDHEADICBC(nxsf,nysf,nxsb,nysb,nz,clat,clon,dssb
     &     , pt,sigf,sigh,sighrev,xplat,xplon,fsb,xmapsb,dmapsb
     &     , xlatsb,xlonsb,dlatsb,dlonsb,zssb,zssdsb,lssb
     &     , mdate0,iin,icbcheadsb,ibyte)
        CALL PARAM(nxsb,nysb,1,1,dssb,clat,clon,xplat,xplon
     &     , xlatsb,xlonsb,vvarmin,vvarmax,xlatsb1d,xlonsb1d
     &     , idim,ndim, plv)
        print*,'            '
        print*,'            '
        print*,vvarmin
        print*,vvarmax
        print*,dssb
        print*,'            '
        print*,'            '
C **** OPEN NetCDF FILE **** C
        ifil = 1
        filrcm(ifil) = trim(rcmdir)//'/SUB.'//trim(rcmext(ifil))
        print*,'INPUT (SUB) FILE: ',filrcm(ifil)
        if (idirect.eq.1) then
          srec = 0
          open(iin,file=filrcm(ifil),status='old'
     &        ,form='unformatted',recl=nysb*nxsb*ibyte
     &        ,access='direct')
        else
          open (iin,file=filrcm(ifil),status='old',form='unformatted')
        endif
C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUSUB(vnamsub,lnamsub,usub,xminsub,xmaxsub
     &     , factsub,offsetsub,nsub2)
C **** COMPUTE VVARMIN AND VVARMAX **** C
        vvarmax(3) = 1050.
        idim(3) = 1
        CALL SETCONST(sigb,1.0,2,1,1,1,1,2,1,1,1)
        if (sub) then
          filsub = 'SUB'//trim(filinfo)//filext
          CALL FEXIST(filsub)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF (filsub,idout,vvarmin,vvarmax,ndim,ierr)
            print*,'OPENING SUB NetCDF FILE',idout
          elseif (iotyp.eq.3) then
            open(un1s,file=filsub,status='unknown',form='unformatted'
     &          ,recl=nxsb*nysb*ibyte,access='direct')
            nr1s = 0
            print*,'OPENING SUB GrADS FILE',un1s
          end if
        end if
        if (subavg .or. subdiur .or. subday) then
          if (subavg) then
            filavgsub = 'SUB'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(filavgsub)
          end if
          if (subdiur) then
            fildiursub = 'SUB'//trim(filinfo)//'DIUR'//filext
            CALL FEXIST(fildiursub)
          end if
          if (subday) then
            fildaysub = 'SUB'//trim(filinfo)//'AVG'//filext
            CALL FEXIST(fildaysub)
            if (iotyp.eq.1 .or. iotyp.eq.2) then
              CALL RCRECDF (fildaysub,idday,vvarmin,vvarmax,ndim,ierr)
            elseif (iotyp.eq.3) then
              open(un2s,file=fildaysub,status='unknown'
     &            ,form='unformatted',recl=nxsb*nysb*ibyte
     &            ,access='direct')
              nr2s = 0
            end if
          end if
          CALL SETCONST(s2davg,0.0,nxsb,nysb,nsub2,nhrsub
     &       , 1,1,nxsb,1,nysb)
          do l=1,nhrsub
            nsubtime(l) = 0
          end do
        end if
        if (mdate0.eq.idate0) then
          idate = idate0
        else
          idate = idate0 + nint(dtsub)
        end if
        idateold = idate
C **** Initialize Max and Min Temperatures
c        do j=1,nysb
c        do i=1,nxsb
c          sfld2d(i,j,nstmax) = sfld2d(i,j,nstanm)
c          sfld2d(i,j,nstmin) = sfld2d(i,j,nstanm)
c        end do
c        end do
        print*,idate,idate0,idate1,idate2
        do while(idate.le.idate2)
C **** READ SUB FILE **** C
          idatenew = idateold + nint(dtsub)
          CALL JULIAN(idatenew,julnc,iyr,imo,idy,ihr)
          CALL RDSUB(vnamsub,lnamsub,usub,idate,iin,srec,idirect,ierr)
          if (ierr.ne.0 .or. idate.gt.idatenew .or.
     &        (idate.gt.ircmext(ifil+1).and.ifil+1.le.nfiles)) then
            print*,'READ ERROR OR END OF FILE REACHED',ierr
            print*,idate,idatenew,idateold,ircmext(ifil+1)
 55         ifil = ifil + 1
            filrcm(ifil) = trim(rcmdir)//'/SUB.'//trim(rcmext(ifil))
            CALL FEXISTNEW(filrcm(ifil),there)
            if (.not.there) go to 94
            if (idateold.ne.ircmext(ifil)) then
              print*,' '
              print*,'INCONSISTENT DATES: ifil=',ifil
              print*,'  idateold=',idateold,'rcmext=',ircmext(ifil)
              go to 55
            end if
            close (iin)
            idate = idatenew
            if (idirect.eq.1) then
              srec = 0
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted',recl=nxsb*nysb*ibyte
     &            ,access='direct')
            else
              open(iin,file=filrcm(ifil),status='old'
     &            ,form='unformatted')
            end if
            CALL RDSUB(vnamsub,lnamsub,usub,idate,iin,srec,idirect
     &         , ierr)
          end if
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtsub)+1
          if (ihr.eq.0) ihr=24/nint(dtsub)
c          CALL TMINMAX(sfld2d,s2davg,nxsb,nysb,nsub2,nhrsub,ihr
c     &       , nstanm,nstmax,nstmin)
c          CALL CALCMSE2D(sfld2d,zssb,nxsb,nysb,nsub2
c     &       , nstanm,nsqanm,nsmsea)
          CALL CALCRH2D(sfld2d,nxsb,nysb,nsub2
     &       , nstanm,nsqanm,nspsrf,nsrha,vmisdat)
c **** AVERAGE DATA **** c
          if (subavg .or. subdiur .or. subday) then
          if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING DATA: ',idate,xhr,ihr
            CALL AVGDATASUB(ihr,vmisdat)
            nsubtime(ihr) = nsubtime(ihr) + 1
            if (subday) then
              ntime = 0
              do l=1,nhrsub
                ntime = ntime + nsubtime(l)
              end do
              CALL JULIAN(idateold,julncx,iyrx,imox,idyx,ihrx)
              if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtsub).and.nday.gt.0)) then
                idatex = iyrx*1000000 + imox*10000 + 1500
                if (nday.eq.-1) then
                  CALL JULIAN(idatex,julncx,iyrx,imox,idyx,ihrx)
                  xhrdy = float(julncx)
                else
                  xhrdy = float(julncx)+dtsub-24.-(nday-1)*24.
                end if
                print*,'WRITING CONTINUAL OUTPUT',xhrdy,idatex
                CALL WRITEAVGSUB(vmisdat,vnamsub,lnamsub,usub,xminsub
     &             , xmaxsub,factsub,offsetsub,vvarmin,vvarmax
     &             , xlatsb1d,xlonsb1d,idim,ndim,xhrdy,nsubtime,idday
     &             , iotyp,un2s,nr2s)
                do l=1,nhrsub
                  nsubtime(l) = 0
                  do ns=1,nsub2
                  do j=1,nysb
                  do i=1,nxsb
                    s2davg(i,j,ns,l) = 0.0
                  end do
                  end do
                  end do
                end do
              end if
            end if
          end if
          end if
C **** WRITE SUB IN COARDS NetCDF CONVENTIONS **** C
          if (sub) then
            if (idate.ge.idate1.and.idate.le.idate2) then
             CALL WRITESUB(vnamsub,lnamsub,usub,xminsub,xmaxsub
     &          , factsub,offsetsub,vvarmin,vvarmax,xlatsb1d,xlonsb1d
     &          , idim,ndim,vmisdat,xhr,ihr,idout,iotyp,un1s,nr1s)
              print*,'SUB DATA WRITTEN:  ', idate
              print*,''
            end if
          end if
C **** INCREMENT TIME **** C
          idateold = idate
          idate = idate + nint(dtsub)
          CALL JULIAN(idate,julnc,iyr,imo,idy,ihr)
        end do
 94     continue
        close (iin)
        if (sub) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un1s)
          end if
          print*,'DONE WRITING SUB DATA!!!'
        end if
        if (subday) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idday,ierr)
          elseif (iotyp.eq.3) then
            close(un2s)
          end if
          print*,'DONE WRITING DAILY SUB DATA!!!'
        end if
        if (subavg .or. subdiur) then
          do ihr=1,nhrsub
            if (nsubtime(ihr).le.0) then
              print*,'Not enough data for average.'
              print*,'nsubtime must have values great than zero.'
              print*,nsubtime
              stop 'SUBAVG-SUBDIUR'
            end if
          end do
        end if
        if (subavg) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(filavgsub,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un3s,file=filavgsub,status='unknown'
     &          ,form='unformatted',recl=nxsb*nysb*ibyte
     &          ,access='direct')
            nr3s = 0
          end if
          CALL WRITEAVGSUB(vmisdat,vnamsub,lnamsub,usub,xminsub
     &       , xmaxsub,factsub,offsetsub,vvarmin,vvarmax,xlatsb1d
     &       , xlonsb1d,idim,ndim,xhrm,nsubtime,idout,iotyp,un3s,nr3s)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un3s)
          end if
          print*,'DONE WRITING AVERAGED SUB DATA'
        end if
        if (subdiur) then
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL RCRECDF(fildiursub,idout,vvarmin,vvarmax,ndim,ierr)
          elseif (iotyp.eq.3) then
            open(un4s,file=fildiursub,status='unknown'
     &          ,form='unformatted',recl=nxsb*nysb*ibyte
     &          ,access='direct')
            nr4s = 0
          end if
          CALL WRITEDIURSUB(vmisdat,vnamsub,lnamsub,usub,xminsub
     &       , xmaxsub,factsub,offsetsub,vvarmin,vvarmax,xlatsb1d
     &       , xlonsb1d,idim,ndim,xhrm,nsubtime,idout,iotyp,un4s,nr4s)
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL CLSCDF(idout,ierr)
          elseif (iotyp.eq.3) then
            close(un4s)
          end if
          print*,'DONE WRITING AVERAGED DIURNAL SUB DATA'
        end if
        print*,'***** if you see <variable not found>, do NOT worry.'
      end if



      print*,'DONE!!!'

      end

      SUBROUTINE XTDOT( px, pd, ni, nj, nk, ni1, nj1 )               

      implicit none
      integer ni,nj,nk,i,j,k,jk,ni1,nj1
      real px(ni,nj,nk), pd(ni,nj,nk)                                 
c                                                                   
c  this routine determines p(.) from p(x) by a 4-point interpolation.  
c  on the x-grid, a p(x) point outside the grid domain is assumed to    
c  satisfy p(0,j)=p(1,j); p(ni,j)=p(ni-1,j); and similarly for the i's.
c                                                                      
      do 5 k=1,nk                                                     

       do j=2,nj1                                                     
       do i=2,ni1                                                     
        pd(i,j,k)=0.25*(px(i,j,k)+px(i-1,j,k)+px(i,j-1,k)+
     &            px(i-1,j-1,k))
       end do
       end do

       do 2 i=2,ni1                                                   
        pd(i,1,k)=0.5*(px(i,1,k)+px(i-1,1,k))                           
2       pd(i,nj,k)=0.5*(px(i,nj1,k)+px(i-1,nj1,k))                     

       do 3 j=2,nj1                                                   
        pd(1,j,k)=0.5*(px(1,j,k)+px(1,j-1,k))                           
3       pd(ni,j,k)=0.5*(px(ni1,j,k)+px(ni1,j-1,k))                     

       pd(1,1,k)=px(1,1,k)                                             
       pd(1,nj,k)=px(1,nj1,k)                                          
       pd(ni,1,k)=px(ni1,1,k)                                         
       pd(ni,nj,k)=px(ni1,nj1,k)                                   

5     continue                                                  

      return                                                    
      end

      SUBROUTINE FEXIST(filnam)
      implicit none

      character filnam*50, yesno*1
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

      SUBROUTINE FEXISTNEW(filnam,there)
      implicit none

      character filnam*50, yesno*1
      logical there
      
 1    inquire(file=filnam,exist=there)
      if (.not.there) then
 2      print*,'FILE CAN NOT BE OPENED BECAUSE IT DOES NOT EXISTS: '
     &        ,filnam
        print*,'DO YOU WANT TO CONTINUE? (y/n)'
        read(*,*) yesno
        if (yesno.eq.'y') then
          print*,'ENTER NEW FILE NAME'
          read(*,*) filnam
          go to 1
        elseif (yesno.eq.'n') then
          return
        else
          print*,'I DO NOT UNDERSTAND YOUR RESPONSE!!!'
          go to 2
        end if
      end if
      print*,'OPEN NEW FILE:',filnam

      return
      end

      SUBROUTINE RDHEADICBC(nxf,nyf,nx,ny,nz,clat,clon,ds
     &         , pt,sigf,sigh,sighrev,xplat,xplon,f,xmap,dmap
     &         , xlat,xlon,dlat,dlon,zs,zssd,ls
     &         , iin,inhead,ibyte)

      implicit none

      integer nxf,nyf,nx,ny,nz,ibyte
      real   sigf(nz+1), f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny)
     &     , xlon(nx,ny), dlat(nx,ny), dlon(nx,ny), zs(nx,ny)
     &     , zssd(nx,ny), ls(nx,ny), sigh(nz), sighrev(nz), grdfac
     &     , ds, clat, clon, pt, xplat, xplon, dto, dtb, dtr
     &     , tmp2d(nxf,nyf)
      integer ibltyp,icup,imoist,iboudy,igrads,ibigend
     &      , iin,ni,nj,nk,i,j,k,kk,ierr
      character proj*6, inhead*70

      open (iin,file=inhead,status='old',form='unformatted'
     &     ,recl=nxf*nyf*ibyte,access='direct')
      read(iin,rec=1,iostat=ierr) ni,nj,nk,ds,clat,clon,xplat,xplon
     &               ,grdfac,proj,sigf,pt,igrads,ibigend
      print*,'ni,nj,nk,ds='
      print*,ni,nj,nk,ds
      print*,'sigf='
      print*,sigf
      print*,'pt,clat,clon,xplat,xplon,proj='
      print*,pt,clat,clon,xplat,xplon,proj
      if (ni.ne.nyf .or. nj.ne.nxf .or. nz.ne.nk) then
        print*,'Grid Dimensions DO NOT MATCH'
        print*,'  nyf=',nyf,'nxf=',nxf,'nz=',nz
        print*,'  ni=',ni,'nj=',nj,'nk=',nk
        print*,'  Also check ibyte in icbc.param: ibyte= ',ibyte
        stop 'BAD DIMENSIONS (SUBROUTINE RDHEADICBC)'
      end if
c     print*,'ZS'
      read (iin,rec=2,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        zs(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'ZSSD'
      read (iin,rec=3,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        zssd(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'LU'
      read (iin,rec=4,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        ls(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'XLAT'
      read (iin,rec=5,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        xlat(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'XLON'
      read (iin,rec=6,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        xlon(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'XMAP'
      read (iin,rec=9,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        xmap(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'DMAP'
      read (iin,rec=10,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        dmap(i,j) = tmp2d(i+1,j+1)
      end do
      end do
c     print*,'F'
      read (iin,rec=11,iostat=ierr) tmp2d
      do j=1,ny
      do i=1,nx
        f(i,j) = tmp2d(i+1,j+1)
      end do
      end do

      if (ierr.ne.0) then
        print*,'END OF FILE REACHED'
        print*,'  Check ibyte in postproc.param: ibyte= ',ibyte
        stop 'EOF (SUBROUTINE RDHEADICBC)'
      end if

      do k=1,nz
        kk = nz - k + 1
        sigh(k) = (sigf(k)+sigf(k+1))/2.
        sighrev(kk) = sigh(k)
      end do

      return
      end

      SUBROUTINE RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev
     &         , xplat,xplon,f,xmap,dmap,xlat,xlon,dlat,dlon
     &         , zs,zssd,ls,mdate0,iin,inhead,idirect)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real   sigf(nz+1), f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny)
     &     , xlon(nx,ny), dlat(nx,ny), dlon(nx,ny), zs(nx,ny)
     &     , zssd(nx,ny), ls(nx,ny), sigh(nz), sighrev(nz), grdfac
     &     , ds, clat, clon, pt, xplat, xplon, dto, dtb, dtr, dtc
     &     , tmp2d(nxf,nyf)
      integer mdate0,ibltyp,icup,imoist,iboudy,igrads,ibigend
     &      , iin,ni,nj,nk,i,j,k,kk,ierr,idirect
      character proj*6, inhead*70
      logical there

      inquire(file=inhead,exist=there)
      if (.not.there) then
        print*,'OUT_HEAD FILE DOES NOT EXIST: ',inhead
        stop 'SUBROUTINE RDHEAD'
      end if
      open (iin,file=inhead,status='old',form='unformatted'
     &     ,recl=nx*ny*ibyte,access='direct')
      read(iin,rec=1,iostat=ierr) mdate0,ibltyp,icup,imoist,iboudy
     &,ni,nj,nk,(sigf(k),k=nz+1,1,-1),ds,pt,clat,clon,xplat,xplon,proj
     &,dto,dtb,dtr,dtc,idirect
      print*,'mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds='
      print*,mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds
      print*,'sigf='
      print*,sigf
      print*,'pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,dtc='
      print*,pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,dtc
      if (ni.ne.nyf .or. nj.ne.nxf .or. nz.ne.nk) then
        print*,'Grid Dimensions DO NOT MATCH'
        print*,'  nyf=',nyf,'nxf=',nxf,'nz=',nz
        print*,'  ni=',ni,'nj=',nj,'nk=',nk
        print*,'  Also check ibyte in postproc.param: ibyte= ',ibyte
        stop 'BAD DIMENSIONS (SUBROUTINE RDHEAD)'
      end if
      if (dto.ne.dtout .or. dtb.ne.dtbat .or. dtr.ne.dtrad
     &    .or. dtc.ne.dtche) then
      print*,'OUTPUT INTERVALS ARE IMPROPERLY DEFINED'
        print*,'dto=',dto,'dtout=',dtout
        print*,'dtb=',dtb,'dtbat=',dtbat
        print*,'dtr=',dtr,'dtrad=',dtrad
        print*,'dtc=',dtc,'dtche=',dtche
        stop 'BAD TIME PARAMETERS (SUBROUTINE RDHEAD)'
      end if
      print*,'Access type= ',idirect
c     print*,'ZS'
      read (iin,rec=2,iostat=ierr) zs
c     print*,'ZSSD'
      read (iin,rec=3,iostat=ierr) zssd
c     print*,'LS'
      read (iin,rec=4,iostat=ierr) ls
c     print*,'SATBRT'
      read (iin,rec=5,iostat=ierr) ls
c     print*,'XLAT'
      read (iin,rec=6,iostat=ierr) xlat
c     print*,'XLON'
      read (iin,rec=7,iostat=ierr) xlon
c     print*,'XMAP'
      read (iin,rec=8,iostat=ierr) xmap
c     print*,'DMAP'
      read (iin,rec=9,iostat=ierr) dmap
c     print*,'F'
      read (iin,rec=10,iostat=ierr) f

      if (ierr.ne.0) then
        print*,'END OF FILE REACHED'
        print*,'  Check ibyte in postproc.param: ibyte= ',ibyte
        stop 'EOF (SUBROUTINE RDHEAD)'
      end if

      do k=1,nz
        kk = nz - k + 1
        sigh(k) = (sigf(k)+sigf(k+1))/2.
        sighrev(kk) = sigh(k)
      end do

      return
      end

      SUBROUTINE RDICBC(vnambc,lnambc,ubc,idate,iin,irec,ierr)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer   nui,  nvi,  nti, nqvi, nmsei,  nrhi, nhgti
     &      , npsi, ntgi
      common /ipoint3d/  nui,  nvi,  nti, nqvi,nmsei, nrhi,nhgti
      common /ipoint2d/ npsi, ntgi
      real ifld2d, ifld3d, i2davg, i3davg
      common /bcflds/ ifld2d(nx,ny,nbc2d), ifld3d(nx,ny,nz,nbc3d)
     &   , i2davg(nx,ny,nbc2d,nhrbc), i3davg(nx,ny,nz,nbc3d,nhrbc)

      real tmp2d(nxf,nyf)
      integer idate, iin, ierr, ni, i,j,k, irec, kk
      character vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13


      print*,''
      ierr=0
      irec = irec + 1
      read (iin,rec=irec,iostat=ierr) idate
      if (ierr.ne.0) return
c     print *,'Reading ICBC:  ',idate
      do ni=1,ni3d
c       print*,'READ VAR:  ',vnambc(ni),lnambc(ni)
        do k=1,nz
          irec = irec + 1
          read (iin,rec=irec,iostat=ierr) tmp2d
          if (ierr.ne.0) return
          kk = nz - k + 1
          if (ni.eq.nui .or. ni.eq.nvi) then
            do j=1,ny
            do i=1,nx
              ifld3d(i,j,kk,ni) = 0.25*(tmp2d(i+2,j+2)+tmp2d(i+1,j+1)
     &                                 +tmp2d(i+1,j+2)+tmp2d(i+2,j+1))
            end do
            end do
          else
            do j=1,ny
            do i=1,nx
              ifld3d(i,j,kk,ni) = tmp2d(i+1,j+1)
c             ifld3d(i,j,k,ni) = tmp2d(i+1,j+1)
            end do
            end do
          end if
        end do
      end do
      do ni=1,ni2d
c       print*,'READ VAR:  ',vnambc(ni+nbc3d),lnambc(ni+nbc3d)
        irec = irec + 1
        read (iin,rec=irec,iostat=ierr) tmp2d
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          if (ni.eq.npsi) then
            ifld2d(i,j,ni) = tmp2d(i+1,j+1)*10.
          else
            ifld2d(i,j,ni) = tmp2d(i+1,j+1)
          end if
        end do
        end do
      end do
      print*,'ICBC DATA READ:',idate

      return
      end

      SUBROUTINE RDATM(vnamout,lnamout,uout,idate,iin,orec,idirect
     &         , ierr)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)

      real tmp2d(nx,ny)
      integer idate, iin, ierr, no, i,j,k,kk,orec,idirect
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13

      integer   nua,  nva, nomega,  nta, nqva, nqca,   nrh, nhgt
     &      , ntha, ntda,nvora,ndiva, npsa,  nrt,ntgb,nsmt,nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca, nrh,nhgt
     &                   ,ntha, ntda, nvora, ndiva
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf, nslp

      print*,' '
      ierr=0
      if (idirect.ne.1) then
        read (iin,iostat=ierr) idate
        if (ierr.ne.0) return
      end if
      print *,'Reading output:  ',idate
      do no=1,no3d
c       print*,'READ VAR:  ',vnamout(no),lnamout(no)
        do k=1,nz
          if (idirect.eq.1) then
            orec = orec + 1
            read (iin,rec=orec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
c         print*,vnamout(no),ierr
          if (ierr.ne.0) return
          kk = nz - k + 1
          do j=1,ny
          do i=1,nx
c           ofld3d(i,j,k,no) = tmp2d(i,j)
            ofld3d(i,j,kk,no) = tmp2d(i,j)
          end do
          end do
        end do
      end do
      do no=1,no2d
c       print*,'READ VAR:  ',vnamout(no+nout3d),lnamout(no+nout3d)
        if (idirect.eq.1) then
          orec = orec + 1
          read (iin,rec=orec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
c       print*,vnamout(no),ierr
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          ofld2d(i,j,no) = tmp2d(i,j)
        end do
        end do
      end do
c     print*,'DONE READING OUTPUT FOR CURRENT TIMESTEP',idate

      return
      end

      SUBROUTINE WRITEBC(vvarmin,vvarmax,vnambc,lnambc,ubc,xmin,xmax
     &        , fact,offset,idim,ndim,xlat1d,xlon1d,sighrev,vmisdat
     &        , idout,xhr,iobctyp,unit,nrec,plv,u_bc)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim, i, j, k, idout, ni, nni, iobctyp, unit, nrec
      real vmin, vmax, vmisdat, misdat
      real xmin(nitot), xmax(nitot), fact(nitot), offset(nitot)
      real sighrev(nz)
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), tmp3d_p(nx,ny,npl)
      real*8 xhr
      integer idim(ndim)
      character vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13

      integer   nui,  nvi,  nqvi, nrhi, nti, ntdi, nthi, nvori, ndivi
     &      , nhgti, npsi, ntgi, nslpi
      common /ipoint3d/  nui,  nvi,  nqvi, nrhi,nti, ntdi, nthi
     &      , nvori, ndivi, nhgti
      common /ipoint2d/ npsi, ntgi, nslpi
      real ifld2d, ifld3d, i2davg, i3davg
      common /bcflds/ ifld2d(nx,ny,nbc2d), ifld3d(nx,ny,nz,nbc3d)
     &   , i2davg(nx,ny,nbc2d,nhrbc), i3davg(nx,ny,nz,nbc3d,nhrbc)
	real ifld3d_p, i3davg_p
      common /bcflds_p/ ifld3d_p(nx,ny,npl,nbc3d)
     &   , i3davg_p(nx,ny,npl,nbc3d,nhrbc)
	logical plv
	integer u_bc(nitot)


c **** WRITE OUT 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do ni=1,nbc3d
	  if(u_bc(ni).eq.1) then
c       print*,ni,vnambc(ni)
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = max(ifld3d(i,j,k,ni),vmisdat)
        end do
        end do
        end do
        if (iobctyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(ni) .or. vmax.gt.xmax(ni)) then
            print*,'Values Out of Range:  FIELD=',vnambc(ni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
            stop 999
          end if
          misdat = xmin(ni)
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(ni),tmp3d,nx,ny,nz,idim,xhr
     &       , lnambc(ni),ubc(ni),fact(ni),offset(ni),vvarmin,vvarmax
     &       , xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	  end if
      end do
	else
	idim(3)=npl
	CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do ni=1,nbc3d
	  if (u_bc(ni).eq.1) then
c       print*,ni,vnambc(ni)
        do k=1,npl
        do j=1,ny1
        do i=1,nx1
          tmp3d_p(i,j,k) = max(ifld3d_p(i,j,k,ni),vmisdat)
        end do
        end do
        end do
        if (iobctyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(ni) .or. vmax.gt.xmax(ni)) then
            print*,'Values Out of Range:  FIELD=',vnambc(ni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
            stop 999
          end if
          misdat = xmin(ni)
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(ni),tmp3d_p,nx,ny,npl,idim,xhr
     &       , lnambc(ni),ubc(ni),fact(ni),offset(ni),vvarmin,vvarmax
     &       , xlat1d,xlon1d,plev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	  end if
      end do
	endif
c **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do ni=1,nbc2d
        nni = ni + nbc3d
	  if(u_bc(nni).eq.1) then
c       print*,ni,nni,vnambc(nni)
        do j=1,ny1
        do i=1,nx1
          tmp2d(i,j) = max(ifld2d(i,j,ni),vmisdat)
        end do
        end do
        if (iobctyp.eq.1) then
          misdat = xmin(ni)
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nni) .or. vmax.gt.xmax(nni)) then
            print*,'Values Out of Range:  FIELD=',vnambc(nni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nni)
            stop 999
          end if
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(nni),tmp2d,nx,ny,1,idim,xhr
     &       , lnambc(nni),ubc(nni),fact(nni),offset(nni)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	  end if
      end do

      return
      end



      SUBROUTINE WRITEOUT(vvarmin,vvarmax,vnamout,lnamout,uout,xmin,xmax
     &        , fact,offset,idim,ndim,xlat1d,xlon1d,sighrev,vmisdat
     &        , idout,xhr,iotyp,unit,nrec,plv,u_out)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
       real ofld3d_p, o3davg_p
      common /outflds_p/ ofld3d_p(nx,ny,npl,nout3d)
     &   , o3davg_p(nx,ny,npl,nout3d,nhrout)
      integer ndim, i, j, k, idout, no, nno, iotyp, unit, nrec
      real vmin, vmax, vmisdat, misdat
      real xmin(notot), xmax(notot), fact(notot), offset(notot)
      real sighrev(nz)
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real tmp2d(nx,ny), tmp3d(nx,ny,nz),tmp3d_p(nx,ny,npl)
      real*8 xhr
      integer idim(ndim)
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13

      integer   nua,  nva, nomega,  nta, nqva, nqca,   nrh, nhgt
     &      , ntha, ntda,nvora, ndiva, npsa, nrt,ntgb,nsmt,nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca, nrh,nhgt
     &                   ,ntha, ntda, nvora, ndiva
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf, nslp
      logical plv
      integer u_out(notot)


c **** WRITE OUT 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,no,vnamout(no)
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = max(ofld3d(i,j,k,no),vmisdat)
        end do
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(no) .or. vmax.gt.xmax(no)) then
            print*,'Values Out of Range:  FIELD=',vnamout(no)
            print*,'MINVAL=',vmin,'XMIN=',xmin(no)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
            stop 999
          end if
          misdat = xmin(no)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(no),tmp3d,nx,ny,nz,idim,xhr
     &       , lnamout(no),uout(no),fact(no),offset(no),vvarmin,vvarmax
     &       , xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	end if
      end do
      else
      idim(3)=npl
            CALL SETCONST(tmp3d,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,no,vnamout(no)
        do k=1,npl
        do j=1,ny1
        do i=1,nx1
          tmp3d_p(i,j,k) = max(ofld3d_p(i,j,k,no),vmisdat)
        end do
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(no) .or. vmax.gt.xmax(no)) then
            print*,'Values Out of Range:  FIELD=',vnamout(no)
            print*,'MINVAL=',vmin,'XMIN=',xmin(no)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
            stop 999
          end if
          misdat = xmin(no)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(no),tmp3d_p,nx,ny,npl,idim,xhr
     &       , lnamout(no),uout(no),fact(no),offset(no),vvarmin,vvarmax
     &       , xlat1d,xlon1d,plev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	end if
      end do
      end if

c **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do no=1,nout2d
        nno = no + nout3d
	if (u_out(nno) .eq. 1) then
c       print*,no,nno,vnamout(nno)
        do j=1,ny1
        do i=1,nx1
          tmp2d(i,j) = max(ofld2d(i,j,no),vmisdat)
        end do
        end do
        if (iotyp.eq.1) then
          misdat = xmin(no)
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nno) .or. vmax.gt.xmax(nno)) then
            print*,'Values Out of Range:  FIELD=',vnamout(nno)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nno)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nno)
            stop 999
          end if
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(nno),tmp2d,nx,ny,1,idim,xhr
     &       , lnamout(nno),uout(nno),fact(nno),offset(nno)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do

      return
      end


      SUBROUTINE RDSRF(vnambat,lnambat,ubat,idate,iin,brec,idirect
     &         , ierr)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real*4 b2d(nx,ny)
      character vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)

      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha

      integer ierr, iin, idate, i, j, nb, brec, idirect

      ierr=0
      if (idirect.ne.1) then
        read(iin,iostat=ierr) idate
        PRINT*,' READING SRF (Sequential):', idate
        if (ierr.ne.0) return
      else
        PRINT*,' READING SRF (GrADS):', idate
      end if
      do nb=1,nbat
c       print*,'FIELD= ',vnambat(nb)
        if (idirect.eq.1) then
          brec = brec + 1
          read(iin,rec=brec,iostat=ierr) b2d
        else
          read(iin,iostat=ierr) b2d
        end if
c       print*,'FIELD= ',vnambat(nb),b2d(1,1),ierr
        do j=1,ny
        do i=1,nx
          bfld2d(i,j,nb) = b2d(i,j)
        end do
        end do
      end do

      return
      end

      SUBROUTINE RDSUB(vnamsub,lnamsub,usub,idate,iin,srec,idirect
     &         , ierr)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real*4 s2d(nxsb,nysb)
      character vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
      integer ierr, iin, idate, i, j, ns, srec, idirect

      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 


      ierr = 0
      if (idirect.ne.1) then
        read(iin,iostat=ierr) idate
        if (ierr.ne.0) return
      end if
      PRINT*,' READING SUB:', idate
      do ns=1,nsub
c       print*,'FIELD= ',vnamsub(ns)
        if (idirect.eq.1) then
          srec = srec + 1
          read(iin,rec=srec,iostat=ierr) s2d
        else
          read(iin,iostat=ierr) s2d
        end if
c       print*,'FIELD= ',vnamsub(ns),s2d(15,85),ierr
        do j=1,nysb
        do i=1,nxsb
          sfld2d(i,j,ns) = s2d(i,j)
        end do
        end do
      end do

      return
      end

      SUBROUTINE WRITEBAT(vnambat,lnambat,ubat,xmin,xmax,fact,offset
     &         , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,vmisdat
     &         , xhr,ihr,idout,iotyp,unit,nrec,u_bat)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim, i, j, nb, ihr, idout, iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real tmp2d(nx,ny), xmin(nbat2), xmax(nbat2), fact(nbat2)
     &   , offset(nbat2), vmin, vmax, sig1(1)
      character vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
      real*8 xhr
      integer idim(ndim)
      real vmisdat, misdat

      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)

      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
	logical plv
	integer u_bat(nbat2)
      sig1(1)=1.
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do nb=1,nbat2
	if (u_bat(nb).eq.1) then
c        if (nb.ne.ntmin .and. nb.ne.ntmax) then
          do i=nb1,nxb
          do j=nb1,nyb
            tmp2d(i,j) = max(bfld2d(i,j,nb),vmisdat)
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nxb,nyb,1,vmin,vmax,vmisdat)
            if (vmin.lt.xmin(nb) .or. vmax.gt.xmax(nb)) then
              print*,'Values Out of Range:  FIELD=',vnambat(nb)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nb)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nb)
c             stop 999
            end if
            misdat = xmin(nb)
          elseif (iotyp.eq.2 .or. iotyp.eq.3) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idout,vnambat(nb),tmp2d,nx,ny,1,idim,xhr
     &         , lnambat(nb),ubat(nb),fact(nb),offset(nb)
     &         , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
          end if
        end if
      end do

      return
      end

      SUBROUTINE WRITESUB(vnamsub,lnamsub,usub,xmin,xmax,fact,offset
     &         , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,vmisdat
     &         , xhr,ihr,idout,iotyp,unit,nrec)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim, i, j, ns, ihr, idout, iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), xlat1d(nysb), xlon1d(nxsb)
      real tmp2d(nxsb,nysb), xmin(nsub2), xmax(nsub2), fact(nsub2)
     &   , offset(nsub2), vmin, vmax, sig1(1)
      character vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
      real*8 xhr
      integer idim(ndim)
      real vmisdat, misdat

      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 

      sig1(1)=1.0
      CALL SETCONST(tmp2d,vmisdat,nxsb,nysb,1,1,1,1,nxsb,1,nysb)
      do ns=1,nsub2
c        if (ns.ne.nstmin .and. ns.ne.nstmax) then
          do i=1,nxsb
          do j=1,nysb
            tmp2d(i,j) = max(sfld2d(i,j,ns),vmisdat)
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nxsb,nysb,1,vmin,vmax,vmisdat)
            if (vmin.lt.xmin(ns) .or. vmax.gt.xmax(ns)) then
              print*,'Values Out of Range:  FIELD=',vnamsub(ns)
              print*,'MINVAL=',vmin,'XMIN=',xmin(ns)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(ns)
c             stop 999
            end if
            misdat = xmin(ns)
          elseif (iotyp.eq.2 .or. iotyp.eq.3) then
            misdat = vmisdat
          end if
c         print*,vnamsub(ns),nrec+1
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idout,vnamsub(ns),tmp2d,nxsb,nysb,1,idim,xhr
     &         , lnamsub(ns),usub(ns),fact(ns),offset(ns)
     &         , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nxsb,nysb,1,nrec)
          end if
c        end if
      end do

      return
      end

      SUBROUTINE RDRAD(vnamrad,lnamrad,iin,idate,rrec,idirect,ierr)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      real tmp2d(nx,ny)
      integer iin,idate,ierr,nr,i,j,k,idatez,rrec,idirect

      character vnamrad(nrtot)*10, lnamrad(nrtot)*20

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp

      ierr=0
      if (idirect.ne.1) then
        read (iin,iostat=ierr) idate
        if (ierr.ne.0) return
      end if
      print *,'READING RADIATION DATA:  ',idate
      do nr=1,nr3d
c       print*,'READ VAR:  ',vnamrad(nr),lnamrad(nr)
        do k=1,nz
          if (idirect.eq.1) then
            rrec = rrec + 1
            read (iin,rec=rrec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
          if (ierr.ne.0) return
          do j=1,ny
          do i=1,nx
            rfld3d(i,j,k,nr) = tmp2d(i,j)
          end do
          end do
        end do
      end do
      do nr=1,nr2d
c       print*,'READ VAR:  ',vnamrad(nr+nr3d),lnamrad(nr+nr3d)
        if (idirect.eq.1) then
          rrec = rrec + 1
          read (iin,rec=rrec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          rfld2d(i,j,nr) = tmp2d(i,j)
        end do
        end do
      end do
c     print*,'DONE READING RADIATION FOR CURRENT TIMESTEP',idate

      return
      end

      SUBROUTINE RDCHE(iin,idate,crec,idirect,ierr)
      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real tmp2d(nx,ny)
      integer iin,idate,ierr,nc,i,j,k,idatez,mdate,crec,idirect

      character vnamche(nctot)*10, lnamche(nctot)*20

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)
c      print*,'idirect is ',idirect
c      print*,'idate in RDCHE is ',idate    

csr changed this line
      if (idirect.ne.1) then
        read (iin,iostat=ierr) idate

      end if
      if (ierr.ne.0) return
      print *,'READING CHEM-TRACER DATA:  ',idate
     
      do nc=1,nc3d
c       print*,'READ VAR: ',vnamche(nc),lnamche(nc)
        do k=1,nz
          if (idirect.eq.1) then
            crec = crec + 1
            read (iin,rec=crec,iostat=ierr) tmp2d
          else
            read (iin,iostat=ierr) tmp2d
          end if
          if (ierr.ne.0) return
          do j=1,ny
          do i=1,nx
            cfld3d(i,j,k,nc) = tmp2d(i,j)
          end do
          end do
        end do
      end do

      do nc=1,nc2d
        if (idirect.eq.1) then
          crec = crec + 1
          read (iin,rec=crec,iostat=ierr) tmp2d
        else
          read (iin,iostat=ierr) tmp2d
        end if
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          cfld2d(i,j,nc) = tmp2d(i,j)
        end do
        end do
      end do

     
      return
      end






      SUBROUTINE TWO2THREE(i,nx,ny,nz,in,out)

      implicit none
      integer nx, ny, nz, i, j, k
      real out(nx,ny,nz), in(ny,nz)

      do k=1,nz
      do j=1,ny
        out(i,j,k) = in(j,k)
      end do
      end do

      return
      end

      SUBROUTINE ONE2TWO(i,nx,ny,in,out)

      implicit none
      integer nx, ny, i, j
      real out(nx,ny), in(ny)

      do j=1,ny
        out(i,j) = in(j)
      end do

      return
      end

      SUBROUTINE WRITERAD(vnamrad,lnamrad,urad,xmin,xmax,fact,offset
     &        , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,sighrev
     &        , vmisdat,idout,xhr,iotyp,unit,nrec)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      real*8 xhr
      character vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
      integer i,j,k,nr,nnr,ndim,idout,iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real sighrev(nz), vmisdat, misdat
     &   , xmin(nrtot), xmax(nrtot), fact(nrtot), offset(nrtot)
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), vmin, vmax
      integer idim(ndim)

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp


c **** WRITE RAD 3-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nr3d
c       print*,nr,vnamrad(nr)
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = rfld3d(i,j,k,nr)
        end do
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nr) .or. vmax.gt.xmax(nr)) then
            print*,'Values Out of Range:  FIELD=',vnamrad(nr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamrad(nr),tmp3d,nx,ny,nz,idim,xhr
     &       , lnamrad(nr),urad(nr),fact(nr),offset(nr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
      end do

c **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do nr=1,nr2d
        nnr = nr + nr3d
c       print*,nr,nnr,vnamrad(nnr)
        do j=1,ny
        do i=2,nx1
          tmp2d(i,j) = rfld2d(i,j,nr)
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nnr) .or. vmax.gt.xmax(nnr)) then
            print*,'Values Out of Range:  FIELD=',vnamrad(nnr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nnr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nnr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamrad(nnr),tmp2d,nx,ny,1,idim,xhr
     &       , lnamrad(nnr),urad(nnr),fact(nnr),offset(nnr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
      end do


      return
      end


      SUBROUTINE WRITECHE(vnamche,lnamche,uche,xmin,xmax,fact,offset
     &        , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,sighrev
     &        , vmisdat,idout,xhr,iotyp,unit,nrec,u_che)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      real*8 xhr
      character vnamche(nctot)*10, lnamche(nctot)*20, uche(nctot)*13
      integer i,j,k,nc,nnc,ndim,idout,iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real sighrev(nz), vmisdat, misdat
     &   , xmin(nctot), xmax(nctot), fact(nctot), offset(nctot)
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), vmin, vmax
      integer idim(ndim)

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)
      integer u_che(nctot)


c **** WRITE RAD 3-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      
      do nc=1,nc3d
	if (u_che(nc).eq.1) then
c       print*,nr,vnamrad(nr)
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = cfld3d(i,j,k,nc)*1.E9
        end do
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nc) .or. vmax.gt.xmax(nc)) then
            print*,'Values Out of Range:  FIELD=',vnamche(nc)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nc)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nc)
            stop 999
          end if
          misdat = xmin(nc)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
      
          CALL WRITECDF(idout,vnamche(nc),tmp3d,nx,ny,nz,idim,xhr
     &       , lnamche(nc),uche(nc),fact(nc),offset(nc)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)

        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	end if
      end do

c **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do nc=1,nc2d
        nnc = nc + nc3d
	if (u_che(nnc).eq.1) then
c       print*,nr,nnr,vnamrad(nnr)
        do j=1,ny
        do i=2,nx1
          tmp2d(i,j) = cfld2d(i,j,nc)
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nnc) .or. vmax.gt.xmax(nnc)) then
            print*,'Values Out of Range:  FIELD=',vnamche(nnc)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nnc)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nnc)
            stop 999
          end if
          misdat = xmin(nc)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
         
          CALL WRITECDF(idout,vnamche(nnc),tmp2d,nx,ny,1,idim,xhr
     &       , lnamche(nnc),uche(nnc),fact(nnc),offset(nnc)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do


      return
      end




      SUBROUTINE AVGDATA3D(favg,f,n1,n2,n3,n4,n5,ihr,vmisdat)

      implicit none
      integer i,j,k,l,n1,n2,n3,n4,n5,ihr
      real favg(n1,n2,n3,n4,n5), f(n1,n2,n3,n4), vmisdat

      do l=1,n4
      do k=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,k,l).gt.vmisdat) then
          favg(i,j,k,l,ihr) = favg(i,j,k,l,ihr) + f(i,j,k,l)
        else
          favg(i,j,k,l,ihr) = vmisdat
        end if
      end do
      end do
      end do
      end do

      return
      end

      SUBROUTINE AVGDATA2D(favg,f,n1,n2,n3,n4,ihr,vmisdat)

      implicit none
      integer i,j,l,n1,n2,n3,n4,ihr
      real favg(n1,n2,n3,n4), f(n1,n2,n3), vmisdat

      do l=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,l).gt.vmisdat) then
          favg(i,j,l,ihr) = favg(i,j,l,ihr) + f(i,j,l)
        else
          favg(i,j,l,ihr) = vmisdat
        end if
      end do
      end do
      end do

      return
      end

      SUBROUTINE AVGDATABAT(ihr,vmisdat)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      real vmisdat,misdat
      integer ihr,i,j,nb

      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)

      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha

      if (vmisdat.gt.0.0) then
        misdat = -1.0*vmisdat
      else
        misdat = vmisdat
      end if
      do nb=1,nbat2
c        if ((nb.eq.ntmin.or.nb.eq.ntmax) .and. ihr.lt.nhrbat) then
c        else
          do j=1,ny
          do i=1,nx
            if (bfld2d(i,j,nb).gt.misdat) then
              b2davg(i,j,nb,ihr) = b2davg(i,j,nb,ihr) + bfld2d(i,j,nb)
            else
              b2davg(i,j,nb,ihr) = vmisdat
            end if
          end do
          end do
c       end if
      end do

      return
      end

      SUBROUTINE AVGDATASUB(ihr,vmisdat)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      real vmisdat,misdat
      integer ihr,i,j,nb

      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 

      if (vmisdat.gt.0.0) then
        misdat = -1.0*vmisdat
      else
        misdat = vmisdat
      end if
      do nb=1,nsub2
c        if ((nb.eq.nstmin.or.nb.eq.nstmax) .and. ihr.lt.nhrsub) then
c        else
          do j=1,nysb
          do i=1,nxsb
            if (sfld2d(i,j,nb).gt.misdat) then
              s2davg(i,j,nb,ihr) = s2davg(i,j,nb,ihr) + sfld2d(i,j,nb)
            else
              s2davg(i,j,nb,ihr) = vmisdat
            end if
          end do
          end do
c       end if
      end do

      return
      end

      SUBROUTINE AVGRAIN(favg,f,f1,nx,ny,nhr,ihr)

      implicit none
	integer nx, ny, nhr, i, j, ihr
      real favg(nx,ny,nhr), f(nx,ny), f1(nx,ny)

      do j=1,ny
      do i=1,nx
        favg(i,j,ihr) = favg(i,j,ihr) + (f(i,j)-f1(i,j))
        f1(i,j) = f(i,j)
      end do
      end do

      return
      end

      SUBROUTINE WRITEAVGBAT(vmisdat,vnambat,lnambat,ubat,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,xhr1,nbattime,idbat,iotyp,unit,nrec,u_bat)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim, i, j, nb, ihr, idbat, iotyp, unit, nrec
      real favgsum(nx,ny,nbat2), sig1(1)
     &   , xmin(nbat2), xmax(nbat2), fact(nbat2), offset(nbat2)
     &   , vvarmin(ndim), vvarmax(ndim), tmp2d(nx,ny), const
     &   , xlat1d(ny), xlon1d(nx), vmisdat, misdat, vmin, vmax, xntimes
      character vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
      real*8 xhr1 ,xhravg
      integer idim(ndim), nbattime(nhrbat)

      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)
      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
	integer u_bat(nbat2)

      idim(3) = 1
      sig1(1)=1.

      print*,'COMPUTING AVERAGE BAT FIELDS:',nbattime

      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      CALL SETCONST(favgsum,0.0,nx,ny,nbat2,1,1,1,nx,1,ny)
      do ihr=1,nhrbat
      do nb=1,nbat2
	if(u_bat(nb).eq.1) then
        const=1.0
c       if (nb.eq.ntmin .or. nb.eq.ntmax) then
c          xntimes = const/float(nbattime(ihr))
c        else
          xntimes = const/float(nhrbat*nbattime(ihr))
c        end if
        do j=nb1,nyb
        do i=nb1,nxb
         if (b2davg(i,j,nb,ihr).gt.vmisdat) then
          favgsum(i,j,nb) = favgsum(i,j,nb) + b2davg(i,j,nb,ihr)*xntimes
         else
          favgsum(i,j,nb) = vmisdat
         end if
        end do
        end do
	end if
        end do
      end do
      xhravg = xhr1
      do nb=1,nbat2
	if (u_bat(nb).eq.1) then
        if (xntimes.le.0) then
          print*,'NOTHING TO AVERAGE -- nbattime = 0'
          stop 999
        end if
        do j=nb1,nyb
        do i=nb1,nxb
          tmp2d(i,j) = favgsum(i,j,nb)
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nxb,nyb,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(nb) .or. vmax.gt.xmax(nb)) then
            print*,'Values Out of Range:  FIELD=',vnambat(nb)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nb)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nb)
            stop 999
          end if
          misdat = xmin(nb)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idbat,vnambat(nb),tmp2d,nx,ny,1,idim
     &       , xhravg,lnambat(nb),ubat(nb),fact(nb),offset(nb)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do

      return
      end

      SUBROUTINE WRITEAVGSUB(vmisdat,vnamsub,lnamsub,usub,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,xhr1,nsubtime,idsub,iotyp,unit,nrec)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim, i, j, ns, ihr, idsub, iotyp, unit, nrec
      real favgsum(nxsb,nysb,nsub2), sig1(1)
     &   , xmin(nsub2), xmax(nsub2), fact(nsub2), offset(nsub2)
     &   , vvarmin(ndim), vvarmax(ndim), tmp2d(nxsb,nysb), const
     &   , xlat1d(nysb), xlon1d(nxsb), vmisdat, misdat, vmin, vmax
     &   , xntimes
      character vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
      real*8 xhr1 ,xhravg
      integer idim(ndim), nsubtime(nhrsub)

      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 

      idim(3) = 1
      sig1(1) = 1.

      print*,'COMPUTING AVERAGE SUB FIELDS:',nsubtime

      CALL SETCONST(tmp2d,vmisdat,nxsb,nysb,1,1,1,1,nxsb,1,nysb)
      CALL SETCONST(favgsum,0.0,nxsb,nysb,nsub2,1,1,1,nxsb,1,nysb)
      do ihr=1,nhrsub
      do ns=1,nsub2
        const=1.0
c        if (ns.eq.nstmin .or. ns.eq.nstmax) then
c          xntimes = const/float(nsubtime(ihr))
c        else
          xntimes = const/float(nhrsub*nsubtime(ihr))
c        end if
        do j=1,nysb
        do i=1,nxsb
         if (s2davg(i,j,ns,ihr).gt.vmisdat) then
          favgsum(i,j,ns) = favgsum(i,j,ns) + s2davg(i,j,ns,ihr)*xntimes
         else
          favgsum(i,j,ns) = vmisdat
         end if
        end do
        end do
        end do
      end do
      xhravg = xhr1
      do ns=1,nsub2
        if (xntimes.le.0) then
          print*,'NOTHING TO AVERAGE -- nsubtime = 0'
          stop 999
        end if
        do j=1,nysb
        do i=1,nxsb
          tmp2d(i,j) = favgsum(i,j,ns)
        end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nxsb,nysb,1,vmin,vmax,vmisdat)
          if (vmin.lt.xmin(ns) .or. vmax.gt.xmax(ns)) then
            print*,'Values Out of Range:  FIELD=',vnamsub(ns)
            print*,'MINVAL=',vmin,'XMIN=',xmin(ns)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(ns)
            stop 999
          end if
          misdat = xmin(ns)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idsub,vnamsub(ns),tmp2d,nxsb,nysb,1,idim
     &       , xhravg,lnamsub(ns),usub(ns),fact(ns),offset(ns)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nxsb,nysb,1,nrec)
        end if
      end do

      return
      end

      SUBROUTINE WRITEDIURBAT(vmisdat,vnambat,lnambat,ubat,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,xhr1,nbattime,idbat,iotyp,unit,nrec,u_bat)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      character vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
      integer ndim, i, j, ihr, nb, idbat, iotyp, unit, nrec
      real xmin(nbat2), xmax(nbat2), fact(nbat2), offset(nbat2)
     &   , tmp2d(nx,ny), vvarmin(ndim), vvarmax(ndim)
     &   , xlat1d(ny), xlon1d(nx), const, xntimes
     &   , vmisdat, misdat, vmin, vmax, sig1(1)
      real*8 xhr1 ,xhravg
      integer idim(ndim), nbattime(nhrbat), u_bat(nbat2)

      real bfld2d, b2davg
      common /batflds/ bfld2d(nx,ny,nbat2), b2davg(nx,ny,nbat2,nhrbat)
      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha

      idim(3) = 1
      sig1(1) = 1.
      CALL SETCONST(tmp2d,vmisdat,nx,ny,1,1,1,1,nx,1,ny)
      do ihr=1,nhrbat
        xhravg = xhr1 + float(ihr-1)*dtbat
        if (nbattime(ihr).le.0) then
          print*,'NOTHING TO AVERAGE -- nbattime = 0'
          stop 999
        end if
        const=1.0
        do nb=1,nbat2
	if (u_bat(nb).eq.1) then
c          if (nb.eq.ntmin .or. nb.eq.ntmax) then
c            xntimes = const/float(nbattime(ihr))
c          else
            xntimes = const/float(nbattime(ihr))
c          end if
c          if ((nb.eq.ntmin.or.nb.eq.ntmax) .and. ihr.lt.nhrbat) then
            do j=nb1,nyb
            do i=nb1,nxb
              tmp2d(i,j) = vmisdat
            end do
            end do
c          else
            do j=nb1,nyb
            do i=nb1,nxb
              if (b2davg(i,j,nb,ihr).gt.vmisdat) then
                tmp2d(i,j) = b2davg(i,j,nb,ihr)*xntimes
              else
                tmp2d(i,j) = vmisdat
              end if
            end do
            end do
c          end if
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nxb,nyb,1,vmin,vmax,vmisdat)
            if (vmin.lt.xmin(nb) .or. vmax.gt.xmax(nb)) then
              print*,'Values Out of Range:  FIELD=',vnambat(nb)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nb)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nb)
              stop 999
            end if
            misdat = xmin(nb)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idbat,vnambat(nb),tmp2d,nx,ny,1,idim
     &         , xhravg,lnambat(nb),ubat(nb),fact(nb),offset(nb)
     &         , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
          end if
	  end if
        end do
      end do

      return
      end

      SUBROUTINE WRITEDIURSUB(vmisdat,vnamsub,lnamsub,usub,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,xhr1,nsubtime,idsub,iotyp,unit,nrec)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      character vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
      integer ndim, i, j, ihr, ns, idsub, iotyp, unit, nrec
      real xmin(nsub2), xmax(nsub2), fact(nsub2), offset(nsub2)
     &   , tmp2d(nxsb,nysb), vvarmin(ndim), vvarmax(ndim)
     &   , xlat1d(nysb), xlon1d(nxsb), const, xntimes
     &   , vmisdat, misdat, vmin, vmax, sig1(1)
      real*8 xhr1 ,xhravg
      integer idim(ndim), nsubtime(nhrbat)

      real sfld2d, s2davg
      common /subflds/ sfld2d(nxsb,nysb,nsub2)
     &               , s2davg(nxsb,nysb,nsub2,nhrsub)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha  

      idim(3) = 1
      sig1(1) = 1.
      CALL SETCONST(tmp2d,vmisdat,nxsb,nysb,1,1,1,1,nxsb,1,nysb)
      do ihr=1,nhrsub
        xhravg = xhr1 + float(ihr-1)*dtsub
        if (nsubtime(ihr).le.0) then
          print*,'NOTHING TO AVERAGE -- nsubtime = 0'
          stop 999
        end if
        const=1.0
        do ns=1,nsub2
c          if (ns.eq.nstmin .or. ns.eq.nstmax) then
c            xntimes = const/float(nsubtime(ihr))
c          else
            xntimes = const/float(nsubtime(ihr))
c          end if
c          if ((ns.eq.nstmin.or.ns.eq.nstmax) .and. ihr.lt.nhrsub) then
            do j=1,nysb
            do i=1,nxsb
              tmp2d(i,j) = vmisdat
            end do
            end do
c          else
            do j=1,nysb
            do i=1,nxsb
              if (s2davg(i,j,ns,ihr).gt.vmisdat) then
                tmp2d(i,j) = s2davg(i,j,ns,ihr)*xntimes
              else
                tmp2d(i,j) = vmisdat
              end if
            end do
            end do
c          end if
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nxsb,nysb,1,vmin,vmax,vmisdat)
            if (vmin.lt.xmin(ns) .or. vmax.gt.xmax(ns)) then
              print*,'Values Out of Range:  FIELD=',vnamsub(ns)
              print*,'MINVAL=',vmin,'XMIN=',xmin(ns)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(ns)
              stop 999
            end if
            misdat = xmin(ns)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idsub,vnamsub(ns),tmp2d,nxsb,nysb,1,idim
     &         , xhravg,lnamsub(ns),usub(ns),fact(ns),offset(ns)
     &         , vvarmin,vvarmax,xlat1d,xlon1d,sig1,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nxsb,nysb,1,nrec)
          end if
        end do
      end do

      return
      end

      SUBROUTINE WRITEAVGBC(sighrev,vnambc,lnambc,ubc,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &         , xhr1,nbctime,idout,vmisdat,iobctyp,unit,nrec,plv
     &         , u_bc)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'


      integer ndim, idim(ndim), nbctime(nhrbc), iobctyp, unit, nrec
     &     ,i,j,k,ihr,idout,ni,nni
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), xlat1d(ny), xlon1d(nx)
     &   , vvarmin(ndim),vvarmax(ndim),sighrev(nz),tmp3d_p(nx,ny,npl)
     &   , xmin(nitot), xmax(nitot), fact(nitot), offset(nitot)
     &   , xntimes, vmisdat, vmax, vmin, misdat
      character vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13
      real*8 xhr1 ,xhravg

      integer   nui,  nvi,  nqvi, nrhi, nti, ntdi, nthi, nvori, ndivi
     &      , nhgti, npsi, ntgi, nslpi
      common /ipoint3d/  nui,  nvi,  nqvi, nrhi,nti, ntdi, nthi
     &      , nvori, ndivi, nhgti
      common /ipoint2d/ npsi, ntgi, nslpi
      real ifld2d, ifld3d, i2davg, i3davg
      common /bcflds/ ifld2d(nx,ny,nbc2d), ifld3d(nx,ny,nz,nbc3d)
     &   , i2davg(nx,ny,nbc2d,nhrbc), i3davg(nx,ny,nz,nbc3d,nhrbc)
	real ifld3d_p, i3davg_p
      common /bcflds_p/ ifld3d_p(nx,ny,npl,nbc3d)
     &   , i3davg_p(nx,ny,npl,nbc3d,nhrbc)
	logical plv
	integer u_bc(nitot)

      print*,'COMPUTING AVERAGE ICBC FIELDS:',nbctime
      xhravg = xhr1
      print*,'xhravg=',xhravg
c **** WRITE ICBC AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
	idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do ni=1,nbc3d
	  if (u_bc(ni).eq.1) then
c       print*,vnambc(ni)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
        do ihr=1,nhrbc
          xntimes = 1./float(nbctime(ihr)*nhrbc)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (i3davg(i,j,k,ni,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + i3davg(i,j,k,ni,ihr)*xntimes
              i3davg(i,j,k,ni,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iobctyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(ni).or.vmax.gt.xmax(ni))then
            print*,'Values Out of Range:  FIELD=',vnambc(ni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
            stop 999
          end if
          misdat = xmin(ni)
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(ni),tmp3d,nx,ny,nz,idim,xhravg
     &       , lnambc(ni),ubc(ni),fact(ni),offset(ni)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	  end if
      end do
	else
	idim(3)=npl
	CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do ni=1,nbc3d
	  if(u_bc(ni).eq.1) then
c       print*,vnambc(ni)
        CALL SETCONST(tmp3d_p,0.0,nx,ny,npl,1,1,1,nx1,1,ny1)
        do ihr=1,nhrbc
          xntimes = 1./float(nbctime(ihr)*nhrbc)
          do k=1,npl
          do j=1,ny
          do i=1,nx
            if (i3davg_p(i,j,k,ni,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k) 
     &         + i3davg_p(i,j,k,ni,ihr)*xntimes
              i3davg_p(i,j,k,ni,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iobctyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(ni).or.vmax.gt.xmax(ni))then
            print*,'Values Out of Range:  FIELD=',vnambc(ni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
            stop 999
          end if
          misdat = xmin(ni)
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(ni),tmp3d_p,nx,ny,npl,idim,xhravg
     &       , lnambc(ni),ubc(ni),fact(ni),offset(ni)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,plev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	  end if
      end do
	endif
c **** WRITE ICBC AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do ni=1,nbc2d
        nni = nbc3d + ni
	  if (u_bc(nni).eq.1) then
c       print*,vnambc(nni)
        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrbc
          xntimes = 1./float(nbctime(ihr)*nhrbc)
          do j=1,ny
          do i=1,nx
            if (i2davg(i,j,ni,ihr).gt.vmisdat) then
              tmp2d(i,j) = tmp2d(i,j) + i2davg(i,j,ni,ihr)*xntimes
              i2davg(i,j,ni,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
        end do
        if (iobctyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nni).or.vmax.gt.xmax(nni))then
            print*,'Values Out of Range:  FIELD=',vnambc(nni)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nni)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nni)
            stop 999
          end if
          misdat = xmin(ni)
        elseif (iobctyp.eq.2) then
          misdat = vmisdat
        end if
        if (iobctyp.eq.1 .or. iobctyp.eq.2) then
          CALL WRITECDF(idout,vnambc(nni),tmp2d,nx,ny,1,idim,xhravg
     &       , lnambc(nni),ubc(nni),fact(nni),offset(nni)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
        else if (iobctyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	  end if
      end do

      return
      end


      SUBROUTINE WRITEAVGOUT(sighrev,vnamout,lnamout,uout,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &         , xhr1,nouttime,idout,vmisdat,iotyp,unit,nrec,plv,u_out)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      real ofld3d_p, o3davg_p
      common /outflds_p/ ofld3d_p(nx,ny,npl,nout3d)
     &   , o3davg_p(nx,ny,npl,nout3d,nhrout)

      integer ndim, idim(ndim), nouttime(nhrout), iotyp, unit, nrec
     &     ,i,j,k,ihr,idout,no,nno
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), xlat1d(ny), xlon1d(nx)
     &   , vvarmin(ndim), vvarmax(ndim), sighrev(nz)
     &   , xmin(notot), xmax(notot), fact(notot), offset(notot)
     &   , xntimes, vmisdat, vmax, vmin, misdat, tmp3d_p(nx,ny,npl)
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
      real*8 xhr1 ,xhravg

      integer   nua,  nva, nomega,  nta, nqva, nqca, nrh, nhgt
     &      , ntha, ntda, nvora,ndiva,npsa, ntgb,nsmt,nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca,  nrh,nhgt
     &                   ,ntha, ntda, nvora, ndiva
      common /opoint2d/ npsa,  ntgb, nsmt,  nbf, nslp
      logical plv
      integer u_out(notot)

      print*,'COMPUTING AVERAGE OUT FIELDS:',nouttime
      xhravg = xhr1
      print*,'xhravg=',xhravg
c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,vnamout(no)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (o3davg(i,j,k,no,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + o3davg(i,j,k,no,ihr)*xntimes
              o3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(no).or.vmax.gt.xmax(no))then
            print*,'Values Out of Range:  FIELD=',vnamout(no)
            print*,'MINVAL=',vmin,'XMIN=',xmin(no)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
            stop 999
          end if
          misdat = xmin(no)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(no),tmp3d,nx,ny,nz,idim,xhravg
     &       , lnamout(no),uout(no),fact(no),offset(no)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	end if
      end do
      else
      idim(3)=npl
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,vnamout(no)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)
          do k=1,npl
          do j=1,ny
          do i=1,nx
            if (o3davg_p(i,j,k,no,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k)
     &		+ o3davg_p(i,j,k,no,ihr)*xntimes
              o3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(no).or.vmax.gt.xmax(no))then
            print*,'Values Out of Range:  FIELD=',vnamout(no)
            print*,'MINVAL=',vmin,'XMIN=',xmin(no)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
            stop 999
          end if
          misdat = xmin(no)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(no),tmp3d_p,nx,ny,npl,idim,xhravg
     &       , lnamout(no),uout(no),fact(no),offset(no)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	end if
      end do
      end if

c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do no=1,nout2d
        nno = nout3d + no
	if (u_out(nno) .eq.1) then
c       print*,vnamout(nno)
        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)
          do j=1,ny
          do i=1,nx
            if (o2davg(i,j,no,ihr).gt.vmisdat) then
              tmp2d(i,j) = tmp2d(i,j) + o2davg(i,j,no,ihr)*xntimes
              o2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nno).or.vmax.gt.xmax(nno))then
            print*,'Values Out of Range:  FIELD=',vnamout(nno)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nno)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nno)
            stop 999
          end if
          misdat = xmin(no)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(nno),tmp2d,nx,ny,1,idim,xhravg
     &       , lnamout(nno),uout(nno),fact(nno),offset(nno)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do


      return
      end

      SUBROUTINE WRITEDIURBC(sighrev,vnambc,lnambc,ubc,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &         , xhr1,nbctime,idout,vmisdat,iobctyp,unit,nrec,plv
     &         , u_bc)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim,idim(ndim),nbctime(nhrbc),iobctyp, unit, nrec
     &      , idout,ni,nni,ihr,i,j,k
      real xmin(nitot), xmax(nitot), fact(nitot), offset(nitot)
     &   , vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
     &   , sighrev(nz), vmisdat, misdat, vmin, vmax, xntimes
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz),tmp3d_p(nx,ny,npl)
      character vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13
      real*8 xhr1 ,xhravg

      integer   nui,  nvi,  nqvi, nrhi, nti, ntdi, nthi, nvori, ndivi
     &      , nhgti, npsi, ntgi, nslpi
      common /ipoint3d/  nui,  nvi,  nqvi, nrhi,nti, ntdi, nthi
     &      , nvori, ndivi, nhgti
      common /ipoint2d/ npsi, ntgi,nslpi
      real ifld2d, ifld3d, i2davg, i3davg
      common /bcflds/ ifld2d(nx,ny,nbc2d), ifld3d(nx,ny,nz,nbc3d)
     &   , i2davg(nx,ny,nbc2d,nhrbc), i3davg(nx,ny,nz,nbc3d,nhrbc)
	real ifld3d_p, i3davg_p
      common /bcflds_p/ ifld3d_p(nx,ny,npl,nbc3d)
     &   , i3davg_p(nx,ny,npl,nbc3d,nhrbc)
	logical plv
	integer u_bc(nitot)


      print*,'COMPUTING AVERAGE FIELDS FOR DIURNAL OUTPUT:',nbctime

c **** WRITE ICBC AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do ni=1,nbc3d
	  if (u_bc(ni).eq.1) then
c       print*,vnambc(ni)
        do ihr=1,nhrbc
          xhravg = xhr1 + float(ihr-1)*dtbc
          xntimes = 1./float(nbctime(ihr))
c         xhravg = float(ihr-1)*dtbc
          if (nbctime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nbctime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (i3davg(i,j,k,ni,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = i3davg(i,j,k,ni,ihr)*xntimes
              i3davg(i,j,k,ni,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iobctyp.eq.1) then
            CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(ni).or.vmax.gt.xmax(ni))then
              print*,'Values Out of Range:  FIELD=',vnambc(ni)
              print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
              stop 999
            end if
            misdat = xmin(ni)
          elseif (iobctyp.eq.2) then
            misdat = vmisdat
          end if
          if (iobctyp.eq.1 .or. iobctyp.eq.2) then
            CALL WRITECDF(idout,vnambc(ni),tmp3d,nx,ny,nz,idim,xhravg
     &         , lnambc(ni),ubc(ni),fact(ni),offset(ni),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
          else if (iobctyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
          end if
        end do
      end if
      end do
      else
         idim(3)=npl
	CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do ni=1,nbc3d
         if (u_bc(ni).eq.1) then
c       print*,vnambc(ni)
        do ihr=1,nhrbc
          xhravg = xhr1 + float(ihr-1)*dtbc
          xntimes = 1./float(nbctime(ihr))
c         xhravg = float(ihr-1)*dtbc
          if (nbctime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nbctime = 0'
            stop 999
          end if
          do k=1,npl
          do j=1,ny1
          do i=1,nx1
            if (i3davg_p(i,j,k,ni,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = i3davg_p(i,j,k,ni,ihr)*xntimes
              i3davg_p(i,j,k,ni,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iobctyp.eq.1) then
            CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(ni).or.vmax.gt.xmax(ni))then
              print*,'Values Out of Range:  FIELD=',vnambc(ni)
              print*,'MINVAL=',vmin,'XMIN=',xmin(ni)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(ni)
              stop 999
            end if
            misdat = xmin(ni)
          elseif (iobctyp.eq.2) then
            misdat = vmisdat
          end if
          if (iobctyp.eq.1 .or. iobctyp.eq.2) then
            CALL WRITECDF(idout,vnambc(ni),tmp3d_p,nx,ny,npl,idim,xhravg
     &         , lnambc(ni),ubc(ni),fact(ni),offset(ni),vvarmin
     &         , vvarmax,xlat1d,xlon1d,plev,0,misdat,iobctyp)
          else if (iobctyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
          end if

       end do
      end if
      end do
      endif

c **** WRITE ICBC AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do ni=1,nbc2d
        nni = nbc3d + ni
	  if (u_bc(nni).eq.1) then
c       print*,vnambc(nni)
        do ihr=1,nhrbc
          xhravg = xhr1 + float(ihr-1)*dtbc
          xntimes = 1./float(nbctime(ihr))
          print*,'nbctime(ihr)=',nbctime(ihr),'xntimes=',xntimes
     &          ,'ihr=',ihr,xhravg
          if (nbctime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nbctime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (i2davg(i,j,ni,ihr).gt.vmisdat) then
              tmp2d(i,j) = i2davg(i,j,ni,ihr)*xntimes
              i2davg(i,j,ni,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
          if (iobctyp.eq.1) then
            CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nni).or.vmax.gt.xmax(nni))then
              print*,'Values Out of Range:  FIELD=',vnambc(nni)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nni)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nni)
              stop 999
            end if
            misdat = xmin(ni)
          elseif (iobctyp.eq.2) then
            misdat = vmisdat
          end if
          CALL WRITECDF(idout,vnambc(nni),tmp2d,nx,ny,1,idim
     &       , xhravg,lnambc(nni),ubc(nni),fact(nni),offset(nni)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iobctyp)
        end do
      end if
      end do

      return
      end


      SUBROUTINE WRITEDIUROUT(sighrev,vnamout,lnamout,uout,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &         , xhr1,nouttime,idout,vmisdat,iotyp,unit,nrec,plv,u_out)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      real ofld3d_p, o3davg_p
      common /outflds_p/ ofld3d_p(nx,ny,npl,nout3d)
     &   , o3davg_p(nx,ny,npl,nout3d,nhrout)
      integer ndim,idim(ndim),nouttime(nhrout),iotyp, unit, nrec
     &      , idout,no,nno,ihr,i,j,k
      real xmin(notot), xmax(notot), fact(notot), offset(notot)
     &   , vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
     &   , sighrev(nz), vmisdat, misdat, vmin, vmax, xntimes
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), tmp3d_p(nx,ny,npl)
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
      real*8 xhr1 ,xhravg

      integer   nua,  nva, nomega,  nta, nqva, nqca,  nrh, nhgt
     &      , ntha, ntda, nvora,ndiva,npsa,nrt,ntgb,nsmt, nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca,  nrh, nhgt
     &                   ,ntha, ntda, nvora, ndiva
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf, nslp
	logical plv
	integer u_out(notot)

      print*,'COMPUTING AVERAGE FIELDS FOR DIURNAL OUTPUT:',nouttime

c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,vnamout(no)
        do ihr=1,nhrout
          xhravg = xhr1 + float(ihr-1)*dtout
          xntimes = 1./float(nouttime(ihr))
c         xhravg = float(ihr-1)*dtout
          if (nouttime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (o3davg(i,j,k,no,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = o3davg(i,j,k,no,ihr)*xntimes
              o3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(no).or.vmax.gt.xmax(no))then
              print*,'Values Out of Range:  FIELD=',vnamout(no)
              print*,'MINVAL=',vmin,'XMIN=',xmin(no)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
              stop 999
            end if
            misdat = xmin(no)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idout,vnamout(no),tmp3d,nx,ny,nz,idim,xhravg
     &         , lnamout(no),uout(no),fact(no),offset(no),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
          end if
        end do
	end if
      end do
      else
      idim(3)=npl
            CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do no=1,nout3d
      if (u_out(no).eq.1) then
c       print*,vnamout(no)
        do ihr=1,nhrout
          xhravg = xhr1 + float(ihr-1)*dtout
          xntimes = 1./float(nouttime(ihr))
c         xhravg = float(ihr-1)*dtout
          if (nouttime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,npl
          do j=1,ny1
          do i=1,nx1
            if (o3davg_p(i,j,k,no,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = o3davg_p(i,j,k,no,ihr)*xntimes
              o3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(no).or.vmax.gt.xmax(no))then
              print*,'Values Out of Range:  FIELD=',vnamout(no)
              print*,'MINVAL=',vmin,'XMIN=',xmin(no)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(no)
              stop 999
            end if
            misdat = xmin(no)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idout,vnamout(no),tmp3d_p,nx,ny,npl,idim,xhravg
     &         , lnamout(no),uout(no),fact(no),offset(no),vvarmin
     &         , vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
          end if
        end do
	end if
      end do
	end if

c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do no=1,nout2d
        nno = nout3d + no
	if (u_out(nno) .eq.1) then
c       print*,vnamout(nno)
        do ihr=1,nhrout
          xhravg = xhr1 + float(ihr-1)*dtout
          xntimes = 1./float(nouttime(ihr))
          print*,'nouttime(ihr)=',nouttime(ihr),'xntimes=',xntimes
     &          ,'ihr=',ihr,xhravg
          if (nouttime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (o2davg(i,j,no,ihr).gt.vmisdat) then
              tmp2d(i,j) = o2davg(i,j,no,ihr)*xntimes
              o2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nno).or.vmax.gt.xmax(nno))then
              print*,'Values Out of Range:  FIELD=',vnamout(nno)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nno)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nno)
              stop 999
            end if
            misdat = xmin(no)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idout,vnamout(nno),tmp2d,nx,ny,1,idim,xhravg
     &         , lnamout(nno),uout(nno),fact(nno),offset(nno),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
          end if
        end do
	end if
      end do

      return
      end

      SUBROUTINE WRITEAVGRAD(xhr1,sighrev,vnamrad,lnamrad,urad,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,vmisdat,nradtime,idrad,iotyp,unit,nrec
     &         , plv,u_rad)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer i,j,k,nr,nnr,ndim,idrad,ihr,iotyp, unit, nrec
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), vmisdat, misdat, vmin, vmax
     &   , vvarmin(ndim), vvarmax(ndim), sighrev(nz), xntimes
     &   , xmin(nrtot), xmax(nrtot), fact(nrtot), offset(nrtot)
     &   , xlat1d(ny), xlon1d(nx), tmp3d_p(nx,ny,npl)
      character vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
      integer idim(ndim), nradtime(nhrrad)
      real*8 xhr1,xhravg

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      real rfld3d_p, r3davg_p
      common /radflds_p/ rfld3d_p(nx,ny,npl,nr3d)
     &   , r3davg_p(nx,ny,npl,nr3d,nhrrad)
      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp
	logical plv
	integer u_rad(nrtot)

      print*,'COMPUTING AVERAGE RAD FIELDS:',nradtime
      xhravg = xhr1
      print*,'nradtime=',nradtime
      print*,'xhravg=',xhravg

c **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nr3d
      if (u_rad(nr).eq.1) then
c       print*,vnamrad(nr)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (r3davg(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + r3davg(i,j,k,nr,ihr)*xntimes
              r3davg(i,j,k,nr,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
            print*,'Values Out of Range:  FIELD=',vnamrad(nr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idrad,vnamrad(nr),tmp3d,nx,ny,nz,idim,xhravg
     &       , lnamrad(nr),urad(nr),fact(nr),offset(nr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	end if
      end do
      else
      idim(3)=npl
            CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nr3d
      if (u_rad(nr).eq.1) then
c       print*,vnamrad(nr)
        CALL SETCONST(tmp3d_p,0.0,nx,ny,npl,1,1,1,nx1,1,ny1)
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do k=1,npl
          do j=1,ny
          do i=1,nx
            if (r3davg_p(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k)
     &		+ r3davg_p(i,j,k,nr,ihr)*xntimes
              r3davg_p(i,j,k,nr,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
            print*,'Values Out of Range:  FIELD=',vnamrad(nr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idrad,vnamrad(nr),tmp3d_p,nx,ny,npl,idim,xhravg
     &       , lnamrad(nr),urad(nr),fact(nr),offset(nr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	end if
      end do
      end if

c **** WRITE RAD AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do nr=1,nr2d
        nnr = nr3d + nr
	if (u_rad(nnr).eq. 1) then
c       print*,vnamrad(nnr)
        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do j=1,ny
          do i=1,nx
            if (r2davg(i,j,nr,ihr).gt.vmisdat) then
              tmp2d(i,j) = tmp2d(i,j) + r2davg(i,j,nr,ihr)*xntimes
              r2davg(i,j,nr,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nnr).or.vmax.gt.xmax(nnr))then
            print*,'Values Out of Range:  FIELD=',vnamrad(nnr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nnr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nnr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idrad,vnamrad(nnr),tmp2d,nx,ny,1,idim,xhravg
     &       , lnamrad(nnr),urad(nnr),fact(nnr),offset(nnr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do

      return
      end
      SUBROUTINE WRITEAVGCHE(xhr1,sighrev,vnamche,lnamche,uche,xmin,xmax
     &         , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &         , idim,ndim,vmisdat,nchetime,idche,iotyp,unit,nrec
     &         , plv,u_che)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      integer i,j,k,nr,nnr,ndim,idche,ihr,iotyp, unit, nrec
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), vmisdat, misdat, vmin, vmax
     &   , vvarmin(ndim), vvarmax(ndim), sighrev(nz), xntimes
     &   , xmin(nctot), xmax(nctot), fact(nctot), offset(nctot)
     &   , xlat1d(ny), xlon1d(nx), tmp3d_p(nx,ny,npl)
      character vnamche(nctot)*10, lnamche(nctot)*20, uche(nctot)*13
      integer idim(ndim), nchetime(nhrche)
      real*8 xhr1,xhravg

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)
      real cfld3d_p, c3davg_p
      common /cheflds_p/ cfld3d_p(nx,ny,npl,nc3d)
     &   , c3davg_p(nx,ny,npl,nc3d,nhrche)
      logical plv
      integer u_che(nctot)


      print*,'COMPUTING AVERAGE CHE FIELDS:',nchetime
      xhravg = xhr1
      print*,'nchetime=',nchetime
      print*,'xhravg=',xhravg

c **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nc3d
      if (u_che(nr).eq.1) then
c       print*,vnamrad(nr)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (c3davg(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + c3davg(i,j,k,nr,ihr)*xntimes
              c3davg(i,j,k,nr,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
            print*,'Values Out of Range:  FIELD=',vnamche(nr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idche,vnamche(nr),tmp3d,nx,ny,nz,idim,xhravg
     &       , lnamche(nr),uche(nr),fact(nr),offset(nr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
        end if
	end if
      end do
      else
      idim(3)=npl
            CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nc3d
      if (u_che(nr).eq.1) then
c       print*,vnamrad(nr)
        CALL SETCONST(tmp3d_p,0.0,nx,ny,npl,1,1,1,nx1,1,ny1)
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do k=1,npl
          do j=1,ny
          do i=1,nx
            if (c3davg_p(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k)
     &		+ c3davg_p(i,j,k,nr,ihr)*xntimes
              c3davg_p(i,j,k,nr,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
            print*,'Values Out of Range:  FIELD=',vnamche(nr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idche,vnamche(nr),tmp3d_p,nx,ny,npl,idim,xhravg
     &       , lnamche(nr),uche(nr),fact(nr),offset(nr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
        end if
	end if
      end do
      end if
      print*,'repere1'
c **** WRITE RAD AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do nr=1,nc2d
        nnr = nc3d + nr
	if (u_che(nnr).eq.1) then
c       print*,vnamrad(nnr)
        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do j=1,ny
          do i=1,nx
            if (c2davg(i,j,nr,ihr).gt.vmisdat) then
              tmp2d(i,j) = tmp2d(i,j) + c2davg(i,j,nr,ihr)*xntimes
              c2davg(i,j,nr,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
        end do
        if (iotyp.eq.1) then
          CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
          if(vmin.lt.xmin(nnr).or.vmax.gt.xmax(nnr))then
            print*,'Values Out of Range:  FIELD=',vnamche(nnr)
            print*,'MINVAL=',vmin,'XMIN=',xmin(nnr)
            print*,'MAXVAL=',vmax,'XMAX=',xmax(nnr)
            stop 999
          end if
          misdat = xmin(nr)
        elseif (iotyp.eq.2) then
          misdat = vmisdat
        end if
        if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idche,vnamche(nnr),tmp2d,nx,ny,1,idim,xhravg
     &       , lnamche(nnr),uche(nnr),fact(nnr),offset(nnr)
     &       , vvarmin,vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
        else if (iotyp.eq.3) then
          CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
        end if
	end if
      end do

      return
      end

      SUBROUTINE WRITEDIURRAD(xhr1,sighrev,vnamrad,lnamrad,urad,xmin
     &        , xmax,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &        , idim,ndim,vmisdat,nradtime,idrad,iotyp,unit,nrec
     &        , plv, u_rad)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      integer ndim,idrad,nr,nnr,ihr,i,j,k,iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), sighrev(nz), vmisdat, misdat
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), xntimes, vmin, vmax
     &   , xlat1d(ny), xlon1d(nx), tmp3d_p(nx,ny,npl)
      integer idim(ndim), nradtime(nhrrad)
      real*8 xhr1,xhravg
      character vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
      real xmin(nrtot), xmax(nrtot), fact(nrtot), offset(nrtot)

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      real rfld3d_p, r3davg_p
      common /radflds_p/ rfld3d_p(nx,ny,npl,nr3d)
     &   , r3davg_p(nx,ny,npl,nr3d,nhrrad)
      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp
      logical plv
      integer u_rad(nrtot)

c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not.plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nr3d
      if (u_rad(nr).eq.1) then
c       print*,vnamrad(nr)
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))
c         xhravg = float(ihr-1)*dtrad
          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nradtime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (r3davg(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = r3davg(i,j,k,nr,ihr)*xntimes
              r3davg(i,j,k,nr,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
              print*,'Values Out of Range:  FIELD=',vnamrad(nr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idrad,vnamrad(nr),tmp3d,nx,ny,nz,idim,xhravg
     &         , lnamrad(nr),urad(nr),fact(nr),offset(nr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
          end if
        end do
	end if
      end do
      else
      idim(3)=npl
            CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do nr=1,nr3d
      if (u_rad(nr).eq.1) then
c       print*,vnamrad(nr)
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))
c         xhravg = float(ihr-1)*dtrad
          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nradtime = 0'
            stop 999
          end if
          do k=1,npl
          do j=1,ny1
          do i=1,nx1
            if (r3davg(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = r3davg_P(i,j,k,nr,ihr)*xntimes
              r3davg_p(i,j,k,nr,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
              print*,'Values Out of Range:  FIELD=',vnamrad(nr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
          CALL WRITECDF(idrad,vnamrad(nr),tmp3d_p,nx,ny,npl,idim,xhravg
     &         , lnamrad(nr),urad(nr),fact(nr),offset(nr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
          end if
        end do
	end if
      end do
      end if


c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do nr=1,nr2d
        nnr = nr3d + nr
	if (u_rad(nnr).eq.1) then
c       print*,vnamrad(nnr)
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))
          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nradtime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (r2davg(i,j,nr,ihr).gt.vmisdat) then
              tmp2d(i,j) = r2davg(i,j,nr,ihr)*xntimes
              r2davg(i,j,nr,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nnr).or.vmax.gt.xmax(nnr))then
              print*,'Values Out of Range:  FIELD=',vnamrad(nnr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nnr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nnr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idrad,vnamrad(nnr),tmp2d,nx,ny,1,idim,xhravg
     &         , lnamrad(nnr),urad(nnr),fact(nnr),offset(nnr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
          end if
        end do
	end if
      end do


      return
      end


      SUBROUTINE WRITEDIURCHE(xhr1,sighrev,vnamche,lnamche,uche,xmin
     &        , xmax,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &        , idim,ndim,vmisdat,nchetime,idche,iotyp,unit,nrec
     &        , plv, u_che)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      integer ndim,idche,nr,nnr,ihr,i,j,k,iotyp, unit, nrec
      real vvarmin(ndim), vvarmax(ndim), sighrev(nz), vmisdat, misdat
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), xntimes, vmin, vmax
     &   , xlat1d(ny), xlon1d(nx), tmp3d_p(nx,ny,npl)
      integer idim(ndim), nchetime(nhrche)
      real*8 xhr1,xhravg
      character vnamche(nctot)*10, lnamche(nctot)*20, uche(nctot)*13
      real xmin(nctot), xmax(nctot), fact(nctot), offset(nctot)

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)
      real cfld3d_p, c3davg_p
      common /cheflds_p/ cfld3d_p(nx,ny,npl,nc3d)
     &   , c3davg_p(nx,ny,npl,nc3d,nhrche)
	logical plv
	integer u_che(nctot)

c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
	if (.not. plv) then
      idim(3) = nz
      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
      do nr=1,nc3d
      if (u_che(nr).eq.1) then
       print*,vnamche(nr)
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))
c         xhravg = float(ihr-1)*dtche
          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nchetime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (c3davg(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d(i,j,k) = c3davg(i,j,k,nr,ihr)*xntimes
              c3davg(i,j,k,nr,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d,nx,ny,nz,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
              print*,'Values Out of Range:  FIELD=',vnamche(nr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idche,vnamche(nr),tmp3d,nx,ny,nz,idim,xhravg
     &         , lnamche(nr),uche(nr),fact(nr),offset(nr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d,nx,ny,nz,nrec)
          end if
        end do
	end if
      end do
      else
      idim(3)= npl
            CALL SETCONST(tmp3d_p,vmisdat,nx,ny,npl,1,1,1,nx,1,ny)
      do nr=1,nc3d
      if (u_che(nr).eq.1) then
       print*,vnamche(nr)
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))
c         xhravg = float(ihr-1)*dtche
          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nchetime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (c3davg_p(i,j,k,nr,ihr).gt.vmisdat) then
              tmp3d_p(i,j,k) = c3davg_p(i,j,k,nr,ihr)*xntimes
              c3davg_p(i,j,k,nr,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdat
            end if
          end do
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp3d_p,nx,ny,npl,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nr).or.vmax.gt.xmax(nr))then
              print*,'Values Out of Range:  FIELD=',vnamche(nr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
           CALL WRITECDF(idche,vnamche(nr),tmp3d_p,nx,ny,npl,idim,xhravg
     &         , lnamche(nr),uche(nr),fact(nr),offset(nr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,plev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp3d_p,nx,ny,npl,nrec)
          end if
        end do
	end if
      end do
      end if


c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      idim(3) = 1
      do nr=1,nc2d
        nnr = nc3d + nr
	if (u_che(nnr).eq.1) then
c       print*,vnamche(nnr)
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))
          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nchetime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (c2davg(i,j,nr,ihr).gt.vmisdat) then
              tmp2d(i,j) = c2davg(i,j,nr,ihr)*xntimes
              c2davg(i,j,nr,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdat
            end if
          end do
          end do
          if (iotyp.eq.1) then
            CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
            if(vmin.lt.xmin(nnr).or.vmax.gt.xmax(nnr))then
              print*,'Values Out of Range:  FIELD=',vnamche(nnr)
              print*,'MINVAL=',vmin,'XMIN=',xmin(nnr)
              print*,'MAXVAL=',vmax,'XMAX=',xmax(nnr)
              stop 999
            end if
            misdat = xmin(nr)
          elseif (iotyp.eq.2) then
            misdat = vmisdat
          end if
          if (iotyp.eq.1 .or. iotyp.eq.2) then
            CALL WRITECDF(idche,vnamche(nnr),tmp2d,nx,ny,1,idim,xhravg
     &         , lnamche(nnr),uche(nnr),fact(nnr),offset(nnr),vvarmin
     &         , vvarmax,xlat1d,xlon1d,sighrev,0,misdat,iotyp)
          else if (iotyp.eq.3) then
            CALL WRITEGRADS(unit,tmp2d,nx,ny,1,nrec)
          end if
        end do
	end if
      end do


      return
      end

      SUBROUTINE MMVLUBC(vnambc,lnambc,ubc,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      character vnambc(nitot)*10, lnambc(nitot)*20, ubc(nitot)*13
      real xmin(nitot), xmax(nitot), fact(nitot), offset(nitot), aaa
      integer l

      integer   nui,  nvi,  nqvi, nrhi, nti, ntdi, nthi, nvori, ndivi
     &      , nhgti, npsi, ntgi, nslpi
      common /ipoint3d/  nui,  nvi,  nqvi, nrhi,nti, ntdi, nthi
     &      , nvori, ndivi, nhgti
      common /ipoint2d/ npsi, ntgi, nslpi


      vnambc(nui)  = 'U'
      vnambc(nvi)  = 'V'
      vnambc(nti)  = 'TK'
      vnambc(nqvi) = 'QD'
      vnambc(nbc3d+npsi) = 'PS'
      vnambc(nbc3d+ntgi) = 'TGRND'
      vnambc(nbc3d+nslpi) = 'SLP'
      vnambc(nrhi)  = 'RH'
      vnambc(nhgti) = 'HGT'
      vnambc(ntdi) =  'TD'
      vnambc(nthi)  = 'TH'
      vnambc(nvori) = 'VOR'
      vnambc(ndivi) = 'DIV'

      lnambc(nui)  = 'Zonal Wind'
      lnambc(nvi)  = 'Meridional Wind'
      lnambc(nti)  = 'Temperature'
      lnambc(nqvi) = 'Mixing Ratio'
      lnambc(nbc3d+npsi) = 'Surface Pressure'
      lnambc(nbc3d+ntgi) = 'Surface Temperature'
      lnambc(nbc3d+nslpi) = 'Sea Level Temperature'
      lnambc(nrhi)  = 'Relative Humidity'
      lnambc(nhgti) = 'Geopotential Height'
      lnambc(ntdi)  = 'Dew Point Temperature'
      lnambc(nthi)  = 'Potential Temperature'
      lnambc(nvori) = 'Vorticity (Vertical Component)'
      lnambc(nvori) = 'Vorticity (Horizontal Compnent)'

      ubc(nui)  = 'm/s'
      ubc(nvi)  = 'm/s'
      ubc(nti)  = 'K'
      ubc(nqvi) = 'kg/kg'
      ubc(nbc3d+npsi) = 'hPa'
      ubc(nbc3d+ntgi) = 'K'
      ubc(nbc3d+nslpi) = 'hPa'
      ubc(nrhi)  = 'fraction'
      ubc(nhgti) = 'm'
      ubc(ntdi)  = 'K'
      ubc(nthi)  = 'K'
      ubc(nvori) = 'm/s'
      ubc(ndivi)   = 'm/s'

      xmax(nui)  = 210.0
      xmax(nvi)  = 210.0
      xmax(nti)  = 350.0
      xmax(nqvi) = 0.1
      xmax(nbc3d+npsi) = 1200.0
      xmax(nbc3d+ntgi) = 350.0
      xmax(nbc3d+nslpi) = 1200.0
      xmax(nrhi)  = 30.0
      xmax(nhgti) = 40000.0
      xmax(ntdi)  = 350.0
      xmax(nthi)  = 350.0
      xmax(nvori) = 210.0
      xmax(ndivi) = 210.0

      xmin(nui)  = -210.0
      xmin(nvi)  = -210.0
      xmin(nti)  = 160.0
      xmin(nqvi) = -0.001
      xmin(nbc3d+npsi) = 300.0
      xmin(nbc3d+ntgi) = 200.0
      xmin(nbc3d+nslpi) = 200.0
      xmin(nrhi)  = -0.5
      xmin(nhgti) = -100.0
      xmin(ntdi)  = 160.0
      xmin(nthi)  = 160.0
      xmin(nvori) = -210.0
      xmin(ndivi) = -210.0


      aaa = 2.**16.-1.
      do l=1,nitot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end


      SUBROUTINE MMVLUOUT(vnamout,lnamout,uout,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
      real xmin(notot), xmax(notot), fact(notot), offset(notot), aaa
      integer l

      integer   nua,  nva, nomega,  nta, nqva, nqca,   nrh, nhgt
     &      , ntha, ntda, nvora, ndiva,npsa, nrt, ntgb, nsmt, nbf,nslp
      common /opoint3d/  nua,  nva,  nomega,nta, nqva, nqca,  nrh, nhgt
     &                   ,ntha, ntda, nvora, ndiva
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf, nslp

      vnamout(nua)  = 'U'
      vnamout(nva)  = 'V'
      vnamout(nomega) = 'OMEGA'
      vnamout(nta)  = 'TK'
      vnamout(nqva) = 'QD'
      vnamout(nqca) = 'QC'
      vnamout(nout3d+npsa) = 'PS'
      vnamout(nout3d+nrt)  = 'RT'
      vnamout(nout3d+ntgb) = 'TGRND'
      vnamout(nout3d+nsmt) = 'SMT'
      vnamout(nout3d+nbf)  = 'RB'
      vnamout(nout3d+nslp) = 'SLP'
      vnamout(nrh)  = 'RH'
      vnamout(nhgt) = 'HGT'
      vnamout(ntda) =  'TD'
      vnamout(ntha)  = 'TH'
      vnamout(nvora) = 'VOR'
      vnamout(ndiva) = 'DIV'

      lnamout(nua)  = 'Zonal Wind'
      lnamout(nva)  = 'Meridional Wind'
      lnamout(nomega) = 'Omega'
      lnamout(nta)  = 'Temperature'
      lnamout(nqva) = 'Mixing Ratio'
      lnamout(nqca) = 'Cloud Mixing Ratio'
      lnamout(nout3d+npsa) = 'Surface Pressure'
      lnamout(nout3d+ntgb) = 'Ground Temperature'
      lnamout(nout3d+nrt)  = 'Total Precip'
      lnamout(nout3d+nsmt) = 'Total Soil Water'
      lnamout(nout3d+nbf)  = 'Base Flow'
      lnamout(nout3d+nslp) = 'Sea Level Temperature'
      lnamout(nrh)  = 'Relative Humidity'
      lnamout(nhgt) = 'Geopotential Height'
      lnamout(ntda)  = 'Dew Point Temperature'
      lnamout(ntha)  = 'Potential Temperature'
      lnamout(nvora) = 'Vorticity (Vertical Component)'
      lnamout(nvora) = 'Vorticity (Horizontal Compnent)'

      uout(nua)  = 'm/s'
      uout(nva)  = 'm/s'
      uout(nomega) = 'hPa'
      uout(nta)  = 'K'
      uout(nqva) = 'kg/kg'
      uout(nqca) = 'kg/kg'
      uout(nout3d+npsa) = 'hPa'
      uout(nout3d+ntgb) = 'K'
      uout(nout3d+nrt)  = 'mm/day'
      uout(nout3d+nsmt) = 'mm'
      uout(nout3d+nbf)  = 'mm/day'
      uout(nout3d+nslp) = 'hPa'
      uout(nrh)  = 'fraction'
      uout(nhgt) = 'm'
      uout(ntda)  = 'K'
      uout(ntha)  = 'K'
      uout(nvora) = 'm/s'
      uout(ndiva)   = 'm/s'

      xmax(nua)  = 210.0
      xmax(nva)  = 210.0
      xmax(nomega) = 0.1
      xmax(nta)  = 350.0
      xmax(nqva) = 0.1
      xmax(nqca) = 0.1
      xmax(nout3d+npsa) = 1200.0
      xmax(nout3d+ntgb) = 350.0
      xmax(nout3d+nrt)  = 2500.0
      xmax(nout3d+nsmt) = 3000.0
      xmax(nout3d+nbf)  = 200.0
      xmax(nout3d+nslp) = 1200.0
      xmax(nrh)  = 30.0
      xmax(nhgt) = 40000.0
      xmax(ntda)  = 350.0
      xmax(ntha)  = 350.0
      xmax(nvora) = 210.0
      xmax(ndiva) = 210.0

      xmin(nua)  = -210.0
      xmin(nva)  = -210.0
      xmin(nomega) = -0.1
      xmin(nta)  = 160.0
      xmin(nqva) = -0.001
      xmin(nqca) = -0.001
      xmin(nout3d+npsa) = 300.0
      xmin(nout3d+ntgb) = 180.0
      xmin(nout3d+nrt)  = -10.0
      xmin(nout3d+nsmt) = 0.0
      xmin(nout3d+nbf)  = -10.0
      xmin(nout3d+nslp) = 200.0
      xmin(nrh)  = -0.5
      xmin(nhgt) = -100.0
      xmin(ntda)  = 160.0
      xmin(ntha)  = 160.0
      xmin(nvora) = -210.0
      xmin(ndiva) = -210.0

      aaa = 2.**16.-1.
      do l=1,notot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end

      SUBROUTINE MMVLUBAT(vnambat,lnambat,ubat,xmin,xmax
     &         , fact,offset,nbat2)

      implicit none
	integer l, nbat2
	real aaa
      character vnambat(nbat2)*10, lnambat(nbat2)*20, ubat(nbat2)*13
      real xmin(nbat2), xmax(nbat2), fact(nbat2), offset(nbat2)
      integer           nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              , ntgmax,  ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha
      common /bpoint/   nux,   nvx, ndrag,   ntg,   ntf, ntanm, nqanm
     &              ,  nsmu,  nsmr,   npt,   net, nrnfs, nsnow,   nsh
     &              ,  nlwn,  nswn,  nlwd,  nswi,  nprc, npsrf, nzpbl
     &              ,  ntgmax, ntgmin, ntamax, ntamin, w10max, psmin
     &              , nrha

      lnambat(nux) = 'Anemom Zonal Winds'
      vnambat(nux) = 'UA'
      ubat(nux) =  'm/s'
      xmax(nux) = 50.0
      xmin(nux) = -50.0

      lnambat(nvx) = 'Anemom Merid Winds'
      vnambat(nvx) = 'VA'
      ubat(nvx) =  'm/s'
      xmax(nvx) = 50.0
      xmin(nvx) = -50.0

      lnambat(ndrag) = 'Surface Drag Stress'
      vnambat(ndrag) =  'DRAG'
      ubat(ndrag) =  'si'
      xmax(ndrag) = 1.0
      xmin(ndrag) = -1.0

      vnambat(ntg)   =  'TG'
      lnambat(ntg)   = 'Ground Temperature'
      ubat(ntg)   =  'K'
      xmax(ntg)   = 350.0
      xmin(ntg)   = 180.0

      vnambat(ntf)   =  'TF'
      lnambat(ntf)   = 'Foliage Temp'
      ubat(ntf)   =  'K'
      xmax(ntf)   = 350.0
      xmin(ntf)   = 180.0

      lnambat(ntanm) = 'Anemom Temp'
      vnambat(ntanm) =  'TA'
      ubat(ntanm) =  'K'
      xmax(ntanm) = 350.0
      xmin(ntanm) = 180.0

      lnambat(nqanm) = 'Anemom Spec Humidity'
      vnambat(nqanm) =  'QA'
      ubat(nqanm) =  'kg/kg'
      xmax(nqanm) = 0.20
      xmin(nqanm) = -1.0E-5

      lnambat(nsmu)  = 'Top Layer Soil Moist'
      vnambat(nsmu)  =  'SMU'
      ubat(nsmu)  =  'mm'
      xmax(nsmu)  = 80.0
      xmin(nsmu)  = -1.0

      lnambat(nsmr)  = 'Root Lay Soil Moist'
      vnambat(nsmr)  =  'SMR'
      ubat(nsmr)  =  'mm'
      xmax(nsmr)  = 1200.0
      xmin(nsmr)  = -1.0

      lnambat(net)   = 'Evapotranspiration'
      vnambat(net)   =  'ET'
      ubat(net)   =  'mm/day'
      xmax(net)   = 150.0
      xmin(net)   = -5.0

      lnambat(nrnfs) = 'Surface Runoff'
      vnambat(nrnfs) =  'RNFS'
      ubat(nrnfs) =  'mm/day'
      xmax(nrnfs) = 2000.0
      xmin(nrnfs) = -200.0

      lnambat(nsnow) = 'Snow Depth'
      vnambat(nsnow) =  'SNOW'
      ubat(nsnow) =  'mm H2O'
      xmax(nsnow) = 1000.0
      xmin(nsnow) = -1.0

      lnambat(nsh)   = 'Sensible Heat'
      vnambat(nsh)   =  'SH'
      ubat(nsh)   =  'W/m2'
      xmax(nsh)   = 1000.0
      xmin(nsh)   = -300.0

      lnambat(nlwn)  = 'Net Longwave'
      vnambat(nlwn)  =  'LWN'
      ubat(nlwn)  =  'W/m2'
      xmax(nlwn)  = 750.0
      xmin(nlwn)  = -300.0

      lnambat(nlwd)  = 'Downward Longwave'
      vnambat(nlwd)  =  'LWD'
      ubat(nlwd)  =  'W/m2'
      xmax(nlwd)  = 750.0
      xmin(nlwd)  = -300.0

      lnambat(nswn)  = 'Net Solar Absorbed'
      vnambat(nswn)  =  'SWN'
      ubat(nswn)  =  'W/m2'
      xmax(nswn)  = 1200.0
      xmin(nswn)  = -1.0

      lnambat(nswi)  = 'Solar Incident'
      vnambat(nswi)  =  'SWI'
      ubat(nswi)  =  'W/m2'
      xmax(nswi)  = 1400.0
      xmin(nswi)  = -1.0

      lnambat(nprc)  = 'Convective Precip'
      vnambat(nprc)  =  'RC'
      ubat(nprc)  =  'mm/day'
      xmax(nprc)  = 1500.0
      xmin(nprc)  = -1.0

      lnambat(npt) = 'Total Precipitation'
      vnambat(npt) =  'RT'
      ubat(npt) =  'mm/day'
      xmax(npt) = 2500.0
      xmin(npt) = -1.0

      lnambat(nzpbl) = 'PBL Height'
      vnambat(nzpbl) =  'ZPBL'
      ubat(nzpbl) =  'm'
      xmax(nzpbl) = 6000.0
      xmin(nzpbl) = -1.0

      lnambat(npsrf) = 'Surface Pressure'
      vnambat(npsrf) =  'PSRF'
      ubat(npsrf) =  'hPa'
      xmax(npsrf) = 1500.0
      xmin(npsrf) = 300.0


      lnambat(nrha)  = 'Relative Humidity'
      vnambat(nrha)  =  'RHA'
      ubat(nrha)  =  'fraction'
      xmax(nrha)  = 5.0
      xmin(nrha)  = -0.1

      lnambat(ntgmax) = 'Max Ground Temp'
      vnambat(ntgmax) =  'TGMAX'
      ubat(ntgmax) =  'K'
      xmax(ntgmax) = 350.0
      xmin(ntgmax) = 200.0

      lnambat(ntgmin) = 'Min Ground Temp'
      vnambat(ntgmin) =  'TGMIN'
      ubat(ntgmin) =  'K'
      xmax(ntgmin) = 350.0
      xmin(ntgmin) = 200.0

      lnambat(ntamax) = 'Max Anemom Temp'
      vnambat(ntamax) =  'TAMAX'
      ubat(ntamax) =  'K'
      xmax(ntamax) = 350.0
      xmin(ntamax) = 200.0

      lnambat(ntamin) = 'Min Anemom Temp'
      vnambat(ntamin) =  'TAMIN'
      ubat(ntamin) =  'K'
      xmax(ntamin) = 350.0
      xmin(ntamin) = 200.0
   
      lnambat(w10max) = 'Max 10m Wind Speed'
      vnambat(w10max) = 'W10MX'
      ubat(w10max) =  'm/s'
      xmax(w10max) = 500.0
      xmin(w10max) = -500.0
      
      lnambat(psmin) = 'Min Surface Pressure'
      vnambat(psmin) =  'PSMIN'
      ubat(psmin) =  'hPa'
      xmax(psmin) = 1500.0
      xmin(psmin) = 300.0

      aaa = 2.**16.-1.
      do l=1,nbat2
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end

      SUBROUTINE MMVLUSUB(vnamsub,lnamsub,usub,xmin,xmax
     &         , fact,offset,nsub2)

      implicit none
      integer l, nsub2
      real aaa
      character vnamsub(nsub2)*10, lnamsub(nsub2)*20, usub(nsub2)*13
      real xmin(nsub2), xmax(nsub2), fact(nsub2), offset(nsub2)
      integer           nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 
      common /spoint/   nsux,   nsvx, nsdrag,   nstg,   nstf, nstanm
     &              , nsqanm,  nssmu,  nssmr,   nspt,   nset, nsrnfs
     &              , nssnow,   nssh,  nsprc, nspsrf, nsrha 

      lnamsub(nsux) = 'Anemom Zonal Winds'
      vnamsub(nsux) = 'UA'
      usub(nsux) =  'm/s'
      xmax(nsux) = 50.0
      xmin(nsux) = -50.0

      lnamsub(nsvx) = 'Anemom Merid Winds'
      vnamsub(nsvx) = 'VA'
      usub(nsvx) =  'm/s'
      xmax(nsvx) = 50.0
      xmin(nsvx) = -50.0

      lnamsub(nsdrag) = 'Surface Drag Stress'
      vnamsub(nsdrag) =  'DRAG'
      usub(nsdrag) =  'si'
      xmax(nsdrag) = 1.0
      xmin(nsdrag) = -1.0

      vnamsub(nstg)   =  'TG'
      lnamsub(nstg)   = 'Ground Temperature'
      usub(nstg)   =  'K'
      xmax(nstg)   = 350.0
      xmin(nstg)   = 180.0

      vnamsub(nstf)   =  'TF'
      lnamsub(nstf)   = 'Foliage Temp'
      usub(nstf)   =  'K'
      xmax(nstf)   = 350.0
      xmin(nstf)   = 180.0

      lnamsub(nstanm) = 'Anemom Temp'
      vnamsub(nstanm) =  'TA'
      usub(nstanm) =  'K'
      xmax(nstanm) = 350.0
      xmin(nstanm) = 180.0

      lnamsub(nsqanm) = 'Anemom Spec Humidity'
      vnamsub(nsqanm) =  'QA'
      usub(nsqanm) =  'kg/kg'
      xmax(nsqanm) = 0.20
      xmin(nsqanm) = -1.0E-5

      lnamsub(nssmu)  = 'Top Layer Soil Moist'
      vnamsub(nssmu)  =  'SMU'
      usub(nssmu)  =  'mm'
      xmax(nssmu)  = 80.0
      xmin(nssmu)  = -1.0

      lnamsub(nssmr)  = 'Root Lay Soil Moist'
      vnamsub(nssmr)  =  'SMR'
      usub(nssmr)  =  'mm'
      xmax(nssmr)  = 1200.0
      xmin(nssmr)  = -1.0

      lnamsub(nset)   = 'Evapotranspiration'
      vnamsub(nset)   =  'ET'
      usub(nset)   =  'mm/day'
      xmax(nset)   = 150.0
      xmin(nset)   = -5.0

      lnamsub(nsrnfs) = 'Surface Runoff'
      vnamsub(nsrnfs) =  'RNFS'
      usub(nsrnfs) =  'mm/day'
      xmax(nsrnfs) = 2000.0
      xmin(nsrnfs) = -200.0

      lnamsub(nssnow) = 'Snow Depth'
      vnamsub(nssnow) =  'SNOW'
      usub(nssnow) =  'mm H2O'
      xmax(nssnow) = 1000.0
      xmin(nssnow) = -1.0

      lnamsub(nssh)   = 'Sensible Heat'
      vnamsub(nssh)   =  'SH'
      usub(nssh)   =  'W/m2'
      xmax(nssh)   = 1000.0
      xmin(nssh)   = -300.0

      lnamsub(nsprc)  = 'Convective Precip'
      vnamsub(nsprc)  =  'RC'
      usub(nsprc)  =  'mm/day'
      xmax(nsprc)  = 1500.0
      xmin(nsprc)  = -1.0

      lnamsub(nspt) = 'Total Precipitation'
      vnamsub(nspt) =  'RT'
      usub(nspt) =  'mm/day'
      xmax(nspt) = 2500.0
      xmin(nspt) = -1.0

      lnamsub(nspsrf) = 'Surface Pressure'
      vnamsub(nspsrf) =  'PSRF'
      usub(nspsrf) =  'hPa'
      xmax(nspsrf) = 1500.0
      xmin(nspsrf) = 300.0

      lnamsub(nsrha)  = 'Relative Humidity'
      vnamsub(nsrha)  =  'RHA'
      usub(nsrha)  =  'fraction'
      xmax(nsrha)  = 5.0
      xmin(nsrha)  = -0.1

      
      aaa = 2.**16.-1.
      do l=1,nsub2
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end

      SUBROUTINE MMVLUCHE(vnamche,lnamche,uche,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      character vnamche(nctot)*10, lnamche(nctot)*20, uche(nctot)*13
      real xmin(nctot), xmax(nctot), fact(nctot), offset(nctot), aaa
      integer l,n,r

      character tracname(ntrac)*10

      tracname(1) = 'TR1'
      tracname(2) = 'TR2'
      tracname(3) = 'TR3'
      tracname(4) = 'TR4'
      tracname(5) = 'TR5'
      tracname(6) = 'TR6'
      tracname(7) = 'TR7'
      tracname(8) = 'TR8'
      tracname(9) = 'TR9'
      tracname(10) = 'TR10'


       r = (nc3d-3)/ntrac
       do n=1,ntrac
       vnamche(r*(n-1)+1) = tracname(n)
       lnamche(r*(n-1)+1 ) = 'MMR_'//tracname(n)
       uche(r*(n-1)+1)    = 'micro-g/Kg'
       end do
       print*,'r is ',r
       lnamche(7) = 'aer mix. ext. coef'
       vnamche(7) = 'aext8'
       uche(7) =  'na'
       xmax(7) = 0.5
       xmin(7) = 0.

       lnamche(8) = 'aer mix. scat. alb'
       vnamche(8) = 'assa9'
       uche(8) =  'na'
       xmax(8) = 0.5
       xmin(8) = 0.

       lnamche(9) = 'aer mix. scat. alb'
       vnamche(9) = 'agfu8'
       uche(9) =  'na'
       xmax(9) = 0.5
       xmin(9) = 0.

       r= (nc2d-2)/ntrac

       do n=1,ntrac
       vnamche(nc3d+r*(n-1)+1)='BURD'//tracname(n)
       vnamche(nc3d+r*(n-1)+2)='WDLS'//tracname(n)
       vnamche(nc3d+r*(n-1)+3)='WDCV'//tracname(n)
       vnamche(nc3d+r*(n-1)+4)='DRDP'//tracname(n)
       vnamche(nc3d+r*(n-1)+5)='XGAZ'//tracname(n)
       vnamche(nc3d+r*(n-1)+6)='XSAQ'//tracname(n)
       vnamche(nc3d+r*(n-1)+7)='EMRA'//tracname(n)

       lnamche(nc3d+r*(n-1)+1)='Int column '//tracname(n)
       lnamche(nc3d+r*(n-1)+2)='Wetdep lsc '//tracname(n)
       lnamche(nc3d+r*(n-1)+3)='Wetdep cvc '//tracname(n)
       lnamche(nc3d+r*(n-1)+4)='Drydep surf'//tracname(n)
       lnamche(nc3d+r*(n-1)+5)='gaz conv'//tracname(n)
       lnamche(nc3d+r*(n-1)+6)='aq conv'//tracname(n)
       lnamche(nc3d+r*(n-1)+7)='Emission rate'//tracname(n)

       uche(nc3d+r*(n-1)+1)='mg/m2'
       uche(nc3d+r*(n-1)+2)='mg/m2'
       uche(nc3d+r*(n-1)+3)='mg/m2'
       uche(nc3d+r*(n-1)+4)='mg/m2'
       uche(nc3d+r*(n-1)+5)='mg/m2'
       uche(nc3d+r*(n-1)+6)='mg/m2'
       uche(nc3d+r*(n-1)+7)='micro-g/m2.s'


c       print*,'1 :',vnamche(1)
c       print*,'2 :',vnamche(2)
c       print*,'3 :',vnamche(3)
c       print*,'4 :',vnamche(4)
c       print*,'5 :',vnamche(5)
c       print*,'6 :',vnamche(6)
c       print*,'7 :',vnamche(7)
c       print*,'8 :',vnamche(8)

       end do

        do n=1,nc3d - 3
           xmax(n) = 10.E6
           xmin(n) = 0.
        end do

C           xmax(5) = 10.E9
C           xmax(6) = 10.E9
        do n=nc3d+1, nc3d+nc2d - 2
           xmax(n) = 10.E10
           xmin(n) = 0.
        end do

c      xmax(ncld)        = 1.1
c      xmax(nclwp)       = 5000.0
c      xmax(nqrs)        = 1.0e-3
c      xmax(nqrl)        = 1.0e-3
c      xmax(nr3d+nfsw)   = 1200.0
c     xmax(nr3d+nflw)   = 500.0
c      xmax(nr3d+nclrst) = 1500.0
c      xmax(nr3d+nclrss) = 1500.0
c      xmax(nr3d+nclrlt) = 1500.0
c      xmax(nr3d+nclrls) = 500.0
c      xmax(nr3d+nsolin) = 1500.0
c      xmax(nr3d+nsabtp) = 1500.0
c      xmax(nr3d+nfirtp) = 500.0

c      xmin(ncld)        = -0.1
c      xmin(nclwp)       = -10.0
c      xmin(nqrs)        = -1.0e-3
c      xmin(nqrl)        = -1.0e-3
c      xmin(nr3d+nfsw)   = -10.0
c      xmin(nr3d+nflw)   = -100.0
c      xmin(nr3d+nclrst) = -10.0
c      xmin(nr3d+nclrss) = -10.0
c      xmin(nr3d+nclrlt) = -10.0
c      xmin(nr3d+nclrls) = -10.0
c      xmin(nr3d+nsolin) = -10.0
c      xmin(nr3d+nsabtp) = -10.0
c      xmin(nr3d+nfirtp) = -10.0

       lnamche(52) = 'TOArad forcing'
       vnamche(52) = 'acsto'
       uche(52) =  'W/m^2'
       xmax(52) = 200.0
       xmin(52) = -200.0

       lnamche(53) = 'TOArad forcing'
       vnamche(53) = 'acsto'
       uche(53) =  'W/m^2'
       xmax(53) = 200.0
       xmin(53) = -200.0

           aaa = 2.**16.-1.
       do l=1,nctot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end

      SUBROUTINE MMVLURAD(vnamrad,lnamrad,urad,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'
      character vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
      real xmin(nrtot), xmax(nrtot), fact(nrtot), offset(nrtot), aaa
      integer l

      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp

      vnamrad(ncld)        = 'FC'
      vnamrad(nclwp)       = 'CLWP'
      vnamrad(nqrs)        = 'QRS'
      vnamrad(nqrl)        = 'QRL'
      vnamrad(nr3d+nfsw)   = 'FSW'
      vnamrad(nr3d+nflw)   = 'FLW'
      vnamrad(nr3d+nclrst) = 'CLRST'
      vnamrad(nr3d+nclrss) = 'CLRSS'
      vnamrad(nr3d+nclrlt) = 'CLRLT'
      vnamrad(nr3d+nclrls) = 'CLRLS'
      vnamrad(nr3d+nsolin) = 'SOLIN'
      vnamrad(nr3d+nsabtp) = 'SABTP'
      vnamrad(nr3d+nfirtp) = 'FIRTP'

      lnamrad(ncld)        = 'Cloud Fraction'
      lnamrad(nclwp)       = 'Cld Liquid H2O Path'
      lnamrad(nqrs)        = 'Solar Heating Rate'
      lnamrad(nqrl)        = 'LW Cooling Rate'
      lnamrad(nr3d+nfsw)   = 'Surface Abs solar'
      lnamrad(nr3d+nflw)   = 'LW Cooling of Surf'
      lnamrad(nr3d+nclrst) = 'Clr Sky Col Abs Sol'
      lnamrad(nr3d+nclrss) = 'Clr Sky Surf Abs Sol'
      lnamrad(nr3d+nclrlt) = 'Clr Sky Net Up Flx'
      lnamrad(nr3d+nclrls) = 'Clr Sky LW Surf Cool'
      lnamrad(nr3d+nsolin) = 'Instant Incid Solar'
      lnamrad(nr3d+nsabtp) = 'Column Abs Solar'
      lnamrad(nr3d+nfirtp) = 'Net Up Flux at Top'

      urad(ncld)        = 'fraction'
      urad(nclwp)       = 'g/m2'
      urad(nqrs)        = 'K/s'
      urad(nqrl)        = 'K/s'
      urad(nr3d+nfsw)   = 'W/m2'
      urad(nr3d+nflw)   = 'W/m2'
      urad(nr3d+nclrst) = 'W/m2'
      urad(nr3d+nclrss) = 'W/m2'
      urad(nr3d+nclrlt) = 'W/m2'
      urad(nr3d+nclrls) = 'W/m2'
      urad(nr3d+nsolin) = 'W/m2'
      urad(nr3d+nsabtp) = 'W/m2'
      urad(nr3d+nfirtp) = 'W/m2'
 
      xmax(ncld)        = 1.1
      xmax(nclwp)       = 5000.0
      xmax(nqrs)        = 1.0e-2
      xmax(nqrl)        = 1.0e-2
      xmax(nr3d+nfsw)   = 1200.0
      xmax(nr3d+nflw)   = 500.0
      xmax(nr3d+nclrst) = 1500.0
      xmax(nr3d+nclrss) = 1500.0
      xmax(nr3d+nclrlt) = 1500.0
      xmax(nr3d+nclrls) = 500.0
      xmax(nr3d+nsolin) = 1500.0
      xmax(nr3d+nsabtp) = 1500.0
      xmax(nr3d+nfirtp) = 500.0

      xmin(ncld)        = -0.1
      xmin(nclwp)       = -10.0
      xmin(nqrs)        = -1.0e-2
      xmin(nqrl)        = -1.0e-2
      xmin(nr3d+nfsw)   = -10.0
      xmin(nr3d+nflw)   = -100.0
      xmin(nr3d+nclrst) = -10.0
      xmin(nr3d+nclrss) = -10.0
      xmin(nr3d+nclrlt) = -10.0
      xmin(nr3d+nclrls) = -10.0
      xmin(nr3d+nsolin) = -10.0
      xmin(nr3d+nsabtp) = -10.0
      xmin(nr3d+nfirtp) = -10.0

      aaa = 2.**16.-1.
      do l=1,nrtot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end


c      SUBROUTINE TMINMAX(f2d,f2davg,nx,ny,nfld,nhr,ihr
c     &         , nta,ntmax,ntmin)
c     
c      implicit none
c      integer nx, ny, nfld, nta, ntmax, ntmin, nhr, i, j, ihr
c      real f2d(nx,ny,nfld), f2davg(nx,ny,nfld,nhr)

c      do j=1,ny
c      do i=1,nx
c        if (ihr.le.1) then
c          f2d(i,j,ntmax) = f2d(i,j,nta)
c          f2d(i,j,ntmin) = f2d(i,j,nta)
c        else
c          if (f2d(i,j,nta) .gt. f2d(i,j,ntmax)) then
c            f2d(i,j,ntmax) = f2d(i,j,nta)
c          elseif (f2d(i,j,nta) .lt. f2d(i,j,ntmin)) then
c            f2d(i,j,ntmin) = f2d(i,j,nta)
c          end if
c        end if
c      end do
c      end do

c      return
c      end

      SUBROUTINE CALCMSE(f2d,f3d,nx,ny,nz,nfld2d,nfld3d
     &       , zs,sigh,pt,nta,nqva,npsa,nmse,nx1,ny1)

      implicit none

      integer nx,ny,nz,nfld2d,nfld3d
      real f2d(nx,ny,nfld2d), f3d(nx,ny,nz,nfld3d)
      integer nta,nqva,npsa,nmse

      real tv(nx,ny,nz), p(nx,ny,nz), z(nx,ny,nz)
      real zs(nx,ny), sigh(nz)
      real g, r, cp, ps, lh, ta, qa, za, t, q, pt, dp, tvbar
      integer i, j, k, k1, nx1, ny1

      parameter (g=9.81, r=287., cp=1004. ,lh=2.5e06)

      do j=1,ny1
      do i=1,nx1
        t = f3d(i,j,nz,nta)
        q = f3d(i,j,nz,nqva)
        ps = f2d(i,j,npsa)
        tv(i,j,nz) = t*(1.+0.608*q)
        p(i,j,nz) = (ps-pt*10.)*sigh(nz) + pt*10.
        dp = log(p(i,j,nz)) - log(ps)
        z(i,j,nz) = zs(i,j) - dp*tv(i,j,nz)*r/g
        f3d(i,j,nz,nmse) = g*z(i,j,nz) + cp*t + lh*q
      end do
      end do
      do k=nz-1,1,-1
      do j=1,ny1
      do i=1,nx1
        k1 = k+1
        t = f3d(i,j,k,nta)
        q = f3d(i,j,k,nqva)
        ps = f2d(i,j,npsa)
        tv(i,j,k) = t*(1.+0.608*q)
        tvbar = 0.5*(tv(i,j,k)+tv(i,j,k1))
        p(i,j,k) = (ps-pt*10.)*sigh(k) + pt*10.
        dp = log(p(i,j,k)) - log(p(i,j,k1))
        z(i,j,k) = z(i,j,k1) - dp*tvbar*r/g
        f3d(i,j,k,nmse) = g*z(i,j,k) + cp*t + lh*q
      end do
      end do
      end do

      return
      end

      SUBROUTINE CALCMSE2D(f2d,zs,nx,ny,nfld,nta,nqa,nmse)

      implicit none
      integer nx, ny, nfld, nta, nqa, nmse, nhr, ihr
      real zs(nx,ny)
      real f2d(nx,ny,nfld)

      real g, cp, lh, ta, qa, za
      integer i, j, k

      parameter (g=9.81, cp=1004. ,lh=2.5e06)

      do j=1,ny
      do i=1,nx
        ta = f2d(i,j,nta)
        qa = f2d(i,j,nqa)
        za = zs(i,j) + 2.
        f2d(i,j,nmse) = g*za + cp*ta + lh*qa
      end do
      end do

      return
      end

      

      SUBROUTINE CALCRH2D(f2d,nx,ny,nfld,nta,nqa,npsrf,nrh,vmisdat)

      implicit none
      integer nx, ny, nfld, nta, nqa, npsrf, nrh, i, j
      real f2d(nx,ny,nfld), vmisdat
      real pres, ta, qa, satvp, qs, svp1, svp2, svp3, ep2
      parameter (svp1=6.112, svp2=17.67, svp3=29.65, ep2=0.622)

      do j=1,ny
      do i=1,nx
        pres = f2d(i,j,npsrf)
        if (pres.gt.0.) then
          ta = f2d(i,j,nta)
          qa = f2d(i,j,nqa)
          if (ta.gt.273.15) then
            satvp = svp1*exp(svp2*(ta-273.15)/(ta-svp3)) ! SAT'N VAP PRES
          else
            satvp=svp1*exp(22.514-6.15e3/ta)
          end if
          qs = ep2*satvp/(pres-satvp)                ! SAT. MIXING RATIO
          f2d(i,j,nrh) = qa/qs
        else
          f2d(i,j,nrh) = vmisdat
        end if
      end do
      end do

      return
      end

      SUBROUTINE CALCHGT(fld2d,fld3d,nx,ny,nz,nfld2d,nfld3d
     &       , zs,sigf,sigh,pt,nta,nqva,npsa,nhgt,nx1,ny1)

      implicit none

      integer nx,ny,nz,nfld2d,nfld3d
      real fld2d(nx,ny,nfld2d), fld3d(nx,ny,nz,nfld3d)
      integer nta,nqva,npsa,nhgt

      integer i,j,k,k1,k2,nx1,ny1
      real zs(nx,ny),sigh(nz),dsig(nz),sigf(nz+1)
      real pt,ep1,rg,g,rgog,ps,pf,tv,t,t1,t2,q,q1,q2,tv1,tv2

      ep1 = 0.608
      rg=287.04
      g = 9.805
      rgog = rg/g

      do k=1,nz
        dsig(k)=sigf(k)-sigf(k+1)
      end do

      do i=1,nx1
      do j=1,ny1
        ps = fld2d(i,j,npsa)/10.
        t = fld3d(i,j,1,nta)
        q = fld3d(i,j,1,nqva)
        pf=pt/(ps-pt)
        tv = t*(1.0+ep1*q)
        fld3d(i,j,nz,nhgt)=zs(i,j)+tv*rgog*log((1.+pf)/(sigh(nz)+pf))
        DO k1=nz-1,1,-1
          k2=k1+1
          t1 = fld3d(i,j,k1,nta)
          t2 = fld3d(i,j,k2,nta)
          q1 = fld3d(i,j,k1,nqva)
          q2 = fld3d(i,j,k2,nqva)
          tv1 = t1*(1.0+ep1*q1)
          tv2 = t2*(1.0+ep1*q2)
          tv=(tv1*dsig(k1)+tv2*dsig(k2))/(dsig(k1)+dsig(k2))
          fld3d(i,j,k1,nhgt)=fld3d(i,j,k2,nhgt)
     &              +tv*rgog*log((sigh(k2)+pf)/(sigh(k1)+pf))
        end do
      end do
      end do
      
      return
      end

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

      SUBROUTINE WRITEBATHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls
     &         , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhro
     &         , iotyp)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      real f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny), xlon(nx,ny)
     &   , dlat(nx,ny), dlon(nx,ny), zs(nx,ny), ls(nx,ny)
      real tmp2d(nx,ny)
      integer ils(nx,ny), itex(nx,ny), icol(nx,ny), iotyp
      real rootf(20), vegc(20), seasf(20), rough(20), displa(20)
     &   , rsmin(20), xla(20), xlai0(20), sai(20), sqrtdi(20)
     &   , fc(20), depuv(20), deprv(20), deptv(20), iexsol(20)
     &   , kolsol(20), albvgs(20), albvgl(20)
      real xmopor(12), xmosuc(12), xmohyd(12), xmowil(12), bee(12)
     &   , skrat(12)
      character varnam*10, lname*20, units*13
      real vmin, vmax, vmisdat, misdat, fact, offset, xmin, xmax, aaa
      integer ndim, idim(ndim), idout, i, j, julnc
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real*8 xhro

c* root f: fraction of roots in root zone
      data rootf/.30,.80,.67,.67,.50,.80,.80,.90,.90,.30,.80,
     &  7*.50,2*0.5/
c* vegc is maximum fractional cover of vegetation
      data vegc/.85,.8,.8,.8,.8,.9,.8,0.0,0.6,0.8,0.35,0.0,.8,2*0.,
     &3*0.8,0.8,0./
c* seasf is the difference between vegc and fractional cover at 269k
      data seasf/0.6,0.1,0.1,0.3,0.3,0.5,0.3,0.,0.2,0.6,0.1,0.0,0.4,2*0.
     &0,.2,.3,.2,0.4,0.0/
c* rough is an aerodynamic roughness length (m) =approx 0.1*veg height
c* also used snow masking depth in subrout albedo
      data rough/0.08,0.05,2*1.0,0.8,2.0,0.1,0.05,0.04,0.06,0.1,0.01,
     &0.03,2*0.0004,2*0.1,0.8,0.3,0.0/
c ******      displacement height (meter)
c ******      if great parts of veg. are covered by snow, use displa=0
c ******      because then the new displa-theory is not valid
      data displa/0.,0.,9.,9.,0.,18.,14*0./
c ******      min stomatl resistance (s/m)
      data rsmin/45.,60.,2*80.,120.,2*60.,200.,80.,45.,150.,200.,45.
     &          ,2*200.,80.,120.,100.,2*120./
c ******      max leaf area index (ratio unit cover per unit ground)
      data xla/6.,2.,5*6.,0.,3*6.,0.,6.,2*0.,3*6.,6.0,0.0/
c ******      min leaf area index **lai depends on temp as veg cover
      data xlai0/0.5,0.5,5.,2*1.,5.,0.5,0.,3*0.5,0.,0.5,2*0.,5.,1.,3.,
     & 0.5, 0.0/
c ******      stem area index (projected area of non-transpiring sfcs)
      data sai/ .5, 4., 5*2., 2*.5, 9*2., 2*2. /
c ******      inverse square root of leaf dimension - used for
c ******      calculating fluxes from foliage
      data sqrtdi/10.,17*5.0,2*5./
c ******      fc = light dependence of stomatal resistance
      data fc/ .02, .02, 4*.06, 11*.02, .06, 2*.02 /

c ******      depuv is depth of upper soil layer (mm)
c ******      deprv is depth of root zone (mm)
c ******      deptv is depth of total soil (mm)
      data depuv/20*100./
      data deprv/2*1000.,2*1500.,2000.,1500.,11*1000.,2000.,2*2000./
      data deptv/18*3000.,2*3000./
 
c ******      iexsol is soil texture type (see subr soilbc)
      data iexsol/6,6,6,6,7,8,6,1,6,6,3,12,6,6,6,6,5,6,6,0/
c ******      kolsol is soil color type (see subr. albedo)
      data kolsol/5,3,4,4,4,4,4,1,3,3,2,1,5,5,5,4,3,4,4,0/
c ******      xmopor is fraction of soil that is voids
      data xmopor/.33,.36,.39,.42,.45,.48,.51,.54,.57,.6,.63,.66/
c ******      xmosuc is the minimum soil suction (mm)
      data xmosuc/3*30.0,9*200./
c ******      xmohyd is the max. hydraulic conductivity (mm/s)
      data xmohyd/0.20e-0,0.80e-1,0.32e-1,0.13e-1,0.89e-2,0.63e-
     &2,0.45e-2,0.32e-2,0.22e-2,0.16e-2,0.11e-2,0.80e-3/
c ******      xmowilt is fraction of water content at which permanent
c                 wilting occurs
      data xmowil/.095,.128,.161,.266,.3,.332,.378,.419,.455,
     &.487,.516,.542/
c ******      bee is the clapp and hornbereger "b" parameter
      data bee/3.5,4.0,4.5,5.0,5.5,6.0,6.8,7.6,8.4,9.2,10.0,10.8/
c ******      bskrat is ratio of soil thermal conduc. to that of loam -
c                 a function of texture
      data skrat/1.7,1.5,1.3,1.2,1.1,1.0,.95,.90,.85,.80,.75,.7/

c ******      albvgs is vegetation albedo for wavelengths < 0.7 microns
      data albvgs/.1,.1,.05,.05,.08,.04,.08,.2,.1,.08,.17,.8,.06,2*.07,
     &.05,.08,.06,2*0.06/
c ******      albvgl is vegetation albedo for wavelengths > 0.7 microns
      data albvgl/.3,.3,.23,.23,.28,.20,.30,.4,.3,.28,.34,.6,.18,2*.2,
     &.23,.28,.24,2*.18/

      idim(3) = 1
      aaa = 2.**16.-1.
c     CALL XTDOT(zs,zs,nx,ny,1,nx-1,ny-1)
c     CALL XTDOT(f,f,nx,ny,1,nx-1,ny-1)
      do j=1,ny
      do i=1,nx
        ils(i,j) = anint(ls(i,j))
        itex(i,j) = anint(iexsol(ils(i,j)))
        icol(i,j) = anint(kolsol(ils(i,j)))
      end do
      end do


      varnam='ZB'
      print*,varnam
      lname='Terrain Elevation'
      units='m'
      xmax=7000.
      xmin=-100.
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = zs(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(zb)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(zb)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &   , 1.0,0,misdat,iotyp)

      varnam='LU'
      print*,varnam
      lname='Land Use Type'
      units='unitless'
      xmax=21.
      xmin=-1.
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = ls(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(lu)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(lu)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='F'
      print*,varnam
      lname='Coriolus'
      units='rad/sec'
      xmax=0.001
      xmin=-0.001
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = f(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(f)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(f)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMAP'
      print*,varnam
      lname='Cross-Grid Map Fact'
      units='unitless'
      xmax=2.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = xmap(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmap)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmap)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DMAP'
      print*,varnam
      lname='Dot Grid Map Factor'
      units='degrees'
      xmax=2.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = dmap(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dmap)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(dmap)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLAT'
      print*,varnam
      lname='Cross Grid Latitude'
      units='degrees'
      xmax=100.0
      xmin=-100.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = xlat(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xlat)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xlat)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLON'
      print*,varnam
      lname='Cross Grid Longitude'
      units='degrees'
      xmax=200.0
      xmin=-200.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = xlon(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xlon)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xlon)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DLAT'
      print*,varnam
      lname='Dot Grid Latitude'
      units='degrees'
      xmax=100.0
      xmin=-100.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = dlat(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dlat)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(dlat)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DLON'
      print*,varnam
      lname='Dot Grid Longitude'
      units='degrees'
      xmax=200.0
      xmin=-200.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=3,ny-2
      do i=3,nx-2
        tmp2d(i-2,j-2) = dlon(i,j)
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dlon)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(dlon)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='land'
      print*,varnam
      lname='Land Mask'
      units='unitless'
      xmax=1.5
      xmin=-0.5
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        if (ils(i,j).eq.15) then
          tmp2d(i,j) = -999.
        else
          tmp2d(i,j) = 1.0
        end if
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(land)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(land)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='ROOTF'
      print*,varnam
      lname='Root Zone Root Frac'
      units='fraction'
      xmax=1.5
      xmin=-0.5
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = rootf(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(rootf)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(rootf)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='VEGC'
      print*,varnam
      lname='Max Vegetation Cover'
      units='fraction'
      xmax=1.5
      xmin=-0.5
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = vegc(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(vegc)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(vegc)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='SEASF'
      print*,varnam
      lname='VEGC(MAX)-VEGC @269K'
      units='fraction'
      xmax=1.5
      xmin=-0.5
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = seasf(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(seasf)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(seasf)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='ROUGH'
      print*,varnam
      lname='Roughness Length'
      units='m'
      xmax=3.5
      xmin=-0.5
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = rough(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(rough)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(rough)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DISPLA'
      print*,varnam
      lname='Displacement Height'
      units='m'
      xmax=30.0
      xmin=-2.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = displa(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(displa)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(displa)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='RSMIN'
      print*,varnam
      lname='Min Stomatl Resist'
      units='s/m'
      xmax=300.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = rsmin(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(rsmin)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(rsmin)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLA'
      print*,varnam
      lname='Max Leaf Area Index'
      units='unitless'
      xmax=10.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xla(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xla)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xla)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLAI0'
      print*,varnam
      lname='Min Leaf Area Index'
      units='unitless'
      xmax=10.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xlai0(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xlai0)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xlai0)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='SAI'
      print*,varnam
      lname='Stem Area Index'
      units='unitless'
      xmax=10.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = sai(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(sai)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(sai)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='SQRTDI'
      print*,varnam
      lname='Inv SQRT Leaf Dim'
      units='m-0.5'
      xmax=50.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = sqrtdi(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(sqrtdi)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(sqrtdi)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='FC'
      print*,varnam
      lname='Light Depend on RS'
      units='unitless'
      xmax=50.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = fc(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(fc)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(fc)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DEPUV'
      print*,varnam
      lname='Top Soil Layer Depth'
      units='mm'
      xmax=200.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = depuv(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(depuv)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(depuv)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DEPRV'
      print*,varnam
      lname='Root Soil Depth'
      units='mm'
      xmax=4000.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = deprv(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(deprv)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(deprv)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DEPTV'
      print*,varnam
      lname='Total Soil Depth'
      units='mm'
      xmax=8000.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = deptv(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(deptv)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(deptv)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='IEXSOL'
      print*,varnam
      lname='Soil Texture'
      units='unitless'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = iexsol(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(iexsol)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(iexsol)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='KOLSOL'
      print*,varnam
      lname='Soil Color'
      units='unitless'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = kolsol(ils(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(kolsol)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(kolsol)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMOPOR'
      print*,varnam
      lname='Soil Porosity'
      units='fraction'
      xmax=1.5
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xmopor(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmopor)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmopor)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMOSUC'
      print*,varnam
      lname='Min Soil Suction'
      units='mm'
      xmax=300.0
      xmin=-10.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xmosuc(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmosuc)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmosuc)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMOHYD'
      print*,varnam
      lname='Sat Soil Conductiv'
      units='mm/s'
      xmax=1.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xmohyd(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmohyd)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmohyd)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMOWIL'
      print*,varnam
      lname='Soil Wilting Point'
      units='Fraction'
      xmax=1.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = xmowil(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmowil)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmowil)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='BEE'
      print*,varnam
      lname='Clapp Hornbereger'
      units='unitless'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = bee(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(bee)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(bee)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='SKRAT'
      print*,varnam
      lname='Thermal Conductivity'
      units='ratio'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = skrat(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(skrat)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(skrat)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='ALBVGS'
      print*,varnam
      lname='Veg Albedo < 0.7 um'
      units='fraction'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = albvgs(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(albvgs)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(albvgs)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='ALBVGL'
      print*,varnam
      lname='Veg Albedo > 0.7 um'
      units='fraction'
      xmax=20.0
      xmin=-1.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      do j=1,ny
      do i=1,nx
        tmp2d(i,j) = albvgl(itex(i,j))
      end do
      end do
      if (iotyp.eq.1) then
        CALL GETMINMAX(tmp2d,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(albvgl)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(albvgl)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,tmp2d,nx,ny,1,idim,xhro,lname
     &   , units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      return 
      end

      SUBROUTINE WRITEHEAD(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls
     &         , vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim,idout,xhro
     &         , iotyp)

      implicit none
      include 'postproc.param'
      include 'postproc1.param'

      character varnam*10, lname*20, units*13
      real f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny), xlon(nx,ny)
     &   , dlat(nx,ny), dlon(nx,ny), zs(nx,ny), zssd(nx,ny), ls(nx,ny)
      integer ndim, idout, iotyp
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real aaa, xmin, xmax, fact, offset, vmisdat, misdat, vmin, vmax
      integer idim(ndim)
      real*8 xhro

      aaa = 2.**16.-1.
c     CALL XTDOT(zs,zs,nx,ny,1,nx-1,ny-1)
c     CALL XTDOT(f,f,nx,ny,1,nx-1,ny-1)

      idim(3) = 1

      varnam='HT'
      lname='Terrain Elevation'
      units='m'
      xmax=7000.
      xmin=-100.
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(zs,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(zs)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(zs)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,zs,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='HTSD'
      lname='Elevation Std Dev'
      units='m'
      xmax=5000.
      xmin=-100.
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(zssd,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(zssd)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(zssd)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,zssd,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='LU'
      lname='Land Use Type'
      units='unitless'
      xmax=21.
      xmin=-1.
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(ls,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(ls)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(ls)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,ls,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='F'
      lname='Coriolus'
      units='rad/sec'
      xmax=0.001
      xmin=-0.001
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(f,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(f)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(f)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,f,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XMAP'
      lname='Cross-Grid Map Fact'
      units='unitless'
      xmax=2.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(xmap,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xmap)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xmap)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,xmap,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DMAP'
      lname='Dot Grid Map Factor'
      units='degrees'
      xmax=2.0
      xmin=0.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(dmap,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dmap)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(dmap)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,dmap,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLAT'
      lname='Cross Grid Latitude'
      units='degrees'
      xmax=100.0
      xmin=-100.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(xlat,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xlat)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xlat)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,xlat,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='XLON'
      lname='Cross Grid Longitude'
      units='degrees'
      xmax=200.0
      xmin=-200.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(xlon,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(xlon)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(xlon)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,xlon,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DLAT'
      lname='Dot Grid Latitude'
      units='degrees'
      xmax=100.0
      xmin=-100.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(dlat,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dlat)=',vmin,'XMIN=',xmin
        print*,'MAXVAL(dlat)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,dlat,nx,ny,1,idim,xhro,lname,units
     &   , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      varnam='DLON'
      lname='Dot Grid Longitude'
      units='degrees'
      xmax=200.0
      xmin=-200.0
      fact=(xmax-xmin)/aaa
      offset=(xmax+xmin)/2.
      if (iotyp.eq.1) then
        CALL GETMINMAX(dlon,nx,ny,1,vmin,vmax,vmisdat)
        if (vmin.lt.xmin .or. vmax.gt.xmax) then
          print*,'Values Out of Range:  FIELD=',varnam
          print*,'MINVAL(dlon)=',vmin,'XMIN=',xmin
          print*,'MAXVAL(dlon)=',vmax,'XMAX=',xmax
          stop 999
        end if
        misdat = xmin
      elseif (iotyp.eq.2) then
        misdat = vmisdat
      end if
      CALL WRITECDF(idout,varnam,dlon,nx,ny,1,idim,xhro,lname,units
     &  , fact,offset,vvarmin,vvarmax,xlat1d,xlon1d
     &    ,1.0,0,misdat,iotyp)

      return
      end

      SUBROUTINE SETCONST(f,const,n1,n2,n3,n4,n5,i1,ni1,i2,ni2)

      implicit none
      integer i,j,k,l,m,n1,n2,n3,n4,n5,i1,ni1,i2,ni2
      real f(n1,n2,n3,n4,n5), const

      do m=1,n5
      do l=1,n4
      do k=1,n3
      do j=i2,ni2
      do i=i1,ni1
        f(i,j,k,l,m) = const
      end do
      end do
      end do
      end do
      end do

      return
      end

      SUBROUTINE GETMINMAX(f,n1,n2,n3,vmin,vmax,vmisdat)

      implicit none
      integer i,j,k,n1,n2,n3
      real f(n1,n2,n3),vmin,vmax,vmisdat,misdat

      if (vmisdat.gt.0.0) then
        misdat = -1.0*vmisdat
      else
        misdat = vmisdat
      end if
      vmin=1.e30
      vmax=-1.e30
      do k=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,k).gt.misdat) then
          if (f(i,j,k).gt.vmax) vmax=f(i,j,k)
          if (f(i,j,k).lt.vmin) vmin=f(i,j,k)
        end if
      end do
      end do
      end do

      return
      end


      SUBROUTINE PARAM(nx,ny,kz,np,ds,clat,clon,xplat,xplon
     &         , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim
     &         , plv)

      implicit none
      integer nx, ny, kz, ndim,np
      real ds, clat, clon, xplat, xplon, xlat(nx,ny), xlon(nx,ny)
      real vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      integer idim(ndim)
      integer i, j
      logical plv

      vvarmin(1) = xlon(1,ny/2)
      vvarmin(2) = xlat(nx/2,1)
      vvarmin(3) = 1050.
      vvarmax(1) = xlon(nx,ny/2)
      vvarmax(2) = xlat(nx/2,ny)
      vvarmax(3) = 0.
      idim(1) = nx
      idim(2) = ny
      if (.not.plv) then
      idim(3) = kz
      else
      idim(3) = np
      endif
      do i=1,nx
        xlon1d(i) = xlon(i,ny/2)
      end do
      do j=1,ny
        xlat1d(j) = xlat(nx/2,j)
      end do

      return
      end


      SUBROUTINE WRITEGRADS(idout,outvar,nx,ny,nk,nrec)

      implicit none
      integer nx, ny, nk, idout, nrec, i, j, k
      real outvar(nx,ny,nk), xlat1d(ny), xlon1d(nx)

      do k=1,nk
        nrec = nrec + 1
        write(idout,rec=nrec) ((outvar(i,j,k),i=1,nx),j=1,ny)
      end do

      return
      end

      SUBROUTINE HTSIG(f3d,f2d,nta,nhga,npsa,zs,sig,pt,im,jm,km,im1,jm1
     & ,n3d,n2d)
      implicit none
      include 'postproc.param'
	include 'postproc1.param'
      INTEGER IM,JM,KM,n2d,n3d
      REAL    f3d(IM,JM,KM,n3d),f2d(IM,JM,n2d)
      REAL    PSTAR(IM,JM),ZS(IM,JM)
      REAL    SIG(KM)
      real    rgas,grav,bltop,tlapse
      data    rgas,grav,bltop,tlapse
     &       /287.04,9.80616,0.96,-6.5E-3/
      INTEGER I,J,K,nta,nhga,npsa,im1,jm1
      REAL    TBAR,pt
	  pstar=f2d(:,:,npsa)
      DO J=1,JM1
      DO I=1,IM1
         f3d(I,J,KM,nhga) = ZS(I,J) + RGAS/GRAV*f3d(I,J,KM,nta)
     &             * LOG(PSTAR(I,J)/((PSTAR(I,J)-PT)*SIG(KM)+PT))
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM1
      DO I=1,IM1
         TBAR = 0.5*( f3d(I,J,K,nta)+f3d(I,J,K+1,nta) )
         f3d(I,J,K,nhga) = f3d(I,J,K+1,nhga) +RGAS/GRAV*TBAR
     &            * LOG(((PSTAR(I,J)-PT)*SIG(K+1)+PT)
     &                 /((PSTAR(I,J)-PT)*SIG(K)+PT))
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
	subroutine intlin(fp,f,f2,npsa,pt,sig,im,jm,km,nv,p
     &           ,      kp,n3d,n2d,im1,jm1)
      implicit none
      integer im,jm,km,im1,jm1,kp,n3d,n2d,nv,npsa
      real    fp(im,jm,kp,n3d),f(im,jm,km,n3d)
      real    pstar(im,jm),f2(im,jm,n2d)
      real    sig(km),p(kp)
      real    rgas,grav,bltop,tlapse
      data    rgas,grav,bltop,tlapse
     &       /287.04,9.80616,0.96,-6.5E-3/
      integer i,j,k,n,pt,skip
      integer k1,k1p
      real    sigp,wp,w1
	pstar = f2(:,:,npsa)
      do j=1,jm1
      do i=1,im1
      do n=1,kp
         sigp = (p(n)-pt) / (pstar(i,j)-pt)
         k1=0
         do k=1,km
            if (sigp.gt.sig(k)) k1=k
         enddo
            if(sigp.le.sig(1)) then
            fp(i,j,n,nv) = f(i,j,1,nv)
         	else if((sigp.gt.sig(1)).and.(sigp.lt.sig(km))) then
            k1p = k1 + 1
            wp  = (sigp-sig(k1))/(sig(k1p)-sig(k1))
            w1  = 1.-wp
            fp(i,j,n,nv)  = w1*f(i,j,k1,nv)+wp*f(i,j,k1p,nv)
         	else if(sigp.ge.sig(km)) then
            fp(i,j,n,nv)  = f(i,j,km,nv)
         endif
      enddo
      enddo
      enddo
      return
      end

	subroutine intlog(fp,f,f2,npsa,pt,sig,im,jm,km,nv,p
     &           ,      kp,n3d,n2d,im1,jm1)
      implicit none
      integer im,jm,im1,jm1,km,kp,n3d,n2d,nv
      real    fp(im,jm,kp,n3d),f(im,jm,km,n3d)
      real    pstar(im,jm),f2(im,jm,n2d)
      real    sig(km),p(kp),pt
      real    rgas,grav,bltop,tlapse
      data    rgas,grav,bltop,tlapse
     &       /287.04,9.80616,0.96,-6.5E-3/
      integer i,j,k,n,npsa
      integer k1,k1p,kbc
      real    sigp,wp,w1
	
	pstar=f2(:,:,npsa)
	
      do k=1,km
        if(sig(k).lt.bltop) kbc = k
      enddo
      do j=1,jm1
      do i=1,im1
      do n=1,kp
         sigp = (p(n)-pt) / (pstar(i,j)-pt)
         k1=0
         do k=1,km
            if (sigp.gt.sig(k)) k1=k
         enddo
            if(sigp.le.sig(1)) then
              fp(i,j,n,nv) = f(i,j,1,nv)
            else if((sigp.gt.sig(1)).and.(sigp.lt.sig(km))) then
              k1p = k1 + 1
              wp  = log(sigp/sig(k1)) / log(sig(k1p)/sig(k1))
              w1  = 1. - wp
              fp(i,j,n,nv)= w1*f(i,j,k1,nv) + wp*f(i,j,k1p,nv)
            else if((sigp.ge.sig(km)).and.(sigp.le.1.))then
              fp(i,j,n,nv)= f(i,j,km,nv)
            else if(sigp.gt.1.) then
              fp(i,j,n,nv) = f(i,j,kbc,nv) 
     &        * exp(-rgas*tlapse*log(sigp/sig(kbc))/grav)
            endif
      enddo
      enddo
      enddo
	
	return
      end
	
	subroutine height(fp,f,f2,nta,npsa,zs,sig,im,jm,km,kp,nhga
     &           ,      p,n3d,n2d,pt,im1,jm1)

      implicit none
      integer im,jm,im1,jm1,km,kp,n3d,n2d,nhga,npsa,nta
      real    t(im,jm,km),f(im,jm,km,n3d),fp(im,jm,kp,n3d)
      real    pstar(im,jm),zs(im,jm),f2(im,jm,n2d)
      real    sig(km),p(kp)
	real    rgas,grav,bltop,tlapse
      data    rgas,grav,bltop,tlapse
     &       /287.04,9.80616,0.96,-6.5E-3/
      real    psig(100),pt
      integer i,j,k,kbc,n,kt,kb
      real    psfc,temp,wt,wb

      do k=1,km
         if(sig(k).lt.bltop) kbc=k
      enddo
	pstar=f2(:,:,npsa)
      do j=1,jm1
      do i=1,im1
      do k=1,km
         psig(k) = sig(k) * (pstar(i,j)-pt) + pt
      enddo
         psfc = pstar(i,j)
         do n = 1,kp
           kt = 1
           do k=1,km
             if (psig(k).lt.p(n)) kt=k
           enddo
           kb = kt + 1
           if(p(n).le.psig(1)) then
             temp = f(i,j,1,nta)
             fp(i,j,n,nhga) =f(i,j,1,nhga)
     &        +rgas*temp*log(psig(1)/p(n))/grav
           else if((p(n).gt.psig(1)).and.(p(n).lt.psig(km))) then
             wt = log(psig(kb)/p(n)) / log(psig(kb)/psig(kt))
             wb = log(p(n)/psig(kt)) / log(psig(kb)/psig(kt))
             temp = wt * f(i,j,kt,nta) + wb * f(i,j,kb,nta)
             temp = ( temp + f(i,j,kb,nta) ) / 2.
             fp(i,j,n,nhga) =f(i,j,kb,nhga)+rgas*temp*log(psig(kb)/p(n))/grav
           else if((p(n).ge.psig(km)).and.(p(n).le.psfc)) then
             temp = f(i,j,km,nta)
             fp(i,j,n,nhga) =zs(i,j)+rgas*temp*log(psfc/p(n))/grav
           else if(p(n).gt.psfc) then
             temp = f(i,j,kbc,nta) - tlapse * (f(i,j,kbc,nta)-zs(i,j))
             fp(i,j,n,nhga) =zs(i,j)-(temp/tlapse)
     &             * ( 1.-exp(-rgas*tlapse*log(p(n)/psfc)/grav))
           endif
      enddo
      enddo
      enddo
      return
      end
      subroutine calcslp(f3,f2,nhga,nta,npsa,pt,zs,slp,sig
     &                  ,im,jm,km,n3d,n2d,nx1,ny1)
      implicit none
      integer im,jm,km
      real    f3(im,jm,km,n3d),f2(im,jm,n2d)
      real    pstar(im,jm),zs(im,jm)
      real    sig(km),pt
      real    rgas,grav,bltop,tlapse
      parameter (rgas=287.04,grav=9.80616,bltop=0.96,tlapse=6.5E-3)
      integer kbc,i,j,k,nhga,nta,n3d,n2d,npsa,slp,nx1,ny1
      real    tsfc
	pstar(:,:)=f2(:,:,npsa)
      do k=1,km
         if(sig(k).lt.bltop) kbc=k
      enddo
      do j=1,ny1
      do i=1,nx1
         tsfc = f3(i,j,kbc,nta)-tlapse*(f3(i,j,kbc,nhga)-zs(i,j))
		 f2(i,j,slp) = pstar(i,j)
     &   * exp( -grav/(rgas*tlapse)*log(1.-zs(i,j)*tlapse/tsfc))
      enddo
      enddo

      return
      end
      SUBROUTINE CALCVD(fld3d,nx,ny,nz,nfld3d,ds,dmap,xmap
     &       , nua,nva, nvor,ndiv,nx1,ny1)

      implicit none

      integer nx,ny,nz,nfld3d
      real fld3d(nx,ny,nz,nfld3d)
      integer nua, nva, nvor, ndiv

      integer   i, j, k, nx1, ny1
      real	    ds,ds2r
      real      u(nx,ny,nz),v(nx,ny,nz),u1,u2,u3,u4
      real      v1,v2,v3,v4
      real      dmap(nx,ny),xmap(nx,ny)
      ds2r=1.0/(2.0*ds)
      
	u(:,:,:) = fld3d(:,:,:,nua)
	v(:,:,:) = fld3d(:,:,:,nva)
      
      do k=1,nz
      do j=1,ny1
      do i=1,nx1
        u1=u(i  ,j  ,k)/dmap(i  ,j  )
        u2=u(i+1,j  ,k)/dmap(i+1,j  )
        u3=u(i  ,j+1,k)/dmap(i  ,j+1)
        u4=u(i+1,j+1,k)/dmap(i+1,j+1)
        v1=v(i  ,j  ,k)/dmap(i  ,j  )
        v2=v(i+1,j  ,k)/dmap(i+1,j  )
        v3=v(i  ,j+1,k)/dmap(i  ,j+1)
        v4=v(i+1,j+1,k)/dmap(i+1,j+1)
        fld3d(i,j,k,nvor)=xmap(i,j)*xmap(i,j)*ds2r*
     &       ((v4-v2+v3-v1)-(u2-u1+u4-u3))
        fld3d(i,j,k,ndiv)=xmap(i,j)*xmap(i,j)*ds2r*
     &       ((u3-u1+u4-u2)+(v2-v1+v4-v3))
        
      enddo
      enddo
      enddo
      return
      end
      SUBROUTINE CALCRH(fld2d,fld3d,nx,ny,nz,nfld2d,nfld3d
     &       , sigh,pt,nta,nqva,npsa,nrh,ntd,nth,nx1,ny1)

      implicit none

      integer nx,ny,nz,nfld2d,nfld3d
      real fld2d(nx,ny,nfld2d), fld3d(nx,ny,nz,nfld3d)
      real sigh(nz)
      integer nta, nqva, npsa, nrh, ntd, nth


      integer i, j, k, nx1, ny1
      real pres, pt, t, q, hl, satvp, qs,x,dpd,tmp

      real svp1, svp2, svp3, ep2, rovcp
      parameter (svp1=6.112, svp2=17.67, svp3=29.65, ep2=0.622)
	rovcp = 287.04/1004.
      do k=1,nz
      do j=1,ny1
      do i=1,nx1
        pres = (fld2d(i,j,npsa)-pt*10.)*sigh(k) + pt*10. ! PRES AT LEVEL K
        t = fld3d(i,j,k,nta)
        q = fld3d(i,j,k,nqva)
        if (t.gt.273.15) then
          satvp = svp1*exp(svp2*(t-273.15)/(t-svp3)) ! SATURATION VAP PRESS.
        else
          satvp=svp1*exp(22.514-6.15e3/t)
        end if
        qs = ep2*satvp/(pres-satvp)                   ! SAT. MIXING RATIO
        fld3d(i,j,k,nrh) = q/qs
	  X = 1.-fld3d(i,j,k,nrh)
	  TMP = t-273.15
	  DPD = (14.55+0.144*TMP)*X
     &       + 2*((2.5+0.007*TMP)*X)**3
     &       + (15.9+0.117*TMP)*(X**14)  
	  fld3d(i,j,k,ntd)= t-DPD	! DEW POINT TEMP
	  fld3d(i,j,k,nth)=t*(1000./pres)**rovcp	! POTENTIAL TEMP

      end do
      end do
      end do
      return
      end

