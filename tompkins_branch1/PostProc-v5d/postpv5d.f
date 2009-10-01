compile on linux: pgf77 -byteswapio -fast -Ktrap=fp -o postpv5d postpv5d.f -L../Commons/env/liblinux  -lnetcdf
compile on SUN: f77 -fast -fns -ftrap=common -o postpv5d postpv5d.f -L../Commons/env/libsun -lnetcdf
      implicit none
      include 'postpv5d.param'
      logical there
      character filout*70, filbat*70, filrad*70, filchem*70
     &        , inout(5)*70, inbat(5)*70, inrad(5)*70, inchem(5)*70
     &        , inhead*70
     &        , vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
     &        , vnamrad(nrtot)*10, lnamrad(nrtot)*20, urad(nrtot)*13
     &        , vnamche(nctot)*10,lnamche(nctot)*20
     &        , uche(nctot)*13 
     &        , vnamv5d(nv5d)*10,lnamv5d(nv5d)*20
     &        , uv5d(nv5d)*13 

      integer  nvarlev(nv5d)

      real xminout(notot), xmaxout(notot), factout(notot)
     &   , offsetout(notot)
     &   , xminrad(nrtot), xmaxrad(nrtot), factrad(nrtot)
     &   , offsetrad(nrtot)
     &   , xminche(nctot), xmaxche(nctot),factche(nctot)
     &   ,offsetche(nctot)


      real vmisdat,vmisdatv5d,clat,clon,ds,pt,xplat,xplon,
     &     ddeg,xday,xdayold
     &   , xtime,xtimeold
     
c   modif v5D

      parameter (vmisdatv5d=1.E30)

      integer iin,ndim,idout,idday,ierr,i,j,ii,l,nb
     &      , iyr,iyr0,iyr1,iyr2,iyrx,imo,imo0,imo1,imo2,imox
     &      , idy,idy0,idy1,idy2,idyx,ihr,ihr0,ihr1,ihr2,ihrx
     &      , mdate,idate,idate0,idate1,idate2,idatex,idateold,idatenew
     &      ,julnc,julnc0,julnc1,julnc2,julncx,julmid
     &      ,julda,julda0,julda1,julda2
     &      , numtimes, numvars,nzv5d,itv5d,itcv5d,itrv5d 
     &      , indfv5d,iio,iic,iir
      parameter(iin=10,ndim=3)

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      
      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)

      
c pressure grid
      logical zpress       
      common /press/zpress

      real ofld3d_p,  o3davg_p
      common /outfld3d_p/ofld3d_p(nx,ny,nzp,nout3d), 
     &           o3davg_p(nx,ny,nzp,nout3d,nhrout)

      real  cfld3d_p,  c3davg_p  
      common /cheflds_p/ cfld3d_p(nx,ny,nzp,nc3d),
     &             c3davg_p(nx,ny,nzp,nc3d,nhrout)

      real  rfld3d_p,  r3davg_p
      common /radflds_p/  rfld3d_p(nx,ny,nzp,nr3d),
     &             r3davg_p(nx,ny,nzp,nr3d,nhrout)


      real f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny)
     &    , xlon(nx,ny), dlat(nx,ny), dlon(nx,ny), zs(nx,ny)
     &    ,zssd(nx,ny), ls(nx,ny)
      real sigh(nz), sighrev(nz), sigf(nz+1), sigb(2)
     &   , vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)

      real*8 xhr, xhr0, xhr1, xhr2, xhrm, xhrdy
      integer idatev5d, idim(ndim), nday, ntime
     &      , nouttime(nhrout), nradtime(nhrrad),nchetime(nhrche)

      integer   nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt, nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf, ntopo,nlus
      common /opoint3d/  nua,  nva,  nta, nqva, nqca,nmse,nrh,nhgt,nsig
      data          nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt, nsig
     &            /    1,    2,    3,    4,    5,    6,    7,  8, 9 /
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf,ntopo,nlus
      data              npsa,  nrt, ntgb, nsmt, nbf, ntopo, nlus
     &            /       1,    2,    3,    4,    5 , 6, 7/

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

      integer naer1, naer2
      common /chepoint3d/ naer1, naer2
      data               naer1, naer2
     &                  /1,      2   / 


      logical out,head,rad,chem,outhead,noavg,avg,dayavg,diuavg
      common /logic/out,rad,chem 



       print*,'read  postp.in file'

      open (11,file=files,status='old')
      read (11,*) idate0
      read (11,*) idate1
      read (11,*) idate2
      read (11,*) zpress
      read (11,*) head
      read (11,*) inhead
      read (11,*) inout(1)
      read (11,*) inout(2)
      read (11,*) inout(3)
      read (11,*) inout(4)
      read (11,*) inout(5)
      read (11,*) inrad(1)
      read (11,*) inrad(2)
      read (11,*) inrad(3)
      read (11,*) inrad(4)
      read (11,*) inrad(5)
      read (11,*) inchem(1)
      read (11,*) inchem(2)
      read (11,*) inchem(3)
      read (11,*) inchem(4)
      read (11,*) inchem(5)
      read (11,*) outhead
      read (11,*)
      read (11,*) out
      read (11,*) chem
      read (11,*) rad
      read (11,*) nday
      read (11,*) noavg
      read (11,*) avg
      read (11,*) dayavg
      read (11,*) diuavg
      read (11,*) filout
      read (11,*)

 
      print*,' '
      print*,'**************************************************'
      print*,'IDATE0 =',idate0,'IDATE1 =',idate1,'IDATE2 =',idate2
      print*,'**************************************************'
      print*,' '
      print*,' '

      if (zpress .and. out ) then
      print*,' '
      print*,'Sigma level interpolated on pressure grid defined'
      print*,'in postpv5d.param file !! '

      else if ( zpress .and. .not.out) then 
      print*,' '
      print*,' need to read atmospheric fields to perform any 
     &         vertical interpolation '
      stop

      end if

      if ( ( noavg .and. (avg .or. dayavg .or. diuavg))
     & .or.( avg .and. (noavg .or. dayavg .or. diuavg)) 
     & .or.( dayavg .and. (avg .or. noavg .or. diuavg)) 
     & .or.( diuavg .and. (avg .or. noavg .or. dayavg))
     &   ) then 
       print*, 'noavg =',noavg, 'avg', avg
       print*, 'dayavg',dayavg,'diuavg',diuavg 
       print*,' cannot be performed at the same time !' 
       stop       
      else if ( .not. noavg .and. .not. avg .and. .not. dayavg 
     & .and. .not. diuavg) then 
       print*, ' choose an averaging option !'
        stop 
    
       end if
     
     

      if (out) then
        print*,' '
        print*,' '
        print*,' '
        print*,'**************************************************'
        print*,'***INPUT***  ATMOS OUTPUT FROM MODEL:  INOUT = ',inout
        print*,'**************************************************'
        print*,'  WRITING ALL OUTPUT DATA'
        print*,'**OUTPUT**  FILOUT = ',filout
     
      end if

  
        if (chem ) then
        print*,' '
        print*,' '
        print*,' '
        print*,'**************************************************'
        print*,'***INPUT***  CHEM OUTPUT FROM MODEL: INCHEM = ',inchem
        print*,'**************************************************'
        print*,'  WRITING ALL CHEM DATA'
        print*,'**OUTPUT**  FILOUT = ',filout
      end if

      if (rad) then
        print*,' '
        print*,' '
        print*,' '
        print*,'**************************************************'
        print*,'***INPUT***  RAD OUTPUT FROM MODEL:  INRAD = ',inrad
        print*,'**************************************************'
        print*,'  WRITING ALL RAD DATA'
        print*,'**OUTPUT**  FILOUT = ',filout       
      end if
      print*,' '
      print*,' '
      print*,'**************************************************'
      print*,' '
      print*,' '


      print*,idate0,idate1,idate2
      CALL JULIAN(idate0,julnc0,julda0,iyr0,imo0,idy0,ihr0)
      xhr0 = float(julnc0)
      CALL JULIAN(idate1,julnc1,julda1,iyr1,imo1,idy1,ihr1)
      xhr1 = float(julnc1)
      CALL JULIAN(idate2,julnc2,julda2,iyr2,imo2,idy2,ihr2)
      xhr2 = float(julnc2)
      print*,'xhr0=',xhr0,'xhr1=',xhr1,'xhr2=',xhr2
      julmid = (julnc1+julnc2)/2
      xhrm = float(julmid-mod(julmid,24))

  

     
C **** READ HEADER FILE **** C
      open (iin,file=inhead,status='old',form='unformatted'
     &    ,recl=nx*ny*4,access='direct')
      CALL RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon
     &   , f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,mdate,iin)
      print*,'HEADER READ IN',clon,clat
c
      close (iin)
c
      
      CALL PARAM(nx,ny,nz,ds,ddeg,clat,clon,xplat,xplon
     &     , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim)


c**** open the data files*******       
        iio = 1
        open (iin,file=inout(iio),status='old',form='unformatted')
        print*,'INPUT (ATM) FILE: ',inout(iio)

        iic = 1
        if (chem) then
        open (iin+1,file=inchem(iic),status='old',form='unformatted')
        print*,'INPUT (CHEM) FILE: ',inchem(iic)
        end if

        iir = 1        
        if (rad) then
        open (iin+2,file=inrad(iir),status='old',form='unformatted')
        print*,'INPUT (RAD) FILE: ',inrad(iir)
        end if


C **** SETUP MIN, MAX, VARNAM, LNAME, UNITS DATA FOR NetCDF ****
        CALL MMVLUOUT(vnamout,lnamout,uout,xminout,xmaxout
     &     , factout,offsetout)


        CALL MMVLUCHE(vnamche,lnamche,uche,xminche,xmaxche
     &     , factche,offsetche)
c

        CALL MMVLURAD (vnamrad,lnamrad,urad,xminrad,xmaxrad
     &    , factrad,offsetrad)
        
       
c***** PREPARE THE V5D OUTPUT FILE **********

       if (noavg) then 
         dtv5d = dtout
         numtimes= int ( (xhr2-xhr1)/dtv5d )+ 1   
         filout = 'ALL_'//filout 
       end if 

       if (avg) then      
         dtv5d =1
         numtimes= 1
         filout = 'AVG_'//filout 
       end if     

       if (dayavg) then      
         dtv5d = dtout 
         numtimes= int ((xhr2-xhr1)/(dtv5d * nhrout) ) 
         filout = 'DAY_'//filout 
       end if    

       if (diuavg) then 
         dtv5d =dtout
         numtimes = nhrout
         filout = 'DIU_'//filout 
       end if 

    
       if (zpress)then
         nzv5d= nzp
       else
         nzv5d= nz
       end if
       
       vnamv5d(1:notot)=vnamout(:)
       lnamv5d(1:notot)=lnamout(:)
       uv5d(1:notot)=uout(:)
       nvarlev(1:nout3d)=nzv5d
       nvarlev(nout3d+1:notot)=1

 
       if ( out .and. chem ) then

       vnamv5d(notot+1:notot+nctot)=vnamche(:) 
       lnamv5d(notot+1:notot+nctot)=lnamche(:) 
       uv5d(notot+1:notot+nctot)= uche(:) 
       nvarlev(notot+1:notot+nc3d)=nzv5d
       nvarlev(notot+nc3d+1:notot+nctot)=1 
       end if

       if ( out .and. rad .and. .not.chem) then

       vnamv5d(notot+1:notot+nrtot)=vnamrad(:) 
       lnamv5d(notot+1:notot+nrtot)=lnamrad(:) 
       uv5d(notot+1:notot+nrtot)= urad(:) 
       nvarlev(notot+1:notot+nr3d)=nzv5d
       nvarlev(notot+nr3d+1:notot+nrtot)=1 
       end if
       print*,vnamv5d

       if ( out .and. chem .and. rad) then
       vnamv5d(notot+nctot +1:notot+nctot+nrtot)=
     &                       vnamrad(:)      

       lnamv5d(notot+nctot +1:notot+nctot+nrtot)=
     &                       lnamrad(:)         

       uv5d(notot+nctot +1:notot+nctot+nrtot)=
     &                       urad(:)  

       nvarlev(notot+nctot+1:notot+nctot+nr3d)=nzv5d
       nvarlev(notot+nctot+nr3d+1:notot+nctot+nrtot)=1

       end if
      print*,nvarlev
       idatev5d = idate1/1000000
       itv5d = 1
       itcv5d= 1 
       itrv5d= 1 

       if (zpress)then

           filout = 'P_'//filout
           CALL FEXIST(filout)
           print*,'OUTPUT V5D FILE: ',filout

           CALL PREP_V5D (ihr1,ihr2,julda1,julda2,idatev5d,dtv5d,
     &                numtimes,
     &                nv5d,
     &                nvarlev,
     &                vnamv5d,lnamv5d,uv5d,clat,clon,xlat,xlon,
     &                ds,nzp, plev*10.,filout)

       else 

           filout = 'S_'//filout
           CALL FEXIST(filout)
           print*,'OUTPUT V5D FILE: ',filout         

           CALL PREP_V5D (ihr1,ihr2,julda1,julda2,idatev5d,
     &                dtv5d,numtimes,
     &                nv5d,
     &                nvarlev,
     &                vnamv5d,lnamv5d,uv5d,clat,clon,xlat,xlon,
     &                ds,nz,sigh,filout)
      
       end if


C **** ZERO OUT AVERAGE ARRAYS **** C
c        if (outavg .or. outdiur .or. outday) then
         if (avg .or. dayavg .or. diuavg ) then
          CALL SETCONST(o3davg,vmisdatv5d,nx,ny,nz,nout2d,
     & nhrout,1,nx,1,ny)
          CALL SETCONST(o2davg,vmisdatv5d,nx,ny,nout2d,
     & nhrout,1,1,nx,1,ny)
          CALL SETCONST(o3davg,0.0,nx,ny,nz,nout2d,nhrout,1,nx1,1,ny1)
          CALL SETCONST(o2davg,0.0,nx,ny,nout2d,nhrout,1,1,nx1,1,ny1)
          do l=1,nhrout
            nouttime(l) = 0
          end do
        
          if (chem) then
           CALL SETCONST(c3davg,vmisdatv5d,nx,ny,nz,nc3d,
     & nhrche,1,nx,1,ny)
           CALL SETCONST(c2davg,vmisdatv5d,nx,ny,nc2d,
     & nhrche,1,1,nx,1,ny)
           CALL SETCONST(c3davg,0.0,nx,ny,nz,nc3d,nhrche,1,nx1,1,ny1)
           CALL SETCONST(c2davg,0.0,nx,ny,nc2d,nhrche,1,1,nx1,1,ny1)
           do l=1,nhrche
            nchetime(l) = 0
           end do
          end if

          if (rad) then
           CALL SETCONST(r3davg,vmisdatv5d,nx,ny,nz,nr3d,
     &  nhrrad,1,nx,1,ny)
           CALL SETCONST(r2davg,vmisdatv5d,nx,ny,nr2d,
     &  nhrrad,1,1,nx,1,ny)
           CALL SETCONST(r3davg,0.0,nx,ny,nz,nr3d,nhrrad,1,nx1,1,ny1)
           CALL SETCONST(r2davg,0.0,nx,ny,nr2d,nhrrad,1,1,nx1,1,ny1)
           do l=1,nhrche
            nradtime(l) = 0
           end do
          end if

        end if



c**********************************      
c**********************************
C **** PROCESS THE DATA ***********  
c**********************************
c**********************************

c************************************
c****** evolving fields **********
c*******************************

        idate = idate0
        idateold=idate       
        do while(idate.le.idate2)


        indfv5d=0

        


c************************
C ***** 1 ATM DATA ******
c************************

          idatenew = idateold + nint(dtout)
          CALL JULIAN(idatenew,julnc,julda,iyr,imo,idy,ihr)
          CALL RDOUT(vnamout,lnamout,uout,idate,iin,ierr)
c test        
     

C **** for end of arm  files **** C
          if (ierr.ne.0 .or. idate.gt.idatenew) then
           print*,'END OF ATM FIL REACHED:File Index =',iio,'ierr=',ierr
            ierr = 0
            iio = iio + 1
            print*,inout(iio)
            CALL FEXISTNEW(inout(iio),there)
            if (.not.there) go to 99
            close (iin)
            open (iin,file=inout(iio),status='old',form='unformatted')
            idate = 0
            do while(idate.lt.idatenew)
              print*,'SEARCHING FOR PROPER DAY:',idate,idatenew
              CALL RDOUT(vnamout,lnamout,uout,idate,iin,ierr)
              if (ierr.ne.0) stop 'READ OUT ERROR'
              if (idate.gt.idatenew) then
                print*,'FILE ERROR (ATM): DATE EXCEEDED'
                print*,idate,idatenew
                stop 999
              end if
            end do
          end if
          
          CALL JULIAN(idate,julnc,julda,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtout)
          if (ihr.eq.0) ihr=24/nint(dtout)
     
c *** more  surface and diagnostic variabls         
          ofld2d(:,:,ntopo)=zs(:,:)        
          where( ls(:,:) == 14 .or. ls(:,:)==15 )
            ofld2d(:,:,nlus)= 0
          elsewhere 
            ofld2d(:,:,nlus)= 1
          end where

          CALL CALCMSE(zs,sigh,pt)
          CALL CALCRH(sigh,pt)
          CALL CALCHGT(zs,sigf,sigh,pt)

c *** vertical interpolation ***

          if (zpress) CALL VERT_INT('out',zs,sighrev,vmisdatv5d)

c **** average atm data  **** c
          
          if ( avg .or.dayavg .or. diuavg ) then 
            if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING ATM DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(o2davg,ofld2d,nx,ny,nout2d,nhrout,ihr
     &         , vmisdatv5d)

          if ( .not. zpress) then
            CALL AVGDATA3D(o3davg,ofld3d,nx,ny,nz,nout3d,nhrout,ihr
     &         , vmisdatv5d)
          else

            CALL AVGDATA3D(o3davg_p,ofld3d_p,nx,ny,nzp,nout3d,nhrout,ihr
     &         , vmisdatv5d)

          end if 
             nouttime(ihr) = nouttime(ihr) + 1

               if(dayavg) then           
                ntime = 0
                do l=1,nhrout
                 ntime = ntime + nouttime(l)
                end do
                print*,'ntime ATM',ntime
                CALL JULIAN(idateold,julncx,julda,iyrx,imox,idyx,ihrx)
                if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtout).and.nday.gt.0)) then
                 idatex = iyr*1000000 + imo*10000 + 1500
                 CALL JULIAN(idatex,julncx,julda,iyrx,imox,idyx,ihrx)
                 xhrdy = float(julncx)
                 CALL WRITEAVGOUT_V5D(nouttime,vmisdatv5d,itv5d,indfv5d)
                 print*,'DAILY DATA WRITTEN: ',xhr,idate
                 itv5d = itv5d+1
                 print*,'write avg day',itv5d
                 do l=1,nhrout
                  nouttime(l) = 0
                 end do                 
                end if  
               end if 
            end if
          end if

C **** WRITE atm data at each time step **** c
       
          if (out .and. noavg) then        
          if (idate.ge.idate1.and.idate.le.idate2) then
                     
            CALL WRITEOUT_V5D(vmisdatv5d,itv5d,indfv5d)
            itv5d=itv5d + 1
          end if
          end if


c**********************************
C**** 2 CHEM DATA *****************
c**********************************
          if (chem) then

            idatenew = idateold + nint(dtche)
            CALL JULIAN(idatenew,julnc,julda,iyr,imo,idy,ihr)
            if (mod(xhr-xhr0,dtche) .eq. 0) then
            CALL RDCHE(iin+1,idate,ierr)
            end if 

C **** for end of che  files **** C
          if (ierr.ne.0 .or. idate.gt.idatenew) then
           print*,'END OF CHE FIL REACHED:File Index =',iic,'ierr=',ierr
            ierr = 0
            iic = iic + 1
            print*,inchem(iic)
            CALL FEXISTNEW(inchem(iic),there)
            if (.not.there) go to 99
            close (iin +1)
           open (iin+1,file=inchem(iic),status='old',form='unformatted')
            idate = 0
            do while(idate.lt.idatenew)
              print*,'SEARCHING FOR PROPER DAY:',idate,idatenew
              CALL RDCHE(iin+1,idate,ierr)
              if (ierr.ne.0) stop 'READ OUT ERROR'
              if (idate.gt.idatenew) then
                print*,'FILE ERROR (CHE): DATE EXCEEDED'
                print*,idate,idatenew
                stop 999
              end if
            end do
          end if

          CALL JULIAN(idate,julnc,julda,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtche)
          if (ihr.eq.0) ihr=24/nint(dtche)


c *** vertical interpolation ***

          if (zpress) CALL VERT_INT('che',zs,sighrev,vmisdatv5d)

c **** average che data  **** c
          
          if ( avg .or.dayavg .or. diuavg ) then 
            if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING CHE DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(c2davg,cfld2d,nx,ny,nc2d,nhrche,ihr
     &         , vmisdatv5d)

          if (.not. zpress) then
            CALL AVGDATA3D(c3davg,cfld3d,nx,ny,nz,nc3d,nhrche,ihr
     &         , vmisdatv5d)
          else
            CALL AVGDATA3D(c3davg_p,cfld3d_p,nx,ny,nzp,nc3d,nhrche,ihr
     &         , vmisdatv5d)
          end if

             nchetime(ihr) = nchetime(ihr) + 1

               if(dayavg) then           
                ntime = 0
                do l=1,nhrche
                 ntime = ntime + nchetime(l)
                end do
                print*,'ntime CHE',ntime
                CALL JULIAN(idateold,julncx,julda,iyrx,imox,idyx,ihrx)
                if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtche).and.nday.gt.0)) then
                 idatex = iyr*1000000 + imo*10000 + 1500
                 CALL JULIAN(idatex,julncx,julda,iyrx,imox,idyx,ihrx)
                 xhrdy = float(julncx)
                 CALL WRITEAVGCHE_V5D(nchetime,vmisdatv5d,
     &                 itcv5d,indfv5d)
                 print*,'DAILY DATA WRITTEN: ',xhr,idate
                 itcv5d = itcv5d+1
                 print*,'write avg day',itv5d
                 do l=1,nhrche
                  nchetime(l) = 0
                 end do                 
                end if  
               end if 
            end if
          end if


C **** WRITE chem data at each time step **** c

           if (chem .and. noavg) then
           if (idate.ge.idate1.and.idate.le.idate2) then            
   
            CALL WRITECHE_V5D(vmisdatv5d,itcv5d,indfv5d)
             itcv5d = itcv5d+1
           end if
           end if

         end if ! fin chem

c**********************************
C**** 3 RAD  DATA *****************
c**********************************
          if (rad) then
 

            idatenew = idateold + nint(dtche)
            CALL JULIAN(idatenew,julnc,julda,iyr,imo,idy,ihr)
            if (mod(xhr-xhr0,dtrad) .eq. 0) then
            CALL RDRAD(iin+2,idate,ierr)
            end if 
C **** for end of che  files **** C
          if (ierr.ne.0 .or. idate.gt.idatenew) then
           print*,'END OF RADFIL REACHED:File Index =',iir,'ierr=',ierr
            ierr = 0
            iir = iir + 1
            print*,inrad(iir)
            CALL FEXISTNEW(inrad(iir),there)
            if (.not.there) go to 99
            close (iin +2)
           open (iin+2,file=inrad(iir),status='old',form='unformatted')
            idate = 0
            do while(idate.lt.idatenew)
              print*,'SEARCHING FOR PROPER DAY:',idate,idatenew
              CALL RDRAD(iin+2,idate,ierr)
              if (ierr.ne.0) stop 'READ OUT ERROR'
              if (idate.gt.idatenew) then
                print*,'FILE ERROR (CHE): DATE EXCEEDED'
                print*,idate,idatenew
                stop 999
              end if
            end do
          end if
  
          CALL JULIAN(idate,julnc,julda,iyr,imo,idy,ihr)
          xhr = float(julnc)
          ihr = ihr/nint(dtrad)
          if (ihr.eq.0) ihr=24/nint(dtrad)


c *** vertical interpolation ***

            if (zpress) CALL VERT_INT('rad',zs,sighrev,vmisdatv5d)


c **** average rad data  **** c
   
          if ( avg .or.dayavg .or. diuavg ) then 
            if (idate.gt.idate1 .and. idate.le.idate2) then
            print*,'AVERAGING RAD DATA: ',idate,xhr,ihr
            CALL AVGDATA2D(r2davg,rfld2d,nx,ny,nr2d,nhrrad,ihr
     &         , vmisdatv5d)

          if (.not. zpress) then 
            CALL AVGDATA3D(r3davg,rfld3d,nx,ny,nz,nr3d,nhrrad,ihr
     &         , vmisdatv5d)
          else
            CALL AVGDATA3D(r3davg_p,rfld3d_p,nx,ny,nzp,nr3d,nhrrad,ihr
     &        ,vmisdatv5d)
          end if

             nradtime(ihr) = nradtime(ihr) + 1
         
               if(dayavg) then           
                ntime = 0
                do l=1,nhrrad
                 ntime = ntime + nradtime(l)
                end do
                print*,'ntime RAD',ntime
                CALL JULIAN(idateold,julncx,julda,iyrx,imox,idyx,ihrx)
                if ((nday.eq.-1.and.imo.ne.imox) .or.
     &            (ntime.ge.nint(nday*24./dtrad).and.nday.gt.0)) then
                 idatex = iyr*1000000 + imo*10000 + 1500
                 CALL JULIAN(idatex,julncx,julda,iyrx,imox,idyx,ihrx)
                 xhrdy = float(julncx)
                 CALL WRITEAVGRAD_V5D(nradtime,vmisdatv5d,
     &                                itrv5d,indfv5d)
                 print*,'DAILY DATA WRITTEN: ',xhr,idate
                 itrv5d = itrv5d+1
                 print*,'write avg day',itrv5d
                 do l=1,nhrrad
                  nradtime(l) = 0
                 end do                 
                end if  
               end if 
            end if
          end if

C **** WRITE rad data at each time step **** c

           if (rad .and. noavg) then
           if (idate.ge.idate1.and.idate.le.idate2) then            
 
            print*,itcv5d        
            CALL WRITERAD_V5D(vmisdatv5d,itrv5d,indfv5d)
            itrv5d = itrv5d+1
           end if
           end if

         end if ! fin rad



C **** INCREMENT TIME **** 
          idateold = idate
          idate = idate + nint(dtout)
          CALL JULIAN(idate,julnc,julda,iyr,imo,idy,ihr)
        end do
C**** 
       

 
c **** CLOSE INPUT (OUT) FILE **** c
 99   close (iin)
      close (iin+1)          
      close (iin+2)   


C **** WRITE OUT AVERAGED ATM  FIELDS **** C
c      if (outavg .or. outdiur) then
       if (out) then

       if (avg .or. diuavg) then
        do ihr=1,nhrout
          print*,nhrout
          if (nouttime(ihr).le.0) then
            print*,'Not enough data for average.'
            print*,'nouttime must have values great than zero.'
            print*,nouttime
            stop 'OUTAVG-OUTDIUR'
          end if
        end do
      end if

      if (avg) then
       itv5d=1 
       CALL WRITEAVGOUT_V5D(nouttime,vmisdatv5d,itv5d,indfv5d)
               
        print*,'DONE WRITING AVERAGED DATA IN V5D!!!'
      end if

      if (diuavg) then 
     
      CALL WRITEDIUROUT_V5D(xhrm,nouttime,vmisdatv5d,indfv5d)
      end if  

      end if

C **** WRITE OUT AVERAGED CHE  FIELDS **** C

      if (chem) then

       if (avg .or. diuavg) then
        do ihr=1,nhrche
          print*,nhrche
          if (nchetime(ihr).le.0) then
            print*,'Not enough data for average.'
            print*,'nouttime must have values great than zero.'
            print*,nouttime
            stop 'OUTAVG-OUTDIUR'
          end if
        end do
      end if

      if (avg) then
       itcv5d=1 
       CALL WRITEAVGCHE_V5D (nchetime,vmisdatv5d,itcv5d,indfv5d )                    
        print*,'DONE WRITING AVERAGED DATA IN V5D!!!'
      end if


      if (diuavg) then 
     
      CALL WRITEDIURCHE_V5D(xhr1,nchetime,vmisdatv5d,indfv5d)
      end if  

      end if


C **** WRITE OUT AVERAGED RAD FIELDS **** C

      if (rad) then

       if (avg .or. diuavg) then
        do ihr=1,nhrrad
          print*,nhrrad
          if (nradtime(ihr).le.0) then
            print*,'Not enough data for average.'
            print*,'nouttime must have values great than zero.'
            print*,nradtime
            stop 'OUTAVG-OUTDIUR'
          end if
        end do
      end if

      if (avg) then
       itrv5d=1 
       CALL WRITEAVGRAD_V5D (nchetime, vmisdatv5d, itrv5d,indfv5d )                    
        print*,'DONE WRITING AVERAGED DATA IN V5D!!!'
      end if


      if (diuavg) then 
     
      CALL WRITEDIURRAD_V5D(xhr1,nchetime,vmisdatv5d,indfv5d)
      end if  

      end if



c **** CLOSE V5D file **** c        




      CALL CLSV5D


      end



c######################################################################
c####################################################################
c#####################################################################

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
      SUBROUTINE RDHEAD(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon
     &         , f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,mdate,iin)

      implicit none
      include 'postpv5d.param'

      real   sigf(nz+1), f(nx,ny), xmap(nx,ny), dmap(nx,ny), xlat(nx,ny)
     &     , xlon(nx,ny), dlat(nx,ny), dlon(nx,ny), zs(nx,ny)
     &     , zssd(nx,ny), ls(nx,ny), sigh(nz), sighrev(nz)
     &     , ds, clat, clon, pt, xplat, xplon, dto, dtb, dtr
      integer iin,mdate,mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,k,kk
     &        ,iotyp
      character proj*6

      rewind (iin)
      read(iin,rec=1) mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,sigf,ds
     &       ,pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,iotyp
      print*,'mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds='
      print*,mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds
      print*,'sigf='
      print*,sigf
      print*,'pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr='
      print*,pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr
      if (ni.ne.nyf .or. nj.ne.nxf .or. nz.ne.nk) then
        print*, 'Grid Dimensions DO NOT MATCH'
        print*, 'nyf=',nyf,'nxf=',nxf,'nz=',nz
        print*, 'ni=',ni,'nj=',nj,'nk=',nk
        stop 'BAD DIMENSIONS'
      end if
      if (dto.ne.dtout .or. dtb.ne.dtbat .or. dtr.ne.dtrad) then
        print*,'OUTPUT INTERVALS ARE IMPROPERLY DEFINED'
        print*,'dto=',dto,'dtout=',dtout
        print*,'dtb=',dtb,'dtbat=',dtbat
        print*,'dtr=',dtr,'dtrad=',dtrad
        stop 'BAD TIME PARAMETERS'
      end if
c     print*,'ZS'
      read (iin,rec=2) zs
c     print*,'ZSSD'
      read (iin,rec=3) zssd
c     print*,'LU'
      read (iin,rec=4) ls
c     print*,'LU'
      read (iin,rec=5) ls
c     print*,'XLAT'
      read (iin,rec=6) xlat
c     print*,'XLON'
      read (iin,rec=7) xlon
c     print*,'XMAP'
      read (iin,rec=8) xmap
c     print*,'DMAP'
      read (iin,rec=9) dmap
c     print*,'F'
      read (iin,rec=10) f

      do k=1,nz
        kk = nz - k + 1
        sigh(k) = (sigf(k)+sigf(k+1))/2.
        sighrev(kk) = sigh(k)
      end do

      return
      end



      SUBROUTINE RDOUT(vnamout,lnamout,uout,mdate,iin,ierr)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)

      real tmp2d(nx,ny), xday
      integer mdate, iin, ierr, no, i,j,k
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13

      integer nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt, nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf
      common /opoint3d/  nua,  nva,  nta, nqva, nqca, nmse,nrh,nhgt,nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf

   
      read (iin,iostat=ierr) mdate
      if (ierr.ne.0) return
      xday = xday/1440.
      print *,'Reading output:  ',mdate
      do no=1,no3d
c       print*,'READ VAR:  ',vnamout(no),lnamout(no)
        do k=1,nz
          read (iin,iostat=ierr) tmp2d
          if (ierr.ne.0) return
          do j=1,ny
          do i=1,nx
            ofld3d(i,j,k,no) = tmp2d(i,j)
          end do
          end do
        end do
      end do
      
      do no=1,no2d
        read (iin,iostat=ierr) tmp2d
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          ofld2d(i,j,no) = tmp2d(i,j)
        end do
        end do
      end do
c   
      return
      end

      SUBROUTINE WRITEOUT_V5D( vmisdatv5d ,itv5d, indfv5d)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)

      logical zpress       
      common /press/zpress

      real ofld3d_p,  o3davg_p
      common /outfld3d_p/ofld3d_p(nx,ny,nzp,nout3d), 
     &           o3davg_p(nx,ny,nzp,nout3d,nhrout)

      integer i, j, k, no, nno, itv5d,indfv5d
      real vmisdatv5d

      real tmp2d(nx,ny), tmp3d(nx,ny,nz)
      real tmp3d_p(nx,ny,nzp)
    
     
      nno = indfv5d
 
      if (.not.zpress)then

      do no=1,nout3d
c       print*,no,vnamout(no)
        nno = nno + 1
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = ofld3d(i,j,k,no)
        end do
        end do
        end do
               
         CALL WRITEV5D(itv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
      end do
c    if pressure grid for v5d
      else
       do no=1,nout3d
          nno = nno + 1
         do k=1,nzp
         do j=1,ny1
         do i=1,nx1
          tmp3d_p(i,j,k) = ofld3d_p(i,j,k,no)
         end do
         end do
         end do        
         CALL WRITEV5D(itv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
       end do

      end if

      indfv5d= indfv5d + nout3d 
      
      nno = indfv5d

      do no=1,nout2d
        nno = nno + 1
c       print*,no,nno,vnamout(nno)
        do j=1,ny1
        do i=1,nx1
          tmp2d(i,j) = ofld2d(i,j,no)
        end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)
      end do

       indfv5d= indfv5d + nout2d 
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITEAVGOUT_V5D( nouttime,vmisdatv5d,itv5d,indfv5d )

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
           
      logical zpress       
      common /press/zpress 
      real ofld3d_p,  o3davg_p
      common /outfld3d_p/ofld3d_p(nx,ny,nzp,nout3d), 
     &           o3davg_p(nx,ny,nzp,nout3d,nhrout)


      integer  nouttime(nhrout)
     &     ,i,j,k,ihr,indfv5d,no,nno,itv5d
      real tmp2d(nx,ny), tmp3d(nx,ny,nz),tmp3d_p(nx,ny,nzp)
     &   , xntimes,vmisdatv5d

      real*8 xhr1 ,xhravg

     

      print*,'COMPUTING AVERAGE OUT FIELDS:',nouttime
 

c **** WRITE OUT AVERAGED 3-D FIELDS **** c

c      CALL SETCONST(tmp3d,vmisdat,nx,ny,nz,1,1,1,nx,1,ny)
       
       nno = indfv5d

      if (.not. zpress) then  

      do no=1,nout3d
c       print*,vnamout(no)
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)
      
        nno = nno +1
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)

          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (o3davg(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + o3davg(i,j,k,no,ihr)*xntimes
              o3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do 
   
        CALL WRITEV5D(itv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
      end do
     
      else
 
       do no=1,nout3d
        CALL SETCONST(tmp3d_p,0.0,nx,ny,nzp,1,1,1,nx1,1,ny1)

        nno = nno + 1
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)

          do k=1,nzp
          do j=1,ny
          do i=1,nx
            if (o3davg_p(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k) + 
     &                         o3davg_p(i,j,k,no,ihr)*xntimes
              o3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do 
   
        CALL WRITEV5D(itv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
      end do

      end if  
      indfv5d= indfv5d + nout3d 
      

c **** WRITE OUT AVERAGED 2-D FIELDS IN VIS5D  FORMAT **** c
      nno = indfv5d
      do no=1,nout2d
        nno = nno + 1
        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrout
          xntimes = 1./float(nouttime(ihr)*nhrout)
          do j=1,ny
          do i=1,nx
            if (o2davg(i,j,no,ihr).ne.vmisdatv5d) then  
              tmp2d(i,j) = tmp2d(i,j) + o2davg(i,j,no,ihr)*xntimes
              o2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
        end do

        CALL WRITEV5D(itv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)

      end do                
      
      indfv5d= indfv5d + nout2d 
      
      return
      end
c
c***************************************************************
      SUBROUTINE WRITEDIUROUT_V5D( xhr1,nouttime,vmisdatv5d,indfv5d)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      
      logical zpress       
      common /press/zpress

      real ofld3d_p,  o3davg_p
      common /outfld3d_p/ofld3d_p(nx,ny,nzp,nout3d), 
     &           o3davg_p(nx,ny,nzp,nout3d,nhrout)

      integer nouttime(nhrout)
     &      ,no,nno,ihr,i,j,k,itv5d,indfv5d
      
      real  vmisdatv5d, xntimes
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz),tmp3d_p(nx,ny,nzp)
   
      real*8 xhr1 ,xhravg


      print*,'COMPUTING AVERAGE FIELDS FOR DIURNAL OUTPUT:',nouttime


c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
      nno = indfv5d
     
      if (.not. zpress) then
      CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)
      do no=1,nout3d
        nno=nno+1
        do ihr=1,nhrout
          xhravg = xhr1 + float(ihr-1)*dtout
          xntimes = 1./float(nouttime(ihr))
          if (nouttime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (o3davg(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = o3davg(i,j,k,no,ihr)*xntimes
              o3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d,nx,ny,nz,vmisdatv5d)

        end do
      end do
      
      else

      do no=1,nout3d
        nno=nno+1
        do ihr=1,nhrout
          xhravg = xhr1 + float(ihr-1)*dtout
          xntimes = 1./float(nouttime(ihr))
          if (nouttime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nzp
          do j=1,ny1
          do i=1,nx1
            if (o3davg_p(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = o3davg_p(i,j,k,no,ihr)*xntimes
              o3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)

        end do
      end do

      endif 

      indfv5d= indfv5d + nout3d 

c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
    
      nno = indfv5d
      do no=1,nout2d
        nno = nno+1
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
            if (o2davg(i,j,no,ihr).ne.vmisdatv5d) then
              tmp2d(i,j) = o2davg(i,j,no,ihr)*xntimes
              o2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
         CALL WRITEV5D(ihr,nno,tmp2d,nx,ny,1,vmisdatv5d)
        end do
      end do

      indfv5d= indfv5d + nout2d 

      return
      end

c
      SUBROUTINE RDCHE(iin,idate,ierr)

      implicit none
      include 'postpv5d.param'
      real tmp2d(nx,ny), xday
      integer iin,idate,ierr,nc,i,j,k,idatez

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nr2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)



      read (iin,iostat=ierr) idate
      if (ierr.ne.0) return
      print *,'READING CHEM-TRACER DATA:  ',idate
      do nc=1,nc3d
        do k=1,nz
          read (iin,iostat=ierr) tmp2d
          if (ierr.ne.0) return
          do j=1,ny
          do i=1,nx
            cfld3d(i,j,k,nc) = tmp2d(i,j)
          end do
          end do
        end do
c          print*,nc,maxval( cfld3d(:,:,:,nc))
      end do
      do nc=1,nc2d
        read (iin,iostat=ierr) tmp2d
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          cfld2d(i,j,nc) = tmp2d(i,j)
        end do
        end do
c        print*,nc,maxval( cfld2d(:,:,nc))
      end do
c     print*,'DONE CHEM FOR CURRENT TIMESTEP',idate
   
      return
      end

      

cccccccccccccccccc

      SUBROUTINE WRITECHE_V5D( vmisdatv5d,itv5d,indfv5d)

      implicit none
      include 'postpv5d.param'
     
      integer i,j,k,no,nno,indfv5d,itv5d
    
      real vmisdat,vmisdatv5d

      real  tmp2d(nx,ny), tmp3d(nx,ny,nz)      
      real tmp3d_p(nx,ny,nzp) 

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)

      logical zpress       
      common /press/zpress

      real  cfld3d_p  
      common /cheflds_p/ cfld3d_p(nx,ny,nzp,nc3d)



c **** WRITE CHE 3-D FIELDS IN NetCDF FORMAT **** c



      nno= indfv5d

      if (.not.zpress)then

       CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)
       do no=1,nc3d
        nno = nno +1
        do k=1,nz
        do j=1,ny1
        do i=1,nx1
          tmp3d(i,j,k) = cfld3d(i,j,k,no)*1e9
        end do
        end do
        end do
   
        CALL WRITEV5D(itv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
       end do

      else 
 
       CALL SETCONST(tmp3d_p,vmisdatv5d,nx,ny,nzp,1,1,1,nx,1,ny)
       do no=1,nc3d
        nno = nno +1

        do k=1,nzp
        do j=1,ny1
        do i=1,nx1
          if (cfld3d_p(i,j,k,no) .ne. vmisdatv5d) then
           tmp3d_p(i,j,k) = cfld3d_p(i,j,k,no)*1.e9
          end if          
        end do
        end do
        end do
       
        CALL WRITEV5D(itv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
    
       end do

       end if
       
      indfv5d= indfv5d + nc3d 

c **** WRITE OUT 2-D FIELDS IN NetCDF FORMAT **** c
      
      CALL SETCONST(tmp2d,vmisdatv5d,nx,ny,1,1,1,1,nx,1,ny)
      nno= indfv5d
      do no=1,nc2d
        nno = nno + 1
        do j=1,ny
        do i=2,nx1
          tmp2d(i,j) = cfld2d(i,j,no)
        end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)
      end do

      indfv5d= indfv5d + nc2d 

      return
      end

cccccccccccccccccccindfv5d
      SUBROUTINE WRITEAVGCHE_V5D(nchetime, vmisdatv5d, itv5d,indfv5d )

      implicit none
      include 'postpv5d.param'

      integer i,j,k,nc,nno,ihr,itv5d,indfv5d
      real tmp2d(nx,ny), tmp3d(nx,ny,nz),tmp3d_p(nx,ny,nzp),
     &    vmisdatv5d, xntimes
     
     
      integer  nchetime(nhrche)
   

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)

      logical zpress       
      common /press/zpress

      real  cfld3d_p,  c3davg_p  
      common /cheflds_p/ cfld3d_p(nx,ny,nzp,nc3d),
     &             c3davg_p(nx,ny,nzp,nc3d,nhrout)
     

      print*,'COMPUTING AVERAGE CHE FIELDS:',nchetime
     

c **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
   
      nno= indfv5d
      if (.not. zpress) then

      CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)
      do nc=1,nc3d
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)

        nno = nno +1
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (c3davg(i,j,k,nc,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + c3davg(i,j,k,nc,ihr)*xntimes
     &                                      *1.E9
              c3davg(i,j,k,nc,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
      end do
  
      else 

      CALL SETCONST(tmp3d_p,vmisdatv5d,nx,ny,nzp,1,1,1,nx,1,ny)
      do nc=1,nc3d
        CALL SETCONST(tmp3d_p,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)

        nno = nno +1
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do k=1,nzp
          do j=1,ny
          do i=1,nx
            if (c3davg_p(i,j,k,nc,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k) + 
     &                         c3davg_p(i,j,k,nc,ihr)*xntimes *1.E9
              c3davg_p(i,j,k,nc,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
      end do

      end if 

      indfv5d= indfv5d + nc3d 

c **** WRITE CHE  AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
   
      nno= indfv5d
      do nc=1,nc2d
        nno = nno +1

        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrche
          xntimes = 1./float(nchetime(ihr)*nhrche)
          do j=1,ny
          do i=1,nx
            if (c2davg(i,j,nc,ihr).ne.vmisdatv5d) then
              tmp2d(i,j) = tmp2d(i,j) + c2davg(i,j,nc,ihr)*xntimes
              c2davg(i,j,nc,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)
      end do

      indfv5d = indfv5d +nc2d

      return
      end
cccccccccccccccccccc

      SUBROUTINE WRITEDIURCHE_V5D(xhr1,nchetime,vmisdatv5d,indfv5d)

      implicit none
      include 'postpv5d.param'

      real cfld2d, cfld3d, c2davg, c3davg
      common /cheflds/ cfld2d(nx,ny,nc2d), cfld3d(nx,ny,nz,nc3d)
     &   , c2davg(nx,ny,nc2d,nhrche), c3davg(nx,ny,nz,nc3d,nhrche)


      logical zpress       
      common /press/zpress

      real  cfld3d_p ,  c3davg_p
      common /cheflds_p/ cfld3d_p(nx,ny,nzp,nc3d),
     &             c3davg_p(nx,ny,nzp,nc3d,nhrout)

      integer nchetime(nhrout),no,nno,ihr,i,j,k,itv5d,indfv5d
      real  sighrev(nz), vmisdatv5d, xntimes
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), tmp3d_p(nx,ny,nzp)
      
      real*8 xhr1 ,xhravg



      print*,'COMPUTING AVERAGE FIELDS FOR DIURNAL CHE OUTPUT:',nchetime

      nno = indfv5d

c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** 

      if (.not. zpress) then 

      CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)

      do no=1,nc3d
        nno=nno+1 
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))

          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (c3davg(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = c3davg(i,j,k,no,ihr)*xntimes
     &                         *1.E9
              c3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d,nx,ny,nz,vmisdatv5d)

        end do
      end do

      else

      CALL SETCONST(tmp3d_p,vmisdatv5d,nx,ny,nzp,1,1,1,nx,1,ny)

      do no=1,nc3d
        nno=nno+1 
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))

          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nzp
          do j=1,ny1
          do i=1,nx1
            if (c3davg_p(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = c3davg_p(i,j,k,no,ihr)*xntimes
     &                         *1.E9
              c3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)

        end do
      end do

      end if

      indfv5d= indfv5d + nc3d 

c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
      nno =  indfv5d     
      do no=1,nc2d        
        nno = nno + 1
        do ihr=1,nhrche
          xhravg = xhr1 + float(ihr-1)*dtche
          xntimes = 1./float(nchetime(ihr))
          print*,'nchetime(ihr)=',nchetime(ihr),'xntimes=',xntimes
     &          ,'ihr=',ihr,xhravg
          if (nchetime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (c2davg(i,j,no,ihr).ne.vmisdatv5d) then
              tmp2d(i,j) = c2davg(i,j,no,ihr)*xntimes
              c2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
         CALL WRITEV5D(ihr,nno,tmp2d,nx,ny,1,vmisdatv5d)
        end do
      end do
      
      indfv5d= indfv5d + nc2d 
      return
      end


ccccccccccccccccccccccc
      SUBROUTINE RDRAD(iin,mdate,ierr)

      implicit none
      include 'postpv5d.param'
      real tmp2d(nx,ny), xday
      integer iin,mdate,ierr,nr,i,j,k

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      integer   ncld,  nclwp,   nqrs,   nqrl
     &     ,   nfsw,   nflw, nclrst, nclrss, nclrlt
     &     , nclrls, nsolin, nsabtp, nfirtp
      common /rpoint3d/   ncld,  nclwp,   nqrs,   nqrl
      common /rpoint2d/   nfsw,   nflw, nclrst, nclrss, nclrlt
     &               , nclrls, nsolin, nsabtp, nfirtp

      read (iin,iostat=ierr) mdate
      if (ierr.ne.0) return
      print *,'reading radiation:  ',mdate
      do nr=1,nr3d
        do k=1,nz
          read (iin,iostat=ierr) tmp2d
          if (ierr.ne.0) return
          do j=1,ny
          do i=1,nx
            rfld3d(i,j,k,nr) = tmp2d(i,j)
          end do
          end do
        end do
   
      end do
      do nr=1,nr2d
        read (iin,iostat=ierr) tmp2d
        if (ierr.ne.0) return
        do j=1,ny
        do i=1,nx
          rfld2d(i,j,nr) = tmp2d(i,j)
        end do
        end do
      end do
     
c     print*,'DONE READING RADIATION FOR CURRENT TIMESTEP',xday

      return
      end

ccccccccccccccccccccccccccccccccccc
        
      SUBROUTINE WRITERAD_V5D( vmisdatv5d, itrv5d, indfv5d )

      implicit none
      include 'postpv5d.param'


      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)
      
      logical zpress       
      common /press/zpress

      real  rfld3d_p
      common /radflds_p/  rfld3d_p(nx,ny,nzp,nr3d)
    
     
      integer i, j, k, no, nno,itrv5d,indfv5d 
      real vmisdatv5d

      real tmp2d(nx,ny), tmp3d(nx,ny,nz)
      real tmp3d_p(nx,ny,nzp)
      logical outf, out,rad,chem 
      common /logic/out,rad,chem 

      

      nno=indfv5d

      if(.not. zpress) then
        do no=1,nr3d
c         print*,no,vnamout(no)
          nno = nno +1
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
           tmp3d(i,j,k) = rfld3d(i,j,k,no)
          end do
          end do
          end do 


           CALL WRITEV5D(itrv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
         end do

      else 

         do no=1,nr3d
          nno = nno + 1
          do k=1,nzp
          do j=1,ny1
          do i=1,nx1
            if (rfld3d_p(i,j,k,no) .ne. vmisdatv5d) then
             tmp3d_p(i,j,k) = rfld3d_p(i,j,k,no) 
            else
             tmp3d_p(i,j,k)=vmisdatv5d
            end if
          end do
          end do
          end do 
         CALL WRITEV5D(itrv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
         end do
          
      end if
      
      indfv5d= indfv5d + nr3d 
      

      nno = indfv5d

      do no=1,nr2d
        nno= nno +1
        do j=1,ny1
        do i=1,nx1
          tmp2d(i,j) = rfld2d(i,j,no)
        end do
        end do       
      CALL WRITEV5D(itrv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)
      end do

      indfv5d= indfv5d + nr2d 

      return
      end 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITEAVGRAD_V5D(nradtime, vmisdatv5d, itv5d,indfv5d )

      implicit none
      include 'postpv5d.param'

      integer i,j,k,nr,nno,ihr,itv5d,indfv5d
      real tmp2d(nx,ny), tmp3d(nx,ny,nz), tmp3d_p(nx,ny,nzp)
     &   , vmisdatv5d, sighrev(nz), xntimes
     
  
      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)   

      
      logical zpress       
      common /press/zpress

      real  rfld3d_p,r3davg_p
      common /radflds_p/  rfld3d_p(nx,ny,nzp,nr3d), 
     &          r3davg_p(nx,ny,nzp,nr3d,nhrout)

      integer  nradtime(nhrche)
 

      logical outf, out,rad,chem 
      common /logic/out,rad,chem 
 

      print*,'COMPUTING AVERAGE RAD FIELDS:',nradtime
     

c **** WRITE RAD AVERAGED 3-D FIELDS IN NetCDF FORMAT **** c
   
c       rad fields after (dynamic +chem)

      
      nno=indfv5d

      if (.not. zpress) then 
      
      CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)

      do nr=1,nr3d
        CALL SETCONST(tmp3d,0.0,nx,ny,nz,1,1,1,nx1,1,ny1)

        nno = nno +1
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            if (r3davg(i,j,k,nr,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = tmp3d(i,j,k) + r3davg(i,j,k,nr,ihr)*xntimes
              r3davg(i,j,k,nr,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp3d,nx,ny,nz,vmisdatv5d)
      end do
  
      else
     
      do nr=1,nr3d
        CALL SETCONST(tmp3d_p,0.0,nx,ny,nzp,1,1,1,nx1,1,ny1)

        nno = nno +1
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do k=1,nzp
          do j=1,ny
          do i=1,nx
            if (r3davg_p(i,j,k,nr,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = tmp3d_p(i,j,k) + 
     &                         r3davg_p(i,j,k,nr,ihr)*xntimes
              r3davg_p(i,j,k,nr,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)
      end do

      end if
      
      indfv5d= indfv5d + nr3d 

c **** WRITE CHE  AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
            
      nno = indfv5d
      
      do nr=1,nr2d
        nno = nno +1

        CALL SETCONST(tmp2d,0.0,nx,ny,1,1,1,1,nx1,1,ny1)
        do ihr=1,nhrrad
          xntimes = 1./float(nradtime(ihr)*nhrrad)
          do j=1,ny
          do i=1,nx
            if (r2davg(i,j,nr,ihr).ne.vmisdatv5d) then
              tmp2d(i,j) = tmp2d(i,j) + r2davg(i,j,nr,ihr)*xntimes
              r2davg(i,j,nr,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
        end do
        CALL WRITEV5D(itv5d,nno,tmp2d,nx,ny,1,vmisdatv5d)
      end do
      indfv5d = indfv5d + nr2d
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITEDIURRAD_V5D(xhr1,nradtime,vmisdatv5d,indfv5d)

      implicit none
      include 'postpv5d.param'


      integer nradtime(nhrout),no,nno,ihr,i,j,k,itv5d,indfv5d
      real  vmisdatv5d, xntimes
     &   , tmp2d(nx,ny), tmp3d(nx,ny,nz), tmp3d_p(nx,ny,nzp)
      
      real*8 xhr1 ,xhravg

  
      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)
     &   , r2davg(nx,ny,nr2d,nhrrad), r3davg(nx,ny,nz,nr3d,nhrrad)   


      logical zpress       
      common /press/zpress

      real  rfld3d_p,r3davg_p
      common /radflds_p/  rfld3d_p(nx,ny,nzp,nr3d), 
     &          r3davg_p(nx,ny,nzp,nr3d,nhrout)

      logical outf, out,rad,chem 
      common /logic/out,rad,chem  


      print*,'COMPUTING AVERAGE FIELDS FOR DIURNAL RAD OUTPUT:',nradtime

c **** WRITE OUT AVERAGED 3-D FIELDS IN NetCDF FORMAT **** 


      nno = indfv5d   

      if (.not. zpress) then 

      CALL SETCONST(tmp3d,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)

      do no=1,nr3d
        nno=nno+1 
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))

          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nz
          do j=1,ny1
          do i=1,nx1
            if (r3davg(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d(i,j,k) = r3davg(i,j,k,no,ihr)*xntimes
              r3davg(i,j,k,no,ihr) = 0.0
            else
              tmp3d(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d,nx,ny,nz,vmisdatv5d)

        end do
      end do


      else 

      CALL SETCONST(tmp3d_p,vmisdatv5d,nx,ny,nz,1,1,1,nx,1,ny)

      do no=1,nr3d
        nno=nno+1 
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))

          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do k=1,nzp
          do j=1,ny1
          do i=1,nx1
            if (r3davg_p(i,j,k,no,ihr).ne.vmisdatv5d) then
              tmp3d_p(i,j,k) = r3davg_p(i,j,k,no,ihr)*xntimes
              r3davg_p(i,j,k,no,ihr) = 0.0
            else
              tmp3d_p(i,j,k) = vmisdatv5d
            end if
          end do
          end do
          end do

         CALL WRITEV5D(ihr,nno,tmp3d_p,nx,ny,nzp,vmisdatv5d)

        end do
      end do

      end if  

       indfv5d= indfv5d + nr3d 


c **** WRITE OUT AVERAGED 2-D FIELDS IN NetCDF FORMAT **** c
       nno = indfv5d
 
      do no=1,nr2d        
        nno = nno + 1
        do ihr=1,nhrrad
          xhravg = xhr1 + float(ihr-1)*dtrad
          xntimes = 1./float(nradtime(ihr))
          print*,'nradtime(ihr)=',nradtime(ihr),'xntimes=',xntimes
     &          ,'ihr=',ihr,xhravg
          if (nradtime(ihr).le.0) then
            print*,'NOTHING TO AVERAGE -- nouttime = 0'
            stop 999
          end if
          do j=1,ny
          do i=1,nx
            if (r2davg(i,j,no,ihr).ne.vmisdatv5d) then
              tmp2d(i,j) = r2davg(i,j,no,ihr)*xntimes
              r2davg(i,j,no,ihr) = 0.0
            else
              tmp2d(i,j) = vmisdatv5d
            end if
          end do
          end do
         CALL WRITEV5D(ihr,nno,tmp2d,nx,ny,1,vmisdatv5d)
        end do
      end do

      indfv5d= indfv5d + nr2d 

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



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


      SUBROUTINE AVGDATA3D(favg,f,n1,n2,n3,n4,n5,ihr,vmisdatv5d)

      implicit none
      integer i,j,k,l,n1,n2,n3,n4,n5,ihr
      real favg(n1,n2,n3,n4,n5), f(n1,n2,n3,n4), vmisdatv5d

      do l=1,n4
      do k=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,k,l).ne.vmisdatv5d) then
          favg(i,j,k,l,ihr) = favg(i,j,k,l,ihr) + f(i,j,k,l)
        else
          favg(i,j,k,l,ihr) = vmisdatv5d
        end if
      end do
      end do
      end do
      end do

      return
      end

      SUBROUTINE AVGDATA2D(favg,f,n1,n2,n3,n4,ihr,vmisdatv5d)

      implicit none
      integer i,j,l,n1,n2,n3,n4,ihr
      real favg(n1,n2,n3,n4), f(n1,n2,n3), vmisdatv5d

      do l=1,n3
      do j=1,n2
      do i=1,n1
        if (f(i,j,l).ne.vmisdatv5d) then
          favg(i,j,l,ihr) = favg(i,j,l,ihr) + f(i,j,l)
        else
          favg(i,j,l,ihr) = vmisdatv5d
        end if
      end do
      end do
      end do

      return
      end



      SUBROUTINE MMVLUOUT(vnamout,lnamout,uout,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postpv5d.param'
      character vnamout(notot)*10, lnamout(notot)*20, uout(notot)*13
      real xmin(notot), xmax(notot), fact(notot), offset(notot), aaa
      integer l

      integer nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt, nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf,ntopo,nlus
      common /opoint3d/nua,  nva,  nta, nqva, nqca, nmse,nrh, nhgt, nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf,ntopo,nlus

      vnamout(nua)  = 'U'
      vnamout(nva)  = 'V'
      vnamout(nta)  = 'TK'
      vnamout(nqva) = 'QD'
      vnamout(nqca) = 'QC'
      vnamout(nout3d+npsa) = 'PS'
      vnamout(nout3d+nrt)  = 'RT'
      vnamout(nout3d+ntgb) = 'TGRND'
      vnamout(nout3d+nsmt) = 'SMT'
      vnamout(nout3d+nbf)  = 'RB'
      vnamout(nout3d+ntopo)  = 'TOPOG'
      vnamout(nout3d+nlus)  = 'LUSE'
      vnamout(nmse) = 'MSE'
      vnamout(nrh)  = 'RH'
      vnamout(nhgt) = 'HGT'
      vnamout(nsig) = 'SIGM'
      

      lnamout(nua)  = 'Zonal Wind'
      lnamout(nva)  = 'Meridional Wind'
      lnamout(nta)  = 'Temperature'
      lnamout(nqva) = 'Mixing Ratio'
      lnamout(nqca) = 'Cloud Mixing Ratio'
      lnamout(nout3d+npsa) = 'Surface Pressure'
      lnamout(nout3d+ntgb) = 'Ground Temperature'
      lnamout(nout3d+nrt)  = 'Total Precip'
      lnamout(nout3d+nsmt) = 'Total Soil Water'
      lnamout(nout3d+nbf)  = 'Base Flow'
      lnamout(nout3d+ntopo)  = 'Model topog'
      lnamout(nout3d+nlus)  = 'Model land use'
      lnamout(nmse) = 'Moist Static Energy'
      lnamout(nrh)  = 'Relative Humidity'
      lnamout(nhgt) = 'Geopotential Height'
      lnamout(nsig) = 'Sigma level'

      uout(nua)  = 'm/s'
      uout(nva)  = 'm/s'
      uout(nta)  = 'K'
      uout(nqva) = 'kg/kg'
      uout(nqca) = 'kg/kg'
      uout(nout3d+npsa) = 'hPa'
      uout(nout3d+ntgb) = 'K'
      uout(nout3d+nrt)  = 'mm/day'
      uout(nout3d+nsmt) = 'mm'
      uout(nout3d+nbf)  = 'mm/day'
      uout(nout3d+ntopo)  = 'm'
      uout(nout3d+nlus)  = ''
      uout(nmse) = 'm2/s2'
      uout(nrh)  = 'fraction'
      uout(nhgt) = 'm'
      uout(nsig)  = 'fraction'

      xmax(nua)  = 210.0
      xmax(nva)  = 210.0
      xmax(nta)  = 350.0
      xmax(nqva) = 0.05
      xmax(nqca) = 0.005
      xmax(nout3d+npsa) = 1200.0
      xmax(nout3d+ntgb) = 350.0
      xmax(nout3d+nrt)  = 1000.0
      xmax(nout3d+nsmt) = 3000.0
      xmax(nout3d+nbf)  = 200.0
      xmax(nmse) = 550000.0
      xmax(nrh)  = 30.0
      xmax(nhgt) = 40000.0

      xmin(nua)  = -210.0
      xmin(nva)  = -210.0
      xmin(nta)  = 160.0
      xmin(nqva) = -0.001
      xmin(nqca) = -0.001
      xmin(nout3d+npsa) = 300.0
      xmin(nout3d+ntgb) = 200.0
      xmin(nout3d+nrt)  = -10.0
      xmin(nout3d+nsmt) = 0.0
      xmin(nout3d+nbf)  = -10.0
      xmin(nmse) = 200000.0
      xmin(nrh)  = -0.5
      xmin(nhgt) = -100.0

      aaa = 2.**16.-1.
      do l=1,notot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end


      SUBROUTINE MMVLURAD(vnamrad,lnamrad,urad,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postpv5d.param'
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
      xmax(nclwp)       = 1500.0
      xmax(nqrs)        = 1.0e-3
      xmax(nqrl)        = 1.0e-3
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
      xmin(nqrs)        = -1.0e-3
      xmin(nqrl)        = -1.0e-3
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


      SUBROUTINE MMVLUCHE(vnamche,lnamche,uche,xmin,xmax
     &         , fact,offset)

      implicit none
      include 'postpv5d.param'
      character vnamche(nctot)*10, lnamche(nctot)*20
     &         , uche(nctot)*13
      real xmin(nctot), xmax(nctot),fact(nctot),offset(nctot),
     &                 aaa

      integer naer1, naer2
      common /chepoint3d/ naer1, naer2
      
      character tracname(ntrac)*10  
      integer l,n,r

      tracname(1) = 'TR1'
      tracname(2) = 'TR2'        
      tracname(3) = 'TR3'  
      tracname(4) = 'TR4'
      tracname(5) = 'TR5'
      tracname(6) = 'TR6'
      
      r = nc3d/ntrac
      do n=1,ntrac             
       vnamche(r*(n-1)+1) = tracname(n)
       lnamche(r*(n-1)+1 ) = 'MMR_'//tracname(n)
       uche(r*(n-1)+1)    = 'micro-g/Kg'
      end do     


      r= nc2d/ntrac
      do n=1,ntrac
       vnamche(nc3d+r*(n-1)+1)='BURD'//tracname(n) 
       vnamche(nc3d+r*(n-1)+2)='WDLS'//tracname(n) 
       vnamche(nc3d+r*(n-1)+3)='WDCV'//tracname(n)  
       vnamche(nc3d+r*(n-1)+4)='DRDP'//tracname(n)  
       vnamche(nc3d+r*(n-1)+5)='EMRA'//tracname(n)  

       lnamche(nc3d+r*(n-1)+1)='Int column '//tracname(n) 
       lnamche(nc3d+r*(n-1)+2)='Wetdep lsc '//tracname(n) 
       lnamche(nc3d+r*(n-1)+3)='Wetdep cvc '//tracname(n)  
       lnamche(nc3d+r*(n-1)+4)='Drydep surf'//tracname(n)  
       lnamche(nc3d+r*(n-1)+5)='Emission rate'//tracname(n)  

       uche(nc3d+r*(n-1)+1)='mg/m2'
       uche(nc3d+r*(n-1)+2)='mg/m2' 
       uche(nc3d+r*(n-1)+3)='mg/m2' 
       uche(nc3d+r*(n-1)+4)='mg/m2'  
       uche(nc3d+r*(n-1)+5)='micro-g/m2.s'  

      end do 
   
      
      aaa = 2.**16.-1.
      do l=1,nctot
        fact(l)=(xmax(l)-xmin(l))/aaa
        offset(l)=(xmax(l)+xmin(l))/2.
      end do

      return
      end


      SUBROUTINE CALCMSE(zs,sigh,pt)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)

      real tv(nx,ny,nz), p(nx,ny,nz), z(nx,ny,nz)
      real zs(nx,ny), sigh(nz)
	real g, r, cp, ps, lh, ta, qa, za, t, q, pt, dp, tvbar
	integer i, j, k, k1

      integer nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt, nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf
      common /opoint3d/  nua,  nva,  nta, nqva, nqca,nmse,nrh,nhgt,nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf

      parameter (g=9.81, r=287., cp=1004. ,lh=2.5e06)

      do j=1,ny1
      do i=1,nx1
        k = 1
        t = ofld3d(i,j,k,nta)
        q = ofld3d(i,j,k,nqva)
        ps = ofld2d(i,j,npsa)
        tv(i,j,k)  = t*(1.+0.608*q)
        p(i,j,k)   = (ps-pt*10.)*sigh(k) + pt*10.
        dp         = log(p(i,j,k)) - log(ps)
        z(i,j,k)   = zs(i,j) - dp*tv(i,j,k)*r/g
        ofld3d(i,j,k,nmse) = g*z(i,j,k) + cp*t + lh*q
      end do
      end do
      do k=2,nz
      do j=1,ny1
      do i=1,nx1
        k1 = k-1
        t = ofld3d(i,j,k,nta)
        q = ofld3d(i,j,k,nqva)
        ps = ofld2d(i,j,npsa)
        tv(i,j,k)  = t*(1.+0.608*q)
        tvbar      = 0.5*(tv(i,j,k)+tv(i,j,k1))
        p(i,j,k)   = (ps-pt*10.)*sigh(k) + pt*10.
        dp         = log(p(i,j,k)) - log(p(i,j,k1))
        z(i,j,k)   = z(i,j,k1) - dp*tvbar*r/g
        ofld3d(i,j,k,nmse) = g*z(i,j,k) + cp*t + lh*q
      end do
      end do
      end do

      return
      end



      SUBROUTINE CALCRH(sigh,pt)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      real sigh(nz)

      integer nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt,nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf
      common /opoint3d/  nua,  nva,  nta, nqva, nqca,nmse,nrh,nhgt,nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf

	integer i, j, k
	real tfrz, tr, pres, pt, t, q, hl, satvp, qs

      parameter (tfrz = 273.16, tr=1.0/tfrz)


      do k=1,nz
      do j=1,ny1
      do i=1,nx1
        pres = (ofld2d(i,j,npsa)-pt*10.)*sigh(k) + pt*10. ! PRES AT LEVEL K
        t = ofld3d(i,j,k,nta)
        q = ofld3d(i,j,k,nqva)
        hl = 597.3 - 0.566*(t-tfrz)                     ! LATENT HEAT OF EVAP.
        satvp = 6.11*exp(9.045*hl*(tr-1.0/t))           ! SATURATION VAP PRESS.
        qs = 0.622*satvp/(pres-satvp)                   ! SAT. MIXING RATIO
        ofld3d(i,j,k,nrh) = q/qs
      end do
      end do
      end do

      return
      end



      SUBROUTINE CALCHGT(zs,sigf,sigh,pt)

      implicit none
      include 'postpv5d.param'

      real ofld2d, ofld3d, o2davg, o3davg
      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)

      integer nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt,nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf
      common /opoint3d/  nua,  nva,  nta, nqva, nqca,nmse,nrh,nhgt,nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf

	integer i, j, k, k1
      real zs(nx,ny), sigh(nz), dsig(nz), sigf(nz+1)
	real pt, ep1, rg, g, rgog, ps, t, q, pf, tv, t1, q1, tvbar, tv2

      ep1 = 0.608
      RG=287.04
      g = 9.805
      rgog = rg/g

      do K=1,nz
        DSIG(K)=SIGF(K)-SIGF(K+1)
      end do

      do i=1,nx1
      do j=1,ny1
        ps = ofld2d(i,j,npsa)
        t = ofld3d(i,j,1,nta)
        q = ofld3d(i,j,1,nqva)
        PF=PT/(PS/10.-PT)
        Tv = T*(1.0+ep1*Q)
        ofld3d(i,j,1,nhgt)=ZS(I,J)+Tv*rgog*LOG((1.+PF)/(sigh(1)+PF))
        DO K=2,nz
          K1=K-1
          t = ofld3d(i,j,k,nta)
          t1 = ofld3d(i,j,k1,nta)
          q = ofld3d(i,j,k,nqva)
          q1 = ofld3d(i,j,k1,nqva)
          Tv = T*(1.0+ep1*Q)
          Tv2 = T1*(1.0+ep1*Q1)
          TvBAR=(Tv*DSIG(K)+Tv2*DSIG(K1))/(DSIG(K)+DSIG(K1))
          ofld3d(i,j,k,nhgt)=ofld3d(i,j,k1,nhgt)
     &              +TvBAR*rgog*LOG((sigh(K1)+PF)/(sigh(K)+PF))
        end do
      end do
      end do
      
      return
      end

      SUBROUTINE JULIAN(idate,julnc,julday,iyr,imo,idy,ihr)

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




      SUBROUTINE PARAM(nx,ny,kz,ds,ddeg,clat,clon,xplat,xplon
     &         , xlat,xlon,vvarmin,vvarmax,xlat1d,xlon1d,idim,ndim)

      real   vvarmin(ndim), vvarmax(ndim), xlat1d(ny), xlon1d(nx)
      real xlat(nx,ny), xlon(nx,ny)
      integer idim(ndim)

      ddeg=ds*180./3.14159265359/6370000.
      print*,'DS=',ds,'DDEG=',ddeg
      vvarmin(1) = xlon(1,ny/2)
      vvarmin(2) = xlat(nx/2,1)
      vvarmin(3) = 1050.
      vvarmax(1) = xlon(nx,ny/2)
      vvarmax(2) = xlat(nx/2,ny)
      vvarmax(3) = 0.
      idim(1) = nx
      idim(2) = ny
      idim(3) = kz
      do i=1,nx
        xlon1d(i) = xlon(i,ny/2)
      end do
      do j=1,ny
        xlat1d(j) = xlat(nx/2,j)
      end do

      return
      end


c------------------------------------------------------------------------
c-----------------------------------------------------------------------    
      SUBROUTINE PREP_V5D(ihr1,ihr2,julda1,julda2,idat,dt,numtimes,
     &                     numvars,
     &                    nl,
     &                    vname,lname,unit,clat,clon,xlat,xlon,
     &                    ds,nzlev,vertgrid, v5dfile)

      implicit none
      include 'postpv5d.param'
      include "v5df.h"
      
        
      integer numtimes, numvars,julda1,julda2,ihr1,ihr2,nf3d 
      integer idat
      integer nzlev
      real clat,clon,ds,dt
      real xlat(nx,ny)
      real xlon(nx,ny)
    
   
      integer projection,vertical


      real latn,lonw,dis  
      integer szproj   
c      integer proj_args(6) 
      real proj_args(4) 
      integer nl(numvars) 
      integer vertgrid(nzlev)
      integer dates(numtimes)
      integer times(numtimes)
c      date begining output    
       integer julday,ihr 
       character vname(numvars)*10, lname(numvars)*20,unit(numvars)*13
       character*7 v5dfile  

       integer compressmode 
       data compressmode / 1 /

      integer n,i
     
      
      ihr=ihr1
      julday=julda1
      
      
        do i=1,numtimes         
           if ( mod(ihr,24).eq. 0)then
             ihr=0
             julday=julday+1
           end if
           times(i)= ihr*10000                    
           dates(i)= idat*1000+julday           
           ihr=ihr + int(dt)           
        end do

         
c projection   

c        projection=2
        projection=5
        latn = xlat(1, ny1)
        lonw = xlon(1, ny1)

        dis  = ds/1000.
                  
        szproj = 6
 
        print*, 'XLAT', xlat(1, ny1),xlat(1,1)
        print*, 'XLAT', xlat(nx1, ny1)
        print*, 'XLON', xlon(nx1,1),xlon(1,1) 
         print*, 'XLON', xlon(nx1,ny1),xlon(1,ny1) 
         proj_args(1)= 51.74
         proj_args(2)= 34.24
         proj_args(3)= (xlat(1, ny1) - xlat(1,1))/ ny1
         proj_args(4)= (xlon(nx1, 1) - xlon(1,1))/ nx1

c         proj_args(3)= (30.35 +7.18)/80.  
c         proj_args(4)= (34.9228  +24.4292)/121.
c        call maklamb (latn, lonw, clon, clat, dis,szproj, proj_args)  

         print*, proj_args
c vertical (cf v5d doc)
       
      vertical = 3

      n = v5dcreate( v5dfile ,numtimes, numvars, ny1, nx1, nl,
     &               vname, times, dates, compressmode,
     &               projection, proj_args, vertical, vertgrid )


      do n=1,numvars            
        call v5dsetunits(n,unit(n))
      end do   

      END 
c ---------------------------------------------------
      subroutine maklamb(latn, lonw, clon, clat, dis,szproj,proj_args)
      implicit none
      real lonw, latn
      real clon, clat
      real dis
      integer  szproj 
      real   proj_args(szproj)

c     local 
      real pi,theta,theta1,theta2,htheta,htheta1,htheta2
      real cone,rad
      real radius, degtorad, alpha
      real stdlat1, stdlat2, ratio, x, y, temp
      data radius /6371./
c
       stdlat1 = 60
       stdlat2 = 30

c  make sure stdlat1 is the larger number

      if (stdlat1 .lt. stdlat2) then
        temp = stdlat1
        stdlat1 = stdlat2
        stdlat2 = temp
      endif

c  put map parameters in as defined by vis5d; 
c    south is positive direction, so negate ypole
 
       pi = 4.0*atan( 1.0 )
       degtorad = pi/180.0
 
       theta  = (90.0-latn) * degtorad
       theta1 = (90.0-stdlat1) * degtorad
       theta2 = (90.0-stdlat2) * degtorad
       htheta = theta /2.0
       htheta1= theta1/2.0
       htheta2= theta2/2.0
       cone=alog(sin(theta1)/sin(theta2))/
     & alog(tan(htheta1)/tan(htheta2))
      
       alpha  = cone * (lonw - clon) * degtorad
           
       rad  = radius * sin(theta1) * 
     &  (tan(htheta)/tan(htheta1))**cone/cone
 
       x = rad * sin(alpha) / dis
       y = rad * cos(alpha) / dis
            
       proj_args(1)=stdlat1
       proj_args(2)=stdlat2
       proj_args(3)=-y
       proj_args(4)=-x
       proj_args(5)=-clon
       proj_args(6)=dis

      end
c-------------------------------------------------------------------

      SUBROUTINE WRITEV5D(it,no,arr,mx,my,mz,vmisdatv5d)
     
      implicit none
      include 'postpv5d.param'
      include "v5df.h"

      integer no,mx,my,mz
      real dt,vmisdatv5d
      real arr(mx,my,mz)
      integer n,it,i,j,k,jj,ii
      

      real G(my,mx,mz)
                  
            do i=1,mx
              ii = (mx) - i +1
              do j=1,my            
                jj = (my) - j +1
                do k=1,mz 
                if (arr(j,ii,k).lt.-1.e+29) arr(j,ii,k)=vmisdatv5d
                G(j,i,k) = arr(i,jj,k)
               end do
              end do
             end do 
            where( G(:,:,:)< 1.E-30 .and. G(:,:,:)>-1.E-30 )
               G=0.
            end where

c      write the 3-D grid to the v5d file
            n = v5dwrite( it, no, G )
            if (n .eq. 0) then
c              error
               call exit(1)
            endif
         
      END
c ---------------------------------------
c       
      SUBROUTINE CLSV5D()

      implicit none
      

      include 'postpv5d.param'
      include "v5df.h"
       
      integer n
c     close the v5d file and exit
      n = v5dclose()
c      if (n .eq. 0) then
c        failed
c         call exit(1)
c      else
c        success
c         call exit(0)
c      endif

      return 
      end
c-------------------------------------------------------------
      subroutine VERT_INT(field,ht,sig,vmisdatv5d)
      

      implicit none 
      include 'postpv5d.param'
      
      
      character*3 field
      real ht      (nx,ny)
      real sig     (nz)     
      integer vmisdatv5d
      real ofld2d, ofld3d, o2davg, o3davg

      common /outflds/ ofld2d(nx,ny,nout2d), ofld3d(nx,ny,nz,nout3d)
     &   , o2davg(nx,ny,nout2d,nhrout), o3davg(nx,ny,nz,nout3d,nhrout)
      
      real ofld3d_p,  o3davg_p
      common /outfld3d_p/ofld3d_p(nx,ny,nzp,nout3d), 
     &           o3davg_p(nx,ny,nzp,nout3d,nhrout)
cccccccccccccccccccccccccccc
      real cfld2d, cfld3d  
      common /cheflds/ cfld2d(nx,ny,nc2d),cfld3d(nx,ny,nz,nc3d)

      real  cfld3d_p  
      common /cheflds_p/ cfld3d_p(nx,ny,nzp,nc3d)
cccccccccccccccccccccccccccc

      real rfld2d, rfld3d, r2davg, r3davg
      common /radflds/ rfld2d(nx,ny,nr2d), rfld3d(nx,ny,nz,nr3d)

      real  rfld3d_P
      common /radflds_p/  rfld3d_p(nx,ny,nzp,nr3d)


      integer   nua,  nva,  nta, nqva, nqca, nmse,  nrh, nhgt,nsig
     &      , npsa,  nrt, ntgb, nsmt,  nbf

      common /opoint3d/  nua,  nva,  nta, nqva, nqca, nmse,nrh,nhgt,nsig
      common /opoint2d/ npsa,  nrt, ntgb, nsmt,  nbf
      
c---------------------------------------------------      
c     local variables
      integer i,j,k,kk,itr
      real ta(nx1,ny1,nz), hsp(nx1,ny1,nz),w3d(nx1,ny1,nz) 

      real psa(nx1,ny1), tga(nx1,ny1), 
     &     slp1(nx1,ny1), slp2(nx1,ny1)  

           
      real hp(nx1,ny1,nzp),w3d_p(nx1,ny1,nzp) 
           
      print*,'INTERPOLATION ', field
     
C pstar for interpolations !  

        psa(:,:)= ofld2d(1:nx1,1:ny1,npsa)/ 10. - 8.

c      a changer avec temperature de surface
        tga(:,:) = ofld3d(1:nx1,1:ny1,1,nta)

c      temperature field used, care must be consistent with sig 
       do k=1,nz   
          kk=nz -k +1
          ta(:,:,k)= ofld3d(1:nx1,1:ny1,kk,nta)
       end do 
             
     
      if ( field .eq.'out') then
 
      CALL HTSIG(ta,hsp,psa,ht,sig,nx1,ny1,nz)
     
C 
C     to calculate Sea-Level Pressure using
C         1. ERRICO's solution described in height
C         2. a simple formulae
C         3. MM5 method
      CALL SLPRES(hsp,ta,psa,ht,tga,slp1,slp2,sig,nx1,ny1,nz)
c
C     to interpolate H,U,V,T,Q and QC
C        1. For Heights
      CALL HEIGHT(hp,hsp,ta,psa,ht,sig,nx1,ny1,nz,plev,nzp)

C        2. For Zonal and Meridional Winds
      do k=1,nz   
         kk=nz -k +1
         w3d(:,:,k)=  ofld3d(1:nx1,1:ny1,kk,nua) 
      end do 
      CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
    
      do k=1,nzp
       ofld3d_p(1:nx1,1:ny1,k,nua)= w3d_p(:,:,k)
      end do
 
      do k=1,nz   
         kk=nz -k +1                  
         w3d(:,:,k)=  ofld3d(1:nx1,1:ny1,kk,nva)
      end do
      CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
      do k=1,nzp
        ofld3d_p(1:nx1,1:ny1,k,nva)= w3d_p(:,:,k)
      end do
     


C        3. For Temperatures
      
      CALL INTLOG(w3d_p,ta,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
      do k=1,nzp
        ofld3d_p(1:nx1,1:ny1,k,nta)= w3d_p(:,:,k)
      end do

c 
C        4. For Moisture qva & qca
      
c      CALL HUMID1(ta,qva,psa,sig,ix,jx,kx)
             
c     relative humidity
      do k=1,nz
        kk=nz -k +1
        w3d(:,:,k)=  ofld3d(1:nx1,1:ny1,kk,nrh)
      end do
      CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
      do k=1,nzp
        ofld3d_p(1:nx1,1:ny1,k,nrh)= w3d_p(:,:,k)
      end do


c      cloud mixing ratio
      do k=1,nz
        kk=nz -k +1
        w3d(:,:,k)=  ofld3d(1:nx1,1:ny1,kk,nqca)
      end do    
      CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
      do k=1,nzp
        kk=nzp -k +1
        ofld3d_p(1:nx1,1:ny1,k,nqca)= w3d_p(:,:,k)      
      end do
       
c      CALL HUMID2(tp,qvp,plev,ix,jx,knp)
c      CALL INTLIN(qcp,qca,psa,sig,ix,jx,kx,plev,knp 
       
c     sigma level on pressure grid 
       do k=1,nzp
                
        ofld3d_p(1:nx1,1:ny1,k,nsig)=(plev(k)-8)
     &                   / psa(1:nx1,1:ny1) 
         
       end do


      end if

c    traceurs

      if (field .eq.'che') then  

      do itr=1,nc3d 
       do k=1,nz   
         kk=nz -k +1
         w3d(:,:,k)=  cfld3d(1:nx1,1:ny1,kk,itr)
       end do 
       
        CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
        
        do k=1,nzp
          kk=nzp -k +1
          cfld3d_p(1:nx1,1:ny1,k,itr)= w3d_p(:,:,k)      
        end do
      end do

      end if

c--------------
      if (field .eq.'rad') then  
      do itr=1,nr3d 
       do k=1,nz   
         kk=nz -k +1
         w3d(:,:,k)=  rfld3d(1:nx1,1:ny1,kk,itr)
       end do 

        CALL INTLIN(w3d_p,w3d,psa,sig,nx1,ny1,nz,plev,nzp,vmisdatv5d)
        
        do k=1,nzp
          kk=nzp -k +1
          rfld3d_p(1:nx1,1:ny1,k,itr)= w3d_p(:,:,k)      
        end do
      end do
      end if

      END



c ----------------------------------------------------------------------------
      SUBROUTINE HTSIG(T,H,PSTAR,HT,SIG,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL  T(IM,JM,KM),H(IM,JM,KM)
      REAL  PSTAR(IM,JM),HT(IM,JM)
      REAL  SIG(KM)
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K
      REAL  TBAR
C
      DO J=1,JM
      DO I=1,IM
         H(I,J,KM) = HT(I,J) + RGAS/GRAV*T(I,J,KM)
     &             * LOG((PSTAR(I,J)+PTOP)/(PSTAR(I,J)*SIG(KM)+PTOP))
      ENDDO
      ENDDO
      DO K=KM-1,1,-1
      DO J=1,JM
      DO I=1,IM
         TBAR = 0.5*( T(I,J,K)+T(I,J,K+1) )
         H(I,J,K) = H(I,J,K+1) +RGAS/GRAV*TBAR
     &            * LOG((PSTAR(I,J)*SIG(K+1)+PTOP)
     &                 /(PSTAR(I,J)*SIG(K)+PTOP))
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END

c -----------------------------------------------------------------------------
      SUBROUTINE SLPRES(H,T,PSTAR,HT,TG,SLP1,SLP2,SIG,IM,JM,KM)
      implicit none
      INTEGER IM,JM,KM
      REAL  T(IM,JM,KM),H(IM,JM,KM)
      REAL  PSTAR(IM,JM),HT(IM,JM),TG(IM,JM)
      REAL  SLP1(IM,JM),SLP2(IM,JM)
      REAL  SIG(KM)
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER KBC,I,J,K
      REAL  TSFC
C
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
      DO J=1,JM
         DO I=1,IM
            TSFC = T(I,J,KBC)-TLAPSE*(H(I,J,KBC)-HT(I,J))
            SLP1(I,J) = (PSTAR(I,J)+PTOP)
     &      * EXP( -GRAV/(RGAS*TLAPSE)*LOG(1.-HT(I,J)*TLAPSE/TSFC))
         ENDDO
      ENDDO

      DO J=1,JM
         DO I=1,IM
            SLP2(I,J) = (PSTAR(I,J)+PTOP)
     &      * EXP( GRAV*HT(I,J)/(RGAS*0.5*(TG(I,J)+288.15)))
         ENDDO
      ENDDO

      RETURN
      END

c-----------------------------------------------------------------------------
      SUBROUTINE INTLIN(FP,F,PSTAR,SIG,IM,JM,KM,P,KP,vmisdatv5d)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL  FP(IM,JM,KP),F(IM,JM,KM)
      REAL  PSTAR(IM,JM)
      REAL  SIG(KM),P(KP)
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P
      REAL  SIGP,WP,W1,vmisdatv5d
C
C  INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
C        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
C        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
    
      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / PSTAR(I,J)
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
c  under surface levels                     
            IF (SIGP .GT. 1.) FP(I,J,N)= vmisdatv5d
           
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
c ------------------------------------------------------------------------------
c ----------------------------------------------------------------------------
      SUBROUTINE INTLOG(FP,F,PSTAR,SIG,IM,JM,KM,P,KP,vmisdatv5d)
      implicit none
      INTEGER IM,JM,KM,KP
      REAL  FP(IM,JM,KP),F(IM,JM,KM)
      REAL  PSTAR(IM,JM)
      REAL  SIG(KM),P(KP)
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      INTEGER I,J,K,N
      INTEGER K1,K1P,KBC
      REAL  SIGP,WP,W1,vmisdatv5d
C
C  INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
C        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
C        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
C        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
C        CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
C        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
C        TWO EXTREME TEMPERATURES IN THE LAYER.

C
C** FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      DO K=1,KM
        IF(SIG(K).LT.BLTOP) KBC = K
      ENDDO
      DO J=1,JM
        DO I=1,IM
          DO N=1,KP
            SIGP = (P(N)-PTOP) / PSTAR(I,J)
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
C                  ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
            ENDIF
c  under surface levels                     
            IF (SIGP .GT. 1.) FP(I,J,N)= vmisdatv5d            
          ENDDO
        ENDDO
      ENDDO
      
      RETURN
      END
c------------------------------------------------------------------------------
c -----------------------------------------------------------------------------
      BLOCK DATA
      implicit none
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      DATA PTOP/8./
      DATA RGAS/287.04/
      DATA GRAV/9.80616/
      DATA BLTOP/.96/
      DATA TLAPSE/-6.5E-3/
      END
c---------------------------------------------------------------------------
      SUBROUTINE HEIGHT(HP,H,T,PSTAR,HT,SIG,IM,JM,KM,P,KP)

C  HEIGHT DETERMINES THE HEIGHT OF PRESSURE LEVELS.
C     ON INPUT:
C        H AND T ARE HEIGHT AND TEMPERATURE ON SIGMA, RESPECTIVELY.
C        PSTAR = SURFACE PRESSURE - MODEL TOP PRESSURE.
C        SIG = SIGMA LEVELS.
C        P = PRESSURE LEVELS DESIRED.
C     ON OUTPUT:
C        ALL FIELDS EXCEPT H ARE UNCHANGED.
C        H HAS HEIGHT FIELDS AT KP PRESSURE LEVELS.
C
C  FOR UPWARD EXTRAPOLATION, T IS CONSIDERED TO HAVE 0 VERITCAL DERIV.
C  FOR DOWNWARD EXTRAPOLATION, T HAS LAPSE RATE OF TLAPSE (K/KM)
C     AND EXTRAPOLATION IS DONE FROM THE LOWEST SIGMA LEVEL ABOVE
C     THE BOUNDARY LAYER (TOP ARBITRARILY TAKEN AT SIGMA = BLTOP).
C     EQUATION USED IS EXACT SOLUTION TO HYDROSTATIC RELATION,
C     GOTTEN FROM R. ERRICO (ALSO USED IN SLPRES ROUTINE):
C      Z = Z0 - (T0/TLAPSE) * (1.-EXP(-R*TLAPSE*LN(P/P0)/G))
C
      implicit none
      INTEGER IM,JM,KM,KP
      REAL  T(IM,JM,KM),H(IM,JM,KM),HP(IM,JM,KP)
      REAL  PSTAR(IM,JM),HT(IM,JM)
      REAL  SIG(KM),P(KP)
      REAL  PTOP,RGAS,GRAV,BLTOP,TLAPSE
      COMMON /CONST/ PTOP,RGAS,GRAV,BLTOP,TLAPSE
      REAL  PSIG(20)
      INTEGER I,J,K,KBC,N,KT,KB
      REAL  PSFC,TEMP,WT,WB
C
      DO K=1,KM
         IF(SIG(K).LT.BLTOP) THEN
           KBC=K
         ENDIF
      ENDDO
      PRINT *,'FIRST SIGMA LEVEL ABOVE BNDY LAYER:', SIG(KBC)
C


      DO J=1,JM
        DO I=1,IM
          DO K=1,KM
            PSIG(K) = SIG(K) * PSTAR(I,J) + PTOP
          ENDDO
          PSFC = PSTAR(I,J) + PTOP
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
C
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
c-----------------------------------------------------------------------------
