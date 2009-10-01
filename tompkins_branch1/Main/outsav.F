ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outsav(iutl)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine writes the model output for restart.            c
c                                                                     c
c     iutl : is the output unit number for large-domain variables.    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer iutl
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/parrad.cb'
      include 'Commons/param1.cb'
#ifdef DIAG
      include 'Commons/diagnosis.cb'
#endif
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
      include 'Commons/iunits.cb'
#ifdef MPP1
      include 'Common2/pmoist.cb'
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
      include 'Common2/bats.cb2'
      include 'Common2/rad.cb'
      include 'Common2/bdycod.cb'
      include 'Common2/trachem.cb'
      include 'Common2/io.cb'
#else
      include 'Commons/pmoist.cb'
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
      include 'Commons/bats.cb2'
      include 'Commons/rad.cb'
      include 'Commons/bdycod.cb'
      include 'Commons/trachem.cb'
#endif
      include 'Commons/date.cb'
      real*8  absnxt,              ! Nearest layer absorptivities
     $     abstot,              ! Non-adjacent layer absorptivites
     $     emstot               ! Total emissivity
#ifdef MPP1
      common /radbuf/ absnxt(plond,plev,4,jxp)
     &      , abstot(plond,plevp,plevp,jxp),emstot(plond,plevp,jxp)
      real*8  absnxt_io
      real*8  abstot_io
      real*8  emstot_io
      common /radbuf_io/absnxt_io(plond,plev,4,mjx-1)
     &                 ,abstot_io(plond,plevp,plevp,mjx-1)
     &                 ,emstot_io(plond,plevp,mjx-1)
#else
      common /radbuf/ absnxt(plond,plev,4,jlx)
     &      , abstot(plond,plevp,plevp,jlx),emstot(plond,plevp,jlx)
#endif
      integer m
      real*8  am,an,dstor,hstor
#ifdef MPP1
      common /dpass/ m(nsplit)
     &             , am(kx,nsplit), an(nsplit)
     &             , dstor(ix,0:jxp+1,nsplit), hstor(ix,0:jxp+1,nsplit)
#else
      common /dpass/ m(nsplit)
     a             , am(kx,nsplit), an(nsplit)
     b             , dstor(ix,jx,nsplit), hstor(ix,jx,nsplit)
#endif
cleexx
c
      integer lcount,iin,iout,numpts
      common /ldata/ lcount, iin, iout, numpts
      integer ilake, jlake, depth, freeze
      real*8  hi, hii, hs, eta, tlake(400)
      integer j,n

c----------------------------------------------------------------------
c-----output large-domain variables:
c
      write(iutl) mdate0
      write(iutl) ktau,xtime,ldatez,lyear,lmonth,lday,lhour,ntime
#ifdef MPP1
      if(EHSO4) then
         write(iutl) UB0_IO,VB0_IO,QB0_IO,TB0_IO,PS0_IO,TS0_IO,SO0_IO
      else
         write(iutl) UB0_IO,VB0_IO,QB0_IO,TB0_IO,PS0_IO,TS0_IO
      endif
      write(iutl) ua_io
      write(iutl) ub_io
      write(iutl) va_io
      write(iutl) vb_io
      write(iutl) ta_io
      write(iutl) tb_io
      write(iutl) qva_io
      write(iutl) qvb_io
      write(iutl) qca_io
      write(iutl) qcb_io
      write(iutl) psa_io,psb_io,satbrt_io,satbrt1_io,f_io
      write(iutl) ht_io,ht1_io,msfx_io,msfd_io,xlat_io,xlong_io
      write(iutl) tga_io,tgb_io,rainc_io,rainnc_io
      if(icup.eq.1) then
         write(iutl) rsheat_io,  rswat_io
      elseif(icup.eq.3) then
         write(iutl) TBASE_io,CLDEFI_io
      elseif(icup.eq.4) then
         write(iutl) cbmf2d_io
      endif
      write(iutl) hfx_io,qfx_io,snowc_io,uvdrag_io
#ifdef DIAG
      write(iutl) tdini,tdadv,tqini,tqadv,tqeva,tqrai
#endif
      write(iutl) absnxt_io, abstot_io, emstot_io
      if(ipptls.eq.1) write(iutl) fcc_io,rh0_io
      write(iutl) sol2d_io, solvd2d_io, solvs2d_io
     &          ,  flw2d_io, flwd2d_io, fsw2d_io
     &          , sabv2d_io, sinc2d_io
      write(iutl) taf2d_io, tlef2d_io, tgbb_io
     &          ,  ssw2d_io,  srw2d_io
     &          ,  tg2d_io,   tgb2d_io,  swt2d_io
     &          ,  scv2d_io, gwet2d_io,  veg2d_io, veg2d1_io
     &          ,  sag2d_io, sice2d_io,  dew2d_io, ircp2d_io
     &          , text2d_io,  col2d_io, ocld2d_io
     &          , heatrt_io, o3prof_io
      write(iutl) pptnc_io, pptc_io, prca2d_io, prnca2d_io 
      if (iocnflx.eq.2) write(iutl) zpbl_io
      if(ichem.eq.1) then
        write(iutl) chia_io
        write(iutl) chib_io
        write(iutl) remlsc_io
        write(iutl) remcvc_io
        write(iutl) remdrd_io
        write(iutl) ssw2da_io
        write(iutl) sdeltk2d_io
        write(iutl) sdelqk2d_io
        write(iutl) sfracv2d_io
        write(iutl) sfracb2d_io
        write(iutl) sfracs2d_io
        write(iutl) svegfrac2d_io
#ifdef DIAG
        write(iutl) tchiad
        write(iutl) tchitb
        write(iutl) tchie
#endif
      endif
#else
      if(EHSO4) then
         write(iutl) UB0,VB0,QB0,TB0,PS0,TS0,SO0
      else
         write(iutl) UB0,VB0,QB0,TB0,PS0,TS0
      endif
      write(iutl) ua
      write(iutl) ub
      write(iutl) va
      write(iutl) vb
      write(iutl) ta
      write(iutl) tb
      write(iutl) qva
      write(iutl) qvb
      write(iutl) qca
      write(iutl) qcb
      write(iutl) psa,psb,satbrt,satbrt1,f
      write(iutl) ht,ht1,msfx,msfd,xlat,xlong
      write(iutl) tga,tgb,rainc,rainnc
      if(icup.eq.1) then
         write(iutl) rsheat,  rswat
      elseif(icup.eq.3) then
         write(iutl) TBASE,CLDEFI
      elseif(icup.eq.4) then
         write(iutl) cbmf2d
      endif
      write(iutl) hfx,qfx,snowc,uvdrag
#ifdef DIAG
      write(iutl) tdini,tdadv,tqini,tqadv,tqeva,tqrai
#endif
      write(iutl) absnxt, abstot, emstot
      if(ipptls.eq.1) write(iutl)  fcc,rh0
      write(iutl) sol2d, solvd2d, solvs2d
     &           ,  flw2d, flwd2d, fsw2d
     &           , sabv2d, sinc2d
      write(iutl)  taf2d, tlef2d, tgbb
     b           ,  ssw2d,  srw2d
     c           ,  tg2d,   tgb2d,  swt2d
     d           ,  scv2d, gwet2d,  veg2d, veg2d1
     e           ,  sag2d, sice2d,  dew2d, ircp2d
     g           , text2d,  col2d, ocld2d
     h           , heatrt, o3prof
      write(iutl)  pptnc, pptc, prca2d, prnca2d
      if (iocnflx.eq.2) write(iutl) zpbl

Cchem2
      if(ichem.eq.1) then
c pronostic concentrartions
        write(iutl) chia
        write(iutl) chib
c cumul removal terms (3d, 2d)
        write(iutl) remlsc
        write(iutl) remcvc
        write(iutl) remdrd
        write(iutl) ssw2da
        write(iutl) sdeltk2d
        write(iutl) sdelqk2d
        write(iutl) sfracv2d
        write(iutl) sfracb2d
        write(iutl) sfracs2d
        write(iutl) svegfrac2d
c cumul ad, dif, emis terms ( scalar)
#ifdef DIAG
        write(iutl) tchiad
        write(iutl) tchitb
        write(iutl) tchie
#endif
      endif
Cchem2_
#endif

c---lake model
c
      if(lakemod.eq.1)then
       write(iutl) numpts
       print *, 'writing lake model restart file. numpts = ', numpts
       do n=1,numpts
        read(iin) ilake, jlake, depth, freeze, hi ,hii, hs
     a          , eta, (tlake(j),j=1,depth)
c       print *,' in outsav ilake jlake = ',ilake,jlake

        write(iutl) ilake, jlake, depth, freeze, hi, hii, hs
     a             , eta, (tlake(j),j=1,depth)
       enddo
       rewind(iin)
      endif
c
#ifdef MPP1
      write(iutl) dstor_io
      write(iutl) hstor_io
      write(iutl) uj1,uj2,ujlx,ujl
      write(iutl) ui1_io,ui2_io,uilx_io,uil_io
      write(iutl) vj1,vj2,vjlx,vjl
      write(iutl) vi1_io,vi2_io,vilx_io,vil_io
#else
      write(iutl) dstor
      write(iutl) hstor
      write(iutl) uj1,uj2,ujlx,ujl
      write(iutl) ui1,ui2,uilx,uil
      write(iutl) vj1,vj2,vjlx,vjl
      write(iutl) vi1,vi2,vilx,vil
#endif
      return
      end
