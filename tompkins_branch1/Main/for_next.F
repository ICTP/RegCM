      subroutine for_next
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
      include 'Commons/date.cb'
      include 'Commons/cupeman.cb'
#ifdef MPP1
      include 'Common2/pmoist.cb'
      include 'Common2/trachem.cb'
#else
      include 'Commons/pmoist.cb'
      include 'Commons/trachem.cb'
#endif
      integer itr,nc

      character*80 chname
      integer isystm,system
      external system

      isystm=system('/bin/mv -f regcm.in regcm0.in')
      OPEN(99,FILE='regcm.in',FORM='FORMATTED',STATUS='unknown')

      write(99,201) '$restartparam'
      ifrest = .true.
      if (ifrest) then
         write(99,101) 'ifrest  = .true. '
      else
         write(99,101) 'ifrest  = .false.'
      endif
      write(99,102) 'idate0  = ',idate0
      write(99,102) 'idate1  = ',idate1
      write(99,102) 'idate2  = ',idate2
      write(99,102) 'nslice  = ',nslice
      write(99,200) '$'

      write(99,202) '$timeparam'
      write(99,103) 'radfrq  = ',radfrq
      write(99,103) 'abemh   = ',abemh
      write(99,103) 'abatm   = ',abatm
      write(99,103) 'dt      = ',dtmin*60.
      write(99,102) 'ibdyfrq = ',ibdyfrq
      write(99,200) '$'

      write(99,203) '$outparam'
      if (ifsave) then
        write(99,101) 'ifsave  = .true. '
      else
        write(99,101) 'ifsave  = .false.'
      endif
      write(99,104) 'savfrq  = ',savfrq
      if (iftape) then
        write(99,101) 'iftape  = .true. '
      else
        write(99,101) 'iftape  = .false.'
      endif
      write(99,104) 'tapfrq  = ',tapfrq
      if (ifrad) then
        write(99,101) 'ifrad   = .true. '
      else
        write(99,101) 'ifrad   = .false.'
      endif
      write(99,104) 'radisp  = ',radisp
      if (ifbat) then
        write(99,101) 'ifbat   = .true. '
      else
        write(99,101) 'ifbat   = .false.'
      endif
      if (ifsub) then
        write(99,101) 'ifsub   = .true. '
      else
        write(99,101) 'ifsub   = .false.'
      endif
      write(99,104) 'batfrq  = ',batfrq
      if (ifprt) then
        write(99,101) 'ifprt   = .true. '
      else
        write(99,101) 'ifprt   = .false.'
      endif
      write(99,104) 'prtfrq  = ',prtfrq
      write(99,105) 'kxout   = ',kxout
      write(99,105) 'jxsex   = ',jxsex
      write(99,102) 'iotyp   = ',iotyp
      write(99,102) 'ibintyp = ',ibintyp
      if (ifchem) then
        write(99,101) 'ifchem  = .true. '
      else
        write(99,101) 'ifchem  = .false.'
      endif
      write(99,104) 'chemfrq = ',chemfrq
      write(99,200) '$'

      write(99,201) '$physicsparam'
      write(99,102) 'iboudy  = ',iboudy
      write(99,102) 'ibltyp  = ',ibltyp
      write(99,102) 'icup    = ',icup
      write(99,105) 'igcc    = ',igcc
      write(99,102) 'ipptls  = ',ipptls
      write(99,102) 'iocnflx = ',iocnflx
      write(99,102) 'ipgf    = ',ipgf
      write(99,102) 'iemiss  = ',iemiss
      write(99,102) 'lakemod = ',lakemod
      write(99,102) 'ichem   = ',ichem
      write(99,200) '$'

      write(99,204) '$subexparam'
      if(ncld.ne.1)          write(99,102) 'ncld    = ',ncld
      if(fcmax.ne.0.80)      write(99,106) 'fcmax   = ',fcmax
      if(qck1land.ne.0.0005) write(99,107) 'qck1land= ',qck1land
      if(qck1oce.ne.0.0005)  write(99,107) 'qck1oce = ',qck1oce
      if(gulland.ne.0.4)     write(99,106) 'gulland = ',gulland
      if(guloce.ne.0.4)      write(99,106) 'guloce  = ',guloce
      if(rhmax.ne.1.01)      write(99,106) 'rhmax   = ',rhmax
      if(rh0oce.ne.0.90)     write(99,106) 'rh0oce  = ',rh0oce
      if(rh0land.ne.0.80)    write(99,106) 'rh0land = ',rh0land
      if(tc0.ne.238.0)       write(99,106) 'tc0     = ',tc0
      if(cevap.ne.0.00002)   write(99,107) 'cevap   = ',cevap
      if(caccr.ne.6.0)       write(99,106) 'caccr   = ',caccr
      write(99,200) '$'

      write(99,204) '$grellparam'
      if(shrmin.ne.0.25)  write(99,106) 'shrmin  = ',shrmin
      if(shrmax.ne.0.50)  write(99,106) 'shrmax  = ',shrmax
      if(edtmin.ne.0.25)  write(99,106) 'edtmin  = ',edtmin
      if(edtmax.ne.1.00)  write(99,106) 'edtmax  = ',edtmax
      if(edtmino.ne.0.00) write(99,106) 'edtmino = ',edtmino
      if(edtmaxo.ne.1.00) write(99,106) 'edtmaxo = ',edtmaxo
      if(edtminx.ne.0.25) write(99,106) 'edtminx = ',edtminx
      if(edtmaxx.ne.1.00) write(99,106) 'edtmaxx = ',edtmaxx
      if(pbcmax.ne.150.)  write(99,106) 'pbcmax  = ',pbcmax
      if(mincld.ne.150.)  write(99,106) 'mincld  = ',mincld
      if(htmin.ne.-250.)  write(99,106) 'htmin   = ',htmin
      if(htmax.ne.500.)   write(99,106) 'htmax   = ',htmax
      if(skbmax.ne.0.4)   write(99,106) 'skbmax  = ',skbmax
      if(dtauc.ne.30.)    write(99,106) 'dtauc   = ',dtauc
      write(99,200) '$'

      write(99,202) '$emanparam'
      if(minsig.ne.0.95)  write(99,106) 'minsig  = ',minsig
      if(elcrit.ne.0.0011)write(99,106) 'elcrit  = ',elcrit
      if(tlcrit.ne.-55.0) write(99,106) 'tlcrit  = ',tlcrit
      if(entp.ne.1.5)     write(99,106) 'entp    = ',entp
      if(sigd.ne.0.05)    write(99,106) 'sigd    = ',sigd
      if(sigs.ne.0.12)    write(99,106) 'sigs    = ',sigs
      if(omtrain.ne.50.0) write(99,106) 'omtrain = ',omtrain
      if(omtsnow.ne.5.5)  write(99,106) 'omtsnow = ',omtsnow
      if(coeffr.ne.1.0)   write(99,106) 'coeffr  = ',coeffr
      if(coeffs.ne.0.8)   write(99,106) 'coeffs  = ',coeffs
      if(cu.ne.0.7)       write(99,106) 'cu      = ',cu
      if(betae.ne.10.0)   write(99,106) 'betae   = ',betae
      if(dtmax.ne.0.9)    write(99,106) 'dtmax   = ',dtmax
      if(alphae.ne.0.2)   write(99,106) 'alphae  = ',alphae
      if(damp.ne.0.1)     write(99,106) 'damp    = ',damp
      write(99,200) '$'

      write(99,204) '$chemparam'
      if(ichremlsc.ne.1) write(99,108) 'ichremlsc = ',ichremlsc
      if(ichremcvc.ne.1) write(99,108) 'ichremcvc = ',ichremcvc
      if(ichdrdepo.ne.1) write(99,108) 'ichdrdepo = ',ichdrdepo
      if(ichcumtra.ne.1) write(99,108) 'ichcumtra = ',ichcumtra
      write(99,108) 'idirect   = ',idirect
      nc=0

      do itr=1,ntr
        if(chtrname(itr).eq.'SO2'.or.chtrname(itr).eq.'SO4')then
          write(chname(nc+1:nc+6),111) chtrname(itr)
          nc=nc+6
        else if(chtrname(itr).eq.'DUST') then
          write(chname(nc+1:nc+7),112) chtrname(itr)
          nc=nc+7
        else
          write(chname(nc+1:nc+8),113) chtrname(itr)
          nc=nc+8
        endif
      enddo
      write(99,109) 'chtrname  = ',chname(1:nc)
      write(99,114) 'chtrsol   = ',chtrsol
      write(99,116) 'chtrdpv   = ',chtrdpv
      write(99,114) 'dustbsiz  = ',dustbsiz
      write(99,200) '$'

      close(99)
 200  format(1x,A1)
 201  format(1x,A13)
 202  format(1x,A10)
 203  format(1x,A9)
 204  format(1x,A11)

 101  format(1x,A17,',')
 102  format(1x,A10,I10,',')
 103  format(1x,A10,F6.0,',')
 104  format(3x,A10,F6.0,',')
 105  format(3x,A10,I10,',')
 106  format(1x,A10,F10.3,',')
 107  format(1x,A10,E10.3,',')
 108  format(1x,A12,I4,',')
 109  format(1x,A12,A75)
 111  format("'",A3,"',")
 112  format("'",A4,"',")
 113  format("'",A5,"',")
 114  format(1x,A12,10(F5.2,','))
 116  format(1x,A12,20(F7.5,','))
      return
      end
