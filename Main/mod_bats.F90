!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_bats
!
! Storage for Surface (BATS and shared by CLM) variables
!
      use mod_memutil
      use mod_runparams
      use mod_bats_param
!
      real(8) , allocatable , dimension(:,:) :: p1d0 , qs1d0 , ts1d0
!
      integer , allocatable , dimension(:,:) :: ldoc1d
      real(8) , allocatable , dimension(:,:) :: delq1d , delt1d ,       &
             drag1d , emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d ,    &
             p1d , pbp1d , prcp1d , q2m1d , qg1d , qs1d , resp1d ,      &
             rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
             sice1d , ssw1d , t2m1d , taf1d , tg1d , tgb1d , tlef1d ,   &
             ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d , aldirs1d ,  &
             aldifs1d
!
      real(8) , allocatable , dimension(:,:) :: bfc , bsw , evmx0 ,     &
             fdry , fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc ,      &
             rnet , texrat , vegt , wiltr , wt , xkmx
!
      real(8) , allocatable , dimension(:,:) :: cdr , cdrn , cdrx ,     &
           & cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg , &
           & etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,    &
           & fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,     &
           & rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,    &
           & tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,&
           & xlsai , xrun , z1log , z2fra , z10fra , zlgocn , zlglnd ,  &
           & zlgsno , zlgveg , zlgdis
!
      real(8) , allocatable , dimension(:,:) :: cn1 , rgr , wta0 ,      &
             wtaq0 , wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq ,  &
             wtgq , wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi ,        &
             wtsqi , df
!
      integer , allocatable , dimension(:,:) :: lveg , oveg
!
      real(8) :: difrat , rdnnsg
      real(4) :: rrnnsg
!
      integer :: ilat , ihis , mhis , ncase
!
      real(8) , allocatable , dimension(:) :: flw1d , fsw1d
      real(8) , allocatable , dimension(:) :: czen , sola , vpdd
      real(8) , allocatable , dimension(:) :: ems , us1d , vs1d
      real(8) , allocatable , dimension(:) :: albdif , albdir , albvl , &
                               & albvs , albvsd , aldifl , aldifs ,     &
                               & aldirl , aldirs , emiss1d , fracd ,    &
                               & sabveg , solis , solvd , solvs , albvld
!
      real(8) , allocatable , dimension(:) :: coszrs
!
      real(8) , pointer , dimension(:,:) :: flw2d , flwa2d ,       &
             flwd2d , flwda2d , fsw2d , fswa2d , pptc , pptnc ,    &
             prca2d , prnca2d , sabv2d , sina2d , sinc2d , sol2d , &
             solvd2d , solvs2d , svga2d
!
      real(8) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
            sdelqk2d , sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
      integer , allocatable , dimension(:,:,:) :: ocld2d , veg2d1
      integer , allocatable , dimension(:,:) :: veg2d , ldmsk
!
      real(8) , pointer , dimension(:,:,:) :: col2d , dew2d ,        &
           emiss2d , evpa2d , gwet2d , ircp2d , rno2d , rnos2d ,     &
           sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d , &
           taf2d , tg2d , tgb2d , tlef2d
!
      real(8) ,allocatable, dimension(:,:,:) :: ht1 , satbrt1 , xlat1 , &
                                             &  xlon1
!
      real(4) , target , allocatable, dimension(:,:,:) :: fbat
!
      real(4) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
            flwd_o , fswa_o , prcv_o , psmn_o , ps_o , q2m_o , rnos_o , &
            rsw_o , scv_o , sena_o , sina_o , ssw_o , t2mn_o , t2mx_o , &
            t2m_o , tgmn_o , tgmx_o , tg_o , tlef_o , tpr_o , u10m_o ,  &
            v10m_o , w10x_o , zpbl_o , aldirs_o , aldifs_o
!
      real(4) , target ,allocatable, dimension(:,:,:,:) :: fsub
!
      real(4) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
             ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
             t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
!
      ! dtskin is difference between skin temp and bulk sst
      real(8) , allocatable , dimension(:,:) :: deltas , tdeltas ,      &
                           &                    dtskin
      logical , allocatable , dimension(:,:) :: firstcall

      contains

        subroutine allocate_mod_bats
        implicit none

        rrnnsg = 1.0/real(nnsg)
        rdnnsg = d_one/dble(nnsg)

        allocate(veg2d(iym1,jxp))
        veg2d = -1
        allocate(ldmsk(iym1,jxp))
        ldmsk = -1

        call getmem2d(flw2d,iym1,jxp,'bats:flw2d')
        call getmem2d(flwa2d,iym1,jxp,'bats:flwa2d')
        call getmem2d(flwd2d,iym1,jxp,'bats:flwd2d')
        call getmem2d(flwda2d,iym1,jxp,'bats:flwda2d')
        call getmem2d(fsw2d,iym1,jxp,'bats:fsw2d')
        call getmem2d(fswa2d,iym1,jxp,'bats:fswa2d')
        call getmem2d(pptc,iym1,jxp,'bats:pptc')
        call getmem2d(pptnc,iym1,jxp,'bats:pptnc')
        call getmem2d(prca2d,iym1,jxp,'bats:prca2d')
        call getmem2d(prnca2d,iym1,jxp,'bats:prnca2d')
        call getmem2d(sabv2d,iym1,jxp,'bats:sabv2d')
        call getmem2d(sina2d,iym1,jxp,'bats:sina2d')
        call getmem2d(sinc2d,iym1,jxp,'bats:sinc2d')
        call getmem2d(sol2d,iym1,jxp,'bats:sol2d')
        call getmem2d(solvd2d,iym1,jxp,'bats:solvd2d')
        call getmem2d(solvs2d,iym1,jxp,'bats:solvs2d')
        call getmem2d(svga2d,iym1,jxp,'bats:svga2d')
        if ( ichem == 1 ) then
          call getmem2d(ssw2da,iym1,jxp,'bats:ssw2da')
          call getmem2d(sdeltk2d,iym1,jxp,'bats:sdeltk2d')
          call getmem2d(sdelqk2d,iym1,jxp,'bats:sdelqk2d')
          call getmem2d(sfracv2d,iym1,jxp,'bats:sfracv2d')
          call getmem2d(sfracb2d,iym1,jxp,'bats:sfracb2d')
          call getmem2d(sfracs2d,iym1,jxp,'bats:sfracs2d')
          call getmem2d(svegfrac2d,iym1,jxp,'bats:svegfrac2d')
        end if

        allocate(ocld2d(nnsg,iym1,jxp))
        ocld2d = -1
        allocate(veg2d1(nnsg,iym1,jxp))
        veg2d1 = -1

        call getmem3d(col2d,nnsg,iym1,jxp,'bats:col2d')
        call getmem3d(dew2d,nnsg,iym1,jxp,'bats:dew2d')
        call getmem3d(emiss2d,nnsg,iym1,jxp,'bats:emiss2d')
        call getmem3d(evpa2d,nnsg,iym1,jxp,'bats:evpa2d')
        call getmem3d(gwet2d,nnsg,iym1,jxp,'bats:gwet2d')
        call getmem3d(ircp2d,nnsg,iym1,jxp,'bats:ircp2d')
        call getmem3d(rno2d,nnsg,iym1,jxp,'bats:rno2d')
        call getmem3d(rnos2d,nnsg,iym1,jxp,'bats:rnos2d')
        call getmem3d(sag2d,nnsg,iym1,jxp,'bats:sag2d')
        call getmem3d(scv2d,nnsg,iym1,jxp,'bats:scv2d')
        call getmem3d(sena2d,nnsg,iym1,jxp,'bats:sena2d')
        call getmem3d(sice2d,nnsg,iym1,jxp,'bats:sice2d')
        call getmem3d(srw2d,nnsg,iym1,jxp,'bats:srw2d')
        call getmem3d(ssw2d,nnsg,iym1,jxp,'bats:ssw2d')
        call getmem3d(swt2d,nnsg,iym1,jxp,'bats:swt2d')
        call getmem3d(taf2d,nnsg,iym1,jxp,'bats:taf2d')
        call getmem3d(tg2d,nnsg,iym1,jxp,'bats:tg2d')
        call getmem3d(tgb2d,nnsg,iym1,jxp,'bats:tgb2d')
        call getmem3d(tlef2d,nnsg,iym1,jxp,'bats:tlef2d')

        allocate(ht1(nnsg,iy,jxp))
        allocate(satbrt1(nnsg,iy,jxp))
        allocate(xlat1(nnsg,iy,jxp))
        allocate(xlon1(nnsg,iy,jxp))
        ht1 = d_zero
        satbrt1 = d_zero
        xlat1 = d_zero
        xlon1 = d_zero

        if (idcsst == 1) then
          allocate(deltas(iy,jxp))
          allocate(tdeltas(iy,jxp))
          allocate(dtskin(iy,jxp))
          allocate(firstcall(iy,jxp))
          deltas = d_zero
          tdeltas = d_zero
          dtskin = d_zero
          firstcall = .false.
        end if

        allocate(ldoc1d(nnsg,iym1))
        ldoc1d = -1

        allocate(p1d0(nnsg,iym1))
        p1d0 = d_zero
        allocate(qs1d0(nnsg,iym1))
        qs1d0 = d_zero
        allocate(ts1d0(nnsg,iym1))
        ts1d0 = d_zero
        allocate(delq1d(nnsg,iym1))
        delq1d = d_zero
        allocate(delt1d(nnsg,iym1))
        delt1d = d_zero
        allocate(drag1d(nnsg,iym1))
        drag1d = d_zero
        allocate(emiss_1d(nnsg,iym1))
        emiss_1d = d_zero
        allocate(evpr1d(nnsg,iym1))
        evpr1d = d_zero
        allocate(gwet1d(nnsg,iym1))
        gwet1d = d_zero
        allocate(ircp1d(nnsg,iym1))
        ircp1d = d_zero
        allocate(ldew1d(nnsg,iym1))
        ldew1d = d_zero
        allocate(p1d(nnsg,iym1))
        p1d = d_zero
        allocate(pbp1d(nnsg,iym1))
        pbp1d = d_zero
        allocate(prcp1d(nnsg,iym1))
        prcp1d = d_zero
        allocate(q2m1d(nnsg,iym1))
        q2m1d = d_zero
        allocate(qg1d(nnsg,iym1))
        qg1d = d_zero
        allocate(qs1d(nnsg,iym1))
        qs1d = d_zero
        allocate(resp1d(nnsg,iym1))
        resp1d = d_zero
        allocate(rhs1d(nnsg,iym1))
        rhs1d = d_zero
        allocate(rno1d(nnsg,iym1))
        rno1d = d_zero
        allocate(rnos1d(nnsg,iym1))
        rnos1d = d_zero
        allocate(rsw1d(nnsg,iym1))
        rsw1d = d_zero
        allocate(sag1d(nnsg,iym1))
        sag1d = d_zero
        allocate(scv1d(nnsg,iym1))
        scv1d = d_zero
        allocate(sent1d(nnsg,iym1))
        sent1d = d_zero
        allocate(sice1d(nnsg,iym1))
        sice1d = d_zero
        allocate(ssw1d(nnsg,iym1))
        ssw1d = d_zero
        allocate(t2m1d(nnsg,iym1))
        t2m1d = d_zero
        allocate(taf1d(nnsg,iym1))
        taf1d = d_zero
        allocate(tg1d(nnsg,iym1))
        tg1d = d_zero
        allocate(tgb1d(nnsg,iym1))
        tgb1d = d_zero
        allocate(tlef1d(nnsg,iym1))
        tlef1d = d_zero
        allocate(ts1d(nnsg,iym1))
        ts1d = d_zero
        allocate(tsw1d(nnsg,iym1))
        tsw1d = d_zero
        allocate(u10m1d(nnsg,iym1))
        u10m1d = d_zero
        allocate(v10m1d(nnsg,iym1))
        v10m1d = d_zero
        allocate(veg1d(nnsg,iym1))
        veg1d = d_zero
        allocate(z1d(nnsg,iym1))
        z1d = d_zero
        allocate(bfc(nnsg,iym1))
        bfc = d_zero
        allocate(bsw(nnsg,iym1))
        bsw = d_zero
        allocate(evmx0(nnsg,iym1))
        evmx0 = d_zero
        allocate(fdry(nnsg,iym1))
        fdry = d_zero
        allocate(fwet(nnsg,iym1))
        fwet = d_zero
        allocate(gwmx0(nnsg,iym1))
        gwmx0 = d_zero
        allocate(gwmx1(nnsg,iym1))
        gwmx1 = d_zero
        allocate(gwmx2(nnsg,iym1))
        gwmx2 = d_zero
        allocate(porsl(nnsg,iym1))
        porsl = d_zero
        allocate(relfc(nnsg,iym1))
        relfc = d_zero
        allocate(rnet(nnsg,iym1))
        rnet = d_zero
        allocate(texrat(nnsg,iym1))
        texrat = d_zero
        allocate(vegt(nnsg,iym1))
        vegt = d_zero
        allocate(wiltr(nnsg,iym1))
        wiltr = d_zero
        allocate(wt(nnsg,iym1))
        wt = d_zero
        allocate(xkmx(nnsg,iym1))
        xkmx = d_zero
        allocate(cdr(nnsg,iym1))
        cdr = d_zero
        allocate(cdrn(nnsg,iym1))
        cdrn = d_zero
        allocate(cdrx(nnsg,iym1))
        cdrx = d_zero
        allocate(cf(nnsg,iym1))
        cf = d_zero
        allocate(cgrnd(nnsg,iym1))
        cgrnd = d_zero
        allocate(cgrndl(nnsg,iym1))
        cgrndl = d_zero
        allocate(cgrnds(nnsg,iym1))
        cgrnds = d_zero
        allocate(clead(nnsg,iym1))
        clead = d_zero
        allocate(densi(nnsg,iym1))
        densi = d_zero
        allocate(efpr(nnsg,iym1))
        efpr = d_zero
        allocate(eg(nnsg,iym1))
        eg = d_zero
        allocate(etr(nnsg,iym1))
        etr = d_zero
        allocate(etrrun(nnsg,iym1))
        etrrun = d_zero
        allocate(evaps(nnsg,iym1))
        evaps = d_zero
        allocate(evapw(nnsg,iym1))
        evapw = d_zero
        allocate(fevpg(nnsg,iym1))
        fevpg = d_zero
        allocate(flnet(nnsg,iym1))
        flnet = d_zero
        allocate(flneto(nnsg,iym1))
        flneto = d_zero
        allocate(fseng(nnsg,iym1))
        fseng = d_zero
        allocate(htvp(nnsg,iym1))
        htvp = d_zero
        allocate(ps(nnsg,iym1))
        ps = d_zero
        allocate(pw(nnsg,iym1))
        pw = d_zero
        allocate(qice(nnsg,iym1))
        qice = d_zero
        allocate(qsatl(nnsg,iym1))
        qsatl = d_zero
        allocate(rhosw(nnsg,iym1))
        rhosw = d_zero
        allocate(ribd(nnsg,iym1))
        ribd = d_zero
        allocate(rlai(nnsg,iym1))
        rlai = d_zero
        allocate(rpp(nnsg,iym1))
        rpp = d_zero
        allocate(scrat(nnsg,iym1))
        scrat = d_zero
        allocate(scvk(nnsg,iym1))
        scvk = d_zero
        allocate(sdrop(nnsg,iym1))
        sdrop = d_zero
        allocate(seasb(nnsg,iym1))
        seasb = d_zero
        allocate(sigf(nnsg,iym1))
        sigf = d_zero
        allocate(sm(nnsg,iym1))
        sm = d_zero
        allocate(tm(nnsg,iym1))
        tm = d_zero
        allocate(uaf(nnsg,iym1))
        uaf = d_zero
        allocate(vspda(nnsg,iym1))
        vspda = d_zero
        allocate(wata(nnsg,iym1))
        wata = d_zero
        allocate(watr(nnsg,iym1))
        watr = d_zero
        allocate(watt(nnsg,iym1))
        watt = d_zero
        allocate(watu(nnsg,iym1))
        watu = d_zero
        allocate(wta(nnsg,iym1))
        wta = d_zero
        allocate(xlai(nnsg,iym1))
        xlai = d_zero
        allocate(xlsai(nnsg,iym1))
        xlsai = d_zero
        allocate(xrun(nnsg,iym1))
        xrun = d_zero
        allocate(cn1(nnsg,iym1))
        cn1 = d_zero
        allocate(df(nnsg,iym1))
        df = d_zero
        allocate(rgr(nnsg,iym1))
        rgr = d_zero
        allocate(wta0(nnsg,iym1))
        wta0 = d_zero
        allocate(wtaq0(nnsg,iym1))
        wtaq0 = d_zero
        allocate(wtg(nnsg,iym1))
        wtg = d_zero
        allocate(wtg0(nnsg,iym1))
        wtg0 = d_zero
        allocate(wtg2(nnsg,iym1))
        wtg2 = d_zero
        allocate(wtga(nnsg,iym1))
        wtga = d_zero
        allocate(wtgaq(nnsg,iym1))
        wtgaq = d_zero
        allocate(wtgl(nnsg,iym1))
        wtgl = d_zero
        allocate(wtglq(nnsg,iym1))
        wtglq = d_zero
        allocate(wtgq(nnsg,iym1))
        wtgq = d_zero
        allocate(wtgq0(nnsg,iym1))
        wtgq0 = d_zero
        allocate(wtl0(nnsg,iym1))
        wtl0 = d_zero
        allocate(wtlh(nnsg,iym1))
        wtlh = d_zero
        allocate(wtlq(nnsg,iym1))
        wtlq = d_zero
        allocate(wtlq0(nnsg,iym1))
        wtlq0 = d_zero
        allocate(wtshi(nnsg,iym1))
        wtshi = d_zero
        allocate(wtsqi(nnsg,iym1))
        wtsqi = d_zero
        allocate(z2fra(nnsg,iym1))
        z2fra = d_zero
        allocate(z10fra(nnsg,iym1))
        z10fra = d_zero
        allocate(z1log(nnsg,iym1))
        z1log = d_zero
        allocate(zlgocn(nnsg,iym1))
        zlgocn = d_zero
        allocate(zlglnd(nnsg,iym1))
        zlglnd = d_zero
        allocate(zlgsno(nnsg,iym1))
        zlgsno = d_zero
        allocate(zlgveg(nnsg,iym1))
        zlgveg = d_zero
        allocate(zlgdis(nnsg,iym1))
        zlgdis = d_zero
        allocate(aldirs1d(nnsg,iym1))
        aldirs1d = d_zero
        allocate(aldifs1d(nnsg,iym1))
        aldifs1d = d_zero
        allocate(lveg(nnsg,iym1))
        lveg = -1
        allocate(oveg(nnsg,iym1))
        oveg = -1
        allocate(coszrs(iy))
        coszrs = d_zero
        allocate(flw1d(iym1))
        flw1d = d_zero
        allocate(fsw1d(iym1))
        fsw1d = d_zero
        allocate(us1d(iym1))
        us1d = d_zero
        allocate(vs1d(iym1))
        vs1d = d_zero
        allocate(czen(iym1))
        czen = d_zero
        allocate(sola(iym1))
        sola = d_zero
        allocate(vpdd(iym1))
        vpdd = d_zero
        allocate(ems(iym1))
        ems = d_zero
        allocate(albdif(iym1))
        albdif = d_zero
        allocate(albdir(iym1))
        albdir = d_zero
        allocate(albvl(iym1))
        albvl = d_zero
        allocate(albvld(iym1))
        albvld = d_zero
        allocate(albvs(iym1))
        albvs = d_zero
        allocate(albvsd(iym1))
        albvsd = d_zero
        allocate(aldifl(iym1))
        aldifl = d_zero
        allocate(aldifs(iym1))
        aldifs = d_zero
        allocate(aldirl(iym1))
        aldirl = d_zero
        allocate(aldirs(iym1))
        aldirs = d_zero
        allocate(emiss1d(iym1))
        emiss1d = d_zero
        allocate(fracd(iym1))
        fracd = d_zero
        allocate(sabveg(iym1))
        sabveg = d_zero
        allocate(solis(iym1))
        solis = d_zero
        allocate(solvd(iym1))
        solvd = d_zero
        allocate(solvs(iym1))
        solvs = d_zero

        allocate(fbat(jxp,iym2,numbat))
        fbat = 0.0
        ps_o   => fbat(:,:,1)
        u10m_o => fbat(:,:,2)
        v10m_o => fbat(:,:,3)
        drag_o => fbat(:,:,4)
        tg_o   => fbat(:,:,5)
        tlef_o => fbat(:,:,6)
        t2m_o  => fbat(:,:,7)
        q2m_o  => fbat(:,:,8)
        ssw_o  => fbat(:,:,9)
        rsw_o  => fbat(:,:,10)
        tpr_o  => fbat(:,:,11)
        evpa_o => fbat(:,:,12)
        rnos_o => fbat(:,:,13)
        scv_o  => fbat(:,:,14)
        sena_o => fbat(:,:,15)
        flwa_o => fbat(:,:,16)
        fswa_o => fbat(:,:,17)
        flwd_o => fbat(:,:,18)
        sina_o => fbat(:,:,19)
        prcv_o => fbat(:,:,20)
        zpbl_o => fbat(:,:,21)
        tgmx_o => fbat(:,:,22)
        tgmn_o => fbat(:,:,23)
        t2mx_o => fbat(:,:,24)
        t2mn_o => fbat(:,:,25)
        w10x_o => fbat(:,:,26)
        psmn_o => fbat(:,:,27)
        aldirs_o => fbat(:,:,28)
        aldifs_o => fbat(:,:,29)

        allocate(fsub(nnsg,jxp,iym2,numsub))
        fsub = 0.0
        ps_s   => fsub(:,:,:,1)
        u10m_s => fsub(:,:,:,2)
        v10m_s => fsub(:,:,:,3)
        drag_s => fsub(:,:,:,4)
        tg_s   => fsub(:,:,:,5)
        tlef_s => fsub(:,:,:,6)
        t2m_s  => fsub(:,:,:,7)
        q2m_s  => fsub(:,:,:,8)
        ssw_s  => fsub(:,:,:,9)
        rsw_s  => fsub(:,:,:,10)
        tpr_s  => fsub(:,:,:,11)
        evpa_s => fsub(:,:,:,12)
        rnos_s => fsub(:,:,:,13)
        scv_s  => fsub(:,:,:,14)
        sena_s => fsub(:,:,:,15)
        prcv_s => fsub(:,:,:,16)
!
      end subroutine allocate_mod_bats 
!
      end module mod_bats
