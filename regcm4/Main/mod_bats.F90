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

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_bats_param
      use mod_main
      use mod_pbldim
      use mod_slice
      use mod_date

      implicit none
      real(8) , allocatable , target , dimension(:,:,:) :: spacebs1d
      real(8) , allocatable , target , dimension(:,:) :: spaceb1d
      private :: spacebs1d , spaceb1d
!
      real(8) , pointer , dimension(:,:) :: p1d0 , qs1d0 , ts1d0
!
      real(8) , pointer , dimension(:,:) :: delq1d , delt1d , drag1d ,  &
           & emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d , ldoc1d ,    &
           & p1d , pbp1d , prcp1d , q2m_1d , qg1d , qs1d , resp1d ,     &
           & rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
           & sice1d , ssw1d , t2m_1d , taf1d , tg1d , tgb1d , tlef1d ,  &
           & ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d
!
      real(8) , pointer , dimension(:,:) :: bfc , bsw , evmx0 , fdry ,  &
           & fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc , rnet ,      &
           & texrat , vegt , wiltr , wt , xkmx
!
      real(8) , pointer , dimension(:,:) :: aarea , cdr , cdrn , cdrx , &
           & cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg , &
           & etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,    &
           & fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,     &
           & rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,    &
           & tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,&
           & xlsai , xrun , z1 , z1log
!
      real(8) , pointer , dimension(:,:) :: cn1 , rgr , wta0 , wtaq0 ,  &
           & wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq , wtgq ,   &
           & wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi , df
!
      integer , allocatable , dimension(:,:) :: imelt , lveg
!
      real(8) :: difrat
!
      integer :: ilat , ihis , mhis , ncase
!
      real(8) , pointer , dimension(:) :: flw1d , fsw1d , us1d , vs1d
      real(8) , pointer , dimension(:) :: czen , sola , vpdd
      real(8) , pointer , dimension(:) :: ems
      real(8) , pointer , dimension(:) :: albdif , albdir , albvl ,     &
                               & albvs , albvsd , aldifl , aldifs ,     &
                               & aldirl , aldirs , emiss1d , fracd ,    &
                               & sabveg , solis , solvd , solvs , albvld
!
      real(8) , allocatable , dimension(:) :: coszrs
!
      real(8) , allocatable, dimension(:,:) :: flw2d , flwa2d , flwd2d ,&
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
!
      real(8) , allocatable, dimension(:,:,:) :: col2d , dew2d ,        &
           & emiss2d , evpa2d , gwet2d , ircp2d , ocld2d , rno2d ,      &
           & rnos2d , sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , &
           & swt2d , taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1 ,  &
           & lkdpth
      real(8) ,allocatable, dimension(:,:,:) :: ht1 , satbrt1 , xlat1 , &
                                             &  xlon1
!
      real(4) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
                                     & flwd_o , fswa_o , prcv_o ,       &
                                     & psmn_o , ps_o , q2m_o , rnos_o , &
                                     & rsw_o , scv_o , sena_o , sina_o ,&
                                     & ssw_o , t2mn_o , t2mx_o , t2m_o ,&
                                     & tgmn_o , tgmx_o , tg_o , tlef_o ,&
                                     & tpr_o , u10m_o , v10m_o ,        &
                                     & w10x_o , zpbl_o
!
      real(4) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
!
      real(4) , target , allocatable, dimension(:,:,:) :: fbat
      real(4) , target ,allocatable, dimension(:,:,:,:) :: fsub
!
#ifdef MPP1
#ifdef CLM
      real(8) , allocatable , target , dimension(:,:,:) :: spaceclm
      private :: spaceclm
      ! Direct solar rad incident on surface (<0.7)
      real(8) , pointer , dimension(:,:) :: sols2d
      ! Direct solar rad incident on surface (>=0.7)
      real(8) , pointer , dimension(:,:) :: soll2d
      ! Diffuse solar rad incident on surface (<0.7)
      real(8) , pointer , dimension(:,:) :: solsd2d
      ! Diffuse solar rad incident on surface (>=0.7)
      real(8) , pointer , dimension(:,:) :: solld2d
      real(8) , pointer , dimension(:,:) :: aldirs2d
      real(8) , pointer , dimension(:,:) :: aldirl2d
      real(8) , pointer , dimension(:,:) :: aldifs2d
      real(8) , pointer , dimension(:,:) :: aldifl2d
      real(8) , pointer , dimension(:,:) :: coszrs2d
      real(8) , pointer , dimension(:,:) :: rs2d
      real(8) , pointer , dimension(:,:) :: ra2d
      ! 2 meter specific humidity
      real(8) , pointer , dimension(:,:) :: q2d
#endif
#endif

#ifdef DCSST
      ! dtskin is difference between skin tem and bulk sst
      real(8) , allocatable , dimension(:,:) :: deltas , tdeltas ,      &
                           &                    dtskin
      logical , allocatable , dimension(:,:) :: firstcall
#endif

      contains

       subroutine allocate_mod_bats 
       implicit none
#ifdef MPP1
        allocate(flw2d(iym1,jxp)) 
        allocate(flwa2d(iym1,jxp))
        allocate(flwd2d(iym1,jxp))
        allocate(flwda2d(iym1,jxp))
        allocate(fsw2d(iym1,jxp))
        allocate(fswa2d(iym1,jxp))
        allocate(pptc(iym1,jxp))  
        allocate(pptnc(iym1,jxp))
        allocate(prca2d(iym1,jxp))
        allocate(prnca2d(iym1,jxp))
        allocate(sabv2d(iym1,jxp))
        allocate(sdelqk2d(iym1,jxp))
        allocate(sdeltk2d(iym1,jxp))
        allocate(sfracb2d (iym1,jxp))
        allocate(sfracs2d(iym1,jxp))
        allocate(sfracv2d(iym1,jxp))
        allocate(sina2d(iym1,jxp))
        allocate(sinc2d(iym1,jxp))
        allocate(sol2d(iym1,jxp))
        allocate(solvd2d(iym1,jxp))
        allocate(solvs2d(iym1,jxp))
        allocate(ssw2da(iym1,jxp))
        allocate(svegfrac2d(iym1,jxp))
        allocate(svga2d(iym1,jxp))
        allocate(veg2d(iym1,jxp))
        allocate(col2d(nnsg,iym1,jxp))
        allocate(dew2d(nnsg,iym1,jxp))
        allocate(emiss2d(nnsg,iym1,jxp))
        allocate(evpa2d(nnsg,iym1,jxp))
        allocate(gwet2d(nnsg,iym1,jxp))
        allocate(ircp2d(nnsg,iym1,jxp))
        allocate(ocld2d(nnsg,iym1,jxp))
        allocate(rno2d(nnsg,iym1,jxp))
        allocate(rnos2d(nnsg,iym1,jxp))
        allocate(sag2d(nnsg,iym1,jxp))
        allocate(scv2d(nnsg,iym1,jxp))
        allocate(sena2d(nnsg,iym1,jxp))
        allocate(sice2d(nnsg,iym1,jxp))
        allocate(srw2d(nnsg,iym1,jxp))
        allocate(ssw2d(nnsg,iym1,jxp))
        allocate(swt2d(nnsg,iym1,jxp))
        allocate(taf2d(nnsg,iym1,jxp))
        allocate(text2d(nnsg,iym1,jxp))
        allocate(tg2d(nnsg,iym1,jxp))
        allocate(tgb2d(nnsg,iym1,jxp))
        allocate(tlef2d(nnsg,iym1,jxp))
        allocate(veg2d1(nnsg,iym1,jxp))
        allocate(lkdpth(nnsg,iym1,jxp))
        allocate(ht1(nnsg,iy,jxp))
        allocate(satbrt1(nnsg,iy,jxp))
        allocate(xlat1(nnsg,iy,jxp))
        allocate(xlon1(nnsg,iy,jxp))
        allocate(fbat(jxp,iym2,numbat))
        allocate(fsub(nnsg,jxp,iym2,numsub))
#ifdef CLM
        allocate(spaceclm(iym1,jxp,12))
        sols2d   => spaceclm(:,:,1)
        soll2d   => spaceclm(:,:,2)
        solsd2d  => spaceclm(:,:,3)
        solld2d  => spaceclm(:,:,4)
        aldirs2d => spaceclm(:,:,5)
        aldirl2d => spaceclm(:,:,6)
        aldifs2d => spaceclm(:,:,7)
        aldifl2d => spaceclm(:,:,8)
        coszrs2d => spaceclm(:,:,9)
        rs2d     => spaceclm(:,:,10)
        ra2d     => spaceclm(:,:,11)
        q2d      => spaceclm(:,:,12)
#endif
#ifdef DCSST
        allocate(deltas(iy,jxp))
        allocate(tdeltas(iy,jxp))
        allocate(dtskin(iy,jxp))
        allocate(firstcall(iy,jxp))
#endif
#else
#ifdef BAND
        allocate(flw2d(iym1,jx)) 
        allocate(flwa2d(iym1,jx))
        allocate(flwd2d(iym1,jx))
        allocate(flwda2d(iym1,jx))
        allocate(fsw2d(iym1,jx))
        allocate(fswa2d(iym1,jx))
        allocate(pptc(iym1,jx))  
        allocate(pptnc(iym1,jx))
        allocate(prca2d(iym1,jx))
        allocate(prnca2d(iym1,jx))
        allocate(sabv2d(iym1,jx))
        allocate(sdelqk2d(iym1,jx))
        allocate(sdeltk2d(iym1,jx))
        allocate(sfracb2d (iym1,jx))
        allocate(sfracs2d(iym1,jx))
        allocate(sfracv2d(iym1,jx))
        allocate(sina2d(iym1,jx))
        allocate(sinc2d(iym1,jx))
        allocate(sol2d(iym1,jx))
        allocate(solvd2d(iym1,jx))
        allocate(solvs2d(iym1,jx))
        allocate(ssw2da(iym1,jx))
        allocate(svegfrac2d(iym1,jx))
        allocate(svga2d(iym1,jx))
        allocate(veg2d(iym1,jx))
        allocate(col2d(nnsg,iym1,jx))
        allocate(dew2d(nnsg,iym1,jx))
        allocate(emiss2d(nnsg,iym1,jx))
        allocate(evpa2d(nnsg,iym1,jx))
        allocate(gwet2d(nnsg,iym1,jx))
        allocate(ircp2d(nnsg,iym1,jx))
        allocate(ocld2d(nnsg,iym1,jx))
        allocate(rno2d(nnsg,iym1,jx))
        allocate(rnos2d(nnsg,iym1,jx))
        allocate(sag2d(nnsg,iym1,jx))
        allocate(scv2d(nnsg,iym1,jx))
        allocate(sena2d(nnsg,iym1,jx))
        allocate(sice2d(nnsg,iym1,jx))
        allocate(srw2d(nnsg,iym1,jx))
        allocate(ssw2d(nnsg,iym1,jx))
        allocate(swt2d(nnsg,iym1,jx))
        allocate(taf2d(nnsg,iym1,jx))
        allocate(text2d(nnsg,iym1,jx))
        allocate(tg2d(nnsg,iym1,jx))
        allocate(tgb2d(nnsg,iym1,jx))
        allocate(tlef2d(nnsg,iym1,jx))
        allocate(veg2d1(nnsg,iym1,jx))
        allocate(lkdpth(nnsg,iym1,jx))
#else
        allocate(flw2d(iym1,jxm1)) 
        allocate(flwa2d(iym1,jxm1))
        allocate(flwd2d(iym1,jxm1))
        allocate(flwda2d(iym1,jxm1))
        allocate(fsw2d(iym1,jxm1))
        allocate(fswa2d(iym1,jxm1))
        allocate(pptc(iym1,jxm1))  
        allocate(pptnc(iym1,jxm1))
        allocate(prca2d(iym1,jxm1))
        allocate(prnca2d(iym1,jxm1))
        allocate(sabv2d(iym1,jxm1))
        allocate(sdelqk2d(iym1,jxm1))
        allocate(sdeltk2d(iym1,jxm1))
        allocate(sfracb2d (iym1,jxm1))
        allocate(sfracs2d(iym1,jxm1))
        allocate(sfracv2d(iym1,jxm1))
        allocate(sina2d(iym1,jxm1))
        allocate(sinc2d(iym1,jxm1))
        allocate(sol2d(iym1,jxm1))
        allocate(solvd2d(iym1,jxm1))
        allocate(solvs2d(iym1,jxm1))
        allocate(ssw2da(iym1,jxm1))
        allocate(svegfrac2d(iym1,jxm1))
        allocate(svga2d(iym1,jxm1))
        allocate(veg2d(iym1,jxm1))
        allocate(col2d(nnsg,iym1,jxm1))
        allocate(dew2d(nnsg,iym1,jxm1))
        allocate(emiss2d(nnsg,iym1,jxm1))
        allocate(evpa2d(nnsg,iym1,jxm1))
        allocate(gwet2d(nnsg,iym1,jxm1))
        allocate(ircp2d(nnsg,iym1,jxm1))
        allocate(ocld2d(nnsg,iym1,jxm1))
        allocate(rno2d(nnsg,iym1,jxm1))
        allocate(rnos2d(nnsg,iym1,jxm1))
        allocate(sag2d(nnsg,iym1,jxm1))
        allocate(scv2d(nnsg,iym1,jxm1))
        allocate(sena2d(nnsg,iym1,jxm1))
        allocate(sice2d(nnsg,iym1,jxm1))
        allocate(srw2d(nnsg,iym1,jxm1))
        allocate(ssw2d(nnsg,iym1,jxm1))
        allocate(swt2d(nnsg,iym1,jxm1))
        allocate(taf2d(nnsg,iym1,jxm1))
        allocate(text2d(nnsg,iym1,jxm1))
        allocate(tg2d(nnsg,iym1,jxm1))
        allocate(tgb2d(nnsg,iym1,jxm1))
        allocate(tlef2d(nnsg,iym1,jxm1))
        allocate(veg2d1(nnsg,iym1,jxm1))
        allocate(lkdpth(nnsg,iym1,jxm1))
#endif
        allocate(ht1(nnsg,iy,jx))
        allocate(satbrt1(nnsg,iy,jx))
#ifdef BAND
        allocate(fbat(jx,iym2,numbat))
        allocate(fsub(nnsg,jx,iym2,numsub))
#else
        allocate(fbat(jxm2,iym2,numbat))
        allocate(fsub(nnsg,jxm2,iym2,numsub))
#endif
#ifdef DCSST
        allocate(deltas(iy,jx))
        allocate(tdeltas(iy,jx))
        allocate(dtskin(iy,jx))
        allocate(firstcall(iy,jx))
#endif
#endif
        allocate(spacebs1d(nnsg,iym1,123))
        p1d0     => spacebs1d(:,:,1)
        qs1d0    => spacebs1d(:,:,2)
        ts1d0    => spacebs1d(:,:,3)
        delq1d   => spacebs1d(:,:,4)
        delt1d   => spacebs1d(:,:,5)
        drag1d   => spacebs1d(:,:,6)
        emiss_1d => spacebs1d(:,:,7)
        evpr1d   => spacebs1d(:,:,8)
        gwet1d   => spacebs1d(:,:,9)
        ircp1d   => spacebs1d(:,:,10)
        ldew1d   => spacebs1d(:,:,11)
        ldoc1d   => spacebs1d(:,:,12)
        p1d      => spacebs1d(:,:,13)
        pbp1d    => spacebs1d(:,:,14)
        prcp1d   => spacebs1d(:,:,15)
        q2m_1d   => spacebs1d(:,:,16)
        qg1d     => spacebs1d(:,:,17)
        qs1d     => spacebs1d(:,:,18)
        resp1d   => spacebs1d(:,:,19)
        rhs1d    => spacebs1d(:,:,20)
        rno1d    => spacebs1d(:,:,21)
        rnos1d   => spacebs1d(:,:,22)
        rsw1d    => spacebs1d(:,:,23)
        sag1d    => spacebs1d(:,:,24)
        scv1d    => spacebs1d(:,:,25)
        sent1d   => spacebs1d(:,:,26)
        sice1d   => spacebs1d(:,:,27)
        ssw1d    => spacebs1d(:,:,28)
        t2m_1d   => spacebs1d(:,:,29)
        taf1d    => spacebs1d(:,:,30)
        tg1d     => spacebs1d(:,:,31)
        tgb1d    => spacebs1d(:,:,32)
        tlef1d   => spacebs1d(:,:,33)
        ts1d     => spacebs1d(:,:,34)
        tsw1d    => spacebs1d(:,:,35)
        u10m1d   => spacebs1d(:,:,36)
        v10m1d   => spacebs1d(:,:,37)
        veg1d    => spacebs1d(:,:,38)
        z1d      => spacebs1d(:,:,39)
        bfc      => spacebs1d(:,:,40)
        bsw      => spacebs1d(:,:,41)
        evmx0    => spacebs1d(:,:,42)
        fdry     => spacebs1d(:,:,43)
        fwet     => spacebs1d(:,:,44)
        gwmx0    => spacebs1d(:,:,45)
        gwmx1    => spacebs1d(:,:,46)
        gwmx2    => spacebs1d(:,:,47)
        porsl    => spacebs1d(:,:,48)
        relfc    => spacebs1d(:,:,49)
        rnet     => spacebs1d(:,:,50)
        texrat   => spacebs1d(:,:,51)
        vegt     => spacebs1d(:,:,52)
        wiltr    => spacebs1d(:,:,53)
        wt       => spacebs1d(:,:,54)
        xkmx     => spacebs1d(:,:,55)
        aarea    => spacebs1d(:,:,56)
        cdr      => spacebs1d(:,:,57)
        cdrn     => spacebs1d(:,:,58)
        cdrx     => spacebs1d(:,:,59)
        cf       => spacebs1d(:,:,60)
        cgrnd    => spacebs1d(:,:,61)
        cgrndl   => spacebs1d(:,:,62)
        cgrnds   => spacebs1d(:,:,63)
        clead    => spacebs1d(:,:,64)
        densi    => spacebs1d(:,:,65)
        efpr     => spacebs1d(:,:,66)
        eg       => spacebs1d(:,:,67)
        etr      => spacebs1d(:,:,68)
        etrrun   => spacebs1d(:,:,69)
        evaps    => spacebs1d(:,:,70)
        evapw    => spacebs1d(:,:,71)
        fevpg    => spacebs1d(:,:,72)
        flnet    => spacebs1d(:,:,73)
        flneto   => spacebs1d(:,:,74)
        fseng    => spacebs1d(:,:,75)
        htvp     => spacebs1d(:,:,76)
        ps       => spacebs1d(:,:,77)
        pw       => spacebs1d(:,:,78)
        qice     => spacebs1d(:,:,79)
        qsatl    => spacebs1d(:,:,80)
        rhosw    => spacebs1d(:,:,81)
        ribd     => spacebs1d(:,:,82)
        rlai     => spacebs1d(:,:,83)
        rpp      => spacebs1d(:,:,84)
        scrat    => spacebs1d(:,:,85)
        scvk     => spacebs1d(:,:,86)
        sdrop    => spacebs1d(:,:,87)
        seasb    => spacebs1d(:,:,88)
        sigf     => spacebs1d(:,:,89)
        sm       => spacebs1d(:,:,90)
        tm       => spacebs1d(:,:,91)
        uaf      => spacebs1d(:,:,92)
        vspda    => spacebs1d(:,:,93)
        wata     => spacebs1d(:,:,94)
        watr     => spacebs1d(:,:,95)
        watt     => spacebs1d(:,:,96)
        watu     => spacebs1d(:,:,97)
        wta      => spacebs1d(:,:,98)
        xlai     => spacebs1d(:,:,99)
        xlsai    => spacebs1d(:,:,100)
        xrun     => spacebs1d(:,:,101)
        z1       => spacebs1d(:,:,102)
        z1log    => spacebs1d(:,:,103)
        cn1      => spacebs1d(:,:,104)
        df       => spacebs1d(:,:,105)
        rgr      => spacebs1d(:,:,106)
        wta0     => spacebs1d(:,:,107)
        wtaq0    => spacebs1d(:,:,108)
        wtg      => spacebs1d(:,:,109)
        wtg0     => spacebs1d(:,:,110)
        wtg2     => spacebs1d(:,:,111)
        wtga     => spacebs1d(:,:,112)
        wtgaq    => spacebs1d(:,:,113)
        wtgl     => spacebs1d(:,:,114)
        wtglq    => spacebs1d(:,:,115)
        wtgq     => spacebs1d(:,:,116)
        wtgq0    => spacebs1d(:,:,117)
        wtl0     => spacebs1d(:,:,118)
        wtlh     => spacebs1d(:,:,119)
        wtlq     => spacebs1d(:,:,120)
        wtlq0    => spacebs1d(:,:,121)
        wtshi    => spacebs1d(:,:,122)
        wtsqi    => spacebs1d(:,:,123)
        allocate(imelt(nnsg,iym1))
        allocate(lveg(nnsg,iym1))
        allocate(spaceb1d(iym1,24))
        flw1d   => spaceb1d(:,1)
        fsw1d   => spaceb1d(:,2)
        us1d    => spaceb1d(:,3)
        vs1d    => spaceb1d(:,4)
        czen    => spaceb1d(:,5)
        sola    => spaceb1d(:,6)
        vpdd    => spaceb1d(:,7)
        ems     => spaceb1d(:,8)
        albdif  => spaceb1d(:,9)
        albdir  => spaceb1d(:,10)
        albvl   => spaceb1d(:,11)
        albvld  => spaceb1d(:,12)
        albvs   => spaceb1d(:,13)
        albvsd  => spaceb1d(:,14)
        aldifl  => spaceb1d(:,15)
        aldifs  => spaceb1d(:,16)
        aldirl  => spaceb1d(:,17)
        aldirs  => spaceb1d(:,18)
        emiss1d => spaceb1d(:,19)
        fracd   => spaceb1d(:,20)
        sabveg  => spaceb1d(:,21)
        solis   => spaceb1d(:,22)
        solvd   => spaceb1d(:,23)
        solvs   => spaceb1d(:,24)
        allocate(coszrs(iy))
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
        u10m_s => fsub(:,:,:,1)
        v10m_s => fsub(:,:,:,2)
        drag_s => fsub(:,:,:,3)
        tg_s   => fsub(:,:,:,4)
        tlef_s => fsub(:,:,:,5)
        t2m_s  => fsub(:,:,:,6)
        q2m_s  => fsub(:,:,:,7)
        ssw_s  => fsub(:,:,:,8)
        rsw_s  => fsub(:,:,:,9)
        tpr_s  => fsub(:,:,:,10)
        evpa_s => fsub(:,:,:,11)
        rnos_s => fsub(:,:,:,12)
        scv_s  => fsub(:,:,:,13)
        sena_s => fsub(:,:,:,14)
        prcv_s => fsub(:,:,:,15)
        ps_s   => fsub(:,:,:,16)
#ifdef DCSST
        firstcall(:,:) = .false.
#endif
!
      end subroutine allocate_mod_bats 
!
      end module mod_bats
