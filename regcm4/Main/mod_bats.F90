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

      use mod_dynparam
      use mod_param1 , only : dtbat , dtmin
      use mod_bats_param

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
           & swt2d , taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
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
        u10m_o => fbat(:,:,1)
        v10m_o => fbat(:,:,2)
        drag_o => fbat(:,:,3)
        tg_o   => fbat(:,:,4)
        tlef_o => fbat(:,:,5)
        t2m_o  => fbat(:,:,6)
        q2m_o  => fbat(:,:,7)
        ssw_o  => fbat(:,:,8)
        rsw_o  => fbat(:,:,9)
        tpr_o  => fbat(:,:,10)
        evpa_o => fbat(:,:,11)
        rnos_o => fbat(:,:,12)
        scv_o  => fbat(:,:,13)
        sena_o => fbat(:,:,14)
        flwa_o => fbat(:,:,15)
        fswa_o => fbat(:,:,16)
        flwd_o => fbat(:,:,17)
        sina_o => fbat(:,:,18)
        prcv_o => fbat(:,:,19)
        ps_o   => fbat(:,:,20)
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
!
!
      subroutine interf(ivers,j,k,istart,iend,ng)

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  this subroutine interfaces mm42d and bats variables
!
!  ivers = 1 ,   regcm2d --> bats
!  ivers = 2 ,   bats --> regcm2d
!
      use mod_param2 , only : iocnflx , kbats
      use mod_param3 , only : r8pt
      use mod_main , only : psb , hfx , qfx , ht , zpbl , uvdrag ,      &
                  &         snowc , tgbb , tga , tgb
      use mod_pbldim , only : thx3d , za
      use mod_slice , only : ubx3d , vbx3d , qvb3d
      use mod_constants , only : tau1 , zlnd , zoce , zsno , rgti ,     &
                  &              rgas , tzero , lh0 , lh1 , lsvp1 ,     &
                  &              lsvp2 , ep2 , lrate
      use mod_date , only : jyear , jyear0 , jyearr , ntime , ktau ,    &
                  &         ktaur
      implicit none
!
! Dummy arguments
!
      integer , intent (in) :: ivers , j , k , istart , iend , ng
!
! Local variables
!
      real(8) :: amxtem , facb , facs , fact , factuv , facv , fracb ,  &
               & fracs , fracv , hl , mmpd , rh0 , satvp , sfac ,       &
               & solvt , wpm2
      integer :: i , n , nnn
      real(4) :: real_4
!
!     ******    fbat contains nummap+1 fields to be written out to
!     iutbat ******    fields are written out in following order:
!     1.  anemon west  wind (41)
!     2.  anemom south wind (42)
!     3.  drag- surface stress, in si (48)
!     4.  ground temp (4)
!     5.  temp of foliage (7)
!     6.  anemom temp (5)
!     7.  anemom spec. humidity (15)
!     8.  upper layer soil water (26)
!     9.  root zone soil water (25)
!     10.  accum precip (19)
!     11.  accum evap (39)
!     12.  accum surf runoff (12)
!     13.  snow depth in mm h2o (30)
!     14.  accum sensible heat (40)
!     15.  accum net ir (38)
!     16.  accum net solar abs (37)
!     17.  accum downward ir
!     18.  accum solar incident at surface
!     19.  convective precipitation
!     20.  surface pressure
!     21.  pbl height
 
      if ( ivers.eq.1 ) then ! regcm2d --> bats

        do i = istart, iend
          do n = 1 , ng
            p1d0(n,i) = (psb(i,j)+r8pt)*1000.
            z1d(n,i) = za(i,k,j)
            ts1d0(n,i) = thx3d(i,k,j)
            qs1d0(n,i) = qvb3d(i,k,j)/(1.+qvb3d(i,k,j))
            qs1d(n,i) = qs1d0(n,i)
 
            hl = lh0 - lh1*(ts1d0(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(1./tzero-1./ts1d0(n,i)))
            rh0 = dmax1(qs1d0(n,i)/(ep2*satvp/(p1d0(n,i)*0.01-satvp)), &
                & 0.D0)
 
            ts1d(n,i) = ts1d0(n,i) - lrate*rgti*(ht1(n,i,j)-ht(i,j))
            p1d(n,i) = p1d0(n,i)*(ts1d(n,i)/ts1d0(n,i))
 
            hl = lh0 - lh1*(ts1d(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(1./tzero-1./ts1d(n,i)))
            qs1d(n,i) = dmax1(rh0*ep2*satvp/(p1d(n,i)*0.01-satvp),0.D0)
 
            tg1d(n,i) = tg2d(n,i,j)
            rhs1d(n,i) = p1d(n,i)/(rgas*ts1d(n,i))
            prcp1d(n,i) = pptnc(i,j) + pptc(i,j)
!
!           quantities stored on 2d surface array for bats use only
!
            tgb1d(n,i) = tgb2d(n,i,j)
            taf1d(n,i) = taf2d(n,i,j)
            tlef1d(n,i) = tlef2d(n,i,j)
            tsw1d(n,i) = swt2d(n,i,j)
            rsw1d(n,i) = srw2d(n,i,j)
            ssw1d(n,i) = ssw2d(n,i,j)
            ldew1d(n,i) = dew2d(n,i,j)
            sag1d(n,i) = sag2d(n,i,j)
            scv1d(n,i) = scv2d(n,i,j)
            sice1d(n,i) = sice2d(n,i,j)
            gwet1d(n,i) = gwet2d(n,i,j)
            sent1d(n,i) = hfx(i,j)
            evpr1d(n,i) = qfx(i,j)
            ldoc1d(n,i) = ocld2d(n,i,j)
            ircp1d(n,i) = ircp2d(n,i,j)
            lveg(n,i) = nint(veg2d1(n,i,j))
            amxtem = dmax1(298.-tgb1d(n,i),0.D0)
            sfac = 1. - dmax1(0.D0,1.-0.0016*amxtem**2)
            if ( lveg(n,i).eq.0 ) then
              veg1d(n,i) = 0.
            else
              veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
            end if
            emiss_1d(n,i) = emiss2d(n,i,j)
          end do
 
          rh0 = 0.0D0
          do n = 1 , ng
            rh0 = rh0 + (qs1d(n,i)-qs1d0(n,i))
          end do
          rh0 = rh0/ng
          do n = 1 , ng
            qs1d(n,i) = dmax1(qs1d(n,i)-rh0,0.0D0)
          end do
 
          us1d(i) = ubx3d(i,k,j)
          vs1d(i) = vbx3d(i,k,j)
          fsw1d(i) = fsw2d(i,j)
          flw1d(i) = flw2d(i,j)
          solis(i) = sol2d(i,j)
          sabveg(i) = sabv2d(i,j)
          solvt = solvd2d(i,j) + solvs2d(i,j)
          if ( solvt.gt.0.0 ) then
            fracd(i) = solvd2d(i,j)/solvt
          else
            fracd(i) = 0.2
          end if
          czen(i) = dmax1(coszrs(i),0.D0)
        end do
 
      else if ( ivers.eq.2 ) then ! bats --> regcm2d
 
        do i = istart, iend
          uvdrag(i,j) = 0.0
          hfx(i,j) = 0.0
          qfx(i,j) = 0.0
          tgb(i,j) = 0.0
          tga(i,j) = 0.0
          tgbb(i,j) = 0.0
!chem2
          ssw2da(i,j) = 0.0
          sdeltk2d(i,j) = 0.0
          sdelqk2d(i,j) = 0.0
          sfracv2d(i,j) = 0.0
          sfracb2d(i,j) = 0.0
          sfracs2d(i,j) = 0.0
          svegfrac2d(i,j) = 0.0
!chem2_
          do n = 1 , ng
            uvdrag(i,j) = uvdrag(i,j) + drag1d(n,i)
            hfx(i,j) = hfx(i,j) + sent1d(n,i)
            qfx(i,j) = qfx(i,j) + evpr1d(n,i)
            tgb(i,j) = tgb(i,j) + tg1d(n,i)
            tga(i,j) = tga(i,j) + tg1d(n,i)
!chem2
            ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
            sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
            sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
            sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
            sfracb2d(i,j) = sfracb2d(i,j) + (1.-veg1d(n,i))             &
                          & *(1.-scvk(n,i))
            sfracs2d(i,j) = sfracs2d(i,j) + veg1d(n,i)*wt(n,i)          &
                          & + (1.-veg1d(n,i))*scvk(n,i)
            svegfrac2d(i,j) = svegfrac2d(i,j) + veg1d(n,i)
!chem2_
            if ( iocnflx.eq.1 .or.                                      &
               & (iocnflx.eq.2 .and. ocld2d(n,i,j).ge.0.5) ) then
              tgbb(i,j) = tgbb(i,j)                                     &
                        & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)    &
                        & *tlef1d(n,i)**4)**0.25
            else
              tgbb(i,j) = tgbb(i,j) + tg1d(n,i)
            end if
            if ( ocld2d(n,i,j).lt.0.5 ) then
              ssw1d(n,i) = -1.E34
              rsw1d(n,i) = -1.E34
              tsw1d(n,i) = -1.E34
              rno1d(n,i) = -1.E34
              rnos1d(n,i) = -1.E34
              scv1d(n,i) = -1.E34
            end if
          end do
          uvdrag(i,j) = uvdrag(i,j)/float(ng)
          hfx(i,j) = hfx(i,j)/float(ng)
          qfx(i,j) = qfx(i,j)/float(ng)
          tgb(i,j) = tgb(i,j)/float(ng)
          tga(i,j) = tga(i,j)/float(ng)
          tgbb(i,j) = tgbb(i,j)/float(ng)
!chem2
          ssw2da(i,j) = ssw2da(i,j)/float(ng)
          sdeltk2d(i,j) = sdeltk2d(i,j)/float(ng)
          sdelqk2d(i,j) = sdelqk2d(i,j)/float(ng)
          sfracv2d(i,j) = sfracv2d(i,j)/float(ng)
          sfracb2d(i,j) = sfracb2d(i,j)/float(ng)
          sfracs2d(i,j) = sfracs2d(i,j)/float(ng)
          svegfrac2d(i,j) = svegfrac2d(i,j)/float(ng)
!chem2_
          do n = 1 , ng
            snowc(n,i,j) = scv1d(n,i)
            tg2d(n,i,j) = tg1d(n,i)
            tgb2d(n,i,j) = tgb1d(n,i)
            taf2d(n,i,j) = taf1d(n,i)
            tlef2d(n,i,j) = tlef1d(n,i)
            swt2d(n,i,j) = tsw1d(n,i)
            srw2d(n,i,j) = rsw1d(n,i)
            ssw2d(n,i,j) = ssw1d(n,i)
            dew2d(n,i,j) = ldew1d(n,i)
            sag2d(n,i,j) = sag1d(n,i)
            scv2d(n,i,j) = scv1d(n,i)
            sice2d(n,i,j) = sice1d(n,i)
            gwet2d(n,i,j) = gwet1d(n,i)
            ocld2d(n,i,j) = ldoc1d(n,i)
            ircp2d(n,i,j) = ircp1d(n,i)
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
            if ( rnos2d(n,i,j).gt.-1.E10 .and. rnos1d(n,i).gt.-1.E10 )  &
               & then
              rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)/tau1*dtbat
            else
              rnos2d(n,i,j) = -1.E34
            end if
            if ( rno2d(n,i,j).gt.-1.E10 .and. rnos1d(n,i)               &
               & .gt.-1.E10 .and. rno1d(n,i).gt.-1.E10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + (rno1d(n,i)-rnos1d(n,i))    &
                           & /tau1*dtbat
            else
              rno2d(n,i,j) = -1.E34
            end if
          end do
!
!         quantities stored on 2d surface array for bats use only
!
          prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
          prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
          if ( prnca2d(i,j) < 1E-30 ) prnca2d(i,j) = 0.0
          if ( prca2d(i,j) < 1E-30 ) prca2d(i,j) = 0.0
          flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
          flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
          fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
          svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
          sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
          pptnc(i,j) = 0.
          pptc(i,j) = 0.
        end do

        do i = istart, iend
#ifdef MPP1
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          do n = 1 , ng
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j,i-1) = tg1d(n,i)
            u10m_s(n,j,i-1) = u10m1d(n,i)
            v10m_s(n,j,i-1) = v10m1d(n,i)
            t2m_s(n,j,i-1) = t2m_1d(n,i)
 
            u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
            v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
            t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
            tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/float(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/float(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/float(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/float(ng)
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**2+    &
                        & v10m_o(j,i-1)**2))
          real_4 = (psb(i,j)+r8pt)*10.
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
#else
#ifdef BAND
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          do n = 1 , ng
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j,i-1) = tg1d(n,i)
            u10m_s(n,j,i-1) = u10m1d(n,i)
            v10m_s(n,j,i-1) = v10m1d(n,i)
            t2m_s(n,j,i-1) = t2m_1d(n,i)

            u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
            v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
            t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
            tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/float(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/float(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/float(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/float(ng)
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**&
                          & 2+v10m_o(j,i-1)**2))
          real_4 = (psb(i,j)+r8pt)*10.
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
#else
          u10m_o(j-1,i-1) = 0.0
          v10m_o(j-1,i-1) = 0.0
          tg_o(j-1,i-1) = 0.0
          t2m_o(j-1,i-1) = 0.0
          do n = 1 , ng
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j-1,i-1) = tg1d(n,i)
            u10m_s(n,j-1,i-1) = u10m1d(n,i)
            v10m_s(n,j-1,i-1) = v10m1d(n,i)
            t2m_s(n,j-1,i-1) = t2m_1d(n,i)

            u10m_o(j-1,i-1) = u10m_o(j-1,i-1) + u10m1d(n,i)
            v10m_o(j-1,i-1) = v10m_o(j-1,i-1) + v10m1d(n,i)
            t2m_o(j-1,i-1) = t2m_o(j-1,i-1) + t2m_1d(n,i)
            tg_o(j-1,i-1) = tg_o(j-1,i-1) + tg1d(n,i)
          end do
          u10m_o(j-1,i-1) = u10m_o(j-1,i-1)/float(ng)
          v10m_o(j-1,i-1) = v10m_o(j-1,i-1)/float(ng)
          t2m_o(j-1,i-1) = t2m_o(j-1,i-1)/float(ng)
          tg_o(j-1,i-1) = tg_o(j-1,i-1)/float(ng)
          tgmx_o(j-1,i-1) = amax1(tgmx_o(j-1,i-1),tg_o(j-1,i-1))
          tgmn_o(j-1,i-1) = amin1(tgmn_o(j-1,i-1),tg_o(j-1,i-1))
          t2mx_o(j-1,i-1) = amax1(t2mx_o(j-1,i-1),t2m_o(j-1,i-1))
          t2mn_o(j-1,i-1) = amin1(t2mn_o(j-1,i-1),t2m_o(j-1,i-1))
          w10x_o(j-1,i-1) = amax1(w10x_o(j-1,i-1),sqrt(u10m_o(j-1,i-1)**&
                          & 2+v10m_o(j-1,i-1)**2))
          real_4 = (psb(i,j)+r8pt)*10.
          psmn_o(j-1,i-1) = amin1(psmn_o(j-1,i-1),real_4)
#endif
#endif
        end do

        if ( mod(ntime+nint(dtmin*60.),kbats).eq.0 .or.                 &
           & (jyear.eq.jyearr .and. ktau.eq.ktaur) ) then
          if ( jyear.eq.jyear0 .and. ktau.le.1 ) then
            mmpd = 86400./dtbat
            wpm2 = 1./dtbat
          else if ( jyear.eq.jyear0 .and. dble(ktau*dtmin)              &
                  & .le.batfrq*60.+0.01 ) then
            mmpd = 24./(batfrq-dtmin/60.)
            wpm2 = 1./((batfrq-dtmin/60.)*3600.)
          else
            mmpd = 24./batfrq
            wpm2 = 1./(batfrq*3600.)
          end if
          do i = istart, iend
#ifdef MPP1
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j,i-1) = q2m_1d(n,i)
              drag_s(n,j,i-1) = drag1d(n,i)
              evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j,i-1) = p1d(n,i)*0.01
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
              drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
              evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
              sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j,i-1) = q2m_o(j,i-1)/float(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/float(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/float(ng)*mmpd
            sena_o(j,i-1) = sena_o(j,i-1)/float(ng)*wpm2
            flwa_o(j,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j,i-1) = flwda2d(i,j)*wpm2
            sina_o(j,i-1) = sina2d(i,j)*wpm2
            prcv_o(j,i-1) = prca2d(i,j)*mmpd
            ps_o(j,i-1) = (psb(i,j)+r8pt)*10.
            zpbl_o(j,i-1) = zpbl(i,j)
 
            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + tlef1d(n,i)
                ssw_o(j,i-1) = ssw_o(j,i-1) + ssw1d(n,i)
                rsw_o(j,i-1) = rsw_o(j,i-1) + rsw1d(n,i)
                rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,i,j)
                if (abs(scv1d(n,i)) > 1E-30) then
                  scv_o(j,i-1) = scv_o(j,i-1) + scv1d(n,i)
                end if
                tlef_s(n,j,i-1) = tlef1d(n,i)
                ssw_s(n,j,i-1) = ssw1d(n,i)
                rsw_s(n,j,i-1) = rsw1d(n,i)
                rnos_s(n,j,i-1) = rnos2d(n,i,j)*mmpd
                if (abs(scv1d(n,i)) > 1E-30) then
                  scv_s(n,j,i-1) = scv1d(n,i)
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = -1.E34
                ssw_s(n,j,i-1) = -1.E34
                rsw_s(n,j,i-1) = -1.E34
                rnos_s(n,j,i-1) = -1.E34
                scv_s(n,j,i-1) = -1.E34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
              scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
            else
              tlef_o(j,i-1) = -1.E34
              ssw_o(j,i-1) = -1.E34
              rsw_o(j,i-1) = -1.E34
              rnos_o(j,i-1) = -1.E34
              scv_o(j,i-1) = -1.E34
            end if
#else
#ifdef BAND
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j,i-1) = q2m_1d(n,i)
              drag_s(n,j,i-1) = drag1d(n,i)
              evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j,i-1) = p1d(n,i)*0.01

              q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
              drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
              evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
              sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j,i-1) = q2m_o(j,i-1)/float(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/float(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/float(ng)*mmpd
            sena_o(j,i-1) = sena_o(j,i-1)/float(ng)*wpm2
            flwa_o(j,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j,i-1) = flwda2d(i,j)*wpm2
            sina_o(j,i-1) = sina2d(i,j)*wpm2
            prcv_o(j,i-1) = prca2d(i,j)*mmpd
            ps_o(j,i-1) = (psb(i,j)+r8pt)*10.
            zpbl_o(j,i-1) = zpbl(i,j)

            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + tlef1d(n,i)
                ssw_o(j,i-1) = ssw_o(j,i-1) + ssw1d(n,i)
                rsw_o(j,i-1) = rsw_o(j,i-1) + rsw1d(n,i)
                rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,i,j)
                scv_o(j,i-1) = scv_o(j,i-1) + scv1d(n,i)
                tlef_s(n,j,i-1) = tlef1d(n,i)
                ssw_s(n,j,i-1) = ssw1d(n,i)
                rsw_s(n,j,i-1) = rsw1d(n,i)
                rnos_s(n,j,i-1) = rnos2d(n,i,j)*mmpd
                scv_s(n,j,i-1) = scv1d(n,i)
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = -1.E34
                ssw_s(n,j,i-1) = -1.E34
                rsw_s(n,j,i-1) = -1.E34
                rnos_s(n,j,i-1) = -1.E34
                scv_s(n,j,i-1) = -1.E34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
              scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
            else
              tlef_o(j,i-1) = -1.E34
              ssw_o(j,i-1) = -1.E34
              rsw_o(j,i-1) = -1.E34
              rnos_o(j,i-1) = -1.E34
              scv_o(j,i-1) = -1.E34
            end if
#else
            drag_o(j-1,i-1) = 0.0
            q2m_o(j-1,i-1) = 0.0
            evpa_o(j-1,i-1) = 0.0
            sena_o(j-1,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j-1,i-1) = q2m_1d(n,i)
              drag_s(n,j-1,i-1) = drag1d(n,i)
              evpa_s(n,j-1,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j-1,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j-1,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j-1,i-1) = p1d(n,i)*0.01

              q2m_o(j-1,i-1) = q2m_o(j-1,i-1) + q2m_1d(n,i)
              drag_o(j-1,i-1) = drag_o(j-1,i-1) + drag1d(n,i)
              evpa_o(j-1,i-1) = evpa_o(j-1,i-1) + evpa2d(n,i,j)
              sena_o(j-1,i-1) = sena_o(j-1,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j-1,i-1) = q2m_o(j-1,i-1)/float(ng)
            drag_o(j-1,i-1) = drag_o(j-1,i-1)/float(ng)
            evpa_o(j-1,i-1) = evpa_o(j-1,i-1)/float(ng)*mmpd
            sena_o(j-1,i-1) = sena_o(j-1,i-1)/float(ng)*wpm2
            flwa_o(j-1,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j-1,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j-1,i-1) = flwda2d(i,j)*wpm2
            sina_o(j-1,i-1) = sina2d(i,j)*wpm2
            prcv_o(j-1,i-1) = prca2d(i,j)*mmpd
            ps_o(j-1,i-1) = (psb(i,j)+r8pt)*10.
            zpbl_o(j-1,i-1) = zpbl(i,j)

            tlef_o(j-1,i-1) = 0.0
            ssw_o(j-1,i-1) = 0.0
            rsw_o(j-1,i-1) = 0.0
            rnos_o(j-1,i-1) = 0.0
            scv_o(j-1,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j-1,i-1) = tlef_o(j-1,i-1) + tlef1d(n,i)
                ssw_o(j-1,i-1) = ssw_o(j-1,i-1) + ssw1d(n,i)
                rsw_o(j-1,i-1) = rsw_o(j-1,i-1) + rsw1d(n,i)
                rnos_o(j-1,i-1) = rnos_o(j-1,i-1) + rnos2d(n,i,j)
                if (abs(scv1d(n,i)) > 1E-34) then
                  scv_o(j-1,i-1) = scv_o(j-1,i-1) + scv1d(n,i)
                end if
                tlef_s(n,j-1,i-1) = tlef1d(n,i)
                ssw_s(n,j-1,i-1) = ssw1d(n,i)
                rsw_s(n,j-1,i-1) = rsw1d(n,i)
                rnos_s(n,j-1,i-1) = rnos2d(n,i,j)*mmpd
                if (abs(scv1d(n,i)) > 1E-34) then
                  scv_s(n,j-1,i-1) = scv1d(n,i)
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j-1,i-1) = -1.E34
                ssw_s(n,j-1,i-1) = -1.E34
                rsw_s(n,j-1,i-1) = -1.E34
                rnos_s(n,j-1,i-1) = -1.E34
                scv_s(n,j-1,i-1) = -1.E34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j-1,i-1) = tlef_o(j-1,i-1)/float(nnn)
              ssw_o(j-1,i-1) = ssw_o(j-1,i-1)/float(nnn)
              rsw_o(j-1,i-1) = rsw_o(j-1,i-1)/float(nnn)
              rnos_o(j-1,i-1) = rnos_o(j-1,i-1)/float(nnn)*mmpd
              scv_o(j-1,i-1) = scv_o(j-1,i-1)/float(nnn)
            else
              tlef_o(j-1,i-1) = -1.E34
              ssw_o(j-1,i-1) = -1.E34
              rsw_o(j-1,i-1) = -1.E34
              rnos_o(j-1,i-1) = -1.E34
              scv_o(j-1,i-1) = -1.E34
            end if
#endif
#endif
 
!           ******    reset accumulation arrays to zero
            do n = 1 , ng
              evpa2d(n,i,j) = 0.
              rnos2d(n,i,j) = 0.
              sena2d(n,i,j) = 0.
            end do
            prnca2d(i,j) = 0.
            prca2d(i,j) = 0.
            flwa2d(i,j) = 0.
            flwda2d(i,j) = 0.
            fswa2d(i,j) = 0.
            svga2d(i,j) = 0.
            sina2d(i,j) = 0.
          end do
        end if
!
      else ! end ivers test
      end if
!
      end subroutine interf
!
!     provides leaf and stem area parameters;
!     depends on climate through subsoil temperatures.
!
      subroutine vcover
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      implicit none
!
! Local variables
!
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) seasb(n,i)                        &
               & = dmax1(0.D0,1.-0.0016*dmax1(298.-tgb1d(n,i),0.D0)**2)
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              xlai(n,i) = xla(lveg(n,i))
              xlai(n,i) = xlai(n,i) + (xlai0(lveg(n,i))-xlai(n,i))      &
                         & *(1.-seasb(n,i))
              rlai(n,i) = xlai(n,i) + sai(lveg(n,i))
              xlsai(n,i) = xlai(n,i) + sai(lveg(n,i))
              vegt(n,i) = sigf(n,i)*xlsai(n,i)
            end if
          end if
        end do
      end do
!
      end subroutine vcover
!
!
!     this subroutine provides root function in terms of maximum
!     transpiration rate plants can sustain depending on soil moisture.
!
!     trsmx0 is a prescribed constant (kg/m**2/s).
!     trsmx is the maximum transpiration rate,
!        including a low temperature correction (=seasb)
!        and a correction for fractional vegetation (=sigf).
!
!     rotf is ratio of moisture extracton from top to total when
!           fully saturated
!     rootf is ratio of roots in upper soil layer
!                    to roots in root soil layer
!     bsw is the b param in clapp and hornberger
!
!     "wlt  " are ratios factors controlling the saturation
!                 cf wilting (see ewing paper)
!     wlttb (total) & wltub (upper) become 1 at the wilting point
!     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
!
!     etrc= max poss transpiration given the soil moisture distributions
!     efpr = the relative contribution of upper soil layer to
!     evapotranspiration - need soil moist. budget (subrout water)
!
      subroutine root
      use mod_ictp01
      use mod_constants , only : trsmx0
      implicit none
!
! Local variables
!
      real(8) :: bneg , rotf , trsmx , wlttb , wltub , wmli
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!             trsmx = trsmx0*sigf(n,i)*seasb(n,i)
              trsmx = trsmx0*sigf(n,i)
              rotf = rootf(lveg(n,i))
              bneg = -bsw(n,i)
              wmli = 1./(wiltr(n,i)**bneg-1.)
              wlttb = (watr(n,i)**bneg-1.)*wmli
              wltub = (watu(n,i)**bneg-1.)*wmli
              wlttb = dmin1(wlttb,1.D0)
              wltub = dmin1(wltub,1.D0)
              etrc(n,i) = trsmx*(1.-(1.-rotf)*wlttb-rotf*wltub)
              efpr(n,i) = trsmx*rotf*(1.-wltub)
              if ( etrc(n,i).lt.1.E-12 ) then
                etrc(n,i) = 1.E-12
                efpr(n,i) = 1.0
              else
                efpr(n,i) = efpr(n,i)/etrc(n,i)
              end if
            end if
          end if
        end do
      end do
!
      end subroutine root
!
!              update soil moisture and runoff
!
!     new algorithms for three soil layers (dickinson & kennedy 8-88)
!     calculate fluxes through air, surface layer, and root layer faces
!
!          b = b of clapp and hornberger
!       est0 = soil water flux, out top
!      gwatr = net input of water to the soil surface
!       ircp = leaf interception
!     wflux1 = soil water flux, 10 cm
!     wflux2 = soil water flux, 1 m
!     rsubss = soil water flux by grav. drainage thru 10 cm interface
!     rsubsr = soil water flux by grav. drainage thru 1 m interface
!     rsubst = soil water flux by grav. drainage thru 10 m interface
!       rsur = surface runoff
!     rno1d(n,i) = total runoff (mm/day)
!     rnos1d(n,i) = surface runoff (mm/day)
!
!     xkmxr and wflux1 determine flow thru upper/root soil interface
!     evmxt, xkmx1, and xkmx2 determine flow thru lower interfaces
!
!     veg type 10 "irrigated crop" is irrigated through reducing
!          the runoff (rsur), i.e., by adding a negative number
!          if the land isn't at least 60% saturated.
!     veg type 13 and 14 are water covered (lake, swamp, rice paddy);
!          negative runoff keeps this land saturated.
!
      subroutine water
!
      use mod_constants , only : drain , tau1 , csoilc , tzero
      implicit none
!
! Local variables
!
      real(8) :: b , bfac , bfac2 , delwat , est0 , evmax , evmxr ,     &
               & evmxt , rap , vakb , wtg2c , xxkb
      real(8) , dimension(nnsg,iym1) :: gwatr , rnof , rsubsr ,         &
           & rsubss , rsubst , rsur , wflux1 , wflux2 , wfluxc , xkmx1 ,&
           & xkmx2 , xkmxr
      integer :: n , i
!
!***********************************************************************
!
 
!=======================================================================
!     1.   define soil water fluxes
!=======================================================================
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
!
!           1.1  reduce infiltration for frozen ground
!
            if ( tgb1d(n,i).gt.tzero ) then
              xkmxr(n,i) = xkmx(n,i)
            else
              xkmxr(n,i) = 0.
            end if
!
!           1.11 permafrost or ice sheet
!
            if ( lveg(n,i).eq.9 .or. lveg(n,i).eq.12 ) then
              xkmx1(n,i) = 0.
              xkmx2(n,i) = 0.
            else
              xkmx1(n,i) = xkmx(n,i)
              xkmx2(n,i) = drain
            end if
!
!           1.2  diffusive fluxes
!
            evmxr = evmx0(n,i)*xkmxr(n,i)/xkmx(n,i)
            evmxt = evmx0(n,i)*xkmx1(n,i)/xkmx(n,i)
            b = bsw(n,i)
            bfac = watr(n,i)**(3.+bfc(n,i))*watu(n,i)**(b-bfc(n,i)-1)
            bfac2 = watt(n,i)**(2.+bfc(n,i))*watr(n,i)**(b-bfc(n,i))
            wfluxc(n,i) = evmxr*(depuv(lveg(n,i))/deprv(lveg(n,i)))     &
                         & **0.4*bfac
            wflux1(n,i) = wfluxc(n,i)*watr(n,i)
            wflux2(n,i) = evmxt*(depuv(lveg(n,i))/deprv(lveg(n,i)))     &
                         & **0.5*bfac2*(watt(n,i)-watr(n,i))
!
!           1.3  gravitational drainage
!
            rsubss(n,i) = xkmxr(n,i)*watr(n,i)**(b+0.5)*watu(n,i)       &
                         & **(b+2.5)
            rsubsr(n,i) = xkmx1(n,i)*watt(n,i)**(b+0.5)*watr(n,i)       &
                         & **(b+2.5)
            rsubst(n,i) = xkmx2(n,i)*watt(n,i)**(2.*b+3.)
!
!           1.32 bog or water
!
            if ( (lveg(n,i).ge.13) .and. (lveg(n,i).le.15) ) then
              rsubst(n,i) = 0.0
              rsubss(n,i) = 0.0
              rsubsr(n,i) = 0.0
            end if
!
!           1.4  fluxes through internal surfaces
!
            wflux1(n,i) = wflux1(n,i) - rsubss(n,i)
            wflux2(n,i) = wflux2(n,i) - rsubsr(n,i)
          end if
        end do
      end do
!
!     1.5  net flux at air interface
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            gwatr(n,i) = pw(n,i) - evapw(n,i) + sm(n,i)                 &
                        & + etrrun(n,i)/dtbat
!
!=======================================================================
!           2.   define runoff terms
!=======================================================================
!
!           2.1  surface runoff
!
            wata(n,i) = 0.5*(watu(n,i)+watr(n,i))
!
!           2.11 increase surface runoff over frozen ground
!
            if ( tg1d(n,i).lt.tzero ) then
              rsur(n,i) = dmin1(1.D0,wata(n,i)**1)*                     &
                         & dmax1(0.D0,gwatr(n,i))
            else
              rsur(n,i) = dmin1(1.D0,wata(n,i)**4)                      &
                         & *dmax1(0.D0,gwatr(n,i))
            end if
!
!           2.12 irrigate cropland
!
            if ( lveg(n,i).eq.10 .and. watr(n,i).lt.relfc(n,i) )        &
               & rsur(n,i) = rsur(n,i) + (rsw1d(n,i)-relfc(n,i)*        &
               &             gwmx1(n,i))/dtbat
!
!           2.13 saturate swamp or rice paddy
!
            if ( (lveg(n,i).ge.13) .and. (lveg(n,i).lt.16) )            &
               & rsur(n,i) = rsur(n,i) + dmin1(0.D0,(rsw1d(n,i)-        &
               &             gwmx1(n,i))/dtbat)
!
!           2.2  total runoff
!
            rnof(n,i) = rsur(n,i) + rsubst(n,i)
!
!=======================================================================
!           3.   increment soil moisture
!=======================================================================
!
!           3.1  update top layer with implicit treatment
!           of flux from below
!
            ssw1d(n,i) = ssw1d(n,i) + dtbat*(gwatr(n,i)-efpr(n,i)*      &
                        & etr(n,i)-rsur(n,i)+wflux1(n,i))
            ssw1d(n,i) = ssw1d(n,i)/(1.+wfluxc(n,i)*dtbat/gwmx0(n,i))
!
!           3.2  update root zone
!
            rsw1d(n,i) = rsw1d(n,i) + dtbat*(gwatr(n,i)-etr(n,i)-       &
                        & rsur(n,i)+wflux2(n,i))
!
!           3.3  update total water
!
            tsw1d(n,i) = tsw1d(n,i) + dtbat*(gwatr(n,i)-etr(n,i)-       &
                        & rnof(n,i))
          end if
        end do
      end do
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
!
!=======================================================================
!           4.   check whether soil water exceeds maximum capacity or
!           becomes negative (should rarely or never happen)
!=======================================================================
!
!           4.1  surface water assumed to move downward into soil
!
            if ( ssw1d(n,i).gt.gwmx0(n,i) ) ssw1d(n,i) = gwmx0(n,i)
!
!           4.2  excess root layer water assumed to move downward
!
            if ( rsw1d(n,i).gt.gwmx1(n,i) ) rsw1d(n,i) = gwmx1(n,i)
!
!           4.3  excess total water assumed to go to subsurface runoff
!
            if ( tsw1d(n,i).gt.gwmx2(n,i) ) then
              delwat = tsw1d(n,i) - gwmx2(n,i)
              tsw1d(n,i) = gwmx2(n,i)
              rsubst(n,i) = rsubst(n,i) + delwat/dtbat
            end if
!
!           4.4  check for negative water in top layer
!
            if ( ssw1d(n,i).le.1.E-2 ) ssw1d(n,i) = 1.E-2
!
!=======================================================================
!           5.   accumulate leaf interception
!=======================================================================
!
            ircp1d(n,i) = ircp1d(n,i) + sigf(n,i)*(dtbat*               &
                         & prcp1d(n,i)) - (sdrop(n,i)+etrrun(n,i))
!
!=======================================================================
!           6.   evaluate runoff (incremented in ccm)
!=======================================================================
!
!*          update total runoff
!
            rnof(n,i) = rsur(n,i) + rsubst(n,i)
            rno1d(n,i) = rnof(n,i)*tau1
            rnos1d(n,i) = rsur(n,i)*tau1
          else                       ! ocean or sea ice
            rnof(n,i) = 0.
            rno1d(n,i) = 0.
            rnos1d(n,i) = 0.
          end if
        end do
      end do
!
!=======================================================================
!     7.   calculate potential evaporation and use mod_to determine
!     wetness factor, allowing for snow being saturated
!=======================================================================
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            xxkb = dmin1(rough(lveg(n,i)),1.D0)
            vakb = (1.-sigf(n,i))*vspda(n,i) + sigf(n,i)                &
                 & *(xxkb*uaf(n,i)+(1.-xxkb)*vspda(n,i))
            wtg2c = (1.-sigf(n,i))*cdrx(n,i)*vakb
            rap = rhs1d(n,i)*(csoilc*uaf(n,i)*sigf(n,i)*(qg1d(n,i)+     &
                & delq1d(n,i)-qs1d(n,i))+wtg2c*(qg1d(n,i)-qs1d(n,i)))
            bfac = watr(n,i)**(3.+bfc(n,i))*watu(n,i)                   &
                 & **(bsw(n,i)-bfc(n,i)-1)
            est0 = evmx0(n,i)*bfac*watu(n,i)
            evmax = dmax1(est0,0.D0)
            gwet1d(n,i) = dmin1(1.D0,evmax/dmax1(1.D-14,rap))
            gwet1d(n,i) = scvk(n,i) + gwet1d(n,i)*(1.0-scvk(n,i))
          end if
        end do
      end do
!
      end subroutine water
! 
!     update snow cover and snow age
!
!     three-part if block:
!       if snow cover < 0, then snow cover and snow age = 0
!       if antarctica, snow age = 0 (katabatic winds keep snow fresh)
!       if elsewhere, snow age follows given formulae
!
!        ps = snow precipitation rate
!     evaps = moisture flux from ground to atmosphere
!        sm = snow melt rate
!     sdrop = snow fallen from vegetation
!
!     aging of snow consists of three factors:
!           age1: snow crystal growth
!           age2: surface melting
!           age3: accumulation  of other particles, soot, etc., which
!                      is small in southern hemisphere
!
      subroutine snow
!
      use mod_param1 , only : dtbat
      use mod_constants , only : tzero
      implicit none
!
! Local variables
!
      real(8) :: age1 , age2 , age3 , arg , arg2 , dela , dela0 , dels ,&
               & sge , tage
      integer :: n , i
      real(8) , dimension(nnsg,iym1) :: sold
!
      age3 = 0.3
 
!=======================================================================
!     1.   partition soil evaporation and precipitation
!     between water and snow
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
 
            evapw(n,i) = fevpg(n,i)
            evaps(n,i) = scvk(n,i)*evapw(n,i)
            if ( ldoc1d(n,i).gt.1.5 ) evaps(n,i) = fevpg(n,i)
            evapw(n,i) = (1.-scvk(n,i))*evapw(n,i)
!
!           ******                tm  is temperature of precipitation
            if ( tm(n,i).ge.tzero ) then
              pw(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
              ps(n,i) = 0.0
            else
!             ******                snowing
              pw(n,i) = 0.0
              ps(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
            end if
          end if
        end do
      end do
!
!=======================================================================
!     2.   update snow cover
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            sold(n,i) = scv1d(n,i)
            scv1d(n,i) = scv1d(n,i) + dtbat                             &
                        & *(ps(n,i)-evaps(n,i)-sm(n,i)) + sdrop(n,i)
            scv1d(n,i) = dmax1(scv1d(n,i),0.D0)
            sag1d(n,i) = dmax1(sag1d(n,i),0.D0)
 
!           ******           snow cover except for antarctica
!=======================================================================
!           3.   increment non-dimensional "age" of snow;
!           10 mm snow restores surface to that of new snow.
!=======================================================================
            if ( scv1d(n,i).gt.0. ) then
              arg = 5.E3*(1./tzero-1./tg1d(n,i))
              age1 = dexp(arg)
              arg2 = dmin1(0.D0,10.*arg)
              age2 = dexp(arg2)
              tage = age1 + age2 + age3
              dela0 = 1.E-6*dtbat
              dela = dela0*tage
              dels = 0.1*dmax1(0.D0,scv1d(n,i)-sold(n,i))
              sge = (sag1d(n,i)+dela)*(1.0-dels)
              sag1d(n,i) = dmax1(0.D0,sge)
            end if
 
!           ******           antarctica
            if ( scv1d(n,i).gt.800. ) sag1d(n,i) = 0.
          end if
        end do
      end do
 
      end subroutine snow
      end module mod_bats
