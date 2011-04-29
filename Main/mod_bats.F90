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
      use mod_runparams
      use mod_bats_param
!
      real(8) , allocatable , target , dimension(:,:) :: spaceb1d
      real(8) , allocatable , target , dimension(:,:,:) :: spacebs1d
      private :: spacebs1d , spaceb1d
!
      real(8) , pointer , dimension(:,:) :: p1d0 , qs1d0 , ts1d0
!
      integer , allocatable , dimension(:,:) :: ldoc1d
      real(8) , pointer , dimension(:,:) :: delq1d , delt1d , drag1d ,  &
           & emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d ,             &
           & p1d , pbp1d , prcp1d , q2m_1d , qg1d , qs1d , resp1d ,     &
           & rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
           & sice1d , ssw1d , t2m_1d , taf1d , tg1d , tgb1d , tlef1d ,  &
           & ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d , aldirs1d , &
           & aldifs1d
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
           & xlsai , xrun , z1log , z2fra , z10fra , zlgocn , zlglnd ,  &
           & zlgsno , zlgveg , zlgdis
!
      real(8) , pointer , dimension(:,:) :: cn1 , rgr , wta0 , wtaq0 ,  &
           & wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq , wtgq ,   &
           & wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi , df
!
      integer , allocatable , dimension(:,:) :: lveg , oveg
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
      real(8) , allocatable , target , dimension(:,:,:) :: spaceb2d
      private :: spaceb2d
!
      real(8) , pointer , dimension(:,:) :: flw2d , flwa2d , flwd2d , &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sina2d , sinc2d , sol2d ,&
                                    & solvd2d , solvs2d , svga2d
!
      real(8) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
                               sdelqk2d , sfracv2d , sfracb2d , &
                               sfracs2d , svegfrac2d
!
      real(8) , allocatable , target , dimension(:,:,:,:) :: spacebs2d
      private :: spacebs2d
!
      integer , allocatable , dimension(:,:,:) :: ocld2d , veg2d1
      integer , allocatable , dimension(:,:) :: veg2d
!
      real(8) , pointer, dimension(:,:,:) :: col2d , dew2d ,        &
           & emiss2d , evpa2d , gwet2d , ircp2d , rno2d ,      &
           & rnos2d , sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , &
           & swt2d , taf2d , tg2d , tgb2d , tlef2d
!
      real(8) ,allocatable, dimension(:,:,:) :: ht1 , satbrt1 , xlat1 , &
                                             &  xlon1
!
      real(4) , target , allocatable, dimension(:,:,:) :: fbat
!
      real(4) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
                                     & flwd_o , fswa_o , prcv_o ,       &
                                     & psmn_o , ps_o , q2m_o , rnos_o , &
                                     & rsw_o , scv_o , sena_o , sina_o ,&
                                     & ssw_o , t2mn_o , t2mx_o , t2m_o ,&
                                     & tgmn_o , tgmx_o , tg_o , tlef_o ,&
                                     & tpr_o , u10m_o , v10m_o ,        &
                                     & w10x_o , zpbl_o , aldirs_o ,  &
                                     & aldifs_o
!
      real(4) , target ,allocatable, dimension(:,:,:,:) :: fsub
!
      real(4) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
!
      ! dtskin is difference between skin temp and bulk sst
      real(8) , allocatable , dimension(:,:) :: deltas , tdeltas ,      &
                           &                    dtskin
      logical , allocatable , dimension(:,:) :: firstcall

      contains

        subroutine allocate_mod_bats(lmpi,lband)
        implicit none
        logical , intent(in) :: lmpi , lband
        integer :: nj , njm1 , njm2

        if (lmpi) then
          nj = jxp
          njm1 = jxp
          njm2 = jxp
        else
          if (lband) then
            nj = jx
            njm1 = jx
            njm2 = jx
          else
            nj = jx
            njm1 = jxm1
            njm2 = jxm2
          end if
        end if

        allocate(veg2d(iym1,njm1))
        veg2d = -1

        if ( ichem == 1 ) then
          allocate(spaceb2d(iym1,njm1,24))
        else
          allocate(spaceb2d(iym1,njm1,17))
        end if
        spaceb2d = d_zero

        flw2d      => spaceb2d(:,:,1)
        flwa2d     => spaceb2d(:,:,2)
        flwd2d     => spaceb2d(:,:,3)
        flwda2d    => spaceb2d(:,:,4)
        fsw2d      => spaceb2d(:,:,5)
        fswa2d     => spaceb2d(:,:,6)
        pptc       => spaceb2d(:,:,7)
        pptnc      => spaceb2d(:,:,8)
        prca2d     => spaceb2d(:,:,9)
        prnca2d    => spaceb2d(:,:,10)
        sabv2d     => spaceb2d(:,:,11)
        sina2d     => spaceb2d(:,:,12)
        sinc2d     => spaceb2d(:,:,13)
        sol2d      => spaceb2d(:,:,14)
        solvd2d    => spaceb2d(:,:,15)
        solvs2d    => spaceb2d(:,:,16)
        svga2d     => spaceb2d(:,:,17)
        if ( ichem == 1 ) then
          ssw2da     => spaceb2d(:,:,18)
          sdeltk2d   => spaceb2d(:,:,19)
          sdelqk2d   => spaceb2d(:,:,20)
          sfracv2d   => spaceb2d(:,:,21)
          sfracb2d   => spaceb2d(:,:,22)
          sfracs2d   => spaceb2d(:,:,23)
          svegfrac2d => spaceb2d(:,:,24)
        end if

        allocate(ocld2d(nnsg,iym1,njm1))
        ocld2d = -1
        allocate(veg2d1(nnsg,iym1,njm1))
        veg2d1 = -1

        allocate(spacebs2d(nnsg,iym1,njm1,19))
        spacebs2d = d_zero
        col2d    => spacebs2d(:,:,:,1)
        dew2d    => spacebs2d(:,:,:,2)
        emiss2d  => spacebs2d(:,:,:,3)
        evpa2d   => spacebs2d(:,:,:,4)
        gwet2d   => spacebs2d(:,:,:,5)
        ircp2d   => spacebs2d(:,:,:,6)
        rno2d    => spacebs2d(:,:,:,7)
        rnos2d   => spacebs2d(:,:,:,8)
        sag2d    => spacebs2d(:,:,:,9)
        scv2d    => spacebs2d(:,:,:,10)
        sena2d   => spacebs2d(:,:,:,11)
        sice2d   => spacebs2d(:,:,:,12)
        srw2d    => spacebs2d(:,:,:,13)
        ssw2d    => spacebs2d(:,:,:,14)
        swt2d    => spacebs2d(:,:,:,15)
        taf2d    => spacebs2d(:,:,:,16)
        tg2d     => spacebs2d(:,:,:,17)
        tgb2d    => spacebs2d(:,:,:,18)
        tlef2d   => spacebs2d(:,:,:,19)

        allocate(ht1(nnsg,iy,nj))
        allocate(satbrt1(nnsg,iy,nj))
        allocate(xlat1(nnsg,iy,nj))
        allocate(xlon1(nnsg,iy,nj))
        ht1 = d_zero
        satbrt1 = d_zero
        xlat1 = d_zero
        xlon1 = d_zero

        if (idcsst == 1) then
          allocate(deltas(iy,nj))
          allocate(tdeltas(iy,nj))
          allocate(dtskin(iy,nj))
          allocate(firstcall(iy,nj))
          deltas = d_zero
          tdeltas = d_zero
          dtskin = d_zero
          firstcall = .false.
        end if

        allocate(ldoc1d(nnsg,iym1))
        ldoc1d = -1

        allocate(spacebs1d(nnsg,iym1,130))
        spacebs1d = d_zero
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
        p1d      => spacebs1d(:,:,12)
        pbp1d    => spacebs1d(:,:,13)
        prcp1d   => spacebs1d(:,:,14)
        q2m_1d   => spacebs1d(:,:,15)
        qg1d     => spacebs1d(:,:,16)
        qs1d     => spacebs1d(:,:,17)
        resp1d   => spacebs1d(:,:,18)
        rhs1d    => spacebs1d(:,:,19)
        rno1d    => spacebs1d(:,:,20)
        rnos1d   => spacebs1d(:,:,21)
        rsw1d    => spacebs1d(:,:,22)
        sag1d    => spacebs1d(:,:,23)
        scv1d    => spacebs1d(:,:,24)
        sent1d   => spacebs1d(:,:,25)
        sice1d   => spacebs1d(:,:,26)
        ssw1d    => spacebs1d(:,:,27)
        t2m_1d   => spacebs1d(:,:,28)
        taf1d    => spacebs1d(:,:,29)
        tg1d     => spacebs1d(:,:,30)
        tgb1d    => spacebs1d(:,:,31)
        tlef1d   => spacebs1d(:,:,32)
        ts1d     => spacebs1d(:,:,33)
        tsw1d    => spacebs1d(:,:,34)
        u10m1d   => spacebs1d(:,:,35)
        v10m1d   => spacebs1d(:,:,36)
        veg1d    => spacebs1d(:,:,37)
        z1d      => spacebs1d(:,:,38)
        bfc      => spacebs1d(:,:,39)
        bsw      => spacebs1d(:,:,40)
        evmx0    => spacebs1d(:,:,41)
        fdry     => spacebs1d(:,:,42)
        fwet     => spacebs1d(:,:,43)
        gwmx0    => spacebs1d(:,:,44)
        gwmx1    => spacebs1d(:,:,45)
        gwmx2    => spacebs1d(:,:,46)
        porsl    => spacebs1d(:,:,47)
        relfc    => spacebs1d(:,:,48)
        rnet     => spacebs1d(:,:,49)
        texrat   => spacebs1d(:,:,50)
        vegt     => spacebs1d(:,:,51)
        wiltr    => spacebs1d(:,:,52)
        wt       => spacebs1d(:,:,53)
        xkmx     => spacebs1d(:,:,54)
        aarea    => spacebs1d(:,:,55)
        cdr      => spacebs1d(:,:,56)
        cdrn     => spacebs1d(:,:,57)
        cdrx     => spacebs1d(:,:,58)
        cf       => spacebs1d(:,:,59)
        cgrnd    => spacebs1d(:,:,60)
        cgrndl   => spacebs1d(:,:,61)
        cgrnds   => spacebs1d(:,:,62)
        clead    => spacebs1d(:,:,63)
        densi    => spacebs1d(:,:,64)
        efpr     => spacebs1d(:,:,65)
        eg       => spacebs1d(:,:,66)
        etr      => spacebs1d(:,:,67)
        etrrun   => spacebs1d(:,:,68)
        evaps    => spacebs1d(:,:,69)
        evapw    => spacebs1d(:,:,70)
        fevpg    => spacebs1d(:,:,71)
        flnet    => spacebs1d(:,:,72)
        flneto   => spacebs1d(:,:,73)
        fseng    => spacebs1d(:,:,74)
        htvp     => spacebs1d(:,:,75)
        ps       => spacebs1d(:,:,76)
        pw       => spacebs1d(:,:,77)
        qice     => spacebs1d(:,:,78)
        qsatl    => spacebs1d(:,:,79)
        rhosw    => spacebs1d(:,:,80)
        ribd     => spacebs1d(:,:,81)
        rlai     => spacebs1d(:,:,82)
        rpp      => spacebs1d(:,:,83)
        scrat    => spacebs1d(:,:,84)
        scvk     => spacebs1d(:,:,85)
        sdrop    => spacebs1d(:,:,86)
        seasb    => spacebs1d(:,:,87)
        sigf     => spacebs1d(:,:,88)
        sm       => spacebs1d(:,:,89)
        tm       => spacebs1d(:,:,90)
        uaf      => spacebs1d(:,:,91)
        vspda    => spacebs1d(:,:,92)
        wata     => spacebs1d(:,:,93)
        watr     => spacebs1d(:,:,94)
        watt     => spacebs1d(:,:,95)
        watu     => spacebs1d(:,:,96)
        wta      => spacebs1d(:,:,97)
        xlai     => spacebs1d(:,:,98)
        xlsai    => spacebs1d(:,:,99)
        xrun     => spacebs1d(:,:,100)
        cn1      => spacebs1d(:,:,101)
        df       => spacebs1d(:,:,102)
        rgr      => spacebs1d(:,:,103)
        wta0     => spacebs1d(:,:,104)
        wtaq0    => spacebs1d(:,:,105)
        wtg      => spacebs1d(:,:,106)
        wtg0     => spacebs1d(:,:,107)
        wtg2     => spacebs1d(:,:,108)
        wtga     => spacebs1d(:,:,109)
        wtgaq    => spacebs1d(:,:,110)
        wtgl     => spacebs1d(:,:,111)
        wtglq    => spacebs1d(:,:,112)
        wtgq     => spacebs1d(:,:,113)
        wtgq0    => spacebs1d(:,:,114)
        wtl0     => spacebs1d(:,:,115)
        wtlh     => spacebs1d(:,:,116)
        wtlq     => spacebs1d(:,:,117)
        wtlq0    => spacebs1d(:,:,118)
        wtshi    => spacebs1d(:,:,119)
        wtsqi    => spacebs1d(:,:,120)
        z2fra    => spacebs1d(:,:,121)
        z10fra   => spacebs1d(:,:,122)
        z1log    => spacebs1d(:,:,123)
        zlgocn   => spacebs1d(:,:,124)
        zlglnd   => spacebs1d(:,:,125)
        zlgsno   => spacebs1d(:,:,126)
        zlgveg   => spacebs1d(:,:,127)
        zlgdis   => spacebs1d(:,:,128)
        aldirs1d => spacebs1d(:,:,129)
        aldifs1d => spacebs1d(:,:,130)
        allocate(lveg(nnsg,iym1))
        lveg = -1
        allocate(oveg(nnsg,iym1))
        oveg = -1
        allocate(coszrs(iy))
        coszrs = d_zero
        allocate(spaceb1d(iym1,24))
        spaceb1d = d_zero
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

        allocate(fbat(njm2,iym2,numbat))
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

        allocate(fsub(nnsg,njm2,iym2,numsub))
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
