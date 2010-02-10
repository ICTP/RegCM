      module mod_bats

      use regcm_param

      implicit none
!
! COMMON /BATCNT/
!
      integer :: npts
!
! COMMON /BATS1D/
!
      real(8) , dimension(nnsg,nbmax) :: delq1d , delt1d , drag1d ,     &
           & emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d , ldoc1d ,    &
           & p1d , pbp1d , prcp1d , q2m_1d , qg1d , qs1d , resp1d ,     &
           & rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
           & sice1d , ssw1d , t2m_1d , taf1d , tg1d , tgb1d , tlef1d ,  &
           & ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d
      real(8) , dimension(nbmax) :: flw1d , fsw1d , us1d , vs1d
!
! COMMON /BATVAR/
!
      real(8) , dimension(nnsg,nbmax) :: bfc , bsw , evmx0 , fdry ,     &
           & fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc , rnet ,      &
           & texrat , vegt , wiltr , wt , xkmx
      real(8) , dimension(nbmax) :: czen , sola , vpdd
      real(8) :: difrat
!
! COMMON /BDCN/
!
      real(8) , dimension(20) :: albvgl , albvgs , crough , deprv ,     &
                               & deptv , depuv , displa , fc , freza ,  &
                               & frezu , rough , rsmin , sai , seasf ,  &
                               & sqrtdi , vegc , xla , xlai0
      real(8) , dimension(12) :: bee , skrat , xmofc , xmohyd , xmopor ,&
                               & xmosuc , xmowil
      real(8) , dimension(130) :: c
      real(8) :: ch2o , cice , csnw , csoilc , cwi , cws , dewmx ,      &
               & dewmxi , drain , pi , rmax0 , tau1 , trsmx0 , vonkar , &
               & zlnd , zoce , zsno
      integer , dimension(20) :: iexsol , kolsol
      integer :: ihis , lat , mhis , ncase
!
! COMMON /BDPRM/
!
      real(8) , dimension(nnsg,nbmax) :: aarea , cdr , cdrn , cdrx ,    &
           & cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg , &
           & etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,    &
           & fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,     &
           & rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,    &
           & tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,&
           & xlsai , xrun , z1 , z1log
      real(8) , dimension(nbmax) :: ems
      integer , dimension(nnsg,nbmax) :: imelt , lveg
!
! COMMON /LEAFWT/
!
      real(8) , dimension(nnsg,nbmax) :: cn1 , df , rgr , wta0 , wtaq0 ,&
           & wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq , wtgq ,   &
           & wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi

!
! COMMON /BAT1D/
!
      real(8) , dimension(ilx) :: albdif , albdir , albvl , albvld ,    &
                                & albvs , albvsd , aldifl , aldifs ,    &
                                & aldirl , aldirs , emiss1d , fracd ,   &
                                & sabveg , solis , solvd , solvs
      real(8) , dimension(ix) :: coszrs
!
! COMMON /BAT2D/
!
#ifdef MPP1
      real(8) , dimension(ilx,jxp) :: flw2d , flwa2d , flwd2d ,         &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
#else
      real(8) , dimension(ilx,jlx) :: flw2d , flwa2d , flwd2d ,         &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
#endif
!
! COMMON /BAT2D1/
!
#ifdef MPP1
      real(8) , dimension(nnsg,ilx,jxp) :: col2d , dew2d , emiss2d ,    &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,ix,jxp) :: ht1 , satbrt1
#else
      real(8) , dimension(nnsg,ilx,jlx) :: col2d , dew2d , emiss2d ,    &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,ix,jx) :: ht1 , satbrt1
#endif
!
! COMMON /BATOUT/
!
#ifdef MPP1
      real(4) , dimension(jxp,ix-2) :: drag_o , evpa_o , flwa_o ,       &
                                     & flwd_o , fswa_o , prcv_o ,       &
                                     & psmn_o , ps_o , q2m_o , rnos_o , &
                                     & rsw_o , scv_o , sena_o , sina_o ,&
                                     & ssw_o , t2mn_o , t2mx_o , t2m_o ,&
                                     & tgmn_o , tgmx_o , tg_o , tlef_o ,&
                                     & tpr_o , u10m_o , v10m_o ,        &
                                     & w10x_o , zpbl_o
#else
      real(4) , dimension(jx-2,ix-2) :: drag_o , evpa_o , flwa_o ,      &
                                      & flwd_o , fswa_o , prcv_o ,      &
                                      & psmn_o , ps_o , q2m_o , rnos_o ,&
                                      & rsw_o , scv_o , sena_o ,        &
                                      & sina_o , ssw_o , t2mn_o ,       &
                                      & t2mx_o , t2m_o , tgmn_o ,       &
                                      & tgmx_o , tg_o , tlef_o , tpr_o ,&
                                      & u10m_o , v10m_o , w10x_o ,      &
                                      & zpbl_o
#endif
!
! COMMON /BATUNIT/
!
      integer :: nrcbat , nrcsub
!
! COMMON /SUBOUT/
!
#ifdef MPP1
      real(4) , dimension(nnsg,jxp,ix-2) :: drag_s , evpa_s , prcv_s ,  &
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#else
      real(4) , dimension(nnsg,jx-2,ix-2) :: drag_s , evpa_s , prcv_s , &
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#endif
!
! COMMON BATOUTio
!
#ifdef MPP1
      real(kind=4) , dimension(jxp,ix-2,numbat) :: fbat
#else
      real(kind=4) , dimension(jx-2,ix-2,numbat) :: fbat
#endif

#ifdef MPP1
      real(kind=4) , dimension(NNSG,jxp,ix-2,numsub) :: fsub
#else
      real(kind=4) , dimension(NNSG,jx-2,ix-2,numsub) :: fsub
#endif
!
! COMMON /BATS1D0/
!
      real(8) , dimension(nnsg,nbmax) :: p1d0 , qs1d0 , ts1d0
!
!*    vegc is maximum fractional cover of vegetation
      data vegc/.85 , .8 , .8 , .8 , .8 , .9 , .8 , 0.0 , 0.6 , 0.8 ,   &
         & 0.35 , 0.0 , .8 , 2*0. , 5*0.8/
!*    seasf is the difference between vegc and fractional cover at 269k
      data seasf/0.6 , 0.1 , 0.1 , 0.3 , 0.3 , 0.5 , 0.3 , 0. , 0.2 ,   &
         & 0.6 , 0.1 , 0.0 , 0.4 , 2*0.0 , .2 , .3 , .2 , 2*.4/
!*    rough is an aerodynamic roughness length (m) =approx 0.1*veg
!*    height also used snow masking depth in subrout albedo
      data rough/0.08 , 0.05 , 2*1.0 , 0.8 , 2.0 , 0.1 , 0.05 , 0.04 ,  &
         & 0.06 , 0.1 , 0.01 , 0.03 , 2*0.0004 , 2*0.1 , 0.8 , 2*0.3/
!     ******      displacement height (meter)
!     ******      if great parts of veg. are covered by snow, use
!     displa=0 ******      because then the new displa-theory is not
!     valid
      data displa/0. , 0. , 9. , 9. , 0. , 18. , 14*0./
!     ******      min stomatl resistance (s/m)
!cc   data rsmin/153.0,4*200.0,150.0,14*200.0/   ! shuttleworth
!     data rsmin/120.0,4*200.0,150.0,14*200.0/   ! bats1e numbers
!Sara
!     data rsmin/45.,60.,2*80.,120.,2*60.,200.,80.,45.,150.,200.,45.
!     &          ,2*200.,80.,120.,100.,2*120./
      data rsmin/20*200./
!Sara_
!     ******      max leaf area index (ratio unit cover per unit ground)
!     ORIGINAL
!     data xla/6.,2.,5*6.,0.,3*6.,0.,6.,2*0.,5*6./
!     Laura 21/04/08
      data xla/4. , 2. , 4*6. , 3. , 0. , 2. , 4. , 1. , 0. , 4. ,      &
         & 2*0. , 4. , 4. , 5. , 4. , 1./
!     ******      min leaf area index **lai depends on temp as veg cover
!     ORIGINAL
!     data xlai0/0.5,0.5,5.,2*1.,5.,0.5,0.,3*0.5,0.,0.5,2*0.,5.,1.,3.,
!     &2*0.5/
!     Laura 21/04/08
      data xlai0/2*0.5 , 5. , 2*1. , 5. , 1 , 0. , 0.5 , 2. , 0.5 , 0. ,&
         & 2. , 2*0. , 3. , 1. , 3. , 0.5 , 1.0/

!     ******      stem area index (projected area of non-transpiring
!     sfcs)
      data sai/.5 , 4. , 5*2. , 2*.5 , 11*2./
!     ******      inverse square root of leaf dimension - used for
!     ******      calculating fluxes from foliage
      data sqrtdi/10. , 19*5.0/
!     ******      fc = light dependence of stomatal resistance
      data fc/.02 , .02 , 4*.06 , 11*.02 , .06 , 2*.02/

!     ******      depuv is depth of upper soil layer (mm)
!     ******      deprv is depth of root zone (mm)
!     ******      deptv is depth of total soil (mm)
      data depuv/20*100./
      data deprv/2*1000. , 2*1500. , 2000. , 1500. , 11*1000. , 2000. , &
         & 2*2000./
      data deptv/20*3000./
!     ******      iexsol is soil texture type (see subr soilbc)
!     ORIGINAL
!     data iexsol/6,6,6,6,7,8,6,3,6,6,5,12,6,6,6,6,5,6,6,6/
!     Laura  04/04/08 changed soil texture for desert: 3->1
      data iexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12 , 6 ,  &
         & 6 , 6 , 6 , 5 , 6 , 6 , 6/
!     ******      kolsol is soil color type (see subr. albedo)
!     Dec. 15, 2008
!     data kolsol/5,3,4,4,4,4,4,1,3,3,2,1,5,5,5,4,3,4,4,4/
      data kolsol/6 , 4 , 5 , 5 , 5 , 5 , 5 , 1 , 4 , 4 , 2 , 1 , 6 ,   &
         & 6 , 6 , 5 , 4 , 5 , 5 , 5/
!     Dec. 15, 2008_
!     ******      xmopor is fraction of soil that is voids
      data xmopor/.33 , .36 , .39 , .42 , .45 , .48 , .51 , .54 , .57 , &
         & .6 , .63 , .66/
!     ******      xmosuc is the minimum soil suction (mm)
      data xmosuc/3*30.0 , 9*200./
!     ******      xmohyd is the max. hydraulic conductivity (mm/s)
      data xmohyd/0.20E-0 , 0.80E-1 , 0.32E-1 , 0.13E-1 , 0.89E-2 ,     &
         & 0.63E-2 , 0.45E-2 , 0.32E-2 , 0.22E-2 , 0.16E-2 , 0.11E-2 ,  &
         & 0.80E-3/
!     ******      xmowilt is fraction of water content at which
!     permanent wilting occurs
      data xmowil/.095 , .128 , .161 , .266 , .3 , .332 , .378 , .419 , &
         & .455 , .487 , .516 , .542/
      data xmofc/.404 , .477 , .547 , .614 , .653 , .688 , .728 , .763 ,&
         & .794 , .820 , .845 , .866/
!     ******      bee is the clapp and hornbereger "b" parameter
      data bee/3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 6.8 , 7.6 , 8.4 ,    &
         & 9.2 , 10.0 , 10.8/
!     ******      bskrat is ratio of soil thermal conduc. to that of
!     loam - a function of texture
      data skrat/1.7 , 1.5 , 1.3 , 1.2 , 1.1 , 1.0 , .95 , .90 , .85 ,  &
         & .80 , .75 , .7/

!     Dec. 15, 2008
!     ******      albvgs is vegetation albedo for wavelengths < 0.7
!     microns data
!     albvgs/.1,.1,.05,.05,.08,.04,.08,.2,.1,.08,.17,.8,.06,2*.07,
!     &.05,.08,.06,2*0.06/
      data albvgs/.1 , .1 , .04 , .04 , .06 , .04 , .08 , .2 , .1 ,     &
         & .08 , .17 , .8 , .06 , 2*.07 , .05 , .08 , .05 , 2*0.06/
!     ******      albvgl is vegetation albedo for wavelengths > 0.7
!     microns data
!     albvgl/.3,.3,.23,.23,.28,.20,.30,.4,.3,.28,.34,.6,.18,2*.2,
!     &.23,.28,.24,2*.18/
      data albvgl/.3 , .3 , .20 , .20 , .26 , .20 , .30 , .4 , .3 ,     &
         & .28 , .34 , .6 , .18 , 2*.2 , .23 , .28 , .23 , 2*.18/
!     Dec. 15, 2008_

      contains
        subroutine fillbat
          implicit none
          fbat(:,:,1)  = u10m_o(:,:)
          fbat(:,:,2)  = v10m_o(:,:)
          fbat(:,:,3)  = drag_o(:,:)
          fbat(:,:,4)  = tg_o(:,:)
          fbat(:,:,5)  = tlef_o(:,:)
          fbat(:,:,6)  = t2m_o(:,:)
          fbat(:,:,7)  = q2m_o(:,:)
          fbat(:,:,8)  = ssw_o(:,:)
          fbat(:,:,9)  = rsw_o(:,:)
          fbat(:,:,10) = tpr_o(:,:)
          fbat(:,:,11) = evpa_o(:,:)
          fbat(:,:,12) = rnos_o(:,:)
          fbat(:,:,13) = scv_o(:,:)
          fbat(:,:,14) = sena_o(:,:)
          fbat(:,:,15) = flwa_o(:,:)
          fbat(:,:,16) = fswa_o(:,:)
          fbat(:,:,17) = flwd_o(:,:)
          fbat(:,:,18) = sina_o(:,:)
          fbat(:,:,19) = prcv_o(:,:)
          fbat(:,:,20) = ps_o(:,:)
          fbat(:,:,21) = zpbl_o(:,:)
          fbat(:,:,22) = tgmx_o(:,:)
          fbat(:,:,23) = tgmn_o(:,:)
          fbat(:,:,24) = t2mx_o(:,:)
          fbat(:,:,25) = t2mn_o(:,:)
          fbat(:,:,26) = w10x_o(:,:)
          fbat(:,:,27) = psmn_o(:,:)
        end subroutine fillbat

        subroutine fillsub
        implicit none
          fsub(:,:,:,1)  = u10m_s(:,:,:)
          fsub(:,:,:,2)  = v10m_s(:,:,:)
          fsub(:,:,:,3)  = drag_s(:,:,:)
          fsub(:,:,:,4)  = tg_s(:,:,:)
          fsub(:,:,:,5)  = tlef_s(:,:,:)
          fsub(:,:,:,6)  = t2m_s(:,:,:)
          fsub(:,:,:,7)  = q2m_s(:,:,:)
          fsub(:,:,:,8)  = ssw_s(:,:,:)
          fsub(:,:,:,9)  = rsw_s(:,:,:)
          fsub(:,:,:,10) = tpr_s(:,:,:)
          fsub(:,:,:,11) = prcv_s(:,:,:)
          fsub(:,:,:,12) = evpa_s(:,:,:)
          fsub(:,:,:,13) = rnos_s(:,:,:)
          fsub(:,:,:,14) = scv_s(:,:,:)
          fsub(:,:,:,15) = sena_s(:,:,:)
          fsub(:,:,:,16) = ps_s(:,:,:)
        end subroutine fillsub

      end module mod_bats
