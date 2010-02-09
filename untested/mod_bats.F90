      module bats

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
      integer , parameter :: numbat = 21 + 6
#ifdef MPP1
      real(kind=4) , dimension(mjx-2,ix-2,numbat) :: fbat_io
      real(kind=4) , dimension(jxp,ix-2,numbat) :: fbat
#else
      real(kind=4) , dimension(jx-2,ix-2,numbat) :: fbat
#endif

#ifdef MPP1
!
! COMMON /BAT2D1IO/
!
      real(8) , dimension(nnsg,ilx,jxbb) :: col2d_io , dew2d_io ,       &
           & evpa2d_io , gwet2d_io , ircp2d_io , ocld2d_io , rno2d_io , &
           & rnos2d_io , sag2d_io , scv2d_io , sena2d_io , sice2d_io ,  &
           & srw2d_io , ssw2d_io , swt2d_io , taf2d_io , text2d_io ,    &
           & tg2d_io , tgb2d_io , tlef2d_io , veg2d1_io
      real(8) , dimension(nnsg,ix,mjx) :: ht1_io , satbrt1_io
!
! COMMON /BAT2DIO/
!
      real(8) , dimension(ilx,jxbb) :: flw2d_io , flwd2d_io , fsw2d_io ,&
                                     & sabv2d_io , sdelqk2d_io ,        &
                                     & sdeltk2d_io , sfracb2d_io ,      &
                                     & sfracs2d_io , sfracv2d_io ,      &
                                     & sinc2d_io , sol2d_io ,           &
                                     & solvd2d_io , solvs2d_io ,        &
                                     & ssw2da_io , svegfrac2d_io ,      &
                                     & veg2d_io
!
! COMMON /SUBOUTIO/
!
      integer , parameter :: numsub = 16
      real(4) , dimension(nnsg,mjx-2,ix-2,numsub) :: fsub_io
      real(kind=4)  fsub(NNSG,jxp,ix-2,numsub)
      equivalence (fsub(1,1,1,1),u10m_s(1,1,1))

#endif

      end module bats
