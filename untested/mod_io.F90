      module io

      use regcm_param

      implicit none

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
! COMMON /BCVARSIO/
!
      real(8) , dimension(ix,mjx) :: ps0_io , ps1_io , ts0_io , ts1_io
      real(8) , dimension(ix,kx,mjx) :: qb0_io , qb1_io , so0_io ,      &
                                      & so1_io , tb0_io , tb1_io ,      &
                                      & ub0_io , ub1_io , vb0_io ,      &
                                      & vb1_io
!
! COMMON /BDYCODIO/
!
      real(8) , dimension(kx,mjx) :: ui1_io , ui2_io , uilx_io ,        &
                                   & uil_io , vi1_io , vi2_io ,         &
                                   & vilx_io , vil_io
!
! COMMON /IO/
!
      real(8) , dimension(ix,mjx,12,ntr) :: chemsrc_io
      real(8) , dimension(ix,mjx,ntr) :: ddsfc_io , dtrace_io ,         &
           & wdcvc_io , wdlsc_io
      real(8) , dimension(ilx,jxbb) :: pptc_io , pptnc_io , prca2d_io , &
                                     & prnca2d_io
!
! COMMON /MAINCMIO/
!
      real(8) , dimension(ix,kx,mjx,ntr) :: chia_io , chib_io
      real(8) , dimension(ix,mjx) :: cldefi_io , f_io , hfx_io ,        &
                                   & htsd_io , ht_io , msfd_io ,        &
                                   & msfx_io , psa_io , psb_io ,        &
                                   & qfx_io , rainc_io , rainnc_io ,    &
                                   & satbrt_io , tga_io , tgbb_io ,     &
                                   & tgb_io , uvdrag_io , xlat_io ,     &
                                   & xlong_io , zpbl_io
      real(8) , dimension(ix,kx,mjx) :: omega_io , qca_io , qcb_io ,    &
                                      & qva_io , qvb_io , ta_io ,       &
                                      & tbase_io , tb_io , ua_io ,      &
                                      & ub_io , va_io , vb_io
      real(8) , dimension(nnsg,ix,mjx) :: snowc_io
!
! COMMON /PMOISTIO/
!
      real(8) , dimension(ix,mjx) :: cbmf2d_io
      real(8) , dimension(ix,kx,mjx) :: fcc_io , rsheat_io , rswat_io
!
! COMMON /RAD1IO/
!
      real(8) , dimension(ilx,kx,mjx-1) :: heatrt_io
      real(8) , dimension(ilx,kxp1,mjx-1) :: o3prof_io
!
! COMMON /SPLITIO/
!
      real(8) , dimension(ix,mjx,nsplit) :: dstor_io , hstor_io
!
! COMMON /TRACEIO/
!
      real(8) , dimension(ix-1,kx,mjx-1) :: aerasp_io , aerext_io ,     &
           & aerssa_io
      real(8) , dimension(ix-1,mjx-1) :: aersrrf_io , aertarf_io
      real(8) , dimension(ix,mjx,ntr) :: cemtrac_io , cemtr_io ,        &
           & wxaq_io , wxsg_io
      real(8) , dimension(ix,mjx) :: dustsotex_io
      real(8) , dimension(ix,kx,mjx,ntr) :: rxsaq1_io , rxsaq2_io ,     &
           & rxsg_io
!
! COMMON /TRACHEMIO/
!
      real(8) , dimension(ix,kx,mjx,ntr) :: remcvc_io , remlsc_io
      real(8) , dimension(ix,mjx,ntr) :: remdrd_io

#endif

      end module io
