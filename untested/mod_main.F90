      module mod_main

      use mod_regcm_param

      implicit none
!
! COMMON /MAINCM/
!
#ifdef MPP1
      real(8) , dimension(ix,jxp) :: cldefi , f , hfx , hgfact , htsd , &
                                   & qfx , rainc , rainnc , tgb , tgbb ,&
                                   & xlat , xlong , zpbl
      real(8) , dimension(ix,0:jxp+1) :: ht
      real(8) , dimension(ix,-1:jxp+2) :: msfd , msfx , pdotb , psa ,   &
           & rainp
      real(8) , dimension(ix,0:jxp+2) :: psb
      real(8) , dimension(ix,kx,-1:jxp+2) :: qca , qcb , qva , qvb ,    &
           & ta , tb , ua , ub , va , vb
      real(8) , dimension(ix,jxp+1) :: satbrt , tga
      real(8) , dimension(nnsg,ix,jxp) :: snowc
      real(8) , dimension(ix,kx,jxp) :: so4 , tbase
      real(8) , dimension(ix,0:jxp) :: uvdrag
#else
      real(8) , dimension(ix,jx) :: cldefi , f , hfx , hgfact , ht ,    &
                                  & htsd , msfd , msfx , pdotb , psa ,  &
                                  & psb , qfx , rainc , rainnc ,        &
                                  & satbrt , tga , tgb , tgbb , xlat ,  &
                                  & xlong , zpbl
      real(8) , dimension(ix,kx,jx) :: qca , qcb , qva , qvb , so4 ,    &
                                     & ta , tb , tbase , ua , ub , va , &
                                     & vb
      real(8) , dimension(nnsg,ix,jx) :: snowc
      real(8) , dimension(ix,jx) :: uvdrag
#endif

#ifdef MPP1
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
! COMMON /TMPSRF/
!
      real(8) , dimension(ix,nnsg*3+8,jxp) :: inisrf0
      real(8) , dimension(ix,nnsg*3+8,mjx) :: inisrf_0
!
! COMMON /TMP_var1/
!
      real(8) , dimension(kx,8) :: var1snd , var1rcv
!
#endif

      end module mod_main
