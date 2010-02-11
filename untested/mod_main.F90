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

      end module mod_main
