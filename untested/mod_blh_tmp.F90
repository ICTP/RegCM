      module blh_tmp
      use regcm_param
      implicit none
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: cgh , kvc , kvh , kvm , kvq
      real(8) , dimension(ix,jxp) :: hfxv , obklen , th10 , ustr ,      &
                                   & xhfx , xqfx
#else
      real(8) , dimension(ix,kx,jlx) :: cgh , kvc , kvh , kvm , kvq
      real(8) , dimension(ix,jx) :: hfxv , obklen , th10 , ustr ,      &
                                   & xhfx , xqfx
#endif
      end module blh_tmp
