      module mod_radbuf
      use mod_regcm_param
      implicit none
      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity
#ifdef MPP1
      real(8) , dimension(ixm1,kx,4,jxp) :: absnxt
      real(8) , dimension(ixm1,kxp1,kx + 1,jxp) :: abstot
      real(8) , dimension(ixm1,kxp1,jxp) :: emstot
#else
      real(8) , dimension(ixm1,kx,4,jxm1) :: absnxt
      real(8) , dimension(ixm1,kxp1,kx + 1,jxm1) :: abstot
      real(8) , dimension(ixm1,kxp1,jxm1) :: emstot
#endif
      end module mod_radbuf

