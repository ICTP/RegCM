      module mod_radbuf
      use mod_regcm_param
      implicit none
      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity
#ifdef MPP1
      real(8) , dimension(ix - 1,kx,4,jxp) :: absnxt
      real(8) , dimension(ix - 1,kx + 1,kx + 1,jxp) :: abstot
      real(8) , dimension(ix - 1,kx + 1,jxp) :: emstot
      real(8) , dimension(ix - 1,kx,4,mjx-1) :: absnxt_io
      real(8) , dimension(ix - 1,kx + 1,kx + 1,mjx-1) :: abstot_io
      real(8) , dimension(ix - 1,kx + 1,mjx-1) :: emstot_io
#else
      real(8) , dimension(ix - 1,kx,4,jlx) :: absnxt
      real(8) , dimension(ix - 1,kx + 1,kx + 1,jlx) :: abstot
      real(8) , dimension(ix - 1,kx + 1,jlx) :: emstot
#endif
      end module mod_radbuf

