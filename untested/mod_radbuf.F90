      module radbuf
      use parrad
      implicit none
      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity
#ifdef MPP1
      real(8) , dimension(plond,plev,4,jxp) :: absnxt
      real(8) , dimension(plond,plevp,plevp,jxp) :: abstot
      real(8) , dimension(plond,plevp,jxp) :: emstot
      real(8) , dimension(plond,plev,4,mjx-1) :: absnxt_io
      real(8) , dimension(plond,plevp,plevp,mjx-1) :: abstot_io
      real(8) , dimension(plond,plevp,mjx-1) :: emstot_io
#else
      real(8) , dimension(plond,plev,4,jlx) :: absnxt
      real(8) , dimension(plond,plevp,plevp,jlx) :: abstot
      real(8) , dimension(plond,plevp,jlx) :: emstot
#endif
      end module radbuf

