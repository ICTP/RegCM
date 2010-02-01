      module o3blk
      use regcm_param
      implicit none
      real(8) , dimension(31) :: o3ann , o3sum , o3win , o3wrk , ppann ,&
                               & ppsum , ppwin , ppwrk
      real(8) , dimension(32) :: ppwrkh
      real(8) , dimension(kxp2) :: prlevh
      end module o3blk
