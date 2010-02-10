      module mod_param1

      use mod_regcm_param

      implicit none
!
! COMMON /I8PARAM1/
!
      integer :: jyear , jyear0 , jyearr , ktau , ktaur , ntime
!
! COMMON /IPARAM1/
!
      integer :: ibdyfrq , ifrabe , klake , nbatst , nslice
!
! COMMON /PARAM1/
!
      real(8) :: abatm , abemh , c200 , c201 , c203 , dt , dt0 , dt2 ,  &
               & dtbat , dtlake , dtmin , dx , dx16 , dx2 , dx4 , dx8 , &
               & dxsq , fnudge , gnudge , xkhmax , xkhz , xtime
      real(8) , dimension(nsplit) :: dtau
      end module mod_param1
