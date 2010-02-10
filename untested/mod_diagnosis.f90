      module mod_diagnosis

      use mod_regcm_param

      implicit none
!
! COMMON /DIAG_ATM/
!
      real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
! COMMON /TRACMASS/
!
      real(8) , dimension(ntr) :: tchiad , tchie , tchitb
      real(8) , dimension(ntr,2) :: tremcvc , tremdrd , tremlsc ,       &
                                  & trxsaq1 , trxsaq2 , trxsg , ttrace
      end module mod_diagnosis
