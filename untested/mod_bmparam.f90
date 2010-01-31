      module bmparam
      implicit none
      integer , parameter :: itb = 100 , jtb = 150
      real(8) :: pl , rdp , rdq , rdth , rdthe , thl
      real(8) , dimension(itb,jtb) :: ptbl
      real(8) , dimension(jtb) :: qs0 , sqs , sthe , the0
      real(8) , dimension(jtb,itb) :: ttbl
      real(8) :: aice , aliq , bice , bliq , cice1 , cliq , dice ,      &
               & dliq , xls0 , xls1
      end module bmparam
