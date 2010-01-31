      module param2
      implicit none
!
! COMMON /IPARAM2/
!
      integer :: ibintyp , ibltyp , iboudy , ichem , icnt , icup ,      &
               & idirect , iemiss , igcc , iocnflx , iotyp , ipgf ,     &
               & ipptls , kbats , kchem , lakemod , maschk , nradisp ,  &
               & ntrad , ntsave , nttape
!
! COMMON /LPARAM2/
!
      logical :: ifbat , ifchem , ifprt , ifrad , ifrest , ifsave ,     &
               & ifsub , iftape , rfstrt
!
! COMMON /PARAM2/
!
      real(8) :: batfrq , bdytim , chemfrq , prtfrq , prttim , radfrq , &
               & radisp , savfrq , savtim , tapfrq , taptim , tbdybe
      integer :: nprtfrq , nsavfrq , ntapfrq
      end module param2
