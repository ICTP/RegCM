      subroutine intgtb(pa,za,tlayer,zrcm,tp,zp,sccm,ni,nj,nlev1)
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj , nlev1
      real , dimension(ni,nj) :: pa , tlayer , za , zrcm
      real , dimension(nlev1) :: sccm
      real , dimension(ni,nj,nlev1) :: tp , zp
      intent (in) ni , nj , nlev1 , sccm , tp , zp , zrcm
      intent (out) pa , za
      intent (inout) tlayer
!
! Local variables
!
      integer :: i , j , k , kb , kt
!
!     INTGTB CALCULATES ALL VARIABLES NEEDED TO COMPUTE P* ON THE RCM
!     TOPOGRAPHY.  THE MEAN TEMPERATURE IN THE LAYER BETWEEN
!     THE TOPOGRAPHY AND THE PRESSURE LEVEL ABOVE IS CALULATED
!     BY LINEARLY INTERPOLATING WITH HEIGHT THE TEMPS ON
!     PRESSURE LEVELS.
!     INPUT:    TP        TEMPS ON ECMWF PRESSURE LEVELS
!     ZP        HEIGHTS OF ECMWF PRESSURE LEVELS
!     ZRCM      RCM TOPOGRAPHY
!     SCCM      ECMWF PRESSURE LEVELS (DIVIDED BY 1000.)
!     OUTPUT:   TLAYER    MEAN LAYER TEMP ABOVE RCM SURFACE
!     PA        PRESSURE AT TOP OF LAYER
!     ZA        HEIGHT AT PRESSURE PA
!
      do i = 1 , ni
        do j = 1 , nj
 
          kt = 0
          do k = 1 , nlev1 - 1
            if ( zrcm(i,j)<=zp(i,j,nlev1+1-k) .and. zrcm(i,j)           &
               & >zp(i,j,nlev1-k) ) kt = k
          end do
          kb = kt + 1
 
          if ( kt/=0 ) then
            tlayer(i,j) = (tp(i,j,nlev1+1-kt)*(zrcm(i,j)-zp(i,j,nlev1+1-&
                        & kb))+tp(i,j,nlev1+1-kb)                       &
                        & *(zp(i,j,nlev1+1-kt)-zrcm(i,j)))              &
                        & /(zp(i,j,nlev1+1-kt)-zp(i,j,nlev1+1-kb))
            tlayer(i,j) = (tp(i,j,nlev1+1-kt)+tlayer(i,j))/2.
            za(i,j) = zp(i,j,nlev1+1-kt)
            pa(i,j) = 100.*sccm(kt)
          else
            tlayer(i,j) = tp(i,j,1)
            za(i,j) = zp(i,j,1)
            pa(i,j) = 100.
          end if
 
        end do
      end do
 
!     PRINT *, 'ZRCM, ZP(6)   =', ZRCM(5,5), ZP(5,5,NLEV1+1-6)
!     PRINT *, '      TP(6)   =',            TP(5,5,NLEV1+1-6)
!     PRINT *, 'TLAYER, ZA, PA =', TLAYER(5,5), ZA(5,5), PA(5,5)
 
      end subroutine intgtb
