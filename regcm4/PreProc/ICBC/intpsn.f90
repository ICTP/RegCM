      subroutine intpsn(psrcm,zrcm,pa,za,tlayer,pt,ni,nj)
      use mod_constants , only : govr
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj
      real :: pt
      real , dimension(ni,nj) :: pa , psrcm , tlayer , za , zrcm
      intent (in) ni , nj , pa , pt , tlayer , za , zrcm
      intent (out) psrcm
!
! Local variables
!
      real :: tb
      integer :: i , j
!
!     EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
!     USE TLAYER CALCULATED IN INTGTB.
!     PSRCM = SURFACE PRESSURE - PTOP
!
      do i = 1 , ni
        do j = 1 , nj
          tb = tlayer(i,j)
          psrcm(i,j) = pa(i,j)*exp(-govr*(zrcm(i,j)-za(i,j))/tb) - pt
        end do
      end do
 
!     PRINT *, 'ZRCM, ZA, PA, PT =', ZRCM(5,5), ZA(5,5), PA(5,5), PT
!     PRINT *, 'TLAYER(5,5), PSRCM(5,5) = ', TLAYER(5,5), PSRCM(5,5)
 
      end subroutine intpsn
