      subroutine intv1(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: psccm = 100.
!
! Dummy arguments
!
      integer :: kccm , krcm , ni , nj
      real :: pt
      real , dimension(ni,nj,kccm) :: fccm
      real , dimension(ni,nj,krcm) :: frcm
      real , dimension(ni,nj) :: psrcm
      real , dimension(kccm) :: sccm
      real , dimension(krcm) :: srcm
      intent (in) fccm , kccm , krcm , ni , nj , psrcm , pt , sccm ,    &
                & srcm
      intent (out) frcm
!
! Local variables
!
      real :: dp1 , pt1 , rc , rc1 , sc
      integer :: i , j , k , k1 , k1p , n
!
!     INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!     HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!     IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!     LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATUES IN THE LAYER.
! 
      do i = 1 , ni
        do j = 1 , nj
          dp1 = psrcm(i,j)/psccm
          pt1 = pt/psccm
          do n = 1 , krcm
            sc = srcm(n)*dp1 + pt1
            k1 = 0
            do k = 1 , kccm
              if ( sc>sccm(k) ) k1 = k
            end do
!
!           CONDITION FOR SC .LT. SCCM(1) FOLLOWS
!
            if ( k1==0 ) then
              frcm(i,j,n) = fccm(i,j,kccm)
!
!             CONDITION FOR SCCM(1) .LT. SC .LT. SCCM(KCCM) FOLLOWS
!
            else if ( k1/=kccm ) then
              k1p = k1 + 1
              rc = (sc-sccm(k1))/(sccm(k1)-sccm(k1p))
              rc1 = rc + 1.
              frcm(i,j,n) = rc1*fccm(i,j,kccm+1-k1)                     &
                          & - rc*fccm(i,j,kccm+1-k1p)
!
!             CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
!
            else
              frcm(i,j,n) = fccm(i,j,1)
!
            end if
          end do
        end do
      end do
 
      end subroutine intv1
