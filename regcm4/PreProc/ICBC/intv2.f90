!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine intv2(frcm,fccm,psrcm,srcm,sccm,pt,ni,nj,krcm,kccm)
      use mod_constants , only : rgas , gti , lrate
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: rgas2 = rgas/2.
      real , parameter :: b1 = -gti/lrate
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
      real :: a1 , dp1 , pt1 , rc , rc1 , sc
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
              rc = log(sc/sccm(k1))/log(sccm(k1)/sccm(k1p))
              rc1 = rc + 1.
              frcm(i,j,n) = rc1*fccm(i,j,kccm+1-k1)                     &
                          & - rc*fccm(i,j,kccm+1-k1p)
!
!             CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
!
            else
              a1 = rgas2*log(sc/sccm(kccm))
              frcm(i,j,n) = fccm(i,j,1)*(b1-a1)/(b1+a1)
!
            end if
          end do
        end do
      end do
 
      end subroutine intv2
