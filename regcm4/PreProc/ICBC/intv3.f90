!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine intv3(fsccm,fccm,psrccm,sccm,ptop,ni,nj,kccm)
      use mod_constants , only : rgas , gti , lrate
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: rgas2 = rgas/2.
      real , parameter :: b1 = -gti/lrate
!
! Dummy arguments
!
      integer :: kccm , ni , nj
      real :: ptop
      real , dimension(ni,nj,kccm) :: fccm
      real , dimension(ni,nj) :: fsccm , psrccm
      real , dimension(kccm) :: sccm
      intent (in) fccm , kccm , ni , nj , psrccm , ptop , sccm
      intent (out) fsccm
!
! Local variables
!
      real :: a1 , rc , rc1 , sc
      integer :: i , j , k , k1 , kp1
!
!**   INTV3 IS FOR VERTICAL INTERPOLATION OF TSCCM.  THE INTERPOLATION
!     IS LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATUES IN THE LAYER.
!
      do i = 1 , ni
        do j = 1 , nj
          sc = (psrccm(i,j)+ptop)/100.
          do k = 1 , kccm - 1
            if ( sc<=sccm(k+1) .and. sc>=sccm(k) ) k1 = k
          end do
 
!
!         CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
!
          if ( sc>sccm(kccm) ) then
            a1 = rgas2*log(sc/sccm(kccm))
            fsccm(i,j) = fccm(i,j,kccm+1-kccm)*(b1-a1)/(b1+a1)
          else
!
!         CONDITION FOR SC .LT. SCCM(KCCM) FOLLOWS
!
            kp1 = k1 + 1
            rc = log(sc/sccm(k1))/log(sccm(k1)/sccm(kp1))
            rc1 = rc + 1.
            fsccm(i,j) = rc1 * fccm(i,j,kccm+1-k1) -                    &
                       & rc  * fccm(i,j,kccm+1-kp1)
          end if
!
        end do
      end do
 
      end subroutine intv3
