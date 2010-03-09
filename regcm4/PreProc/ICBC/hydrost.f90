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

      subroutine hydrost(h,t,phis,ps,pt,sigmaf,sigmah,dsigma,ni,nj,nk)
      use mod_constants , only : rovg
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj , nk
      real :: pt
      real , dimension(nk) :: dsigma , sigmah
      real , dimension(ni,nj,nk) :: h , t
      real , dimension(ni,nj) :: phis , ps
      real , dimension(nk+1) :: sigmaf
      intent (in) ni , nj , nk , phis , ps , pt , sigmaf , sigmah , t
      intent (inout) dsigma , h
!
! Local variables
!
      integer :: i , j , k , k1 , k2
      real :: pf , tbar
!
!     ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
!     THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
!     HEIGHT IS COMPUTED IN THE RCM MODEL.
!
      do k = 1 , nk
        dsigma(k) = sigmaf(k+1) - sigmaf(k)
      end do
!
!     SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
!     DEFINED AT CROSS POINTS AND AT HALF LEVELS.
!
      do i = 1 , ni
        do j = 1 , nj
          pf = pt/ps(i,j)
          h(i,j,nk) = phis(i,j) + rovg*t(i,j,nk)                        &
                    & *log((1.+pf)/(sigmah(nk)+pf))
          do k2 = 1 , nk - 1
            k = nk - k2
            k1 = k + 1
            tbar = (t(i,j,k)*dsigma(k)+t(i,j,k1)*dsigma(k1))            &
                 & /(dsigma(k)+dsigma(k1))
            h(i,j,k) = h(i,j,k1)                                        &
                     & + rovg*tbar*log((sigmah(k1)+pf)/(sigmah(k)+pf))
          end do
        end do
      end do
      end subroutine hydrost
