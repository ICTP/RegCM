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

      subroutine htsig(f3d,f2d,nta,nhga,npsa,zs,sig,pt,im,jm,km,im1,jm1,&
                     & n3d,n2d)
      use mod_constants , only : rgas , rgti
      implicit none
!
! Dummy arguments
!
      integer :: im , im1 , jm , jm1 , km , n2d , n3d , nhga , npsa ,   &
               & nta
      real(4) :: pt
      real(4) , dimension(im,jm,n2d) :: f2d
      real(4) , dimension(im,jm,km,n3d) :: f3d
      real(4) , dimension(km) :: sig
      real(4) , dimension(im,jm) :: zs
      intent (in) f2d , im , im1 , jm , jm1 , km , n2d , n3d , nhga ,   &
                & npsa , nta , pt , sig , zs
      intent (inout) f3d
!
! Local variables
!
      real(4) :: tbar
      integer :: i , j , k
      real(4) , dimension(im,jm) :: pstar
!
      pstar = f2d(:,:,npsa)
      do j = 1 , jm1
        do i = 1 , im1
          f3d(i,j,km,nhga) = zs(i,j) + rgas*rgti*f3d(i,j,km,nta)        &
                           & *log(pstar(i,j)                            &
                           & /((pstar(i,j)-pt)*sig(km)+pt))
        end do
      end do
      do k = km - 1 , 1 , -1
        do j = 1 , jm1
          do i = 1 , im1
            tbar = 0.5*(f3d(i,j,k,nta)+f3d(i,j,k+1,nta))
            f3d(i,j,k,nhga) = f3d(i,j,k+1,nhga)                         &
                            & + rgas*rgti*tbar*log(((pstar(i,j)-pt)     &
                            & *sig(k+1)+pt)/((pstar(i,j)-pt)*sig(k)+pt))
          end do
        end do
      end do
      end subroutine htsig
