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

      subroutine calcmse(f2d,f3d,nx,ny,kz,nfld2d,nfld3d,zs,sigh,pt,nta, &
                       & nqva,npsa,nmse,nx1,ny1)
 
      use mod_constants , only : gti , rovg , wlhv , cpd
      implicit none
!
! Dummy arguments
!
      integer :: nfld2d , nfld3d , nmse , npsa , nqva , nta , nx , nx1 ,&
               & ny , ny1 , kz
      real(4) :: pt
      real(4) , dimension(nx,ny,nfld2d) :: f2d
      real(4) , dimension(nx,ny,kz,nfld3d) :: f3d
      real(4) , dimension(kz) :: sigh
      real(4) , dimension(nx,ny) :: zs
      intent (in) f2d , nfld2d , nfld3d , nmse , npsa , nqva , nta ,    &
                & nx , nx1 , ny , ny1 , kz , pt , sigh , zs
      intent (inout) f3d
!
! Local variables
!
      real(4) :: dp , ps , q , t , tvbar
      integer :: i , j , k , k1
      real(4) , dimension(nx,ny,kz) :: p , tv , z
!
      do j = 1 , ny1
        do i = 1 , nx1
          t = f3d(i,j,kz,nta)
          q = f3d(i,j,kz,nqva)
          ps = f2d(i,j,npsa)
          tv(i,j,kz) = t*(1.+0.608*q)
          p(i,j,kz) = (ps-pt*10.)*sigh(kz) + pt*10.
          dp = log(p(i,j,kz)) - log(ps)
          z(i,j,kz) = zs(i,j) - dp*tv(i,j,kz)*rovg
          f3d(i,j,kz,nmse) = gti*z(i,j,kz) + cpd*t + wlhv*q
        end do
      end do
      do k = kz - 1 , 1 , -1
        do j = 1 , ny1
          do i = 1 , nx1
            k1 = k + 1
            t = f3d(i,j,k,nta)
            q = f3d(i,j,k,nqva)
            ps = f2d(i,j,npsa)
            tv(i,j,k) = t*(1.+0.608*q)
            tvbar = 0.5*(tv(i,j,k)+tv(i,j,k1))
            p(i,j,k) = (ps-pt*10.)*sigh(k) + pt*10.
            dp = log(p(i,j,k)) - log(p(i,j,k1))
            z(i,j,k) = z(i,j,k1) - dp*tvbar*rovg
            f3d(i,j,k,nmse) = gti*z(i,j,k) + cpd*t + wlhv*q
          end do
        end do
      end do
 
      end subroutine calcmse
