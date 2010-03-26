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

      subroutine intlog(fp,f,f2,npsa,pt,sig,im,jm,km,nv,p,kp,n3d,n2d,   &
                      & im1,jm1)

      use mod_constants , only : rgas , rgti , lrate , bltop
      implicit none
!
! Dummy arguments
!
      integer :: im , im1 , jm , jm1 , km , kp , n2d , n3d , npsa , nv
      real(4) :: pt
      real(4) , dimension(im,jm,km,n3d) :: f
      real(4) , dimension(im,jm,n2d) :: f2
      real(4) , dimension(im,jm,kp,n3d) :: fp
      real(4) , dimension(kp) :: p
      real(4) , dimension(km) :: sig
      intent (in) f , f2 , im , im1 , jm , jm1 , km , kp , n2d , n3d ,  &
                & npsa , nv , p , pt , sig
      intent (out) fp
!
! Local variables
!
      real(4) :: sigp , w1 , wp
      integer :: i , j , k , k1 , k1p , kbc , n
      real(4) , dimension(im,jm) :: pstar
!
      pstar = f2(:,:,npsa)
!
      do k = 1 , km
        if ( sig(k)<bltop ) kbc = k
      end do
      do j = 1 , jm1
        do i = 1 , im1
          do n = 1 , kp
            sigp = (p(n)-pt)/(pstar(i,j)-pt)
            k1 = 0
            do k = 1 , km
              if ( sigp>sig(k) ) k1 = k
            end do
            if ( sigp<=sig(1) ) then
              fp(i,j,n,nv) = f(i,j,1,nv)
            else if ( (sigp>sig(1)) .and. (sigp<sig(km)) ) then
              k1p = k1 + 1
              wp = log(sigp/sig(k1))/log(sig(k1p)/sig(k1))
              w1 = 1. - wp
              fp(i,j,n,nv) = w1*f(i,j,k1,nv) + wp*f(i,j,k1p,nv)
            else if ( (sigp>=sig(km)) .and. (sigp<=1.) ) then
              fp(i,j,n,nv) = f(i,j,km,nv)
            else if ( sigp>1. ) then
              fp(i,j,n,nv) = f(i,j,kbc,nv)                              &
                           & *exp(-rgas*lrate*log(sigp/sig(kbc))*rgti)
            else
            end if
          end do
        end do
      end do
!
      end subroutine intlog
