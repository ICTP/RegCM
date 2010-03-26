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

      subroutine height(fp,f,f2,nta,npsa,zs,sig,im,jm,km,kp,nhga,p,n3d, &
                      & n2d,pt,im1,jm1)
 
      use mod_constants , only : rgti , rgas , lrate , bltop
      implicit none
!
! Dummy arguments
!
      integer :: im , im1 , jm , jm1 , km , kp , n2d , n3d , nhga ,     &
               & npsa , nta
      real(4) :: pt
      real(4) , dimension(im,jm,km,n3d) :: f
      real(4) , dimension(im,jm,n2d) :: f2
      real(4) , dimension(im,jm,kp,n3d) :: fp
      real(4) , dimension(kp) :: p
      real(4) , dimension(km) :: sig
      real(4) , dimension(im,jm) :: zs
      intent (in) f , f2 , im , im1 , jm , jm1 , km , kp , n2d , n3d ,  &
                & nhga , npsa , nta , p , pt , sig , zs
      intent (out) fp
!
! Local variables
!
      real(4) :: psfc , temp , wb , wt
      integer :: i , j , k , kb , kbc , kt , n
      real(4) , dimension(km+1) :: psig
      real(4) , dimension(im,jm) :: pstar
!
      do k = 1 , km
        if ( sig(k)<bltop ) kbc = k
      end do
      pstar = f2(:,:,npsa)
      do j = 1 , jm1
        do i = 1 , im1
          do k = 1 , km
            psig(k) = sig(k)*(pstar(i,j)-pt) + pt
          end do
          psfc = pstar(i,j)
          do n = 1 , kp
            kt = 1
            do k = 1 , km
              if ( psig(k)<p(n) ) kt = k
            end do
            kb = kt + 1
            if ( p(n)<=psig(1) ) then
              temp = f(i,j,1,nta)
              fp(i,j,n,nhga) = f(i,j,1,nhga)                            &
                             & + rgas*temp*log(psig(1)/p(n))*rgti
            else if ( (p(n)>psig(1)) .and. (p(n)<psig(km)) ) then
              wt = log(psig(kb)/p(n))/log(psig(kb)/psig(kt))
              wb = log(p(n)/psig(kt))/log(psig(kb)/psig(kt))
              temp = wt*f(i,j,kt,nta) + wb*f(i,j,kb,nta)
              temp = (temp+f(i,j,kb,nta))/2.
              fp(i,j,n,nhga) = f(i,j,kb,nhga)                           &
                             & + rgas*temp*log(psig(kb)/p(n))
            else if ( (p(n)>=psig(km)) .and. (p(n)<=psfc) ) then
              temp = f(i,j,km,nta)
              fp(i,j,n,nhga) = zs(i,j) + rgas*temp*log(psfc/p(n))*rgti
            else if ( p(n)>psfc ) then
              temp = f(i,j,kbc,nta) - lrate*(f(i,j,kbc,nta)-zs(i,j))
              fp(i,j,n,nhga) = zs(i,j) - (temp/lrate)                   &
                             & *(1.-exp(-rgas*lrate*log(p(n)/psfc)      &
                             & *rgti))
            else
            end if
          end do
        end do
      end do
      end subroutine height
