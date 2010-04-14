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

      subroutine intlin(fp,f,pstar,sig,im,jm,km,p,kp)
      use mod_preproc_param , only : ptop
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km , kp
      real(4) , dimension(im,jm,km) :: f
      real(4) , dimension(im,jm,kp) :: fp
      real(4) , dimension(kp) :: p
      real(4) , dimension(im,jm) :: pstar
      real(4) , dimension(km) :: sig
      intent (in) f , im , jm , km , kp , p , pstar , sig
      intent (out) fp
!
! Local variables
!
      integer :: i , j , k , k1 , k1p , n
      real(4) :: sigp , w1 , wp
!
!     INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!     HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!     IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
 
      do j = 1 , jm
        do i = 1 , im
          do n = 1 , kp
            sigp = (p(n)-ptop)/(pstar(i,j)-ptop)
            k1 = 0
            do k = 1 , km
              if ( sigp>sig(k) ) k1 = k
            end do
            if ( sigp<=sig(1) ) then
              fp(i,j,n) = f(i,j,1)
            else if ( (sigp>sig(1)) .and. (sigp<sig(km)) ) then
              k1p = k1 + 1
              wp = (sigp-sig(k1))/(sig(k1p)-sig(k1))
              w1 = 1. - wp
              fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,k1p)
            else if ( sigp>=sig(km) ) then
              fp(i,j,n) = f(i,j,km)
            else
            end if
          end do
        end do
      end do
      end subroutine intlin
