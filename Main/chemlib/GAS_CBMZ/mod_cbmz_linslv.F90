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

module mod_cbmz_linslv

  use mod_intkinds
  use mod_realkinds

  private

  public :: linslv

  contains

! linslv.f
!    from bnrchemv5.f
!
! This contains the subroutines LINSLV and RESOLV
!   for inverting a matrix.
!   (Dates from Prather, 1988, and Numerical Recipes)
!
! It solves for X in the matrix equation AX=B.
!
!     call LINSLV(A, B, X, N) for AX=B, used dimension N
!          A is updated to reduced form  by LINSLV.
!          NOTE: declared dimension of A, B, and X are set at 100 (not N
!
!     call RESOLV(A, B, X, N) for AX=B, with A reduced from LINSLV.
!          (if error in last solution is small)
!
! THE ORIGINAL CALL:
!       err=abs(xoo(i)/xr( 1 ,ic))
!      if(errxo.lt.ermat) lnumat=.false.
!      if(lnumat) call LINSLV(fxo,xoo,nchem)
!      if(.not.lnumat) call RESOLV(fxo,xoo,nchem)
!
! MODIFICATIONS FROM bnrchemv5.f:
!   include 'commq5.f' and 'commb5.f' are removed.
!   matrix A is parameter of subroutine.
!
!
! ------------------------------------------------------
 
    subroutine linslv(a,b,x,n)

!  *** SUB -LINSLV- SOLVES THE FOLLOWING MATRIX EQUATION--
!  ***     A(N,N)*X(N) = B(N)  BY REDUCING THE A-MATRIX IN PLACE.
!  *** THE ENTRY -RESOLV- ASSUMES THAT THE A-MATRIX HAS BEEN PROPERLY
!  ***     REDUCED AND JUST SOLVES FOR X(N).  THIS OPTION SAVES TIME
!  ***     WHEN THE SYSTEM IS TO BE RESOLVED WITH A NEW B-VECTOR..  (Pra
   
      implicit none
!
      integer(ik4) :: n
      real(rk8) , dimension(100,100) :: a
      real(rk8) , dimension(100) :: b , x
      intent (inout) a
!
      real(rk8) :: div , smax
      integer(ik4) :: i , j , jp , jp1 , k , kr , krm1 , krmax , krp1
      integer(ik4) , dimension(100) :: ipa
      real(rk8) , dimension(100) :: s
!
      do kr = 1 , n
        do k = 1 , n
          s(k) = a(k,kr)
        end do
        if ( kr /= 1 ) then
          krm1 = kr - 1
          do j = 1 , krm1
            jp = ipa(j)
            a(j,kr) = s(jp)
            s(jp) = s(j)
            jp1 = j + 1
            do i = jp1 , n
              s(i) = s(i) - a(i,j)*a(j,kr)
            end do
          end do
        end if
        krmax = kr
        smax = dabs(s(kr))
        do i = kr , n
          if ( dabs(s(i)) > smax ) then
            krmax = i
            smax = dabs(s(i))
          end if
        end do
        ipa(kr) = krmax
        a(kr,kr) = s(krmax)
        div = 1.0D0/s(krmax)
        s(krmax) = s(kr)
        if ( kr /= n ) then
          krp1 = kr + 1
          do i = krp1 , n
            a(i,kr) = s(i)*div
          end do
        end if
      end do
      call resolv(a,b,x,ipa,n)
    end subroutine linslv

    subroutine resolv(a,b,x,ipa,n)

      implicit none
!
      integer(ik4) :: n
      real(rk8) , dimension(100,100) :: a
      real(rk8) , dimension(100) :: b , x
      integer(ik4) , dimension(100) :: ipa
      intent (in) a , b , ipa , n
      intent (inout) x
!
      integer(ik4) :: i , ii , iip1 , ip , ip1 , j
      real(rk8) , dimension(100) :: s
      real(rk8) :: summ
!
      do i = 1 , n
        s(i) = b(i)
      end do
      do i = 1 , n
        ip = ipa(i)
        x(i) = s(ip)
        s(ip) = s(i)
        if ( i == n ) exit
        ip1 = i + 1
        do j = ip1 , n
          s(j) = s(j) - a(j,i)*x(i)
        end do
      end do
 
      do i = 1 , n
        ii = n + 1 - i
        summ = x(ii)
!       NOTE POSSIBLE ERROR IN THIS NEXT CHANGED LINE.
        if ( ii < n ) then
          iip1 = ii + 1
          do j = iip1 , n
            summ = summ - a(ii,j)*x(j)
          end do
        end if
        x(ii) = summ/a(ii,ii)
      end do
 
    end subroutine resolv

end module mod_cbmz_linslv
