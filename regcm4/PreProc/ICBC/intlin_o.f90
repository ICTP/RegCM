      subroutine intlin_o(fp,f,pstar,sig,ptop,im,jm,km,p,kp)
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km , kp
      real :: ptop
      real , dimension(im,jm,km) :: f
      real , dimension(im,jm,kp) :: fp
      real , dimension(kp) :: p
      real , dimension(im,jm) :: pstar
      real , dimension(km) :: sig
      intent (in) f , im , jm , km , kp , p , pstar , ptop , sig
      intent (out) fp
!
! Local variables
!
      integer :: i , j , k , k1 , k1p , n
      real :: sigp , w1 , wp
!
!     INTLIN IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE
!     HUMIDITY. THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION
!     IS NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!
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
      end subroutine intlin_o
