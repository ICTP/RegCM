      subroutine intlog(fp,f,ps,p3d,im,jm,km,p,kp)
      use mod_constants , only : rgas , rgti , lrate
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km , kp
      real , dimension(im,jm,km) :: f , p3d
      real , dimension(im,jm,kp) :: fp
      real , dimension(kp) :: p
      real , dimension(im,jm) :: ps
      intent (in) f , im , jm , km , kp , p , p3d , ps
      intent (out) fp
!
! Local variables
!
      real :: sigp , w1 , wp
      integer :: i , j , k , k1 , k1p , kbc , n
      real , dimension(61) :: sig
      real , parameter :: bltop = 0.96
!
!     INTLOG IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
!     LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF TLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATURES IN THE LAYER.
 
!
!**   FIND FIRST SIGMA LEVEL ABOVE BOUNDARY LAYER (LESS THAN SIG=BLTOP)
      do j = 1 , jm
        do i = 1 , im
          if ( ps(i,j)>-9995.0 ) then
            do k = 1 , km
              sig(k) = p3d(i,j,k)/ps(i,j)
              if ( sig(k)<bltop ) kbc = k
            end do
            do n = 1 , kp
              sigp = p(n)/ps(i,j)
              k1 = 0
              do k = 1 , km
                if ( sigp>sig(k) ) k1 = k
              end do
              if ( sigp<=sig(1) ) then
                fp(i,j,n) = f(i,j,1)
              else if ( (sigp>sig(1)) .and. (sigp<sig(km)) ) then
                k1p = k1 + 1
                wp = log(sigp/sig(k1))/log(sig(k1p)/sig(k1))
                w1 = 1. - wp
                fp(i,j,n) = w1*f(i,j,k1) + wp*f(i,j,k1p)
              else if ( (sigp>=sig(km)) .and. (sigp<=1.) ) then
                fp(i,j,n) = f(i,j,km)
              else if ( sigp>1. ) then
                fp(i,j,n) = f(i,j,kbc)                                  &
                          & *exp(+rgas*lrate*log(sigp/sig(kbc))*rgti)
!               ***** FROM R. ERRICO, SEE ROUTINE HEIGHT *****
              else
              end if
            end do
          else
            do n = 1 , kp
              fp(i,j,n) = -9999.0
            end do
          end if
        end do
      end do
 
      end subroutine intlog
