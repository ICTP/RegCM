      subroutine calcslp(f3,f2,nhga,nta,npsa,zs,slp,sig,im,jm,km,n3d,   &
                       & n2d,nx1,ny1)
      use mod_constants , only : rgas , gti , bltop , lrate
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km , n2d , n3d , nhga , npsa , nta , nx1 ,   &
               & ny1 , slp
      real(4) , dimension(im,jm,n2d) :: f2
      real(4) , dimension(im,jm,km,n3d) :: f3
      real(4) , dimension(km) :: sig
      real(4) , dimension(im,jm) :: zs
      intent (in) f3 , im , jm , km , n2d , n3d , nhga , npsa , nta ,   &
                & nx1 , ny1 , sig , slp , zs
      intent (inout) f2
!
! Local variables
!
      integer :: i , j , k , kbc
      real(4) , dimension(im,jm) :: pstar
      real(4) :: tsfc
!
      pstar(:,:) = f2(:,:,npsa)
      do k = 1 , km
        if ( sig(k)<bltop ) kbc = k
      end do
      do j = 1 , ny1
        do i = 1 , nx1
          tsfc = f3(i,j,kbc,nta) - lrate*(f3(i,j,kbc,nhga)-zs(i,j))
          f2(i,j,slp) = pstar(i,j)                                      &
                      & *exp(-gti/(rgas*lrate)*log(1.-zs(i,j)           &
                      & *lrate/tsfc))
        end do
      end do
 
      end subroutine calcslp
