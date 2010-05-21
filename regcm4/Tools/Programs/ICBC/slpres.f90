      subroutine slpres(h,t,pstar,ht,tg,slp1,slp2,sig,im,jm,km)
      use mod_preproc_param , only : ptop
      use mod_constants , only : bltop , gti , rgas , lrate , stdt
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km
      real(4) , dimension(im,jm,km) :: h , t
      real(4) , dimension(im,jm) :: ht , pstar , slp1 , slp2 , tg
      real(4) , dimension(km) :: sig
      intent (in) h , ht , im , jm , km , pstar , sig , t , tg
      intent (out) slp1 , slp2
!
! Local variables
!
      integer :: i , j , k , kbc
      real(4) :: tsfc
!
      do k = 1 , km
        if ( sig(k)<bltop ) kbc = k
      end do
      do j = 1 , jm
        do i = 1 , im
          tsfc = t(i,j,kbc) - lrate*(h(i,j,kbc)-ht(i,j))
          slp1(i,j) = pstar(i,j)                                        &
                    & *exp(-gti/(rgas*lrate)*log(1.-ht(i,j)*lrate/tsfc))
        end do
      end do
 
      do j = 1 , jm
        do i = 1 , im
          slp2(i,j) = pstar(i,j)                                        &
                    & *exp(gti*ht(i,j)/(rgas*0.5*(tg(i,j)+stdt)))
        end do
      end do
 
      end subroutine slpres
