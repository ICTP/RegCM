      subroutine htsig_o(t,h,pstar,ht,sig,ptop,im,jm,km)
      use mod_constants, only : rgas , rgti
      implicit none
!
! Dummy arguments
!
      integer :: im , jm , km
      real :: ptop
      real , dimension(im,jm,km) :: h , t
      real , dimension(im,jm) :: ht , pstar
      real , dimension(km) :: sig
      intent (in) ht , im , jm , km , pstar , ptop , sig , t
      intent (inout) h
!
! Local variables
!
      real :: tbar
      integer :: i , j , k
!
      do j = 1 , jm
        do i = 1 , im
          h(i,j,km) = ht(i,j) + rgas*rgti*t(i,j,km)                     &
                    & *log(pstar(i,j)/((pstar(i,j)-ptop)*sig(km)+ptop))
        end do
      end do
      do k = km - 1 , 1 , -1
        do j = 1 , jm
          do i = 1 , im
            tbar = 0.5*(t(i,j,k)+t(i,j,k+1))
            h(i,j,k) = h(i,j,k+1)                                       &
                     & + rgas*rgti*tbar*log(((pstar(i,j)-ptop)*sig(k+1) &
                     & +ptop)/((pstar(i,j)-ptop)*sig(k)+ptop))
          end do
        end do
      end do
      end subroutine htsig_o
