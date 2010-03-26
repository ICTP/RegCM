      subroutine writegrads(idout,outvar,nx,ny,nk,nrec)
 
      implicit none
!
! Dummy arguments
!
      integer :: idout , nk , nrec , nx , ny
      real(4) , dimension(nx,ny,nk) :: outvar
      intent (in) idout , nk , nx , ny , outvar
      intent (inout) nrec
!
! Local variables
!
      integer :: i , j , k
!
      do k = 1 , nk
        nrec = nrec + 1
        write (idout,rec=nrec) ((outvar(i,j,k),i=1,nx),j=1,ny)
      end do
 
      end subroutine writegrads
