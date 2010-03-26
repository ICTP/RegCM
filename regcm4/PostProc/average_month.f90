      subroutine average_month(im,ny,in_fn,outfn,ndays)
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nvar = 27
!
! Dummy arguments
!
      integer :: im , ny
      character(256) :: outfn
      character(256) , dimension(12) :: in_fn
      integer , dimension(12) :: ndays
      intent (in) im , in_fn , ndays , ny , outfn
!
! Local variables
!
      real(4) , dimension(im,ny,nvar) :: b , c
      integer :: i , iopen , j , month , mrec , n , nday , nrec ,       &
               & nrecord
!
      open (20,file=outfn,form='unformatted',recl=im*ny*4,              &
           &access='direct')
      mrec = 0
      iopen = 0
      do month = 1 , 12
        iopen = iopen + 1
        open (10,file=in_fn(iopen),form='unformatted',recl=im*ny*4,     &
            & access='direct')
        nrecord = ndays(iopen)
        nrec = 0
        print * , 'processing file' , in_fn(iopen) , ndays(iopen)
        do n = 1 , nvar
          do j = 1 , ny
            do i = 1 , im
              c(i,j,n) = 0.0
            end do
          end do
        end do
        do nday = 1 , nrecord
          do n = 1 , nvar
            nrec = nrec + 1
            read (10,rec=nrec) ((b(i,j,n),i=1,im),j=1,ny)
            do j = 1 , ny
              do i = 1 , im
                c(i,j,n) = c(i,j,n) + b(i,j,n)
              end do
            end do
          end do
        end do
        do n = 1 , nvar
          do j = 1 , ny
            do i = 1 , im
              c(i,j,n) = c(i,j,n)/float(nrecord)
            end do
          end do
          mrec = mrec + 1
          write (20,rec=mrec) ((c(i,j,n),i=1,im),j=1,ny)
        end do
        close (10)
      end do
      close (20)
      end subroutine average_month
