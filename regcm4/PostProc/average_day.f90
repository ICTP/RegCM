      subroutine average_day(ix,iy,filetoread,filetowrite,ndays)
!
! this routine averages a month over the number of days of the month
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nvar = 27
!
! Dummy arguments
!
      character(256) :: filetoread
      character(256) :: filetowrite
      integer :: ix , iy , ndays
      intent (in) filetoread , filetowrite , ix , iy , ndays
!
! Local variables
!
      real(4) , dimension(ix,iy,nvar) :: b , c
      integer :: i , iopen , j , mday , mrec , n , nday , nrec
!
      iopen = 0
      iopen = iopen + 1
      open (10,file=filetoread,form='unformatted',recl=ix*iy*4,         &
           &access='direct')
      mrec = 0
      nrec = 0
      open (20,file=filetowrite,form='unformatted',recl=ix*iy*4,        &
           &access='direct')
      do nday = 1 , ndays
        do n = 1 , 21
          do j = 1 , iy
            do i = 1 , ix
              c(i,j,n) = 0.0
            end do
          end do
        end do
        do j = 1 , iy
          do i = 1 , ix
            c(i,j,nvar-5) = -1.E20
            c(i,j,nvar-4) = 1.E20
            c(i,j,nvar-3) = -1.E20
            c(i,j,nvar-2) = 1.E20
            c(i,j,nvar-1) = -1.E20
            c(i,j,nvar) = 1.E20
          end do
        end do
        do mday = 1 , 8
          do n = 1 , nvar
            nrec = nrec + 1
            read (10,rec=nrec) ((b(i,j,n),i=1,ix),j=1,iy)
          end do
          do n = 1 , nvar - 6
            do j = 1 , iy
              do i = 1 , ix
                c(i,j,n) = c(i,j,n) + b(i,j,n)
              end do
            end do
          end do
          do j = 1 , iy
            do i = 1 , ix
              if ( c(i,j,nvar-5)<b(i,j,nvar-5) ) c(i,j,nvar-5)          &
                 & = b(i,j,nvar-5)
              if ( c(i,j,nvar-4)>b(i,j,nvar-4) ) c(i,j,nvar-4)          &
                 & = b(i,j,nvar-4)
              if ( c(i,j,nvar-3)<b(i,j,nvar-3) ) c(i,j,nvar-3)          &
                 & = b(i,j,nvar-3)
              if ( c(i,j,nvar-2)>b(i,j,nvar-2) ) c(i,j,nvar-2)          &
                 & = b(i,j,nvar-2)
              if ( c(i,j,nvar-1)<b(i,j,nvar-1) ) c(i,j,nvar-1)          &
                 & = b(i,j,nvar-1)
              if ( c(i,j,nvar)>b(i,j,nvar) ) c(i,j,nvar) = b(i,j,nvar)
            end do
          end do
        end do
        do n = 1 , nvar - 6
          do j = 1 , iy
            do i = 1 , ix
              c(i,j,n) = c(i,j,n)/8.
            end do
          end do
        end do
        do n = 1 , nvar
          mrec = mrec + 1
          write (20,rec=mrec) ((c(i,j,n),i=1,ix),j=1,iy)
        end do
      end do
      close (10)
      close (20)
      end subroutine average_day
