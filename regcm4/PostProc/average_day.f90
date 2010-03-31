!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine average_day(iy,jx,filetoread,filetowrite,ndays)
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
      integer :: iy , jx , ndays
      intent (in) filetoread , filetowrite , iy , jx , ndays
!
! Local variables
!
      real(4) , dimension(iy,jx,nvar) :: b , c
      integer :: i , iopen , j , mday , mrec , n , nday , nrec
!
      iopen = 0
      iopen = iopen + 1
      open (10,file=filetoread,form='unformatted',recl=iy*jx*4,         &
           &access='direct')
      mrec = 0
      nrec = 0
      open (20,file=filetowrite,form='unformatted',recl=iy*jx*4,        &
           &access='direct')
      do nday = 1 , ndays
        do n = 1 , 21
          do j = 1 , jx
            do i = 1 , iy
              c(i,j,n) = 0.0
            end do
          end do
        end do
        do j = 1 , jx
          do i = 1 , iy
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
            read (10,rec=nrec) ((b(i,j,n),i=1,iy),j=1,jx)
          end do
          do n = 1 , nvar - 6
            do j = 1 , jx
              do i = 1 , iy
                c(i,j,n) = c(i,j,n) + b(i,j,n)
              end do
            end do
          end do
          do j = 1 , jx
            do i = 1 , iy
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
          do j = 1 , jx
            do i = 1 , iy
              c(i,j,n) = c(i,j,n)/8.
            end do
          end do
        end do
        do n = 1 , nvar
          mrec = mrec + 1
          write (20,rec=mrec) ((c(i,j,n),i=1,iy),j=1,jx)
        end do
      end do
      close (10)
      close (20)
      end subroutine average_day
