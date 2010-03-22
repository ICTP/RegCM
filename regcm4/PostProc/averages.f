      subroutine average_day(ix,iy,filetoread,filetowrite,ndays)
! this routine averages a month over the number of days of the month
      implicit none
      integer ix,iy,nvar,iopen,ndays
      parameter (nvar=27)
      real*4  b(ix,iy,nvar),c(ix,iy,nvar)
      character*14 filetoread
Cf2py intent(in) ix 
Cf2py intent(in) iy 
Cf2py intent(in) filetoread 
Cf2py intent(in) filetowrite 
Cf2py intent(in) ndays
      character*23 filetowrite
      integer i,j,n,mrec,nrec,nday,mday,month,nrecord
      iopen=0
         iopen=iopen+1
         open(10,file=filetoread,form='unformatted'
     &          ,recl=ix*iy*4,access='direct')
         mrec=0
         nrec=0
         open(20,file=filetowrite,form='unformatted'
     &          ,recl=ix*iy*4,access='direct')
         do nday=1,ndays
            do n=1,21
               do j=1,iy
               do i=1,ix
                  c(i,j,n)=0.0
               enddo
               enddo
            enddo
            do j=1,iy
            do i=1,ix
               c(i,j,nvar-5)= -1.e20
               c(i,j,nvar-4)=  1.e20
               c(i,j,nvar-3)= -1.e20
               c(i,j,nvar-2)=  1.e20
               c(i,j,nvar-1)= -1.e20
               c(i,j,nvar  )=  1.e20
            enddo
            enddo
            do mday=1,8
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) ((b(i,j,n),i=1,ix),j=1,iy)
               enddo
               do n=1,nvar-6
                  do j=1,iy
                  do i=1,ix
                     c(i,j,n)=c(i,j,n)+b(i,j,n)
                  enddo
                  enddo
               enddo
               do j=1,iy
               do i=1,ix
         if(c(i,j,nvar-5).lt.b(i,j,nvar-5)) c(i,j,nvar-5)=b(i,j,nvar-5)
         if(c(i,j,nvar-4).gt.b(i,j,nvar-4)) c(i,j,nvar-4)=b(i,j,nvar-4)
         if(c(i,j,nvar-3).lt.b(i,j,nvar-3)) c(i,j,nvar-3)=b(i,j,nvar-3)
         if(c(i,j,nvar-2).gt.b(i,j,nvar-2)) c(i,j,nvar-2)=b(i,j,nvar-2)
         if(c(i,j,nvar-1).lt.b(i,j,nvar-1)) c(i,j,nvar-1)=b(i,j,nvar-1)
         if(c(i,j,nvar  ).gt.b(i,j,nvar  )) c(i,j,nvar  )=b(i,j,nvar  )
               enddo
               enddo
            enddo
            do n=1,nvar-6
               do j=1,iy
               do i=1,ix
                  c(i,j,n)=c(i,j,n)/8.
               enddo
               enddo
            enddo
            do n=1,nvar
               mrec=mrec+1
               write(20,rec=mrec) ((c(i,j,n),i=1,ix),j=1,iy)
            enddo
         enddo
         close(10)
         close(20)
      end subroutine average_day
!
!
!
! 
      subroutine month_average(im,ny,in_fn,outfn,ndays)
      implicit none
      integer, intent(in) :: im
      integer, intent(in) :: ny 
      character (len=22), intent(in) ::   outfn
      integer :: ndays(12)
      integer :: iopen
      integer, parameter :: nvar=27
Cf2py intent(in) im 
Cf2py intent(in) ny 
Cf2py character in_fn(12*17)
      character (len=17), intent(in) ::  in_fn(12)
!      character*23 in_fn(12)
Cf2py intent(in) outfn 
Cf2py intent(in) ndays(12)
      real*4  b(im,ny,nvar),c(im,ny,nvar)
      integer i,j,n,mrec,nrec,nday,month,nrecord
      open(20,file=outfn,form='unformatted',recl=im*ny*4
     &                  ,access='direct')
      mrec=0
      iopen=0
      do month=1,12
         iopen=iopen+1
         open(10,file=in_fn(iopen),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         nrecord=ndays(iopen)
         nrec=0
         print*, 'processing file', in_fn(iopen),ndays(iopen)
         do n=1,nvar
            do j=1,ny
            do i=1,im
               c(i,j,n)=0.0
            enddo
            enddo
         enddo
         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) ((b(i,j,n),i=1,im),j=1,ny)
               do j=1,ny
               do i=1,im
                  c(i,j,n)=c(i,j,n)+b(i,j,n)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do j=1,ny
            do i=1,im
               c(i,j,n)=c(i,j,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(i,j,n),i=1,im),j=1,ny)
         enddo
         close(10)
      enddo
      close(20)
      end subroutine month_average
