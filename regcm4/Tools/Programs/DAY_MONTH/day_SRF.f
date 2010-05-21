      implicit none
      integer im,ny,nvar
      parameter (im=158,ny=107,nvar=27)
      real*4  b(im,ny,nvar),c(im,ny,nvar)
      character*14 in_fn(12)
      data in_fn/
     & 'SRF.1991010100','SRF.1991020100','SRF.1991030100',
     & 'SRF.1991040100','SRF.1991050100','SRF.1991060100',
     & 'SRF.1991070100','SRF.1991080100','SRF.1991090100',
     & 'SRF.1991100100','SRF.1991110100','SRF.1991120100'/
      character*12 outfn(12)
      data outfn/
     & 'SRF.19910101','SRF.19910201','SRF.19910301',
     & 'SRF.19910401','SRF.19910501','SRF.19910601',
     & 'SRF.19910701','SRF.19910801','SRF.19910901',
     & 'SRF.19911001','SRF.19911101','SRF.19911201'/
      integer i,j,n,mrec,nrec,nday,mday,month,nrecord
      do month=1,12
         open(20,file=outfn(month),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         mrec=0
         open(10,file=in_fn(month),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         nrec=0
         if(month.eq.1) nrec=nvar
         if(month.eq.1.or.month.eq.3.or.month.eq.5 .or.
     &      month.eq.7.or.month.eq.8.or.month.eq.10.or.
     &      month.eq.12) then
            nrecord=31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or.
     &           month.eq.11) then
            nrecord=30
         else
Cnormal Year     nrecord=3
            nrecord=28

C Leap  Year     nrecord=4
!           nrecord=29
         endif
         do nday=1,nrecord
            do n=1,21
               do j=1,ny
               do i=1,im
                  c(i,j,n)=0.0
               enddo
               enddo
            enddo
            do j=1,ny
            do i=1,im
               c(i,j,nvar-5)= -1.e20
               c(i,j,nvar-4)=  1.e20
               c(i,j,nvar-3)= -1.e20
               c(i,j,nvar-2)=  1.e20
               c(i,j,nvar-1)= -1.e20
               c(i,j,nvar  )=  1.e20
            enddo
            enddo
            do mday=1,24
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) ((b(i,j,n),i=1,im),j=1,ny)
               enddo
               do n=1,nvar-6
                  do j=1,ny
                  do i=1,im
                     c(i,j,n)=c(i,j,n)+b(i,j,n)
                  enddo
                  enddo
               enddo
               do j=1,ny
               do i=1,im
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
               do j=1,ny
               do i=1,im
                  c(i,j,n)=c(i,j,n)/24.
               enddo
               enddo
            enddo
            do n=1,nvar
               mrec=mrec+1
               write(20,rec=mrec) ((c(i,j,n),i=1,im),j=1,ny)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      stop
      end
