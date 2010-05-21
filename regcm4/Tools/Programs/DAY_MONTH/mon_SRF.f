      implicit none
      integer im,ny,nvar
      parameter (im=158,ny=107,nvar=27)
      character*12 in_fn(12)
      data in_fn/
     & 'SRF.19910101','SRF.19910201','SRF.19910301',
     & 'SRF.19910401','SRF.19910501','SRF.19910601',
     & 'SRF.19910701','SRF.19910801','SRF.19910901',
     & 'SRF.19911001','SRF.19911101','SRF.19911201'/
      character*11 outfn
      data outfn/'SRF.1991mon'/
      real*4  b(im,ny,nvar),c(im,ny,nvar)
      integer i,j,n,mrec,nrec,nday,month,nrecord
      open(20,file=outfn,form='unformatted',recl=im*ny*4
     &                  ,access='direct')
      mrec=0
      do month=1,12
         open(10,file=in_fn(month),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         nrec=0
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
      stop
      end
