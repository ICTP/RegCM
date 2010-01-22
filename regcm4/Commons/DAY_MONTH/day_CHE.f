      implicit none
      integer im,ny,nvar
      parameter (im=158,ny=107,nvar=156)
      real*4  b(im,ny),c(im,ny,nvar)
      character*14 in_fn(12)
      data in_fn/
     & 'CHE.1991010100','CHE.1991020100','CHE.1991030100',
     & 'CHE.1991040100','CHE.1991050100','CHE.1991060100',
     & 'CHE.1991070100','CHE.1991080100','CHE.1991090100',
     & 'CHE.1991100100','CHE.1991110100','CHE.1991120100'/
      character*12 outfn(12)
      data outfn/
     & 'CHE.19910101','CHE.19910201','CHE.19910301',
     & 'CHE.19910401','CHE.19910501','CHE.19910601',
     & 'CHE.19910701','CHE.19910801','CHE.19910901',
     & 'CHE.19911001','CHE.19911101','CHE.19911201'/
      integer i,j,n,mrec,nrec,nday,mday,month,nrecord
      do month=1,12
         open(20,file=outfn(month),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         mrec=0
         open(10,file=in_fn(month),form='unformatted'
     &          ,recl=im*ny*4,access='direct')
         nrec=0
         if(month.eq.1) nrec=nvar
         nrecord=30
         if(month.eq.1.or.month.eq.3.or.month.eq.5 .or.
     &      month.eq.7.or.month.eq.8.or.month.eq.10.or.
     &      month.eq.12) then
            nrecord=31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or.
     &           month.eq.11) then
            nrecord=30
         else
Cnormal Year     nrecord=28
            nrecord=28

C Leap  Year     nrecord=29
!           nrecord=29
         endif
         do nday=1,nrecord
            do n=1,nvar
               do j=1,ny
               do i=1,im
                  c(i,j,n)=0.0
               enddo
               enddo
            enddo
            do mday=1,4
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) b
                  do j=1,ny
                  do i=1,im
                     c(i,j,n)=c(i,j,n)+b(i,j)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do j=1,ny
               do i=1,im
                  c(i,j,n)=c(i,j,n)/4.
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(i,j,n),i=1,im),j=1,ny)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      stop
      end
