      implicit none
      character*14 in_fn(12)
      data in_fn/
     & 'SRF.1990010100','SRF.1990020100','SRF.1990030100',
     & 'SRF.1990040100','SRF.1990050100','SRF.1990060100',
     & 'SRF.1990070100','SRF.1990080100','SRF.1990090100',
     & 'SRF.1990100100','SRF.1990110100','SRF.1990120100'/
      character*12 outfn(12)
      data outfn/
     & 'SRF.19900101','SRF.19900201','SRF.19900301',
     & 'SRF.19900401','SRF.19900501','SRF.19900601',
     & 'SRF.19900701','SRF.19900801','SRF.19900901',
     & 'SRF.19901001','SRF.19901101','SRF.19901201'/
      real*4  b(119,98,27),c(119,98,27)
      integer i,j,n,mrec,nrec,nday,mday,month,nrecord
      do month=1,12
         open(10,file=in_fn(month),form='unformatted'
     &          ,recl=119*98*4,access='direct')
         nrec=0
         if(month.eq.1) nrec=27
         open(20,file=outfn(month),form='unformatted'
     &          ,recl=119*98*4,access='direct')
         mrec=0
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
            do n=1,21
               do j=1,98
               do i=1,119
                  c(i,j,n)=0.0
               enddo
               enddo
            enddo
            do j=1,98
            do i=1,119
               c(i,j,22)= -1.e20
               c(i,j,23)=  1.e20
               c(i,j,24)= -1.e20
               c(i,j,25)=  1.e20
               c(i,j,26)= -1.e20
               c(i,j,27)=  1.e20
            enddo
            enddo
            do mday=1,24
               do n=1,27
                  nrec=nrec+1
                  read(10,rec=nrec) ((b(i,j,n),i=1,119),j=1,98)
               enddo
               do n=1,21
                  do j=1,98
                  do i=1,119
                     c(i,j,n)=c(i,j,n)+b(i,j,n)
                  enddo
                  enddo
               enddo
               do j=1,98
               do i=1,119
                  if(c(i,j,22).lt.b(i,j,22)) c(i,j,22)=b(i,j,22)
                  if(c(i,j,23).gt.b(i,j,23)) c(i,j,23)=b(i,j,23)
                  if(c(i,j,24).lt.b(i,j,24)) c(i,j,24)=b(i,j,24)
                  if(c(i,j,25).gt.b(i,j,25)) c(i,j,25)=b(i,j,25)
                  if(c(i,j,26).lt.b(i,j,26)) c(i,j,26)=b(i,j,26)
                  if(c(i,j,27).gt.b(i,j,27)) c(i,j,27)=b(i,j,27)
               enddo
               enddo
            enddo
            do n=1,21
               do j=1,98
               do i=1,119
                  c(i,j,n)=c(i,j,n)/24.
               enddo
               enddo
            enddo
            do n=1,27
               mrec=mrec+1
               write(20,rec=mrec) ((c(i,j,n),i=1,119),j=1,98)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      stop
      end
