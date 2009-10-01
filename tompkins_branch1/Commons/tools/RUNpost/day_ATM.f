      implicit none
      character*14 in_fn(72)
      data in_fn/
     & 'ATM.1990010100','ATM.1990010600','ATM.1990011100',
     & 'ATM.1990011600','ATM.1990012100','ATM.1990012600',
     & 'ATM.1990020100','ATM.1990020600','ATM.1990021100',
     & 'ATM.1990021600','ATM.1990022100','ATM.1990022600',
     & 'ATM.1990030100','ATM.1990030600','ATM.1990031100',
     & 'ATM.1990031600','ATM.1990032100','ATM.1990032600',
     & 'ATM.1990040100','ATM.1990040600','ATM.1990041100',
     & 'ATM.1990041600','ATM.1990042100','ATM.1990042600',
     & 'ATM.1990050100','ATM.1990050600','ATM.1990051100',
     & 'ATM.1990051600','ATM.1990052100','ATM.1990052600',
     & 'ATM.1990060100','ATM.1990060600','ATM.1990061100',
     & 'ATM.1990061600','ATM.1990062100','ATM.1990062600',
     & 'ATM.1990070100','ATM.1990070600','ATM.1990071100',
     & 'ATM.1990071600','ATM.1990072100','ATM.1990072600',
     & 'ATM.1990080100','ATM.1990080600','ATM.1990081100',
     & 'ATM.1990081600','ATM.1990082100','ATM.1990082600',
     & 'ATM.1990090100','ATM.1990090600','ATM.1990091100',
     & 'ATM.1990091600','ATM.1990092100','ATM.1990092600',
     & 'ATM.1990100100','ATM.1990100600','ATM.1990101100',
     & 'ATM.1990101600','ATM.1990102100','ATM.1990102600',
     & 'ATM.1990110100','ATM.1990110600','ATM.1990111100',
     & 'ATM.1990111600','ATM.1990112100','ATM.1990112600',
     & 'ATM.1990120100','ATM.1990120600','ATM.1990121100',
     & 'ATM.1990121600','ATM.1990122100','ATM.1990122600'/
      character*12 outfn(12)
      data outfn/
     & 'ATM.19900101','ATM.19900201','ATM.19900301',
     & 'ATM.19900401','ATM.19900501','ATM.19900601',
     & 'ATM.19900701','ATM.19900801','ATM.19900901',
     & 'ATM.19901001','ATM.19901101','ATM.19901201'/
      real*4  b(398,330),c(398,330,113)
      integer i,j,n,mrec,nrec,nday,mday,month,nrecord,n5
      do month=11,12
         open(20,file=outfn(month),form='unformatted'
     &          ,recl=398*330*4,access='direct')
         mrec=0
         do n5=1,6
            open(10,file=in_fn((month-1)*6+n5),form='unformatted'
     &          ,recl=398*330*4,access='direct')
            nrec=0
            if(month.eq.1.and.n5.eq.1) nrec=113
            nrecord=5
            if(n5.eq.6) then
               if(month.eq.1.or.month.eq.3.or.month.eq.5 .or.
     &            month.eq.7.or.month.eq.8.or.month.eq.10.or.
     &            month.eq.12) then
                  nrecord=6
               else if(month.eq.4.or.month.eq.6.or.month.eq.9.or.
     &                 month.eq.11) then
                  nrecord=5
               else
Cnormal Year     nrecord=3
                  nrecord=3

C Leap  Year     nrecord=4
!                 nrecord=4
               endif
            endif
            do nday=1,nrecord
               do n=1,113
                  do j=1,330
                  do i=1,398
                     c(i,j,n)=0.0
                  enddo
                  enddo
               enddo
               do mday=1,4
                  do n=1,113
                     nrec=nrec+1
                     read(10,rec=nrec) b
                     do j=1,330
                     do i=1,398
                        c(i,j,n)=c(i,j,n)+b(i,j)
                     enddo
                     enddo
                  enddo
               enddo
               do n=1,113
                  do j=1,330
                  do i=1,398
                     c(i,j,n)=c(i,j,n)/4.
                  enddo
                  enddo
                  mrec=mrec+1
                  write(20,rec=mrec) ((c(i,j,n),i=1,398),j=1,330)
               enddo
            enddo
            close(10)
         enddo
         close(20)
      enddo
      stop
      end
