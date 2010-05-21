      implicit none
      integer i,j,n
      integer*2 ii2(240,121)
      real*4  ss4(240,121)
      real*4  ssmax,ssmin
      real*8  offset,scale

      open(10,file='tskinI2.dat',form='unformatted'
     &       ,recl=240*121*2,access='direct')
      open(21,file='tskin_00.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(22,file='tskin_06.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(23,file='tskin_12.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(24,file='tskin_18.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
!     ssmax= -1.E20
!     ssmin=  1.E20
      ssmax= 344.2406
      ssmin= 168.7359
      offset=0.5*(ssmax+ssmin)
      scale =(ssmax-offset)/32766.
      write(*,*) offset,scale
      do n=1,7456
         read(21,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               ii2(i,j) = nint((ss4(i,j)-offset)/scale)
            else
               ii2(i,j) = -32767
            endif
         enddo
         enddo
         write(10,rec=(n-1)*4+1) ii2
         read(22,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               ii2(i,j) = nint((ss4(i,j)-offset)/scale)
            else
               ii2(i,j) = -32767
            endif
         enddo
         enddo
         write(10,rec=(n-1)*4+2) ii2
         read(23,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               ii2(i,j) = nint((ss4(i,j)-offset)/scale)
            else
               ii2(i,j) = -32767
            endif
         enddo
         enddo
         write(10,rec=(n-1)*4+3) ii2
         read(24,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               ii2(i,j) = nint((ss4(i,j)-offset)/scale)
            else
               ii2(i,j) = -32767
            endif
         enddo
         enddo
         write(10,rec=(n-1)*4+4) ii2
      enddo
      close(10)
      close(21)
      close(22)
      close(23)
      close(23)
      stop
      end
