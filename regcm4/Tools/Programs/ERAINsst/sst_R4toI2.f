      implicit none
      integer i,j,n
      integer*2 ii2(240,121)
      real*4  ss4(240,121)
      real*4  ssmax,ssmin

      open(21,file='sst_00.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(22,file='sst_06.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(23,file='sst_12.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      open(24,file='sst_18.dat',form='unformatted'
     &       ,recl=240*121*4,access='direct')
      ssmax= -1.E20
      ssmin=  1.E20
      do n=1,7456
         read(21,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               if(ss4(i,j).gt.ssmax) ssmax= ss4(i,j)
               if(ss4(i,j).lt.ssmin) ssmin= ss4(i,j)
            endif
         enddo
         enddo
         read(22,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               if(ss4(i,j).gt.ssmax) ssmax= ss4(i,j)
               if(ss4(i,j).lt.ssmin) ssmin= ss4(i,j)
            endif
         enddo
         enddo
         read(23,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               if(ss4(i,j).gt.ssmax) ssmax= ss4(i,j)
               if(ss4(i,j).lt.ssmin) ssmin= ss4(i,j)
            endif
         enddo
         enddo
         read(24,rec=n) ss4
         do j=1,121
         do i=1,240
            if(ss4(i,j).gt.-9990.) then
               if(ss4(i,j).gt.ssmax) ssmax= ss4(i,j)
               if(ss4(i,j).lt.ssmin) ssmin= ss4(i,j)
            endif
         enddo
         enddo
      enddo
      write(*,*) ssmax,ssmin
      stop
      end
