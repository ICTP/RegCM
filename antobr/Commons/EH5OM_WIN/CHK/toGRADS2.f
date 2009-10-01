      implicit none
      integer*2 ia(192,96)
      real*4  a(192,96)
      real*8  offset,scale
      integer i,j,k
      open(10,file=
     & '/home/RAID2-D13/EH5OM/RF/1961/EHgRF1961JAN'
     &       ,form='unformatted'
     &       ,recl=192*96*2+16,access='direct')
      open(20,file='JAN61g.dat',form='unformatted'
     &       ,recl=192*96*4,access='direct')
      do k=1,105
         read(10,rec=k) offset,scale,ia
         do j=1,96
         do i=1,192
            a(i,j)=ia(i,97-j)*scale+offset
         enddo
         enddo
         write(20,rec=k) a
      enddo
      stop
      end
