      implicit none
      integer*2 ia(44,44)
      real*4  a(44,44)
      real*8  offset,scale
      integer i,j,k
      open(10,file='EH_RF1961JAN',form='unformatted'
     &       ,recl=44*44*2+16,access='direct')
      open(20,file='JAN61.dat',form='unformatted'
     &       ,recl=44*44*4,access='direct')
      do k=1,105
         read(10,rec=k) offset,scale,ia
         do j=1,44
         do i=1,44
            a(i,j)=ia(i,j)*scale+offset
         enddo
         enddo
         write(20,rec=k) a
      enddo
      stop
      end
