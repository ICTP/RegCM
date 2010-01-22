      implicit none
      real*4  s(50)
      real*4  a(46,22,41)
      integer i,j,k
      do j=1,22
      do i=1,46
         read(48) s
         do k=1,41
            a(i,j,k)=s(k)
         enddo
      enddo
      enddo
      open(10,file='CHK48.dat',form='unformatted'
     &       ,recl=46*22,access='direct')
      do k=1,41
         write(10,rec=k) ((a(i,j,k),i=1,46),j=1,22)
      enddo
      stop
      end
