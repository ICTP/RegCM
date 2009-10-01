c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine minimi(array,ix,kx,ks,kend,kt,istart,iend)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer ix,kx,kend,istart,iend
      real*8  array(ix,kx)
      integer kt(ix),ks(ix)
      integer i,k
      real*8  x
c
      do 200 i=istart,iend
       kt(i)=ks(i)
       x=array(i,ks(i))
       do 100 k=ks(i)+1,kend
        if(array(i,k).lt.x)then
         x=array(i,k)
         kt(i)=k
        endif
100    continue
200   continue
c
      return
      end
