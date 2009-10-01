ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine maximi(array,ix,kx,ks,ke,max,istart,iend)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer ix,kx,ks,ke,istart,iend
      real*8  array(ix,kx)
      integer max(ix)
      integer i,k
      real*8  x,xar

      do 200 i=istart,iend
       max(i)=ks
       x=array(i,ks)
       do 100 k=ks,ke
        xar=array(i,k)
        if(xar.ge.x) then
         x=xar
         max(i)=k
        endif
100    continue
200   continue
c
      return
      end
