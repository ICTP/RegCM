
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wheneq(n,array,inc,target,index,nval)
c
c $Id: wheneq.F,v 1.2 1995/02/28 22:19:52 ccmproc2 Exp $
c $Author: ccmproc2 $
c
      implicit none
      integer n,inc,nval
      integer index(*),
     $        array(*),
     $        target
      integer ina,i
c
      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do i=1,n
         if(array(ina) .eq. target) then
           nval=nval+1
           index(nval)=i
         end if
         ina=ina+inc
      enddo
      return
      end
