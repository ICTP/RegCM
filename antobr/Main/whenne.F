cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine whenne(n,array,inc,target,index,nval)
c
c $Id: whenne.F,v 1.1.1.1 1995/02/09 23:26:18 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer n,inc
      integer array(*),target
      integer index(*),nval
c
      integer i,ina
c
      ina=1
      nval=0
      if(inc.lt.0)ina=(-inc)*(n-1)+1
      do i=1,n
         if (array(ina).ne.target) then
            nval=nval+1
            index(nval)=i
         end if
         ina = ina+inc
      end do
      return
      end
