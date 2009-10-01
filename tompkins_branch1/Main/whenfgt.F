
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine whenfgt(n,array,inc,target,index,nval)
c
c $Id: whenfgt.F,v 1.1.1.1 1995/02/09 23:26:17 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer n,inc,nval
      integer index(*)
      real*8  array(*),target
      integer ina,i
c
      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do 100 i=1,n
          if(array(ina) .gt. target) then
          nval=nval+1
          index(nval)=i
          end if
          ina=ina+inc
 100    continue
      return
      end
