      integer function isrchfgt(n, array, inc, target)
      implicit none
c
c $Id: isrchfgt.F,v 1.1.1.1 1995/02/09 23:26:15 ccm2 Exp $
c $Author: ccm2 $
c
      integer n,inc
      real*8  array(*),target
      integer i,ind
c
      if (n.le.0) then
         isrchfgt = 0
         return
      end if
      ind = 1
      do i=1,n
         if (array(ind).gt.target) then
            isrchfgt = i
            return
         else
            ind = ind + inc
         end if
      end do
      isrchfgt = n + 1
      return
      end
