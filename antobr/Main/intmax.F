
      integer function intmax(n,ix,inc)
c
c $Id: intmax.F,v 1.1.1.1 1995/02/09 23:26:15 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer n,inc
      integer ix(*)
c
      integer i,mx
c
      mx = ix(1)
      intmax = 1
      do i=1+inc,inc*n,inc
         if (ix(i).gt.mx) then
            mx = ix(i)
            intmax = i
         end if
      end do
      return
      end
