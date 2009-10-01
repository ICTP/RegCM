c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vcheki (ier,numerr,aname)
      implicit none
      integer ier,numerr
      character*8 aname
c
c  flag of detected errors in linear algebra routines
c
      if (ier.ne.0) then
         numerr=numerr+1
         print 1,aname,ier
      endif
    1 format('0 error in determination of ',a8,' using library routine
     a   ier=',i4)
c
      return
      end
