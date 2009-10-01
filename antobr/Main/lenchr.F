c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function lenchr(chrstg)
c----------------------------------------------------------------------
c
cl         return position of right-most non-blank, non-null character
cl         in chrstg.
c
c----------------------------------------------------------------------
      implicit none
      character*(*) chrstg
      lenchr = len(chrstg)
      if (chrstg(lenchr:lenchr).eq.' '.or.
     1    chrstg(lenchr:lenchr).eq.char(0)) then
10       continue
         lenchr = lenchr - 1
         if (chrstg(lenchr:lenchr).eq.' '.or.
     1       chrstg(lenchr:lenchr).eq.char(0)) then
            if (lenchr.gt.0) go to 10
         endif
      endif
      return
      end
