cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printsp( array, ix, jx, kx, iinc, jinc, kref, name)
c
      implicit none
      integer ix,jx,kx,iinc,jinc,kref
      real*8  array(ix,jx,kx)
      character*18 cform
      character*08 name
      integer i,j,ineg,inum,jnum,jcheck
c
      ineg = - iinc
      inum = (ix-2) / iinc + 1
      jnum = (jx-2) / jinc + 1
c
      jcheck = 11 * jnum
      if(jcheck.gt.127) print *,'too many values printed.'
c
      write(cform,5) jnum
5     format('(1x,i2,2x,',i2,'g11.3)')
c
      print *,'  field = ', name
      write(6,cform)(i,(array(i,j,kref),j=1,jx-1,jinc),i=ix-1,1,ineg)
c
      return
      end
