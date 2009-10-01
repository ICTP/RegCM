c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine ljchrn(chrv,nlshf)
c
c****     routine left justifies a character variable null filling on
c****     the right.  chrv is the character variable.  nlshf tells
c****     how many characters the variable was shifted left
c
      implicit none
      integer nlshf,lnchrv
      character *(*) chrv
      character *1 blank,null
      data blank /' '/
c
      null = char(0)
      lnchrv = len(chrv)
      nlshf = 0
10    continue
      if (chrv(1:1).eq.blank.or.chrv(1:1).eq.null) then
         call lshfch(chrv,lnchrv,1,1)
         nlshf = nlshf+1
         if (nlshf.ne.lnchrv) go to 10
      endif
      return
      end
