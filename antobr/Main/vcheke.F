c
c  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
c
      subroutine vcheke (er,ei,nk,numerr,aname)
      implicit none
      integer nk,numerr
      real*8  tol
      parameter (tol=1.e-9)
      character*8 aname
      real*8  er(nk), ei(nk)
      integer numneg,n,nimag
      real*8  emax
c
c  check that eigenvalues are real and positive valued.
c
      numneg=0
      emax=0.
      do 1 n=1,nk
      if (er(n).le.0.) numneg=numneg+1
      if (er(n).gt.emax) emax=er(n)
    1 continue
c
      nimag=0
      do 2 n=1,nk
      if ( ei(n)/emax .gt. tol ) nimag=nimag+1
    2 continue
c
      if ( numneg+nimag .eq. 0 ) return
c
      numerr=numerr+1
      print 3,aname,numneg,nimag
    3 format('0 problem with equivalent depths determined from ',a8,/,
     a10x,i3,' depths are nonpositive valued',10x,i3,' are not real')
c
      return
      end
