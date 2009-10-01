
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine resetr(pa      ,kdim    ,pvalue  )
C-----------------------------------------------------------------------
C
C Reset array pa(kdim) to pvalue
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Bath, June 1992
C
C-----------------------------------------------------------------------
c
c $Id: resetr.F,v 1.1.1.1 1995/02/09 23:27:04 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer kdim      ! Dimension of array pa
      real*8  pvalue       ! Value to store in pa
C
C Output arguments
C
      real*8  pa(kdim)     ! Array to reset
C
C---------------------------Local variable------------------------------
C
      integer j         ! Loop index
C
C-----------------------------------------------------------------------
C
      do j=1,kdim
         pa(j) = pvalue
      end do
C
      return
      end
