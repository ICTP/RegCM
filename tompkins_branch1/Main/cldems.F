cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cldems(clwp, fice, rei, emis)
C-----------------------------------------------------------------------
C
C Compute cloud emissivity using cloud liquid water path (g/m**2)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, J. Kiehl, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: cldems.F,v 1.2 1995/02/17 21:28:09 jhack Exp $
c $Author: jhack $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
c
      include 'Commons/parrad.cb'
C
      real*8  kabs                   ! longwave absorption coeff (m**2/g)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  clwp(plond,plev),      ! cloud liquid water path (g/m**2)
     $     rei(plond,plev),       ! ice effective drop size (microns)
     $     fice(plond,plev)       ! fractional ice content within cloud
C
C Output arguments
C
      real*8  emis(plond,plev)       ! cloud emissivity (fraction)
C
C---------------------------Local workspace-----------------------------
C
      integer i,k                 ! longitude, level indices
      real*8  kabsl,                 ! longwave absorption coeff (m**2/g)
     $     kabsi                  ! ice absorption coefficient
      parameter (kabsl = 0.090361)
C
C-----------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            kabsi = 0.005 + 1./rei(i,k)
            kabs = kabsl*(1.-fice(i,k))+kabsi*fice(i,k)
            emis(i,k) = 1. - dexp(-1.66*kabs*clwp(i,k))
         end do
      end do
C
      return
      end
