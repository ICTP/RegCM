
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cldefr(ioro, t, rel, rei, fice, ps, pmid, jslc)
C-----------------------------------------------------------------------
C
C Compute cloud drop size
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl, January 1993
C
C-----------------------------------------------------------------------
c
c $Id: cldefr.F,v 1.2 1995/02/17 21:28:07 jhack Exp $
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
      include 'Commons/param2.cb'
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ioro(plond)       ! nint(oro(i))
      real*8  t(plond,plev)        ! Temperature
      real*8  ps(plond),           ! surface pressure
     $     pmid(plond,plev)     ! midpoint pressures
      integer jslc
C
C Output arguments
C
      real*8  rel(plond,plev),     ! liquid effective drop size (microns)
     $     rei(plond,plev),     ! ice effective drop size (microns)
     $     fice(plond,plev)     ! fractional ice content within cloud
      real*8  pirnge,              ! nrmlzd pres range for ice particle changes
     $     picemn,              ! normalized pressure below which rei=reimax
     $     rirnge,              ! range of ice radii (reimax - 10 microns)
     $     reimax,              ! maximum ice effective radius
     $     pnrml,               ! normalized pressure
     $     weight               ! coef. for determining rei as fn of P/PS
C
C---------------------------Local workspace-----------------------------
C
      integer i,k               ! longitude, level indices
      real*8  rliq                 ! temporary liquid drop size
      real*8  pi                   ! pi
C
C-----------------------------------------------------------------------
C
      pi = 4. * datan(1.d0)
      do k=1,plev
         do i=1,plon
c
c Define liquid drop size
c
            if(ioro(i).ne.1) then
c
c     Effective liquid radius over ocean and sea ice
c
              rliq = 10.0
            else
c
c     Effective liquid radius over land
c
              rliq = 5.0 +
     $               5.0*dmin1(1.d0,dmax1(0.d0,(263.16-t(i,k))*0.05))
            endif
c
            rel(i,k) = rliq
cfil
c  test radius = 10.0
c           rel(i,k) = 10.0
cfil
C+            rei(i,k) = 30.0
c
c     Determine rei as function of normalized pressure
c
            reimax   = 30.0
            rirnge   = 20.0
            pirnge   = 0.4
            picemn   = 0.4
c
            pnrml    = pmid(i,k)/ps(i)
            weight   = dmax1(dmin1((pnrml-picemn)/pirnge,1.d0),0.d0)
            rei(i,k) = reimax - rirnge*weight
c
c Define fractional amount of cloud that is ice
c
c if warmer than -10 degrees C then water phase
c
             if(t(i,k).gt.263.16) fice(i,k) = 0.0
c
c if colder than -10 degrees C but warmer than -30 C mixed phase
c
             if(t(i,k).le.263.16.and.t(i,k).ge.243.16) then
                   fice(i,k) =(263.16-t(i,k)) / 20.0
             endif
c
c if colder than -30 degrees C then ice phase
c
             if(t(i,k).lt.243.16) fice(i,k) = 1.0
c
c Turn off ice radiative properties by setting fice = 0.0
c
cfil   no-ice test
c             fice(i,k) = 0.0
c
         end do
      end do
C
      return
      end
