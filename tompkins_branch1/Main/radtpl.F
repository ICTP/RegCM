
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radtpl(tnm    ,ts     ,qnm    ,pnm    ,plh2o  ,
     $                  tplnka ,s2c    ,s2t    ,w      ,tplnke ,
     $                  tint   ,tint4  ,tlayr  ,tlayr4 ,pmln   ,
     $                  piln   )
C-----------------------------------------------------------------------
C
C Compute temperatures and path lengths for longwave radiation
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radtpl.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
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
C------------------------------Commons----------------------------------
      include 'Commons/crdcon.cb'
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  tnm(plond,plev),     ! Model level temperatures
     $     ts(plond),           ! Surface skin temperature
     $     qnm(plond,plev),     ! Model level specific humidity
     $     pnm(plond,plevp),    ! Pressure at model interfaces (dynes/cm2)
     $     plh2o(plond,plevp)   ! Pressure weighted h2o path
C
C Output arguments
C
      real*8  tplnka(plond,plevp), ! Level temperature from interface temperatures
     $     s2c(plond,plevp),    ! H2o continuum path length
     $     s2t(plond,plevp),    ! H2o tmp and prs wghtd path length
     $     w(plond,plevp),      ! H2o path length
     $     tplnke(plond),       ! Equal to tplnka
     $     tint(plond,plevp),   ! Layer interface temperature
     $     tint4(plond,plevp),  ! Tint to the 4th power
     $     tlayr(plond,plevp),  ! K-1 level temperature
     $     tlayr4(plond,plevp), ! Tlayr to the 4th power
     $     pmln(plond,plev),    ! Ln(pmidm1)
     $     piln(plond,plevp)    ! Ln(pintm1)
C
C---------------------------Local variables-----------------------------
C
      integer   i,              ! Longitude index
     $          k               ! Level index
      real*8  r296,                ! Inverse stand temp for h2o continuum
     $     repsil,              ! Inver ratio mol weight h2o to dry air
     $         dy,              ! Thickness of layer for tmp interp
     $       dpnm,              ! Pressure thickness of layer
     $     dpnmsq,              ! Prs squared difference across layer
     $       rtnm               ! Inverse level temperature
C
C-----------------------------------------------------------------------
C
      r296   = 1./296.d0
      repsil = 1./epsilo
C
C Set the top and bottom intermediate level temperatures,
C top level planck temperature and top layer temp**4.
C
C Tint is lower interface temperature
C (not available for bottom layer, so use ground temperature)
C
      do i=1,plon
         tint(i,plevp)  = ts(i)
         tint4(i,plevp) = tint(i,plevp)**4
         tplnka(i,1)    = tnm(i,1)
         tint(i,1)      = tplnka(i,1)
         tlayr4(i,1)    = tplnka(i,1)**4
         tint4(i,1)     = tlayr4(i,1)
      end do
C
C Intermediate level temperatures are computed using temperature
C at the full level below less dy*delta t,between the full level
C
      do k=2,plev
         do i=1,plon
            dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
            tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
            tint4(i,k) = tint(i,k)**4
         end do
      end do
C
C Now set the layer temp=full level temperatures and establish a
C planck temperature for absorption (tplnka) which is the average
C the intermediate level temperatures.  Note that tplnka is not
C equal to the full level temperatures.
C
      do k=2,plevp
         do i=1,plon
            tlayr(i,k)  = tnm(i,k-1)
            tlayr4(i,k) = tlayr(i,k)**4
            tplnka(i,k) = .5*(tint(i,k) + tint(i,k-1))
         end do
      end do
C
C Calculate tplank for emissivity calculation.
C Assume isothermal tplnke i.e. all levels=ttop.
C
      do i=1,plon
         tplnke(i)  = tplnka(i,1)
         tlayr(i,1) = tint(i,1)
      end do
C
C Now compute h2o path fields:
C
      do i=1,plon
         s2t(i,1) = plh2o(i,1) * tnm(i,1)
C ccm3.2
C        w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
C        s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil

C ccm3.6.6
         w(i,1)   = sslp * (plh2o(i,1)*2.) / pnm(i,1)
         rtnm     = 1./tnm(i,1)
         s2c(i,1) = plh2o(i,1)*
     $            exp(1800.*(rtnm - r296))*qnm(i,1)*repsil
      end do
      do k=1,plev
         do i=1,plon
            dpnm       = pnm(i,k+1) - pnm(i,k)
            dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
            rtnm       = 1./tnm(i,k)
            s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
            w(i,k+1)   = w(i,k)   + rga*qnm(i,k)*dpnm
            s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)*
     $             dexp(1800.*(rtnm - r296))*qnm(i,k)*repsil
         end do
      end do
C
      return
      end
