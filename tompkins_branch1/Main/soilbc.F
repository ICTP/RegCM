
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine soilbc
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   this subrout overwrites many of the soil constants
c   as a function of location(jlon,jlat)
c

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'
      integer n,np,itex
      real*8  phi0,dmax,dmin,dmnor,tweak1,ck

c ================================================================
c   new soils data as a fn of texture make porosity, soil suction,
c   hydraul conduc, wilting frac variables rather than consts
c
c   for explanation of params set here see block data
c
c  relfc is the ratio of field capacity to saturated water content,
c           defined so the rate of gravitational drainage at field
c           capacity is assumed to be 2 mm/day (baver et al., 1972)
c  ===============================================================
c

      do 90 np=np1,npts
      do n=1,NNSG

       if(lveg(n,np).ne.0) then

c **********            lveg is set in subr. interf
        freza(lveg(n,np)) = 0.15*deprv(lveg(n,np))
        frezu(lveg(n,np)) = 0.15*depuv(lveg(n,np))
        itex = iexsol(lveg(n,np))
        texrat(n,np) = skrat(itex)
        porsl(n,np)  = xmopor(itex)
        xkmx(n,np)   = xmohyd(itex)
        bsw(n,np)    = bee(itex)
        bfc(n,np)    = 5.8 - bsw(n,np)
     a           * ( 0.8+0.12*(bsw(n,np)-4.)*dlog10(1.e2*xkmx(n,np)) )

        phi0 = xmosuc(itex)
        dmax = bsw(n,np)*phi0*xkmx(n,np) / porsl(n,np)
        dmin = 1.e-3
        dmnor = 1550.*dmin/dmax
        tweak1 = ( bsw(n,np)*(bsw(n,np)-6.)+10.3 )
     a       / ( bsw(n,np)* bsw(n,np)+40.*bsw(n,np) )
        ck = (1.+dmnor)*tweak1*0.23/0.02356
        evmx0(n,np) = 1.02*dmax*ck/dsqrt(depuv(lveg(n,np))
     &                                  *deprv(lveg(n,np)))
        gwmx0(n,np) = depuv(lveg(n,np)) * porsl(n,np)
        gwmx1(n,np) = deprv(lveg(n,np)) * porsl(n,np)
        gwmx2(n,np) = deptv(lveg(n,np)) * porsl(n,np)
        wiltr(n,np) = xmowil(itex)
c **********            force irrigated crop to be at field capacity
        relfc(n,np) = xmofc(itex)

       endif

      enddo
90    continue

      return
      end
