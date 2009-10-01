c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine snow
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     update snow cover and snow age
c
c     three-part if block:
c       if snow cover < 0, then snow cover and snow age = 0
c       if antarctica, snow age = 0 (katabatic winds keep snow fresh)
c       if elsewhere, snow age follows given formulae
c
c        ps = snow precipitation rate
c     evaps = moisture flux from ground to atmosphere
c        sm = snow melt rate
c     sdrop = snow fallen from vegetation
c
c     aging of snow consists of three factors:
c           age1: snow crystal growth
c           age2: surface melting
c           age3: accumulation  of other particles, soot, etc., which
c                      is small in southern hemisphere
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'

      real*8  sold(NNSG,nbmax)
      integer n,np
      real*8  age3,arg,age1,arg2,age2,tage,dela0,dela,dels,sge

      age3=0.3

c=======================================================================
cl       1.   partition soil evaporation and precipitation
cl                   between water and snow
c=======================================================================
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then

        evapw(n,np) = fevpg(n,np)
        evaps(n,np) = scvk(n,np) * evapw(n,np)
        if(ldoc1d(n,np).gt.1.5) evaps(n,np) = fevpg(n,np)
        evapw(n,np) = (1.-scvk(n,np)) * evapw(n,np)
c
c ******                tm  is temperature of precipitation
        if(tm(n,np) .ge. c(67)) then
         pw(n,np) = prcp1d(n,np) * (1.-sigf(n,np))
         ps(n,np) = 0.0
        else
c ******                snowing
         pw(n,np) = 0.0
         ps(n,np) = prcp1d(n,np) * (1.-sigf(n,np))
        endif
       endif
      enddo
90    continue
c
c=======================================================================
cl       2.   update snow cover
c=======================================================================
      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        sold(n,np) = scv1d(n,np)
        scv1d(n,np) = scv1d(n,np)
     a            + c(4)*(ps(n,np)-evaps(n,np)-sm(n,np)) + sdrop(n,np)
        scv1d(n,np) = dmax1( scv1d(n,np), 0.d0 )
        sag1d(n,np) = dmax1( sag1d(n,np), 0.d0 )

c ******           snow cover except for antarctica
c=======================================================================
cl       3.   increment non-dimensional "age" of snow;
cl            10 mm snow restores surface to that of new snow.
c=======================================================================
        if(scv1d(n,np).gt.0.) then
         arg=5.e3*(1./c(67)-1./tg1d(n,np))
         age1=dexp(arg)
         arg2=dmin1(0.d0,10.*arg)
         age2=dexp(arg2)
         tage=age1+age2+age3
         dela0=1.e-6*c(4)
         dela=dela0*tage
         dels=0.1*dmax1( 0.d0, scv1d(n,np)-sold(n,np) )
         sge=(sag1d(n,np)+dela)*(1.0-dels)
         sag1d(n,np) = dmax1( 0.d0, sge )
        endif

c ******           antarctica
        if (scv1d(n,np).gt.800.) sag1d(n,np)=0.
       endif
      enddo
91    continue

      return
      end
