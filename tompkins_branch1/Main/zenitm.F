c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine zenitm(coszrs,ivmx,jslc)
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   this subroutine calculates the cosine of the solar zenith angle
c   for all longitude points of the mm42d domain. it needs as inputs
c   the longitude and latitude of the points, the initial date of the
c   simulation and the gmt. all these quantities are specified
c   in the initialization procedure of mm42d
c
      implicit none
      integer ivmx,jslc
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param3.cb'
#ifdef MPP1
      include 'Common2/main.cb'
#else
      include 'Commons/main.cb'
#endif
      include 'Commons/date.cb'
c
      real*8  coszrs(ix)
      integer ill
      real*8  xt24,tlocap,omega,xxlat
c
      xt24 = dmod(lhour*60.+xtime,1440.d0)
      do 10 ill = 1,ivmx
       tlocap = xt24/60. + xlong(ill,jslc)/15.
       tlocap = dmod(tlocap+24., 24.d0)
       omega  = 15.*(tlocap-12.)*degrad
       xxlat  = xlat(ill,jslc)*degrad
c      coszrs = cosine of solar zenith angle
       coszrs(ill)=dsin(declin)*dsin(xxlat) +
     1            dcos(declin)*dcos(xxlat)*dcos(omega)
       coszrs(ill) = dmax1(0.d0,coszrs(ill))
10    continue

      return
      end
