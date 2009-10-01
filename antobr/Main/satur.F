
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine satur( qsat, t, p, npin )
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ****  calculates saturation vapor pressure (eg)
c           and qsat = saturated specific humidity (dimensionless)
c
c           uses tetens formula (1930) (ref. riegel,1974,jam,p606
c                                                 equation 1)
c
      implicit none
      integer npin
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'

      real*8  a,b
      common/ab/ a(NNSG,nbmax), b(NNSG,nbmax)
      real*8  qsat(NNSG,nbmax), t(NNSG,nbmax), p(NNSG,nbmax)
      integer n,np

      do 90 np=np1,npts
      do n=1,NNSG
       if(t(n,np)-c(67)) 1, 1, 2
1      a(n,np) = c(70)
       b(n,np) = c(71)
       go to 3
2      a(n,np) = c(72)
       b(n,np) = c(73)
3      continue
       eg(n,np)   = c(74) * dexp( a(n,np)*(t(n,np)-c(67)) 
     &                                   /(t(n,np)-b(n,np)) )
       qsat(n,np) = c(75)*eg(n,np) / (p(n,np)-c(76)*eg(n,np))
      enddo
90    continue

      return
      end
