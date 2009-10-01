
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solar1(xtime)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the solar declination angle from the   c
c     julian date.                                                    c
c                                                                     c
c     xtime  : forecast time in minutes.                              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8  xtime
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param3.cb'
      include 'Commons/date.cb'
#ifdef MPP1
      include 'mpiregcm.h'
#endif
      real*8  pie,dayspy,calday,theta,delta,decdeg
c----------------------------------------------------------------------
c
c-----obecl : obliquity = 23.5 degree.
c
cKN      obecl=23.5*degrad
cKN      sinob=dsin(obecl)
c
c-----calculate longitude of the sun from vernal equinox:
c
cKN      julian=julday+nint((xtime/60.+gmt)/24.)
cKN      if (julian .ge. 81) xlong=dpd*dble(julian-81)
cKN      if (julian .lt. 81) xlong=dpd*dble(julian+284)
cKN      xlong=xlong*degrad
cKN      arg=sinob*dsin(xlong)
cKN      declin=dasin(arg)

cKN   added below
cKN   to take the earth's orbital eccentricity into consideration
C
      pie = 3.141592
!     dayspy = 365.
      dayspy = 365.24
      calday = dble(julday)+(NNNNNN-NSTRT0)/4.
     &       + (xtime/60.+gmt)/24.
      theta = 2.*pie*calday/dayspy
C
C Solar declination in radians:
C
      delta = .006918 - .399912*dcos(theta) + .070257*dsin(theta) -
     $        .006758*dcos(2.*theta) + .000907*dsin(2.*theta) -
     $        .002697*dcos(3.*theta) + .001480*dsin(3.*theta)
C
      declin = delta
C
cKN   added above

      decdeg=declin/degrad
c
#ifdef MPP1
      if (myid .eq. 0) print 10,decdeg
#else
      print 10,decdeg
#endif
10    format(11x,'*** solar declination angle = ',f6.2,' degrees.')
c
      return
      end
