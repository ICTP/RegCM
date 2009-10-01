******************BEGIN SUBROUTINE TEMP********************
*             COMPUTES TEMPERATURE PROFILE                *
***********************************************************
      SUBROUTINE TEMP(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de
     a               ,eta,depmax,depth)
      implicit none
      real*8  dt,surf,dz,qe,qh,eta

      real*8  cp
      PARAMETER(cp=4179.98)

      real*8  sw,lnet
      integer depmax, depth
      real*8  t(depmax,2), de(depmax), dnsty(depmax)
      integer k
      real*8  t1,t2,top,bot,tdiff

      surf = 1.0

C******    solve differential equations of heat transfer
      do 10 k= 1,depth
        t(k,2)=t(k,1)
   10 continue

         k = 1
           t1 = sw*( 1.-dexp(-eta*surf) ) / (surf*dnsty(k)*cp)
     a        + (Lnet+qe+qh)             / (surf*dnsty(k)*cp)
           t2 = -de(k)*(t(k,1)-t(k+1,1)) / surf
           t(k,2)=t(k,2)+(t1+t2)*dt

         do 30 k = 2, depth-1
           top = (surf+(k-2)*dz)
           bot = (surf+(k-1)*dz)
           t1 = sw*( dexp(-eta*top)-dexp(-eta*bot) ) / (dz*dnsty(k)*cp)
           t2 = ( de(k-1)*(t(k-1,1)-t(k,1))-de(k)*(t(k,1)-t(k+1,1)))/dz
           t(k,2)=t(k,2)+(t1+t2)*dt
   30    continue

           k=depth
           top = (surf+(k-2)*dz)
           t1 = sw*dexp(-eta*top) / (dz*dnsty(k)*cp)
           t2 = de(k-1)*(t(depth-1,1)-t(depth,1)) / dz
           t(k,2)=t(k,2)+(t1+t2)*dt


       tdiff = 0.
       do 40 k = 1, depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if(k.eq.1) tdiff = tdiff * surf
        t(k,1) = t(k,2)
        dnsty(k)=1000.0*(1.0-1.9549E-05*(dabs((t(k,2)+273.15)
     &          -277.0))**1.68)
   40  continue
csb       print *, 'TEMP: Total temp change = ', tdiff


      return
      end
