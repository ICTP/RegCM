
*****************BEGIN SUBROUTINE MIXER*********************
*              SIMULATES CONVECTIVE MIXING                 *
************************************************************
      SUBROUTINE MIXER ( kmin, surf, dz, T, dnsty, depmax, depth )
      implicit none
      integer kmin
      real*8  surf,dz
      INTEGER depmax, depth
      real*8  T(depmax,2), dnsty(depmax)
      integer k,m,k2
      real*8  avet,avev,vol,tav,tdiff

      do 10 k=kmin,depth-1
        avet=0.0
        avev=0.0

      if (dnsty(k).gt.dnsty(k+1)) then

        do 20 m = kmin, k+1
           if (m .eq. 1) then
              vol = surf
           else
              vol = dz
           endif
           avet=avet+T(m,2)*vol
           avev=avev+vol
   20   continue

        Tav=avet/avev
        do 30 k2 = kmin, k+1
           T(k2,2)=Tav
           dnsty(k2)=1000.0*(1.0-1.9549E-05*(dabs((T(k2,2)+273.15)
     +              -277.0))**1.68)
   30   continue
      endif

   10 continue

      tdiff = 0.
      do 40 k = kmin, depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if(k.eq.1) tdiff = tdiff * surf
        T(k,1)=T(k,2)
   40 continue

      return
      end
