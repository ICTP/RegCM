********************END MAIN PROGRAM************************



*******************BEGIN SUBROUTINE EDDY********************
* COMPUTES DENSITY,EDDY DIFFUSIVITY AND VARIABLE TIME STEP *
************************************************************
      SUBROUTINE EDDY (dt,surf,dz,kv,u2,T,dnsty,de,depth)

      implicit none
      real*8  dm
      parameter(dm=1.38889E-07)
ccc   PARAMETER(dm=5.148e-04,       ! value used in prev simulations

      real*8  dt, surf, dz, kv, u2
      real*8  ks, N2
      integer depth
      real*8  T(depth,2), de(depth), dnsty(depth)
      real*8  demax,rimax,ws,po,dpdz,z,rad,ri
      integer k,kmax

      demax = .5*dz**2 / dt
      demax = .99 * demax
      Rimax=0.0
      do 10 k=1,depth
       dnsty(k)=1000.0*(1.0-1.9549E-05*(dabs((T(k,1)+273.15)
     +          -277.0))**1.68)
   10 continue

      if(u2.lt.0.5) u2 = 0.5

C******     compute eddy diffusion profile
c N2 Brunt-Vaisala number
c Ri gradient Richardson number
C dm molecular diffusion of water
      ks = 0.745 * u2**(-1.84)
      ws = 0.0012 * u2
      Po = 1.0

      do 20 k= 1,depth-1
        dpdz=(dnsty(k+1)-dnsty(k))/dz    ! gtb removed /2.0
        N2 =  dpdz/dnsty(k)*9.8          ! gtb removed minus
        z= surf + dble(k-1)             ! gtb: k was k-1
        rad = 1. + 40. * N2 * (kv*z*dexp(ks*z)/ws)**2
        if (rad .lt. 0) then
          rad=0.0
        endif
        Ri=(-1.0+dsqrt(rad))/20.0
        de(k)=dm+kv*ws*z*Po*dexp(-ks*z)
     +       /(1.0+37.0*Ri**2)
        if (de(k) .gt. demax) then
          de(k) = demax
        endif
        if (dabs(Ri) .gt. Rimax) then
          kmax  = k
          Rimax = dabs(Ri)
        endif
   20 continue
      de(depth)=0.

csb      print *, 'kmax, Rimax = ', kmax, Rimax

      return
      end
