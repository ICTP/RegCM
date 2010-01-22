C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc lake model routines cccccccccccccccccccccccccc

      subroutine lake( iutlak, day, dt, ta, ua, za, q
     a         , sw, lnet, hsen, hlat, ts, evap, prec, hice, hsnow)

      implicit none
      integer iutlak
      real(kind=8)  day,dt,ta,ua,za,q,hsen,hlat,ts,evap,prec,hice,hsnow
      real(kind=8)  sw,lnet,ld,lu,kv
      integer depth, freeze

C***        delta: Stephan-Boltzmann const
C***        cp:    heat capacity
C***        Le:    latent heat of evaporation
      real(kind=8)  delta,cp,le
      parameter ( delta=5.67E-08, cp=4179.98, Le=2.47E06 )

      integer depmax
      parameter ( depmax=400 )
      real(kind=8)  t(depmax,2), de(depmax), dnsty(depmax)

      integer lcount,iin,iout,numpts
      common /ldata/ lcount, iin, iout, numpts

      real(kind=8)  surf,dz,zo,z2,u2,hi,eta,tac,tk,qe,qh,ea,tcutoff,hs
      integer ilake,jlake,j,k,kmin

C***        dt:  time step in seconds
C***        surf:surface thickness
C***        dz:  vertical grid spacing in m
C***        kv:  von Karman const
C***        zo:  surface roughness length
      surf = 1.0
      dz = 1.0
      kv= 0.4
      zo= 0.001
      z2= 2.0

C interpolate winds at z1 m to 2m via log wind profile
      u2  = ua*dlog(z2/zo) / dlog(za/zo)

C******    depth: 1-m slices of lake depth
      read(iin) ilake, jlake, depth, freeze, hi, hice, hsnow,
     a       eta, (t(j,1),j=1,depth)
      do k=1,depth
        t(k,2) = t(k,1)
      enddo

      taC= ta - 273.15
      tK=273.15+t(1,1)
      Lu= -0.97*delta*tK**4
      Ld = Lnet - Lu
      Qe = hlat * Le
      Qh = hsen

C convert mixing ratio to air vapor pressure
      ea=q*88.0/(0.622+0.378*q)

c ******    Check if conditions exist for lake ice
      tcutoff = -0.001
      if ((hice .eq. 0.0) .and. (t(1,1) .gt. tcutoff)) then

c ******    Calculate eddy diffusivities
        CALL EDDY (dt,surf,dz,kv,u2,t,dnsty,de,depth)

c ******    Lake temperature calc using BATS sensible and latent heats
        CALL TEMP(dt,surf,dz,t,sw,lnet,Qe,Qh,dnsty,de,eta,depmax,depth)

c ******    Convective mixer
        kmin = 1
        call mixer ( kmin, surf, dz, t, dnsty, depmax, depth )

      else

        call ice(sw,ld,taC,u2,ea,hs,hi,hice,evap,t,depth,prec)
        if (freeze.eq.0) then
          t(1,1)=t(1,2)
        endif

      endif

      write(iout) ilake, jlake, depth, freeze, hi, hice, hsnow, eta
     &          , (t(j,1),j=1,depth)

      write(iutlak) day, ilake, jlake, depth, evap, hi, hice, hsnow
     &            , (t(j,1),j=1,depth)

      ts = t(1,1) + 273.15
      evap = evap / 3600.        !  convert evap from mm/hr to mm/sec
      hice = hice * 1000.        !  convert ice  from m to mm
      hsnow = hsnow  * 100.      !  convert snow from m depth to mm h20

      lcount = lcount + 1
      if(lcount.eq.numpts) then
         lcount = 0
         iin  = 83 - iin
         iout = 83 - iout
         rewind(iin)
         rewind(iout)
      endif

      return
      end
