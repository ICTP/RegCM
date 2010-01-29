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

      subroutine getdat(jslc,h2ommr,
     &     clat,
     &     cld,
     &     clwp,
     &     coslat,
     &     loctim,
     &     o3mmr,
     &     o3vmr,
     &     pilnm1,
     &     pintm1,
     &     pmidm1,
     &     pmlnm1,
     &     ps,
     &     qm1,
     &     tm1,
     &     ts)
c
c-----------------------------------------------------------------------
c
c interface routine for column model that both initializes
c certain constants and reads external data:
c
c o3 mass mixing ratios are read in, but the model also requires the
c path lengths; they are computed here
c
c also, from the cloud input (fraction and liquid water path), the
c cloud longwave emissivity must be computed; this is done here
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer jslc
c
c------------------------------Parameters-------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
C
c------------------------------Commons----------------------------------
c
cKN   added below
c
#     include "include/param1.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/rad.cb"
c
      real(kind=8)  deltaz(ilx,kx)
      real(kind=8)  h2ommr(plond,plev)
c
cKN   added above
#     include "include/comtim.cb"
#     include "include/crdcon.cb"
#     include "include/date.cb"
C
c----------------------------------------------------------------------
c
c output arguments
c
c
      real(kind=8)  clat(plond),                ! model latitude in radians
     &     cld(plond,plevp),    ! cloud fraction
     &     clwp(plond,plev),    ! cloud liquid water path (g/m**2)
     &     coslat(plond),              ! cosine latitude
     &     loctim(plond),       ! local time of solar computation
     &     o3mmr(plond,plev),   ! o3 mass mixing ratio
     &     o3vmr(plond,plev),   ! o3 volume mixing ratio
     &     pilnm1(plond,plevp), ! ln(pintm1)
     &     pintm1(plond,plevp), ! pressure at model interfaces
     &     pmidm1(plond,plev),  ! pressure at model mid-levels
     &     pmlnm1(plond,plev),  ! ln(pmidm1)
     &     ps(plond),           ! model surface pressure field
     &     qm1(plond,plev),     ! moisture field
     &     tm1(plond,plev),     ! atmospheric temperature
     &     ts(plond)            ! surface (air)  temperature
c
c     local workspace
c
      real(kind=8)  rlat(plond)       ! latitude input
c
      integer
     &     i,                   ! longitude index
     &     k                    ! level  index
c
cKN      character*80 label
c
c-----------------------------------------------------------------------
c
cKN   adopted from regcm2 below
c
      real(kind=8)   amd,  ! effective molecular weight of dry air (g/mol)
     &        amo,  ! molecular weight of ozone (g/mol)
     &       vmmr   ! ozone volume mixing ratio
c
      data amd   /  28.9644   /
      data amo   /  48.0000   /
c
cKN   adopted from regcm2 above
c
      real(kind=8)  gravx,rx
      integer n,nll,kj,ncldm1
      real(kind=8)  ccvtem,clwtem
c-----------------------------------------------------------------------
c
cKN   adopted from regcm2 below
c
cKN   instead of reading the data for cloud water and ozone from files,
cKN   their amounts should be calculated in the followings
c
c set fundamental constants (mks):
c
      gravx   =   9.80616
      rx=287.
c
c begin read of data:
c-----
c-----surface pressure and scaled pressure, from which level pressures
c-----are computed
      do 10 n=1,plond
       ps(n)=(psb(n,jslc)+ptop)*10.
      do 10 nll=1,plev
       pmidm1(n,nll)=(psb(n,jslc)*a(nll)+ptop)*10.
cKN       sclpr(nll)=pmidm1(n,nll)/ps(n)
10    continue
c
c.......... convert pressures from mb to pascals and define
c.......... interface pressures:
c
      do 20 i=1,plond
         ps(i) = ps(i) * 100.
         do 20 k=1,plev
c
           pmidm1(i,k) = pmidm1(i,k) * 100.
           pmlnm1(i,k) = dlog(pmidm1(i,k))
c
20       continue
      do 30 k=1,plevp
        do 30 i=1,plond
          pintm1(i,k)=(psb(i,jslc)*sigma(k)+ptop)*1000.
          pilnm1(i,k) = dlog(pintm1(i,k))
30    continue
c
c-----
c-----air temperatures
c-----
      do 40 nll=1,plev
      do 40 n=1,plond
       tm1(n,nll)=tb(n,nll,jslc)/psb(n,jslc)
40    continue
c-----
c-----surface air temperature
c-----
c-----
c-----h2o mass mixing ratio
c-----
      do 60 nll=1,plev
      do 60 n=1,plond
       h2ommr(n,nll)=dmax1(1.d-7,qvb(n,nll,jslc)/psb(n,jslc))
       qm1(n,nll)=h2ommr(n,nll)
60    continue
c-----
c-----o3 mass mixing ratio
c-----
      do 70 nll=1,plev
      do 70 n=1,plond
       kj=plev+1-nll
       o3mmr(n,nll)=o3prof(n,kj,jslc)
70    continue
c-----
c-----fractional cloud cover (dependent on relative humidity)
c-----
cqc = gary's mods for clouds/radiation tie-in to exmois
      do 80 nll=1,plev
      do 80 n=1,plond

       ccvtem = 0.      !cqc mod
cKN       cldfrc(n,nll)=dmax1(cldfra(n,nll)*0.9999999,ccvtem)
       cld(n,nll)=dmax1(cldfra(n,nll)*0.9999999,ccvtem)
cKN       cldfrc(n,nll)=dmin1(cldfrc(n,nll),0.9999999)
       cld(n,nll)=dmin1(cld(n,nll),0.9999999d0)
c
c   implement here the new formula then multiply by 10e6
cqc       if (tm1(n,nll).gt.t0max) clwtem=clwmax
cqc       if (tm1(n,nll).ge.t0st .and. tm1(n,nll).le.t0max)
cqc     1     clwtem=clw0st+((tm1(n,nll)-t0st)/(t0max-t0st))**2
cqc     1     *(clwmax-clw0st)
cqc       if (tm1(n,nll).ge.t0min .and. tm1(n,nll).lt.t0st)
cqc     1     clwtem=clw0st+(tm1(n,nll)-t0st)/(t0min-t0st)
cqc     1     *(clwmin-clw0st)
cqc       if (tm1(n,nll).lt.t0min) clwtem=clwmin
cqc       clwtem=clwtem*1.e6
c
c   convert liquid water content into liquid water path, i.e. multiply b
c   deltaz
       clwtem = cldlwc(n,nll)  !cqc mod
       deltaz(n,nll)=rx*tm1(n,nll)*(pintm1(n,nll+1)-
     1          pintm1(n,nll))/(gravx*pmidm1(n,nll))
       clwp(n,nll)=clwtem*deltaz(n,nll)
cKN       if (cldfrc(n,nll).eq.0.) clwp(n,nll)=0.
       if (cld(n,nll).eq.0.) clwp(n,nll)=0.
80    continue

c  only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
c      do 89 nll=1,plev
c      do 89 n=1,plon
c         if (pintm1(n,nll+1) .lt. 40000. ) then
c           cld(n,nll)=dmin1(cld(n,nll),0.25d0)
c
c         else
c           cld(n,nll)=dmin1(cld(n,nll),0.7d0)
c
c         end if
c  89  continue

c
c     set cloud fractional cover at top model level = 0
      do 90 n=1,plon
         cld(n,1)=0.
         clwp(n,1)=0.
         cld(n,2)=0.        !yhuang, 8/97 two-level
         clwp(n,2)=0.
90    continue
c
c     set cloud fractional cover at bottom (ncld) model levels = 0
c
      ncldm1=ncld-1
      do 100 nll=plev-ncldm1,plev
      do 100 n=1,plon
cKN       cldfrc(n,nll)=0.
       cld(n,nll)=0.
       clwp(n,nll)=0.
100   continue
c
c-----
c-----ground temperature
c-----
      do 101 n=1,plond
c      tg(n)=tgb(n,jlsc)
c   when using bats calculate an equivalent ground (skin) temperature
c   by averaging over vegetated and non-vegetated areas
cjsp       tg(n)=((1.-vgfrac(n))*tgb(n,jslc)**4.+vgfrac(n)*
cjsp     1   tlef2d(n,jslc)**4.)**0.25
cjsp        tg(n)=tgbb(n,jslc)
        ts(n)=tgbb(n,jslc)
101   continue
c
c cloud cover at surface interface always zero
c
      do 500 i=1,plon
cKN         effcld(i,plevp) = 0.
cKN         cldfrc(i,plevp) = 0.
         cld(i,plevp) = 0.
500   continue
c
cKN   adopted from regcm2 above
c
c----------------------------------------------------------------------
c
      do 105 i=1,plon
c
         do 200 k=1,plev
            if( cld(i,k) .gt. 0.999 ) cld(i,k) = .999
 200     continue
c
cKN   added below
      rlat(i) = xlat(i,jslc)
      calday  = dble(julday)+(NNNNNN-NSTRT0)/4.+(xtime/60.+gmt)/24.
cKN   added above
c
      loctim(i) = (calday-aint(calday))*24.
      clat(i)   = rlat(i)*(pie/180.)
      coslat(i) = dcos(clat(i))
C
 105  continue
c
C     Convert ozone mass mixing ratio to ozone volume mixing ratio:
C
      vmmr = amo/amd
      do k=1,plev
         do i=1,plon
c     o3mmr(i,k) = vmmr*o3vmr(i,k)
            o3vmr(i,k) = o3mmr(i,k)/vmmr
         end do
      end do
c
      return
      end subroutine getdat
