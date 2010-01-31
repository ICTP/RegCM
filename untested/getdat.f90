!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine getdat(jslc,h2ommr,clat,cld,clwp,coslat,loctim,o3mmr,  &
                      & o3vmr,pilnm1,pintm1,pmidm1,pmlnm1,ps,qm1,tm1,ts)
!
!-----------------------------------------------------------------------
!
! interface routine for column model that both initializes
! certain constants and reads external data:
!
! o3 mass mixing ratios are read in, but the model also requires the
! path lengths; they are computed here
!
! also, from the cloud input (fraction and liquid water path), the
! cloud longwave emissivity must be computed; this is done here
!
      use regcm_param
      use parrad
      use param1
      use param3
      use main
      use rad
      use comtim
      use date
      use crdcon
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(plond) :: clat , coslat , loctim , ps , ts
      real(8) , dimension(plond,plevp) :: cld , pilnm1 , pintm1
      real(8) , dimension(plond,plev) :: clwp , h2ommr , o3mmr , o3vmr ,&
           & pmidm1 , pmlnm1 , qm1 , tm1
      intent (in) jslc
      intent (out) clwp , coslat , loctim , o3vmr , pilnm1 , pmlnm1 ,   &
                 & qm1 , ts
      intent (inout) clat , cld , h2ommr , o3mmr , pintm1 , pmidm1 ,    &
                   & ps , tm1
!
! Local variables
!
      real(8) :: amd , amo , ccvtem , clwtem , gravx , rx , vmmr
      real(8) , dimension(ilx,kx) :: deltaz
      integer :: i , k , kj , n , ncldm1 , nll
      real(8) , dimension(plond) :: rlat
!
!     output arguments
!
! model latitude in radians
! cloud fraction
! cloud liquid water path (g/m**2)
! cosine latitude
! local time of solar computation
! o3 mass mixing ratio
! o3 volume mixing ratio
! ln(pintm1)
! pressure at model interfaces
! pressure at model mid-levels
! ln(pmidm1)
! model surface pressure field
! moisture field
! atmospheric temperature
! surface (air)  temperature
!
!-----------------------------------------------------------------------
!
      data amd/28.9644/
      data amo/48.0000/
!
!KN   instead of reading the data for cloud water and ozone from files,
!KN   their amounts should be calculated in the followings
!
!     set fundamental constants (mks):
!
      gravx = 9.80616
      rx = 287.
!
!     begin read of data:
!-----
!-----surface pressure and scaled pressure, from which level pressures
!-----are computed
      do n = 1 , plond
        ps(n) = (psb(n,jslc)+ptop)*10.
        do nll = 1 , plev
          pmidm1(n,nll) = (psb(n,jslc)*a(nll)+ptop)*10.
!KN       sclpr(nll)=pmidm1(n,nll)/ps(n)
        end do
      end do
!
!.......... convert pressures from mb to pascals and define
!.......... interface pressures:
!
      do i = 1 , plond
        ps(i) = ps(i)*100.
        do k = 1 , plev
!
          pmidm1(i,k) = pmidm1(i,k)*100.
          pmlnm1(i,k) = dlog(pmidm1(i,k))
!
        end do
      end do
      do k = 1 , plevp
        do i = 1 , plond
          pintm1(i,k) = (psb(i,jslc)*sigma(k)+ptop)*1000.
          pilnm1(i,k) = dlog(pintm1(i,k))
        end do
      end do
!
!-----
!-----air temperatures
!-----
      do nll = 1 , plev
        do n = 1 , plond
          tm1(n,nll) = tb(n,nll,jslc)/psb(n,jslc)
        end do
      end do
!-----
!-----surface air temperature
!-----
!-----
!-----h2o mass mixing ratio
!-----
      do nll = 1 , plev
        do n = 1 , plond
          h2ommr(n,nll) = dmax1(1.D-7,qvb(n,nll,jslc)/psb(n,jslc))
          qm1(n,nll) = h2ommr(n,nll)
        end do
      end do
!-----
!-----o3 mass mixing ratio
!-----
      do nll = 1 , plev
        do n = 1 , plond
          kj = plev + 1 - nll
          o3mmr(n,nll) = o3prof(n,kj,jslc)
        end do
      end do
!-----
!-----fractional cloud cover (dependent on relative humidity)
!-----
!qc   = gary's mods for clouds/radiation tie-in to exmois
      do nll = 1 , plev
        do n = 1 , plond
 
          ccvtem = 0.   !cqc mod
!KN       cldfrc(n,nll)=dmax1(cldfra(n,nll)*0.9999999,ccvtem)
          cld(n,nll) = dmax1(cldfra(n,nll)*0.9999999,ccvtem)
!KN       cldfrc(n,nll)=dmin1(cldfrc(n,nll),0.9999999)
          cld(n,nll) = dmin1(cld(n,nll),0.9999999D0)
!
!         implement here the new formula then multiply by 10e6
!qc       if (tm1(n,nll).gt.t0max) clwtem=clwmax
!qc       if (tm1(n,nll).ge.t0st .and. tm1(n,nll).le.t0max)
!qc       1     clwtem=clw0st+((tm1(n,nll)-t0st)/(t0max-t0st))**2
!qc       1     *(clwmax-clw0st)
!qc       if (tm1(n,nll).ge.t0min .and. tm1(n,nll).lt.t0st)
!qc       1     clwtem=clw0st+(tm1(n,nll)-t0st)/(t0min-t0st)
!qc       1     *(clwmin-clw0st)
!qc       if (tm1(n,nll).lt.t0min) clwtem=clwmin
!qc       clwtem=clwtem*1.e6
!
!         convert liquid water content into liquid water path, i.e.
!         multiply b deltaz
          clwtem = cldlwc(n,nll)
                               !cqc mod
          deltaz(n,nll) = rx*tm1(n,nll)*(pintm1(n,nll+1)-pintm1(n,nll)) &
                        & /(gravx*pmidm1(n,nll))
          clwp(n,nll) = clwtem*deltaz(n,nll)
!KN       if (cldfrc(n,nll).eq.0.) clwp(n,nll)=0.
          if ( cld(n,nll).eq.0. ) clwp(n,nll) = 0.
        end do
      end do
 
!     only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
!     do 89 nll=1,plev
!     do 89 n=1,plon
!     if (pintm1(n,nll+1) .lt. 40000. ) then
!     cld(n,nll)=dmin1(cld(n,nll),0.25d0)
!
!     else
!     cld(n,nll)=dmin1(cld(n,nll),0.7d0)
!
!     end if
!     89  continue
 
!
!     set cloud fractional cover at top model level = 0
      do n = 1 , plon
        cld(n,1) = 0.
        clwp(n,1) = 0.
        cld(n,2) = 0.       !yhuang, 8/97 two-level
        clwp(n,2) = 0.
      end do
!
!     set cloud fractional cover at bottom (ncld) model levels = 0
!
      ncldm1 = ncld - 1
      do nll = plev - ncldm1 , plev
        do n = 1 , plon
!KN       cldfrc(n,nll)=0.
          cld(n,nll) = 0.
          clwp(n,nll) = 0.
        end do
      end do
!
!-----
!-----ground temperature
!-----
      do n = 1 , plond
!       tg(n)=tgb(n,jlsc)
!       when using bats calculate an equivalent ground (skin)
!       temperature by averaging over vegetated and non-vegetated areas
!jsp    tg(n)=((1.-vgfrac(n))*tgb(n,jslc)**4.+vgfrac(n)*
!jsp    1   tlef2d(n,jslc)**4.)**0.25
!jsp    tg(n)=tgbb(n,jslc)
        ts(n) = tgbb(n,jslc)
      end do
!
!     cloud cover at surface interface always zero
!
      do i = 1 , plon
!KN     effcld(i,plevp) = 0.
!KN     cldfrc(i,plevp) = 0.
        cld(i,plevp) = 0.
      end do
!
!KN   adopted from regcm2 above
!
!----------------------------------------------------------------------
!
      do i = 1 , plon
!
        do k = 1 , plev
          if ( cld(i,k).gt.0.999 ) cld(i,k) = .999
        end do
!
!KN     added below
        rlat(i) = xlat(i,jslc)
        calday = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
!KN     added above
!
        loctim(i) = (calday-aint(calday))*24.
        clat(i) = rlat(i)*(pie/180.)
        coslat(i) = dcos(clat(i))
!
      end do
!
!     Convert ozone mass mixing ratio to ozone volume mixing ratio:
!
      vmmr = amo/amd
      do k = 1 , plev
        do i = 1 , plon
!         o3mmr(i,k) = vmmr*o3vmr(i,k)
          o3vmr(i,k) = o3mmr(i,k)/vmmr
        end do
      end do
!
      end subroutine getdat
