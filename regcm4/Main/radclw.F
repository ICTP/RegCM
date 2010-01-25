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

      subroutine radclw(jslc, ts   ,tnm     ,qnm     ,o3vmr   ,
     $                  pmid    ,pint    ,pmln    ,piln    ,plco2   ,
     $                  plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $                  cld     ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $                  flnsc   ,flntc   ,flwds   ,emiss1d )
C-----------------------------------------------------------------------
C
C Compute longwave radiation heating rates and boundary fluxes
C
C Uses broad band absorptivity/emissivity method to compute clear sky;
C assumes randomly overlapped clouds with variable cloud emissivity to
C include effects of clouds.
C
C Computes clear sky absorptivity/emissivity at lower frequency (in
C general) than the model radiation frequency; uses previously computed
C and stored values for efficiency
C
C Note: This subroutine contains vertical indexing which proceeds
C       from bottom to top rather than the top to bottom indexing
C       used in the rest of the model.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radclw.F,v 1.3 1995/03/03 17:44:43 bonan Exp $
c $Author: bonan $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer jslc
c
C------------------------------Parameters-------------------------------
#     include "parameter.inc"
#     include "include/parrad.cb"
#     include "include/param1.cb"
c
cleexy
#     include "include/param2.cb"
C
C-----------------------------------------------------------------------
      integer plevp2,plevp3,plevp4
      parameter (plevp2=plev+2,plevp3=plev+3,plevp4=plev+4)
C------------------------------Commons----------------------------------
#     include "include/comlun.cb"
#     include "include/comtim.cb"
#     include "include/crdcon.cb"
#     include "include/comctl.cb"
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(kind=8)  ts(plond)            ! Ground (skin) temperature
      real(kind=8)  emiss1d(plond)       ! Emissivity of surface
C
C Input arguments which are only passed to other routines
C
      real(kind=8)  tnm(plond,plev),     ! Level temperature
     $     qnm(plond,plev),        ! Level moisture field
     $     o3vmr(plond,plev),      ! ozone volume mixing ratio
     $     pmid(plond,plev),       ! Level pressure
     $     pint(plond,plevp),      ! Model interface pressure
     $     pmln(plond,plev),       ! Ln(pmid)
     $     piln(plond,plevp),      ! Ln(pint)
     $     plco2(plond,plevp),     ! Path length co2
     $     plh2o(plond,plevp)      ! Path length h2o
      real(kind=8)  n2o(plond,plev),     ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
C
C Input/Output arguments
C
      real(kind=8)  cld(plond,plevp),       ! Cloud cover
     $     tclrsf(plond,plevp)     ! Clear sky fraction
C
C Output arguments
C
      real(kind=8)  qrl(plond,plev),        ! Longwave heating rate
     $     flns(plond),            ! Surface cooling flux
     $     flnt(plond),            ! Net outgoing flux
     $     flnsc(plond),           ! Clear sky surface cooing
     $     flntc(plond),           ! Net clear sky outgoing flux
     $     flwds(plond)            ! Down longwave flux at surface
C
C---------------------------Local variables-----------------------------
C
      integer     i,    ! Longitude index
     $            k,    ! Level index
     $           k1,    ! Level index
     $           k2,    ! Level index
     $           k3,    ! Level index
     $           km,    ! Level index
     $          km1,    ! Level index
     $          km2,    ! Level index
     $          km3,    ! Level index
     $          km4     ! Level index
C
      real(kind=8)  tmp(plond),
     $     tmp1,                ! Temporary 1
     $     absbt(plond)         ! Downward emission at model top
      real(kind=8)  plol(plond,plevp),   ! O3 pressure wghted path length
     $     plos(plond,plevp)    ! O3 path length
C
      real(kind=8)  co2em(plond,plevp),  ! Layer co2 normalized planck funct. derivative
     $     co2eml(plond,plev),  ! Interface co2 normalized planck funct. deriv.
     $     delt(plond),         ! Diff t**4 mid layer to top interface
     $     delt1(plond),        ! Diff t**4 lower intrfc to mid layer
     $     bk1(plond),          ! Absrptvty for vertical quadrature
     $     bk2(plond),          ! Absrptvty for vertical quadrature
     $     ful(plond,plevp),    ! Total upwards longwave flux
     $     fsul(plond,plevp),   ! Clear sky upwards longwave flux
     $     fdl(plond,plevp),    ! Total downwards longwave flux
     $     fsdl(plond,plevp),   ! Clear sky downwards longwv flux
     $     fclb4(plond,plev),   ! Sig t**4 for cld bottom interfc
     $     fclt4(plond,plev),   ! Sig t**4 for cloud top interfc
     $     s(plond,plevp,plevp) ! Flx integral sum
      real(kind=8)  tplnka(plond,plevp), ! Planck fnctn temperature
     $     s2c(plond,plevp),    ! H2o cont amount
     $     s2t(plond,plevp),    ! H2o cont temperature
     $     w(plond,plevp),      ! H2o path
     $     tplnke(plond)        ! Planck fnctn temperature
      real(kind=8)  h2otr(plond,plevp),  ! H2o trnmsn for o3 overlap
     $     co2t(plond,plevp),   ! Prs wghted temperature path
     $     tint(plond,plevp),   ! Interface temperature
     $     tint4(plond,plevp),  ! Interface temperature**4
     $     tlayr(plond,plevp),  ! Level temperature
     $     tlayr4(plond,plevp)  ! Level temperature**4
      real(kind=8)  rtclrsf(plond,plevp),! 1./tclrsf(i,k)
     $     gocp                 ! gravit/cpair
      integer klov(plond),      ! Cloud lowest level index
     $        khiv(plond),      ! Cloud highest level index
     $        khivm(plond)      ! khiv(i) - 1
      integer indx(plond),
     $        nptsc,
     $        ii,
     $        khighest
      logical done(plond),
     $        start(plond)

      real(kind=8)  absnxt,              ! Nearest layer absorptivities
     $     abstot,              ! Non-adjacent layer absorptivites
     $     emstot               ! Total emissivity
#ifdef MPP1
      common /radbuf/ absnxt(plond,plev,4,jxp)
     &      , abstot(plond,plevp,plevp,jxp),emstot(plond,plevp,jxp)
#else
      common /radbuf/ absnxt(plond,plev,4,jlx)
     &      , abstot(plond,plevp,plevp,jlx),emstot(plond,plevp,jlx)
#endif
c
c Trace gas variables
c
      real(kind=8)  ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp),   ! pressure factor for ch4
     $     uptype(plond,plevp)  ! p-type continuum path length
      real(kind=8)  abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
C
C
C------------------------------Externals--------------------------------
C
      integer intmax
      external intmax
c      external radtpl,           ! Compute path lengths
c     $         radems,           ! H2o,co2,o3 emissivity
c     $         radabs            ! H2o,co2,o3 absorptivity
C-----------------------------------------------------------------------
C
C Initialize and recompute the tclrsf array
C
      do i=1,plon
         rtclrsf(i,1) = 1.0/tclrsf(i,1)
      end do
C
      do k=1,plev
         do i=1,plon
            fclb4(i,k) = 0.
            fclt4(i,k) = 0.
            tclrsf(i,k+1) = tclrsf(i,k)*(1. - cld(i,k+1))
            rtclrsf(i,k+1) = 1./tclrsf(i,k+1)
         end do
      end do
C
C Calculate some temperatures needed to derive absorptivity and
C emissivity, as well as some h2o path lengths
C
      call radtpl(tnm     ,ts      ,qnm     ,pint    ,plh2o   ,
     $            tplnka  ,s2c     ,s2t     ,w       ,tplnke  ,
     $            tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    ,
     $            piln    )

c do emissivity and absorptivity calculations
c only if abs/ems computation
c
      IF((jyear.eq.jyear0.and.ktau.eq.0) .or.
     &   (mod(ktau+1,ifrabe).eq.0) ) then

C
C Compute ozone path lengths at frequency of a/e calculation.
C
         call radoz2(o3vmr   ,pint    ,plol    ,plos    )
c
c Compute trace gas path lengths
c
         call trcpth(tnm, pint, cfc11, cfc12, n2o, ch4, qnm,
     $               ucfc11, ucfc12, un2o0,  un2o1,  uch4,
     $               uco211, uco212, uco213, uco221, uco222,
     $               uco223, bn2o0,  bn2o1,  bch4,   uptype)
C
C
C Compute total emissivity:
C
         call radems(s2c     ,s2t     ,w       ,tplnke  ,plh2o   ,
     $               pint    ,plco2   ,tint    ,tint4   ,tlayr   ,
     $               tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  ,
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  ,
     $               co2t    ,h2otr   ,abplnk1 ,abplnk2 ,jslc  )
C
C Compute total absorptivity:
C
         call radabs(pmid    ,pint    ,co2em   ,co2eml  ,tplnka  ,
     $               s2c     ,s2t     ,w       ,h2otr   ,plco2   ,
     $               plh2o   ,co2t    ,tint    ,tlayr   ,plol    ,
     $               plos    ,pmln    ,piln    ,ucfc11  ,ucfc12  ,
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,abplnk1 ,abplnk2 ,
     $               jslc  )

      end if
C
C Find the lowest and highest level cloud for each grid point
C Note: Vertical indexing here proceeds from bottom to top
C
      do i=1,plon
         klov(i) = 0
         done(i) = .false.
      end do
      do k=1,plev
         do i=1,plon
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               klov(i) = k
            end if
         end do
      end do
      call whenne(plon,klov,1,0,indx,nptsc)
      do i=1,plon
         khiv(i) = klov(i)
         done(i) = .false.
      end do
      do k=plev,1,-1
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               khiv(i) = k
            end if
         end do
      end do
      do i=1,plon
         khivm(i) = khiv(i) - 1
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do ii=1,nptsc
         i=indx(ii)
         do k=klov(i),khiv(i)
            fclt4(i,plevp-k) = stebol*tint4(i,plevp2-k)
            fclb4(i,plevp-k) = stebol*tint4(i,plevp3-k)
         end do
      end do
C
C Compute sums used in integrals (all longitude points)
C
C Definition of bk1 & bk2 depends on finite differencing.  for
C trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
C layers only.
C
C delt=t**4 in layer above current sigma level km.
C delt1=t**4 in layer below current sigma level km.
C
      do i=1,plon
         delt(i) = tint4(i,plev) - tlayr4(i,plevp)
         delt1(i) = tlayr4(i,plevp) - tint4(i,plevp)
         s(i,plevp,plevp) = stebol*(delt1(i)*absnxt(i,plev,1,jslc) +
     $                              delt (i)*absnxt(i,plev,4,jslc))
         s(i,plev,plevp)  = stebol*(delt (i)*absnxt(i,plev,2,jslc) +
     $                              delt1(i)*absnxt(i,plev,3,jslc))
      end do
      do k=1,plev-1
         do i=1,plon
            bk2(i) = (abstot(i,k,plev,jslc) +
     &                abstot(i,k,plevp,jslc))*0.5
            bk1(i) = bk2(i)
            s(i,k,plevp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
         end do
      end do
C
C All k, km>1
C
      do km=plev,2,-1
         do i=1,plon
            delt(i)  = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
         end do
         do k=plevp,1,-1
            if (k.eq.km) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,4,jslc)
                  bk1(i) = absnxt(i,km-1,1,jslc)
               end do
            else if(k.eq.km-1) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,2,jslc)
                  bk1(i) = absnxt(i,km-1,3,jslc)
               end do
            else
               do i=1,plon
                  bk2(i) = (abstot(i,k,km-1,jslc) +
     &                      abstot(i,k,km,jslc))*0.5
                  bk1(i) = bk2(i)
               end do
            end if
            do i=1,plon
               s(i,k,km) = s(i,k,km+1) + stebol*
     $                    (bk2(i)*delt(i) + bk1(i)*delt1(i))
            end do
         end do
      end do
C
C Computation of clear sky fluxes always set first level of fsul
C
      do i=1,plon
         if(iemiss.eq.1) then
            fsul(i,plevp) = emiss1d(i)*(stebol*(ts(i)**4))
         else
            fsul(i,plevp) = stebol*(ts(i)**4)
         end if
      end do
C
C Downward clear sky fluxes store intermediate quantities in down flux
C Initialize fluxes to clear sky values.
C
      do i=1,plon
         tmp(i) = fsul(i,plevp) - stebol*tint4(i,plevp)
         fsul(i,1) = fsul(i,plevp)-abstot(i,1,plevp,jslc)*tmp(i)
     &             + s(i,1,2)
         fsdl(i,1) = stebol*(tplnke(i)**4)*emstot(i,1,jslc)
         ful(i,1) = fsul(i,1)
         fdl(i,1) = fsdl(i,1)
      end do
C
C fsdl(i,plevp) assumes isothermal layer
C
      do k=2,plev
         do i=1,plon
            fsul(i,k) = fsul(i,plevp)-abstot(i,k,plevp,jslc)*tmp(i)
     $                + s(i,k,k+1)
            ful(i,k) = fsul(i,k)
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k,jslc) -
     $                  (s(i,k,2) - s(i,k,k+1))
            fdl(i,k) = fsdl(i,k)
         end do
      end do
C
C Store the downward emission from level 1 = total gas emission * sigma
C t**4.  fsdl does not yet include all terms
C
      do i=1,plon
         ful(i,plevp) = fsul(i,plevp)
         absbt(i) = stebol*(tplnke(i)**4)*emstot(i,plevp,jslc)
         fsdl(i,plevp) = absbt(i) - s(i,plevp,2)
         fdl(i,plevp) = fsdl(i,plevp)
      end do
C
C Modifications for clouds
C
C Further qualify longitude subset for computations.  Select only those
C locations where there are clouds (total cloud fraction <= 1.e-3 treated
C as clear)
C
      call whenflt(plon,tclrsf(1,plevp),1,0.999d0,indx,nptsc)
C
C Compute downflux at level 1 for cloudy sky
C
      do ii=1,nptsc
         i=indx(ii)
C
C First clear sky flux plus flux from cloud at level 1
C
         fdl(i,plevp) = fsdl(i,plevp)*tclrsf(i,plev)*
     $         rtclrsf(i,plevp-khiv(i)) + fclb4(i,plev-1)*cld(i,plev)
      end do
C
C Flux emitted by other layers
C Note: Vertical indexing here proceeds from bottom to top
C
      khighest = khiv(intmax(plon,khiv,1))
      do km=3,khighest
         km1 = plevp - km
         km2 = plevp2 - km
         km4 = plevp4 - km
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (km.le.khiv(i)) then
               tmp1 = cld(i,km2)*tclrsf(i,plev)*rtclrsf(i,km2)
               fdl(i,plevp) = fdl(i,plevp) +
     $                       (fclb4(i,km1) - s(i,plevp,km4))*tmp1
            end if
         end do
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do k=1,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.ge.klov(i) .and. k.le.khivm(i)) then
               ful(i,k2) = fsul(i,k2)*(tclrsf(i,plevp)*rtclrsf(i,k1))
            end if
         end do
         do km=1,k
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
CCDIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if (k.le.khivm(i) .and. km.ge.klov(i) .and.
     $             km.le.khivm(i)) then
C
                  ful(i,k2) = ful(i,k2) +
     $                 (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))*
     $                 cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
               end if
            end do
         end do            ! km=1,k
      end do               ! k=1,khighest-1
C
      do k=1,plevp
         k2 = plevp2 - k
         k3 = plevp3 - k
         do i=1,plon
            start(i) = .false.
         end do
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.ge.khiv(i)) then
               start(i) = .true.
               ful(i,k2) = fsul(i,k2)*tclrsf(i,plevp)*
     $                     rtclrsf(i,plevp-khiv(i))
            end if
         end do
         do km=1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
CCDIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if(start(i) .and. km.ge.klov(i) .and. km.le.khiv(i)) then
                  ful(i,k2) = ful(i,k2)  +
     $              (cld(i,km2)*tclrsf(i,km1)*rtclrsf(i,plevp-khiv(i)))*
     $              (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))
               end if
            end do
         end do         ! km=1,khighest
      end do            ! k=1,plevp
C
C Computation of the downward fluxes
C
      do k=2,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.le.khivm(i)) fdl(i,k2) = 0.
         end do
         do km=k+1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km4 = plevp4 - km
CCDIR$ IVDEP
            do ii=1,nptsc
               i=indx(ii)
               if (k.le.khiv(i) .and. km.ge.max0(k+1,klov(i)) .and.
     $             km.le.khiv(i)) then
C
                  fdl(i,k2) = fdl(i,k2) +
     $                 (cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2))*
     $                 (fclb4(i,km1) - s(i,k2,km4) + s(i,k2,k3))
               end if
            end do
         end do            ! km=k+1,khighest
CCDIR$ IVDEP
         do ii=1,nptsc
            i=indx(ii)
            if (k.le.khivm(i)) then
               fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)*
     $              (tclrsf(i,k1)*rtclrsf(i,plevp-khiv(i)))
            end if
         end do
      end do               ! k=1,khighest-1
C
C End cloud modification loops
C
C All longitudes: store history tape quantities
C
      do i=1,plon
C
C Downward longwave flux
C
         flwds(i) = fdl(i,plevp)
C
C Net flux
C
         flns(i) = ful(i,plevp) - fdl(i,plevp)
C
C Clear sky flux at top of atmosphere
C
         flntc(i) = fsul(i,1)
         flnsc(i) = fsul(i,plevp) - fsdl(i,plevp)
C
C Outgoing ir
C
         flnt(i) = ful(i,1) - fdl(i,1)
      end do
C
C Computation of longwave heating (k per sec)
C
      gocp = gravit/cpair
      do k=1,plev
         do i=1,plon
            qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))*
     $                 gocp/((pint(i,k) - pint(i,k+1)))
         end do
      end do
C
      return
      end
