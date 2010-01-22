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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine lftemp(iemiss)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     calculate leaf temperature, leaf fluxes, and net transpiration.
c     documented in ncar tech note, dickinson et al., 1986.
c     modifications by klaus blumel, 1988.
c
c        f l o w   d i a g r a m   f o r   l e f t e m
c
c              lftemp ===> stomat
c                          frawat
c                            root
c                           satur
c                          lfdrag
c                          condch
c                          condcq
c                           deriv
c
c     cf = heat transfer coeff. from leaves; assumes same for moisture
c          (see e.g. d. gates' book) for laminar flow past leaf;
c          from ewing paper, cf has dimensions t**1/2 l**-1;
c          so to make it dimensionless, premultiply by 0.01 (si).
c     cgrnd = deriv. of soil energy flux with respect to soil temp.
c                             (used in tgrund)
c     delt = ts-taf
c     ef = transpiration rate
c     efpot = potential evaporation rate (kg/m**2/s)
c     fevpg = evaporative heat flux from ground
c     flnet = (temp gradient)*(d/dt(sig t**4)); hence, t**3 term
c     fseng = sensible heat flux from ground
c     iter = leaf temp iteration counter; runs from 1 to max of 100
c     ra = leaf aerodynamic resistance factor
c     taf = air temperature within foliage canopy
c     tbef = leaf temp at time step before current one
c     ts = air temperature of lowest model layer
c     uaf = mean wind within canopy
c
c     taf1d(n,np) = temperature of air in canopy
c     delt1d(n,np)= difference between temperature of overlying air
c                          and that in canopy
c     delq1d(n,np)= difference between humidity of overlying air
c                          and that in canopy
c
c     convergence of leaf temperature calculation is declared if
c     enough iterations (itmin) and change of temp small enough and
c     change of latent heat fluxes small enough, or if
c     maximum iteration reached (itmax).
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"
#     include "include/bats.cb2"

      integer iemiss

      real(kind=8)  a,b
      common /ab/ a(NNSG,nbmax), b(NNSG,nbmax)

      real(kind=8)  cdrd
      Common /aba/ cdrd(NNSG,nbmax)
      real(kind=8)  rppq, efe, qsatld, dcd
      real(kind=8)  vpdc
      Common /abb/ vpdc(NNSG,nbmax)
      Common /abc/ rppq(NNSG,nbmax), efe(NNSG,nbmax)
     &         ,qsatld(NNSG,nbmax), dcd(NNSG,nbmax)
      real(kind=8)  etrc
      Common /abd/ etrc(NNSG,nbmax)

      real(kind=8)  rs
      real(kind=8)  ra
      common /ictp01/ rs(NNSG,nbmax), ra(NNSG,nbmax)

      real(kind=8)  efpot(NNSG,nbmax),tbef(NNSG,nbmax), dels(NNSG,nbmax)
      integer n,np,iter,itmax,itfull,itmin
      real(kind=8)  sgtg3,del,efeb,delmax,dtmin,dlemin,rppdry,eg1,epss
      real(kind=8)  dcn,sf1,sf2,qcan,xxkb,vakb,fbare,qbare,qsatdg
c
c=======================================================================
cl            1.   setup information
c=======================================================================
c
cl            1.1  get stress-free stomatal resistance
c                    (1st guess at vapor pressure deficit)
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         vpdc(n,np) = 10.
         if(iemiss.eq.1) then
            sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
         else
            sgtg3 = c(83)*tg1d(n,np)**3
         endif
         flneto(n,np) = 4.0*sgtg3 * (tlef1d(n,np)-tg1d(n,np))
        endif
       endif
      enddo
90    continue
      call stomat
c
cl            1.3  determine fraction of total and green canopy surface
cl                 covered by water
      call frawat
c
cl            1.4  establish root function in terms of etrc = maximum
cl                 sustainable transpiration rate
c                    (routine also returns efpr, used in subr. water to
c                     define upper soil layer transpiration)
      call root

cl            1.5  saturation specific humidity of leaf
      call satur( qsatl, tlef1d, p1d, npts )

c=======================================================================
cl            2.   begin iteration for leaf temperature calculation
c=======================================================================
      del    = 0.
      iter  = 0
      efeb   = 0.
      delmax = 1.
      itmax  = 10
      itfull = itmax
c     itmax = 40
c     itfull = 20
      itmin  =  2
      dtmin  = 0.01
      dlemin = 0.10

      do 40 iter=0,itmax
c
cl            2.1  recalc stability dependent canopy & leaf drag coeffs
       if (iter.eq.0) call condch
       call lfdrag
       call condch

       do 91 np=np1,npts
       do n=1,NNSG
        if (ldoc1d(n,np).gt.0.5) then
         if (sigf(n,np).gt.0.001) then
          ra(n,np)  = 1. / (cf(n,np)*uaf(n,np))
          cn1(n,np) = wtlh(n,np) * rhs1d(n,np)
          df(n,np)  =  cn1(n,np) * c(58)

cl            2.2  decrease foliage conductance for stomatal resistance
          rppdry = ra(n,np) * fdry(n,np) / (rs(n,np)+ra(n,np))
          rpp(n,np) = rppdry + fwet(n,np)

cl            2.3  recalculate saturation vapor pressure
          eg1 = eg(n,np)
          eg(n,np) = c(74) * dexp( a(n,np)*(tlef1d(n,np)-c(67))
     a       / (tlef1d(n,np)-b(n,np)) )
          qsatl(n,np) = qsatl(n,np) * eg(n,np)/eg1
         endif
        endif
       enddo
91     continue

cl            2.4  canopy evapotranspiration
       if (iter.eq.0) call condcq

       epss = 1.e-10
       do 92 np=np1,npts
       do n=1,NNSG
        if (ldoc1d(n,np).gt.0.5) then
         if (sigf(n,np).gt.0.001) then
         efpot(n,np)=cn1(n,np)*( wtgaq(n,np)*qsatl(n,np)
     &                          -wtgq0(n,np)* qg1d(n,np)
     &                          -wtaq0(n,np)* qs1d(n,np) )

cas      if(efpot(n,np).ge.0.) then     !if 0 rpp could have floating pt
         if(efpot(n,np).gt.0.) then
          etr(n,np)=efpot(n,np)*ra(n,np)*fdry(n,np)/(rs(n,np)+ra(n,np))
          rpp(n,np)=dmin1( rpp(n,np),
     a            (etr(n,np)+ldew1d(n,np)/c(4))/efpot(n,np)-epss )
         else
          etr(n,np) = 0.
          rpp(n,np) = 1.
         endif

         if((efpot(n,np).ge.0.) .and. (etr(n,np).ge.etrc(n,np)))then
c* transpiration demand exceeds supply, stomat adjust demand
          rppdry  = ra(n,np) * fdry(n,np) / (rs(n,np)+ra(n,np))
          rppdry  = rppdry  / ( etr(n,np) / etrc(n,np) )
          etr(n,np) = etrc(n,np)
c* recalculate stomatl resistance and rpp
          rs(n,np)  = ra(n,np) * (fdry(n,np)/rppdry-1.)
          rpp(n,np) = rppdry + fwet(n,np)
          rpp(n,np) = dmin1( rpp(n,np),
     a            (etr(n,np)+ldew1d(n,np)/c(4))/efpot(n,np)-epss )
         endif

         rppq(n,np) = c(125) * rpp(n,np)
         efe(n,np) = rppq(n,np)* efpot(n,np)
         if (efe(n,np)*efeb .lt. 0.0) efe(n,np) = 0.1*efe(n,np)
         endif
        endif
       enddo
92     continue
c=======================================================================
cl            3.   solve for leaf temperature
c=======================================================================
cl            3.1  update conductances for changes in rpp and cdr
       call condcq
c
cl            3.2  derivatives of energy fluxes with respect to leaf
cl                   temperature for newton-raphson calculation of
cl                   leaf temperature.
cl                     subr.  input: rs,ra,cdrd,rppq,efe.
cl                     subr. output: qsatld,dcd.
       if(iter.le.itfull) call deriv
c
cl            3.3  compute dcn from dcd, output from subr. deriv
       do 93 np=np1,npts
       do n=1,NNSG
        if (ldoc1d(n,np).gt.0.5) then
         if (sigf(n,np).gt.0.001) then
          dcn = dcd(n,np)*tlef1d(n,np)
c
cl            1.2  radiative forcing for leaf temperature calculation
          if(iemiss.eq.1) then
             sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
          else
             sgtg3 = c(83)*tg1d(n,np)**3
          endif
          sf1 = sigf(n,np) * ( sabveg(np)-flw1d(np)
     a         - (1.-sigf(n,np))*flneto(n,np)+4.0*sgtg3*tg1d(n,np) )
          sf2 = 4.*sigf(n,np)*sgtg3 + df(n,np)*wtga(n,np) + dcd(n,np)

cl            3.4  iterative leaf temperature calculation
          tbef(n,np) = tlef1d(n,np)
          tlef1d(n,np) = ( sf1 + df(n,np) * ( wta0(n,np)*ts1d(n,np)
     a             + wtg0(n,np)*tg1d(n,np) ) - efe(n,np)+dcn ) / sf2
c
cl            3.5  chk magnitude of change; limit to max allowed value
          dels(n,np) = tlef1d(n,np) - tbef(n,np)
          if(dabs(dels(n,np)).gt.delmax)
     a       tlef1d(n,np)=tbef(n,np)+delmax*dels(n,np)/dabs(dels(n,np))

cl            3.6  update dependence of stomatal resistance
cl                 on vapor pressure deficit
          qcan = wtlq0(n,np)*qsatl(n,np) + qg1d(n,np)*wtgq0(n,np)
     a                           + qs1d(n,np)*wtaq0(n,np)
          vpdc(n,np) = (1.-rpp(n,np)) * (qsatl(n,np)-qcan)*1.e3/c(75)
         endif
        endif
       enddo
93     continue

       call stomat

cl            3.8  end iteration

40    continue

      do 94 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
c=======================================================================
cl            4.   update dew accumulation (kg/m**2/s)
c=======================================================================
         ldew1d(n,np)=ldew1d(n,np)+(etr(n,np)-efe(n,np)/c(125))*c(4)

c=======================================================================
cl            5.   collect parameters needed to evaluate
cl                   sensible and latent fluxes
c=======================================================================

cl            5.1  canopy properties
         taf1d(n,np)=wtg0(n,np)*tg1d(n,np) + wta0(n,np)*ts1d(n,np)
     a                               + wtl0(n,np)*tlef1d(n,np)
         delt1d(n,np) =  wtgl(n,np)*ts1d(n,np)
     a              - (wtl0(n,np)*tlef1d(n,np)+wtg0(n,np)*tg1d(n,np))
         delq1d(n,np) =  wtglq(n,np)*qs1d(n,np)
     a            - (wtlq0(n,np)*qsatl(n,np)+wtgq0(n,np)*qg1d(n,np))
         if(iemiss.eq.1) then
            sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
         else
            sgtg3 = c(83)*tg1d(n,np)**3
         endif
         flnet(n,np) = sgtg3 * (tlef1d(n,np)-tg1d(n,np))*4.0
         xxkb = dmin1( rough(lveg(n,np)), 1.d0 )
         vakb = (1.-sigf(n,np))*vspda(n,np)
     a        +     sigf(n,np) *(xxkb*uaf(n,np)+(1.-xxkb)*vspda(n,np))
         wtg2(n,np) = (1.-sigf(n,np))*cdr(n,np)*vakb
         fbare = wtg2(n,np) * (tg1d(n,np)-ts1d(n,np))
         qbare = wtg2(n,np) * (qg1d(n,np)-qs1d(n,np))

cl            5.2  fluxes from soil
         fseng(n,np) = c(58) * rhs1d(n,np) *
     1          ( wtg(n,np) * ( (wta0(n,np)+wtl0(n,np))*tg1d(n,np) -
     2            wta0(n,np)*ts1d(n,np) - wtl0(n,np)*tlef1d(n,np))
     &          + fbare )
         fevpg(n,np) = rhs1d(n,np) * rgr(n,np) *
     1          ( wtg(n,np) * ( (wtaq0(n,np)+wtlq0(n,np))*qg1d(n,np)-
     2            wtaq0(n,np)*qs1d(n,np) - wtlq0(n,np)*qsatl(n,np))
     &          + qbare )

cl            5.3  deriv of soil energy flux with respect to soil temp
         qsatdg = qg1d(n,np) * rgr(n,np) * a(n,np) * (c(67)-b(n,np))
     a       * ( 1. / (tg1d(n,np)-b(n,np)))**2
         cgrnds(n,np) = rhs1d(n,np) * c(58)
     &           * ( wtg(n,np) * (wta0(n,np)+wtl0(n,np)) + wtg2(n,np) )
         cgrndl(n,np) = rhs1d(n,np) * qsatdg
     &           * ( (wta(n,np)+wtlq(n,np))*wtg(n,np)*wtsqi(n,np)
     &             + wtg2(n,np) )
         cgrnd (n,np) = cgrnds(n,np) + cgrndl(n,np)*htvp(n,np)

cl            5.4  reinitialize cdrx
c!! shuttleworth mods #3 removed here !!!!!!
         cdrx(n,np) = cdr(n,np)
c
cl            5.5  fluxes from canopy and soil to overlying air
         fbare = wtg2(n,np) * (tg1d(n,np)-ts1d(n,np))
         qbare = wtg2(n,np) * (qg1d(n,np)-qs1d(n,np))
         sent1d(n,np)=c(58)*rhs1d(n,np)*(-wta(n,np)*delt1d(n,np)+fbare)
         evpr1d(n,np)=rhs1d(n,np)
     &               *(-wta(n,np)*delq1d(n,np)+rgr(n,np)*qbare)
        endif
       endif
      enddo
94    continue

      return
      end
