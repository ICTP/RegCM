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
 
      subroutine lftemp(iemiss)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     calculate leaf temperature, leaf fluxes, and net transpiration.
!     documented in ncar tech note, dickinson et al., 1986.
!     modifications by klaus blumel, 1988.
!
!        f l o w   d i a g r a m   f o r   l e f t e m
!
!              lftemp ===> stomat
!                          frawat
!                            root
!                           satur
!                          lfdrag
!                          condch
!                          condcq
!                           deriv
!
!     cf = heat transfer coeff. from leaves; assumes same for moisture
!          (see e.g. d. gates' book) for laminar flow past leaf;
!          from ewing paper, cf has dimensions t**1/2 l**-1;
!          so to make it dimensionless, premultiply by 0.01 (si).
!     cgrnd = deriv. of soil energy flux with respect to soil temp.
!                             (used in tgrund)
!     delt = ts-taf
!     ef = transpiration rate
!     efpot = potential evaporation rate (kg/m**2/s)
!     fevpg = evaporative heat flux from ground
!     flnet = (temp gradient)*(d/dt(sig t**4)); hence, t**3 term
!     fseng = sensible heat flux from ground
!     iter = leaf temp iteration counter; runs from 1 to max of 100
!     ra = leaf aerodynamic resistance factor
!     taf = air temperature within foliage canopy
!     tbef = leaf temp at time step before current one
!     ts = air temperature of lowest model layer
!     uaf = mean wind within canopy
!
!     taf1d(n,i) = temperature of air in canopy
!     delt1d(n,i)= difference between temperature of overlying air
!                          and that in canopy
!     delq1d(n,i)= difference between humidity of overlying air
!                          and that in canopy
!
!     convergence of leaf temperature calculation is declared if
!     enough iterations (itmin) and change of temp small enough and
!     change of latent heat fluxes small enough, or if
!     maximum iteration reached (itmax).
!
      use mod_regcm_param
      use mod_param1 , only : dtbat
      use mod_bats
      use mod_ictp01
      use mod_constants , only : tzero , c1es , sigm , wlhv , cpd , ep2
      implicit none
!
! Dummy arguments
!
      integer :: iemiss
      intent (in) iemiss
!
! Local variables
!
      real(8) :: dcn , delmax , efeb , eg1 , epss , fbare , qbare ,     &
               & qcan , qsatdg , rppdry , sf1 , sf2 , sgtg3 , vakb ,    &
               & xxkb
      real(8) , dimension(nnsg,ixm1) :: dels , efpot , tbef
      integer :: iter , itfull , itmax , n , i
!
!=======================================================================
!l    1.   setup information
!=======================================================================
!
!l    1.1  get stress-free stomatal resistance
!     (1st guess at vapor pressure deficit)
      do i = 2 , ixm1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              vpdc(n,i) = 10.
              if ( iemiss.eq.1 ) then
                sgtg3 = emiss_1d(n,i)*(sigm*tg1d(n,i)**3)
              else
                sgtg3 = sigm*tg1d(n,i)**3
              end if
              flneto(n,i) = 4.0*sgtg3*(tlef1d(n,i)-tg1d(n,i))
            end if
          end if
        end do
      end do
      call stomat
!
!l    1.3  determine fraction of total and green canopy surface
!l    covered by water
      call frawat
!
!l    1.4  establish root function in terms of etrc = maximum
!l    sustainable transpiration rate
!     (routine also returns efpr, used in subr. water to
!     define upper soil layer transpiration)
      call root
 
!l    1.5  saturation specific humidity of leaf
      call satur(qsatl,tlef1d,p1d)
 
!=======================================================================
!l    2.   begin iteration for leaf temperature calculation
!=======================================================================
      iter = 0
      efeb = 0.
      delmax = 1.
      itmax = 10
      itfull = itmax
!     itmax = 40
!     itfull = 20
 
      do iter = 0 , itmax
!
!l      2.1  recalc stability dependent canopy & leaf drag coeffs
        if ( iter.eq.0 ) call condch
        call lfdrag
        call condch
 
        do i = 2 , ixm1
          do n = 1 , nnsg
            if ( ldoc1d(n,i).gt.0.5 ) then
              if ( sigf(n,i).gt.0.001 ) then
                ra(n,i) = 1./(cf(n,i)*uaf(n,i))
                cn1(n,i) = wtlh(n,i)*rhs1d(n,i)
                df(n,i) = cn1(n,i)*cpd
 
!l              2.2  decrease foliage conductance for stomatal
!               resistance
                rppdry = ra(n,i)*fdry(n,i)/(rs(n,i)+ra(n,i))
                rpp(n,i) = rppdry + fwet(n,i)
 
!l              2.3  recalculate saturation vapor pressure
                eg1 = eg(n,i)
                eg(n,i) = c1es*dexp(a(n,i)*(tlef1d(n,i)-tzero)/         &
                         & (tlef1d(n,i)-b(n,i)))
                qsatl(n,i) = qsatl(n,i)*eg(n,i)/eg1
              end if
            end if
          end do
        end do
 
!l      2.4  canopy evapotranspiration
        if ( iter.eq.0 ) call condcq
 
        epss = 1.E-10
        do i = 2 , ixm1
          do n = 1 , nnsg
            if ( ldoc1d(n,i).gt.0.5 ) then
              if ( sigf(n,i).gt.0.001 ) then
                efpot(n,i) = cn1(n,i)*(wtgaq(n,i)*qsatl(n,i)-wtgq0(n,i) &
                            & *qg1d(n,i)-wtaq0(n,i)*qs1d(n,i))
 
!as             if(efpot(n,i).ge.0.) then     !if 0 rpp could have
!               floating pt
                if ( efpot(n,i).gt.0. ) then
                  etr(n,i) = efpot(n,i)*ra(n,i)*fdry(n,i)/(rs(n,i)+     &
                           & ra(n,i))
                  rpp(n,i) = dmin1(rpp(n,i),(etr(n,i)+ldew1d(n,i)/      &
                            & dtbat)/efpot(n,i)-epss)
                else
                  etr(n,i) = 0.
                  rpp(n,i) = 1.
                end if
 
                if ( (efpot(n,i).ge.0.) .and. (etr(n,i).ge.etrc(n,i)))  &
                   & then
!*                transpiration demand exceeds supply, stomat adjust
!                 demand
                  rppdry = ra(n,i)*fdry(n,i)/(rs(n,i)+ra(n,i))
                  rppdry = rppdry/(etr(n,i)/etrc(n,i))
                  etr(n,i) = etrc(n,i)
!*                recalculate stomatl resistance and rpp
                  rs(n,i) = ra(n,i)*(fdry(n,i)/rppdry-1.)
                  rpp(n,i) = rppdry + fwet(n,i)
                  rpp(n,i) = dmin1(rpp(n,i),(etr(n,i)+ldew1d(n,i)/      &
                            & dtbat)/efpot(n,i)-epss)
                end if
 
                rppq(n,i) = wlhv*rpp(n,i)
                efe(n,i) = rppq(n,i)*efpot(n,i)
                if ( efe(n,i)*efeb.lt.0.0 ) efe(n,i) = 0.1*efe(n,i)
              end if
            end if
          end do
        end do
!=======================================================================
!l      3.   solve for leaf temperature
!=======================================================================
!l      3.1  update conductances for changes in rpp and cdr
        call condcq
!
!l      3.2  derivatives of energy fluxes with respect to leaf
!l      temperature for newton-raphson calculation of
!l      leaf temperature.
!l      subr.  ii: rs,ra,cdrd,rppq,efe.
!l      subr. output: qsatld,dcd.
        if ( iter.le.itfull ) call deriv
!
!l      3.3  compute dcn from dcd, output from subr. deriv
        do i = 2 , ixm1
          do n = 1 , nnsg
            if ( ldoc1d(n,i).gt.0.5 ) then
              if ( sigf(n,i).gt.0.001 ) then
                dcn = dcd(n,i)*tlef1d(n,i)
!
!l              1.2  radiative forcing for leaf temperature calculation
                if ( iemiss.eq.1 ) then
                  sgtg3 = emiss_1d(n,i)*(sigm*tg1d(n,i)**3)
                else
                  sgtg3 = sigm*tg1d(n,i)**3
                end if
                sf1 = sigf(n,i)*(sabveg(i)-flw1d(i)-(1.-sigf(n,i))*     &
                    & flneto(n,i)+4.0*sgtg3*tg1d(n,i))
                sf2 = 4.*sigf(n,i)*sgtg3 + df(n,i)*wtga(n,i) + dcd(n,i)
 
!l              3.4  iterative leaf temperature calculation
                tbef(n,i) = tlef1d(n,i)
                tlef1d(n,i) = (sf1+df(n,i)*(wta0(n,i)*ts1d(n,i)+        &
                       & wtg0(n,i)*tg1d(n,i))-efe(n,i)+dcn)/sf2
!
!l              3.5  chk magnitude of change; limit to max allowed value
                dels(n,i) = tlef1d(n,i) - tbef(n,i)
                if ( dabs(dels(n,i)).gt.delmax )                        &
                  &   tlef1d(n,i) = tbef(n,i) + delmax*dels(n,i)/       &
                  &   dabs(dels(n,i))
 
!l              3.6  update dependence of stomatal resistance
!l              on vapor pressure deficit
                qcan = wtlq0(n,i)*qsatl(n,i) + qg1d(n,i)*wtgq0(n,i)     &
                     & + qs1d(n,i)*wtaq0(n,i)
                vpdc(n,i) = (1.-rpp(n,i))*(qsatl(n,i)-qcan)*1.E3/ep2
              end if
            end if
          end do
        end do
 
        call stomat
 
!l      3.8  end iteration
 
      end do
 
      do i = 2 , ixm1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!=======================================================================
!l            4.   update dew accumulation (kg/m**2/s)
!=======================================================================
              ldew1d(n,i) = ldew1d(n,i) + (etr(n,i)-efe(n,i)/wlhv)*dtbat
 
!=======================================================================
!l            5.   collect parameters needed to evaluate
!l            sensible and latent fluxes
!=======================================================================
 
!l            5.1  canopy properties
              taf1d(n,i) = wtg0(n,i)*tg1d(n,i) + wta0(n,i)              &
                          & *ts1d(n,i) + wtl0(n,i)*tlef1d(n,i)
              delt1d(n,i) = wtgl(n,i)*ts1d(n,i) - (wtl0(n,i)*           &
                           & tlef1d(n,i)+wtg0(n,i)*tg1d(n,i))
              delq1d(n,i) = wtglq(n,i)*qs1d(n,i) - (wtlq0(n,i)*         &
                           & qsatl(n,i)+wtgq0(n,i)*qg1d(n,i))
              if ( iemiss.eq.1 ) then
                sgtg3 = emiss_1d(n,i)*(sigm*tg1d(n,i)**3)
              else
                sgtg3 = sigm*tg1d(n,i)**3
              end if
              flnet(n,i) = sgtg3*(tlef1d(n,i)-tg1d(n,i))*4.0
              xxkb = dmin1(rough(lveg(n,i)),1.D0)
              vakb = (1.-sigf(n,i))*vspda(n,i) + sigf(n,i)              &
                   & *(xxkb*uaf(n,i)+(1.-xxkb)*vspda(n,i))
              wtg2(n,i) = (1.-sigf(n,i))*cdr(n,i)*vakb
              fbare = wtg2(n,i)*(tg1d(n,i)-ts1d(n,i))
              qbare = wtg2(n,i)*(qg1d(n,i)-qs1d(n,i))
 
!l            5.2  fluxes from soil
              fseng(n,i) = cpd*rhs1d(n,i)*(wtg(n,i)*((wta0(n,i)+        &
                          & wtl0(n,i))*tg1d(n,i)-wta0(n,i)*ts1d(n,i)-   &
                          & wtl0(n,i)*tlef1d(n,i))+fbare)
              fevpg(n,i) = rhs1d(n,i)*rgr(n,i)*(wtg(n,i)*((wtaq0(n,i)+  &
                          & wtlq0(n,i))*qg1d(n,i)-wtaq0(n,i)*qs1d(n,i)- &
                          & wtlq0(n,i)*qsatl(n,i))+qbare)
 
!l            5.3  deriv of soil energy flux with respect to soil temp
              qsatdg = qg1d(n,i)*rgr(n,i)*a(n,i)*(tzero-b(n,i))         &
                     & *(1./(tg1d(n,i)-b(n,i)))**2
              cgrnds(n,i) = rhs1d(n,i)*cpd*(wtg(n,i)*(wta0(n,i)+        &
                     & wtl0(n,i))+wtg2(n,i))
              cgrndl(n,i) = rhs1d(n,i)*qsatdg*((wta(n,i)+wtlq(n,i))*    &
                     & wtg(n,i)*wtsqi(n,i)+wtg2(n,i))
              cgrnd(n,i) = cgrnds(n,i) + cgrndl(n,i)*htvp(n,i)
 
!l            5.4  reinitialize cdrx
!!!           shuttleworth mods #3 removed here !!!!!!
              cdrx(n,i) = cdr(n,i)
!
!l            5.5  fluxes from canopy and soil to overlying air
              fbare = wtg2(n,i)*(tg1d(n,i)-ts1d(n,i))
              qbare = wtg2(n,i)*(qg1d(n,i)-qs1d(n,i))
              sent1d(n,i) = cpd*rhs1d(n,i)*(-wta(n,i)*delt1d(n,i)+fbare)
              evpr1d(n,i) = rhs1d(n,i)*(-wta(n,i)*delq1d(n,i)+          &
                     & rgr(n,i)*qbare)
            end if
          end if
        end do
      end do
 
      end subroutine lftemp
