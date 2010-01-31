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
!     taf1d(n,np) = temperature of air in canopy
!     delt1d(n,np)= difference between temperature of overlying air
!                          and that in canopy
!     delq1d(n,np)= difference between humidity of overlying air
!                          and that in canopy
!
!     convergence of leaf temperature calculation is declared if
!     enough iterations (itmin) and change of temp small enough and
!     change of latent heat fluxes small enough, or if
!     maximum iteration reached (itmax).
!
      use regcm_param
      use bats
      use ictp01
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
      real(8) , dimension(nnsg,nbmax) :: dels , efpot , tbef
      integer :: iter , itfull , itmax , n , np
!
!=======================================================================
!l    1.   setup information
!=======================================================================
!
!l    1.1  get stress-free stomatal resistance
!     (1st guess at vapor pressure deficit)
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              vpdc(n,np) = 10.
              if ( iemiss.eq.1 ) then
                sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
              else
                sgtg3 = c(83)*tg1d(n,np)**3
              end if
              flneto(n,np) = 4.0*sgtg3*(tlef1d(n,np)-tg1d(n,np))
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
      call satur(qsatl,tlef1d,p1d,npts)
 
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
 
        do np = np1 , npts
          do n = 1 , nnsg
            if ( ldoc1d(n,np).gt.0.5 ) then
              if ( sigf(n,np).gt.0.001 ) then
                ra(n,np) = 1./(cf(n,np)*uaf(n,np))
                cn1(n,np) = wtlh(n,np)*rhs1d(n,np)
                df(n,np) = cn1(n,np)*c(58)
 
!l              2.2  decrease foliage conductance for stomatal
!               resistance
                rppdry = ra(n,np)*fdry(n,np)/(rs(n,np)+ra(n,np))
                rpp(n,np) = rppdry + fwet(n,np)
 
!l              2.3  recalculate saturation vapor pressure
                eg1 = eg(n,np)
                eg(n,np) = c(74)                                        &
                         & *dexp(a(n,np)*(tlef1d(n,np)-c(67))/(tlef1d(n,&
                         & np)-b(n,np)))
                qsatl(n,np) = qsatl(n,np)*eg(n,np)/eg1
              end if
            end if
          end do
        end do
 
!l      2.4  canopy evapotranspiration
        if ( iter.eq.0 ) call condcq
 
        epss = 1.E-10
        do np = np1 , npts
          do n = 1 , nnsg
            if ( ldoc1d(n,np).gt.0.5 ) then
              if ( sigf(n,np).gt.0.001 ) then
                efpot(n,np) = cn1(n,np)                                 &
                            & *(wtgaq(n,np)*qsatl(n,np)-wtgq0(n,np)     &
                            & *qg1d(n,np)-wtaq0(n,np)*qs1d(n,np))
 
!as             if(efpot(n,np).ge.0.) then     !if 0 rpp could have
!               floating pt
                if ( efpot(n,np).gt.0. ) then
                  etr(n,np) = efpot(n,np)*ra(n,np)*fdry(n,np)           &
                            & /(rs(n,np)+ra(n,np))
                  rpp(n,np) = dmin1(rpp(n,np),(etr(n,np)+ldew1d(n,np)/c(&
                            & 4))/efpot(n,np)-epss)
                else
                  etr(n,np) = 0.
                  rpp(n,np) = 1.
                end if
 
                if ( (efpot(n,np).ge.0.) .and. (etr(n,np).ge.etrc(n,np))&
                   & ) then
!*                transpiration demand exceeds supply, stomat adjust
!                 demand
                  rppdry = ra(n,np)*fdry(n,np)/(rs(n,np)+ra(n,np))
                  rppdry = rppdry/(etr(n,np)/etrc(n,np))
                  etr(n,np) = etrc(n,np)
!*                recalculate stomatl resistance and rpp
                  rs(n,np) = ra(n,np)*(fdry(n,np)/rppdry-1.)
                  rpp(n,np) = rppdry + fwet(n,np)
                  rpp(n,np) = dmin1(rpp(n,np),(etr(n,np)+ldew1d(n,np)/c(&
                            & 4))/efpot(n,np)-epss)
                end if
 
                rppq(n,np) = c(125)*rpp(n,np)
                efe(n,np) = rppq(n,np)*efpot(n,np)
                if ( efe(n,np)*efeb.lt.0.0 ) efe(n,np) = 0.1*efe(n,np)
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
!l      subr.  input: rs,ra,cdrd,rppq,efe.
!l      subr. output: qsatld,dcd.
        if ( iter.le.itfull ) call deriv
!
!l      3.3  compute dcn from dcd, output from subr. deriv
        do np = np1 , npts
          do n = 1 , nnsg
            if ( ldoc1d(n,np).gt.0.5 ) then
              if ( sigf(n,np).gt.0.001 ) then
                dcn = dcd(n,np)*tlef1d(n,np)
!
!l              1.2  radiative forcing for leaf temperature calculation
                if ( iemiss.eq.1 ) then
                  sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
                else
                  sgtg3 = c(83)*tg1d(n,np)**3
                end if
                sf1 = sigf(n,np)                                        &
                    & *(sabveg(np)-flw1d(np)-(1.-sigf(n,np))*flneto(n,  &
                    & np)+4.0*sgtg3*tg1d(n,np))
                sf2 = 4.*sigf(n,np)*sgtg3 + df(n,np)*wtga(n,np)         &
                    & + dcd(n,np)
 
!l              3.4  iterative leaf temperature calculation
                tbef(n,np) = tlef1d(n,np)
                tlef1d(n,np) = (sf1+df(n,np)*(wta0(n,np)*ts1d(n,np)+wtg0&
                             & (n,np)*tg1d(n,np))-efe(n,np)+dcn)/sf2
!
!l              3.5  chk magnitude of change; limit to max allowed value
                dels(n,np) = tlef1d(n,np) - tbef(n,np)
                if ( dabs(dels(n,np)).gt.delmax ) tlef1d(n,np)          &
                   & = tbef(n,np) + delmax*dels(n,np)/dabs(dels(n,np))
 
!l              3.6  update dependence of stomatal resistance
!l              on vapor pressure deficit
                qcan = wtlq0(n,np)*qsatl(n,np) + qg1d(n,np)*wtgq0(n,np) &
                     & + qs1d(n,np)*wtaq0(n,np)
                vpdc(n,np) = (1.-rpp(n,np))*(qsatl(n,np)-qcan)          &
                           & *1.E3/c(75)
              end if
            end if
          end do
        end do
 
        call stomat
 
!l      3.8  end iteration
 
      end do
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
!=======================================================================
!l            4.   update dew accumulation (kg/m**2/s)
!=======================================================================
              ldew1d(n,np) = ldew1d(n,np) + (etr(n,np)-efe(n,np)/c(125))&
                           & *c(4)
 
!=======================================================================
!l            5.   collect parameters needed to evaluate
!l            sensible and latent fluxes
!=======================================================================
 
!l            5.1  canopy properties
              taf1d(n,np) = wtg0(n,np)*tg1d(n,np) + wta0(n,np)          &
                          & *ts1d(n,np) + wtl0(n,np)*tlef1d(n,np)
              delt1d(n,np) = wtgl(n,np)*ts1d(n,np)                      &
                           & - (wtl0(n,np)*tlef1d(n,np)+wtg0(n,np)      &
                           & *tg1d(n,np))
              delq1d(n,np) = wtglq(n,np)*qs1d(n,np)                     &
                           & - (wtlq0(n,np)*qsatl(n,np)+wtgq0(n,np)     &
                           & *qg1d(n,np))
              if ( iemiss.eq.1 ) then
                sgtg3 = emiss_1d(n,np)*(c(83)*tg1d(n,np)**3)
              else
                sgtg3 = c(83)*tg1d(n,np)**3
              end if
              flnet(n,np) = sgtg3*(tlef1d(n,np)-tg1d(n,np))*4.0
              xxkb = dmin1(rough(lveg(n,np)),1.D0)
              vakb = (1.-sigf(n,np))*vspda(n,np) + sigf(n,np)           &
                   & *(xxkb*uaf(n,np)+(1.-xxkb)*vspda(n,np))
              wtg2(n,np) = (1.-sigf(n,np))*cdr(n,np)*vakb
              fbare = wtg2(n,np)*(tg1d(n,np)-ts1d(n,np))
              qbare = wtg2(n,np)*(qg1d(n,np)-qs1d(n,np))
 
!l            5.2  fluxes from soil
              fseng(n,np) = c(58)*rhs1d(n,np)                           &
                          & *(wtg(n,np)*((wta0(n,np)+wtl0(n,np))        &
                          & *tg1d(n,np)-wta0(n,np)*ts1d(n,np)-wtl0(n,np)&
                          & *tlef1d(n,np))+fbare)
              fevpg(n,np) = rhs1d(n,np)*rgr(n,np)                       &
                          & *(wtg(n,np)*((wtaq0(n,np)+wtlq0(n,np))      &
                          & *qg1d(n,np)-wtaq0(n,np)*qs1d(n,np)          &
                          & -wtlq0(n,np)*qsatl(n,np))+qbare)
 
!l            5.3  deriv of soil energy flux with respect to soil temp
              qsatdg = qg1d(n,np)*rgr(n,np)*a(n,np)*(c(67)-b(n,np))     &
                     & *(1./(tg1d(n,np)-b(n,np)))**2
              cgrnds(n,np) = rhs1d(n,np)*c(58)                          &
                           & *(wtg(n,np)*(wta0(n,np)+wtl0(n,np))        &
                           & +wtg2(n,np))
              cgrndl(n,np) = rhs1d(n,np)                                &
                           & *qsatdg*((wta(n,np)+wtlq(n,np))*wtg(n,np)  &
                           & *wtsqi(n,np)+wtg2(n,np))
              cgrnd(n,np) = cgrnds(n,np) + cgrndl(n,np)*htvp(n,np)
 
!l            5.4  reinitialize cdrx
!!!           shuttleworth mods #3 removed here !!!!!!
              cdrx(n,np) = cdr(n,np)
!
!l            5.5  fluxes from canopy and soil to overlying air
              fbare = wtg2(n,np)*(tg1d(n,np)-ts1d(n,np))
              qbare = wtg2(n,np)*(qg1d(n,np)-qs1d(n,np))
              sent1d(n,np) = c(58)*rhs1d(n,np)                          &
                           & *(-wta(n,np)*delt1d(n,np)+fbare)
              evpr1d(n,np) = rhs1d(n,np)                                &
                           & *(-wta(n,np)*delq1d(n,np)+rgr(n,np)*qbare)
            end if
          end if
        end do
      end do
 
      end subroutine lftemp
