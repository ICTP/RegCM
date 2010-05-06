!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine tgrund
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! 
!     present version of diurnal and seasonal force restore (red 7-88)
!     based on dickinson (1988) force restore paper in j. climate.
!     in particular, shows that a 0.1 m or thicker surface layer will
!     dominate thermal conductivity for surface temperature over the
!     diurnal cycle, and the first 1 m gives appropriate conductivity
!     over annual cycle.
!
!     for snow cover, use weighted average of soil and snow properties.
!     asymptotes to snow values for snow depth large compared to
!     daily and annual temperature waves, respectively.
!
!     the term fct2 provides latent heat for the freezing or
!     thawing of ground over temperatures between -4 and 0 deg c;
!     this range may be changed if the 0.25 in the arithmetic fcn
!     fct1=1/range and the range for subsoil temperature freezing
!     are correspondingly changed.
!
!       tau1 = seconds in a day
!       wlhv = latent heat of vaporization of water
!       wlhs = latent heat of sublimation of water
!       wlhf = latent heat of freezing of water
!       htvp =  wlhv or = wlhs if snow or tg<zero
!     depsnw = depth of snow
!     dtime  = nondimensional diurnal time step
!     dtimea = nondimensional annual time step
!     fsk(x) = soil conductivity fcn (x = v(h2o)/v(soil+pores))
!     fsc(x) = soil heat capacity fcn (x = v(h2o)/v(soil+pores))
!         hs = net energy flux into the surface
!         sg = solar flux absorbed by bare ground
!     skd,ska= diurnal and annual thermal diffusivities
!
      use mod_dynparam
      use mod_bats
      use mod_param1 , only : dtbat
      use mod_constants , only : mathpi , csnw , cws , tau1 , tzero ,   &
                               & wlhf
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,iym1) :: bb , bcoef , cc , depann ,      &
           & depdiu , deprat , fct2 , hs , rscsa , rscsd , ska , skd ,  &
           & sks , swtrta , swtrtd
      real(8) :: bcoefd , bcoefs , c31 , c3t , c41 , c4t , cder , depr ,&
               & depu , dt2 , dtime , dtimea , froze2 , frozen , rscss ,&
               & t3 , tbef , tg , tinc , wtas , wtax , wtd , wtds , x , &
               & xkperi , xnu , xnua
      real(8) :: fct1 , fsc , fsk
      real(8) :: dtbat2 , rdtbat2
      integer :: n , i
 
      fsk(x) = (2.9E-7*x+4.E-9)/(((1.-0.6*x)*x+0.09)*(0.23+x))
      fsc(x) = (0.23+x)*4.186E6
      fct1(x) = wlhf*0.25*1.414/x
 
      dtbat2 = dtbat*2.
      rdtbat2 = 1./dtbat2

!=======================================================================
!l    1.   define thermal conductivity, heat capacity,
!l    and other force restore parameters
!=======================================================================
      xnu = 2.*mathpi/tau1
      xnua = xnu/365.
      dtime = dtbat*xnu
      dtimea = dtbat*xnua
      dt2 = 0.5*dtime
      xkperi = 1.4E-6
 
!l    3.4  permafrost temperature
      t3 = 271.
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
 
!l          1.1  frozen ground values using 44 m2/yr for frozen ground
!l          thermal diffusion coefficient, based on the values of
!l          50 and 38 quoted by osterkamp; ice contribution to
!l          specific heat only o.49 that of water
 
            swtrtd(n,i) = watu(n,i)*porsl(n,i)
            if ( tg1d(n,i).lt.tzero ) then
              frozen = 0.85*dmin1(1.D0,.25*(tzero-tg1d(n,i)))
              skd(n,i) = xkperi
              rscsd(n,i) = fsc(swtrtd(n,i)*(1.-0.51*frozen))
            else
              skd(n,i) = fsk(swtrtd(n,i))*texrat(n,i)
              rscsd(n,i) = fsc(swtrtd(n,i))
            end if
            swtrta(n,i) = watr(n,i)*porsl(n,i)
            if ( tgb1d(n,i).lt.tzero ) then
              froze2 = 0.85*dmin1(1.D0,.25*(tzero-tgb1d(n,i)))
              ska(n,i) = xkperi
              rscsa(n,i) = fsc(swtrta(n,i)*(1.-0.51*froze2))
            else
              ska(n,i) = fsk(swtrta(n,i))*texrat(n,i)
              rscsa(n,i) = fsc(swtrta(n,i))
            end if
 
!l          1.2  correct for snow cover, if significant
            depdiu(n,i) = dsqrt(2.*skd(n,i)/xnu)
            bcoef(n,i) = dtime*depdiu(n,i)/(rscsd(n,i)*skd(n,i))
            if ( scv1d(n,i).gt.1.0 ) then
              wtd = dexp(-2.*scrat(n,i)/depdiu(n,i))
              rscss = csnw*rhosw(n,i)
              sks(n,i) = 7.E-7*cws*rhosw(n,i)
              bcoefs = dsqrt(2.*sks(n,i)/xnu)/(rscss*sks(n,i))
              wtds = (1.-wtd)*scvk(n,i)
              bcoefd = dsqrt(2.*skd(n,i)/xnu)/(rscsd(n,i)*skd(n,i))
              bcoef(n,i) = dtime*(wtds*bcoefs+(1.-wtds)*bcoefd)
              depdiu(n,i) = wtds*dsqrt(2.*sks(n,i)/xnu) + (1.-wtds)     &
                           & *depdiu(n,i)
            end if
            depann(n,i) = dsqrt(2.*ska(n,i)/xnua)
            if ( scv1d(n,i).gt.20. ) then
              wtax = dexp(-2.*scrat(n,i)/depann(n,i))
              wtas = (1.-wtax)*scvk(n,i)
              depann(n,i) = wtas*dsqrt(2.*sks(n,i)/xnua) + (1.-wtas)    &
                           & *depann(n,i)
            end if
            deprat(n,i) = depann(n,i)/depdiu(n,i)
 
!=======================================================================
!l          2.   collect force restore terms
!=======================================================================
            cc(n,i) = 1.0
            fct2(n,i) = 0.
!
!l          2.1  add freezing thermal inertia
            if ( (tg1d(n,i).lt.tzero) .and. (tg1d(n,i).gt.(tzero-4.))   &
               & .and. (sice1d(n,i).le.1.E-22) ) then
              depu = depuv(lveg(n,i))*1.E-3
              cc(n,i) = 1. + dmax1(ssw1d(n,i)-frezu(lveg(n,i)),0.D0)    &
                       & *fct1(depu*rscsd(n,i))
            end if
            if ( (tgb1d(n,i).lt.tzero) .and.                            &
               & (tgb1d(n,i).gt.(tzero-4.)) .and.                       &
               & (sice1d(n,i).le.1.E-22) ) then
              depr = deprv(lveg(n,i))*1.E-3
              fct2(n,i) = dmax1(rsw1d(n,i)-freza(lveg(n,i)),0.D0)       &
                         & *fct1(depr*rscsa(n,i))
            end if
 
!l          2.2  large thermal inertial for permanent ice cap
            if ( lveg(n,i).eq.12 ) fct2(n,i) = 1.E3*fct2(n,i)
 
!l          2.3  collect energy flux terms
            rnet(n,i) = fsw1d(i) - sigf(n,i)*(sabveg(i)-flnet(n,i))     &
                       & - (1.-sigf(n,i))                               &
                       & *(flw1d(i)-sigf(n,i)*flneto(n,i))
            hs(n,i) = rnet(n,i) - fseng(n,i) - fevpg(n,i)*htvp(n,i)
            bb(n,i) = bcoef(n,i)*hs(n,i) + dtime*tgb1d(n,i)
 
!l          2.4  add in snowmelt (melt enough snow to reach freezing
!           temp)
            sm(n,i) = 0.0
            if ( scv1d(n,i).gt.0.0 ) then
              cder = bcoef(n,i)*cgrnd(n,i)
              sm(n,i) = (bb(n,i)+(cc(n,i)-dt2+cder)*tg1d(n,i)-tzero     &
                       & *(cc(n,i)+dt2+cder))/(bcoef(n,i)*wlhf)
!             **********              snow melt always between 0 and
!             total snow
              sm(n,i) = dmax1(0.D0,dmin1(sm(n,i),scv1d(n,i)*2.*         &
                       & rdtbat2))
              bb(n,i) = bb(n,i) - bcoef(n,i)*wlhf*sm(n,i)
            end if
          end if
        end do
      end do
 
!=======================================================================
!l    3.   update soil temperatures
!=======================================================================
!l    3.1  update surface soil temperature
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            tbef = tg1d(n,i)
            cder = bcoef(n,i)*cgrnd(n,i)
            tg = (bb(n,i)+(cc(n,i)-dt2+cder)*tg1d(n,i))/(cc(n,i)+       &
                   & dt2+cder)
            tg1d(n,i) = tg
 
!l          3.2  put brakes on large temperature excursions
            tg1d(n,i) = dmin1(tbef+10.,tg1d(n,i))
            tg1d(n,i) = dmax1(tbef-10.,tg1d(n,i))
 
!l          3.3  correct fluxes to present soil temperature
            tinc = tg1d(n,i) - tbef
            sent1d(n,i) = sent1d(n,i) + tinc*cgrnds(n,i)
            evpr1d(n,i) = evpr1d(n,i) + tinc*cgrndl(n,i)
            fseng(n,i) = fseng(n,i) + tinc*cgrnds(n,i)
            fevpg(n,i) = fevpg(n,i) + tinc*cgrndl(n,i)
!
!l          3.5  couple to deep temperature in permafrost
!l          3.6  update subsoil temperature
            if ( lveg(n,i).eq.9 .or. lveg(n,i).eq.12 ) then
              c31 = 0.5*dtimea*(1.+deprat(n,i))
              c41 = dtimea*deprat(n,i)
              tgb1d(n,i) = ((1.-c31+fct2(n,i))*tgb1d(n,i)+c41*tg1d(n,i)+&
                    & dtimea*t3)/(1.+c31+fct2(n,i))
            else
              c3t = 0.5*dtimea*deprat(n,i)
              c4t = dtimea*deprat(n,i)
              tgb1d(n,i) = ((1.-c3t+fct2(n,i))*tgb1d(n,i)+c4t*tg1d(n,i))&
                     & /(1.+c3t+fct2(n,i))
            end if
          end if
        end do
      end do
 
      end subroutine tgrund
