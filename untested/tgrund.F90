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
!     c(125) = latent heat of vaporization of water
!     c(126) = latent heat of sublimation of water
!     c(127) = latent heat of freezing of water
!       htvp = c(125); or = c(126) if snow or tg<zero
!     depsnw = depth of snow
!     dtime  = nondimensional diurnal time step
!     dtimea = nondimensional annual time step
!     fsk(x) = soil conductivity fcn (x = v(h2o)/v(soil+pores))
!     fsc(x) = soil heat capacity fcn (x = v(h2o)/v(soil+pores))
!         hs = net energy flux into the surface
!         sg = solar flux absorbed by bare ground
!     skd,ska= diurnal and annual thermal diffusivities
!
      use mod_regcm_param
      use mod_bats
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,nbmax) :: bb , bcoef , cc , depann ,     &
           & depdiu , deprat , fct2 , hs , rscsa , rscsd , ska , skd ,  &
           & sks , swtrta , swtrtd
      real(8) :: bcoefd , bcoefs , c31 , c3t , c41 , c4t , cder , depr ,&
               & depu , dt2 , dtime , dtimea , froze2 , frozen , rscss ,&
               & t3 , tbef , tg , tinc , wtas , wtax , wtd , wtds , x , &
               & xkperi , xnu , xnua
      real(8) :: fct1 , fsc , fsk
      integer :: n , np
 
      fsk(x) = (2.9E-7*x+4.E-9)/(((1.-0.6*x)*x+0.09)*(0.23+x))
      fsc(x) = (0.23+x)*4.186E6
      fct1(x) = c(127)*0.25*1.414/x
 
!=======================================================================
!l    1.   define thermal conductivity, heat capacity,
!l    and other force restore parameters
!=======================================================================
      xnu = 2.*pi/tau1
      xnua = xnu/365.
      dtime = c(4)*xnu
      dtimea = c(4)*xnua
      dt2 = 0.5*dtime
      xkperi = 1.4E-6
 
!l    3.4  permafrost temperature
      t3 = 271.
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 .and. ldoc1d(n,np).lt.1.5 ) then
 
!l          1.1  frozen ground values using 44 m2/yr for frozen ground
!l          thermal diffusion coefficient, based on the values of
!l          50 and 38 quoted by osterkamp; ice contribution to
!l          specific heat only o.49 that of water
 
            swtrtd(n,np) = watu(n,np)*porsl(n,np)
            if ( tg1d(n,np).lt.c(67) ) then
              frozen = 0.85*dmin1(1.D0,.25*(c(67)-tg1d(n,np)))
              skd(n,np) = xkperi
              rscsd(n,np) = fsc(swtrtd(n,np)*(1.-0.51*frozen))
            else
              skd(n,np) = fsk(swtrtd(n,np))*texrat(n,np)
              rscsd(n,np) = fsc(swtrtd(n,np))
            end if
            swtrta(n,np) = watr(n,np)*porsl(n,np)
            if ( tgb1d(n,np).lt.c(67) ) then
              froze2 = 0.85*dmin1(1.D0,.25*(c(67)-tgb1d(n,np)))
              ska(n,np) = xkperi
              rscsa(n,np) = fsc(swtrta(n,np)*(1.-0.51*froze2))
            else
              ska(n,np) = fsk(swtrta(n,np))*texrat(n,np)
              rscsa(n,np) = fsc(swtrta(n,np))
            end if
 
!l          1.2  correct for snow cover, if significant
            depdiu(n,np) = dsqrt(2.*skd(n,np)/xnu)
            bcoef(n,np) = dtime*depdiu(n,np)/(rscsd(n,np)*skd(n,np))
            if ( scv1d(n,np).gt.1.0 ) then
              wtd = dexp(-2.*scrat(n,np)/depdiu(n,np))
              rscss = csnw*rhosw(n,np)
              sks(n,np) = 7.E-7*cws*rhosw(n,np)
              bcoefs = dsqrt(2.*sks(n,np)/xnu)/(rscss*sks(n,np))
              wtds = (1.-wtd)*scvk(n,np)
              bcoefd = dsqrt(2.*skd(n,np)/xnu)/(rscsd(n,np)*skd(n,np))
              bcoef(n,np) = dtime*(wtds*bcoefs+(1.-wtds)*bcoefd)
              depdiu(n,np) = wtds*dsqrt(2.*sks(n,np)/xnu) + (1.-wtds)   &
                           & *depdiu(n,np)
            end if
            depann(n,np) = dsqrt(2.*ska(n,np)/xnua)
            if ( scv1d(n,np).gt.20. ) then
              wtax = dexp(-2.*scrat(n,np)/depann(n,np))
              wtas = (1.-wtax)*scvk(n,np)
              depann(n,np) = wtas*dsqrt(2.*sks(n,np)/xnua) + (1.-wtas)  &
                           & *depann(n,np)
            end if
            deprat(n,np) = depann(n,np)/depdiu(n,np)
 
!=======================================================================
!l          2.   collect force restore terms
!=======================================================================
            cc(n,np) = 1.0
            fct2(n,np) = 0.
!
!l          2.1  add freezing thermal inertia
            if ( (tg1d(n,np).lt.c(67)) .and. (tg1d(n,np).gt.(c(67)-4.)) &
               & .and. (sice1d(n,np).le.1.E-22) ) then
              depu = depuv(lveg(n,np))*1.E-3
              cc(n,np) = 1. + dmax1(ssw1d(n,np)-frezu(lveg(n,np)),0.D0) &
                       & *fct1(depu*rscsd(n,np))
            end if
            if ( (tgb1d(n,np).lt.c(67)) .and.                           &
               & (tgb1d(n,np).gt.(c(67)-4.)) .and.                      &
               & (sice1d(n,np).le.1.E-22) ) then
              depr = deprv(lveg(n,np))*1.E-3
              fct2(n,np) = dmax1(rsw1d(n,np)-freza(lveg(n,np)),0.D0)    &
                         & *fct1(depr*rscsa(n,np))
            end if
 
!l          2.2  large thermal inertial for permanent ice cap
            if ( lveg(n,np).eq.12 ) fct2(n,np) = 1.E3*fct2(n,np)
 
!l          2.3  collect energy flux terms
            rnet(n,np) = fsw1d(np) - sigf(n,np)*(sabveg(np)-flnet(n,np))&
                       & - (1.-sigf(n,np))                              &
                       & *(flw1d(np)-sigf(n,np)*flneto(n,np))
            hs(n,np) = rnet(n,np) - fseng(n,np) - fevpg(n,np)*htvp(n,np)
            bb(n,np) = bcoef(n,np)*hs(n,np) + dtime*tgb1d(n,np)
 
!l          2.4  add in snowmelt (melt enough snow to reach freezing
!           temp)
            sm(n,np) = 0.0
            if ( scv1d(n,np).gt.0.0 ) then
              cder = bcoef(n,np)*cgrnd(n,np)
              sm(n,np) = (bb(n,np)+(cc(n,np)-dt2+cder)*tg1d(n,np)-c(67) &
                       & *(cc(n,np)+dt2+cder))/(bcoef(n,np)*c(127))
!             **********              snow melt always between 0 and
!             total snow
              sm(n,np) = dmax1(0.D0,dmin1(sm(n,np),scv1d(n,np)*2.*c(7)))
              bb(n,np) = bb(n,np) - bcoef(n,np)*c(127)*sm(n,np)
            end if
          end if
        end do
      end do
 
!=======================================================================
!l    3.   update soil temperatures
!=======================================================================
!l    3.1  update surface soil temperature
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 .and. ldoc1d(n,np).lt.1.5 ) then
            tbef = tg1d(n,np)
            cder = bcoef(n,np)*cgrnd(n,np)
            tg = (bb(n,np)+(cc(n,np)-dt2+cder)*tg1d(n,np))              &
               & /(cc(n,np)+dt2+cder)
            tg1d(n,np) = tg
 
!l          3.2  put brakes on large temperature excursions
            tg1d(n,np) = dmin1(tbef+10.,tg1d(n,np))
            tg1d(n,np) = dmax1(tbef-10.,tg1d(n,np))
 
!l          3.3  correct fluxes to present soil temperature
            tinc = tg1d(n,np) - tbef
            sent1d(n,np) = sent1d(n,np) + tinc*cgrnds(n,np)
            evpr1d(n,np) = evpr1d(n,np) + tinc*cgrndl(n,np)
            fseng(n,np) = fseng(n,np) + tinc*cgrnds(n,np)
            fevpg(n,np) = fevpg(n,np) + tinc*cgrndl(n,np)
!
!l          3.5  couple to deep temperature in permafrost
!l          3.6  update subsoil temperature
            if ( lveg(n,np).eq.9 .or. lveg(n,np).eq.12 ) then
              c31 = 0.5*dtimea*(1.+deprat(n,np))
              c41 = dtimea*deprat(n,np)
              tgb1d(n,np) = ((1.-c31+fct2(n,np))*tgb1d(n,np)+c41*tg1d(n,&
                          & np)+dtimea*t3)/(1.+c31+fct2(n,np))
            else
              c3t = 0.5*dtimea*deprat(n,np)
              c4t = dtimea*deprat(n,np)
              tgb1d(n,np) = ((1.-c3t+fct2(n,np))*tgb1d(n,np)+c4t*tg1d(n,&
                          & np))/(1.+c3t+fct2(n,np))
            end if
          end if
        end do
      end do
 
      end subroutine tgrund
