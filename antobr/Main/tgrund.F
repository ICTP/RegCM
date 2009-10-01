c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine tgrund
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     present version of diurnal and seasonal force restore (red 7-88)
c     based on dickinson (1988) force restore paper in j. climate.
c     in particular, shows that a 0.1 m or thicker surface layer will
c     dominate thermal conductivity for surface temperature over the
c     diurnal cycle, and the first 1 m gives appropriate conductivity
c     over annual cycle.
c
c     for snow cover, use weighted average of soil and snow properties.
c     asymptotes to snow values for snow depth large compared to
c     daily and annual temperature waves, respectively.
c
c     the term fct2 provides latent heat for the freezing or
c     thawing of ground over temperatures between -4 and 0 deg c;
c     this range may be changed if the 0.25 in the arithmetic fcn
c     fct1=1/range and the range for subsoil temperature freezing
c     are correspondingly changed.
c
c       tau1 = seconds in a day
c     c(125) = latent heat of vaporization of water
c     c(126) = latent heat of sublimation of water
c     c(127) = latent heat of freezing of water
c       htvp = c(125); or = c(126) if snow or tg<zero
c     depsnw = depth of snow
c     dtime  = nondimensional diurnal time step
c     dtimea = nondimensional annual time step
c     fsk(x) = soil conductivity fcn (x = v(h2o)/v(soil+pores))
c     fsc(x) = soil heat capacity fcn (x = v(h2o)/v(soil+pores))
c         hs = net energy flux into the surface
c         sg = solar flux absorbed by bare ground
c     skd,ska= diurnal and annual thermal diffusivities
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/bats.cb'
#ifdef MPP1
      include 'Common2/bats.cb2'
#else
      include 'Commons/bats.cb2'
#endif
      real*8  hs(NNSG,nbmax),bb(NNSG,nbmax),cc(NNSG,nbmax)
      real*8  fct2(NNSG,nbmax),bcoef(NNSG,nbmax),deprat(NNSG,nbmax)
      real*8  ska(NNSG,nbmax),swtrta(NNSG,nbmax),rscsa(NNSG,nbmax)
      real*8  skd(NNSG,nbmax),swtrtd(NNSG,nbmax),rscsd(NNSG,nbmax)
      real*8  sks(NNSG,nbmax),depann(NNSG,nbmax),depdiu(NNSG,nbmax)
      real*8  xnu,xnua,dtime,dtimea,dt2,cu,xkperi,t3,frozen
      real*8  froze2,wtd,rscss,bcoefs,wtds,bcoefd,wtax,wtas
      real*8  depu,depr,cder,tbef,tg,tinc,c31,c41,c3t,c4t
      integer n,np

      real*8  x,fsk,fsc,fct1
      fsk(x)=(2.9e-7*x+4.e-9)/(((1.-0.6*x)*x+0.09)*(0.23+x))
      fsc(x)=(0.23+x)*4.186e6
      fct1(x)=c(127)*0.25*1.414/x

c=======================================================================
cl       1.   define thermal conductivity, heat capacity,
cl                 and other force restore parameters
c=======================================================================
      xnu = 2.*pi/tau1
      xnua = xnu/365.
      dtime = c(4)*xnu
      dtimea = c(4)*xnua
      dt2 = 0.5*dtime
      cu=0.
      xkperi=1.4e-6

cl       3.4  permafrost temperature
      t3=271.

      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then

cl       1.1  frozen ground values using 44 m2/yr for frozen ground
cl            thermal diffusion coefficient, based on the values of
cl            50 and 38 quoted by osterkamp; ice contribution to
cl            specific heat only o.49 that of water

        swtrtd(n,np) = watu(n,np) * porsl(n,np)
        if (tg1d(n,np) .lt. c(67)) then
         frozen=0.85*dmin1(1.d0,.25*(c(67)-tg1d(n,np)))
         skd(n,np) = xkperi
         rscsd(n,np) = fsc(swtrtd(n,np)*(1.-0.51*frozen))
        else
         skd(n,np) = fsk(swtrtd(n,np))*texrat(n,np)
         rscsd(n,np)=fsc(swtrtd(n,np))
        endif
        swtrta(n,np) = watr(n,np) * porsl(n,np)
        if (tgb1d(n,np) .lt. c(67)) then
         froze2 = 0.85*dmin1(1.d0,.25*(c(67)-tgb1d(n,np)))
         ska(n,np) = xkperi
         rscsa(n,np) = fsc(swtrta(n,np)*(1.-0.51*froze2))
        else
         ska(n,np) = fsk(swtrta(n,np))*texrat(n,np)
         rscsa(n,np)=fsc(swtrta(n,np))
        endif

cl       1.2  correct for snow cover, if significant
        depdiu(n,np) = dsqrt(2.*skd(n,np)/xnu)
        bcoef(n,np) = dtime*depdiu(n,np)/(rscsd(n,np)*skd(n,np))
        if(scv1d(n,np).gt.1.0) then
         wtd = dexp(-2.*scrat(n,np)/depdiu(n,np))
         rscss = csnw*rhosw(n,np)
         sks(n,np) = 7.e-7*cws*rhosw(n,np)
         bcoefs = dsqrt(2.*sks(n,np)/xnu)/(rscss*sks(n,np))
         wtds = (1.-wtd)*scvk(n,np)
         bcoefd = dsqrt(2.*skd(n,np)/xnu)/(rscsd(n,np)*skd(n,np))
         bcoef(n,np) = dtime * (wtds*bcoefs + (1.-wtds)*bcoefd)
         depdiu(n,np) = wtds*dsqrt(2.*sks(n,np)/xnu)
     &                + (1.-wtds)*depdiu(n,np)
        endif
        depann(n,np) = dsqrt(2.*ska(n,np)/xnua)
        if(scv1d(n,np).gt.20.) then
         wtax = dexp(-2.*scrat(n,np)/depann(n,np))
         wtas = (1.-wtax)*scvk(n,np)
         depann(n,np)= wtas*dsqrt(2.*sks(n,np)/xnua)
     &               + (1.-wtas)*depann(n,np)
        endif
        deprat(n,np) = depann(n,np)/depdiu(n,np)

c=======================================================================
cl       2.   collect force restore terms
c=======================================================================
        cc(n,np) = 1.0
        fct2(n,np) = 0.
c
cl       2.1  add freezing thermal inertia
        if((tg1d(n,np).lt.c(67)).and.(tg1d(n,np).gt.(c(67)-4.))
     &     .and.(sice1d(n,np).le.1.e-22) ) then
         depu = depuv(lveg(n,np))*1.e-3
         cc(n,np) = 1.+dmax1(ssw1d(n,np)-frezu(lveg(n,np)),0.d0)
     a          * fct1(depu*rscsd(n,np))
        endif
        if((tgb1d(n,np).lt.c(67)).and.(tgb1d(n,np).gt.(c(67)-4.))
     &     .and.(sice1d(n,np).le.1.e-22))then
         depr = deprv(lveg(n,np))*1.e-3
         fct2(n,np) =  dmax1(rsw1d(n,np)-freza(lveg(n,np)),0.d0)
     &            * fct1(depr*rscsa(n,np))
        endif

cl       2.2  large thermal inertial for permanent ice cap
        if (lveg(n,np).eq.12) fct2(n,np) = 1.e3*fct2(n,np)

cl       2.3  collect energy flux terms
        rnet(n,np) = fsw1d(np)-sigf(n,np)*(sabveg(np)-flnet(n,np))
     &         - (1.-sigf(n,np))*(flw1d(np)-sigf(n,np)*flneto(n,np))
        hs(n,np) = rnet(n,np) - fseng(n,np) - fevpg(n,np)*htvp(n,np)
        bb(n,np) = bcoef(n,np)*hs(n,np) + dtime*tgb1d(n,np)

cl       2.4  add in snowmelt (melt enough snow to reach freezing temp)
        sm(n,np) = 0.0
        if (scv1d(n,np).gt.0.0) then
         cder = bcoef(n,np)*cgrnd(n,np)
         sm(n,np) = ( bb(n,np) + (cc(n,np)-dt2+cder) * tg1d(n,np)
     &          - c(67) * (cc(n,np)+dt2+cder) ) / (bcoef(n,np)*c(127))
c **********              snow melt always between 0 and total snow
         sm(n,np) = dmax1( 0.d0, dmin1(sm(n,np),scv1d(n,np)*2.*c(7)) )
         bb(n,np) = bb(n,np) - bcoef(n,np)*c(127)*sm(n,np)
        endif
       endif
      enddo
90    continue

c=======================================================================
cl       3.   update soil temperatures
c=======================================================================
cl       3.1  update surface soil temperature
      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5.and.ldoc1d(n,np).lt.1.5) then
        tbef = tg1d(n,np)
        cder = bcoef(n,np)*cgrnd(n,np)
        tg=(bb(n,np)+(cc(n,np)-dt2+cder)*tg1d(n,np))/(cc(n,np)+dt2+cder)
        tg1d(n,np)=tg

cl       3.2  put brakes on large temperature excursions
        tg1d(n,np) = dmin1(tbef+10.,tg1d(n,np))
        tg1d(n,np) = dmax1(tbef-10.,tg1d(n,np))

cl       3.3  correct fluxes to present soil temperature
        tinc = tg1d(n,np) - tbef
        sent1d(n,np) = sent1d(n,np) + tinc*cgrnds(n,np)
        evpr1d(n,np) = evpr1d(n,np) + tinc*cgrndl(n,np)
        fseng(n,np) = fseng(n,np) + tinc*cgrnds(n,np)
        fevpg(n,np) = fevpg(n,np) + tinc*cgrndl(n,np)
c
cl       3.5  couple to deep temperature in permafrost
cl       3.6  update subsoil temperature
        if(lveg(n,np).eq.9.or.lveg(n,np).eq.12) then
         c31 = 0.5*dtimea*(1.+deprat(n,np))
         c41 = dtimea*deprat(n,np)
         tgb1d(n,np) = ( (1.-c31+fct2(n,np))*tgb1d(n,np)
     &                  + c41*tg1d(n,np)
     a               + dtimea*t3 ) / (1.+c31+fct2(n,np))
        else
         c3t = 0.5*dtimea*deprat(n,np)
         c4t = dtimea*deprat(n,np)
         tgb1d(n,np) = ( (1.-c3t+fct2(n,np))*tgb1d(n,np)
     &                  + c4t*tg1d(n,np)) / (1.+c3t+fct2(n,np))
        endif
       endif
      enddo
91    continue

      return
      end
