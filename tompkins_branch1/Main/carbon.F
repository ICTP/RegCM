c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      function carbon(vf,t,rm,tg,xlai,xlsai)
      implicit none
      real*8  carbon,vf,t,rm,tg,xlai,xlsai
      real*8  p,rt,alphtl,betatl,cco2,cco2i,xk,b,xkb,al,gt
      real*8  pm,res,wp,wd,w,ac,bc,ab,ccold
      integer it
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ***  a model for whole leaf photosynthesis based on that derived
c   *** by tenhunen as summarized in the text of d gates
c
c     function g represents temperature dependence of photosynthesis
c======================================================================
c  still to improve:
c  parameter for different vegetation types, now wheat only
c  correct initial values for resp must be set in initb
c  parameterisation of resps for cultivated land ok ?
c  betatl = rai/lai  (must be set reasonably)
c====================================================================
c this fn is not in si - vf (radn) and rm (res) handed in are in si
c    result of fn (carbon) is handed back in si (kg/m**2/s)
c======================================================================
      real*8  g,tmx,sl,r,e,xl,a,pml
      g(t,tmx,sl)=dexp(sl*(1./tmx-1./t))
     &           /(1.+(dexp(sl*( 1./tmx-1./t)*6)))*5.e-3*t
c     ***  temperature dependence of dark respiration
      r(t)=dexp(30.-9.e3/t)
c     ****  light dependence of photosynthesis
      e(xl,a,pml)=a*xl/(1.+(a*xl/pml)**2)**0.5

      p=0
      rt=r(t)
c     ****    alphtl=sai/lai     betatl=rai/lai  (rai=root area index)
      alphtl=(xlsai-xlai)/xlai
      betatl=0.5
      if(vf.lt.2.) go to 1
c      **** convert lambda less than 0.7 micron solar into photon units
c    *****  light intensity (e/m2)
      xl=4.6e-3*vf
c    *** co2 external concentration(mm/m3)
c     ****    335ppm/v
      cco2=13.5
c     ***initial guess for co2 concentration inside chloroplast
      cco2i=cco2
c     ****  co2 half max in absence of oxygen  (mm/m3)
      xk=0.5
c    ****  oxygen inhibition factor
      b=3.56
      xkb=xk*b
c     ****  maximum temperature optimum light saturated photosynthesis
      pml=0.050
c     ****    quantum yield(mm/e)
      al=0.05
      gt=g(t,320.d0,4.d3)
c     ****  maximum photosynthesis
      pm=    e(xl,al,pml*gt)
c    ******   iterate
      res=rt

      do 4 it=1,30
c      ****  photorespiration
       wp=pm/(1.+0.4*(1.+cco2i/xk))
c      ****    total respiration
c      ****dark respiration within daytime leaves
       wd=3.e-4*rt+0.14*p
       w=wp+wd
c      *****    carbon uptake factors
       ac=cco2+xkb+rm*(pm-w)
       bc=4.*rm*(cco2*(pm-w)-xkb*w)
       ab = dsqrt(ac**2-bc)
       p = 0.5 * (ac-ab) / rm
       ccold=cco2i
       cco2i=cco2-p*rm
       if(dabs(cco2i-ccold).le.0.05.and.it.gt.9) goto 55
4     continue
55    continue
c     **** respiration outside leaves
      carbon=1.2e-5*((1.-0.14*(alphtl+betatl))*p
     1       -3.0e-4*(alphtl*rt+betatl*r(tg)))
      return

1     continue
c     ***   nighttime maintenance respiration only
      carbon=-0.36e-8*((1.+alphtl)*0.877*rt+betatl*(r(tg)-0.123*rt))

      end
