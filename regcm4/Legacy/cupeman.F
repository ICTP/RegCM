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

      subroutine cupeman
     *    (t,   q,    qs,     u,    v,      tra,    p,    ph,
     *     nd,  na,   nl,   ntra,   delt, iflag,  ft,     fq,   fu,
     *     fv,  ftra, precip, wd,   tprime, qprime, cbmf, icb, inb )
c
***************************************************************************
*****                       subroutine cupeman (formerly convect)     *****
*****                          version 4.3c                           *****
*****                          20 may, 2002                           *****
*****                          kerry emanuel                          *****
***************************************************************************
c
c-----------------------------------------------------------------------------
c    *** on input:      ***
c
c     t:   array of absolute temperature (k) of dimension nd, with first
c           index corresponding to lowest model level. note that this array
c           will be altered by the subroutine if dry convective adjustment
c           occurs and if ipbl is not equal to 0.
c
c     q:   array of specific humidity (gm/gm) of dimension nd, with first
c            index corresponding to lowest model level. must be defined
c            at same grid levels as t. note that this array will be altered
c            if dry convective adjustment occurs and if ipbl is not equal to 0.
c
c     qs:  array of saturation specific humidity of dimension nd, with first
c            index corresponding to lowest model level. must be defined
c            at same grid levels as t. note that this array will be altered
c            if dry convective adjustment occurs and if ipbl is not equal to 0.
c
c     u:   array of zonal wind velocity (m/s) of dimension nd, witth first
c            index corresponding with the lowest model level. defined at
c            same levels as t. note that this array will be altered if
c            dry convective adjustment occurs and if ipbl is not equal to 0.
c
c     v:   same as u but for meridional velocity.
c
c     tra: array of passive tracer mixing ratio, of dimensions (nd,ntra),
c            where ntra is the number of different tracers. if no
c            convective tracer transport is needed, define a dummy
c            input array of dimension (nd,1). tracers are defined at
c            same vertical levels as t. note that this array will be altered
c            if dry convective adjustment occurs and if ipbl is not equal to 0.
c
c     p:   array of pressure (mb) of dimension nd, with first
c            index corresponding to lowest model level. must be defined
c            at same grid levels as t.
c
c     ph:  array of pressure (mb) of dimension nd+1, with first index
c            corresponding to lowest level. these pressures are defined at
c            levels intermediate between those of p, t, q and qs. the first
c            value of ph should be greater than (i.e. at a lower level than)
c            the first value of the array p.
c
c     nd:  the dimension of the arrays t,q,qs,p,ph,ft and fq
c
c     nl:  the maximum number of levels to which convection can
c            penetrate, plus 1.
c            nl must be less than or equal to nd-1.
c
c     ntra:the number of different tracers. if no tracer transport
c            is needed, set this equal to 1. (on most compilers, setting
c            ntra to 0 will bypass tracer calculation, saving some cpu.)  
c
c     delt: the model time step (sec) between calls to convect
c
c----------------------------------------------------------------------------
c    ***   on output:         ***
c
c     iflag: an output integer whose value denotes the following:
c
c                value                        interpretation
c                -----                        --------------
c                  0               no moist convection; atmosphere is not
c                                  unstable, or surface temperature is less
c                                  than 250 k or surface specific humidity
c                                  is non-positive.
c
c                  1               moist convection occurs.
c
c                  2               no moist convection: lifted condensation
c                                  level is above the 200 mb level.
c
c                  3               no moist convection: cloud base is higher
c                                  then the level nl-1.
c
c                  4               moist convection occurs, but a cfl condition
c                                  on the subsidence warming is violated. this
c                                  does not cause the scheme to terminate.
c
c     ft:   array of temperature tendency (k/s) of dimension nd, defined at same
c             grid levels as t, q, qs and p.
c
c     fq:   array of specific humidity tendencies ((gm/gm)/s) of dimension nd,
c             defined at same grid levels as t, q, qs and p.
c
c     fu:   array of forcing of zonal velocity (m/s^2) of dimension nd,
c             defined at same grid levels as t.
c
c     fv:   same as fu, but for forcing of meridional velocity.
c
c     ftra: array of forcing of tracer content, in tracer mixing ratio per
c             second, defined at same levels as t. dimensioned (nd,ntra).
c
c     precip: scalar convective precipitation rate (mm/s).
c
c     wd:    a convective downdraft velocity scale. for use in surface
c             flux parameterizations. see convect.ps file for details.
c
c     tprime: a convective downdraft temperature perturbation scale (k).
c              for use in surface flux parameterizations. see convect.ps
c              file for details.
c
c     qprime: a convective downdraft specific humidity
c              perturbation scale (gm/gm).
c              for use in surface flux parameterizations. see convect.ps
c              file for details.
c
c     cbmf:   the cloud base mass flux ((kg/m**2)/s). this scalar value must
c              be stored by the calling program and returned to convect at
c              its next call. that is, the value of cbmf must be "remembered"
c              by the calling program between calls to convect.
c
c------------------------------------------------------------------------------
c
c    ***  the parameter na should in general be greater than   ***
c    ***                or equal to  nd + 1                    ***
c
c------------------------------------------------------------------------------
c modifications for regcm:
c   1. units for precipitation were change from mm/day to mm/s
c   2. the thermodynamic constants were made consistent with those
c      of regcm.
c   3. dependence of latent heat of vaporization on tempertature
c      removed because it is neglected in regcm. this is done by
c      setting cpv equal to cl and setting lv0 to the regcm value.
c   4. added cloud base (icb) and cloud top (inb) to the output to
c      compute the cloud fraction and cloud liquid water content.
c   5. each variable is now explicitly declared.  that is, the
c      "implicit none" option was added.
c   6. the value minorig is increased because the thickness of the
c      lowest layer(s) is(are) too small. if the thickness is too
c      small for the given timestep, the mass of the layer is likely
c      to be evacuated.
c   7. a maximum value to the cloud base mass flux has been added.
c

c ****  regcm add-ins to declare every variable explicitly **** c
      implicit none

#     include "include/cupeman.cb"

      integer na, nd, ntra, ipbl, i, j, k, nl, iflag, jc, jn
     &      , ihmin, nk, icb, inb, inb1, jtt
      real(kind=8)  rowl, cpvmcl, eps, epsi, ginv, delti, delt, rdcp
     &      , precip, wd, tprime, qprime, asum, thbar, ahm, rm
     &      , um, vm, dphinv, a2, x, tc, alv, tnew, alvnew
     &      , qnew, ahmin, tvx, tvy, ahmax, cbmf, rh, chi
     &      , plcl, tca, elacrit, cape, capem, by, byp, defrac
     &      , frac, dbosum, tvpplcl, tvaplcl, dtpbl, dtmin
     &      , dtma, cbmfold, delt0, damps, dbo, qti, bf2, anum
     &      , denom, dei, altem, qp1, scrit, alt, asij, smin
     &      , smid, sjmax, sjmin, delp, delm, bsum, wdtrain
     &      , awat, coeff, qsm, afac, sigt, b6, c6, revap
     &      , dhdp, fac, qstm, rat, dpinv, am, cpinv
     &      , amp1, ad, fqold, ftold, fuold, fvold, ftraold
     &      , ents, epmax, cwat, stemp, uav, vav, traav, cpd
     &      , cpv, cl, rv, g, rd
c     real(kind=8) cbmfmax
c ****  end regcm add-ins to declare every variable explicitly **** c

c     parameter (na=48)           
c
      integer nent(na)
c     real(kind=8)  t(nd),q(nd),qs(nd),u(nd),v(nd),tra(nd,ntra),p(nd),ph(nd)
      real(kind=8)  t(nd),q(nd),qs(nd),u(nd),v(nd),tra(nd,1)
      real(kind=8)  p(nd),ph(nd)
c     real(kind=8)  ft(nd),fq(nd),fu(nd),fv(nd),ftra(nd,ntra)
      real(kind=8)  ft(nd),fq(nd),fu(nd),fv(nd),ftra(nd,1)   
      real(kind=8)  uent(na,na),vent(na,na),traent(na,na,ntra),tratm(na)
      real(kind=8)  up(na),vp(na),trap(na,ntra)
      real(kind=8)  m(na),mp(na),ment(na,na),qent(na,na),elij(na,na)
      real(kind=8)  sij(na,na),tvp(na),tv(na),water(na)
      real(kind=8)  qp(na),ep(na),th(na),wt(na),evap(na),clw(na)
      real(kind=8)  sigp(na),tp(na),told(na),cpn(na)
      real(kind=8)  lv(na),lvcp(na),lv0,h(na),hp(na),gz(na),hm(na)
c
c -----------------------------------------------------------------------
c
c   ***                     specify switches                         ***
c
c   ***   ipbl: set to zero to bypass dry adiabatic adjustment       ***
c   ***    any other value results in dry adiabatic adjustment       ***
c   ***     (zero value recommended for use in models with           ***
c   ***                   boundary layer schemes)                    ***
c
c   ***   minorig: lowest level from which convection may originate  ***
c   ***     (should be first model level at which t is defined       ***
c   ***      for models using bulk pbl schemes; otherwise, it should ***
c   ***      be the first model level at which t is defined above    ***
c   ***                      the surface layer)                      ***
c
        ipbl=0
crcm    minorig=2
c
c------------------------------------------------------------------------------
c
c   ***                    specify parameters                        ***
c
c   *** elcrit is the autoconversion thershold water content (gm/gm) ***
c   ***  tlcrit is critical temperature below which the auto-        ***
c   ***       conversion threshold is assumed to be zero             ***
c   ***     (the autoconversion threshold varies linearly            ***
c   ***               between 0 c and tlcrit)                        ***
c   ***   entp is the coefficient of mixing in the entrainment       ***
c   ***                       formulation                            ***
c   ***  sigd is the fractional area covered by unsaturated dndraft  ***
c   ***  sigs is the fraction of precipitation falling outside       ***
c   ***                        of cloud                              ***
c   ***        omtrain is the assumed fall speed (p/s) of rain       ***
c   ***     omtsnow is the assumed fall speed (p/s) of snow          ***
c   ***  coeffr is a coefficient governing the rate of evaporation   ***
c   ***                          of rain                             ***
c   ***  coeffs is a coefficient governing the rate of evaporation   ***
c   ***                          of snow                             ***
c   ***     cu is the coefficient governing convective momentum      ***
c   ***                         transport                            ***
c   ***    dtmax is the maximum negative temperature perturbation    ***
c   ***        a lifted parcel is allowed to have below its lfc      ***
c   ***    alphae and damp are parameters that control the rate of   ***
c   ***                 approach to quasi-equilibrium                ***
c   ***   (their standard values are  0.20 and 0.1, respectively)    ***
c   ***                   (damp must be less than 1)                 ***
c
crcm addin
c the following parameters in now in param.f
c       elcrit=.0011
c       tlcrit=-55.0
c       entp=1.5
c       sigd=0.05
c       sigs=0.12
c       omtrain=50.0
c       omtsnow=5.5 
c       coeffr=1.0
c       coeffs=0.8
c       cu=0.7
c       betae=10.0
c       dtmax=0.9
c       alphae=0.2
c       damp=0.1
crcm end addin
c
c   ***        assign values of thermodynamic constants,        ***
c   ***            gravity, and liquid water density.           ***
c   ***             these should be consistent with             ***
c   ***              those used in calling program              ***
c   ***     note: these are also specified in subroutine tlift  ***
c
crcm  cpd=1005.7
      cpd=1004.
      cpv=1870.0
crcm  cpv=cpd
      cl=2500.0
      rv=461.5
      rd=287.04
      lv0=2.501e6
crcm  g=9.8
      g=9.805
      rowl=1000.0
c
      cpvmcl=cl-cpv 
      eps=rd/rv
      epsi=1./eps
      ginv=1.0/g
      delti=1.0/delt
c
c           ***  initialize output arrays and parameters  ***
c
      do 5 i=1,nd
        ft(i)=0.0
        fq(i)=0.0
        fu(i)=0.0
        fv(i)=0.0
        do 4 j=1,ntra
          ftra(i,j)=0.0
    4   continue
    5 continue
      do 7 i=1,nl+1
        rdcp=(rd*(1.-q(i))+q(i)*rv)/(cpd*(1.-q(i))+q(i)*cpv)
        th(i)=t(i)*(1000.0/p(i))**rdcp
    7 continue
      precip=0.0
      wd=0.0
      tprime=0.0
      qprime=0.0
      iflag=0
c
      if (ipbl.ne.0) then
c
c     ***            perform dry adiabatic adjustment            ***
c
        jc=0
        do 30 i=nl-1,1,-1
          jn=0
          asum=th(i)*(1.+q(i)*epsi-q(i))
          do 10 j=i+1,nl
            asum=asum+th(j)*(1.+q(j)*epsi-q(j))
            thbar=asum/float(j+1-i)
            if ((th(j)*(1.+q(j)*epsi-q(j))).lt.thbar)jn=j
   10     continue
          if (i.eq.1)jn=max(jn,2)
          if (jn.eq.0)goto 30
   12     continue
          ahm=0.0
          rm=0.0
          um=0.0
          vm=0.0
          do k=1,ntra
            tratm(k)=0.0
          end do
          do 15 j=i,jn
            ahm=ahm+(cpd*(1.-q(j))+q(j)*cpv)*t(j)*(ph(j)-ph(j+1))
            rm=rm+q(j)*(ph(j)-ph(j+1))
            um=um+u(j)*(ph(j)-ph(j+1))
            vm=vm+v(j)*(ph(j)-ph(j+1))
            do k=1,ntra
              tratm(k)=tratm(k)+tra(j,k)*(ph(j)-ph(j+1))
            end do
   15     continue
          dphinv=1./(ph(i)-ph(jn+1))
          rm=rm*dphinv
          um=um*dphinv
          vm=vm*dphinv
          do k=1,ntra
            tratm(k)=tratm(k)*dphinv
          end do
          a2=0.0
          do 20 j=i,jn
            q(j)=rm
            u(j)=um
            v(j)=vm
            do k=1,ntra
              tra(j,k)=tratm(k)
            end do
            rdcp=(rd*(1.-q(j))+q(j)*rv)/(cpd*(1.-q(j))+q(j)*cpv)  
            x=(0.001*p(j))**rdcp
            told(j)=t(j)
            t(j)=x
            a2=a2+(cpd*(1.-q(j))+q(j)*cpv)*x*(ph(j)-ph(j+1))
   20     continue
          do 25 j=i,jn
            th(j)=ahm/a2
            t(j)=t(j)*th(j)
            tc=told(j)-273.15
            alv=lv0-cpvmcl*tc
            qs(j)=qs(j)+qs(j)*(1.+qs(j)*(epsi-1.))*alv*(t(j)-
     1      told(j))/(rv*told(j)*told(j))
   25     continue
          if (((th(jn+1)*(1.+q(jn+1)*epsi-q(jn+1))).lt.
     1       (th(jn)*(1.+q(jn)*epsi-q(jn))))) then
            jn=jn+1
            goto 12
          end if
          if (i.eq.1)jc=jn 
   30   continue
c
c   ***   remove any supersaturation that results from adjustment ***
c
        if (jc.gt.1) then
          do 38 j=1,jc
            if (qs(j).lt.q(j)) then 
              alv=lv0-cpvmcl*(t(j)-273.15)  
              tnew=t(j)+alv*(q(j)-qs(j))/(cpd*(1.-q(j))+
     1         cl*q(j)+qs(j)*(cpv-cl+alv*alv/(rv*t(j)*t(j))))
              alvnew=lv0-cpvmcl*(tnew-273.15)
              qnew=(alv*q(j)-(tnew-t(j))*(cpd*(1.-q(j))+cl*q(j)))/alvnew
crcm          precip=precip+24.*3600.*1.0e5*(ph(j)-ph(j+1))*  ! mm/d
              precip=precip+1.0e5*(ph(j)-ph(j+1))*            ! mm/s
     1         (q(j)-qnew)/(g*delt*rowl)
              t(j)=tnew
              q(j)=qnew
              qs(j)=qnew
            end if     
   38     continue  
        end if
c
      end if
c
c  *** calculate arrays of geopotential, heat capacity and static energy
c  
      gz(1)=0.0
      cpn(1)=cpd*(1.-q(1))+q(1)*cpv
      h(1)=t(1)*cpn(1)
      lv(1)=lv0-cpvmcl*(t(1)-273.15)
      hm(1)=lv(1)*q(1)
      tv(1)=t(1)*(1.+q(1)*epsi-q(1))
      ahmin=1.0e12
      ihmin=nl
      do 40 i=2,nl+1
        tvx=t(i)*(1.+q(i)*epsi-q(i))
        tvy=t(i-1)*(1.+q(i-1)*epsi-q(i-1))
        gz(i)=gz(i-1)+0.5*rd*(tvx+tvy)*(p(i-1)-p(i))/ph(i)
        cpn(i)=cpd*(1.-q(i))+cpv*q(i)
        h(i)=t(i)*cpn(i)+gz(i)
        lv(i)=lv0-cpvmcl*(t(i)-273.15)
        hm(i)=(cpd*(1.-q(i))+cl*q(i))*(t(i)-t(1))+lv(i)*q(i)+gz(i)
        tv(i)=t(i)*(1.+q(i)*epsi-q(i))
c
c  ***  find level of minimum moist static energy    ***
c
        if (i.ge.minorig.and.hm(i).lt.ahmin.and.hm(i).lt.hm(i-1)) then
          ahmin=hm(i)
          ihmin=i
        end if
   40 continue
      ihmin=min(ihmin, nl-1)
c
c  ***     find that model level below the level of minimum moist       ***
c  ***  static energy that has the maximum value of moist static energy ***
c
      ahmax=0.0
      do 42 i=minorig,ihmin
        if (hm(i).gt.ahmax) then
          nk=i
          ahmax=hm(i)
        end if
   42 continue
c
c  ***  check whether parcel level temperature and specific humidity   ***
c  ***                          are reasonable                         ***
c  ***      skip convection if hm increases monotonically upward       ***
c
      if (t(nk).lt.250.0.or.q(nk).le.0.0.or.ihmin.eq.(nl-1)) then
        iflag=0
        cbmf=0.0
        return
      end if
c
c   ***  calculate lifted condensation level of air at parcel origin level ***
c   ***       (within 0.2% of formula of bolton, mon. wea. rev.,1980)      ***
c
      rh=q(nk)/qs(nk)
      chi=t(nk)/(1669.0-122.0*rh-t(nk))
      plcl=p(nk)*(rh**chi)
      if (plcl.lt.200.0.or.plcl.ge.2000.0) then
        iflag=2
        cbmf=0.0
        return
      end if
c
c   ***  calculate first level above lcl (=icb)  ***
c
      icb=nl-1
      do 50 i=nk+1,nl
        if (p(i).lt.plcl) then
          icb=min(icb,i)
        end if
   50 continue
      if (icb.ge.(nl-1)) then
        iflag=3
        cbmf=0.0
        return
      end if
c
c   *** find temperature up through icb and test for instability           ***
c
c   *** subroutine tlift calculates part of the lifted parcel virtual      ***
c   ***  temperature, the actual temperature and the adiabatic             ***
c   ***                   liquid water content                             ***
c
      call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,1)
      do 54 i=nk,icb
        tvp(i)=tvp(i)-tp(i)*q(nk)
   54 continue
c
c   ***  if there was no convection at last time step and parcel    ***
c   ***       is stable at icb then skip rest of calculation        ***
c
      if (cbmf.eq.0.0.and.tvp(icb).le.(tv(icb)-dtmax)) then
        iflag=0
        return
      end if
c
c   ***  if this point is reached, moist convective adjustment is necessary ***
c
      if (iflag.ne.4)iflag=1
c
c   ***  find the rest of the lifted parcel temperatures          ***
c
      call tlift(p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,2)
c
c   ***  set the precipitation efficiencies and the fraction of   ***
c   ***          precipitation falling outside of cloud           ***
c   ***      these may be functions of tp(i), p(i) and clw(i)     ***
c                 
      do 57 i=1,nk
        ep(i)=0.0
        sigp(i)=sigs
   57 continue
      do 60 i=nk+1,nl
        tca=tp(i)-273.15
        if (tca.ge.0.0) then
          elacrit=elcrit
        else
          elacrit=elcrit*(1.0-tca/tlcrit)
        end if
        elacrit=max(elacrit,0.0d0)
        epmax=0.999
        ep(i)=epmax*(1.0-elacrit/max(clw(i),1.0d-8))
        ep(i)=max(ep(i),0.0d0)
        ep(i)=min(ep(i),epmax)
        sigp(i)=sigs
   60 continue
c
c   ***       calculate virtual temperature and lifted parcel     ***
c   ***                    virtual temperature                    ***
c
      do 64 i=icb+1,nl
        tvp(i)=tvp(i)-tp(i)*q(nk)
   64 continue
      tvp(nl+1)=tvp(nl)-(gz(nl+1)-gz(nl))/cpd
c
c   ***        now initialize various arrays used in the computations       ***
c
      do 70 i=1,nl+1
        hp(i)=h(i)
        nent(i)=0
        water(i)=0.0
        evap(i)=0.0
        wt(i)=omtsnow
        mp(i)=0.0
        m(i)=0.0
        lvcp(i)=lv(i)/cpn(i)
        do 70 j=1,nl+1
          qent(i,j)=q(j)
          elij(i,j)=0.0
          ment(i,j)=0.0
          sij(i,j)=0.0
          uent(i,j)=u(j)
          vent(i,j)=v(j)
          do 70 k=1,ntra
            traent(i,j,k)=tra(j,k)
   70 continue
      qp(1)=q(1)
      up(1)=u(1)
      vp(1)=v(1)
      do 71 i=1,ntra
        trap(1,i)=tra(1,i)
   71 continue
      do 72 i=2,nl+1
        qp(i)=q(i-1)
        up(i)=u(i-1)
        vp(i)=v(i-1)
        do 72 j=1,ntra
          trap(i,j)=tra(i-1,j)
   72 continue
c
c  ***  find the first model level (inb1) above the parcel's      ***
c  ***          highest level of neutral buoyancy                 ***
c  ***     and the highest level of positive cape (inb)           ***
c
      cape=0.0
      capem=0.0
      inb=icb+1
      inb1=inb
      byp=0.0
      do 82 i=icb+1,nl-1
        by=(tvp(i)-tv(i))*(ph(i)-ph(i+1))/p(i)
        cape=cape+by
        if (by.ge.0.0)inb1=i+1
        if (cape.gt.0.0) then
          inb=i+1
          byp=(tvp(i+1)-tv(i+1))*(ph(i+1)-ph(i+2))/p(i+1)
          capem=cape
        end if
   82 continue
      inb=max(inb,inb1)
      cape=capem+byp
      defrac=capem-cape
      defrac=max(defrac,0.001d0)
      frac=-cape/defrac
      frac=min(frac,1.0d0)
      frac=max(frac,0.0d0)
c
c   ***   calculate liquid water static energy of lifted parcel   ***
c
      do 95 i=icb,inb
        hp(i)=h(nk)+(lv(i)+(cpd-cpv)*t(i))*ep(i)*clw(i)
   95 continue                  
c
c   ***  calculate cloud base mass flux and rates of mixing, m(i),  ***
c   ***                   at each model level                       ***
c
      dbosum=0.0
c   
c   ***     interpolate difference between lifted parcel and      ***
c   ***  environmental temperatures to lifted condensation level  ***
c
      tvpplcl=tvp(icb-1)-rd*tvp(icb-1)*(p(icb-1)-plcl)/
     1    (cpn(icb-1)*p(icb-1))
      tvaplcl=tv(icb)+(tvp(icb)-tvp(icb+1))*(plcl-p(icb))/
     1    (p(icb)-p(icb+1))
      dtpbl=0.0
      do 96 i=nk,icb-1
        dtpbl=dtpbl+(tvp(i)-tv(i))*(ph(i)-ph(i+1))
   96 continue
      dtpbl=dtpbl/(ph(nk)-ph(icb))
      dtmin=tvpplcl-tvaplcl+dtmax+dtpbl
      dtma=dtmin
c
c   ***  adjust cloud base mass flux   ***
c
      cbmfold=cbmf
      delt0=300.0
      damps=damp*delt/delt0 
      cbmf=(1.-damps)*cbmf+0.1*alphae*dtma 
      cbmf=max(cbmf,0.0d0)
crcm addin
c     cbmfmax = delti*ginv*(ph(minorig)-ph(minorig+1))/(0.01*2.)
c     cbmfmax = delti*ginv*(ph(1)-ph(2))/(0.01*2.)
c     cbmf=min(cbmf,0.75*cbmfmax)
crcm end addin
c
c   *** if cloud base mass flux is zero, skip rest of calculation  ***
c
      if (cbmf.eq.0.0.and.cbmfold.eq.0.0) then
        return
      end if
c
c   ***   calculate rates of mixing,  m(i)   ***
c
      m(icb)=0.0
      do 103 i=icb+1,inb
        k=min(i,inb1)
        dbo=abs(tv(k)-tvp(k))+entp*0.02*(ph(k)-ph(k+1))
        dbosum=dbosum+dbo
        m(i)=cbmf*dbo
  103 continue
      do 110 i=icb+1,inb
        m(i)=m(i)/dbosum  
  110 continue     
c
c   ***  calculate entrained air mass flux (ment), total water mixing  ***
c   ***     ratio (qent), total condensed water (elij), and mixing     ***
c   ***                        fraction (sij)                          ***
c
      do 170 i=icb+1,inb
        qti=q(nk)-ep(i)*clw(i)
        do 160 j=icb,inb
          bf2=1.+lv(j)*lv(j)*qs(j)/(rv*t(j)*t(j)*cpd)
          anum=h(j)-hp(i)+(cpv-cpd)*t(j)*(qti-q(j))
          denom=h(i)-hp(i)+(cpd-cpv)*(q(i)-qti)*t(j)
          dei=denom
          if (abs(dei).lt.0.01)dei=0.01
          sij(i,j)=anum/dei
          sij(i,i)=1.0
          altem=sij(i,j)*q(i)+(1.-sij(i,j))*qti-qs(j)
          altem=altem/bf2
          cwat=clw(j)*(1.-ep(j))
          stemp=sij(i,j)
          if ((stemp.lt.0.0.or.stemp.gt.1.0.or.altem.gt.cwat).and.
     &        j.gt.i) then
            anum=anum-lv(j)*(qti-qs(j)-cwat*bf2)
            denom=denom+lv(j)*(q(i)-qti)
            if (abs(denom).lt.0.01)denom=0.01
            sij(i,j)=anum/denom
            altem=sij(i,j)*q(i)+(1.-sij(i,j))*qti-qs(j)
            altem=altem-(bf2-1.)*cwat
          end if
          if (sij(i,j).gt.0.0.and.sij(i,j).lt.0.9) then
            qent(i,j)=sij(i,j)*q(i)+(1.-sij(i,j))*qti
            uent(i,j)=sij(i,j)*u(i)+(1.-sij(i,j))*u(nk)
            vent(i,j)=sij(i,j)*v(i)+(1.-sij(i,j))*v(nk)
            do k=1,ntra
              traent(i,j,k)=sij(i,j)*tra(i,k)+(1.-sij(i,j))*tra(nk,k)
            end do
            elij(i,j)=altem
            elij(i,j)=max(0.0d0,elij(i,j))
            ment(i,j)=m(i)/(1.-sij(i,j))
            nent(i)=nent(i)+1
          end if
          sij(i,j)=max(0.0d0,sij(i,j))
          sij(i,j)=min(1.0d0,sij(i,j))
  160   continue
c
c   ***   if no air can entrain at level i assume that updraft detrains  ***
c   ***   at that level and calculate detrained air flux and properties  ***
c
        if (nent(i).eq.0) then
          ment(i,i)=m(i)
          qent(i,i)=q(nk)-ep(i)*clw(i)
          uent(i,i)=u(nk)
          vent(i,i)=v(nk)
          do j=1,ntra
            traent(i,i,j)=tra(nk,j)
          end do
          elij(i,i)=clw(i)
          sij(i,i)=1.0
        end if 
  170 continue
      sij(inb,inb)=1.0
c
c   ***  normalize entrained air mass fluxes to represent equal  ***
c   ***              probabilities of mixing                     ***
c
      do 200 i=icb+1,inb
        if (nent(i).ne.0) then
          qp1=q(nk)-ep(i)*clw(i)
          anum=h(i)-hp(i)-lv(i)*(qp1-qs(i))
          denom=h(i)-hp(i)+lv(i)*(q(i)-qp1)
          if (abs(denom).lt.0.01)denom=0.01
          scrit=anum/denom
          alt=qp1-qs(i)+scrit*(q(i)-qp1)
          if (alt.lt.0.0)scrit=1.0
          scrit=max(scrit,0.0d0)
          asij=0.0
          smin=1.0
          do 175 j=icb,inb
            if (sij(i,j).gt.0.0.and.sij(i,j).lt.0.9) then
              if (j.gt.i) then
                smid=min(sij(i,j),scrit)
                sjmax=smid
                sjmin=smid
                if (smid.lt.smin.and.sij(i,j+1).lt.smid) then
                  smin=smid
                  sjmax=min(sij(i,j+1),sij(i,j),scrit)
                  sjmin=max(sij(i,j-1),sij(i,j))
                  sjmin=min(sjmin,scrit)
                end if
              else
                sjmax=max(sij(i,j+1),scrit)
                smid=max(sij(i,j),scrit)
                sjmin=0.0
                if (j.gt.1)sjmin=sij(i,j-1)
                sjmin=max(sjmin,scrit)
              end if
              delp=abs(sjmax-smid)
              delm=abs(sjmin-smid)
              asij=asij+(delp+delm)*(ph(j)-ph(j+1))
              ment(i,j)=ment(i,j)*(delp+delm)*(ph(j)-ph(j+1))
            end if
  175     continue
          asij=max(1.0d-21,asij)
          asij=1.0/asij
          do 180 j=icb,inb
            ment(i,j)=ment(i,j)*asij
  180     continue
          bsum=0.0
          do 190 j=icb,inb
            bsum=bsum+ment(i,j)
  190     continue
          if (bsum.lt.1.0e-18) then
            nent(i)=0
            ment(i,i)=m(i)
            qent(i,i)=q(nk)-ep(i)*clw(i)
            uent(i,i)=u(nk)
            vent(i,i)=v(nk)
            do j=1,ntra
              traent(i,i,j)=tra(nk,j)
            end do
            elij(i,i)=clw(i)
            sij(i,i)=1.0
          end if
        end if
  200 continue
c
c   ***  check whether ep(inb)=0, if so, skip precipitating    ***
c   ***             downdraft calculation                      ***
c
      if (ep(inb).lt.0.0001)goto 405
c
c   ***  integrate liquid water equation to find condensed water   ***
c   ***                and condensed water flux                    ***
c
      jtt=2
c
c    ***                    begin downdraft loop                    ***
c
      do 400 i=inb,1,-1
c
c    ***              calculate detrained precipitation             ***
c
        wdtrain=g*ep(i)*m(i)*clw(i)
        if (i.gt.1) then
          do 320 j=1,i-1
            awat=elij(j,i)-(1.-ep(i))*clw(i)
            awat=max(0.0d0,awat)
  320       wdtrain=wdtrain+g*awat*ment(j,i)
        end if
c
c    ***    find rain water and evaporation using provisional   ***
c    ***              estimates of qp(i)and qp(i-1)             ***
c     
c
c  ***  value of terminal velocity and coefficient of evaporation for snow   ***
c 
        coeff=coeffs
        wt(i)=omtsnow
c      
c  ***  value of terminal velocity and coefficient of evaporation for rain   ***
c
        if (t(i).gt.273.0) then
          coeff=coeffr
          wt(i)=omtrain
        end if
        qsm=0.5*(q(i)+qp(i+1))
        afac=coeff*ph(i)*(qs(i)-qsm)/(1.0e4+2.0e3*ph(i)*qs(i))
        afac=max(afac,0.0d0)
        sigt=sigp(i)
        sigt=max(0.0d0,sigt)
        sigt=min(1.0d0,sigt)
        b6=100.*(ph(i)-ph(i+1))*sigt*afac/wt(i)
        c6=(water(i+1)*wt(i+1)+wdtrain/sigd)/wt(i)
        revap=0.5*(-b6+sqrt(b6*b6+4.*c6))
        evap(i)=sigt*afac*revap
        water(i)=revap*revap
c
c    ***  calculate precipitating downdraft mass flux under     ***
c    ***              hydrostatic approximation                 ***
c   
        if (i.eq.1)goto 360
        dhdp=(h(i)-h(i-1))/(p(i-1)-p(i))
        dhdp=max(dhdp,10.0d0)
        mp(i)=100.*ginv*lv(i)*sigd*evap(i)/dhdp
        mp(i)=max(mp(i),0.0d0)
c
c   ***   add small amount of inertia to downdraft              ***
c
        fac=20.0/(ph(i-1)-ph(i))
        mp(i)=(fac*mp(i+1)+mp(i))/(1.+fac)
c   
c    ***      force mp to decrease linearly to zero                 ***
c    ***      between about 950 mb and the surface                  ***
c
        if (p(i).gt.(0.949*p(1))) then
          jtt=max(jtt,i)
          mp(i)=mp(jtt)*(p(1)-p(i))/(p(1)-p(jtt))
        end if              
  360   continue
c
c    ***       find mixing ratio of precipitating downdraft     ***
c
        if (i.eq.inb)goto 400
        if (i.eq.1) then
          qstm=qs(1)
        else
          qstm=qs(i-1)
        end if
        if (mp(i).gt.mp(i+1)) then
          rat=mp(i+1)/mp(i)
          qp(i)=qp(i+1)*rat+q(i)*(1.0-rat)+100.*ginv*
     1       sigd*(ph(i)-ph(i+1))*(evap(i)/mp(i))
          up(i)=up(i+1)*rat+u(i)*(1.-rat)
          vp(i)=vp(i+1)*rat+v(i)*(1.-rat)
          do j=1,ntra
            trap(i,j)=trap(i+1,j)*rat+trap(i,j)*(1.-rat)
          end do
        else
          if (mp(i+1).gt.0.0) then
            qp(i)=(gz(i+1)-gz(i)+qp(i+1)*(lv(i+1)+t(i+1)*(
     1        cl-cpd))+cpd*(t(i+1)-t(i)))/(lv(i)+t(i)*(cl-cpd))
            up(i)=up(i+1)
            vp(i)=vp(i+1)
            do j=1,ntra
              trap(i,j)=trap(i+1,j)
            end do
          end if
        end if
        qp(i)=min(qp(i),qstm)
        qp(i)=max(qp(i),0.0d0)
  400 continue
c
c   ***  calculate surface precipitation in mm/s     ***
c
crcm  precip=precip+wt(1)*sigd*water(1)*3600.*24000./(rowl*g)  ! mm/d
      precip=precip+wt(1)*sigd*water(1)*1000./(rowl*g)  ! mm/s
c
  405 continue
c
c   ***  calculate downdraft velocity scale and surface temperature and  ***
c   ***                    water vapor fluctuations                      ***
c
      wd=betae*abs(mp(icb))*0.01*rd*t(icb)/(sigd*p(icb))
      qprime=0.5*(qp(1)-q(1))
      tprime=lv0*qprime/cpd
c
c   ***  calculate tendencies of lowest level potential temperature  ***
c   ***                      and mixing ratio                        ***
c
      dpinv=0.01/(ph(1)-ph(2))
      am=0.0
      if (nk.eq.1) then
        do 410 k=2,inb
  410     am=am+m(k)
      end if
      if ((2.*g*dpinv*am).ge.delti) iflag=4
      ft(1)=ft(1)+g*dpinv*am*(t(2)-t(1)+(gz(2)-gz(1))/cpn(1))
      ft(1)=ft(1)-lvcp(1)*sigd*evap(1)
      ft(1)=ft(1)+sigd*wt(2)*(cl-cpd)*water(2)*(t(2)-t(1))*dpinv/cpn(1)
      fq(1)=fq(1)+g*mp(2)*(qp(2)-q(1))*dpinv+sigd*evap(1)
      fq(1)=fq(1)+g*am*(q(2)-q(1))*dpinv
      fu(1)=fu(1)+g*dpinv*(mp(2)*(up(2)-u(1))+am*(u(2)-u(1)))
      fv(1)=fv(1)+g*dpinv*(mp(2)*(vp(2)-v(1))+am*(v(2)-v(1)))
      do j=1,ntra
        ftra(1,j)=ftra(1,j)+g*dpinv*(mp(2)*(trap(2,j)-tra(1,j))+
     1    am*(tra(2,j)-tra(1,j)))
      end do
      do 415 j=2,inb
        fq(1)=fq(1)+g*dpinv*ment(j,1)*(qent(j,1)-q(1))
        fu(1)=fu(1)+g*dpinv*ment(j,1)*(uent(j,1)-u(1))
        fv(1)=fv(1)+g*dpinv*ment(j,1)*(vent(j,1)-v(1))
        do k=1,ntra
          ftra(1,k)=ftra(1,k)+g*dpinv*ment(j,1)*(traent(j,1,k)-
     1     tra(1,k))
        end do
  415 continue
c
c   ***  calculate tendencies of potential temperature and mixing ratio  ***
c   ***               at levels above the lowest level                   ***
c
c   ***  first find the net saturated updraft and downdraft mass fluxes  ***
c   ***                      through each level                          ***
c
      do 500 i=2,inb
        dpinv=0.01/(ph(i)-ph(i+1))
        cpinv=1.0/cpn(i)
        amp1=0.0
        ad=0.0
        if (i.ge.nk) then
          do 440 k=i+1,inb+1
  440       amp1=amp1+m(k)
        end if
        do 450 k=1,i
        do 450 j=i+1,inb+1
          amp1=amp1+ment(k,j)
  450   continue
        if ((2.*g*dpinv*amp1).ge.delti) iflag=4
        do 470 k=1,i-1
        do 470 j=i,inb
          ad=ad+ment(j,k)
  470   continue
        ft(i)=ft(i)+g*dpinv*(amp1*(t(i+1)-t(i)+(gz(i+1)-gz(i))*
     1   cpinv)-ad*(t(i)-t(i-1)+(gz(i)-gz(i-1))*cpinv))
     2   -sigd*lvcp(i)*evap(i)
        ft(i)=ft(i)+g*dpinv*ment(i,i)*(hp(i)-h(i)+
     1    t(i)*(cpv-cpd)*(q(i)-qent(i,i)))*cpinv
        ft(i)=ft(i)+sigd*wt(i+1)*(cl-cpd)*water(i+1)*
     1    (t(i+1)-t(i))*dpinv*cpinv
        fq(i)=fq(i)+g*dpinv*(amp1*(q(i+1)-q(i))-ad*(q(i)-q(i-1)))
        fu(i)=fu(i)+g*dpinv*(amp1*(u(i+1)-u(i))-ad*(u(i)-u(i-1)))
        fv(i)=fv(i)+g*dpinv*(amp1*(v(i+1)-v(i))-ad*(v(i)-v(i-1)))
        do k=1,ntra
          ftra(i,k)=ftra(i,k)+g*dpinv*(amp1*(tra(i+1,k)-
     1               tra(i,k))-ad*(tra(i,k)-tra(i-1,k)))
        end do
        do 480 k=1,i-1
          awat=elij(k,i)-(1.-ep(i))*clw(i)
          awat=max(awat,0.0d0)
          fq(i)=fq(i)+g*dpinv*ment(k,i)*(qent(k,i)-awat-q(i))
          fu(i)=fu(i)+g*dpinv*ment(k,i)*(uent(k,i)-u(i))
          fv(i)=fv(i)+g*dpinv*ment(k,i)*(vent(k,i)-v(i))
          do j=1,ntra
            ftra(i,j)=ftra(i,j)+g*dpinv*ment(k,i)*(traent(k,i,j)-
     1     tra(i,j))
          end do
  480   continue
        do 490 k=i,inb
          fq(i)=fq(i)+g*dpinv*ment(k,i)*(qent(k,i)-q(i))
          fu(i)=fu(i)+g*dpinv*ment(k,i)*(uent(k,i)-u(i))
          fv(i)=fv(i)+g*dpinv*ment(k,i)*(vent(k,i)-v(i))
          do j=1,ntra
            ftra(i,j)=ftra(i,j)+g*dpinv*ment(k,i)*(traent(k,i,j)-
     1     tra(i,j))
          end do
  490   continue
        fq(i)=fq(i)+sigd*evap(i)+g*(mp(i+1)*
     1    (qp(i+1)-q(i))-mp(i)*(qp(i)-q(i-1)))*dpinv
        fu(i)=fu(i)+g*(mp(i+1)*(up(i+1)-u(i))-mp(i)*
     1    (up(i)-u(i-1)))*dpinv
        fv(i)=fv(i)+g*(mp(i+1)*(vp(i+1)-v(i))-mp(i)*
     1    (vp(i)-v(i-1)))*dpinv
        do j=1,ntra
          ftra(i,j)=ftra(i,j)+g*dpinv*(mp(i+1)*(trap(i+1,j)-tra(i,j))-
     1    mp(i)*(trap(i,j)-trap(i-1,j)))
        end do
  500 continue
c
c   *** adjust tendencies at top of convection layer to reflect  ***
c   ***       actual position of the level zero cape             ***
c
      fqold=fq(inb)
      fq(inb)=fq(inb)*(1.-frac)
      fq(inb-1)=fq(inb-1)+frac*fqold*((ph(inb)-ph(inb+1))/
     1   (ph(inb-1)-ph(inb)))*lv(inb)/lv(inb-1)
      ftold=ft(inb)
      ft(inb)=ft(inb)*(1.-frac)
      ft(inb-1)=ft(inb-1)+frac*ftold*((ph(inb)-ph(inb+1))/
     1   (ph(inb-1)-ph(inb)))*cpn(inb)/cpn(inb-1)
      fuold=fu(inb)
      fu(inb)=fu(inb)*(1.-frac)
      fu(inb-1)=fu(inb-1)+frac*fuold*((ph(inb)-ph(inb+1))/
     1   (ph(inb-1)-ph(inb)))
      fvold=fv(inb)
      fv(inb)=fv(inb)*(1.-frac)
      fv(inb-1)=fv(inb-1)+frac*fvold*((ph(inb)-ph(inb+1))/
     1   (ph(inb-1)-ph(inb)))
      do k=1,ntra
        ftraold=ftra(inb,k)
        ftra(inb,k)=ftra(inb,k)*(1.-frac)
        ftra(inb-1,k)=ftra(inb-1,k)+frac*ftraold*(ph(inb)-ph(inb+1))/
     1    (ph(inb-1)-ph(inb))
      end do
c
c   ***   very slightly adjust tendencies to force exact   ***
c   ***     enthalpy, momentum and tracer conservation     ***
c
      ents=0.0
      uav=0.0
      vav=0.0
      do 680 i=1,inb
        ents=ents+(cpn(i)*ft(i)+lv(i)*fq(i))*(ph(i)-ph(i+1))
        uav=uav+fu(i)*(ph(i)-ph(i+1))
        vav=vav+fv(i)*(ph(i)-ph(i+1))
  680 continue
      ents=ents/(ph(1)-ph(inb+1))
      uav=uav/(ph(1)-ph(inb+1))
      vav=vav/(ph(1)-ph(inb+1))
      do 640 i=1,inb
        ft(i)=ft(i)-ents/cpn(i)
        fu(i)=(1.-cu)*(fu(i)-uav)
        fv(i)=(1.-cu)*(fv(i)-vav)
  640 continue
      do 700 k=1,ntra
        traav=0.0
        do 690 i=1,inb
          traav=traav+ftra(i,k)*(ph(i)-ph(i+1))
  690   continue
        traav=traav/(ph(1)-ph(inb+1))
        do 695 i=1,inb
          ftra(i,k)=ftra(i,k)-traav
  695   continue
  700 continue
c
c   ***           return           ***
c
      return
      end subroutine cupeman
c
c ---------------------------------------------------------------------------
c
      subroutine tlift(p,t,q,qs,gz,icb,nk,tvp,tpk,clw,nd,nl,kk)
c ****  regcm add-ins to declare every variable explicitly **** c
      implicit none
      integer i, j, nst, nsb, icb, nk, nd, nl, kk
      real(kind=8)  cpd, cpv, cl, rv, rd, cpvmcl, eps, epsi, ah0, cpp
     &       , cpinv, tg, qg, alv, s, ahg, tc, denom, es, rg
c ****  end of regcm add-ins to declare every variable explicitly **** c
      real(kind=8)  gz(nd),tpk(nd),clw(nd),p(nd)
      real(kind=8)  t(nd),q(nd),qs(nd),tvp(nd),lv0
c
c   ***   assign values of thermodynamic constants     ***
c
crcm  cpd=1005.7
      cpd=1004.
      cpv=1870.0
crcm  cpv=cpd
      cl=2500.0
      rv=461.5
      rd=287.04
      lv0=2.501e6
c
      cpvmcl=cl-cpv
      eps=rd/rv
      epsi=1./eps
c
c   ***  calculate certain parcel quantities, including static energy   ***
c
      ah0=(cpd*(1.-q(nk))+cl*q(nk))*t(nk)+q(nk)*(lv0-cpvmcl*(
     1   t(nk)-273.15))+gz(nk)
      cpp=cpd*(1.-q(nk))+q(nk)*cpv
      cpinv=1./cpp
c
      if (kk.eq.1) then
c
c   ***   calculate lifted parcel quantities below cloud base   ***
c
        do 50 i=1,icb-1
          clw(i)=0.0
   50   continue
        do 100 i=nk,icb-1
          tpk(i)=t(nk)-(gz(i)-gz(nk))*cpinv
          tvp(i)=tpk(i)*(1.+q(nk)*epsi)
  100   continue
      end if
c
c    ***  find lifted parcel quantities above cloud base    ***
c
      nst=icb
      nsb=icb
      if (kk.eq.2) then  
        nst=nl
        nsb=icb+1
      end if
      do 300 i=nsb,nst
        tg=t(i)
        qg=qs(i)
        alv=lv0-cpvmcl*(t(i)-273.15)
        do 200 j=1,2
          s=cpd+alv*alv*qg/(rv*t(i)*t(i))
          s=1./s
          ahg=cpd*tg+(cl-cpd)*q(nk)*t(i)+alv*qg+gz(i)
          tg=tg+s*(ah0-ahg)
          tg=max(tg,35.0d0)
          tc=tg-273.15
          denom=243.5+tc
          if (tc.ge.0.0) then  
            es=6.112*exp(17.67*tc/denom)
          else  
            es=exp(23.33086-6111.72784/tg+0.15215*log(tg))
          end if  
          qg=eps*es/(p(i)-es*(1.-eps))
  200   continue
        alv=lv0-cpvmcl*(t(i)-273.15)
        tpk(i)=(ah0-(cl-cpd)*q(nk)*t(i)-gz(i)-alv*qg)/cpd
        clw(i)=q(nk)-qg
        clw(i)=max(0.0d0,clw(i))
        rg=qg/(1.-q(nk))
        tvp(i)=tpk(i)*(1.+rg*epsi)
  300 continue
      return
      end subroutine tlift
