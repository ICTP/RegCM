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
 
      subroutine bmpara(tten,qten,j)
!
! modified by jack kain of penn state to replace the look-up tabl
!  by calculations.
!
!*****************************************************************
!                                                                *
!  convective adjustment for deep or shallow convection          *
!                                                                *
!  references:                                                   *
!                                                                *
!  betts, a.k., 1986:  a new convective adjustment scheme.       *
!    part i: observational and theoretical basis.  quart. j. r.  *
!    met. soc., 112, 677-691.                                    *
!                                                                *
!  betts, a.k., and m.j. miller, 1986:  a new convective         *
!    adjustment scheme.  part ii: single column tests using      *
!    gate wave, bomex, atex and arctic air mass data sets.       *
!    quart. j. r. met. soc., 112, 693-709.                       *
!                                                                *
!  n.b.  part of the code is scalar.  in global models           *
!  convection occurs in less than 30/100 points.  with           *
!  simulataneous vector processing for both deep and shallow     *
!  convection, there would be a lot of redundant vector          *
!  computations.  if vector processing is 10 times faster        *
!  than scalar, one might hope that the cpu time will be about   *
!  the same for both scalar and vector code.                     *
!                                                                *
!*****************************************************************
! *** warning: this subroutine will not work if kx.lt.12;
!
      use regcm_param
      use param1
      use param2
      use param3
      use main
      use pmoist
      use rad
      use mod_bats , only : pptc , veg2d
      use bmparam
      use trachem
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: h1 = 1.E0 , h3000 = 3000.E0 ,              &
                           & h10e5 = 100000.E0 , d00 = 0.E0 ,           &
                           & d608 = 0.608E0 , dm2859 = -287.04/1004.E0 ,&
                           & elwv = 2.50E6 , row = 1.E3 ,               &
                           & epsq = 2.E-12 , a2 = 17.2693882E0 ,        &
                           & a3 = 273.16E0 , a4 = 35.86E0 ,             &
                           & t0 = 273.16E0 , t1 = 274.16E0 ,            &
                           & pq0 = 379.90516E0 , stresh = 1.10E0 ,      &
                           & stabs = 1.0E0 , stabd = 0.90E0 ,           &
                           & rhf = 0.20 , pmn = 6500.0 , epsdn = 1.05 , &
                           & epsth = 6.0 , pbm = 30000.0 ,              &
                           & pqm = 20000.0 , pone = 2500.0 ,            &
                           & pfrz = 15000.0 , pshu = 45000.0 ,          &
                           & zno = 750.0 , zsh = 3999.0
      logical , parameter :: unis = .false. , unil = .true. ,           &
                           & oct90 = .true.
      real(8) , parameter :: fss = 0.60E0 , efimn = 0.20E0 ,            &
                           & efmnt = 0.70E0 , fcc1 = 0.50 ,             &
                           & fcp = h1 - fcc1 , dspbfl = -3875.E0 ,      &
                           & dsp0fl = -5875.E0 , dsptfl = -1875.E0 ,    &
                           & fsl = 1.0E0 , dspbfs = -3875.E0 ,          &
                           & dsp0fs = -5875.E0 , dsptfs = -1875.E0 ,    &
                           & dspbsl = dspbfl*fsl , dsp0sl = dsp0fl*fsl ,&
                           & dsptsl = dsptfl*fsl , dspbss = dspbfs*fss ,&
                           & dsp0ss = dsp0fs*fss , dsptss = dsptfs*fss ,&
                           & epsntp = 0.0010E0 , efifc = 5.0E0 ,        &
                           & avgefi = (efimn+1.E0)*.5E0 ,               &
                           & dspc = -3000.E0 , epsp = 1.E-7 ,           &
                           & stefi = avgefi , slopbl = (dspbfl-dspbsl)  &
                           & /(h1-efimn) , slop0l = (dsp0fl-dsp0sl)     &
                           & /(h1-efimn) , sloptl = (dsptfl-dsptsl)     &
                           & /(h1-efimn) , slopbs = (dspbfs-dspbss)     &
                           & /(h1-efimn) , slop0s = (dsp0fs-dsp0ss)     &
                           & /(h1-efimn) , slopts = (dsptfs-dsptss)     &
                           & /(h1-efimn) , slope = (h1-efmnt)/(h1-efimn)&
                           & , a23m4l = a2*(a3-a4)*elwv , d273 = 1./t0 ,&
                           & cporng = 1./dm2859 , elocp = elwv/1004. ,  &
                           & cprlg = 1004./(row*9.81*elwv) ,            &
                           & rcp = h1/1004.
      integer , parameter :: lp1 = kx + 1 , lm1 = kx - 1
!
! Dummy arguments
!
      integer :: j
      real(8) , dimension(ix,kx) :: qten , tten
      intent (in) j
      intent (inout) qten , tten
!
! Local variables
!
      real(8) :: ak , akclth , apekl , aprdiv , avrgt , avrgtl , cell , &
               & cthrs , den , dentpy , dhdt , difql , diftl , dpkl ,   &
               & dpmix , dqref , drheat , dsp , dsp0k , dspbk , dsptk , &
               & dst , dstq , dtdeta , dthem , ee , efi , es , fefi ,   &
               & fptk , hcorr , otsum , pdiff , pdiffk , pflag , pk0 ,  &
               & pkb , pkl , pkt , potsum , pppk , prainx , preck ,     &
               & psfck , psum , pthrs , ptpk , qkl , qnew , qotsum ,    &
               & qrfkl , qrftp , qs , qsum , qu , rdp0t , rdpsum , rhh ,&
               & rhl , rotsum , rtbar , smix , stabdl , sumde , sumdp , &
               & sumdt , tauk , tcorr , tdpt , thskl , thtpk , thvmkl , &
               & tkl , tlcl , trfkl , tskl , ztop
      real(8) , dimension(ix,kx) :: ape , q , qqmod , t , tmod , tref , &
                                  & z0
      real(8) , dimension(kx) :: apek , apesk , difq , dift , dzq ,     &
                               & fpk , pdp , pk , psk , qk , qrefk ,    &
                               & qsatk , therk , thsk , thvref , tk ,   &
                               & trefk
      real(8) , dimension(ix) :: cldhgt , dsp0 , dspb , dspt , p ,      &
                               & pbot , prtop , psp , xsm , thbt ,      &
                               & thesp , ths , tthbt , tthes
      integer :: i , icond , iconss , iter , ivi , k , kb , kbaseb ,    &
               & kclth , khdeep , khshal , kk , l , l0 , l0m1 , lb ,    &
               & lbm1 , lbtk , lcor , lqm , lshu , ltp1 , ltpk , ltsh , &
               & n , ndeep , ndepth , ndstn , ndstp , nshal , nswap , ll
      integer , dimension(ix) :: ifbuoy , ip300 , kdeep , kshal , lbot ,&
                               & ltop , ml
      integer , dimension(kx) :: kdp , nbotd , nbots , ndpthd , ndpths ,&
                               & ntopd , ntops
      real(8), external :: tpfc

!-----------------------------------------------------------------------

      do k = 1 , kx
        do i = 1 , ilx
          cldlwc(i,k) = 0.
          cldfra(i,k) = 0.
        end do
      end do
      if ( ichem.eq.1 ) then
!
!       icumtop = top level of cumulus clouds
!       icumbot = bottom level of cumulus clouds
!       (calculated in cupara and stored for tractend)
!       before do 100 put
        do i = 2 , ilxm
          icumtop(i,j) = 0
          icumbot(i,j) = 0
        end do
      end if
      icond = 0
      iconss = 0
      tauk = dt2/trel
      cthrs = (0.006350/86400.)*dt2/cprlg
!-----------------------------------------------------------------------
!
!...  xsm is surface mask: =1 water; =0 land
      do i = 2 , ilxm
        if ( veg2d(i,j).ge.0.002 ) then
          xsm(i) = 0.
        else
          xsm(i) = 1.
        end if
      end do
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        do i = 2 , ilxm
          cldefi(i,j) = avgefi*xsm(i) + stefi*(h1-xsm(i))
        end do
      end if
!...lb is currently set to kx-1
      lb = kx - 1
      do k = 1 , kx
        ntopd(k) = 0
        nbotd(k) = 0
        ntops(k) = 0
        nbots(k) = 0
        ndpths(k) = 0
        ndpthd(k) = 0
      end do
!...find melting level...
      do i = 2 , ilxm
        ml(i) = kxp1
      end do
      do i = 2 , ilxm
        do k = 1 , kx
          t(i,k) = tb(i,k,j)/psb(i,j)
          if ( t(i,k).gt.t0 .and. ml(i).eq.kxp1 ) ml(i) = k
          q(i,k) = qvb(i,k,j)/psb(i,j)
          pppk = (a(k)*psb(i,j)+ptop)*1000.
          ape(i,k) = (pppk/h10e5)**dm2859
        end do
        lbot(i) = kx
        thesp(i) = d00
        thbt(i) = d00
        psp(i) = 9.5E4
        tref(i,1) = t(i,1)
!...ifbuoy = 0 means no positive buoyancy; ifbuoy(i) means yes...
!...ip300 is the highest model level in the lowest 300 mb...
        ifbuoy(i) = 0
        ip300(i) = 0
        cell = ptop/psb(i,j)
        do k = 1 , kx
          dzq(k) = r/g*tbase(i,k,j)                                     &
                 & *dlog((sigma(k+1)+cell)/(sigma(k)+cell))
        end do
        z0(i,kx) = 0.5*dzq(kx)
        do k = kx - 1 , 1 , -1
          z0(i,k) = z0(i,k+1) + 0.5*(dzq(k)+dzq(k+1))
        end do
      end do
!--------------padding specific humidity if too small-------------------
      do k = 1 , kx
        do i = 2 , ilxm
          if ( q(i,k).lt.epsq ) q(i,k) = epsq
          pdiff = (1.-a(k))*psb(i,j)
          if ( pdiff.lt.30. .and. ip300(i).eq.0 ) ip300(i) = k
        end do
      end do
!--------------search for maximum buoyancy level------------------------
      do kb = 1 , kx
        do i = 2 , ilxm
          pkl = (a(kb)*psb(i,j)+ptop)*1000.
          psfck = (a(kx)*psb(i,j)+ptop)*1000.
          if ( pkl.ge.psfck-pbm ) then
            tthbt(i) = t(i,kb)*ape(i,kb)
            ee = pkl*q(i,kb)/(0.622+q(i,kb))
            tdpt = 1./(d273-rv/elwv*dlog(ee/611.))
            tdpt = dmin1(tdpt,t(i,kb))
            tlcl = tdpt - (.212+1.571E-3*(tdpt-t0)-4.36E-4*(t(i,kb)-t0))&
                 & *(t(i,kb)-tdpt)
            tthes(i) = tthbt(i)*exp(elocp*q(i,kb)/tlcl)
!--------------check for maximum buoyancy-------------------------------
            if ( tthes(i).gt.thesp(i) ) then
              psp(i) = h10e5*(tthbt(i)/tlcl)**cporng
              thbt(i) = tthbt(i)
              thesp(i) = tthes(i)
            end if
          end if
!-----------------------------------------------------------------------
        end do
      end do
!---------choose cloud base as model level just below psp--------------
      do k = 1 , lm1
        ak = a(k)
        do i = 2 , ilxm
          p(i) = (ak*psb(i,j)+ptop)*1000.
!         cloud bottom cannot be above 200 mb
          if ( p(i).lt.psp(i) .and. p(i).ge.pqm ) lbot(i) = k + 1
        end do
      end do
!***  warning: lbot must not be gt kx-1 in shallow convection
!***  make sure the cloud base is at least 25 mb above the surface
      do i = 2 , ilxm
        pbot(i) = (a(lbot(i))*psb(i,j)+ptop)*1000.
        psfck = (a(kx)*psb(i,j)+ptop)*1000.
        if ( pbot(i).ge.psfck-pone .or. lbot(i).ge.kx ) then
!***      cloud bottom is at the surface so recalculate cloud bottom
          do k = 1 , lm1
            p(i) = (a(kx)*psb(i,j)+ptop)*1000.
            if ( p(i).lt.psfck-pone ) lbot(i) = k
          end do
          pbot(i) = (a(lbot(i))*psb(i,j)+ptop)*1000.
        end if
      end do
!--------------cloud top computation------------------------------------
      do i = 2 , ilxm
        prtop(i) = pbot(i)
        ltop(i) = lbot(i)
      end do
      do ivi = 1 , kx
        l = lp1 - ivi
!--------------find environmental saturation equiv pot temp...
        do i = 2 , ilxm
          p(i) = (a(l)*psb(i,j)+ptop)*1000.
          es = aliq*exp((bliq*t(i,l)-cliq)/(t(i,l)-dliq))
          qs = 0.622*es/(p(i)-es)
          ths(i) = t(i,l)*ape(i,l)*exp(elocp*qs/t(i,l))
        end do
!--------------buoyancy check-------------------------------------------
        do i = 2 , ilxm
          if ( l.le.lbot(i) ) then
            if ( thesp(i).gt.ths(i) ) ifbuoy(i) = 1
            if ( thesp(i).gt.ths(i)-1.5 .and. ifbuoy(i).eq.1 ) ltop(i)  &
               & = l + 1
          end if
        end do
!------------------------------------------------
      end do
!--------------cloud top pressure---------------------------------------
      do i = 2 , ilxm
!       if(kf(i).eq.1) goto 275
        prtop(i) = (a(ltop(i))*psb(i,j)+ptop)*1000.
      end do
!-----------------------------------------------------------------------
!--------------define and smooth dsps and cldefi------------------------
      if ( unis ) then
        do i = 2 , ilxm
          efi = cldefi(i,j)
          dspb(i) = (efi-efimn)*slopbs + dspbss
          dsp0(i) = (efi-efimn)*slop0s + dsp0ss
          dspt(i) = (efi-efimn)*slopts + dsptss
        end do
      else if ( .not.unil ) then
        do i = 2 , ilxm
          efi = cldefi(i,j)
          dspb(i) = ((efi-efimn)*slopbs+dspbss)*xsm(i)                  &
                  & + ((efi-efimn)*slopbl+dspbsl)*(h1-xsm(i))
          dsp0(i) = ((efi-efimn)*slop0s+dsp0ss)*xsm(i)                  &
                  & + ((efi-efimn)*slop0l+dsp0sl)*(h1-xsm(i))
          dspt(i) = ((efi-efimn)*slopts+dsptss)*xsm(i)                  &
                  & + ((efi-efimn)*sloptl+dsptsl)*(h1-xsm(i))
        end do
      else
        do i = 2 , ilxm
          efi = cldefi(i,j)
          dspb(i) = ((efi-efimn)*slopbl+dspbsl)
          dsp0(i) = ((efi-efimn)*slop0l+dsp0sl)
          dspt(i) = ((efi-efimn)*sloptl+dsptsl)
        end do
      end if
!--------------initialize changes of t and q due to convection----------
      do k = 1 , kx
        do i = 2 , ilxm
          tmod(i,k) = d00
          qqmod(i,k) = d00
        end do
      end do
!--------------clean up and gather deep convection points---------------
      khdeep = 0
      nswap = 0
      do i = 2 , ilxm
        if ( ltop(i).gt.lbot(i) ) then
          ltop(i) = lbot(i)
          prtop(i) = pbot(i)
        end if
        cldhgt(i) = z0(i,ltop(i)) - z0(i,lbot(i))
!       cloud is less than 90 mb deep or less than 3 sigma layers deep
        if ( cldhgt(i).lt.zno ) cldefi(i,j) = avgefi*xsm(i)             &
           & + stefi*(h1-xsm(i))
!       cloud has to be at least 290 mb deep
        if ( cldhgt(i).ge.zsh ) then
          khdeep = khdeep + 1
          kdeep(khdeep) = i
        end if
      end do
!************* horizontal loop for deep convection *********************
      do n = 1 , khdeep
        i = kdeep(n)
        dentpy = d00
        avrgt = d00
        preck = d00
        ltpk = ltop(i)
        lbtk = lbot(i)
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdc  deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
        efi = cldefi(i,j)
        dspbk = dspb(i)
        dsp0k = dsp0(i)
        dsptk = dspt(i)
!--------------initialize variables in the convective column------------
        do k = 1 , kx
          dift(k) = d00
          difq(k) = d00
          tkl = t(i,k)
          tk(k) = tkl
          trefk(k) = tkl
          qkl = q(i,k)
          qk(k) = qkl
          qrefk(k) = qkl
          pkl = (a(k)*psb(i,j)+ptop)*1000.
!**************
          tref(i,k) = tpfc(pkl,thesp(i),t(i,k),d273,elwv,qu,ape(i,k))
!***************
          pk(k) = pkl
          psk(k) = pkl
          apekl = ape(i,k)
          apek(k) = apekl
          therk(k) = tref(i,k)*apekl
        end do
!--------------deep convection reference temperature profile------------
        ltp1 = ltpk + 1
        lbm1 = lb - 1
        pkb = pk(lb)
        pkt = pk(ltpk)
!--------------temperature reference profile below freezing level-------
        l0 = lb
        pk0 = pk(lb)
        do l = ltpk , lbm1
          ivi = ltpk + lbm1 - l
          if ( trefk(ivi+1).le.t1 ) then
!--------------temperature reference profile above freezing level-------
            l0m1 = l0 - 1
            rdp0t = h1/(pk0-pkt)
            dthem = therk(l0) - trefk(l0)*apek(l0)
            do ll = ltpk , l0m1
              trefk(l) = (therk(l)-(pk(l)-pkt)*dthem*rdp0t)/apek(l)
            end do
            go to 50
          else
            stabdl = stabd
            trefk(ivi) = ((therk(ivi)-therk(ivi+1))*stabdl+trefk(ivi+1) &
                       & *apek(ivi+1))/apek(ivi)
            l0 = ivi
            pk0 = pk(l0)
          end if
        end do
!--------------freezing level at or above the cloud top-----------------
        l0m1 = l0 - 1
!--------------deep convection reference humidity profile---------------
 50     continue
        do l = ltpk , lb
!--------------saturation pressure difference---------------------------
          if ( pkb-pk0.lt.pfrz ) then
            dsp = dspc
          else if ( l.lt.l0 ) then
            dsp = ((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
          else
            dsp = ((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
          end if
!--------------humidity profile-----------------------------------------
          if ( pk(l).gt.pqm ) then
!           pressure must be below 200 mb
            psk(l) = pk(l) + dsp
            apesk(l) = (psk(l)/h10e5)**dm2859
            thsk(l) = trefk(l)*apek(l)
            qrefk(l) = pq0/psk(l)                                       &
                     & *exp(a2*(thsk(l)-a3*apesk(l))/(thsk(l)-a4*apesk  &
                     & (l)))
          else
            qrefk(l) = q(i,l)
          end if
        end do
!--------------enthalpy conservation integral--------------------------
        do iter = 1 , 2
!-----------------------------------------------------------------------
          sumde = d00
          sumdp = d00
          do l = ltpk , lb
            sumde = ((tk(l)-trefk(l))*cp+(qk(l)-qrefk(l))*elwv)         &
                  & *dsigma(l) + sumde
            sumdp = sumdp + dsigma(l)
          end do
          hcorr = sumde/(sumdp-dsigma(ltpk))
          lcor = ltpk + 1
!--------------find lqm-------------------------------------------------
          do l = 1 , lb
            if ( pk(l).le.pqm ) lqm = l
          end do
!--------------above lqm correct temperature only-----------------------
          if ( lcor.le.lqm ) then
            do l = lcor , lqm
              trefk(l) = trefk(l) + hcorr*rcp
            end do
            lcor = lqm + 1
          end if
!--------------below lqm correct both temperature and moisture----------
          do l = lcor , lb
            tskl = trefk(l)*apek(l)/apesk(l)
            dhdt = qrefk(l)*a23m4l/(tskl-a4)**2 + cp
            trefk(l) = hcorr/dhdt + trefk(l)
            thskl = trefk(l)*apek(l)
            qrefk(l) = pq0/psk(l)                                       &
                     & *exp(a2*(thskl-a3*apesk(l))/(thskl-a4*apesk(l)))
          end do
!-----------------------------------------------------------------------
        end do
        do l = 1 , kx
          thvref(l) = trefk(l)*apek(l)*(qrefk(l)*d608+h1)
        end do
!--------------heating, moistening, precipitation-----------------------
        do l = ltpk , lb
          tkl = tk(l)
          diftl = (trefk(l)-tkl)*tauk
          difql = (qrefk(l)-qk(l))*tauk
          avrgtl = (tkl+tkl+diftl)
          dentpy = (diftl*cp+difql*elwv)*dsigma(l)/avrgtl + dentpy
          avrgt = avrgtl*dsigma(l) + avrgt
          preck = dsigma(l)*diftl + preck
          dift(l) = diftl
          difq(l) = difql
        end do
        dentpy = dentpy + dentpy
        avrgt = avrgt/(sumdp+sumdp)
        if ( dentpy.lt.epsntp .or. preck.le.d00 ) then
          if ( oct90 ) then
            cldefi(i,j) = efimn
          else
            cldefi(i,j) = efimn*xsm(i) + stefi*(h1-xsm(i))
          end if
          ztop = z0(i,lbot(i)) + zsh - 0.000001
          do l = 1 , lb
            if ( z0(i,l).ge.ztop ) ltop(i) = l + 1
          end do
          prtop(i) = pk(ltop(i))
!------------cloud must be at least 2 layers thick---------------------
          if ( lbot(i)-ltop(i).lt.2 ) ltop(i) = lbot(i) - 2
          prtop(i) = pk(ltop(i))
          cldhgt(i) = z0(i,ltop(i)) - z0(i,lbot(i))
          nswap = nswap + 1
          cycle
        end if
!--------------... deep convection otherwise----------------------------
        icond = icond + 1
!***    keep the land value of efi equal to 1 until precip surpasses
!***    a threshold value, currently set to 0.25 inches per 24 hrs
        pthrs = cthrs/psb(i,j)
        drheat = (preck*xsm(i)+dmax1(epsp,preck-pthrs)*(h1-xsm(i)))     &
               & *cp/avrgt
        efi = efifc*dentpy/drheat
!vvvvv  unified or separate land/sea conv.
        if ( .not.(oct90) ) then
          efi = cldefi(i,j)*fcp + efi*fcc1
        else if ( unis ) then
          efi = cldefi(i,j)*fcp + efi*fcc1
        else if ( .not.unil ) then
          efi = (cldefi(i,j)*fcp+efi*fcc1)*xsm(i) + h1 - xsm(i)
        else
          efi = h1
        end if
!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        if ( efi.gt.h1 ) efi = h1
        if ( efi.lt.efimn ) efi = efimn
        cldefi(i,j) = efi
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        fefi = efmnt + slope*(efi-efimn)
!aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        preck = preck*fefi
!--------------update precipitation, temperature & moisture-------------
        prainx = 0.5*((psb(i,j)*1000.*preck*cprlg)*100.)
        rainc(i,j) = prainx + rainc(i,j)
!.....................precipitation rate for bats (mm/s)
        aprdiv = dble(nbatst)
        if ( jyear.eq.jyear0 .and. ktau.eq.0 ) aprdiv = 1.
        pptc(i,j) = pptc(i,j) + prainx/(dtmin*60.)/aprdiv
        do l = ltpk , lb
          tmod(i,l) = dift(l)*fefi/dt2
          qqmod(i,l) = difq(l)*fefi/dt2
        end do
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdc  end of deep convection  dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!-----------------------------------------------------------------------
      end do
      ndeep = 0
      do i = 2 , ilxm
        ltpk = ltop(i)
        lbtk = lbot(i)
        ptpk = prtop(i)
        if ( cldhgt(i).ge.zsh ) then
          ndeep = ndeep + 1
          ndepth = lb - ltpk
          ntopd(ltpk) = ntopd(ltpk) + 1
          nbotd(lb) = nbotd(lb) + 1
          if ( ndepth.gt.0 ) ndpthd(ndepth) = ndpthd(ndepth) + 1
        end if
      end do
!--------------gather shallow convection points-------------------------
      khshal = 0
      ndstn = 0
      ndstp = 0
      do i = 2 , ilxm
        if ( cldhgt(i).ge.zno .and. ltop(i).le.lbot(i)-2 ) then
          if ( cldhgt(i).lt.zsh ) then
            khshal = khshal + 1
            kshal(khshal) = i
          end if
        end if
      end do
!************* horizontal loop for shallow convection ******************
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!scscscscscsc  shallow convection  cscscscscscscscscscscscscscscscscscsc
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!-----------------------------------------------------------------------
      do n = 1 , khshal
        i = kshal(n)
        do k = 1 , kx
          tkl = t(i,k)
          tk(k) = tkl
          trefk(k) = tkl
          qkl = q(i,k)
          qk(k) = qkl
          qrefk(k) = qkl
          qsatk(k) = qkl
          pkl = (a(k)*psb(i,j)+ptop)*1000.
          pk(k) = pkl
          apekl = ape(i,k)
          apek(k) = apekl
          thvmkl = tkl*apekl*(qkl*d608+h1)
          thvref(k) = thvmkl
          pdp(k) = pk(k) - pmn
        end do
!
!...find kdp...kdp(k) is the model level closest to 65 mb (pmn) above k;
!...this is the depth over which relative humidity drop is measured to
!...estimate shallow cloud top... see do 545...
!
        do kk = kx , 1 , -1
          pflag = abs(pk(kx)-pdp(kk))
          do k = kx - 1 , 1 , -1
            pdiffk = abs(pk(k)-pdp(kk))
            if ( pdiffk.lt.pflag ) then
              pflag = pdiffk
              if ( kk.eq.k ) then
                kdp(kk) = k - 1
              else
                kdp(kk) = k
              end if
            end if
          end do
          kdp(kk) = max(1,kdp(kk))
        end do
!--------------search for shallow cloud top-----------------------------
        lbtk = lbot(i)
        ltsh = lbtk
        lbm1 = lbtk - 1
        ztop = z0(i,lbot(i)) + zsh - 0.000001
!--------------cloud top is level just above pbtk-psh ------------------
        do l = 1 , kx
          if ( z0(i,l).ge.ztop ) ltpk = l
        end do
        ptpk = pk(ltpk)
!--------------highest level allowed is level just below pshu-----------
        if ( ptpk.le.pshu ) then
          do l = 1 , kx
            if ( pk(l).le.pshu ) lshu = l + 1
          end do
          ltpk = lshu
          ptpk = pk(ltpk)
        end if
        ltp1 = ltpk + 1
!-----------------------------------------------------------------------
        do l = ltpk , lbtk
          if ( l.ge.ml(i) ) then
            es = aliq*exp((bliq*tk(l)-cliq)/(tk(l)-dliq))
          else
            es = aice*exp((bice*tk(l)-cice1)/(tk(l)-dice))
          end if
          qsatk(l) = 0.622*es/(pk(l)-es)
        end do
!-----------------------------------------------------------------------
        do l = ltp1 , lbm1
          rhl = qk(l)/qsatk(l)
          rhh = qk(kdp(l))/qsatk(kdp(l))
          if ( rhh+rhf.lt.rhl ) ltsh = l
        end do
!
        ltop(i) = ltsh
        prtop(i) = pk(ltsh)
        ltp1 = ltsh
        ltpk = ltsh - 1
        cldhgt(i) = z0(i,ltop(i)) - z0(i,lbot(i))
!       if cloud is not at least 90 mb or 3 sigma layers deep, then no
!       cloud
        if ( cldhgt(i).lt.zno .or. ltop(i).gt.lbot(i)-2 ) then
          ltop(i) = lbot(i)
          prtop(i) = pbot(i)
          cycle
        end if
!--------------scaling potential temperature & table index at top-------
        thtpk = t(i,ltp1)*ape(i,ltp1)
        pkl = (a(ltp1)*psb(i,j)+ptop)*1000.
        ee = pkl*q(i,ltp1)/(0.622+q(i,ltp1))
        tdpt = 1./(d273-rv/elwv*dlog(ee/611.))
        tdpt = dmin1(tdpt,t(i,ltp1))
        tlcl = tdpt - (.212+1.571E-3*(tdpt-t0)-4.36E-4*(t(i,ltp1)-t0))  &
             & *(t(i,ltp1)-tdpt)
        ptpk = h10e5*(thtpk/tlcl)**cporng
        dpmix = ptpk - psp(i)
        if ( abs(dpmix).lt.h3000 ) dpmix = -h3000
!--------------temperature propfile slope-------------------------------
        smix = (thtpk-thbt(i))/dpmix*stabs
        do l = ltp1 , lbtk
          ivi = ltp1 + lbtk - l
          trefk(ivi) = ((pk(ivi)-pk(ivi+1))*smix+trefk(ivi+1)           &
                     & *apek(ivi+1))/apek(ivi)
        end do
!--------------temperature reference profile correction-----------------
        sumdt = d00
        sumdp = d00
        do l = ltp1 , lbtk
          sumdt = (tk(l)-trefk(l))*dsigma(l) + sumdt
          sumdp = sumdp + dsigma(l)
        end do
!
        rdpsum = 1./sumdp
        fpk(lbtk) = trefk(lbtk)
        tcorr = sumdt*rdpsum
        do l = ltp1 , lbtk
          trfkl = trefk(l) + tcorr
          trefk(l) = trfkl
          fpk(l) = trfkl
        end do
!--------------humidity profile equations-------------------------------
        psum = 0.0
        qsum = 0.0
        potsum = 0.0
        qotsum = 0.0
        otsum = 0.0
        dst = 0.0
        fptk = fpk(ltp1)
        do l = ltp1 , lbtk
          dpkl = fpk(l) - fptk
          psum = dpkl*dsigma(l) + psum
          qsum = qk(l)*dsigma(l) + qsum
          rtbar = 2./(trefk(l)+tk(l))
          otsum = dsigma(l)*rtbar + otsum
          potsum = dpkl*rtbar*dsigma(l) + potsum
          qotsum = qk(l)*rtbar*dsigma(l) + qotsum
          dst = (trefk(l)-tk(l))*rtbar*dsigma(l) + dst
        end do
!
        psum = psum*rdpsum
        qsum = qsum*rdpsum
        rotsum = 1./otsum
        potsum = potsum*rotsum
        qotsum = qotsum*rotsum
        dst = dst*rotsum*cp/elwv
!--------------ensure positive entropy change---------------------------
        if ( dst.gt.0. ) then
          prtop(i) = pbot(i)
          ltop(i) = lbot(i)
          ndstp = ndstp + 1
          cycle
        else
          dstq = dst*epsdn
        end if
!--------------check for isothermal atmosphere--------------------------
        den = potsum - psum
        if ( -den/psum.lt.0.00005 ) then
          ltop(i) = lbot(i)
          prtop(i) = pbot(i)
          cycle
        else
!--------------slope of the reference humidity profile------------------
          dqref = (qotsum-dstq-qsum)/den
        end if
!--------------humidity doesn`t increase with height--------------------
        if ( dqref.lt.0.0 ) then
          ltop(i) = lbot(i)
          prtop(i) = pbot(i)
          cycle
        end if
!--------------humidity at the cloud top--------------------------------
        qrftp = qsum - dqref*psum
!--------------humidity profile-----------------------------------------
        do l = ltp1 , lbtk
          qrfkl = (fpk(l)-fptk)*dqref + qrftp
!--------------supersaturation not allowed------------------------------
          qnew = (qrfkl-qk(l))*tauk + qk(l)
          if ( qnew.gt.qsatk(l)*stresh ) then
            ltop(i) = lbot(i)
            prtop(i) = pbot(i)
            go to 100
          end if
!-----------------------------------------------------------------------
          thvref(l) = trefk(l)*apek(l)*(qrfkl*d608+h1)
          qrefk(l) = qrfkl
        end do
!--------------eliminate impossible slopes (betts, dtheta/dq)-----------
        do l = ltp1 , lbtk
          dtdeta = (thvref(l-1)-thvref(l))/(a(l)-a(l-1))
          if ( dtdeta.lt.epsth ) then
            ltop(i) = lbot(i)
            prtop(i) = pbot(i)
            go to 100
          end if
        end do
        if ( dst.gt.0. ) then
          ndstp = ndstp + 1
        else
          ndstn = ndstn + 1
        end if
        dentpy = d00
        do l = ltp1 , lbtk
          dentpy = ((trefk(l)-tk(l))*cp+(qrefk(l)-qk(l))*elwv)          &
                 & /(tk(l)+trefk(l))*dsigma(l) + dentpy
        end do
!--------------relaxation towards reference profiles--------------------
        iconss = iconss + 1
        do l = ltp1 , lbtk
          tmod(i,l) = (trefk(l)-tk(l))/trel
          qqmod(i,l) = (qrefk(l)-qk(l))/trel
        end do
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!scscscscscsc  end of shallow convection   scscscscscscscscscscscscscscs
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!-----------------------------------------------------------------------
      end do
 100  continue
      nshal = 0
      do i = 2 , ilxm
        ltpk = ltop(i)
        lbtk = lbot(i)
        ptpk = prtop(i)
!       no shallow convection if cloud is not at least 90 mb or 3 sigma
!       layers deep
        if ( cldhgt(i).ge.zno ) then
          if ( cldhgt(i).lt.zsh ) then
            nshal = nshal + 1
            ntops(ltpk) = ntops(ltpk) + 1
            nbots(lbtk) = nbots(lbtk) + 1
            ndepth = lbtk - ltpk
            if ( ndepth.gt.0 ) ndpths(ndepth) = ndpths(ndepth) + 1
          end if
!         find cloud fractional cover and liquid water content
          kbaseb = min0(lbtk,kx-2)
          if ( ltpk.le.kbaseb ) then
            kclth = kbaseb - ltpk + 1
            akclth = 1./dble(kclth)
            do k = ltpk , kbaseb
              cldlwc(i,k) = cllwcv
              cldfra(i,k) = 1. - (1.-clfrcv)**akclth
            end do
          end if
          if ( ichem.eq.1 ) then
            icumtop(i,j) = ltpk
            icumbot(i,j) = kbaseb
          end if
        end if
      end do
!-----------------------------------------------------------------------
      do k = 1 , kx
        do i = 2 , ilxm
          tten(i,k) = tten(i,k) + tmod(i,k)*psb(i,j)
          qten(i,k) = qten(i,k) + qqmod(i,k)*psb(i,j)
        end do
      end do
      icon(j) = icond

      end subroutine bmpara
