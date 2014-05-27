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
module module_cu_kf

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_stdio
  use mod_mpmessage

  implicit none

  private

!  public :: kfcps

  real(rk8) , parameter :: rad = 1500.0D0

  contains

#ifdef peacewise
  subroutine kfcps(                                         &
              ids,ide, jds,jde, kds,kde                     &
             ,ims,ime, jms,jme, kms,kme                     &
             ,its,ite, jts,jte, kts,kte                     &
             ,dt,ktau,dx,cudt,adapt_step_flag               &
             ,rho                                           &
             ,raincv,pratec,nca                             &
             ,u,v,th,t,w,qv,dz8w,pcps,pi                    &
             ,w0avg,xlv0,xlv1,xls0,xls1,cp,r,g,ep1          &
             ,ep2,svp1,svp2,svp3,svpt0                      &
             ,stepcu,cu_act_flag,warm_rain                  &
           ! optional arguments
             ,f_qv    ,f_qc    ,f_qr    ,f_qi    ,f_qs      &
             ,rqvcuten,rqccuten,rqrcuten,rqicuten,rqscuten  &
             ,rthcuten                                      &
                                                            )
   implicit none
   integer,      intent(in   ) ::                            &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

   integer,      intent(in   ) :: stepcu
   logical,      intent(in   ) :: warm_rain

   real,         intent(in   ) :: xlv0,xlv1,xls0,xls1
   real,         intent(in   ) :: cp,r,g,ep1,ep2
   real,         intent(in   ) :: svp1,svp2,svp3,svpt0

   integer,      intent(in   ) :: ktau

   real,  dimension( ims:ime , kms:kme , jms:jme )         , &
          intent(in   ) ::                                   &
                                                          u, &
                                                          v, &
                                                          w, &
                                                         th, &
                                                         qv, &
                                                          t, &
                                                       dz8w, &
                                                       pcps, &
                                                        rho, &
                                                         pi
!
   real,  intent(in   ) :: dt, dx
   real,  intent(in   ) :: cudt
   logical,optional, intent(in  ) :: adapt_step_flag

   real, dimension( ims:ime , jms:jme ),                     &
          intent(inout) ::                                   &
                                                     raincv  &
                                                    ,pratec  &
                                                    ,   nca

   real, dimension( ims:ime , kms:kme , jms:jme ),           &
          intent(inout) ::                                   &
                                                      w0avg

   logical, dimension( ims:ime , jms:jme ),                  &
          intent(inout) :: cu_act_flag

!
! optional arguments
!

   real, dimension( ims:ime , kms:kme , jms:jme ),           &
         optional,                                           &
         intent(inout) ::                                    &
                                                   rthcuten  &
                                                  ,rqvcuten  &
                                                  ,rqccuten  &
                                                  ,rqrcuten  &
                                                  ,rqicuten  &
                                                  ,rqscuten

!
! flags relating to the optional tendency arrays declared above
! models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!

   logical, optional ::                                      &
                                                   f_qv      &
                                                  ,f_qc      &
                                                  ,f_qr      &
                                                  ,f_qi      &
                                                  ,f_qs



! local vars

   real, dimension( kts:kte ) ::                             &
                                                        u1d, &
                                                        v1d, &
                                                        t1d, &
                                                       dz1d, &
                                                       qv1d, &
                                                        p1d, &
                                                      rho1d, &
                                                    w0avg1d

   real, dimension( kts:kte )::                              &
                                                       dqdt, &
                                                      dqidt, &
                                                      dqcdt, &
                                                      dqrdt, &
                                                      dqsdt, &
                                                       dtdt

   real    ::         tst,tv,prs,rhoe,w0,scr1,dxsq,tmp

   integer :: i,j,k,ntst,icldck

   logical :: qi_flag , qs_flag
! adjustable time step changes
   real    :: lastdt = -1.0
   real    :: w0avgfctr, w0fctr, w0den

!----------------------------------------------------------------------

!--- call cumulus parameterization
!
!...tst is the number of time steps in 10 minutes...w0avg is close to a
!...running mean vertical velocity...note that if you change tst, it wil
!...change the frequency of the convective intitiation check (see below)
!...note that the ordering of vertical layers must be reversed for w0avg
!...because the ordering is reversed in kfpara...
!
   dxsq=dx*dx
   qi_flag = .false.
   qs_flag = .false.
   if ( present( f_qi ) ) qi_flag = f_qi
   if ( present( f_qs ) ) qs_flag = f_qs

!----------------------
   ntst=stepcu
   tst=float(ntst*2)
!----------------------
!  ntst=nint(1200./(dt*2.))
!  tst=float(ntst)
!  ntst=nint(0.5*tst)
!  ntst=max0(ntst,1)
!----------------------
!  icldck=mod(ktau,ntst)
!----------------------
!  write(0,*) 'dt = ',dt,'  ktau = ',ktau,' dx = ',dx
!  write(0,*) 'cudt = ',cudt
!  write(0,*) 'adapt_step_flag = ',adapt_step_flag,' ids = ',ids
!  write(0,*) 'stepcu = ',stepcu,' warm_rain = ',warm_rain
!  write(0,*) 'f_qv = ',f_qv,' f_qc = ',f_qv
!  write(0,*) 'f_qi = ',f_qi,' f_qs = ',f_qs
!  write(0,*) 'f_qr = ',f_qr
!  stop
   if (lastdt < 0) then
      lastdt = dt
   end if

   if (adapt_step_flag) then
      w0avgfctr = 2 * max(cudt*60,dt) - dt
      w0fctr = dt
      w0den = 2 * max(cudt*60,dt)
   else
      w0avgfctr = (tst-1.)
      w0fctr = 1.
      w0den = tst
   end if

   do j = jts,jte
      do k=kts,kte
         do i= its,ite
!           scr1=-5.0e-4*g*rho(i,k,j)*(w(i,k,j)+w(i,k+1,j))
!           tv=t(i,k,j)*(1.+ep1*qv(i,k,j))
!           rhoe=pcps(i,k,j)/(r*tv)
!           w0=-101.9368*scr1/rhoe
            w0=0.5*(w(i,k,j)+w(i,k+1,j))

!           old:
!
!           w0avg(i,k,j)=(w0avg(i,k,j)*(tst-1.)+w0)/tst
!           new, to support adaptive time step:
!
            w0avg(i,k,j) = ( w0avg(i,k,j) * w0avgfctr + w0 * w0fctr ) / w0den

         end do
      end do
   end do
   lastdt = dt

     do j = jts,jte
     do i= its,ite
        cu_act_flag(i,j) = .true.
     end do
     end do

     do j = jts,jte
        do i=its,ite
! if (i.eq. 110 .and. j .eq. 59 ) then
!   write(0,*) 'nca = ',nca(i,j),' cu_act_flag = ',cu_act_flag(i,j)
!   write(0,*) 'dt = ',dt,' adapt_step_flag = ',adapt_step_flag
! end if
!        if ( nint(nca(i,j)) .gt. 0 ) then
         if ( nca(i,j) .gt. 0.5*dt ) then
            cu_act_flag(i,j) = .false.
         else

            do k=kts,kte
               dqdt(k)=0.
               dqidt(k)=0.
               dqcdt(k)=0.
               dqrdt(k)=0.
               dqsdt(k)=0.
               dtdt(k)=0.
            end do
            raincv(i,j)=0.
            pratec(i,j)=0.
!
! assign vars from 3d to 1d

            do k=kts,kte
               u1d(k) =u(i,k,j)
               v1d(k) =v(i,k,j)
               t1d(k) =t(i,k,j)
               rho1d(k) =rho(i,k,j)
               qv1d(k)=qv(i,k,j)
               p1d(k) =pcps(i,k,j)
               w0avg1d(k) =w0avg(i,k,j)
               dz1d(k)=dz8w(i,k,j)
            end do

!
            call kfpara(i, j,                       &
                 u1d,v1d,t1d,qv1d,p1d,dz1d,         &
                 w0avg1d,dt,dx,dxsq,rho1d,          &
                 xlv0,xlv1,xls0,xls1,cp,r,g,        &
                 ep2,svp1,svp2,svp3,svpt0,          &
                 dqdt,dqidt,dqcdt,dqrdt,dqsdt,dtdt, &
                 raincv,pratec,nca,                 &
                 warm_rain,qi_flag,qs_flag,         &
                 ids,ide, jds,jde, kds,kde,         &
                 ims,ime, jms,jme, kms,kme,         &
                 its,ite, jts,jte, kts,kte          )

            if ( present( rthcuten ) .and. present( rqvcuten ) ) then
              do k=kts,kte
                 rthcuten(i,k,j)=dtdt(k)/pi(i,k,j)
                 rqvcuten(i,k,j)=dqdt(k)
              end do
            end if

            if( present(rqrcuten) .and. present(rqccuten) .and. &
                present(f_qr)                                    ) then
              if ( f_qr            ) then
                 do k=kts,kte
                    rqrcuten(i,k,j)=dqrdt(k)
                    rqccuten(i,k,j)=dqcdt(k)
                 end do
               else
! this is the case for eta microphysics without 3d rain field
                 do k=kts,kte
                    rqrcuten(i,k,j)=0.
                    rqccuten(i,k,j)=dqrdt(k)+dqcdt(k)
                 end do
              end if
            end if

!......     qsten stores graupel tendency if it exists, otherise snow (v2)

            if( present( rqicuten ) .and. qi_flag )then
              do k=kts,kte
                 rqicuten(i,k,j)=dqidt(k)
              end do
            end if

            if( present ( rqscuten ) .and. qs_flag )then
              do k=kts,kte
                 rqscuten(i,k,j)=dqsdt(k)
              end do
            end if
!
         end if
       end do
     end do

   end subroutine kfcps

!-----------------------------------------------------------
   subroutine kfpara (i, j,                                &
                      u0,v0,t0,qv0,p0,dzq,w0avg1d,         &
                      dt,dx,dxsq,rho,                      &
                      xlv0,xlv1,xls0,xls1,cp,r,g,          &
                      ep2,svp1,svp2,svp3,svpt0,            &
                      dqdt,dqidt,dqcdt,dqrdt,dqsdt,dtdt,   &
                      raincv,pratec,nca,                   &
                      warm_rain,qi_flag,qs_flag,           &
                      ids,ide, jds,jde, kds,kde,           &
                      ims,ime, jms,jme, kms,kme,           &
                      its,ite, jts,jte, kts,kte            )
!-----------------------------------------------------------
      implicit none
!-----------------------------------------------------------
      integer, intent(in   ) :: ids,ide, jds,jde, kds,kde, &
                                ims,ime, jms,jme, kms,kme, &
                                its,ite, jts,jte, kts,kte, &
                                i,j
      logical, intent(in   ) :: warm_rain
      logical           :: qi_flag, qs_flag

!
      real, dimension( kts:kte ),                          &
            intent(in   ) ::                           u0, &
                                                       v0, &
                                                       t0, &
                                                      qv0, &
                                                       p0, &
                                                      rho, &
                                                      dzq, &
                                                  w0avg1d
!
      real,  intent(in   ) :: dt,dx,dxsq
!

      real,  intent(in   ) :: xlv0,xlv1,xls0,xls1,cp,r,g
      real,  intent(in   ) :: ep2,svp1,svp2,svp3,svpt0
!
      real, dimension( kts:kte ), intent(inout) ::         &
                                                     dqdt, &
                                                    dqidt, &
                                                    dqcdt, &
                                                    dqrdt, &
                                                    dqsdt, &
                                                     dtdt

      real, dimension( ims:ime , jms:jme ),                &
            intent(inout) ::                       raincv, &
                                                   pratec, &
                                                      nca
!
!...define local variables...
!
      real, dimension( kts:kte ) ::                        &
            q0,z0,tv0,tu,tvu,qu,tz,tvd,                    &
            qd,qes,thtes,tg,tvg,qg,wu,wd,w0,ems,emsd,      &
            umf,uer,udr,dmf,der,ddr,umf2,uer2,             &
            udr2,dmf2,der2,ddr2,dza,thta0,thetee,          &
            thtau,theteu,thtad,theted,qliq,qice,           &
            qlqout,qicout,pptliq,pptice,detlq,detic,       &
            detlq2,detic2,ratio,ratio2

      real, dimension( kts:kte ) ::                        &
            domgdp,exn,rhoe,tvqu,dp,rh,eqfrc,wspd,         &
            qdt,fxm,thtag,thtesg,thpa,thfxtop,             &
            thfxbot,qpa,qfxtop,qfxbot,qlpa,qlfxin,         &
            qlfxout,qipa,qifxin,qifxout,qrpa,              &
            qrfxin,qrfxout,qspa,qsfxin,qsfxout,            &
            ql0,qlg,qi0,qig,qr0,qrg,qs0,qsg

      real, dimension( kts:kte+1 ) :: omg
      real, dimension( kts:kte ) :: rainfb,snowfb

! local vars

      real    :: p00,t00,cv,b61,rlf,rhic,rhbc,pie,         &
                 ttfrz,tbfrz,c5,rate
      real    :: gdry,rocp,aliq,bliq,                      &
                 cliq,dliq,aice,bice,cice,dice
      real    :: fbfrc,p300,dpthmx,thmix,qmix,zmix,pmix,   &
                 rocpq,tmix,emix,tlog,tdpt,tlcl,tvlcl,     &
                 cporq,plcl,es,dlp,tenv,qenv,tven,tvbar,   &
                 zlcl,wkl,wabs,trppt,wsigne,dtlcl,gdt,wlcl,&
                 tvavg,qese,wtw,rholcl,au0,vmflcl,upold,   &
                 upnew,abe,wklcl,thtudl,tudl,ttemp,frc1,   &
                 qnewic,rl,r1,qnwfrz,effq,be,boterm,enterm,&
                 dzz,wsq,udlbe,rei,ee2,ud2,ttmp,f1,f2,     &
                 thttmp,qtmp,tmpliq,tmpice,tu95,tu10,ee1,  &
                 ud1,cldhgt,dptt,qnewlq,dumfdp,ee,tsat,    &
                 thta,p150,usr,vconv,timec,shsign,vws,pef, &
                 cbh,rcbh,pefcbh,peff,peff2,tder,thtmin,   &
                 dtmltd,qs,tadvec,dpdd,frc,dpt,rdd,a1,     &
                 dssdt,dtmp,t1rh,qsrh,pptflx,cpr,cndtnf,   &
                 updinc,aincm2,devdmf,ppr,rced,dpptdf,     &
                 dmflfs,dmflfs2,rced2,ddinc,aincmx,aincm1, &
                 ainc,tder2,pptfl2,fabe,stab,dtt,dtt1,     &
                 dtime,tma,tmb,tmm,bcoeff,acoeff,qvdiff,   &
                 topomg,cpm,dq,abeg,dabe,dfda,frc2,dr,     &
                 udfrc,tuc,qgs,rh0,rhg,qinit,qfnl,err2,    &
                 relerr,rlc,rls,rnc,fabeold,aincold,uefrc, &
                 ddfrc,tdc,defrc

      integer :: kx,k,kl
!
      integer :: istop,ml,l5,l4,kmix,low,                  &
                 lc,mxlayr,llfc,nlayrs,nk,                 &
                 kpbl,klcl,lcl,let,iflag,                  &
                 kfrz,nk1,ltop,nj,ltop1,                   &
                 ltopm1,lvf,kstart,kmin,lfs,               &
                 nd,nic,ldb,ldt,nd1,ndk,                   &
                 nm,lmax,ncount,noitr,                     &
                 nstep,ntc
!
      data p00,t00/1.e5,273.16/
      data cv,b61,rlf/717.,0.608,3.339e5/
      data rhic,rhbc/1.,0.90/
      data pie,ttfrz,tbfrz,c5/3.141592654,268.16,248.16,1.0723e-3/
      data rate/0.01/
!-----------------------------------------------------------
      gdry=-g/cp
      rocp=r/cp
      kl=kte
      kx=kte
!
!     aliq = 613.3
!     bliq = 17.502
!     cliq = 4780.8
!     dliq = 32.19
      aliq = svp1*1000.
      bliq = svp2
      cliq = svp2*svpt0
      dliq = svp3
      aice = 613.2
      bice = 22.452
      cice = 6133.0
      dice = 0.61
!

!...option to feed convectively generated rainwater
!...into grid-resolved rainwater (or snow/graupel)
!...field.  'fbfrc' is the fraction of available
!...precipitation to be fed back (0.0 - 1.0)...
!
      fbfrc=0.0
!
!...scheme is called once  on each north-south slice, the loop below
!...checks for the possibility of initiating parameterized
!...convection at each point within the slice
!
!...see if it is necessary to check for convective triggering at this
!...grid point. if nca>0, convection is already active at this point,
!...just feed back the tendencies saved from the time when convection
!...was initiated.  if nca<0, convection is not active
!...and you may want to check to see if it can be activated for the
!...current conditions.  in previous aplications of this scheme,
!...the variable icldck was used below to save time by only checking
!...for the possibility of convective initiation every 5 or 10
!...minutes...
!

!  10 continue
!sue  p300=1000.*(psb(i,j)*a(kl)+ptop-30.)+pp3d(i,j,kl)

      p300=p0(1)-30000.
!
!...pressure perturbation term is only defined at mid-point of
!...vertical layers...since total pressure is needed at the top and
!...bottom of layers below, do an interpolation...
!
!...input a vertical sounding ... note that model layers are numbered
!...from bottom-up in the kf scheme...
!
      ml=0
!sue  tmprpsb=1./psb(i,j)
!sue  cell=ptop*tmprpsb

      do 15 k=1,kx
!sue     p0(k)=1.e3*(a(nk)*psb(i,j)+ptop)+pp3d(i,j,nk)
!
!...if q0 is above saturation value, reduce it to saturation level...
!
         es=aliq*exp((bliq*t0(k)-cliq)/(t0(k)-dliq))
         qes(k)=ep2*es/(p0(k)-es)
         q0(k)=amin1(qes(k),qv0(k))
         q0(k)=amax1(0.000001,q0(k))
         ql0(k)=0.
         qi0(k)=0.
         qr0(k)=0.
         qs0(k)=0.

         tv0(k)=t0(k)*(1.+b61*q0(k))
         rhoe(k)=p0(k)/(r*tv0(k))

         dp(k)=rho(k)*g*dzq(k)
!
!...dzq is dz between sigma surfaces, dza is dz between model half level
!   dp is the pressure interval between full sigma levels...
!
         if(p0(k).ge.500e2)l5=k
         if(p0(k).ge.400e2)l4=k
         if(p0(k).ge.p300)llfc=k
         if(t0(k).gt.t00)ml=k
   15   continue

        z0(1)=.5*dzq(1)
        do 20 k=2,kl
          z0(k)=z0(k-1)+.5*(dzq(k)+dzq(k-1))
          dza(k-1)=z0(k)-z0(k-1)
   20   continue
        dza(kl)=0.
        kmix=1
   25   low=kmix

        if(low.gt.llfc)goto 325

        lc=low
        mxlayr=0
!
!...assume that in order to support a deep updraft you need a layer of
!...unstable air 50 to 100 mb deep...to approximate this, isolate a
!...group of adjacent individual model layers, with the base at level
!...lc, such that the combined depth of these layers is at least 60 mb..
!
        nlayrs=0
        dpthmx=0.
        do 63 nk=lc,kx
          dpthmx=dpthmx+dp(nk)
          nlayrs=nlayrs+1
   63   if(dpthmx.gt.6.e3)goto 64
        goto 325
   64   kpbl=lc+nlayrs-1
        kmix=lc+1
   18   thmix=0.
        qmix=0.
        zmix=0.
        pmix=0.
        dpthmx=0.
!
!...find the thermodynamic characteristics of the layer by
!...mass-weighting the characteristics of the individual model
!...layers...
!
        do 17 nk=lc,kpbl
          dpthmx=dpthmx+dp(nk)
          rocpq=0.2854*(1.-0.28*q0(nk))
          thmix=thmix+dp(nk)*t0(nk)*(p00/p0(nk))**rocpq
          qmix=qmix+dp(nk)*q0(nk)
          zmix=zmix+dp(nk)*z0(nk)
   17   pmix=pmix+dp(nk)*p0(nk)
        thmix=thmix/dpthmx
        qmix=qmix/dpthmx
        zmix=zmix/dpthmx
        pmix=pmix/dpthmx
        rocpq=0.2854*(1.-0.28*qmix)
        tmix=thmix*(pmix/p00)**rocpq
        emix=qmix*pmix/(ep2+qmix)
!
!...find the temperature of the mixture at its lcl, pressure
!...level of lcl...
!
        tlog=alog(emix/aliq)
        tdpt=(cliq-dliq*tlog)/(bliq-tlog)
        tlcl=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))*(tmix-  &
             tdpt)
        tlcl=amin1(tlcl,tmix)
        tvlcl=tlcl*(1.+0.608*qmix)
        cporq=1./rocpq
        plcl=p00*(tlcl/thmix)**cporq
        do 29 nk=lc,kl
          klcl=nk
          if(plcl.ge.p0(nk))goto 35
   29   continue
        goto 325
   35   k=klcl-1
        dlp=alog(plcl/p0(k))/alog(p0(klcl)/p0(k))
!
!...estimate environmental temperature and mixing ratio at the lcl...
!
        tenv=t0(k)+(t0(klcl)-t0(k))*dlp
        qenv=q0(k)+(q0(klcl)-q0(k))*dlp
        tven=tenv*(1.+0.608*qenv)
        tvbar=0.5*(tv0(k)+tven)
!        zlcl=z0(k)+r*tvbar*alog(p0(k)/plcl)/g
        zlcl=z0(k)+(z0(klcl)-z0(k))*dlp
!
!...check to see if cloud is buoyant using fritsch-chappell trigger
!...function described in kain and fritsch (1992)...w0avg is an
!...aproximate value for the running-mean grid-scale vertical
!...velocity, which gives smoother fields of convective initiation
!...than the instantaneous value...formula relating temperature
!...perturbation to vertical velocity has been used with the most
!...success at grid lengths near 25 km.  for different grid-lengths,
!...adjust vertical velocity to equivalent value for 25 km grid
!...length, assuming linear dependence of w on grid length...
!
        wklcl=0.02*zlcl/2.5e3
        wkl=(w0avg1d(k)+(w0avg1d(klcl)-w0avg1d(k))*dlp)*dx/25.e3- &
            wklcl
        wabs=abs(wkl)+1.e-10
        wsigne=wkl/wabs
        dtlcl=4.64*wsigne*wabs**0.33
        gdt=g*dtlcl*(zlcl-z0(lc))/(tv0(lc)+tven)
        wlcl=1.+.5*wsigne*sqrt(abs(gdt)+1.e-10)
        if(tlcl+dtlcl.gt.tenv)goto 45
        if(kpbl.ge.llfc)goto 325
        goto 25
!
!...convective triggering criteria has been satisfied...compute
!...equivalent potential temperature
!...(theteu) and vertical velocity of the rising parcel at the lcl...
!
   45   theteu(k)=tmix*(1.e5/pmix)**(0.2854*(1.-0.28*qmix))* &
                  exp((3374.6525/tlcl-2.5403)*qmix*(1.+0.81*qmix))
        es=aliq*exp((tenv*bliq-cliq)/(tenv-dliq))
        tvavg=0.5*(tv0(klcl)+tenv*(1.+0.608*qenv))
        plcl=p0(klcl)*exp(g/(r*tvavg)*(z0(klcl)-zlcl))
        qese=ep2*es/(plcl-es)
        gdt=g*dtlcl*(zlcl-z0(lc))/(tv0(lc)+tven)
        wlcl=1.+.5*wsigne*sqrt(abs(gdt)+1.e-10)
        thtes(k)=tenv*(1.e5/plcl)**(0.2854*(1.-0.28*qese))* &
                 exp((3374.6525/tenv-2.5403)*qese*(1.+0.81*qese))
        wtw=wlcl*wlcl
        if(wlcl.lt.0.)goto 25
        tvlcl=tlcl*(1.+0.608*qmix)
        rholcl=plcl/(r*tvlcl)
!
        lcl=klcl
        let=lcl
!
!*******************************************************************
!                                                                  *
!                 compute updraft properties                       *
!                                                                  *
!*******************************************************************
!
!
!...estimate initial updraft mass flux (umf(k))...
!
        wu(k)=wlcl
        au0=pie*rad*rad
        umf(k)=rholcl*au0
        vmflcl=umf(k)
        upold=vmflcl
        upnew=upold
!
!...ratio2 is the degree of glaciation in the cloud (0 to 1),
!...uer is the envir entrainment rate, abe is available buoyant energy,
!   trppt is the total rate of precipitation production...
!
        ratio2(k)=0.
        uer(k)=0.
        abe=0.
        trppt=0.
        tu(k)=tlcl
        tvu(k)=tvlcl
        qu(k)=qmix
        eqfrc(k)=1.
        qliq(k)=0.
        qice(k)=0.
        qlqout(k)=0.
        qicout(k)=0.
        detlq(k)=0.
        detic(k)=0.
        pptliq(k)=0.
        pptice(k)=0.
        iflag=0
        kfrz=lc
!
!...the amount of conv avail pot energy (cape) is calculated with
!   respect to undilute parcel ascent; eq pot temp of undilute
!   parcel is thtudl, undilute temperature is given by tudl...
!
        thtudl=theteu(k)
        tudl=tlcl
!
!...ttemp is used during calculation of the linear glaciation
!   process; it is initially set to the temperature at which
!   freezing is specified to begin.  within the glaciation
!   interval, it is set equal to the updraft temp at the
!   previous model level...
!
        ttemp=ttfrz
!
!...enter the loop for updraft calculations...calculate updraft temp,
!   mixing ratio, vertical mass flux, lateral detrainment of mass and
!   moisture, precipitation rates at each model level...
!
        do 60 nk=k,kl-1
          nk1=nk+1
          ratio2(nk1)=ratio2(nk)
!
!...update updraft properties at the next model lvl to reflect
!   entrainment of environmental air...
!
          frc1=0.
          tu(nk1)=t0(nk1)
          theteu(nk1)=theteu(nk)
          qu(nk1)=qu(nk)
          qliq(nk1)=qliq(nk)
          qice(nk1)=qice(nk)

          call tpmix(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1), &
               qice(nk1),qnewlq,qnewic,ratio2(nk1),rl,xlv0,xlv1,xls0, &
               xls1,ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          tvu(nk1)=tu(nk1)*(1.+0.608*qu(nk1))
!
!...check to see if updraft temp is within the freezing interval,
!   if it is, calculate the fractional conversion to glaciation
!   and adjust qnewlq to reflect the gradual change in thetau
!   since the last model level...the glaciation effects will be
!   determined after the amount of condensate available after
!   precip fallout is determined...ttfrz is the temp at which
!   glaciation begins, tbfrz the temp at which it ends...
!
          if(tu(nk1).le.ttfrz.and.iflag.lt.1)then
            if(tu(nk1).gt.tbfrz)then
              if(ttemp.gt.ttfrz)ttemp=ttfrz
              frc1=(ttemp-tu(nk1))/(ttfrz-tbfrz)
              r1=(ttemp-tu(nk1))/(ttemp-tbfrz)
            else
              frc1=(ttemp-tbfrz)/(ttfrz-tbfrz)
              r1=1.
              iflag=1
            end if
            qnwfrz=qnewlq
            qnewic=qnewic+qnewlq*r1*0.5
            qnewlq=qnewlq-qnewlq*r1*0.5
            effq=(ttfrz-tbfrz)/(ttemp-tbfrz)
            ttemp=tu(nk1)
          end if
!
!  calculate updraft vertical velocity and precipitation fallout...
!
          if(nk.eq.k)then
            be=(tvlcl+tvu(nk1))/(tven+tv0(nk1))-1.
            boterm=2.*(z0(nk1)-zlcl)*g*be/1.5
            enterm=0.
            dzz=z0(nk1)-zlcl
          else
            be=(tvu(nk)+tvu(nk1))/(tv0(nk)+tv0(nk1))-1.
            boterm=2.*dza(nk)*g*be/1.5
            enterm=2.*uer(nk)*wtw/upold
            dzz=dza(nk)
          end if
          wsq=wtw
          call condload(qliq(nk1),qice(nk1),wtw,dzz,boterm,enterm,rate, &
               qnewlq,qnewic,qlqout(nk1),qicout(nk1), g)

!...if vert velocity is less than zero, exit the updraft loop and,
!   if cloud is tall enough, finalize updraft calculations...
!
          if(wtw.le.0.)goto 65
          wabs=sqrt(abs(wtw))
          wu(nk1)=wtw/wabs
!
!  update the abe for undilute ascent...
!
          thtes(nk1)=t0(nk1)*(1.e5/p0(nk1))**(0.2854*(1.-0.28*qes(nk1))) &
                     *                                                   &
                     exp((3374.6525/t0(nk1)-2.5403)*qes(nk1)*(1.+0.81*   &
                     qes(nk1)))
          udlbe=((2.*thtudl)/(thtes(nk)+thtes(nk1))-1.)*dzz
          if(udlbe.gt.0.)abe=abe+udlbe*g
!
!  determine the effects of cloud glaciation if within the specified
!  temp interval...
!
          if(frc1.gt.1.e-6)then
            call dtfrznew(tu(nk1),p0(nk1),theteu(nk1),qu(nk1),qliq(nk1), &
                 qice(nk1),ratio2(nk1),ttfrz,tbfrz,qnwfrz,rl,frc1,effq,  &
                 iflag,xlv0,xlv1,xls0,xls1,ep2,aliq,bliq,cliq,dliq,aice,bice &
                 ,cice,dice)
          end if
!
!  call subroutine to calculate environmental equivalent potential temp.
!  within glaciation interval, thetae must be calculated with respect to
!  same degree of glaciation for all entraining air...
!
          call envirtht(p0(nk1),t0(nk1),q0(nk1),thetee(nk1),ratio2(nk1), &
               rl,ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)

!...rei is the rate of environmental inflow...
!
          rei=vmflcl*dp(nk1)*0.03/rad
          tvqu(nk1)=tu(nk1)*(1.+0.608*qu(nk1)-qliq(nk1)-qice(nk1))
!
!...if cloud parcels are virtually colder than the environment, no
!   entrainment is allowed at this level...
!
          if(tvqu(nk1).le.tv0(nk1))then
            uer(nk1)=0.0
            udr(nk1)=rei
            ee2=0.
            ud2=1.
            eqfrc(nk1)=0.
            goto 55
          end if
          let=nk1
          ttmp=tvqu(nk1)
!
!...determine the critical mixed fraction of updraft and environmental
!   air for estimation of entrainment and detrainment rates...
!
          f1=0.95
          f2=1.-f1
          thttmp=f1*thetee(nk1)+f2*theteu(nk1)
          qtmp=f1*q0(nk1)+f2*qu(nk1)
          tmpliq=f2*qliq(nk1)
          tmpice=f2*qice(nk1)
          call tpmix(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,qnewlq,      &
               qnewic,ratio2(nk1),rl,xlv0,xlv1,xls0,xls1,ep2,aliq,bliq,cliq, &
               dliq,aice,bice,cice,dice)
          tu95=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
          if(tu95.gt.tv0(nk1))then
            ee2=1.
            ud2=0.
            eqfrc(nk1)=1.0
            goto 50
          end if
          f1=0.10
          f2=1.-f1
          thttmp=f1*thetee(nk1)+f2*theteu(nk1)
          qtmp=f1*q0(nk1)+f2*qu(nk1)
          tmpliq=f2*qliq(nk1)
          tmpice=f2*qice(nk1)
          call tpmix(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,qnewlq,      &
               qnewic,ratio2(nk1),rl,xlv0,xlv1,xls0,xls1,ep2,aliq,bliq,cliq, &
               dliq,aice,bice,cice,dice)
          tu10=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
          if(tu10.eq.tvqu(nk1))then
            ee2=1.
            ud2=0.
            eqfrc(nk1)=1.0
            goto 50
          end if
          eqfrc(nk1)=(tv0(nk1)-tvqu(nk1))*f1/(tu10-tvqu(nk1))
          eqfrc(nk1)=amax1(0.0,eqfrc(nk1))
          eqfrc(nk1)=amin1(1.0,eqfrc(nk1))
          if(eqfrc(nk1).eq.1)then
            ee2=1.
            ud2=0.
            goto 50
          elseif(eqfrc(nk1).eq.0.)then
            ee2=0.
            ud2=1.
            goto 50
          else
!
!...subroutine prof5 integrates over the gaussian dist to determine the
!   fractional entrainment and detrainment rates...
!
            call prof5(eqfrc(nk1),ee2,ud2)
          end if
!
   50     if(nk.eq.k)then
            ee1=1.
            ud1=0.
          end if
!
!...net entrainment and detrainment rates are given by the average
!   fractional values in the layer...
!
          uer(nk1)=0.5*rei*(ee1+ee2)
          udr(nk1)=0.5*rei*(ud1+ud2)
!
!...if the calculated updraft detrainment rate is greater than the total
!   updraft mass flux, all cloud mass detrains, exit updraft calculation
!
   55     if(umf(nk)-udr(nk1).lt.10.)then
!
!...if the calculated detrained mass flux is greater than the total
!   updraft flux, impose total detrainment of updraft mass at the
!   previous model
!
            if(udlbe.gt.0.)abe=abe-udlbe*g
            let=nk
!         write(98,1015)p0(nk1)/100.
            goto 65
          end if
          ee1=ee2
          ud1=ud2
          upold=umf(nk)-udr(nk1)
          upnew=upold+uer(nk1)
          umf(nk1)=upnew
!
!...detlq and detic are the rates of detrainment of liquid and ice in
!   the detraining updraft mass...
!
          detlq(nk1)=qliq(nk1)*udr(nk1)
          detic(nk1)=qice(nk1)*udr(nk1)
          qdt(nk1)=qu(nk1)
          qu(nk1)=(upold*qu(nk1)+uer(nk1)*q0(nk1))/upnew
          theteu(nk1)=(theteu(nk1)*upold+thetee(nk1)*uer(nk1))/upnew
          qliq(nk1)=qliq(nk1)*upold/upnew
          qice(nk1)=qice(nk1)*upold/upnew
!
!...kfrz is the highest model level at which liquid condensate is
!   generating pptliq is the rate of generation (fallout) of liquid
!   precip at a giving model lvl, pptice the same for ice, trppt is
!   the total rate of production of precip up to the current model level
!
          if(abs(ratio2(nk1)-1.).gt.1.e-6)kfrz=nk1
          pptliq(nk1)=qlqout(nk1)*(umf(nk)-udr(nk1))
          pptice(nk1)=qicout(nk1)*(umf(nk)-udr(nk1))
          trppt=trppt+pptliq(nk1)+pptice(nk1)
          if(nk1.le.kpbl)uer(nk1)=uer(nk1)+vmflcl*dp(nk1)/dpthmx
   60   continue
!
!...check cloud depth...if cloud is tall enough, estimate the equilibriu
!   temperature level (let) and adjust mass flux profile at cloud top so
!   that mass flux decreases to zero as a linear function of pressure
!   between the let and cloud top...
!
!...ltop is the model level just below the level at which vertical
!   velocity first becomes negative...
!
   65   ltop=nk
        cldhgt=z0(ltop)-zlcl
!
!...if cloud top hgt is less than specified minimum height, go back and
!   the next highest 60mb layer to see if a bigger cloud can be obtained
!   that source air...
!
!       if(cldhgt.lt.4.e3.or.abe.lt.1.)then
        if(cldhgt.lt.3.e3.or.abe.lt.1.)then
          do 70 nk=k,ltop
            umf(nk)=0.
            udr(nk)=0.
            uer(nk)=0.
            detlq(nk)=0.
            detic(nk)=0.
            pptliq(nk)=0.
   70     pptice(nk)=0.
          goto 25
        end if
!
!...if the let and ltop are the same, detrain all of the updraft mass
!   flux this level...
!
        if(let.eq.ltop)then
          udr(ltop)=umf(ltop)+udr(ltop)-uer(ltop)
          detlq(ltop)=qliq(ltop)*udr(ltop)*upnew/upold
          detic(ltop)=qice(ltop)*udr(ltop)*upnew/upold
          trppt=trppt-(pptliq(ltop)+pptice(ltop))
          uer(ltop)=0.
          umf(ltop)=0.
          goto 85
        end if
!
!   begin total detrainment at the level above the let...
!
        dptt=0.
        do 71 nj=let+1,ltop
   71   dptt=dptt+dp(nj)
        dumfdp=umf(let)/dptt
!
!...adjust mass flux profiles, detrainment rates, and precipitation fall
!   rates to reflect the linear decrease in mass flx between the let and
!   ptop
!
        do 75 nk=let+1,ltop
          udr(nk)=dp(nk)*dumfdp
          umf(nk)=umf(nk-1)-udr(nk)
          detlq(nk)=qliq(nk)*udr(nk)
          detic(nk)=qice(nk)*udr(nk)
          trppt=trppt-pptliq(nk)-pptice(nk)
          pptliq(nk)=(umf(nk-1)-udr(nk))*qlqout(nk)
          pptice(nk)=(umf(nk-1)-udr(nk))*qicout(nk)
          trppt=trppt+pptliq(nk)+pptice(nk)
   75   continue
!
!...send updraft characteristics to output files...
!
   85   continue
!
!...extend the updraft mass flux profile down to the source layer for
!   the updraft air...also, define thetae for levels below the lcl...
!
        do 90 nk=1,k
          if(nk.ge.lc)then
            if(nk.eq.lc)then
              umf(nk)=vmflcl*dp(nk)/dpthmx
              uer(nk)=vmflcl*dp(nk)/dpthmx
            elseif(nk.le.kpbl)then
              uer(nk)=vmflcl*dp(nk)/dpthmx
              umf(nk)=umf(nk-1)+uer(nk)
            else
              umf(nk)=vmflcl
              uer(nk)=0.
            end if
            tu(nk)=tmix+(z0(nk)-zmix)*gdry
            qu(nk)=qmix
            wu(nk)=wlcl
          else
            tu(nk)=0.
            qu(nk)=0.
            umf(nk)=0.
            wu(nk)=0.
            uer(nk)=0.
          end if
          udr(nk)=0.
          qdt(nk)=0.
          qliq(nk)=0.
          qice(nk)=0.
          qlqout(nk)=0.
          qicout(nk)=0.
          pptliq(nk)=0.
          pptice(nk)=0.
          detlq(nk)=0.
          detic(nk)=0.
          ratio2(nk)=0.
          ee=q0(nk)*p0(nk)/(ep2+q0(nk))
          tlog=alog(ee/aliq)
          tdpt=(cliq-dliq*tlog)/(bliq-tlog)
          tsat=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(t0(nk)-t00))*( &
               t0(nk)-tdpt)
          thta=t0(nk)*(1.e5/p0(nk))**(0.2854*(1.-0.28*q0(nk)))
          thetee(nk)=thta*                                               &
                     exp((3374.6525/tsat-2.5403)*q0(nk)*(1.+0.81*q0(nk)) &
                     )
          thtes(nk)=thta*                                                &
                    exp((3374.6525/t0(nk)-2.5403)*qes(nk)*(1.+0.81*      &
                    qes(nk)))
          eqfrc(nk)=1.0
   90   continue
!
        ltop1=ltop+1
        ltopm1=ltop-1
!
!...define variables above cloud top...
!
        do 95 nk=ltop1,kx
          umf(nk)=0.
          udr(nk)=0.
          uer(nk)=0.
          qdt(nk)=0.
          qliq(nk)=0.
          qice(nk)=0.
          qlqout(nk)=0.
          qicout(nk)=0.
          detlq(nk)=0.
          detic(nk)=0.
          pptliq(nk)=0.
          pptice(nk)=0.
          if(nk.gt.ltop1)then
            tu(nk)=0.
            qu(nk)=0.
            wu(nk)=0.
          end if
          thta0(nk)=0.
          thtau(nk)=0.
          ems(nk)=dp(nk)*dxsq/g
          emsd(nk)=1./ems(nk)
          tg(nk)=t0(nk)
          qg(nk)=q0(nk)
          qlg(nk)=0.
          qig(nk)=0.
          qrg(nk)=0.
          qsg(nk)=0.
   95   omg(nk)=0.
        omg(kl+1)=0.
        p150=p0(klcl)-1.50e4
        do 100 nk=1,ltop
          thtad(nk)=0.
          ems(nk)=dp(nk)*dxsq/g
          emsd(nk)=1./ems(nk)
!
!...initialize some variables to be used later in the vert advection
!   scheme
!
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qdt(nk)))
          thtau(nk)=tu(nk)*exn(nk)
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*q0(nk)))
          thta0(nk)=t0(nk)*exn(nk)
!
!...lvf is the level at which moisture flux is estimated as the basis
!...for precipitation efficiency calculations...
!
          if(p0(nk).gt.p150)lvf=nk
  100   omg(nk)=0.
        lvf=min0(lvf,let)
        usr=umf(lvf+1)*(qu(lvf+1)+qliq(lvf+1)+qice(lvf+1))
        usr=amin1(usr,trppt)
        if(usr.lt.1.e-8)usr=trppt
!
!     write(98,1025)klcl,zlcl,dtlcl,ltop,p0(ltop),iflag,
!    * tmix-t00,pmix,qmix,abe
!     write(98,1030)p0(let)/100.,p0(ltop)/100.,vmflcl,plcl/100.,
!    * wlcl,cldhgt
!
!...compute convective time scale(timec). the mean wind at the lcl
!...and midtroposphere is used.
!
        wspd(klcl)=sqrt(u0(klcl)*u0(klcl)+v0(klcl)*v0(klcl))
        wspd(l5)=sqrt(u0(l5)*u0(l5)+v0(l5)*v0(l5))
        wspd(ltop)=sqrt(u0(ltop)*u0(ltop)+v0(ltop)*v0(ltop))
        vconv=.5*(wspd(klcl)+wspd(l5))
        if (vconv .gt. 0.) then
           timec=dx/vconv
        else
           timec=3600.
        end if
!       timec=dx/vconv
        tadvec=timec
        timec=amax1(1800.,timec)
        timec=amin1(3600.,timec)
        nic=nint(timec/dt)
        timec=float(nic)*dt
!
!...compute wind shear and precipitation efficiency.
!
!        shsign = cvmgt(1.,-1.,wspd(ltop).gt.wspd(klcl))
        if(wspd(ltop).gt.wspd(klcl))then
          shsign=1.
        else
          shsign=-1.
        end if
        vws=(u0(ltop)-u0(klcl))*(u0(ltop)-u0(klcl))+(v0(ltop)-v0(klcl))* &
            (v0(ltop)-v0(klcl))
        vws=1.e3*shsign*sqrt(vws)/(z0(ltop)-z0(lcl))
        pef=1.591+vws*(-.639+vws*(9.53e-2-vws*4.96e-3))
        pef=amax1(pef,.2)
        pef=amin1(pef,.9)
!
!...precipitation efficiency is a function of the height of cloud base.
!
        cbh=(zlcl-z0(1))*3.281e-3
        if(cbh.lt.3.)then
          rcbh=.02
        else
          rcbh=.96729352+cbh*(-.70034167+cbh*(.162179896+cbh*(-  &
               1.2569798e-2+cbh*(4.2772e-4-cbh*5.44e-6))))
        end if
        if(cbh.gt.25)rcbh=2.4
        pefcbh=1./(1.+rcbh)
        pefcbh=amin1(pefcbh,.9)
!
!... mean pef. is used to compute rainfall.
!
        peff=.5*(pef+pefcbh)
        peff2=peff
!        write(98,1035)pef,pefcbh,lc,let,wkl,vws
!
!*****************************************************************
!                                                                *
!                  compute downdraft properties                  *
!                                                                *
!*****************************************************************
!
!...let downdraft originate at the level of minimum saturation equivalen
!...potential temperature (seqt) in the cloud layer, extend downward to
!...surface, or to the layer below cloud base at which envir seqt is les
!...than min seqt in the cloud layer...let downdraft detrain over a laye
!...of specified pressure-depth (dpdd)...
!
        tder=0.
        kstart=max0(kpbl,klcl)
        thtmin=thtes(kstart+1)
        kmin=kstart+1
        do 104 nk=kstart+2,ltop-1
          thtmin=amin1(thtmin,thtes(nk))
          if(thtmin.eq.thtes(nk))kmin=nk
  104   continue
        lfs=kmin
        if(ratio2(lfs).gt.0.)call envirtht(p0(lfs),t0(lfs),q0(lfs),     &
          thetee(lfs),0.,rl,ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
        eqfrc(lfs)=(thtes(lfs)-theteu(lfs))/(thetee(lfs)-theteu(lfs))
        eqfrc(lfs)=amax1(eqfrc(lfs),0.)
        eqfrc(lfs)=amin1(eqfrc(lfs),1.)
        theted(lfs)=thtes(lfs)
!
!...estimate the effect of melting precipitation in the downdraft...
!
        if(ml.gt.0)then
          dtmltd=0.5*(qu(klcl)-qu(ltop))*rlf/cp
        else
          dtmltd=0.
        end if
        tz(lfs)=t0(lfs)-dtmltd
        es=aliq*exp((tz(lfs)*bliq-cliq)/(tz(lfs)-dliq))
        qs=ep2*es/(p0(lfs)-es)
        qd(lfs)=eqfrc(lfs)*q0(lfs)+(1.-eqfrc(lfs))*qu(lfs)
        thtad(lfs)=tz(lfs)*(p00/p0(lfs))**(0.2854*(1.-0.28*qd(lfs)))
        if(qd(lfs).ge.qs)then
          theted(lfs)=thtad(lfs)*                                       &
                      exp((3374.6525/tz(lfs)-2.5403)*qs*(1.+0.81*qs))
        else
          call envirtht(p0(lfs),tz(lfs),qd(lfs),theted(lfs),0.,rl,ep2,aliq, &
               bliq,cliq,dliq,aice,bice,cice,dice)
        end if
        do 107 nk=1,lfs
          nd=lfs-nk
          if(theted(lfs).gt.thtes(nd).or.nd.eq.1)then
            ldb=nd
!
!...if downdraft never becomes negatively buoyant or if it
!...is shallower 50 mb, don't allow it to occur at all...
!
            if(nk.eq.1.or.(p0(ldb)-p0(lfs)).lt.50.e2)goto 141
! testing ---- no downdraft
!           goto 141
            goto 110
          end if
  107   continue
!
!...allow downdraft to detrain in a single layer, but with downdraft air
!...typically flushed up into higher layers as allowed in the total
!...vertical advection calculations farther down in the code...
!
  110   dpdd=dp(ldb)
        ldt=ldb
        frc=1.
        dpt=0.
!      do 115 nk=ldb,lfs
!        dpt=dpt+dp(nk)
!        if(dpt.gt.dpdd)then
!          ldt=nk
!          frc=(dpdd+dp(nk)-dpt)/dp(nk)
!          goto 120
!        end if
!        if(nk.eq.lfs-1)then
!         ldt=nk
!        frc=1.
!        dpdd=dpt
!        goto 120
!        end if
!115   continue
  120   continue
!
!...take a first guess at the initial downdraft mass flux..
!
        tvd(lfs)=t0(lfs)*(1.+0.608*qes(lfs))
        rdd=p0(lfs)/(r*tvd(lfs))
        a1=(1.-peff)*au0
        dmf(lfs)=-a1*rdd
        der(lfs)=eqfrc(lfs)*dmf(lfs)
        ddr(lfs)=0.
        do 140 nd=lfs-1,ldb,-1
          nd1=nd+1
          if(nd.le.ldt)then
            der(nd)=0.
            ddr(nd)=-dmf(ldt+1)*dp(nd)*frc/dpdd
            dmf(nd)=dmf(nd1)+ddr(nd)
            frc=1.
            theted(nd)=theted(nd1)
            qd(nd)=qd(nd1)
          else
            der(nd)=dmf(lfs)*0.03*dp(nd)/rad
            ddr(nd)=0.
            dmf(nd)=dmf(nd1)+der(nd)
            if(ratio2(nd).gt.0.)call envirtht(p0(nd),t0(nd),q0(nd),      &
              thetee(nd),0.,rl,ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
            theted(nd)=(theted(nd1)*dmf(nd1)+thetee(nd)*der(nd))/dmf(nd)
            qd(nd)=(qd(nd1)*dmf(nd1)+q0(nd)*der(nd))/dmf(nd)
          end if
  140   continue
        tder=0.
!
!...calculation an evaporation rate for given mass flux...
!
        do 135 nd=ldb,ldt
          tz(nd)=                                                        &
                 tpdd(p0(nd),theted(ldt),t0(nd),qs,qd(nd),1.0,xlv0,xlv1, &
                 ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
          es=aliq*exp((tz(nd)*bliq-cliq)/(tz(nd)-dliq))
          qs=ep2*es/(p0(nd)-es)
          dssdt=(cliq-bliq*dliq)/((tz(nd)-dliq)*(tz(nd)-dliq))
          rl=xlv0-xlv1*tz(nd)
          dtmp=rl*qs*(1.-rhbc)/(cp+rl*rhbc*qs*dssdt)
          t1rh=tz(nd)+dtmp
          es=rhbc*aliq*exp((bliq*t1rh-cliq)/(t1rh-dliq))
          qsrh=ep2*es/(p0(nd)-es)
!
!...check to see if mixing ratio at specified rh is less than actual
!...mixing ratio...if so, adjust to give zero evaporation...
!
          if(qsrh.lt.qd(nd))then
            qsrh=qd(nd)
!          t1rh=t1+(qs-qsrh)*rl/cp
            t1rh=tz(nd)
          end if
          tz(nd)=t1rh
          qs=qsrh
          tder=tder+(qs-qd(nd))*ddr(nd)
          qd(nd)=qs
  135   thtad(nd)=tz(nd)*(p00/p0(nd))**(0.2854*(1.-0.28*qd(nd)))
!
!...if downdraft does not evaporate any water for specified relative
!...humidity, no downdraft is allowed...
!
  141   if(tder.lt.1.)then
!          write(98,3004)i,j
 3004       format(' ','i=',i3,2x,'j=',i3)
          pptflx=trppt
          cpr=trppt
          tder=0.
          cndtnf=0.
          updinc=1.
          ldb=lfs
          do 117 ndk=1,ltop
            dmf(ndk)=0.
            der(ndk)=0.
            ddr(ndk)=0.
            thtad(ndk)=0.
            wd(ndk)=0.
            tz(ndk)=0.
  117     qd(ndk)=0.
          aincm2=100.
          goto 165
        end if
!
!...adjust downdraft mass flux so that evaporation rate in downdraft is
!...consistent with precipitation efficiency relationship...
!
        devdmf=tder/dmf(lfs)
        ppr=0.
        pptflx=peff*usr
        rced=trppt-pptflx
!
!...ppr is the total amount of precipitation that falls  out of the
!...updraft from cloud base to the lfs...updraft mass flux will be
!...increased up to the lfs to account for updraft air mixing with
!...environmental air to the updraft, so ppr will increase
!...proportionately...
!
        do 132 nm=klcl,lfs
  132   ppr=ppr+pptliq(nm)+pptice(nm)
        if(lfs.ge.klcl)then
          dpptdf=(1.-peff)*ppr*(1.-eqfrc(lfs))/umf(lfs)
        else
          dpptdf=0.
        end if
!
!...cndtnf is the amount of condensate transferred along with updraft
!...mass the downdraft at the lfs...
!
        cndtnf=(qliq(lfs)+qice(lfs))*(1.-eqfrc(lfs))
        dmflfs=rced/(devdmf+dpptdf+cndtnf)
        if(dmflfs.gt.0.)then
          tder=0.
          goto 141
        end if
!
!...ddinc is the factor by which to increase the first-guess downdraft
!...mass flux to satisfy the precip efficiency relationship, updinc is t
!...which to increase the updraft mass flux below the lfs to account for
!...transfer of mass from updraft to downdraft...
!
!       ddinc=dmflfs/dmf(lfs)
        if(lfs.ge.klcl)then
          updinc=(umf(lfs)-(1.-eqfrc(lfs))*dmflfs)/umf(lfs)
!
!...limit updinc to less than or equal to 1.5...
!
          if(updinc.gt.1.5)then
            updinc=1.5
            dmflfs2=umf(lfs)*(updinc-1.)/(eqfrc(lfs)-1.)
            rced2=dmflfs2*(devdmf+dpptdf+cndtnf)
            pptflx=pptflx+(rced-rced2)
            peff2=pptflx/usr
            rced=rced2
            dmflfs=dmflfs2
          end if
        else
          updinc=1.
        end if
        ddinc=dmflfs/dmf(lfs)
        do 149 nk=ldb,lfs
          dmf(nk)=dmf(nk)*ddinc
          der(nk)=der(nk)*ddinc
          ddr(nk)=ddr(nk)*ddinc
  149   continue
        cpr=trppt+ppr*(updinc-1.)
        pptflx=pptflx+peff*ppr*(updinc-1.)
        peff=peff2
        tder=tder*ddinc
!
!...adjust updraft mass flux, mass detrainment rate, and liquid water an
!   detrainment rates to be consistent with the transfer of the estimate
!   from the updraft to the downdraft at the lfs...
!
        do 155 nk=lc,lfs
          umf(nk)=umf(nk)*updinc
          udr(nk)=udr(nk)*updinc
          uer(nk)=uer(nk)*updinc
          pptliq(nk)=pptliq(nk)*updinc
          pptice(nk)=pptice(nk)*updinc
          detlq(nk)=detlq(nk)*updinc
  155   detic(nk)=detic(nk)*updinc
!
!...zero out the arrays for downdraft data at levels above and below the
!...downdraft...
!
        if(ldb.gt.1)then
          do 156 nk=1,ldb-1
            dmf(nk)=0.
            der(nk)=0.
            ddr(nk)=0.
            wd(nk)=0.
            tz(nk)=0.
            qd(nk)=0.
            thtad(nk)=0.
  156     continue
        end if
        do 157 nk=lfs+1,kx
          dmf(nk)=0.
          der(nk)=0.
          ddr(nk)=0.
          wd(nk)=0.
          tz(nk)=0.
          qd(nk)=0.
          thtad(nk)=0.
  157   continue
        do 158 nk=ldt+1,lfs-1
          tz(nk)=0.
          qd(nk)=0.
  158   continue
!
!
!...set limits on the updraft and downdraft mass fluxes so that the
!   inflow into convective drafts from a given layer is no more than
!   is available in that layer initially...
!
  165   aincmx=1000.
        lmax=max0(klcl,lfs)
        do 166 nk=lc,lmax
          if((uer(nk)-der(nk)).gt.0.)aincm1=ems(nk)/((uer(nk)-der(nk))* &
            timec)
          aincmx=amin1(aincmx,aincm1)
  166   continue
        ainc=1.
        if(aincmx.lt.ainc)ainc=aincmx
!
!...save the relevent variables for a unit updrft and downdrft...they
!...will iteratively adjusted by the factor ainc to satisfy the
!...stabilization closure...
!
        ncount=0
        tder2=tder
        pptfl2=pptflx
        do 170 nk=1,ltop
          detlq2(nk)=detlq(nk)
          detic2(nk)=detic(nk)
          udr2(nk)=udr(nk)
          uer2(nk)=uer(nk)
          ddr2(nk)=ddr(nk)
          der2(nk)=der(nk)
          umf2(nk)=umf(nk)
          dmf2(nk)=dmf(nk)
  170   continue
        fabe=1.
        stab=0.95
        noitr=0
        if(ainc/aincmx.gt.0.999)then
          ncount=0
          goto 255
        end if
        istop=0
  175   ncount=ncount+1
!
!*****************************************************************
!                                                                *
!           compute properties for compensational subsidence     *
!                                                                *
!*****************************************************************
!
!...determine omega value necessary at top and bottom of each layer to
!...satisfy mass continuity...
!
  185   continue
        dtt=timec
        do 200 nk=1,ltop
          domgdp(nk)=-(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
          if(nk.gt.1)then
            omg(nk)=omg(nk-1)-dp(nk-1)*domgdp(nk-1)
            dtt1=0.75*dp(nk-1)/(abs(omg(nk))+1.e-10)
            dtt=amin1(dtt,dtt1)
          end if
  200   continue
        do 488 nk=1,ltop
          thpa(nk)=thta0(nk)
          qpa(nk)=q0(nk)
          nstep=nint(timec/dtt+1)
          dtime=timec/float(nstep)
          fxm(nk)=omg(nk)*dxsq/g
  488   continue
!
!...do an upstream/forward-in-time advection of theta, qv...
!
        do 495 ntc=1,nstep
!
!...assign theta and q values at the top and bottom of each layer based
!...sign of omega...
!
          do 493 nk=1,ltop
            thfxtop(nk)=0.
            thfxbot(nk)=0.
            qfxtop(nk)=0.
  493     qfxbot(nk)=0.
          do 494 nk=2,ltop
            if(omg(nk).le.0.)then
              thfxbot(nk)=-fxm(nk)*thpa(nk-1)
              qfxbot(nk)=-fxm(nk)*qpa(nk-1)
              thfxtop(nk-1)=thfxtop(nk-1)-thfxbot(nk)
              qfxtop(nk-1)=qfxtop(nk-1)-qfxbot(nk)
            else
              thfxbot(nk)=-fxm(nk)*thpa(nk)
              qfxbot(nk)=-fxm(nk)*qpa(nk)
              thfxtop(nk-1)=thfxtop(nk-1)-thfxbot(nk)
              qfxtop(nk-1)=qfxtop(nk-1)-qfxbot(nk)
            end if
  494     continue
!
!...update the theta and qv values at each level..
!
          do 492 nk=1,ltop
            thpa(nk)=thpa(nk)+(thfxbot(nk)+udr(nk)*thtau(nk)+ddr(nk)*     &
                     thtad(nk)+thfxtop(nk)-(uer(nk)-der(nk))*thta0(nk))* &
                     dtime*emsd(nk)
            qpa(nk)=qpa(nk)+(qfxbot(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk)+   &
                    qfxtop(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)

  492     continue
  495   continue
        do 498 nk=1,ltop
          thtag(nk)=thpa(nk)
          qg(nk)=qpa(nk)
  498   continue
!
!...check to see if mixing ratio dips below zero anywhere;  if so,
!...borrow moisture from adjacent layers to bring it back up above zero.
!
        do 499 nk=1,ltop
          if(qg(nk).lt.0.)then
            if(nk.eq.1)then
              call wrf_error_fatal ( 'module_cu_kf.f: problem with kf scheme:  qg = 0 at the surface' )
            end if
            nk1=nk+1
            if(nk.eq.ltop)nk1=klcl
            tma=qg(nk1)*ems(nk1)
            tmb=qg(nk-1)*ems(nk-1)
            tmm=(qg(nk)-1.e-9)*ems(nk)
            bcoeff=-tmm/((tma*tma)/tmb+tmb)
            acoeff=bcoeff*tma/tmb
            tmb=tmb*(1.-bcoeff)
            tma=tma*(1.-acoeff)
            if(nk.eq.ltop)then
              qvdiff=(qg(nk1)-tma*emsd(nk1))*100./qg(nk1)
              if(abs(qvdiff).gt.1.)then
            print *,'--warning-- cloud base water vapor changes by ',    &
            qvdiff,                                                      &
             ' percent when moisture is borrowed to prevent neg values', &
             ' in kain-fritsch'
              end if
            end if
            qg(nk)=1.e-9
            qg(nk1)=tma*emsd(nk1)
            qg(nk-1)=tmb*emsd(nk-1)
          end if
  499   continue
        topomg=(udr(ltop)-uer(ltop))*dp(ltop)*emsd(ltop)
        if(abs(topomg-omg(ltop)).gt.1.e-3)then
!       write(98,*)'error:  mass does not balance in kf scheme;'
!    * ,'topomg, omg =',topomg,omg(ltop)
        write(6,*)'error:  mass does not balance in kf scheme;'        &
       ,'topomg, omg =',topomg,omg(ltop)
          istop=1
          goto 265
        end if
!
!...convert theta to t...
!
! pay attention ...
!
        do 230 nk=1,ltop
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qg(nk)))
          tg(nk)=thtag(nk)/exn(nk)
          tvg(nk)=tg(nk)*(1.+0.608*qg(nk))
  230   continue
!
!*******************************************************************
!                                                                  *
!     compute new cloud and change in available buoyant energy.    *
!                                                                  *
!*******************************************************************
!
!...the following computations are similar to that for updraft
!
        thmix=0.
        qmix=0.
        pmix=0.
        do 217 nk=lc,kpbl
          rocpq=0.2854*(1.-0.28*qg(nk))
          thmix=thmix+dp(nk)*tg(nk)*(p00/p0(nk))**rocpq
          qmix=qmix+dp(nk)*qg(nk)
  217   pmix=pmix+dp(nk)*p0(nk)
        thmix=thmix/dpthmx
        qmix=qmix/dpthmx
        pmix=pmix/dpthmx
        rocpq=0.2854*(1.-0.28*qmix)
        tmix=thmix*(pmix/p00)**rocpq
        es=aliq*exp((tmix*bliq-cliq)/(tmix-dliq))
        qs=ep2*es/(pmix-es)
!
!...remove supersaturation for diagnostic purposes, if necessary...
!
        if(qmix.gt.qs)then
          rl=xlv0-xlv1*tmix
          cpm=cp*(1.+0.887*qmix)
          dssdt=qs*(cliq-bliq*dliq)/((tmix-dliq)*(tmix-dliq))
          dq=(qmix-qs)/(1.+rl*dssdt/cpm)
          tmix=tmix+rl/cp*dq
          qmix=qmix-dq
          rocpq=0.2854*(1.-0.28*qmix)
          thmix=tmix*(p00/pmix)**rocpq
          tlcl=tmix
          plcl=pmix
        else
          qmix=amax1(qmix,0.)
          emix=qmix*pmix/(ep2+qmix)
          tlog=alog(emix/aliq)
          tdpt=(cliq-dliq*tlog)/(bliq-tlog)
          tlcl=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))*(tmix-  &
               tdpt)
          tlcl=amin1(tlcl,tmix)
          cporq=1./rocpq
          plcl=p00*(tlcl/thmix)**cporq
        end if
        tvlcl=tlcl*(1.+0.608*qmix)
        do 235 nk=lc,kl
          klcl=nk
  235   if(plcl.ge.p0(nk))goto 240
  240   k=klcl-1
        dlp=alog(plcl/p0(k))/alog(p0(klcl)/p0(k))
!
!...estimate environmental temperature and mixing ratio at the lcl...
!
        tenv=tg(k)+(tg(klcl)-tg(k))*dlp
        qenv=qg(k)+(qg(klcl)-qg(k))*dlp
        tven=tenv*(1.+0.608*qenv)
        tvbar=0.5*(tvg(k)+tven)
!        zlcl=z0(k)+r*tvbar*alog(p0(k)/plcl)/g
        zlcl=z0(k)+(z0(klcl)-z0(k))*dlp
        tvavg=0.5*(tven+tg(klcl)*(1.+0.608*qg(klcl)))
        plcl=p0(klcl)*exp(g/(r*tvavg)*(z0(klcl)-zlcl))
        theteu(k)=tmix*(1.e5/pmix)**(0.2854*(1.-0.28*qmix))*            &
                  exp((3374.6525/tlcl-2.5403)*qmix*(1.+0.81*qmix))
        es=aliq*exp((tenv*bliq-cliq)/(tenv-dliq))
        qese=ep2*es/(plcl-es)
        thtesg(k)=tenv*(1.e5/plcl)**(0.2854*(1.-0.28*qese))*            &
                  exp((3374.6525/tenv-2.5403)*qese*(1.+0.81*qese))
!
!...compute adjusted abe(abeg).
!
        abeg=0.
        thtudl=theteu(k)
        do 245 nk=k,ltopm1
          nk1=nk+1
          es=aliq*exp((tg(nk1)*bliq-cliq)/(tg(nk1)-dliq))
          qese=ep2*es/(p0(nk1)-es)
          thtesg(nk1)=tg(nk1)*(1.e5/p0(nk1))**(0.2854*(1.-0.28*qese))*    &
                      exp((3374.6525/tg(nk1)-2.5403)*qese*(1.+0.81*qese)  &
                      )
!         dzz=cvmgt(z0(klcl)-zlcl,dza(nk),nk.eq.k)
          if(nk.eq.k)then
            dzz=z0(klcl)-zlcl
          else
            dzz=dza(nk)
          end if
          be=((2.*thtudl)/(thtesg(nk1)+thtesg(nk))-1.)*dzz
  245   if(be.gt.0.)abeg=abeg+be*g
!
!...assume at least 90% of cape (abe) is removed by convection during
!...the period timec...
!
        if(noitr.eq.1)then
!        write(98,1060)fabe
          goto 265
        end if
        dabe=amax1(abe-abeg,0.1*abe)
        fabe=abeg/(abe+1.e-8)
        if(fabe.gt.1.)then
!       write(98,*)'updraft/downdraft couplet increases cape at this '
!    *,'grid point; no convection allowed!'
          goto 325
        end if
        if(ncount.ne.1)then
          dfda=(fabe-fabeold)/(ainc-aincold)
          if(dfda.gt.0.)then
            noitr=1
            ainc=aincold
            goto 255
          end if
        end if
        aincold=ainc
        fabeold=fabe
        if(ainc/aincmx.gt.0.999.and.fabe.gt.1.05-stab)then
!      write(98,1055)fabe
          goto 265
        end if
        if(fabe.le.1.05-stab.and.fabe.ge.0.95-stab)goto 265
        if(ncount.gt.10)then
!       write(98,1060)fabe
          goto 265
        end if
!
!...if more than 10% of the original cape remains, increase the
!...convective mass flux by the factor ainc:
!
        if(fabe.eq.0.)then
          ainc=ainc*0.5
        else
          ainc=ainc*stab*abe/(dabe+1.e-8)
        end if
  255   ainc=amin1(aincmx,ainc)
!...if ainc becomes very small, effects of convection
!...will be minimal so just ignore it...
        if(ainc.lt.0.05)goto 325
!       ainc=amax1(ainc,0.05)
        tder=tder2*ainc
        pptflx=pptfl2*ainc
!     write(98,1080)lfs,ldb,ldt,timec,nstep,ncount,fabeold,aincold
        do 260 nk=1,ltop
          umf(nk)=umf2(nk)*ainc
          dmf(nk)=dmf2(nk)*ainc
          detlq(nk)=detlq2(nk)*ainc
          detic(nk)=detic2(nk)*ainc
          udr(nk)=udr2(nk)*ainc
          uer(nk)=uer2(nk)*ainc
          der(nk)=der2(nk)*ainc
          ddr(nk)=ddr2(nk)*ainc
  260   continue
!
!...go back up for another iteration...
!
        goto 175
  265   continue
!
!...clean things up, calculate convective feedback tendencies for this
!...grid point...
!
!...compute hydrometeor tendencies as is done for t, qv...
!
!...frc2 is the fraction of total condensate
!...generated that goes into precipitiation
        frc2=pptflx/(cpr*ainc)
        do 270 nk=1,ltop
          qlpa(nk)=ql0(nk)
          qipa(nk)=qi0(nk)
          qrpa(nk)=qr0(nk)
          qspa(nk)=qs0(nk)
          rainfb(nk)=pptliq(nk)*ainc*fbfrc*frc2
          snowfb(nk)=pptice(nk)*ainc*fbfrc*frc2
  270   continue
        do 290 ntc=1,nstep
!
!...assign hydrometeors concentrations at the top and bottom of each
!...layer based on the sign of omega...
!
          do 275 nk=1,ltop
            qlfxin(nk)=0.
            qlfxout(nk)=0.
            qifxin(nk)=0.
            qifxout(nk)=0.
            qrfxin(nk)=0.
            qrfxout(nk)=0.
            qsfxin(nk)=0.
            qsfxout(nk)=0.
  275     continue
          do 280 nk=2,ltop
            if(omg(nk).le.0.)then
              qlfxin(nk)=-fxm(nk)*qlpa(nk-1)
              qifxin(nk)=-fxm(nk)*qipa(nk-1)
              qrfxin(nk)=-fxm(nk)*qrpa(nk-1)
              qsfxin(nk)=-fxm(nk)*qspa(nk-1)
              qlfxout(nk-1)=qlfxout(nk-1)+qlfxin(nk)
              qifxout(nk-1)=qifxout(nk-1)+qifxin(nk)
              qrfxout(nk-1)=qrfxout(nk-1)+qrfxin(nk)
              qsfxout(nk-1)=qsfxout(nk-1)+qsfxin(nk)
            else
              qlfxout(nk)=fxm(nk)*qlpa(nk)
              qifxout(nk)=fxm(nk)*qipa(nk)
              qrfxout(nk)=fxm(nk)*qrpa(nk)
              qsfxout(nk)=fxm(nk)*qspa(nk)
              qlfxin(nk-1)=qlfxin(nk-1)+qlfxout(nk)
              qifxin(nk-1)=qifxin(nk-1)+qifxout(nk)
              qrfxin(nk-1)=qrfxin(nk-1)+qrfxout(nk)
              qsfxin(nk-1)=qsfxin(nk-1)+qsfxout(nk)
            end if
  280     continue
!
!...update the hydrometeor concentration values at each level...
!
          do 285 nk=1,ltop
            qlpa(nk)=qlpa(nk)+(qlfxin(nk)+detlq(nk)-qlfxout(nk))*dtime*  &
                     emsd(nk)
            qipa(nk)=qipa(nk)+(qifxin(nk)+detic(nk)-qifxout(nk))*dtime*  &
                     emsd(nk)
            qrpa(nk)=qrpa(nk)+(qrfxin(nk)+qlqout(nk)*udr(nk)-qrfxout(nk) &
                     +rainfb(nk))*dtime*emsd(nk)
            qspa(nk)=qspa(nk)+(qsfxin(nk)+qicout(nk)*udr(nk)-qsfxout(nk) &
                     +snowfb(nk))*dtime*emsd(nk)
  285     continue
  290   continue
        do 295 nk=1,ltop
          qlg(nk)=qlpa(nk)
          qig(nk)=qipa(nk)
          qrg(nk)=qrpa(nk)
          qsg(nk)=qspa(nk)
  295   continue
!     write(98,1080)lfs,ldb,ldt,timec,nstep,ncount,fabe,ainc
!
!...send final parameterized values to output files...
!
        if(istop.eq.1)then
        write(6,1070)'  p  ','   dp ',' dt k/d ',' dr k/d ','   omg  ',  &
      ' domgdp ','   umf  ','   uer  ','   udr  ','   dmf  ','   der  '  &
      ,'   ddr  ','   ems  ','    w0  ','  detlq ',' detic '
          do 300 k=ltop,1,-1
            dtt=(tg(k)-t0(k))*86400./timec
            rl=xlv0-xlv1*tg(k)
            dr=-(qg(k)-q0(k))*rl*86400./(timec*cp)
            udfrc=udr(k)*timec*emsd(k)
            uefrc=uer(k)*timec*emsd(k)
            ddfrc=ddr(k)*timec*emsd(k)
            defrc=-der(k)*timec*emsd(k)
            write (6,1075)p0(k)/100.,dp(k)/100.,dtt,dr,omg(k),domgdp(k)* &
                          1.e4,umf(k)/1.e6,uefrc,udfrc,dmf(k)/1.e6,defrc &
                          ,ddfrc,ems(k)/1.e11,w0avg1d(k)*1.e2,detlq(k) &
                          *timec*emsd(k)*1.e3,detic(k)*timec*emsd(k)*    &
                          1.e3
  300     continue
        write(6,1085)'k','p','z','t0','tg','dt','tu','td','q0','qg',     &
                  'dq','qu','qd','qlg','qig','qrg','qsg','rh0','rhg'
          do 305 k=kx,1,-1
            dtt=tg(k)-t0(k)
            tuc=tu(k)-t00
            if(k.lt.lc.or.k.gt.ltop)tuc=0.
            tdc=tz(k)-t00
            if((k.lt.ldb.or.k.gt.ldt).and.k.ne.lfs)tdc=0.
            es=aliq*exp((bliq*tg(k)-cliq)/(tg(k)-dliq))
            qgs=es*ep2/(p0(k)-es)
            rh0=q0(k)/qes(k)
            rhg=qg(k)/qgs
            write (6,1090)k,p0(k)/100.,z0(k),t0(k)-t00,tg(k)-t00,dtt,tuc &
                          ,tdc,q0(k)*1000.,qg(k)*1000.,(qg(k)-q0(k))*    &
                          1000.,qu(k)*1000.,qd(k)*1000.,qlg(k)*1000.,    &
                          qig(k)*1000.,qrg(k)*1000.,qsg(k)*1000.,rh0,rhg
  305     continue
!
!...if calculations above show an error in the mass budget, print out a
!...to be used later for diagnostic purposes, then abort run...
!
          if(istop.eq.1)then
            do 310 k=1,kx
              write ( wrf_err_message , 1115 )                         &
                            z0(k),p0(k)/100.,t0(k)-273.16,q0(k)*1000.,   &
                            u0(k),v0(k),dp(k)/100.,w0avg1d(k)
              call wrf_message ( trim( wrf_err_message ) )
  310       continue
            call wrf_error_fatal ( 'module_cu_kf.f: kain-fritsch' )
          end if
        end if
        cndtnf=(1.-eqfrc(lfs))*(qliq(lfs)+qice(lfs))*dmf(lfs)
!     write(98,1095)cpr*ainc,tder+pptflx+cndtnf
!
!  evaluate moisture budget...
!
        qinit=0.
        qfnl=0.
        dpt=0.
        do 315 nk=1,ltop
          dpt=dpt+dp(nk)
          qinit=qinit+q0(nk)*ems(nk)
          qfnl=qfnl+qg(nk)*ems(nk)
          qfnl=qfnl+(qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
  315   continue
        qfnl=qfnl+pptflx*timec*(1.-fbfrc)
        err2=(qfnl-qinit)*100./qinit
!     write(98,1110)qinit,qfnl,err2
!        if(abs(err2).gt.0.05)stop 'qverr'
        if(abs(err2).gt.0.05)call wrf_error_fatal( 'module_cu_kf.f: qverr' )
        relerr=err2*qinit/(pptflx*timec+1.e-10)
!     write(98,1120)relerr
!     write(98,*)'tder, cpr, usr, trppt =',
!    *tder,cpr*ainc,usr*ainc,trppt*ainc
!
!...feedback to resolvable scale tendencies.
!
!...if the advective time period (tadvec) is less than specified minimum
!...timec, allow feedback to occur only during tadvec...
!
        if(tadvec.lt.timec)nic=nint(tadvec/dt)
        nca(i,j)=float(nic)*dt
        do 320 k=1,kx
!         if(imoist.ne.2)then
!
!...if hydrometeors are not allowed, they must be evaporated or sublimat
!...and fed back as vapor, along with associated changes in temperature.
!...note:  this will introduce changes in the convective temperature and
!...water vapor feedback tendencies and may lead to supersaturated value
!...of qg...
!
!           rlc=xlv0-xlv1*tg(k)
!           rls=xls0-xls1*tg(k)
!           cpm=cp*(1.+0.887*qg(k))
!           tg(k)=tg(k)-(rlc*(qlg(k)+qrg(k))+rls*(qig(k)+qsg(k)))/cpm
!           qg(k)=qg(k)+(qlg(k)+qrg(k)+qig(k)+qsg(k))
!           dqcdt(k)=0.
!           dqidt(k)=0.
!           dqrdt(k)=0.
!           dqsdt(k)=0.
!         else
            if(.not. qi_flag .and. warm_rain)then
!
!...if ice phase is not allowed, melt all frozen hydrometeors...
!
              cpm=cp*(1.+0.887*qg(k))
              tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
              dqcdt(k)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
              dqidt(k)=0.
              dqrdt(k)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
              dqsdt(k)=0.
            elseif(.not. qi_flag .and. .not. warm_rain)then
!
!...if ice phase is allowed, but mixed phase is not, melt frozen hydrome
!...below the melting level, freeze liquid water above the melting level
!
              cpm=cp*(1.+0.887*qg(k))
              if(k.le.ml)then
                tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
              elseif(k.gt.ml)then
                tg(k)=tg(k)+(qlg(k)+qrg(k))*rlf/cpm
              end if
              dqcdt(k)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
              dqidt(k)=0.
              dqrdt(k)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
              dqsdt(k)=0.
            elseif(qi_flag) then
!
!...if mixed phase hydrometeors are allowed, feed back convective
!...tendency of hydrometeors directly...
!
              dqcdt(k)=(qlg(k)-ql0(k))/timec
              dqidt(k)=(qig(k)-qi0(k))/timec
              dqrdt(k)=(qrg(k)-qr0(k))/timec
              if (qs_flag ) then
                 dqsdt(k)=(qsg(k)-qs0(k))/timec
              else
                 dqidt(k)=dqidt(k)+(qsg(k)-qs0(k))/timec
              end if
            else
              call wrf_error_fatal ( 'module_cu_kf: this combination of imoist,  iice not allowed' )
            end if
!         end if
          dtdt(k)=(tg(k)-t0(k))/timec
          dqdt(k)=(qg(k)-q0(k))/timec
  320   continue

! raincv is in the unit of mm

        pratec(i,j)=pptflx*(1.-fbfrc)/dxsq
        raincv(i,j)=dt*pratec(i,j)
        rnc=raincv(i,j)*nic
!        write(98,909)rnc
 909     format(' convective rainfall =',f8.4,' cm')

  325 continue

1000  format(' ',10a8)
1005  format(' ',f6.0,2x,f6.4,2x,f7.3,1x,f6.4,2x,4(f6.3,2x),2(f7.3,1x))
1010  format(' ',' vertical velocity is negative at ',f4.0,' mb')
1015   format(' ','all remaining mass detrains below ',f4.0,' mb')
1025   format(5x,' klcl=',i2,' zlcl=',f7.1,'m',                          &
        ' dtlcl=',f5.2,' ltop=',i2,' p0(ltop)=',-2pf5.1,'mb frz lv=',    &
        i2,' tmix=',0pf4.1,1x,'pmix=',-2pf6.1,' qmix=',3pf5.1,           &
        ' cape=',0pf7.1)
1030   format(' ',' p0(let) = ',f6.1,' p0(ltop) = ',f6.1,' vmflcl =',    &
      e12.3,' plcl =',f6.1,' wlcl =',f6.3,' cldhgt =',                   &
      f8.1)
1035  format(1x,'pef(ws)=',f4.2,'(cb)=',f4.2,'lc,let=',2i3,'wkl='        &
      ,f6.3,'vws=',f5.2)
1040          format(' ','precip eff = 100%, envir cannot support downd' &
      ,'rafts')
!1045  format('number of downdraft iterations exceeds 10...pptflx'       &
!      ' is different from that given by precip eff relation')
! flic has trouble with this one.
1045  format('number of downdraft iterations exceeds 10')
1050  format(' ','lcount= ',i3,' pptflx/cpr, peff= ',f5.3,1x,f5.3,       &
      'dmf(lfs)/umf(lcl)= ',f5.3)
1055     format(/'*** degree of stabilization =',f5.3,', no more mass f' &
      ,'lux is allowed')
!1060     format(/' iteration does not converge to give the specified '  &
!      'degree of stabilization!  fabe= ',f6.4)
1060  format(/' iteration does not converge.  fabe= ',f6.4)
 1070 format (16a8)
 1075 format (f8.2,3(f8.2),2(f8.3),f8.2,2f8.3,f8.2,6f8.3)
1080   format(2x,'lfs,ldb,ldt =',3i3,' timec, nstep=',f5.0,i3,           &
      'ncount, fabe, ainc=',i2,1x,f5.3,f6.2)
 1085 format (a3,16a7,2a8)
 1090 format (i3,f7.2,f7.0,10f7.2,4f7.3,2f8.3)
1095   format(' ','  ppt production rate= ',f10.0,' total evap+ppt= ',   &
       f10.0)
1105   format(' ','net latent heat release =',e12.5,' actual heating =', &
       e12.5,' j/kg-s, difference = ',f9.3,'percent')
1110   format(' ','initial water =',e12.5,' final water =',e12.5,        &
       ' total water change =',f8.2,'percent')
 1115 format (2x,f6.0,2x,f7.2,2x,f5.1,2x,f6.3,2(2x,f5.1),2x,f7.2,2x,f7.4 &
             )
1120   format(' ','moisture error as function of total ppt =',f9.3,      &
            'percent')

   end subroutine kfpara

!-----------------------------------------------------------------------
   subroutine condload(qliq,qice,wtw,dz,boterm,enterm,rate,qnewlq,     &
                       qnewic,qlqout,qicout,g)
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!  9/18/88...this precipitation fallout scheme is based on the scheme us
!  by ogura and cho (1973).  liquid water fallout from a parcel is cal-
!  culated using the equation dq=-rate*q*dt, but to simulate a quasi-
!  continuous process, and to eliminate a dependency on vertical
!  resolution this is expressed as q=q*exp(-rate*dz).

      real, intent(in   )   :: g
      real, intent(in   )   :: dz,boterm,enterm,rate
      real, intent(inout)   :: qlqout,qicout,wtw,qliq,qice,qnewlq,qnewic

      real :: qtot,qnew,qest,g1,wavg,conv,ratio3,oldq,ratio4,dq,pptdrg

      qtot=qliq+qice
      qnew=qnewlq+qnewic
!
!  estimate the vertical velocity so that an average vertical velocity c
!  be calculated to estimate the time required for ascent between model
!  levels...
!
      qest=0.5*(qtot+qnew)
      g1=wtw+boterm-enterm-2.*g*dz*qest/1.5
      if(g1.lt.0.0)g1=0.
      wavg=(sqrt(wtw)+sqrt(g1))/2.
      conv=rate*dz/wavg
!
!  ratio3 is the fraction of liquid water in fresh condensate, ratio4 is
!  the fraction of liquid water in the total amount of condensate involv
!  in the precipitation process - note that only 60% of the fresh conden
!  sate is is allowed to participate in the conversion process...
!
      ratio3=qnewlq/(qnew+1.e-10)
!     oldq=qtot
      qtot=qtot+0.6*qnew
      oldq=qtot
      ratio4=(0.6*qnewlq+qliq)/(qtot+1.e-10)
      qtot=qtot*exp(-conv)
!
!  determine the amount of precipitation that falls out of the updraft
!  parcel at this level...
!
      dq=oldq-qtot
      qlqout=ratio4*dq
      qicout=(1.-ratio4)*dq
!
!  estimate the mean load of condensate on the updraft in the layer, cal
!  late vertical velocity
!
      pptdrg=0.5*(oldq+qtot-0.2*qnew)
      wtw=wtw+boterm-enterm-2.*g*dz*pptdrg/1.5
!
!  determine the new liquid water and ice concentrations including losse
!  due to precipitation and gains from condensation...
!
      qliq=ratio4*qtot+ratio3*0.4*qnew
      qice=(1.-ratio4)*qtot+(1.-ratio3)*0.4*qnew
      qnewlq=0.
      qnewic=0.

   end subroutine condload

!-----------------------------------------------------------------------
   subroutine dtfrznew(tu,p,thteu,qvap,qliq,qice,ratio2,ttfrz,tbfrz,   &
                       qnwfrz,rl,frc1,effq,iflag,xlv0,xlv1,xls0,xls1,  &
                       ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice     )
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   real,         intent(in   )   :: xlv0,xlv1
   real,         intent(in   )   :: p,ttfrz,tbfrz,effq,xls0,xls1,ep2,aliq, &
                                    bliq,cliq,dliq,aice,bice,cice,dice
   real,         intent(inout)   :: tu,thteu,qvap,qliq,qice,ratio2,    &
                                    frc1,rl,qnwfrz
   integer,      intent(inout)   :: iflag

   real ::       ccp,rv,c5,qlqfrz,qnew,esliq,esice,rlc,rls,pi,es,rlf,a, &
                 b,c,dqvap,dtfrz,tu1,qvap1
!-----------------------------------------------------------------------
!
!...allow glaciation of the updraft to occur as an approximately linear
!   function of temperature in the temperature range ttfrz to tbfrz...
!

      rv=461.5
      c5=1.0723e-3
!
!...adjust the liquid water concentrations from fresh condensate and tha
!   brought up from lower levels to an amount that would be present if n
!   liquid water had frozen thus far...this is necessary because the
!   expression for temp change is multiplied by the fraction equal to th
!   parcel temp decrease since the last model level divided by the total
!   glaciation interval, so that effectively this approximately allows a
!   amount of liquid water to freeze which is equal to this same fractio
!   of the liquid water that was present before the glaciation process w
!   initiated...also, to allow thetau to convert approximately linearly
!   its value with respect to ice, we need to allow a portion of the fre
!   condensate to contribute to the glaciation process; the fractional
!   amount that applies to this portion is 1/2 of the fractional amount
!   frozen of the "old" condensate because this fresh condensate is only
!   produced gradually over the layer...note that in terms of the dynami
!   of the precipitation process, ie. precipitation fallout, this fracti
!   amnt of fresh condensate has already been included in the ice catego
!
      qlqfrz=qliq*effq
      qnew=qnwfrz*effq*0.5
      esliq=aliq*exp((bliq*tu-cliq)/(tu-dliq))
      esice=aice*exp((bice*tu-cice)/(tu-dice))
      rlc=2.5e6-2369.276*(tu-273.16)
      rls=2833922.-259.532*(tu-273.16)
      rlf=rls-rlc
      ccp=1005.7*(1.+0.89*qvap)
!
!  a = d(es)/dt is that calculated from buck`s (1981) empirical formulas
!  for saturation vapor pressure...
!
      a=(cice-bice*dice)/((tu-dice)*(tu-dice))
      b=rls*ep2/p
      c=a*b*esice/ccp
      dqvap=b*(esliq-esice)/(rls+rls*c)-rlf*(qlqfrz+qnew)/(rls+rls/c)
      dtfrz=(rlf*(qlqfrz+qnew)+b*(esliq-esice))/(ccp+a*b*esice)
      tu1=tu
      qvap1=qvap
      tu=tu+frc1*dtfrz
      qvap=qvap-frc1*dqvap
      es=qvap*p/(ep2+qvap)
      esliq=aliq*exp((bliq*tu-cliq)/(tu-dliq))
      esice=aice*exp((bice*tu-cice)/(tu-dice))
      ratio2=(esliq-es)/(esliq-esice)
!
!  typically, ratio2 is very close to (ttfrz-tu)/(ttfrz-tbfrz), usually
!  within 1% (using tu before galciation effects are applied);  if the
!  initial updraft temp is below tbfrz and ratio2 is still less than 1,
!  an adjustment to frc1 and ratio2 is introduced so that glaciation
!  effects are not underestimated; conversely, if ratio2 is greater than
!  frc1 is adjusted so that glaciation effects are not overestimated...
!
      if(iflag.gt.0.and.ratio2.lt.1)then
        frc1=frc1+(1.-ratio2)
        tu=tu1+frc1*dtfrz
        qvap=qvap1-frc1*dqvap
        ratio2=1.
        iflag=1
        goto 20
      end if
      if(ratio2.gt.1.)then
        frc1=frc1-(ratio2-1.)
        frc1=amax1(0.0,frc1)
        tu=tu1+frc1*dtfrz
        qvap=qvap1-frc1*dqvap
        ratio2=1.
        iflag=1
      end if
!
!  calculate a hybrid value of thetau, assuming that the latent heat of
!  vaporization/sublimation can be estimated using the same weighting
!  function as that used to calculate saturation vapor pressure, calcu-
!  late new liquid water and ice concentrations...
!
   20 rlc=xlv0-xlv1*tu
      rls=xls0-xls1*tu
      rl=ratio2*rls+(1.-ratio2)*rlc
      pi=(1.e5/p)**(0.2854*(1.-0.28*qvap))
      thteu=tu*pi*exp(rl*qvap*c5/tu*(1.+0.81*qvap))
      if(iflag.eq.1)then
        qice=qice+frc1*dqvap+qliq
        qliq=0.
      else
        qice=qice+frc1*(dqvap+qlqfrz)
        qliq=qliq-frc1*qlqfrz
      end if
      qnwfrz=0.

   end subroutine dtfrznew

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  this subroutine integrates the area under the curve in the gaussian
!  distribution...the numerical approximation to the integral is taken f
!   handbook of mathematical functions with formulas, graphs and mathema
!  tables  ed. by abramowitz and stegun, nat l bureau of standards appli
!  mathematics series.  june, 1964., may, 1968.
!                                     jack kain
!                                     7/6/89
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!***********************************************************************
!*****    gaussian type mixing profile....******************************
   subroutine prof5(eq,ee,ud)
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   real,         intent(in   )   :: eq
   real,         intent(inout)   :: ee,ud
   real ::       sqrt2p,a1,a2,a3,p,sigma,fe,x,y,ey,e45,t1,t2,c1,c2

      data sqrt2p,a1,a2,a3,p,sigma,fe/2.506628,0.4361836,-0.1201676,    &
      0.9372980,0.33267,0.166666667,0.202765151/
      x=(eq-0.5)/sigma
      y=6.*eq-3.
      ey=exp(y*y/(-2))
      e45=exp(-4.5)
      t2=1./(1.+p*abs(y))
      t1=0.500498
      c1=a1*t1+a2*t1*t1+a3*t1*t1*t1
      c2=a1*t2+a2*t2*t2+a3*t2*t2*t2
      if(y.ge.0.)then
        ee=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*eq*eq/2.
        ud=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*(0.5+eq*eq/2.-  &
           eq)
      else
        ee=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*eq*eq/2.
        ud=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*(0.5+eq* &
           eq/2.-eq)
      end if
      ee=ee/fe
      ud=ud/fe

  end subroutine prof5

#endif
  subroutine tpmix(p,thtu,tu,qu,qliq,qice,qnewlq,qnewic,ratio2,rl,    &
                   xlv0,xlv1,xls0,xls1,                               &
                   ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice        )
    implicit none
    real(rk8) , intent(in) :: xlv0 , xlv1
    real(rk8) , intent(in) :: p , thtu , ratio2 , rl , xls0 , xls1
    real(rk8) , intent(in) :: ep2 , aliq , bliq , cliq , dliq
    real(rk8) , intent(in) :: aice , bice , cice , dice
    real(rk8) , intent(inout) :: qu , qliq , qice , tu , qnewlq , qnewic
    real(rk8) :: es , qs , pi , thtgs , f0 , t1 , t0
    real(rk8) :: esliq , esice , f1 , dt , qnew
    real(rk8) :: dq , qtot , dqice , dqliq , rll , ccp
    real(rk8) , parameter :: c5 = 1.0723D-3
    real(rk8) , parameter :: rv = 461.5D0
    integer(ik4) :: itcnt
    !
    ! This subroutine iteratively extracts wet-bulb temperature from equiv
    ! potential temperature, then checks to see if sufficient moisture is
    ! available to achieve saturation...if not, temperature is adjusted
    ! accordingly, if so, the residual liquid water/ice concentration is
    ! determined...
    !
    ! Iterate to find wet bulb temperature as a function of equivalent pot
    ! temp and prs, assuming saturation vapor pressure...ratio2 is the deg
    ! of glaciation...
    !
    if ( ratio2 < 1.D-6 ) then
      es = aliq*exp((bliq*tu - cliq) / (tu-dliq))
      qs = ep2*es/(p-es)
      pi = (1.D5/p)**(0.2854D0*(d_one - 0.28D0*qs))
      thtgs = tu * pi * exp((3374.6525D0/tu - 2.5403D0)*qs*(d_one + 0.81D0*qs))
    else if ( abs(ratio2-d_one) < 1.D-6 ) then
      es = aice * exp((bice*tu - cice)/(tu-dice))
      qs = ep2 * es / (p-es)
      pi = (1.D5/p)**(0.2854D0*(d_one - 0.28D0*qs))
      thtgs = tu * pi * exp((3114.834D0/tu - 0.278296D0)*qs*(d_one+0.81D0*qs))
    else
      esliq = aliq*exp((bliq*tu - cliq) / (tu-dliq))
      esice = aice*exp((bice*tu - cice) / (tu-dice))
      es = (d_one-ratio2)*esliq + ratio2*esice
      qs = ep2*es / (p-es)
      pi = (1.D5/p)**(0.2854D0*(d_one-0.28D0*qs))
      thtgs = tu * pi * exp(rl * qs * c5/tu * (d_one+0.81D0*qs))
    end if
    f0 = thtgs - thtu
    t1 = tu-0.5D0*f0
    t0 = tu
    itcnt = 0
    do
      if ( ratio2 < 1.D-6 ) then
        es = aliq * exp((bliq*t1 - cliq) / (t1-dliq))
        qs = ep2*es/(p-es)
        pi = (1.D5/p)**(0.2854D0*(d_one-0.28D0*qs))
        thtgs = t1 * pi * exp((3374.6525D0/t1 - 2.5403D0)*qs*(d_one+0.81D0*qs))
      else if ( abs(ratio2-d_one) < 1.D-6 ) then
        es = aice * exp((bice*t1 - cice) / (t1-dice))
        qs = ep2*es / (p-es)
        pi = (1.D5/p)**(0.2854D0*(d_one-0.28D0*qs))
        thtgs = t1 * pi * exp((3114.834D0/t1-0.278296D0)*qs*(d_one+0.81D0*qs))
      else
        esliq = aliq * exp((bliq*t1 - cliq) / (t1-dliq))
        esice = aice * exp((bice*t1 - cice) / (t1-dice))
        es = (d_one-ratio2) * esliq + ratio2*esice
        qs = ep2 * es / (p-es)
        pi = (1.D5/p)**(0.2854D0*(d_one-0.28D0*qs))
        thtgs = t1 * pi * exp(rl * qs * c5/t1 * (d_one+0.81D0*qs))
      end if
      f1 = thtgs - thtu
      if ( abs(f1) < 0.01D0 ) exit
      itcnt = itcnt + 1
      if ( itcnt > 10 ) exit
      dt = f1 * (t1-t0) / (f1-f0)
      t0 = t1
      f0 = f1
      t1 = t1 - dt
    end do
    !
    !If the parcel is supersaturated, calculate concentration of fresh
    !condensate...
    !
    do
      if ( qs <= qu ) then
        qnew = qu - qs
        qu = qs
        exit
      end if
      !
      ! if the parcel is subsaturated, temperature and mixing ratio must be
      ! adjusted...if liquid water or ice is present, it is allowed to evapo
      ! sublimate.
      !
      qnew = d_zero
      dq = qs - qu
      qtot = qliq + qice
      !
      ! If there is enough liquid or ice to saturate the parcel, temp stays
      ! wet bulb value, vapor mixing ratio is at saturated level, and the mi
      ! ratios of liquid and ice are adjusted to make up the original satura
      ! deficit... otherwise, any available liq or ice vaporizes and appropr
      ! adjustments to parcel temp; vapor, liquid, and ice mixing ratios are
      !
      ! Note that the liq and ice may be present in proportions slightly dif
      ! than suggested by the value of ratio2...check to make sure that liq
      ! ice concentrations are not reduced to below zero when evaporation/
      ! sublimation occurs...
      !
      if ( qtot >= dq ) then
        dqice = d_zero
        dqliq = d_zero
        qliq = qliq - (d_one - ratio2)*dq
        if ( qliq < d_zero ) then
          dqice = d_zero-qliq
          qliq = d_zero
        end if
        qice = qice - ratio2*dq + dqice
        if ( qice < d_zero ) then
          dqliq = d_zero - qice
          qice = d_zero
        end if
        qliq = qliq + dqliq
        qu = qs
        exit
      else
        if ( ratio2 < 1.D-6 ) then
          rll = xlv0 - xlv1*t1
        else if( abs(ratio2-d_one) < 1.D-6 ) then
          rll = xls0 - xls1*t1
        else
          rll=rl
        end if
        ccp = cpd * (d_one + 0.89D0*qu)
        if ( qtot < 1.D-10 ) then
          !
          ! If no liquid water or ice is available, temperature is given by:
          t1 = t1 + rll * (dq/(d_one+dq))/ccp
          exit
        else
          !
          ! If some liq water/ice is available, but not enough to achieve satura
          ! the temperature is given by:
          t1 = t1 + rll * ((dq-qtot)/(d_one + dq - qtot))/ccp
          qu = qu + qtot
          qtot = d_zero
        end if
        qliq = d_zero
        qice = d_zero
      end if
      exit
    end do
    tu = t1
    qnewlq = (d_one - ratio2)*qnew
    qnewic = ratio2 * qnew
    if ( itcnt > 10) then
      write(stderr,*) '***** number of iterations in tpmix =', itcnt
    end if
  end subroutine tpmix

  subroutine envirtht(p1,t1,q1,tht1,r1,rl,ep2,aliq,bliq,cliq,dliq, &
                      aice,bice,cice,dice)
    implicit none
    real(rk8) , intent(in) :: p1 , t1 , q1 , r1 , rl , ep2
    real(rk8) , intent(in) :: aliq , bliq , cliq , dliq
    real(rk8) , intent(in) :: aice , bice , cice , dice
    real(rk8) ,  intent(inout) :: tht1
    real(rk8) , parameter :: p00 = 1.0D5
    real(rk8) , parameter :: c1 = 3374.6525
    real(rk8) , parameter :: c2 = 2.5403
    real(rk8) , parameter :: c3 = 3114.834
    real(rk8) , parameter :: c4 = 0.278296
    real(rk8) , parameter :: c5 = 1.0723D-3

    real(rk8) :: ee , tlog , tdpt , tsat , tht , tfpt , tlogic
    real(rk8) :: tsatlq , tsatic

    !  calculate environmental equivalent potential temperature...

    if ( r1 < 1.D-6 ) then
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aliq)
      tdpt = (cliq - dliq*tlog) / (bliq-tlog)
      tsat = tdpt - (0.212D0 + 1.571D-3*(tdpt-tzero) - &
                               4.360D-4*(t1-tzero)) * (t1-tdpt)
      tht = t1 * (p00/p1)**(0.2854D0*(d_one - 0.28D0*q1))
      tht1 = tht * exp((c1/tsat-c2)*q1*(d_one + 0.81D0*q1))
    else if ( abs(r1-d_one) < 1.D-6 ) then
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aice)
      tfpt = (cice - dice*tlog) / (bice - tlog)
      tht = t1*(p00/p1)**(0.2854D0*(d_one - 0.28D0*q1))
      tsat = tfpt - (0.182D0 + 1.13D-3*(tfpt-tzero) - &
                               3.58D-4*(t1-tzero)) * (t1-tfpt)
      tht1 = tht * exp((c3/tsat-c4)*q1*(d_one + 0.81D0*q1))
    else
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aliq)
      tdpt = (cliq-dliq*tlog) / (bliq-tlog)
      tlogic = log(ee/aice)
      tfpt = (cice - dice*tlogic) / (bice-tlogic)
      tht = t1 * (p00/p1)**(0.2854D0*(d_one - 0.28D0*q1))
      tsatlq = tdpt - (0.212D0 + 1.571D-3*(tdpt-tzero) - &
                                 4.36D-4*(t1-tzero)) * (t1-tdpt)
      tsatic = tfpt - (0.182D0 + 1.13D-3*(tfpt-tzero) - &
                                 3.58D-4*(t1-tzero)) * (t1-tfpt)
      tsat = r1 * tsatic + (d_one-r1)*tsatlq
      tht1 = tht * exp(rl*q1*c5/tsat*(d_one+0.81D0*q1))
    end if
  end subroutine envirtht

  ! This subroutine iteratively extracts temperature from equivalent
  ! potential temp.  it is designed for use with downdraft calculations.
  ! if relative humidity is specified to be less than 100%, parcel
  ! temp, specific humidity, and liquid water content are iteratively
  ! calculated.
  real(rk8) function tpdd(p,thted,tgs,rs,rd,rh,xlv0,xlv1, &
                          ep2,aliq,bliq,cliq,dliq,aice,bice,cice,dice)
    implicit none
    real(rk8) , intent(in) :: xlv0 , xlv1
    real(rk8) , intent(in) :: p , thted , tgs , rd , rh , ep2
    real(rk8) , intent(in) :: aliq , bliq , cliq , dliq
    real(rk8) , intent(in) :: aice , bice , cice , dice
    real(rk8) , intent(inout) :: rs
    real(rk8) :: es , pi , thtgs , f0 , t1 , t0 , ccp , f1 , dt
    real(rk8) :: rl , dssdt , t1rh , rsrh
    integer(ik4) :: itcnt

    es = aliq * exp( (bliq*tgs - cliq) / (tgs - dliq) )
    rs = ep2 * es / (p - es)
    pi = (1.D5/p)**(0.2854D0*(1.0D0-0.28D0*rs))
    thtgs = tgs * pi * exp((3374.6525D0/tgs - 2.5403D0) * rs*(1.0D0+0.81D0*rs))
    f0 = thtgs - thted
    t1 = tgs - 0.5D0*f0
    t0 = tgs
    ccp = cpd
    !
    ! Iterate to find wet-bulb temperature...
    !
    itcnt = 0
    iterate_loop: &
    do
      es = aliq * exp((bliq*t1 - cliq)/(t1 - dliq))
      rs = ep2*es/(p-es)
      pi = (1.D5/p)**(0.2854D0*(1.0D0-0.28D0*rs))
      thtgs = t1 * pi * exp((3374.6525D0/t1-2.5403D0)*rs*(1.0D0+0.81D0*rs))
      f1 = thtgs-thted
      if ( abs(f1) < 0.05D0 ) exit iterate_loop
      itcnt = itcnt + 1
      if (itcnt > 10 ) exit iterate_loop
      dt = f1 * (t1-t0) / (f1-f0)
      t0 = t1
      f0 = f1
      t1 = t1-dt
    end do iterate_loop
    rl = xlv0 - xlv1*t1
    !
    ! If relative humidity is specified to be less than 100%, estimate the
    ! temperature and mixing ratio which will yield the appropriate value.
    !
    if ( (rh - 1.0D0) < dlowval ) then
      tpdd = t1
    else
      dssdt = (cliq-bliq * dliq) / ((t1-dliq)*(t1-dliq))
      dt = rl*rs*(1.0D0-rh) / (ccp+rl*rh*rs*dssdt)
      t1rh = t1+dt
      es = rh*aliq*exp((bliq*t1rh - cliq)/(t1rh-dliq))
      rsrh = ep2*es/(p-es)
      !
      ! Check to see if mixing ratio at specified rh is less than actual
      ! mixing ratio...if so, adjust to give zero evaporation...
      !
      if ( rsrh < rd ) then
        rsrh = rd
        t1rh = t1 + (rs-rsrh)*rl/ccp
      end if
      t1 = t1rh
      rs = rsrh
      tpdd = t1
    end if
    if ( itcnt < 10 ) then
      write(stderr,*) '***** number of iterations in tpdd = ', itcnt
    end if
  end function tpdd

  subroutine kfinit( rthcuten,rqvcuten,rqccuten,rqrcuten,           &
                     rqicuten,rqscuten,nca,w0avg,p_qi,p_qs,         &
                     p_first_scalar,restart,                        &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte )
    implicit none
    logical , intent(in) ::  restart
    integer(ik4) , intent(in) ::  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte
    integer(ik4) , intent(in) ::  p_qi , p_qs , p_first_scalar

    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rthcuten
    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rqvcuten
    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rqccuten
    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rqrcuten
    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rqicuten
    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: rqscuten

    real(rk8) , dimension(ims:ime,kms:kme,jms:jme) , intent(out) :: w0avg

    real(rk8) , dimension(ims:ime,jms:jme), intent(inout) :: nca

    integer(ik4) :: i , j , k , itf , jtf , ktf

    jtf = min(jte,jde-1)
    ktf = min(kte,kde-1)
    itf = min(ite,ide-1)

    if ( .not. restart ) then
      do j = jts , jtf
        do k = kts , ktf
          do i = its , itf
            rthcuten(i,k,j) = d_zero
            rqvcuten(i,k,j) = d_zero
            rqccuten(i,k,j) = d_zero
            rqrcuten(i,k,j) = d_zero
          end do
        end do
      end do
      if ( p_qi >= p_first_scalar ) then
        do j = jts , jtf
          do k = kts , ktf
            do i = its , itf
              rqicuten(i,k,j) = d_zero
            end do
          end do
        end do
      end if
      if ( p_qs >= p_first_scalar ) then
        do j = jts , jtf
          do k = kts , ktf
            do i = its , itf
              rqscuten(i,k,j) = d_zero
            end do
          end do
        end do
      end if
      do j = jts , jtf
        do i = its , itf
          nca(i,j) = -d_100
        end do
      end do
      do j = jts , jtf
        do k = kts , ktf
          do i = its , itf
            w0avg(i,k,j) = d_zero
          end do
        end do
      end do
    end if
  end subroutine kfinit

end module module_cu_kf
