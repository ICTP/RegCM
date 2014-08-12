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
module mod_cu_kf

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_memutil
  use mod_stdio
  use mod_regcm_types
  use mod_mpmessage
  use mod_cu_common
  use mod_runparams , only : dx , dxsq , ipptls , dtsec
  use mod_runparams , only : iqv , iqr , iqi , iqs , iqc
  use mod_service

  implicit none

  private

  public :: allocate_mod_cu_kf , kfdrv

  real(rk8) , parameter :: rad = 1500.0D0
  real(rk8) , parameter :: rv = 461.5D0
  real(rk8) , parameter :: p00 = 1.0D5
  real(rk8) , parameter :: c1 = 3374.6525D0
  real(rk8) , parameter :: c2 = 2.5403D0
  real(rk8) , parameter :: c3 = 3114.834D0
  real(rk8) , parameter :: c4 = 0.278296D0
  real(rk8) , parameter :: c5 = 1.0723D-3
  real(rk8) , parameter :: c6 = 0.2854D0
  real(rk8) , parameter :: c7 = 0.28D0
  real(rk8) , parameter :: c8 = 0.81D0
  real(rk8) , parameter :: c9 = 0.89D0

  real(rk8) , parameter :: tl1 = 0.212D0
  real(rk8) , parameter :: tl2 = 1.571D-3
  real(rk8) , parameter :: tl3 = 4.360D-4
  real(rk8) , parameter :: ti1 = 0.182D0
  real(rk8) , parameter :: ti2 = 1.13D-3
  real(rk8) , parameter :: ti3 = 3.58D-4

  real(rk8) , parameter :: cv = 717.0D0
  real(rk8) , parameter :: b61 = 0.608D0
  real(rk8) , parameter :: rlf = 3.339D5
  real(rk8) , parameter :: rhic = 1.0D0
  real(rk8) , parameter :: rhbc = 0.90D0
  real(rk8) , parameter :: ttfrz = 268.16D0
  real(rk8) , parameter :: tbfrz = 248.16D0
  real(rk8) , parameter :: rate = 0.01D0
  !
  ! Option to feed convectively generated rainwater into grid-resolved 
  ! rainwater (or snow/graupel) field.  'fbfrc' is the fraction of
  ! available precipitation to be fed back (0.0 - 1.0)
  !
  real(rk8) , parameter :: fbfrc = d_zero

  integer(ik4) , parameter :: maxiter = 10

  integer(ik4) :: nipoi

  integer(ik4) , dimension(:) , pointer :: imap , jmap
  real(rk8) , dimension(:,:) , pointer :: u0 , v0 , z0 , t0 , qv0 , p0
  real(rk8) , dimension(:,:) , pointer :: rho , dzq , w0avg
  real(rk8) , dimension(:) , pointer :: raincv , pratec
  integer(ik4) , dimension(:) , pointer :: ktop , kbot

  real(rk8) , dimension(:,:) , pointer :: dqdt , dtdt , dqcdt

  ! IPPTLS == 2
  real(rk8) , dimension(:,:) , pointer :: dqidt , dqrdt , dqsdt

  contains

  subroutine allocate_mod_cu_kf
    implicit none
    integer(ik4) :: ii , i , j
    nipoi = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( cuscheme(j,i) == 6 ) then
          nipoi = nipoi + 1
        end if
      end do
    end do
    if ( nipoi == 0 ) return
    call getmem1d(imap,1,nipoi,'mod_cu_kf:imap')
    call getmem1d(jmap,1,nipoi,'mod_cu_kf:jmap')
    ii = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( cuscheme(j,i) == 6 ) then
          imap(ii) = i
          jmap(ii) = j
          ii = ii + 1
        end if
      end do
    end do
    call getmem2d(u0,1,nipoi,1,kz,'mod_cu_kf:u0')
    call getmem2d(v0,1,nipoi,1,kz,'mod_cu_kf:v0')
    call getmem2d(t0,1,nipoi,1,kz,'mod_cu_kf:t0')
    call getmem2d(z0,1,nipoi,1,kz,'mod_cu_kf:t0')
    call getmem2d(qv0,1,nipoi,1,kz,'mod_cu_kf:qv0')
    call getmem2d(p0,1,nipoi,1,kz,'mod_cu_kf:p0')
    call getmem2d(rho,1,nipoi,1,kz,'mod_cu_kf:rho')
    call getmem2d(dzq,1,nipoi,1,kz,'mod_cu_kf:dzq')
    call getmem2d(w0avg,1,nipoi,1,kz,'mod_cu_kf:w0avg')
    call getmem2d(dqdt,1,nipoi,1,kz,'mod_cu_kf:dqdt')
    call getmem2d(dqidt,1,nipoi,1,kz,'mod_cu_kf:dqidt')
    call getmem2d(dqcdt,1,nipoi,1,kz,'mod_cu_kf:dqcdt')
    call getmem2d(dqrdt,1,nipoi,1,kz,'mod_cu_kf:dqrdt')
    call getmem2d(dqsdt,1,nipoi,1,kz,'mod_cu_kf:dqsdt')
    call getmem2d(dtdt,1,nipoi,1,kz,'mod_cu_kf:dtdt')
    call getmem1d(raincv,1,nipoi,'mod_cu_kf:raincv')
    call getmem1d(pratec,1,nipoi,'mod_cu_kf:pratec')
    call getmem1d(ktop,1,nipoi,'mod_cu_kf:ktop')
    call getmem1d(kbot,1,nipoi,'mod_cu_kf:kbot')
  end subroutine allocate_mod_cu_kf

  subroutine kfdrv(m2c,c2m)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    type(cum_2_mod) , intent(inout) :: c2m
    integer :: i , j ,  k , kk , np
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'kfdrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nipoi == 0 ) return

    do k = 1 , kz
      kk = kz - k + 1
      do np = 1 , nipoi
        i = imap(np)
        j = jmap(np)
        u0(np,k) = m2c%uas(j,i,kk)
        v0(np,k) = m2c%vas(j,i,kk)
        t0(np,k) = m2c%tas(j,i,kk)
        z0(np,k) = m2c%zas(j,i,kk)
        qv0(np,k) = m2c%qxas(j,i,kk,iqv)
        p0(np,k) = m2c%pas(j,i,kk)
        rho(np,k) = m2c%rhoas(j,i,kk)
        dzq(np,k) = m2c%dzq(j,i,kk)
        w0avg(np,k) = - d_half * (m2c%qdot(j,i,kk+1)-m2c%qdot(j,i,kk)) / &
                      (egrav*rho(np,k))
      end do
    end do

     !do k = 1 , kz
     !  do np = 1 , nipoi
     !    if ( np == 1 .and. myid == 1 ) then
     !      write(stdout,'(i2,7f8.0)') k, sqrt(u0(np,k)**2+v0(np,k)**2) , &
     !              t0(np,k) , 10000.*qv0(np,k) , &
     !              p0(np,k), 100*rho(np,k) , dzq(np,k) , &
     !              100000000*w0avg(np,k)
     !    end if
     !  end do
     !end do

    dtdt = d_zero
    dqdt = d_zero
    dqidt = d_zero
    dqcdt = d_zero
    dqrdt = d_zero
    dqsdt = d_zero
    raincv = d_zero
    pratec = d_zero
    ktop = -1
    kbot = -1

    if ( ipptls == 2 ) then
      call kfpara(1,kz,1,nipoi,.false.,.true.,.true.)
    else
      call kfpara(1,kz,1,nipoi,.true.,.false.,.false.)
    end if

    do k = 1 , kz
      kk = kz - k + 1
      do np = 1 , nipoi
        i = imap(np)
        j = jmap(np)
        c2m%tten(j,i,kk) = c2m%tten(j,i,kk) + dtdt(np,k)*m2c%psb(j,i)
        c2m%qxten(j,i,kk,iqv) = c2m%qxten(j,i,kk,iqv) + &
                dqdt(np,k)*m2c%psb(j,i)
        c2m%qxten(j,i,kk,iqc) = c2m%qxten(j,i,kk,iqc) + &
                dqcdt(np,k)*m2c%psb(j,i)
      end do
    end do

    if ( ipptls == 2 ) then
      do k = 1 , kz
        kk = kz - k + 1
        do np = 1 , nipoi
          i = imap(np)
          j = jmap(np)
          c2m%qxten(j,i,kk,iqr) = c2m%qxten(j,i,kk,iqr) + &
                  dqrdt(np,k)*m2c%psb(j,i)
          c2m%qxten(j,i,kk,iqi) = c2m%qxten(j,i,kk,iqi) + &
                  dqidt(np,k)*m2c%psb(j,i)
          c2m%qxten(j,i,kk,iqs) = c2m%qxten(j,i,kk,iqs) + &
                  dqsdt(np,k)*m2c%psb(j,i)
        end do
      end do
    end if

    do np = 1 , nipoi
      i = imap(np)
      j = jmap(np)
      if ( pratec(np) > dlowval ) then
        c2m%kcumtop(j,i) = kk-ktop(np)+1
        c2m%kcumbot(j,i) = kk-kbot(np)+1
        c2m%rainc(j,i) = c2m%rainc(j,i) + raincv(np)
        c2m%pcratec(j,i)= c2m%pcratec(j,i) + pratec(np)
        total_precip_points = total_precip_points + 1
      end if
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine kfdrv

  subroutine kfpara(kts,kte,its,ite,warm_rain,qi_flag,qs_flag)
    implicit none
    integer(ik4) , intent(in) :: kts , kte , its , ite
    logical , intent(in) :: warm_rain
    logical , intent(in) :: qi_flag , qs_flag

    real(rk8) , dimension(kts:kte) :: q0 , tv0 , tu , tvu , qu ,      &
            tz , tvd , qd , qes , thtes , tg , tvg , qg , wu , wd ,   &
            ems , emsd , umf , uer , udr , dmf , der , ddr ,          &
            umf2 , uer2 , udr2 , dmf2 , der2 , ddr2 , dza , thta0 ,   &
            thetee , thtau , theteu , thtad , theted , qliq , qice ,  &
            qlqout , qicout , pptliq , pptice , detlq , detic ,       &
            detlq2 , detic2 , ratio2

    real(rk8) , dimension(kts:kte) :: domgdp , exn , rhoe , tvqu ,    &
            dp , eqfrc , wspd , qdt , fxm , thtag , thtesg , thpa ,   &
            thfxtop , thfxbot , qpa , qfxtop , qfxbot , qlpa ,        &
            qlfxin , qlfxout , qipa , qifxin , qifxout , qrpa ,       &
            qrfxin , qrfxout , qspa , qsfxin , qsfxout , ql0 , qlg ,  &
            qi0 , qig , qr0 , qrg , qs0 , qsg

    real(rk8) , dimension(kts:kte+1) :: omg
    real(rk8) , dimension(kts:kte) :: rainfb , snowfb

    real(rk8) :: p300 , dpthmx , thmix , qmix , zmix , pmix ,          &
            rocpq , tmix , emix , tlog , tdpt , tlcl , tvlcl , cporq , &
            plcl , es , dlp , tenv , qenv , tven , tvbar , zlcl ,      &
            wkl , wabs , trppt , wsigne , dtlcl , gdt , wlcl , tvavg , &
            qese , wtw , rholcl , au0 , vmflcl , upold , upnew , abe , &
            wklcl , thtudl , tudl , ttemp , frc1 , qnewic , rl , r1 ,  &
            qnwfrz , effq , be , boterm , enterm , dzz , wsq , udlbe , &
            rei , ee2 , ud2 , ttmp , f1 , f2 , thttmp , qtmp ,         &
            tmpliq , tmpice , tu95 , tu10 , ee1 , ud1 , cldhgt ,       &
            dptt , qnewlq , dumfdp , ee , tsat , thta , p150 , usr ,   &
            vconv , timec , shsign , vws , pef , cbh , rcbh , pefcbh , &
            peff , peff2 , tder , thtmin , dtmltd , qs , tadvec ,      &
            dpdd , frc , dpt , rdd , a1 , dssdt , dtmp , t1rh , qsrh , &
            pptflx , cpr , cndtnf , updinc , aincm2 , devdmf , ppr ,   &
            rced , dpptdf , dmflfs , dmflfs2 , rced2 , ddinc ,         &
            aincmx , aincm1 , ainc , tder2 , pptfl2 , fabe , stab ,    &
            dtt , dtt1 , dtime , tma , tmb , tmm , bcoeff , acoeff ,   &
            qvdiff , topomg , cpm , dq , abeg , dabe , dfda , frc2 ,   &
            dr , udfrc , tuc , qgs , rh0 , rhg , qinit , qfnl , err2 , &
            relerr , rnc , fabeold , aincold , uefrc , ddfrc , tdc ,   &
            defrc

    integer(ik4) :: np , kx , k , kl
    integer(ik4) :: istop , ml , l5 , kmix , low , lc , mxlayr ,      &
            llfc , nlayrs , nk , kpbl , klcl , lcl , let , iflag ,    &
            kfrz , nk1 , ltop , nj , ltop1 , ltopm1 , lvf , kstart ,  &
            kmin , lfs , nd , nic , ldb , ldt , nd1 , ndk , nm ,      &
            lmax , ncount , noitr , nstep , ntc
    logical :: ddraft , doainc

    kl = kte
    kx = kte

    !
    ! checks for the possibility of initiating parameterized
    ! convection at each point
    !
    ! see if it is necessary to check for convective triggering at this
    ! grid point.
    !
    pointloop: &
    do np = its , ite
      p300 = p0(np,1) - 300.0D2
      !
      ! pressure perturbation term is only defined at mid-point of
      ! vertical layers. Since total pressure is needed at the top and
      ! bottom of layers below, do an interpolation
      !
      ! input a vertical sounding. Note that model layers are numbered
      ! from bottom-up in the kf scheme
      !
      ml = 0
      l5 = 0
      llfc = 0
      lvf = 0
      ldt = 0
      ldb = 0
      klcl = 0

      do k = 1 , kx
        ! if q0 is above saturation value, reduce it to saturation level
        es = aliq * exp( (bliq*t0(np,k) - cliq) / (t0(np,k)-dliq) )
        qes(k) = ep2 * es / (p0(np,k)-es)
        q0(k) = min(qes(k),qv0(np,k))
        q0(k) = max(minqx,q0(k))
        ql0(k) = d_zero
        qi0(k) = d_zero
        qr0(k) = d_zero
        qs0(k) = d_zero
        tv0(k) = t0(np,k)*(d_one + b61*q0(k))
        rhoe(k) = p0(np,k) / (rdry*tv0(k))
        dp(k) = rho(np,k) * egrav * dzq(np,k)
        !
        ! dzq is dz between sigma surfaces, dza is dz between model half level
        ! dp is the pressure interval between full sigma levels
        !
        if (p0(np,k) >= 500.0D2) l5 = k
        if (p0(np,k) >= p300)  llfc = k
        if (t0(np,k) >= tzero)   ml = k
      end do

      do k = 2 , kl
        dza(k-1) = z0(np,k) - z0(np,k-1)
      end do
      dza(kl) = d_zero
      kmix = 1
      kfmainloop: &
      do
        low = kmix
        if ( low > llfc ) then
          cycle pointloop
        end if

        lc = low
        mxlayr = 0
        !
        ! assume that in order to support a deep updraft you need a layer of
        ! unstable air 50 to 100 mb deep. To approximate this, isolate a
        ! group of adjacent individual model layers, with the base at level
        ! lc, such that the combined depth of these layers is at least 60 mb..
        !
        nlayrs = 0
        dpthmx = d_zero
        do nk = lc , kx
          dpthmx = dpthmx + dp(nk)
          nlayrs = nlayrs + 1
          if ( dpthmx > 6.0D3 ) exit
        end do
        if ( dpthmx < 6.0D3 ) then
          cycle pointloop
        end if
        kpbl = lc + nlayrs - 1
        kmix = lc + 1
        thmix = d_zero
        qmix = d_zero
        zmix = d_zero
        pmix = d_zero
        dpthmx = d_zero
        !
        ! find the thermodynamic characteristics of the layer by
        ! mass-weighting the characteristics of the individual model
        ! layers
        !
        do nk = lc , kpbl
          dpthmx = dpthmx + dp(nk)
          rocpq = c6 * (d_one - c7*q0(nk))
          thmix = thmix + dp(nk)*t0(np,nk) * (p00/p0(np,nk))**rocpq
          qmix = qmix + dp(nk)*q0(nk)
          zmix = zmix + dp(nk)*z0(np,nk)
          pmix = pmix + dp(nk)*p0(np,nk)
        end do
        thmix = thmix/dpthmx
        qmix = qmix/dpthmx
        zmix = zmix/dpthmx
        pmix = pmix/dpthmx
        rocpq = c6*(d_one - c7*qmix)
        tmix = thmix * (pmix/p00)**rocpq
        emix = qmix * pmix/(ep2+qmix)
        !
        ! find the temperature of the mixture at its lcl, pressure
        ! level of lcl
        !
        tlog = log(emix/aliq)
        tdpt = (cliq - dliq*tlog) / (bliq-tlog)
        tlcl = tdpt - &
                (tl1 + tl2*(tdpt-tzero) - tl3*(tmix-tzero))*(tmix-tdpt)
        tlcl = min(tlcl,tmix)
        tvlcl = tlcl*(d_one + b61*qmix)
        cporq = d_one / rocpq

        plcl = p00 * (tlcl/thmix)**cporq
        find_klcl: &
        do nk = lc , kl
          klcl = nk
          if ( plcl >= p0(np,nk) ) exit find_klcl
        end do find_klcl
        !if ( klcl == kl ) then
        !  print *, 'Exit at klcl == kl',klcl,kl,plcl,p0(np,klcl)
        !  cycle pointloop
        !end if
        k = klcl - 1
        dlp = log(plcl/p0(np,k)) / log(p0(np,klcl)/p0(np,k))
        !
        ! estimate environmental temperature and mixing ratio at the lcl
        !
        tenv = t0(np,k) + (t0(np,klcl)-t0(np,k))*dlp
        qenv = q0(k) + (q0(klcl)-q0(k))*dlp
        tven = tenv * (d_one + b61*qenv)
        tvbar = d_half * (tv0(k)+tven)
        zlcl = z0(np,k) + (z0(np,klcl)-z0(np,k))*dlp
        !
        ! check to see if cloud is buoyant using fritsch-chappell trigger
        ! function described in kain and fritsch (1992). W0avg is an
        ! aproximate value for the running-mean grid-scale vertical
        ! velocity, which gives smoother fields of convective initiation
        ! than the instantaneous value. Formula relating temperature
        ! perturbation to vertical velocity has been used with the most
        ! success at grid lengths near 25 km.  for different grid-lengths,
        ! adjust vertical velocity to equivalent value for 25 km grid
        ! length, assuming linear dependence of w on grid length
        !
        wklcl = 0.02D0 * zlcl/2.5D3
        wkl = (w0avg(np,k) + (w0avg(np,klcl) - w0avg(np,k))*dlp) * &
                (dx/25.D3) - wklcl
        wabs = abs(wkl) + 1.D-10
        wsigne = wkl/wabs
        dtlcl = 4.64D0 * wsigne * wabs**0.33D0
        gdt = egrav * dtlcl * (zlcl-z0(np,lc)) / (tv0(lc)+tven)
        wlcl = d_one + d_half * wsigne * sqrt(abs(gdt)+1.D-10)
        if ( tlcl+dtlcl <= tenv ) cycle kfmainloop
        if ( kpbl >= llfc ) then
          print *, 'Exit at kpbl >= llfc'
          cycle pointloop
        end if
        !
        ! convective triggering criteria has been satisfied. Compute
        ! equivalent potential temperature
        ! (theteu) and vertical velocity of the rising parcel at the lcl
        !
        theteu(k) = tmix * (p00/pmix)**(c6*(d_one - c7*qmix)) * &
                     exp((c1/tlcl - c2)*qmix*(d_one + c8*qmix))
        es = aliq * exp((tenv*bliq - cliq) / (tenv-dliq))
        tvavg = d_half * (tv0(klcl)+tenv*(d_one + b61*qenv))
        plcl = p0(np,klcl) * exp(egrav / (rdry*tvavg) * (z0(np,klcl)-zlcl))
        qese = ep2 * es / (plcl-es)
        gdt = egrav * dtlcl * (zlcl-z0(np,lc)) / (tv0(lc)+tven)
        wlcl = d_one + d_half * wsigne * sqrt(abs(gdt) + 1.D-10)
        thtes(k) = tenv*(p00/plcl)**(c6*(d_one - c7*qese)) * &
                   exp((c1/tenv - c2) * qese * (d_one + c8*qese))
        wtw = wlcl * wlcl
        if ( wlcl < d_zero ) cycle kfmainloop
        tvlcl = tlcl * (d_one + b61*qmix)
        rholcl = plcl / (rdry*tvlcl)
        lcl = klcl
        let = lcl
        !
        !*******************************************************************
        !                                                                  *
        !                 compute updraft properties                       *
        !                                                                  *
        !*******************************************************************
        !
        !
        ! estimate initial updraft mass flux (umf(k))
        !
        wu(k) = wlcl
        au0 = mathpi*rad*rad
        umf(k) = rholcl*au0
        vmflcl = umf(k)
        upold = vmflcl
        upnew = upold
        !
        ! ratio2 is the degree of glaciation in the cloud (0 to 1),
        ! uer is the envir entrainment rate, abe is available buoyant energy,
        ! trppt is the total rate of precipitation production
        !
        ratio2(k) = d_zero
        uer(k) = d_zero
        abe = d_zero
        trppt = d_zero
        tu(k) = tlcl
        tvu(k) = tvlcl
        qu(k) = qmix
        eqfrc(k) = d_one
        qliq(k) = d_zero
        qice(k) = d_zero
        qlqout(k) = d_zero
        qicout(k) = d_zero
        detlq(k) = d_zero
        detic(k) = d_zero
        pptliq(k) = d_zero
        pptice(k) = d_zero
        iflag = 0
        kfrz = lc
        !
        ! the amount of conv avail pot energy (cape) is calculated with
        ! respect to undilute parcel ascent; eq pot temp of undilute
        ! parcel is thtudl, undilute temperature is given by tudl
        !
        thtudl = theteu(k)
        tudl = tlcl
        !
        ! ttemp is used during calculation of the linear glaciation
        ! process; it is initially set to the temperature at which
        ! freezing is specified to begin.  within the glaciation
        ! interval, it is set equal to the updraft temp at the
        ! previous model level
        !
        ttemp = ttfrz
        !
        ! enter the loop for updraft calculations. Calculate updraft temp,
        ! mixing ratio, vertical mass flux, lateral detrainment of mass and
        ! moisture, precipitation rates at each model level
        !
        ud1 = d_zero
        ee1 = d_one
        updraft_loop: &
        do nk = k , kl-1
          nk1 = nk + 1
          ratio2(nk1) = ratio2(nk)
          !
          ! update updraft properties at the next model lvl to reflect
          ! entrainment of environmental air
          !
          frc1 = d_zero
          tu(nk1) = t0(np,nk1)
          theteu(nk1) = theteu(nk)
          qu(nk1) = qu(nk)
          qliq(nk1) = qliq(nk)
          qice(nk1) = qice(nk)

          call tpmix(p0(np,nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1), &
                     qice(nk1),qnewlq,qnewic,ratio2(nk1),rl)
          tvu(nk1) = tu(nk1) * (d_one + b61*qu(nk1))
          !
          ! check to see if updraft temp is within the freezing interval,
          ! if it is, calculate the fractional conversion to glaciation
          ! and adjust qnewlq to reflect the gradual change in thetau
          ! since the last model level. The glaciation effects will be
          ! determined after the amount of condensate available after
          ! precip fallout is determined. Ttfrz is the temp at which
          ! glaciation begins, tbfrz the temp at which it ends
          !
          if ( tu(nk1) <= ttfrz .and. iflag < 1 ) then
            if ( tu(nk1) > tbfrz ) then
              if ( ttemp > ttfrz ) ttemp = ttfrz
              frc1 = (ttemp-tu(nk1)) / (ttfrz-tbfrz)
              r1 = (ttemp-tu(nk1)) / (ttemp-tbfrz)
            else
              frc1 = (ttemp-tbfrz) / (ttfrz-tbfrz)
              r1 = d_one
              iflag = 1
            end if
            qnwfrz = qnewlq
            qnewic = qnewic + qnewlq*r1*d_half
            qnewlq = qnewlq - qnewlq*r1*d_half
            effq = (ttfrz-tbfrz) / (ttemp-tbfrz)
            ttemp = tu(nk1)
          end if
          !
          ! calculate updraft vertical velocity and precipitation fallout
          !
          if ( nk == k ) then
            be = (tvlcl+tvu(nk1)) / (tven+tv0(nk1)) - d_one
            boterm = d_two * (z0(np,nk1)-zlcl) * egrav * be/1.5D0
            enterm = d_zero
            dzz = z0(np,nk1) - zlcl
          else
            be = (tvu(nk)+tvu(nk1)) / (tv0(nk)+tv0(nk1)) - d_one
            boterm = d_two * dza(nk) *egrav * be/1.5D0
            enterm = d_two * uer(nk) * wtw/upold
            dzz = dza(nk)
          end if
          wsq = wtw
          call condload(qliq(nk1),qice(nk1),wtw,dzz,boterm,enterm,rate, &
                        qnewlq,qnewic,qlqout(nk1),qicout(nk1))
          !
          ! if vert velocity is less than zero, exit the updraft loop and,
          ! if cloud is tall enough, finalize updraft calculations
          !
          if ( wtw <= d_zero ) exit updraft_loop
          wabs = sqrt(abs(wtw))
          wu(nk1) = wtw / wabs
          !
          ! update the abe for undilute ascent
          !
          thtes(nk1) = t0(np,nk1) * &
                  (p00/p0(np,nk1))**(c6*(d_one-c7*qes(nk1))) * &
                   exp((c1/t0(np,nk1) - c2)*qes(nk1) * (d_one + c8*qes(nk1)))
          udlbe = ((d_two * thtudl) / (thtes(nk) + thtes(nk1))-d_one) * dzz
          if ( udlbe > d_zero ) abe = abe + udlbe*egrav
          !
          ! determine the effects of cloud glaciation if within the specified
          ! temp interval
          !
          if ( frc1 > 1.D-6 ) then
            call dtfrznew(tu(nk1),p0(np,nk1),theteu(nk1),qu(nk1),qliq(nk1), &
                          qice(nk1),ratio2(nk1),qnwfrz,rl,frc1,effq,iflag)
          end if
          !
          ! call subroutine to calculate environmental equivalent potential
          ! temp. within glaciation interval, thetae must be calculated with
          ! respect to same degree of glaciation for all entraining air
          !
          call envirtht(p0(np,nk1),t0(np,nk1),q0(nk1), &
                        thetee(nk1),ratio2(nk1),rl)
          !
          ! rei is the rate of environmental inflow
          !
          rei = vmflcl * dp(nk1) * 0.03D0/rad
          tvqu(nk1) = tu(nk1) * (d_one + b61 * qu(nk1) - qliq(nk1) - qice(nk1))
          !
          ! if cloud parcels are virtually colder than the environment, no
          ! entrainment is allowed at this level
          !
          if ( tvqu(nk1) <= tv0(nk1) ) then
            uer(nk1) = d_zero
            udr(nk1) = rei
            ee2 = d_zero
            ud2 = d_one
            eqfrc(nk1) = d_zero
          else
            let = nk1
            ttmp = tvqu(nk1)
            !
            ! determine the critical mixed fraction of updraft and environmental
            ! air for estimation of entrainment and detrainment rates
            f1 = 0.95D0
            f2 = d_one - f1
            thttmp = f1 * thetee(nk1) + f2 * theteu(nk1)
            qtmp = f1 * q0(nk1) + f2*qu(nk1)
            tmpliq = f2*qliq(nk1)
            tmpice = f2*qice(nk1)
            call tpmix(p0(np,nk1),thttmp,ttmp,qtmp,tmpliq, &
                       tmpice,qnewlq,qnewic,ratio2(nk1),rl)
            tu95 = ttmp * (d_one + b61*qtmp - tmpliq - tmpice)
            mixed_fraction_compute: &
            do 
              if ( tu95 > tv0(nk1) ) then
                ee2 = d_one
                ud2 = d_zero
                eqfrc(nk1) = d_one
                exit mixed_fraction_compute
              end if
              f1 = 0.10D0
              f2 = d_one - f1
              thttmp = f1 * thetee(nk1) + f2 * theteu(nk1)
              qtmp = f1 * q0(nk1) + f2 * qu(nk1)
              tmpliq = f2 * qliq(nk1)
              tmpice = f2 * qice(nk1)
              call tpmix(p0(np,nk1),thttmp,ttmp,qtmp,tmpliq,tmpice, &
                         qnewlq,qnewic,ratio2(nk1),rl)
              tu10 = ttmp * (d_one + b61*qtmp - tmpliq - tmpice)
              if ( abs(tu10-tvqu(nk1)) < dlowval ) then
                ee2 = d_one
                ud2 = d_zero
                eqfrc(nk1) = d_one
                exit mixed_fraction_compute
              end if
              eqfrc(nk1) = (tv0(nk1)-tvqu(nk1)) * f1 / (tu10-tvqu(nk1))
              eqfrc(nk1) = max(d_zero,eqfrc(nk1))
              eqfrc(nk1) = min(d_one,eqfrc(nk1))
              if ( d_one - eqfrc(nk1) < dlowval ) then
                ee2 = d_one
                ud2 = d_zero
                exit mixed_fraction_compute
              else if ( eqfrc(nk1) < dlowval ) then
                ee2 = d_zero
                ud2 = d_one
                exit mixed_fraction_compute
              else
                !
                ! subroutine prof5 integrates over the gaussian dist to
                ! determine the fractional entrainment and detrainment rates
                !
                call prof5(eqfrc(nk1),ee2,ud2)
              end if
              exit mixed_fraction_compute
            end do mixed_fraction_compute
            !
            ! net entrainment and detrainment rates are given by the average
            ! fractional values in the layer
            !
            if ( nk == k ) then
              ee1 = d_one
              ud1 = d_zero
            end if
            uer(nk1) = d_half*rei*(ee1+ee2)
            udr(nk1) = d_half*rei*(ud1+ud2)
          end if
          !
          ! if the calculated updraft detrainment rate is greater than the total
          ! updraft mass flux, all cloud mass detrains, exit updraft calculation
          !
          if ( umf(nk) - udr(nk1) < 10.0D0 ) then
            !
            ! If the calculated detrained mass flux is greater than the total
            ! updraft flux, impose total detrainment of updraft mass at the
            ! previous model
            !
            if ( udlbe > d_zero ) abe = abe - udlbe*egrav
            let = nk
            exit updraft_loop
          end if
          ee1 = ee2
          ud1 = ud2
          upold = umf(nk) - udr(nk1)
          upnew = upold + uer(nk1)
          umf(nk1) = upnew
          !
          ! detlq and detic are the rates of detrainment of liquid and ice in
          ! the detraining updraft mass
          !
          detlq(nk1) = qliq(nk1)*udr(nk1)
          detic(nk1) = qice(nk1)*udr(nk1)
          qdt(nk1) = qu(nk1)
          qu(nk1) = (upold*qu(nk1) + uer(nk1)*q0(nk1))/upnew
          theteu(nk1) = (theteu(nk1)*upold + thetee(nk1)*uer(nk1))/upnew
          qliq(nk1) = qliq(nk1) * upold/upnew
          qice(nk1) = qice(nk1) * upold/upnew
          !
          ! kfrz is the highest model level at which liquid condensate is
          ! generating pptliq is the rate of generation (fallout) of liquid
          ! precip at a giving model lvl, pptice the same for ice, trppt is
          ! the total rate of production of precip up to the current model level
          !
          if ( abs(ratio2(nk1)-d_one) > 1.D-6 ) kfrz = nk1
          pptliq(nk1) = qlqout(nk1) * (umf(nk) - udr(nk1))
          pptice(nk1) = qicout(nk1) * (umf(nk) - udr(nk1))
          trppt = trppt + pptliq(nk1) + pptice(nk1)
          if ( nk1 <= kpbl ) uer(nk1) = uer(nk1) + vmflcl * dp(nk1)/dpthmx
        end do updraft_loop
        !
        ! check cloud depth. If cloud is tall enough, estimate the equilibriu
        ! temperature level (let) and adjust mass flux profile at cloud top so
        ! that mass flux decreases to zero as a linear function of pressure
        ! between the let and cloud top
        !
        ! ltop is the model level just below the level at which vertical
        ! velocity first becomes negative
        !
        ltop = nk
        cldhgt = z0(np,ltop) - zlcl
        !
        ! if cloud top hgt is less than specified minimum height, go back and
        ! the next highest 60mb layer to see if a bigger cloud can be obtained
        ! that source air
        !
        if ( cldhgt < 3.D3 .or. abe < d_one ) then
          do nk = k , ltop
            umf(nk) = d_zero
            udr(nk) = d_zero
            uer(nk) = d_zero
            detlq(nk) = d_zero
            detic(nk) = d_zero
            pptliq(nk) = d_zero
            pptice(nk) = d_zero
          end do
          cycle kfmainloop
        end if
        !
        ! if the let and ltop are the same, detrain all of the updraft mass
        ! flux this level
        !
        if ( let == ltop ) then
          udr(ltop) = umf(ltop) + udr(ltop) - uer(ltop)
          detlq(ltop) = qliq(ltop) * udr(ltop) * upnew/upold
          detic(ltop) = qice(ltop) * udr(ltop) * upnew/upold
          trppt = trppt - (pptliq(ltop) + pptice(ltop))
          uer(ltop) = d_zero
          umf(ltop) = d_zero
        else
          !
          ! begin total detrainment at the level above the let
          !
          dptt = d_zero
          do nj = let+1 , ltop
            dptt = dptt + dp(nj)
          end do
          dumfdp = umf(let)/dptt
          !
          ! adjust mass flux profiles, detrainment rates, and precipitation fall
          ! rates to reflect the linear decrease in mass flx between the let and
          ! ptop
          !
          do nk = let+1 , ltop
            udr(nk) = dp(nk) * dumfdp
            umf(nk) = umf(nk-1) - udr(nk)
            detlq(nk) = qliq(nk) * udr(nk)
            detic(nk) = qice(nk) * udr(nk)
            trppt = trppt - pptliq(nk) - pptice(nk)
            pptliq(nk) = (umf(nk-1) - udr(nk)) * qlqout(nk)
            pptice(nk) = (umf(nk-1) - udr(nk)) * qicout(nk)
            trppt = trppt + pptliq(nk) + pptice(nk)
          end do
        end if
        !
        ! extend the updraft mass flux profile down to the source layer for
        ! the updraft air. Also, define thetae for levels below the lcl
        !
        do nk = 1 , k
          if ( nk >= lc ) then
            if ( nk == lc ) then
              umf(nk) = vmflcl * dp(nk)/dpthmx
              uer(nk) = vmflcl * dp(nk)/dpthmx
            else if ( nk <= kpbl ) then
              uer(nk) = vmflcl * dp(nk)/dpthmx
              umf(nk) = umf(nk-1) + uer(nk)
            else
              umf(nk) = vmflcl
              uer(nk) = d_zero
            end if
            tu(nk) = tmix + (z0(np,nk)-zmix)*gdry
            qu(nk) = qmix
            wu(nk) = wlcl
          else
            tu(nk) = d_zero
            qu(nk) = d_zero
            umf(nk) = d_zero
            wu(nk) = d_zero
            uer(nk) = d_zero
          end if
          udr(nk) = d_zero
          qdt(nk) = d_zero
          qliq(nk) = d_zero
          qice(nk) = d_zero
          qlqout(nk) = d_zero
          qicout(nk) = d_zero
          pptliq(nk) = d_zero
          pptice(nk) = d_zero
          detlq(nk) = d_zero
          detic(nk) = d_zero
          ratio2(nk) = d_zero
          ee = q0(nk) * p0(np,nk) / (ep2+q0(nk))
          tlog = log(ee/aliq)
          tdpt = (cliq - dliq*tlog) / (bliq-tlog)
          tsat = tdpt - (tl1 + tl2*(tdpt-tzero) - &
                         tl3 * (t0(np,nk)-tzero)) * (t0(np,nk)-tdpt)
          thta = t0(np,nk) * (p00/p0(np,nk))**(c6*(d_one - c7*q0(nk)))
          thetee(nk) = thta * exp((c1/tsat - c2)*q0(nk)*(d_one + c8*q0(nk)))
          thtes(nk) = thta * exp((c1/t0(np,nk) - c2)*qes(nk) * &
                      (d_one + c8*qes(nk)))
          eqfrc(nk) = d_one
        end do
        ltop1 = ltop + 1
        ltopm1 = ltop - 1
        !
        ! define variables above cloud top
        !
        do nk = ltop1 , kx
          umf(nk) = d_zero
          udr(nk) = d_zero
          uer(nk) = d_zero
          qdt(nk) = d_zero
          qliq(nk) = d_zero
          qice(nk) = d_zero
          qlqout(nk) = d_zero
          qicout(nk) = d_zero
          detlq(nk) = d_zero
          detic(nk) = d_zero
          pptliq(nk) = d_zero
          pptice(nk) = d_zero
          if ( nk > ltop1 ) then
            tu(nk) = d_zero
            qu(nk) = d_zero
            wu(nk) = d_zero
          end if
          thta0(nk) = d_zero
          thtau(nk) = d_zero
          ems(nk) = dp(nk) * dxsq * regrav
          emsd(nk) = d_one / ems(nk)
          tg(nk) = t0(np,nk)
          qg(nk) = q0(nk)
          qlg(nk) = d_zero
          qig(nk) = d_zero
          qrg(nk) = d_zero
          qsg(nk) = d_zero
          omg(nk) = d_zero
        end do
        omg(kl+1) = d_zero
        p150 = p0(np,klcl) - 1.50D4
        do nk = 1 , ltop
          thtad(nk) = d_zero
          ems(nk) = dp(nk) * dxsq * regrav
          emsd(nk) = d_one / ems(nk)
          !
          ! initialize some variables to be used later in the vert advection
          ! scheme
          !
          exn(nk) = (p00/p0(np,nk))**(c6*(d_one - c7*qdt(nk)))
          thtau(nk) = tu(nk)*exn(nk)
          exn(nk) = (p00/p0(np,nk))**(c6*(d_one - c7*q0(nk)))
          thta0(nk) = t0(np,nk)*exn(nk)
          !
          ! lvf is the level at which moisture flux is estimated as the basis
          ! for precipitation efficiency calculations
          !
          if ( p0(np,nk) > p150 ) lvf = nk
          omg(nk) = d_zero
        end do
        lvf = min(lvf,let)
        usr = umf(lvf+1) * (qu(lvf+1) + qliq(lvf+1) + qice(lvf+1))
        usr = min(usr,trppt)
        if (usr < 1.D-8 ) usr = trppt
        !
        ! compute convective time scale(timec). the mean wind at the lcl
        ! and midtroposphere is used.
        !
        wspd(klcl) = sqrt(u0(np,klcl)*u0(np,klcl) + v0(np,klcl)*v0(np,klcl))
        wspd(l5) = sqrt(u0(np,l5)*u0(np,l5) + v0(np,l5)*v0(np,l5))
        wspd(ltop) = sqrt(u0(np,ltop)*u0(np,ltop) + v0(np,ltop)*v0(np,ltop))
        vconv = d_half * (wspd(klcl) + wspd(l5))
        if ( vconv > d_zero ) then
          timec = dx/vconv
        else
          timec = 3600.D0
        end if
        tadvec = timec
        timec = max(1800.0D0,timec)
        timec = min(3600.0D0,timec)
        nic = nint(timec/dtsec)
        timec = dble(nic)*dtsec
        !
        ! compute wind shear and precipitation efficiency.
        !
        if ( wspd(ltop) > wspd(klcl) ) then
          shsign = d_one
        else
          shsign = -d_one
        end if
        vws = (u0(np,ltop)-u0(np,klcl)) * (u0(np,ltop)-u0(np,klcl)) + &
              (v0(np,ltop)-v0(np,klcl)) * (v0(np,ltop)-v0(np,klcl))
        vws = 1.0D3 * shsign * sqrt(vws) / (z0(np,ltop)-z0(np,lcl))
        pef = 1.591D0 + vws*(-0.639D0 + vws * (9.53D-2 - vws*4.96D-3))
        pef = max(pef,0.2D0)
        pef = min(pef,0.9D0)
        !
        ! precipitation efficiency is a function of the height of cloud base.
        !
        cbh = (zlcl-z0(np,1)) * 3.281D-3
        if ( cbh < 3.0D0 ) then
          rcbh = 0.02D0
        else
          rcbh = 0.96729352D0 + cbh*(-0.70034167D0 + &
                                cbh*(0.162179896D0 + &
                                cbh*(-1.2569798D-2 + &
                                cbh*(4.2772D-4 - &
                                cbh*5.44D-6))))
        end if
        if ( cbh > 25.0D0 ) rcbh = 2.4D0
        pefcbh = d_one / (d_one+rcbh)
        pefcbh = min(pefcbh,0.9D0)
        !
        ! mean pef. is used to compute rainfall.
        !
        peff = d_half * (pef+pefcbh)
        peff2 = peff
        !
        !*****************************************************************
        !                                                                *
        !                  compute downdraft properties                  *
        !                                                                *
        !*****************************************************************
        !
        ! let downdraft originate at the level of minimum saturation equivalen
        ! potential temperature (seqt) in the cloud layer, extend downward to
        ! surface, or to the layer below cloud base at which envir seqt is les
        ! than min seqt in the cloud layer
        ! let downdraft detrain over a laye of specified pressure-depth (dpdd)
        !
        tder = d_zero
        kstart = max(kpbl,klcl)
        thtmin = thtes(kstart+1)
        kmin = kstart+1
        do nk = kstart+2 , ltop-1
          thtmin = min(thtmin,thtes(nk))
          if ( abs(thtmin-thtes(nk)) < dlowval ) kmin = nk
        end do
        lfs = kmin
        if ( ratio2(lfs) > d_zero ) then
          call envirtht(p0(np,lfs),t0(np,lfs),q0(lfs),thetee(lfs),d_zero,rl)
        end if
        eqfrc(lfs) = (thtes(lfs)-theteu(lfs)) / (thetee(lfs)-theteu(lfs))
        eqfrc(lfs) = max(eqfrc(lfs),d_zero)
        eqfrc(lfs) = min(eqfrc(lfs),d_one)
        theted(lfs) = thtes(lfs)
        !
        ! estimate the effect of melting precipitation in the downdraft
        !
        if ( ml > 0 ) then
          dtmltd = d_half * (qu(klcl)-qu(ltop)) * rlf/cpd
        else
          dtmltd = d_zero
        end if
        tz(lfs) = t0(np,lfs) - dtmltd
        es = aliq*exp((tz(lfs)*bliq - cliq) / (tz(lfs)-dliq))
        qs = ep2 * es / (p0(np,lfs)-es)
        qd(lfs) = eqfrc(lfs) * q0(lfs) + (d_one-eqfrc(lfs)) * qu(lfs)
        thtad(lfs) = tz(lfs) * (p00/p0(np,lfs))**(c6*(d_one - c7*qd(lfs)))
        if ( qd(lfs) >= qs ) then
          theted(lfs) = thtad(lfs)*exp((c1/tz(lfs) - c2)*qs*(d_one + c8*qs))
        else
          call envirtht(p0(np,lfs),tz(lfs),qd(lfs),theted(lfs),d_zero,rl)
        end if
        ddraft = .true.
        ddraft_find: &
        do nk = 1 , lfs
          nd = lfs - nk
          if ( theted(lfs) > thtes(nd) .or. nd == 1 ) then
            ldb = nd
            !
            ! if downdraft never becomes negatively buoyant or if it
            ! is shallower 50 mb, don't allow it to occur at all
            !
            if ( nk == 1 .or. (p0(np,ldb)-p0(np,lfs)) < 50.D2 ) then
              ddraft = .false.
            end if
            exit ddraft_find
          end if
        end do ddraft_find
        !
        ! allow downdraft to detrain in a single layer, but with downdraft air
        ! typically flushed up into higher layers as allowed in the total
        ! vertical advection calculations farther down in the code
        !
        if ( ddraft ) then
          dpdd = dp(ldb)
          ldt = ldb
          frc = d_one
          dpt = d_zero
          !
          ! take a first guess at the initial downdraft mass flux..
          !
          tvd(lfs) = t0(np,lfs) * (d_one + b61*qes(lfs))
          rdd = p0(np,lfs) / (rdry*tvd(lfs))
          a1 = (d_one-peff) * au0
          dmf(lfs) = -a1*rdd
          der(lfs) = eqfrc(lfs)*dmf(lfs)
          ddr(lfs) = d_zero
          do nd = lfs-1 , ldb,-1
            nd1 = nd + 1
            if ( nd <= ldt ) then
              der(nd) = d_zero
              ddr(nd) = -dmf(ldt+1) * dp(nd) * frc/dpdd
              dmf(nd) = dmf(nd1) + ddr(nd)
              frc = d_one
              theted(nd) = theted(nd1)
              qd(nd) = qd(nd1)
            else
              der(nd) = dmf(lfs) * 0.03D0 * dp(nd)/rad
              ddr(nd) = d_zero
              dmf(nd) = dmf(nd1) + der(nd)
              if ( ratio2(nd) > d_zero ) then
                call envirtht(p0(np,nd),t0(np,nd),q0(nd),thetee(nd),d_zero,rl)
              end if
              theted(nd) = (theted(nd1)*dmf(nd1) + thetee(nd) * der(nd))/dmf(nd)
              qd(nd) = (qd(nd1) * dmf(nd1) + q0(nd)*der(nd)) / dmf(nd)
            end if
          end do
          tder = d_zero
          !
          ! calculation an evaporation rate for given mass flux
          !
          do nd = ldb , ldt
            tz(nd) = tpdd(p0(np,nd),theted(ldt),t0(np,nd),qs,qd(nd),d_one)
            es = aliq * exp((tz(nd)*bliq - cliq) / (tz(nd)-dliq))
            qs = ep2 * es / (p0(np,nd)-es)
            dssdt = (cliq - bliq*dliq) / ((tz(nd)-dliq)*(tz(nd)-dliq))
            rl = xlv0 - xlv1*tz(nd)
            dtmp = rl * qs * (d_one-rhbc) / (cpd + rl*rhbc*qs*dssdt)
            t1rh = tz(nd) + dtmp
            es = rhbc * aliq * exp((bliq*t1rh - cliq) / (t1rh-dliq))
            qsrh = ep2 * es / (p0(np,nd)-es)
            !
            ! check to see if mixing ratio at specified rh is less than actual
            ! mixing ratio. if so, adjust to give zero evaporation
            !
            if ( qsrh < qd(nd) ) then
              qsrh = qd(nd)
              t1rh = tz(nd)
            end if
            tz(nd) = t1rh
            qs = qsrh
            tder = tder + (qs-qd(nd)) * ddr(nd)
            qd(nd) = qs
            thtad(nd) = tz(nd)*(p00/p0(np,nd))**(c6*(d_one-c7*qd(nd)))
          end do
        end if
        !
        ! if downdraft does not evaporate any water for specified relative
        ! humidity, no downdraft is allowed
        !
        evaploop: &
        do
          if ( tder < d_one ) then
            pptflx = trppt
            cpr = trppt
            tder = d_zero
            cndtnf = d_zero
            updinc = d_one
            ldb = lfs
            do ndk = 1 , ltop
              dmf(ndk) = d_zero
              der(ndk) = d_zero
              ddr(ndk) = d_zero
              thtad(ndk) = d_zero
              wd(ndk) = d_zero
              tz(ndk) = d_zero
              qd(ndk) = d_zero
            end do
            aincm2 = 100.0D0
            exit evaploop
          end if
          !
          ! adjust downdraft mass flux so that evaporation rate in downdraft
          ! is consistent with precipitation efficiency relationship
          !
          devdmf = tder / dmf(lfs)
          ppr = d_zero
          pptflx = peff * usr
          rced = trppt - pptflx
          !
          ! ppr is the total amount of precipitation that falls  out of the
          ! updraft from cloud base to the lfs. updraft mass flux will be
          ! increased up to the lfs to account for updraft air mixing with
          ! environmental air to the updraft, so ppr will increase
          ! proportionately
          !
          do nm = klcl , lfs
            ppr = ppr + pptliq(nm) + pptice(nm)
          end do
          if ( lfs >= klcl ) then
            dpptdf = (d_one-peff) * ppr * (d_one-eqfrc(lfs)) / umf(lfs)
          else
            dpptdf = d_zero
          end if
          !
          ! cndtnf is the amount of condensate transferred along with updraft
          ! mass the downdraft at the lfs
          !
          cndtnf = (qliq(lfs) + qice(lfs)) * (d_one - eqfrc(lfs))
          dmflfs = rced / (devdmf + dpptdf + cndtnf)
          if ( dmflfs > d_zero ) then
            tder = d_zero
            cycle evaploop
          end if
          !
          ! ddinc is the factor by which to increase the first-guess
          ! downdraft mass flux to satisfy the precip efficiency relationship,
          ! updinc is t which to increase the updraft mass flux below the lfs
          ! to account for transfer of mass from updraft to downdraft
          !
          ! ddinc=dmflfs/dmf(lfs)
          if ( lfs >= klcl ) then
            updinc = (umf(lfs) - (d_one - eqfrc(lfs)) * dmflfs) / umf(lfs)
            !
            ! limit updinc to less than or equal to 1.5
            !
            if ( updinc > 1.5D0 ) then
              updinc = 1.5D0
              dmflfs2 = umf(lfs) * (updinc-d_one) / (eqfrc(lfs)-d_one)
              rced2 = dmflfs2 * (devdmf + dpptdf + cndtnf)
              pptflx = pptflx + (rced-rced2)
              peff2 = pptflx/usr
              rced = rced2
              dmflfs = dmflfs2
            end if
          else
            updinc = d_one
          end if
          ddinc = dmflfs / dmf(lfs)
          do nk = ldb , lfs
            dmf(nk) = dmf(nk)*ddinc
            der(nk) = der(nk)*ddinc
            ddr(nk) = ddr(nk)*ddinc
          end do
          cpr = trppt + ppr*(updinc-d_one)
          pptflx = pptflx + peff*ppr*(updinc-d_one)
          peff = peff2
          tder = tder*ddinc
          !
          ! adjust updraft mass flux, mass detrainment rate, and liquid
          ! water and detrainment rates to be consistent with the transfer
          ! of the estimate from the updraft to the downdraft at the lfs
          !
          do nk = lc , lfs
            umf(nk) = umf(nk) * updinc
            udr(nk) = udr(nk) * updinc
            uer(nk) = uer(nk) * updinc
            pptliq(nk) = pptliq(nk) * updinc
            pptice(nk) = pptice(nk) * updinc
            detlq(nk) = detlq(nk) * updinc
            detic(nk) = detic(nk) * updinc
          end do
          !
          ! zero out the arrays for downdraft data at levels above and
          ! below the downdraft
          !
          if ( ldb > 1 ) then
            do nk = 1 , ldb-1
              dmf(nk) = d_zero
              der(nk) = d_zero
              ddr(nk) = d_zero
              wd(nk) = d_zero
              tz(nk) = d_zero
              qd(nk) = d_zero
              thtad(nk) = d_zero
            end do
          end if
          do nk = lfs + 1 , kx
            dmf(nk) = d_zero
            der(nk) = d_zero
            ddr(nk) = d_zero
            wd(nk) = d_zero
            tz(nk) = d_zero
            qd(nk) = d_zero
            thtad(nk) = d_zero
          end do
          do nk = ldt + 1 , lfs - 1
            tz(nk) = d_zero
            qd(nk) = d_zero
          end do
          exit evaploop
        end do evaploop
        !
        !
        ! set limits on the updraft and downdraft mass fluxes so that the
        ! inflow into convective drafts from a given layer is no more than
        ! is available in that layer initially
        !
        aincmx = 1000.0D0
        aincm1 = ems(lc)
        lmax = max(klcl,lfs)
        do nk = lc , lmax
          if ( (uer(nk)-der(nk)) > d_zero ) then
            aincm1 = ems(nk) / ((uer(nk)-der(nk)) * timec)
          end if
          aincmx = min(aincmx,aincm1)
        end do
        ainc = d_one
        if ( aincmx < ainc ) ainc = aincmx
        !
        ! save the relevent variables for a unit updrft and downdrft. they
        ! will iteratively adjusted by the factor ainc to satisfy the
        ! stabilization closure
        !
        ncount = 0
        tder2 = tder
        pptfl2 = pptflx
        do nk = 1 , ltop
          detlq2(nk) = detlq(nk)
          detic2(nk) = detic(nk)
          udr2(nk) = udr(nk)
          uer2(nk) = uer(nk)
          ddr2(nk) = ddr(nk)
          der2(nk) = der(nk)
          umf2(nk) = umf(nk)
          dmf2(nk) = dmf(nk)
        end do
        fabe = d_one
        fabeold = d_zero
        stab = 0.95D0
        noitr = 0
        istop = 0
        nstep = 0
        aincold = d_zero
        closure_loop: &
        do
          if ( ainc/aincmx > 0.999D0 ) then
            ncount = 0
          else
            ncount = ncount + 1
            !
            !*****************************************************************
            !                                                                *
            !           compute properties for compensational subsidence     *
            !                                                                *
            !*****************************************************************
            !
            ! determine omega value necessary at top and bottom of each layer
            ! to satisfy mass continuity
            !
            dtt = timec
            do nk = 1 , ltop
              domgdp(nk) = -(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
              if ( nk > 1 ) then
                omg(nk) = omg(nk-1) - dp(nk-1)*domgdp(nk-1)
                dtt1 = 0.75D0 * dp(nk-1) / (abs(omg(nk))+1.D-10)
                dtt = min(dtt,dtt1)
              end if
            end do
            do nk = 1 , ltop
              thpa(nk) = thta0(nk)
              qpa(nk) = q0(nk)
              nstep = nint(timec/dtt + d_one)
              dtime = timec/dble(nstep)
              fxm(nk) = omg(nk)*dxsq*regrav
            end do
            !
            ! do an upstream/forward-in-time advection of theta, qv
            !
            do ntc = 1 , nstep
              !
              ! assign theta and q values at the top and bottom of each
              ! layer based sign of omega
              !
              do nk = 1 , ltop
                thfxtop(nk) = d_zero
                thfxbot(nk) = d_zero
                qfxtop(nk) = d_zero
                qfxbot(nk) = d_zero
              end do
              do nk = 2 , ltop
                if ( omg(nk) <= d_zero ) then
                  thfxbot(nk) = -fxm(nk)*thpa(nk-1)
                  qfxbot(nk) = -fxm(nk)*qpa(nk-1)
                  thfxtop(nk-1) = thfxtop(nk-1)-thfxbot(nk)
                  qfxtop(nk-1) = qfxtop(nk-1)-qfxbot(nk)
                else
                  thfxbot(nk) = -fxm(nk)*thpa(nk)
                  qfxbot(nk) = -fxm(nk)*qpa(nk)
                  thfxtop(nk-1) = thfxtop(nk-1)-thfxbot(nk)
                  qfxtop(nk-1) = qfxtop(nk-1)-qfxbot(nk)
                end if
              end do
              !
              ! update the theta and qv values at each level..
              !
              do nk = 1 , ltop
                thpa(nk) = thpa(nk) + &
                     (thfxbot(nk)+udr(nk)*thtau(nk)+ddr(nk)*thtad(nk) + &
                      thfxtop(nk)-(uer(nk)-der(nk))*thta0(nk))*dtime*emsd(nk)
                qpa(nk) = qpa(nk) + &
                     (qfxbot(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk) + &
                      qfxtop(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)
              end do
            end do
            do nk = 1 , ltop
              thtag(nk) = thpa(nk)
              qg(nk) = qpa(nk)
            end do
            !
            ! check to see if mixing ratio dips below zero anywhere;  if so,
            ! borrow moisture from adjacent layers to bring it back up
            ! above zero.
            !
            do nk = 1 , ltop
              if ( qg(nk) < d_zero ) then
                if ( nk == 1 ) then
                  call fatal(__FILE__,__LINE__, &
                          'problem with kf scheme:  qg = 0 at the surface')
                end if
                nk1 = nk + 1
                if ( nk == ltop ) nk1 = klcl
                tma = qg(nk1)*ems(nk1)
                tmb = qg(nk-1)*ems(nk-1)
                tmm = (qg(nk)-1.D-9)*ems(nk)
                bcoeff = -tmm / ((tma*tma) / tmb+tmb)
                acoeff = bcoeff * tma/tmb
                tmb = tmb * (d_one-bcoeff)
                tma = tma * (d_one-acoeff)
                if ( nk == ltop ) then
                  qvdiff = (qg(nk1) - tma*emsd(nk1)) * d_100/qg(nk1)
                  if ( abs(qvdiff) > d_one ) then
                    write(stderr,*) 'Warning from KF Cumulus:'
                    write(stderr,*) 'Cloud base water vapor changes by ', &
                            qvdiff, ' percent'
                    write(stderr,*) 'when moisture is borrowed to prevent&
                           & neg values'
                  end if
                end if
                qg(nk) = 1.D-9
                qg(nk1) = tma*emsd(nk1)
                qg(nk-1) = tmb*emsd(nk-1)
              end if
            end do
            topomg = (udr(ltop) - uer(ltop)) * dp(ltop) * emsd(ltop)
            if ( abs(topomg-omg(ltop)) > 1.D-3 ) then
              write(stderr,*) 'Error from KF Cumulus:'
              write(stderr,*) 'Mass does not balance in scheme : ', &
                              'topomg, omg = ',topomg,omg(ltop)
              istop = 1
              exit closure_loop
            end if
            !
            ! convert theta to t
            !
            ! pay attention
            !
            do nk = 1 , ltop
              exn(nk) = (p00/p0(np,nk))**(c6*(d_one - c7*qg(nk)))
              tg(nk) = thtag(nk) / exn(nk)
              tvg(nk) = tg(nk) * (d_one + b61*qg(nk))
            end do
            !
            !*****************************************************************
            !                                                                *
            !   compute new cloud and change in available buoyant energy.    *
            !                                                                *
            !*****************************************************************
            !
            ! the following computations are similar to that for updraft
            !
            thmix = d_zero
            qmix = d_zero
            pmix = d_zero
            do nk = lc , kpbl
              rocpq = c6*(d_one - c7*qg(nk))
              thmix = thmix + dp(nk)*tg(nk)*(p00/p0(np,nk))**rocpq
              qmix = qmix + dp(nk)*qg(nk)
              pmix = pmix + dp(nk)*p0(np,nk)
            end do
            thmix = thmix/dpthmx
            qmix = qmix/dpthmx
            pmix = pmix/dpthmx
            rocpq = c6*(d_one - c7*qmix)
            tmix = thmix * (pmix/p00)**rocpq
            es = aliq * exp((tmix*bliq - cliq) / (tmix - dliq))
            qs = ep2 * es / (pmix-es)
            !
            ! remove supersaturation for diagnostic purposes, if necessary
            !
            if ( qmix > qs ) then
              rl = xlv0 - xlv1*tmix
              cpm = cpd*(d_one + 0.887D0*qmix)
              dssdt = qs * (cliq - bliq*dliq) / ((tmix-dliq)*(tmix-dliq))
              dq = (qmix-qs) / (d_one + rl*dssdt/cpm)
              tmix = tmix + rl/cpd*dq
              qmix = qmix - dq
              rocpq = c6*(d_one - c7*qmix)
              thmix = tmix * (p00/pmix)**rocpq
              tlcl = tmix
              plcl = pmix
            else
              qmix = max(qmix,d_zero)
              emix = qmix * pmix / (ep2+qmix)
              tlog = log(emix/aliq)
              tdpt = (cliq - dliq*tlog) / (bliq-tlog)
              tlcl = tdpt - (tl1 + tl2 * (tdpt-tzero) - &
                             tl3 * (tmix-tzero)) * (tmix-tdpt)
              tlcl = min(tlcl,tmix)
              cporq = d_one/rocpq
              plcl = p00 * (tlcl/thmix)**cporq
            end if
            tvlcl = tlcl*(d_one + b61*qmix)
            internal_find_plcl: &
            do nk = lc , kl
              klcl = nk
              if ( plcl >= p0(np,nk) ) exit internal_find_plcl
            end do internal_find_plcl
            k = klcl - 1
            dlp = log(plcl/p0(np,k)) / log(p0(np,klcl)/p0(np,k))
            !
            ! estimate environmental temperature and mixing ratio at the lcl
            !
            tenv = tg(k) + (tg(klcl)-tg(k))*dlp
            qenv = qg(k) + (qg(klcl)-qg(k))*dlp
            tven = tenv * (d_one + b61*qenv)
            tvbar = d_half * (tvg(k)+tven)
            zlcl = z0(np,k) + (z0(np,klcl)-z0(np,k))*dlp
            tvavg = d_half * (tven + tg(klcl)*(d_one + b61*qg(klcl)))
            plcl = p0(np,klcl) * exp(egrav / (rdry * tvavg)*(z0(np,klcl)-zlcl))
            theteu(k) = tmix * (p00/pmix)**(c6*(d_one - c7*qmix)) * &
                        exp((c1/tlcl - c2)*qmix*(d_one + c8*qmix))
            es = aliq * exp((tenv*bliq - cliq) / (tenv-dliq))
            qese = ep2 * es / (plcl-es)
            thtesg(k) = tenv*(p00/plcl)**(c6*(d_one - c7*qese)) * &
                        exp((c1/tenv - c2)*qese*(d_one + c8*qese))
            !
            ! compute adjusted abe(abeg).
            !
            abeg = d_zero
            thtudl = theteu(k)
            do nk = k , ltopm1
              nk1 = nk + 1
              es = aliq * exp((tg(nk1)*bliq - cliq) / (tg(nk1)-dliq))
              qese = ep2 * es / (p0(np,nk1)-es)
              thtesg(nk1) = tg(nk1)*(p00/p0(np,nk1))**(c6*(d_one - c7*qese)) * &
                            exp((c1/tg(nk1) - c2)*qese*(d_one + c8*qese))
              if ( nk == k ) then
                dzz = z0(np,klcl) - zlcl
              else
                dzz = dza(nk)
              end if
              be = ((d_two*thtudl) / (thtesg(nk1)+thtesg(nk))-d_one)*dzz
              if ( be > d_zero ) abeg = abeg + be*egrav
            end do
            !
            ! assume at least 90% of cape (abe) is removed by convection during
            ! the period timec
            !
            if ( noitr == 1 ) then
              exit closure_loop
            end if
            dabe = max(abe-abeg,0.1D0*abe)
            fabe = abeg / (abe+1.D-8)
            if ( fabe > d_one ) then
              print *, 'Exit at fabe > d_one'
              cycle pointloop ! No convection
            end if
            doainc = .true.
            if ( ncount /= 1 ) then
              dfda = (fabe-fabeold) / (ainc-aincold)
              if ( dfda > d_zero ) then
                noitr = 1
                ainc = aincold
                doainc = .false.
              end if
            end if
            if ( doainc ) then
              aincold = ainc
              fabeold = fabe
              if ( ainc/aincmx > 0.999D0 .and. fabe > 1.05D0 - stab ) then
                exit closure_loop
              end if
              if ( fabe <= 1.05D0 - stab .and. &
                   fabe >= 0.95D0 - stab) exit closure_loop
              if ( ncount > 10 ) then
                exit closure_loop
              end if
              !
              ! if more than 10% of the original cape remains, increase the
              ! convective mass flux by the factor ainc:
              !
              if ( dabs(fabe) < dlowval ) then
                ainc = ainc*d_half
              else
                ainc = ainc * stab * abe / (dabe+1.D-8)
              end if
            end if
          end if
          ainc = min(aincmx,ainc)
          ! if ainc becomes very small, effects of convection
          ! will be minimal so just ignore it
          if ( ainc < 0.05D0 ) then
              print *, 'Exit at ainc < 0.05D0'
            cycle pointloop
          end if
          tder = tder2 * ainc
          pptflx = pptfl2 * ainc
          do nk = 1 , ltop
            umf(nk) = umf2(nk) * ainc
            dmf(nk) = dmf2(nk) * ainc
            detlq(nk) = detlq2(nk) * ainc
            detic(nk) = detic2(nk) * ainc
            udr(nk) = udr2(nk) * ainc
            uer(nk) = uer2(nk) * ainc
            der(nk) = der2(nk) * ainc
            ddr(nk) = ddr2(nk) * ainc
          end do
          !
          ! go back up for another iteration
          !
        end do closure_loop
        !
        ! clean things up, calculate convective feedback tendencies for this
        ! grid point
        !
        ! compute hydrometeor tendencies as is done for t, qv
        !
        ! frc2 is the fraction of total condensate
        ! generated that goes into precipitiation
        frc2 = pptflx / (cpr*ainc)
        do nk = 1 , ltop
          qlpa(nk) = ql0(nk)
          qipa(nk) = qi0(nk)
          qrpa(nk) = qr0(nk)
          qspa(nk) = qs0(nk)
          rainfb(nk) = pptliq(nk) * ainc*fbfrc*frc2
          snowfb(nk) = pptice(nk) * ainc*fbfrc*frc2
        end do
        do ntc = 1 , nstep
          !
          ! assign hydrometeors concentrations at the top and bottom of each
          ! layer based on the sign of omega
          !
          do nk = 1 , ltop
            qlfxin(nk) = d_zero
            qlfxout(nk) = d_zero
            qifxin(nk) = d_zero
            qifxout(nk) = d_zero
            qrfxin(nk) = d_zero
            qrfxout(nk) = d_zero
            qsfxin(nk) = d_zero
            qsfxout(nk) = d_zero
          end do
          do nk = 2 , ltop
            if ( omg(nk) <= d_zero ) then
              qlfxin(nk) = -fxm(nk)*qlpa(nk-1)
              qifxin(nk) = -fxm(nk)*qipa(nk-1)
              qrfxin(nk) = -fxm(nk)*qrpa(nk-1)
              qsfxin(nk) = -fxm(nk)*qspa(nk-1)
              qlfxout(nk-1) = qlfxout(nk-1)+qlfxin(nk)
              qifxout(nk-1) = qifxout(nk-1)+qifxin(nk)
              qrfxout(nk-1) = qrfxout(nk-1)+qrfxin(nk)
              qsfxout(nk-1) = qsfxout(nk-1)+qsfxin(nk)
            else
              qlfxout(nk) = fxm(nk)*qlpa(nk)
              qifxout(nk) = fxm(nk)*qipa(nk)
              qrfxout(nk) = fxm(nk)*qrpa(nk)
              qsfxout(nk) = fxm(nk)*qspa(nk)
              qlfxin(nk-1) = qlfxin(nk-1)+qlfxout(nk)
              qifxin(nk-1) = qifxin(nk-1)+qifxout(nk)
              qrfxin(nk-1) = qrfxin(nk-1)+qrfxout(nk)
              qsfxin(nk-1) = qsfxin(nk-1)+qsfxout(nk)
            end if
          end do
          !
          ! update the hydrometeor concentration values at each level
          !
          do nk = 1 , ltop
            qlpa(nk) = qlpa(nk) + &
                    (qlfxin(nk)+detlq(nk)-qlfxout(nk))*dtime*emsd(nk)
            qipa(nk) = qipa(nk) + &
                    (qifxin(nk)+detic(nk)-qifxout(nk))*dtime*emsd(nk)
            qrpa(nk) = qrpa(nk) + &
                    (qrfxin(nk)+qlqout(nk)*udr(nk) - &
                     qrfxout(nk)+rainfb(nk))*dtime*emsd(nk)
            qspa(nk) = qspa(nk) + &
                    (qsfxin(nk)+qicout(nk)*udr(nk) - &
                     qsfxout(nk)+snowfb(nk))*dtime*emsd(nk)
          end do
        end do
        do nk = 1 , ltop
          qlg(nk) = qlpa(nk)
          qig(nk) = qipa(nk)
          qrg(nk) = qrpa(nk)
          qsg(nk) = qspa(nk)
        end do
        !
        ! send final parameterized values to output files
        !
        if ( istop == 1 ) then
          write(6,1070) '  p  ','   dp ',' dt k/d ',' dr k/d ', &
                        '   omg  ', ' domgdp ','   umf  ',      &
                        '   uer  ','   udr  ','   dmf  ',       &
                        '   der  ','   ddr  ','   ems  ',       &
                        '    w0  ','  detlq ',' detic '
          do k = ltop , 1 , -1
            dtt = (tg(k)-t0(np,k)) * 86400.0D0 / timec
            rl = xlv0 - xlv1*tg(k)
            dr = -(qg(k)-q0(k)) * rl * 86400.0D0 / (timec*cpd)
            udfrc = udr(k)*timec*emsd(k)
            uefrc = uer(k)*timec*emsd(k)
            ddfrc = ddr(k)*timec*emsd(k)
            defrc = -der(k)*timec*emsd(k)
            write (6,1075) p0(np,k)/100.0D0,            &
                           dp(k)/100.0D0,               &
                           dtt,                         &
                           dr,                          &
                           omg(k),                      &
                           domgdp(k)*1.D4,              &
                           umf(k)/1.D6,                 &
                           uefrc,                       &
                           udfrc,                       &
                           dmf(k)/1.D6,                 &
                           defrc,                       &
                           ddfrc,                       &
                           ems(k)/1.D11,                &
                           w0avg(np,k)*1.D2,            &
                           detlq(k)*timec*emsd(k)*1.D3, &
                           detic(k)*timec*emsd(k)*1.D3
          end do
          write(6,1085) 'k','p','z','t0','tg','dt','tu','td','q0','qg', &
                        'dq','qu','qd','qlg','qig','qrg','qsg','rh0','rhg'
          do k = kx , 1 , -1
            dtt = tg(k) - t0(np,k)
            tuc = tu(k) - tzero
            if ( k < lc .or. k > ltop ) tuc = d_zero
            tdc = tz(k) - tzero
            if ( (k < ldb .or. k > ldt) .and. k /= lfs ) tdc = d_zero
            es = aliq*exp((bliq*tg(k) - cliq) / (tg(k)-dliq))
            qgs = es * ep2 / (p0(np,k)-es)
            rh0 = q0(k)/qes(k)
            rhg = qg(k)/qgs
            write (6,1090) k, &
                           p0(np,k)/100.0D0,       &
                           z0(np,k),               &
                           t0(np,k)-tzero,         &
                           tg(k)-tzero,            &
                           dtt,                    &
                           tuc,                    &
                           tdc,                    &
                           q0(k)*1000.0D0,         &
                           qg(k)*1000.0D0,         &
                           (qg(k)-q0(k))*1000.0D0, &
                           qu(k)*1000.0D0,         &
                           qd(k)*1000.0D0,         &
                           qlg(k)*1000.0D0,        &
                           qig(k)*1000.0D0,        &
                           qrg(k)*1000.0D0,        &
                           qsg(k)*1000.0D0,        &
                           rh0,                    &
                           rhg
          end do
          !
          ! if calculations above show an error in the mass budget, print out a
          ! to be used later for diagnostic purposes, then abort run
          !
          call fatal(__FILE__,__LINE__,'kain-fritsch')
        end if
        cndtnf = (d_one - eqfrc(lfs)) * (qliq(lfs) + qice(lfs))*dmf(lfs)
        !
        !  evaluate moisture budget
        !
        qinit = d_zero
        qfnl = d_zero
        dpt = d_zero
        do nk = 1 , ltop
          dpt = dpt + dp(nk)
          qinit = qinit + q0(nk)*ems(nk)
          qfnl = qfnl + qg(nk)*ems(nk)
          qfnl = qfnl + (qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
        end do
        qfnl = qfnl + pptflx*timec*(d_one-fbfrc)
        err2 = (qfnl-qinit) * 100.0D0/qinit
        if ( abs(err2) > 0.05D0 ) then
          call fatal(__FILE__,__LINE__,'qverr' )
        end if
        relerr = err2*qinit / (pptflx*timec+1.D-10)
        !
        ! feedback to resolvable scale tendencies.
        !
        ! if the advective time period (tadvec) is less than specified minimum
        ! timec, allow feedback to occur only during tadvec
        !
        if ( tadvec < timec ) nic = nint(tadvec/dtsec)
        do k = 1 , kx
          !
          ! if hydrometeors are not allowed, they must be evaporated or
          ! sublimated and fed back as vapor, along with associated changes
          ! in temperature.
          ! note:  this will introduce changes in the convective temperature
          ! and water vapor feedback tendencies and may lead to supersaturated
          ! value of qg
          !
          if ( .not. qi_flag .and. warm_rain ) then
            !
            ! if ice phase is not allowed, melt all frozen hydrometeors
            !
            cpm = cpd * (d_one + 0.887*qg(k))
            tg(k) = tg(k) - (qig(k)+qsg(k)) * rlf/cpm
            dqcdt(np,k) = (qlg(k)+qig(k)-ql0(k)-qi0(k)) / timec
            dqidt(np,k) = d_zero
            dqrdt(np,k) = (qrg(k)+qsg(k)-qr0(k)-qs0(k)) / timec
            dqsdt(np,k) = d_zero
          else if ( .not. qi_flag .and. .not. warm_rain ) then
            !
            ! if ice phase is allowed, but mixed phase is not, melt frozen
            ! hydrome below the melting level, freeze liquid water above
            ! the melting level
            !
            cpm = cpd * (d_one + 0.887*qg(k))
            if ( k <= ml ) then
              tg(k) = tg(k) - (qig(k)+qsg(k)) * rlf/cpm
            else if ( k > ml ) then
              tg(k) = tg(k) + (qlg(k)+qrg(k)) * rlf/cpm
            end if
            dqcdt(np,k) = (qlg(k)+qig(k)-ql0(k)-qi0(k)) / timec
            dqidt(np,k) = d_zero
            dqrdt(np,k) = (qrg(k)+qsg(k)-qr0(k)-qs0(k)) / timec
            dqsdt(np,k) = d_zero
          else if ( qi_flag ) then
            !
            ! if mixed phase hydrometeors are allowed, feed back convective
            ! tendency of hydrometeors directly
            !
            dqcdt(np,k) = (qlg(k)-ql0(k))/timec
            dqidt(np,k) = (qig(k)-qi0(k))/timec
            dqrdt(np,k) = (qrg(k)-qr0(k))/timec
            if ( qs_flag ) then
              dqsdt(np,k) = (qsg(k)-qs0(k))/timec
            else
              dqidt(np,k) = dqidt(np,k)+(qsg(k)-qs0(k))/timec
            end if
          else
            call fatal(__FILE__,__LINE__, &
                    'this combination of imoist, iice not allowed')
          end if
          dtdt(np,k) = (tg(k)-t0(np,k)) / timec
          dqdt(np,k) = (qg(k)-q0(k)) / timec
          ktop(np) = ldt
          kbot(np) = ldb
        end do

        ! raincv is in the unit of mm

        pratec(np) = pptflx*(d_one-fbfrc)/dxsq
        raincv(np) = dtsec*pratec(np)
        rnc = raincv(np)*nic
        exit kfmainloop
      end do kfmainloop
    end do pointloop

 1070 format (16a8)
 1075 format (f8.2,3(f8.2),2(f8.3),f8.2,2f8.3,f8.2,6f8.3)
 1085 format (a3,16a7,2a8)
 1090 format (i3,f7.2,f7.0,10f7.2,4f7.3,2f8.3)

  end subroutine kfpara

  subroutine condload(qliq,qice,wtw,dz,boterm,enterm,rate,qnewlq,     &
                       qnewic,qlqout,qicout)
    implicit none
    !-----------------------------------------------------------------------
    ! 9/18/88 this precipitation fallout scheme is based on the scheme us
    ! by ogura and cho (1973).  liquid water fallout from a parcel is cal-
    ! culated using the equation dq=-rate*q*dt, but to simulate a quasi-
    ! continuous process, and to eliminate a dependency on vertical
    ! resolution this is expressed as q=q*exp(-rate*dz).

    real(rk8) , intent(in) :: dz , boterm , enterm , rate
    real(rk8) , intent(inout) :: qlqout , qicout , wtw , qliq , qice
    real(rk8) , intent(inout) :: qnewlq , qnewic
    real(rk8) :: qtot , qnew , qest , g1 , wavg , conv , ratio3
    real(rk8) :: oldq , ratio4 , dq , pptdrg

    qtot = qliq + qice
    qnew = qnewlq + qnewic
    !
    ! estimate the vertical velocity so that an average vertical velocity c
    ! be calculated to estimate the time required for ascent between model
    ! levels
    !
    qest = d_half * (qtot + qnew)
    g1 = wtw + boterm - enterm - d_two * egrav * dz * qest/1.5D0
    g1 = max(d_zero,g1)
    wavg = d_half * (sqrt(wtw) + sqrt(g1))
    conv = rate * dz/wavg
    !
    ! ratio3 is the fraction of liquid water in fresh condensate, ratio4 is
    ! the fraction of liquid water in the total amount of condensate involv
    ! in the precipitation process - note that only 60% of the fresh conden
    ! sate is is allowed to participate in the conversion process
    !
    ratio3 = qnewlq / (qnew + 1.D-10)
    qtot = qtot + 0.6D0*qnew
    oldq = qtot
    ratio4 = (0.6D0*qnewlq + qliq) / (qtot + 1.D-10)
    qtot = qtot * exp(-conv)
    !
    ! determine the amount of precipitation that falls out of the updraft
    ! parcel at this level
    !
    dq = oldq - qtot
    qlqout = ratio4*dq
    qicout = (d_one - ratio4) * dq
    !
    ! estimate the mean load of condensate on the updraft in the layer, cal
    ! late vertical velocity
    !
    pptdrg = d_half * (oldq + qtot - 0.2D0*qnew)
    wtw = wtw + boterm - enterm - d_two * egrav * dz * pptdrg/1.5D0
    !
    ! determine the new liquid water and ice concentrations including losse
    ! due to precipitation and gains from condensation
    !
    qliq = ratio4*qtot + ratio3 * 0.4D0 * qnew
    qice = (d_one-ratio4) * qtot + (d_one-ratio3) * 0.4D0 * qnew
    qnewlq = d_zero
    qnewic = d_zero
  end subroutine condload

  subroutine dtfrznew(tu,p,thteu,qvap,qliq,qice,ratio2,qnwfrz,rl, &
                      frc1,effq,iflag)
    implicit none
    real(rk8) , intent(in) :: p , effq
    real(rk8) , intent(inout) :: tu , thteu , qvap , qliq , qice
    real(rk8) , intent(inout) :: ratio2 , frc1 , rl , qnwfrz
    integer(ik4) , intent(inout) :: iflag
    real(rk8) :: ccp , qlqfrz , qnew , esliq , esice , rlc , rls , &
            pi , es , rlf , a , b , c , dqvap , dtfrz , tu1 , qvap1
    !
    ! allow glaciation of the updraft to occur as an approximately linear
    ! function of temperature in the temperature range ttfrz to tbfrz
    !
    ! adjust the liquid water concentrations from fresh condensate and tha
    ! brought up from lower levels to an amount that would be present if n
    ! liquid water had frozen thus far. This is necessary because the
    ! expression for temp change is multiplied by the fraction equal to th
    ! parcel temp decrease since the last model level divided by the total
    ! glaciation interval, so that effectively this approximately allows a
    ! amount of liquid water to freeze which is equal to this same fractio
    ! of the liquid water that was present before the glaciation process w
    ! initiateds. Also, to allow thetau to convert approximately linearly
    ! its value with respect to ice, we need to allow a portion of the fre
    ! condensate to contribute to the glaciation process; the fractional
    ! amount that applies to this portion is 1/2 of the fractional amount
    ! frozen of the "old" condensate because this fresh condensate is only
    ! produced gradually over the layer. Note that in terms of the dynami
    ! of the precipitation process, ie. precipitation fallout, this fracti
    ! amnt of fresh condensate has already been included in the ice catego
    !
    qlqfrz = qliq*effq
    qnew = qnwfrz*effq*d_half
    esliq = aliq * exp( (bliq*tu - cliq)/(tu-dliq) )
    esice = aice * exp( (bice*tu - cice)/(tu-dice) )
    rlc = 2.5D6 - 2369.276D0 * (tu-tzero)
    rls = 2833922.0D0 - 259.532*(tu-tzero)
    rlf = rls - rlc
    ccp = cpd * (d_one + c9*qvap)
    !
    ! a = d(es)/dt is that calculated from buck`s (1981) empirical formulas
    ! for saturation vapor pressure
    !
    a = (cice - bice*dice) / ((tu-dice)*(tu-dice))
    b = rls * ep2/p
    c = a * b * esice/ccp
    dqvap = b * (esliq-esice) / (rls + rls*c) - rlf*(qlqfrz+qnew)/(rls + rls/c)
    dtfrz = (rlf * (qlqfrz+qnew) + b*(esliq-esice))/(ccp + a*b*esice)
    tu1 = tu
    qvap1 = qvap
    tu = tu + frc1*dtfrz
    qvap = qvap - frc1*dqvap
    es = qvap * p/(ep2+qvap)
    esliq = aliq * exp((bliq*tu - cliq)/(tu-dliq))
    esice = aice * exp((bice*tu - cice)/(tu-dice))
    ratio2 = (esliq-es) / (esliq-esice)
    !
    ! typically, ratio2 is very close to (ttfrz-tu)/(ttfrz-tbfrz), usually
    ! within 1% (using tu before galciation effects are applied);  if the
    ! initial updraft temp is below tbfrz and ratio2 is still less than 1,
    ! an adjustment to frc1 and ratio2 is introduced so that glaciation
    ! effects are not underestimated; conversely, if ratio2 is greater than
    ! frc1 is adjusted so that glaciation effects are not overestimated
    !
    if ( iflag > 0 .and. ratio2 < d_one ) then
      frc1 = frc1 + (d_one-ratio2)
      tu = tu1 + frc1*dtfrz
      qvap = qvap1 - frc1*dqvap
      ratio2 = d_one
      iflag = 1
    else
      if ( ratio2 > d_one ) then
        frc1 = frc1 - (ratio2 - d_one)
        frc1 = max(d_zero,frc1)
        tu = tu1 + frc1*dtfrz
        qvap = qvap1 - frc1*dqvap
        ratio2 = d_one
        iflag = 1
      end if
    end if
    !
    !  calculate a hybrid value of thetau, assuming that the latent heat of
    !  vaporization/sublimation can be estimated using the same weighting
    !  function as that used to calculate saturation vapor pressure, calcu-
    !  late new liquid water and ice concentrations
    !
    rlc = xlv0 - xlv1*tu
    rls = xls0 - xls1*tu
    rl = ratio2 * rls + (d_one-ratio2)*rlc
    pi = (p00/p)**(c6*(d_one - c7*qvap))
    thteu = tu * pi * exp(rl * qvap * c5/tu * (d_one + c8*qvap))
    if ( iflag == 1 ) then
      qice = qice + frc1*dqvap + qliq
      qliq = d_zero
    else
      qice = qice + frc1*(dqvap+qlqfrz)
      qliq = qliq - frc1*qlqfrz
    end if
    qnwfrz = d_zero
  end subroutine dtfrznew

  !-----------------------------------------------------------------------
  !  This subroutine integrates the area under the curve in the gaussian
  !  distribution. The numerical approximation to the integral is taken f
  !  handbook of mathematical functions with formulas, graphs and mathema
  !  tables  ed. by abramowitz and stegun, nat l bureau of standards appli
  !  mathematics series.  June, 1964., May, 1968.
  !       Jack Kain
  !       7/6/89
  !*****    gaussian type mixing profile   ******************************
  subroutine prof5(eq,ee,ud)
    implicit none
    real(rk8) , intent(in) :: eq
    real(rk8) , intent(inout) :: ee , ud
    real(rk8) , parameter :: sqrt2p = 2.506628D0
    real(rk8) , parameter :: a1 = 0.4361836D0
    real(rk8) , parameter :: a2 = -0.1201676D0
    real(rk8) , parameter :: a3 = 0.9372980D0
    real(rk8) , parameter :: p = 0.33267D0
    real(rk8) , parameter :: sigma = 0.166666667D0
    real(rk8) , parameter :: fe = 0.202765151D0
    real(rk8) :: x , y , ey , e45 , t1 , t2 , c1 , c2

    x = (eq - 0.5D0)/sigma
    y = 6.0D0*eq - 3.0D0
    ey = exp(y*y/(-2.0D0))
    e45 = exp(-4.5D0)
    t2 = d_one/(d_one + p*abs(y))
    t1 = 0.500498D0
    c1 = a1*t1 + a2*t1*t1 + a3*t1*t1*t1
    c2 = a1*t2 + a2*t2*t2 + a3*t2*t2*t2
    if ( y >= d_zero ) then
      ee = sigma * (0.5D0*(sqrt2p - e45*c1 - ey*c2) + sigma*(e45-ey)) - &
           e45 * eq * eq/2.0D0
      ud = sigma * (0.5D0*(ey*c2 - e45*c1) + sigma*(e45-ey)) - &
           e45 * (0.5D0 + eq*eq/2.0D0 - eq)
    else
      ee = sigma * (0.5D0*(ey*c2 - e45*c1) + sigma*(e45-ey)) - e45*eq*eq/2.0D0
      ud = sigma * (0.5D0*(sqrt2p - e45*c1 - ey*c2) + sigma*(e45-ey)) - &
                e45*(0.5D0 + eq*eq/2.0D0 - eq)
    end if
    ee = ee/fe
    ud = ud/fe
  end subroutine prof5

  subroutine tpmix(p,thtu,tu,qu,qliq,qice,qnewlq,qnewic,ratio2,rl)
    implicit none
    real(rk8) , intent(in) :: p , thtu , ratio2 , rl
    real(rk8) , intent(inout) :: qu , qliq , qice , tu , qnewlq , qnewic
    real(rk8) :: es , qs , pi , thtgs , f0 , t1 , t0
    real(rk8) :: esliq , esice , f1 , dt , qnew
    real(rk8) :: dq , qtot , dqice , dqliq , rll , ccp
    real(rk8) , parameter :: diffconv = 0.01D0
    integer(ik4) :: itcnt
    !
    ! This subroutine iteratively extracts wet-bulb temperature from equiv
    ! potential temperature, then checks to see if sufficient moisture is
    ! available to achieve saturation. If not, temperature is adjusted
    ! accordingly, if so, the residual liquid water/ice concentration is
    ! determined
    !
    ! Iterate to find wet bulb temperature as a function of equivalent pot
    ! temp and prs, assuming saturation vapor pressure. Ratio2 is the deg
    ! of glaciation
    !
    if ( ratio2 < 1.D-6 ) then
      es = aliq*exp((bliq*tu - cliq) / (tu-dliq))
      qs = ep2 * es / (p-es)
      pi = (p00/p)**(c6*(d_one - c7*qs))
      thtgs = tu * pi * exp((c1/tu - c2)*qs*(d_one + c8*qs))
    else if ( abs(ratio2-d_one) < 1.D-6 ) then
      es = aice * exp((bice*tu - cice)/(tu-dice))
      qs = ep2 * es / (p-es)
      pi = (p00/p)**(c6*(d_one - c7*qs))
      thtgs = tu * pi * exp((c3/tu - c4)*qs*(d_one+c8*qs))
    else
      esliq = aliq*exp((bliq*tu - cliq) / (tu-dliq))
      esice = aice*exp((bice*tu - cice) / (tu-dice))
      es = (d_one-ratio2)*esliq + ratio2*esice
      qs = ep2 * es / (p-es)
      pi = (p00/p)**(c6*(d_one-c7*qs))
      thtgs = tu * pi * exp(rl * qs * c5/tu * (d_one+c8*qs))
    end if
    f0 = thtgs - thtu
    t1 = tu - d_half*f0
    t0 = tu
    itcnt = 0
    do
      if ( ratio2 < 1.D-6 ) then
        es = aliq * exp((bliq*t1 - cliq) / (t1-dliq))
        qs = ep2 * es / (p-es)
        pi = (p00/p)**(c6*(d_one-c7*qs))
        thtgs = t1 * pi * exp((c1/t1 - c2)*qs*(d_one+c8*qs))
      else if ( abs(ratio2-d_one) < 1.D-6 ) then
        es = aice * exp((bice*t1 - cice) / (t1-dice))
        qs = ep2 * es / (p-es)
        pi = (p00/p)**(c6*(d_one-c7*qs))
        thtgs = t1 * pi * exp((c3/t1-c4)*qs*(d_one+c8*qs))
      else
        esliq = aliq * exp((bliq*t1 - cliq) / (t1-dliq))
        esice = aice * exp((bice*t1 - cice) / (t1-dice))
        es = (d_one-ratio2) * esliq + ratio2*esice
        qs = ep2 * es / (p-es)
        pi = (p00/p)**(c6*(d_one-c7*qs))
        thtgs = t1 * pi * exp(rl * qs * c5/t1 * (d_one+c8*qs))
      end if
      f1 = thtgs - thtu
      if ( abs(f1) < diffconv ) exit
      itcnt = itcnt + 1
      if ( itcnt > maxiter ) exit
      dt = f1 * (t1-t0) / (f1-f0)
      t0 = t1
      f0 = f1
      t1 = t1 - dt
    end do
    !
    !If the parcel is supersaturated, calculate concentration of fresh
    !condensate
    !
    do
      if ( qs <= qu ) then
        qnew = qu - qs
        qu = qs
        exit
      end if
      !
      ! if the parcel is subsaturated, temperature and mixing ratio must be
      ! adjusted. If liquid water or ice is present, it is allowed to evapo
      ! sublimate.
      !
      qnew = d_zero
      dq = qs - qu
      qtot = qliq + qice
      !
      ! If there is enough liquid or ice to saturate the parcel, temp stays
      ! wet bulb value, vapor mixing ratio is at saturated level, and the mi
      ! ratios of liquid and ice are adjusted to make up the original satura
      ! deficit. Otherwise, any available liq or ice vaporizes and appropr
      ! adjustments to parcel temp; vapor, liquid, and ice mixing ratios are
      !
      ! Note that the liq and ice may be present in proportions slightly dif
      ! than suggested by the value of ratio2. Check to make sure that liq
      ! ice concentrations are not reduced to below zero when evaporation/
      ! sublimation occurs
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
        ccp = cpd * (d_one + c9*qu)
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
    if ( itcnt > maxiter) then
      write(stderr,*) '***** number of iterations in tpmix =', itcnt
    end if
  end subroutine tpmix

  subroutine envirtht(p1,t1,q1,tht1,r1,rl)
    implicit none
    real(rk8) , intent(in) :: p1 , t1 , q1 , r1 , rl
    real(rk8) ,  intent(inout) :: tht1
    real(rk8) :: ee , tlog , tdpt , tsat , tht , tfpt , tlogic
    real(rk8) :: tsatlq , tsatic

    !  calculate environmental equivalent potential temperature

    if ( r1 < 1.D-6 ) then
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aliq)
      tdpt = (cliq - dliq*tlog) / (bliq-tlog)
      tsat = tdpt - (tl1 + tl2*(tdpt-tzero) - tl3*(t1-tzero)) * (t1-tdpt)
      tht = t1 * (p00/p1)**(c6*(d_one - c7*q1))
      tht1 = tht * exp((c1/tsat-c2)*q1*(d_one + c8*q1))
    else if ( abs(r1-d_one) < 1.D-6 ) then
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aice)
      tfpt = (cice - dice*tlog) / (bice - tlog)
      tht = t1*(p00/p1)**(c6*(d_one - c7*q1))
      tsat = tfpt - (ti1 + ti2*(tfpt-tzero) - ti3*(t1-tzero)) * (t1-tfpt)
      tht1 = tht * exp((c3/tsat-c4)*q1*(d_one + c8*q1))
    else
      ee = q1 * p1 / (ep2+q1)
      tlog = log(ee/aliq)
      tdpt = (cliq-dliq*tlog) / (bliq-tlog)
      tlogic = log(ee/aice)
      tfpt = (cice - dice*tlogic) / (bice-tlogic)
      tht = t1 * (p00/p1)**(c6*(d_one - c7*q1))
      tsatlq = tdpt - (tl1 + tl2*(tdpt-tzero) - tl3*(t1-tzero)) * (t1-tdpt)
      tsatic = tfpt - (ti1 + ti2*(tfpt-tzero) - ti3*(t1-tzero)) * (t1-tfpt)
      tsat = r1 * tsatic + (d_one-r1)*tsatlq
      tht1 = tht * exp(rl*q1*c5/tsat*(d_one+c8*q1))
    end if
  end subroutine envirtht

  ! This subroutine iteratively extracts temperature from equivalent
  ! potential temp.  it is designed for use with downdraft calculations.
  ! if relative humidity is specified to be less than 100%, parcel
  ! temp, specific humidity, and liquid water content are iteratively
  ! calculated.
  real(rk8) function tpdd(p,thted,tgs,rs,rd,rh)
    implicit none
    real(rk8) , intent(in) :: p , thted , tgs , rd , rh
    real(rk8) , intent(inout) :: rs
    real(rk8) :: es , pi , thtgs , f0 , t1 , t0 , ccp , f1 , dt
    real(rk8) :: rl , dssdt , t1rh , rsrh
    real(rk8) , parameter :: diffconv = 0.05D0
    integer(ik4) :: itcnt

    es = aliq * exp( (bliq*tgs - cliq) / (tgs - dliq) )
    rs = ep2 * es / (p - es)
    pi = (p00/p)**(c6*(d_one - c7*rs))
    thtgs = tgs * pi * exp( (c1/tgs - c2) * rs * (d_one + c8*rs) )
    f0 = thtgs - thted
    t1 = tgs - d_half*f0
    t0 = tgs
    ccp = cpd
    !
    ! Iterate to find wet-bulb temperature
    !
    itcnt = 0
    iterate_loop: &
    do
      es = aliq * exp((bliq*t1 - cliq)/(t1 - dliq))
      rs = ep2 * es/(p - es)
      pi = (p00/p)**(c6*(d_one - c7*rs))
      thtgs = t1 * pi * exp((c1/t1 - c2)*rs*(d_one + c8*rs))
      f1 = thtgs - thted
      if ( abs(f1) < diffconv ) exit iterate_loop
      itcnt = itcnt + 1
      if (itcnt > maxiter ) exit iterate_loop
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
    if ( (rh - d_one) < dlowval ) then
      tpdd = t1
    else
      dssdt = (cliq-bliq * dliq) / ((t1-dliq)*(t1-dliq))
      dt = rl * rs * (d_one-rh) / (ccp + rl*rh*rs*dssdt)
      t1rh = t1 + dt
      es = rh * aliq * exp( (bliq*t1rh - cliq) / (t1rh-dliq) )
      rsrh = ep2 * es/(p-es)
      !
      ! Check to see if mixing ratio at specified rh is less than actual
      ! mixing ratio. If so, adjust to give zero evaporation
      !
      if ( rsrh < rd ) then
        rsrh = rd
        t1rh = t1 + (rs-rsrh) * rl/ccp
      end if
      t1 = t1rh
      rs = rsrh
      tpdd = t1
    end if
    if ( itcnt < maxiter ) then
      write(stderr,*) '***** number of iterations in tpdd = ', itcnt
    end if
  end function tpdd

end module mod_cu_kf
