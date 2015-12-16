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
  use mod_runparams , only : dx , dxsq , ipptls , ibltyp , dt
  use mod_runparams , only : iqv , iqr , iqi , iqs , iqc
  use mod_runparams , only : kf_entrate , kf_min_pef , kf_max_pef
  use mod_runparams , only : kf_dpp , kf_min_dtcape , kf_max_dtcape
  use mod_runparams , only : kf_tkemax
  use mod_runparams , only : ichem , clfrcv
  use mod_service

  implicit none

  private
  !
  !  RegCM KF code
  !
  ! The code below is adapted from kfeta code from WRF 3.6.1 codebase.
  ! It is missing the kf_trigger == 2 option, and we use by default the
  ! kf_trigger == 1 option.
  !
  ! Trieste , August 2014
  !
  !             Graziano Giuliani
  !
  integer , parameter :: kf_trigger = 1
  !
  public :: allocate_mod_cu_kf , kfdrv , kf_lutab , kfwavg
  !
  !  V3.3: A new trigger function is added based Ma and Tan (2009):
  !   Ma, L.-M. and Z.-M. Tan, 2009: Improving the behavior of
  !   the cumulus parameterization for tropical cyclone prediction:
  !   Convection trigger. Atmospheric Research, 92, 190 - 211.
  !
  !  WRF v3.5 with diagnosed deep and shallow KF cloud fraction using
  !  CAM3-CAM5 methodology, along with captured liquid and ice condensates.
  !    JAH & KA (U.S. EPA) -- May 2013
  !
  integer(ik4) , parameter :: kfnt = 250
  integer(ik4) , parameter :: kfnp = 220
  integer(ik4) , parameter :: kfna = 200
  real(rk8) , dimension(kfnt,kfnp) , private , save :: ttab , qstab
  real(rk8) , dimension(kfnp) , private , save :: the0k
  real(rk8) , dimension(kfna) , private , save :: alu
  real(rk8) , private , save :: rdpr , rdthk , plutop

  integer(ik4) :: nipoi

  integer(ik4) , dimension(:) , pointer :: imap , jmap
  real(rk8) , dimension(:,:) , pointer :: u0 , v0 , z0 , t0 , qv0 , p0
  real(rk8) , dimension(:,:) , pointer :: rho , dzq , w0avg , tke
  real(rk8) , dimension(:) , pointer :: raincv , pratec
  real(rk8) , dimension(:,:) , pointer :: qc_kf , qi_kf
  integer(ik4) , dimension(:) , pointer :: ktop , kbot
  real(rk8) , dimension(:,:) , pointer :: cldfra_dp_kf , cldfra_sh_kf
  real(rk8) , dimension(:,:) , pointer :: pptliq , pptice

  real(rk8) , dimension(:,:) , pointer :: dqdt , dtdt , dqcdt
  real(rk8) , dimension(:,:) , pointer :: tpart_h , tpart_v

  real(rk8) , dimension(:,:,:) , pointer :: kfwavg

  ! IPPTLS == 2
  real(rk8) , dimension(:,:) , pointer :: dqidt , dqrdt , dqsdt

  real(rk8) , parameter :: t00 = tzero
  real(rk8) , parameter :: p00 = 1.0D5

  real(rk8) , parameter :: astrt = 1.0D-3
  real(rk8) , parameter :: aincb = 0.050D0

  real(rk8) , parameter :: c1 = 3374.6525D0
  real(rk8) , parameter :: c2 = 2.5403D0
  real(rk8) , parameter :: c4 = 0.810D0
  real(rk8) , parameter :: dpmin = 5.0D3
  real(rk8) , parameter :: ttfrz = tzero - 5.0D0
  real(rk8) , parameter :: tbfrz = tzero - 25.0D0

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
    call getmem2d(z0,1,nipoi,1,kz,'mod_cu_kf:z0')
    call getmem2d(qv0,1,nipoi,1,kz,'mod_cu_kf:qv0')
    call getmem2d(p0,1,nipoi,1,kz,'mod_cu_kf:p0')
    call getmem2d(rho,1,nipoi,1,kz,'mod_cu_kf:rho')
    call getmem2d(tke,1,nipoi,1,kz,'mod_cu_kf:tke')
    call getmem2d(dzq,1,nipoi,1,kz,'mod_cu_kf:dzq')
    call getmem2d(w0avg,1,nipoi,1,kz,'mod_cu_kf:w0avg')
    call getmem2d(dqdt,1,nipoi,1,kz,'mod_cu_kf:dqdt')
    call getmem2d(dqidt,1,nipoi,1,kz,'mod_cu_kf:dqidt')
    call getmem2d(dqcdt,1,nipoi,1,kz,'mod_cu_kf:dqcdt')
    call getmem2d(dqrdt,1,nipoi,1,kz,'mod_cu_kf:dqrdt')
    call getmem2d(dqsdt,1,nipoi,1,kz,'mod_cu_kf:dqsdt')
    call getmem2d(dtdt,1,nipoi,1,kz,'mod_cu_kf:dtdt')
    call getmem2d(cldfra_dp_kf,1,nipoi,1,kz,'mod_cu_kf:cldfra_dp_kf')
    call getmem2d(pptliq,1,nipoi,1,kz,'mod_cu_kf:pptliq')
    call getmem2d(pptice,1,nipoi,1,kz,'mod_cu_kf:pptice')
    call getmem2d(cldfra_sh_kf,1,nipoi,1,kz,'mod_cu_kf:cldfra_sh_kf')
    call getmem2d(qc_kf,1,nipoi,1,kz,'mod_cu_kf:qc_kf')
    call getmem2d(qi_kf,1,nipoi,1,kz,'mod_cu_kf:qi_kf')
    call getmem1d(raincv,1,nipoi,'mod_cu_kf:raincv')
    call getmem1d(pratec,1,nipoi,'mod_cu_kf:pratec')
    call getmem1d(ktop,1,nipoi,'mod_cu_kf:ktop')
    call getmem1d(kbot,1,nipoi,'mod_cu_kf:kbot')
    if ( kf_trigger == 2 ) then
      call getmem2d(tpart_v,1,nipoi,1,kz,'mod_cu_kf:tpart_v')
      call getmem2d(tpart_h,1,nipoi,1,kz,'mod_cu_kf:tpart_h')
    end if
    call getmem3d(kfwavg,jci1,jci2,ici1,ici2,1,kz,'mod_cu_kf:kfwavg')
  end subroutine allocate_mod_cu_kf

  subroutine kfdrv(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    integer :: i , j ,  k , kk , np
    real(rk8) :: w0
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'kfdrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nipoi == 0 ) return

    do k = 1 , kz
      kk = kzp1 - k
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
        w0 = d_half * (m2c%was(j,i,kk)+m2c%was(j,i,kk+1))
        ! Average over four timesteps.
        w0avg(np,k) = (kfwavg(j,i,kk)*3.0D0+w0)*d_rfour
      end do
    end do

    if ( ibltyp == 2 ) then
      do k = 1 , kz
        kk = kzp1 - k
        do np = 1 , nipoi
          i = imap(np)
          j = jmap(np)
          tke(np,k) = d_half * (m2c%tkeas(j,i,kk)+m2c%tkeas(j,i,kk+1))
        end do
      end do
    else
      tke(:,:) = kf_tkemax
    end if

    dtdt(:,:) = d_zero
    dqdt(:,:) = d_zero
    dqidt(:,:) = d_zero
    dqcdt(:,:) = d_zero
    dqrdt(:,:) = d_zero
    dqsdt(:,:) = d_zero
    raincv(:) = d_zero
    pratec(:) = d_zero
    pptliq(:,:) = d_zero
    pptice(:,:) = d_zero
    ktop(:) = 0
    kbot(:) = 0

    if ( kf_trigger == 2 ) then
      call fatal(__FILE__,__LINE__,'Not implemented kf_trigger == 2')
    end if

    if ( ipptls == 2 ) then
      call kfpara(1,kz,1,nipoi,.true.,.true.,.false.)
    else
      call kfpara(1,kz,1,nipoi,.false.,.false.,.true.)
    end if

    do k = 1 , kz
      kk = kz - k + 1
      do np = 1 , nipoi
        i = imap(np)
        j = jmap(np)
        cu_tten(j,i,kk) = dtdt(np,k)
        cu_qten(j,i,kk,iqv) = dqdt(np,k)
        cu_qten(j,i,kk,iqc) = dqcdt(np,k)
        cu_cldfrc(j,i,k) = max(cldfra_sh_kf(np,k),cldfra_dp_kf(np,k))
        kfwavg(j,i,kk) = w0avg(np,k)
      end do
    end do

    if ( ipptls == 2 ) then
      do k = 1 , kz
        kk = kz - k + 1
        do np = 1 , nipoi
          i = imap(np)
          j = jmap(np)
          cu_qten(j,i,kk,iqr) = dqrdt(np,k)
          cu_qten(j,i,kk,iqi) = dqidt(np,k)
          cu_qten(j,i,kk,iqs) = dqsdt(np,k)
        end do
      end do
    end if

    ! Build for chemistry 3d table of constant precipitation rate
    ! from the surface to the top of the convection
    if ( ichem == 1 ) then
      do k = 1 , kz
        kk = kz - k + 1
        do np = 1 , nipoi
          i = imap(np)
          j = jmap(np)
          cu_convpr(j,i,kk) = pptliq(np,k) + pptice(np,k)
        end do
      end do
    end if
    do np = 1 , nipoi
      i = imap(np)
      j = jmap(np)
      if ( pratec(np) > dlowval ) then
        cu_ktop(j,i) = kz-ktop(np)+1
        cu_kbot(j,i) = kz-kbot(np)+1
        cu_prate(j,i)= cu_prate(j,i) + pratec(np)
        total_precip_points = total_precip_points + 1
      end if
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine kfdrv
  !
  ! The KF scheme that is currently used in experimental runs of EMCs
  ! Eta model  jsk 8/00
  !
  subroutine kfpara(kts,kte,its,ite,f_qi,f_qs,warm_rain)
    implicit none
    integer(ik4) , intent(in) :: kts , kte , its , ite
    logical , intent(in) :: f_qi , f_qs
    logical , intent(in) :: warm_rain

    real(rk8) , dimension(kts:kte) :: q0 , tv0 , tu , tvu , qu , tz ,  &
            tvd , qd , qes , tg , tvg , qg , wu , wd , ems , emsd ,  &
            umf , uer , udr , dmf , der , ddr , umf2 , uer2 , udr2 , &
            dmf2 , der2 , ddr2 , dza , thta0 , thetee , thtau ,      &
            theteu , thtad , theted , qliq , qice , qlqout ,         &
            qicout , detlq , detic , detlq2 ,      &
            detic2 , ratio2
    real(rk8) , dimension(kts:kte) :: domgdp , exn , tvqu , dp , rh ,  &
            eqfrc , wspd , qdt , fxm , thtag , thpa , thfxout ,      &
            thfxin , qpa , qfxout , qfxin , qlpa , qlfxin , qlfxout ,&
            qipa , qifxin , qifxout , qrpa , qrfxin , qrfxout ,      &
            qspa , qsfxin , qsfxout , ql0 , qlg , qi0 , qig , qr0 ,  &
            qrg , qs0 , qsg

    real(rk8) , dimension(kts:kte+1) :: omg
    real(rk8) , dimension(kts:kte) :: rainfb , snowfb
    real(rk8) , dimension(kts:kte) :: cldhgt , qsd , dilfrc , ddilfrc , &
            tgu , qgu , thteeg

    real(rk8) :: fbfrc , p300 , dpthmx , qmix , zmix , pmix , tmix , emix , &
            tlog , tdpt , tlcl , tvlcl , plcl , es , dlp , tenv , qenv ,  &
            tven , zlcl , wkl , trppt , dtlcl , gdt , wlcl , wtw ,        &
            rholcl , au0 , vmflcl , upold , upnew , abe , wklcl , ttemp , &
            frc1 , qnewic , rl , be , boterm , enterm , dzz , rei , ee2 , &
            ud2 , ttmp , f1 , f2 , thttmp , qtmp , tmpliq , tmpice ,      &
            tu95 , tu10 , ee1 , ud1 , dptt , qnewlq , dumfdp , vconv ,    &
            timec , shsign , vws , pef , cbh , rcbh , pefcbh , peff ,     &
            peff2 , tder , tadvec , dpdd , rdd , a1 , dssdt ,             &
            dtmp , t1rh , qsrh , pptflx , cpr , cndtnf , updinc ,         &
            aincm2 , ddinc , aincmx , aincm1 , ainc , tder2 , pptfl2 ,    &
            fabe , stab , dtt , dtt1 , dtime , tma , tmb , tmm , bcoeff , &
            acoeff , topomg , cpm , dq , abeg , dabe , dfda ,             &
            frc2 , dr , udfrc , tuc , qgs , rh0 , rhg , qinit , qfnl ,    &
            err2 , relerr , fabeold , aincold , uefrc , ddfrc , tdc ,     &
            defrc , rhbar , dmffrc , dilbe
    real(rk8) :: tp , avalue , aintrp , qfrz , qss , pptmlt , dtmelt , &
            rhh , evac , binc
    integer(ik4) :: indlu , nu , nuchm , nnn , klfs
    real(rk8) :: chmin , pm15 , chmax , dtrh , rad , dppp
    real(rk8) :: tvdiff , dttot , absomg , absomgtc , frdp
    real(rk8) :: xcldfra , umf_new
    integer(ik4) :: kx , k , kl
    integer(ik4) :: ncheck
    integer(ik4) , dimension(kts:kte) :: kcheck
    integer(ik4) :: istop , ml , l5 , kmix , low , lc , llfc , nlayrs , &
            nk , kpbl , klcl , lcl , let , iflag , nk1 , ltop , nj ,  &
            ltop1 , ltopm1 , kstart , lfs , nd , nic , ldb , ldt ,    &
            nd1 , ndk , lmax , ncount , noitr , nstep , ntc , ishall , np
    logical :: iprnt
    real(rk8) :: qslcl , rhlcl , dqssdt    !jfb
    integer , parameter :: maxiter = 100

    kl = kte
    kx = kte

    modelpoints: &
    do np = its , ite
      iprnt = .false.
      !                                                    ! PPT FB MODS
      ! OPTION TO FEED CONVECTIVELY GENERATED RAINWATER    ! PPT FB MODS
      ! INTO GRID-RESOLVED RAINWATER (OR SNOW/GRAUPEL)     ! PPT FB MODS
      ! FIELD.  "FBFRC" IS THE FRACTION OF AVAILABLE       ! PPT FB MODS
      ! PRECIPITATION TO BE FED BACK (0.0 - 1.0)...        ! PPT FB MODS
      fbfrc = d_zero                                       ! PPT FB MODS
      ! mods to allow shallow convection...
      ishall = 0
      llfc = 1
      l5 = 1
      p300 = p0(np,1) - 30000.0D0
      !
      ! Pressure perturbation term is only defined at mid-point of
      ! vertical layers. Since total pressure is needed at the top and
      ! bottom of layers below, do an interpolation.
      !
      ! Input a vertical sounding. Note that model layers are numbered
      ! from bottom-up in the kf scheme.
      !
      ml = 0
      do k = 1 , kx
        !
        ! Saturation vapor pressure (ES) is calculated following Buck (1981)
        ! If q0 is above saturation value, reduce it to saturation level.
        !
        es = aliq * exp((bliq*t0(np,k)-cliq)/(t0(np,k)-dliq))
        qes(k) = ep2 * es/(p0(np,k)-es)
        q0(k) = min(qes(k),qv0(np,k))
        q0(k) = max(0.000001D0,q0(k))
        ql0(k) = d_zero
        qi0(k) = d_zero
        qr0(k) = d_zero
        qs0(k) = d_zero
        rh(k) = q0(k) / qes(k)
        dilfrc(k) = d_one
        tv0(k) = t0(np,k) * (d_one + ep1*q0(k))
        ! dp is the pressure interval between full sigma levels
        dp(k) = rho(np,k)*egrav*dzq(np,k)
        ! If Turbulent Kinetic Energy (TKE) is available from turbulent
        ! mixing scheme use it for shallow convection.
        ! For now, assume it is not available
        cldhgt(k) = d_zero
        if ( p0(np,k) >= d_half*p0(np,1) ) l5 = k
        if ( p0(np,k) >= p300) llfc = k
      end do
      ! dzq is dz between sigma surfaces, dza is dz between model half level
      do k = 2 , kl
        dza(k-1) = z0(np,k)-z0(np,k-1)
      end do
      dza(kl) = d_zero
      !
      !
      ! To save time, specify a pressure interval to move up in sequential
      ! check of different ~50 mb deep groups of adjacent model layers in
      ! the process of identifying updraft source layer (USL).  Note that
      ! this search is terminated as soon as a buoyant parcel is found and
      ! this parcel can produce a cloud greater than specifed minimum depth
      ! (CHMIN)...For now, set interval at 15 mb...
      !
      kcheck(:) = 0
      kcheck(1) = 1
      ncheck = 1
      pm15 = p0(np,1) - 15.0D2
      do k = 2 , llfc
        if ( p0(np,k) < pm15 ) then
          ncheck = ncheck+1
          kcheck(ncheck) = k
          pm15 = pm15 - 15.0D2
        end if
      end do

      nu = 0
      nuchm = 0
      usl: &
      do
        nu = nu+1
        if ( nu > ncheck ) then
          if ( ishall == 1 ) then
            chmax = d_zero
            do nk = 1 , ncheck
              nnn = kcheck(nk)
              if ( cldhgt(nnn) > chmax ) then
                nuchm = nk
                chmax = cldhgt(nnn)
              end if
            end do
            nu = nuchm-1
            fbfrc = d_one
            cycle usl
          else
            cycle modelpoints
          end if
        end if
        kmix = kcheck(nu)
        low = kmix
        lc = low
        !
        ! Assume that in order to support a deep updraft you need a layer of
        ! unstable air at least 50 mb deep. To approximate this, isolate a
        ! group of adjacent individual model layers, with the base at level
        ! lc, such that the combined depth of these layers is at least 50 mb.
        !
        nlayrs = 0
        dpthmx = d_zero
        nk = lc-1
        if ( nk+1 < kts ) then
          write(stderr,*) 'WOULD GO OFF BOTTOM: KF_ETA_PARA NK', NK
          write(stderr,*) 'AT I = ',imap(np), ', J = ', jmap(np)
          call fatal(__FILE__,__LINE__,'KF FATAL ERROR')
        else
          calcdpth: &
          do
            nk = nk + 1
            if ( nk > kte ) then
              write(stderr,*) 'WOULD GO OFF TOP: KF_ETA_PARA.'
              write(stderr,*) 'AT I = ',imap(np), ', J = ', jmap(np)
              write(stderr,*) 'DPTHMX = ', dpthmx
              write(stderr,*) 'DPMIN  = ', dpmin
              call fatal(__FILE__,__LINE__,'KF FATAL ERROR')
            end if
            dpthmx = dpthmx+dp(nk)
            nlayrs = nlayrs+1
            if ( dpthmx > dpmin ) then
              exit calcdpth
            end if
          end do calcdpth
        end if
        if ( dpthmx < dpmin ) then
          cycle modelpoints
        end if
        kpbl = lc+nlayrs-1
        !
        ! For computational simplicity without much loss in accuracy,
        ! mix temperature instead of theta for evaluating convective
        ! initiation (triggering) potential.
        !
        tmix = d_zero
        qmix = d_zero
        zmix = d_zero
        pmix = d_zero
        !
        ! Find the thermodynamic characteristics of the layer by
        ! mass-weighting the characteristics of the individual model
        ! layers.
        !
        do nk = lc , kpbl
          tmix = tmix + dp(nk)*t0(np,nk)
          qmix = qmix + dp(nk)*q0(nk)
          zmix = zmix + dp(nk)*z0(np,nk)
          pmix = pmix + dp(nk)*p0(np,nk)
        end do
        tmix = tmix/dpthmx
        qmix = qmix/dpthmx
        zmix = zmix/dpthmx
        pmix = pmix/dpthmx
        emix = qmix*pmix/(ep2+qmix)
        !
        ! Find the temperature of the mixture at its lcl.
        ! Calculate dewpoint using lookup table.
        !
        ainc = aincb
        a1 = emix/aliq
        tp = (a1-astrt)/ainc
        indlu = max(1, min(kfna-1,int(tp)+1))
        avalue = (indlu-1)*ainc + astrt
        aintrp = (a1-avalue)/ainc
        tlog = aintrp*alu(indlu+1) + (d_one-aintrp)*alu(indlu)
        tdpt = (cliq-dliq*tlog) / (bliq-tlog)
        tlcl = tdpt - (0.212D0 + 1.571D-3*(tdpt-t00) - &
                       4.36D-4*(tmix-t00))*(tmix-tdpt)
        tlcl = min(tlcl,tmix)
        tvlcl = tlcl*(d_one + ep1*qmix)
        zlcl = zmix + (tlcl-tmix)/gdry
        findklcl1: &
        do nk = lc , kl
          klcl = nk
          if ( zlcl <= z0(np,nk) ) exit findklcl1
        end do findklcl1
        if ( zlcl > z0(np,kl) ) cycle modelpoints

        k = max(1,klcl-1)
        ! Calculate DLP using Z instead of log(P)
        dlp = (zlcl-z0(np,k)) / (z0(np,klcl)-z0(np,k))
        !
        ! Estimate environmental temperature and mixing ratio at the lcl.
        !
        tenv = t0(np,k) + (t0(np,klcl)-t0(np,k))*dlp
        qenv = q0(k) + (q0(klcl)-q0(k))*dlp
        tven = tenv*(d_one + ep1*qenv)
        !
        ! Check to see if cloud is buoyant using fritsch-chappell trigger
        ! function described in kain and fritsch (1992). w0 is an
        ! aproximate value for the running-mean grid-scale vertical
        ! velocity, which gives smoother fields of convective initiation
        ! than the instantaneous value. Formula relating temperature
        ! perturbation to vertical velocity has been used with the most
        ! success at grid lengths near 25 km. For different grid-lengths,
        ! adjust vertical velocity to equivalent value for 25 km grid
        ! length, assuming linear dependence of w on grid length.
        !
        if ( zlcl < 2.0D3 ) then  ! Kain (2004) Eq. 2
          wklcl = 0.02D0 * zlcl/2.0D3
        else
          wklcl = 0.02D0          ! units of m/s
        end if
        wkl = (w0avg(np,k) + &
                (w0avg(np,klcl) - w0avg(np,k))*dlp)*dx/25.0D3 - wklcl
        if ( wkl < 0.0001D0 ) then
          dtlcl = d_zero
        else
          dtlcl = 4.64D0*wkl**0.33D0  ! Kain (2004) Eq. 1
        end if

        if ( kf_trigger == 2 ) then
          dtlcl = max(tpart_h(np,klcl) + tpart_v(np,klcl), d_zero)
        end if

        dtrh = d_zero
        if ( kf_trigger == 3 ) then
          !
          ! for ETA model, give parcel an extra temperature perturbation based
          ! the threshold RH for condensation (U00).
          ! as described in Narita and Ohmori (2007, 12th Mesoscale Conf.)
          ! for now, just assume U00 = 0.75.
          ! !!!!!! for MM5, SET DTRH = 0. !!!!!!!!
          qslcl = qes(k) + (qes(klcl)-qes(k))*dlp
          rhlcl = qenv/qslcl
          dqssdt = qmix*(cliq-bliq*dliq)/((tlcl-dliq)*(tlcl-dliq))
          if ( rhlcl >= 0.75D0 .and. rhlcl <= 0.95D0 ) then
            dtrh = 0.25D0*(rhlcl-0.75D0)*qmix/dqssdt
          else if ( rhlcl > 0.95D0 ) then
            dtrh = (d_one/rhlcl-d_one)*qmix/dqssdt
          else
            dtrh = d_zero
          end if
        end if   ! kf_trigger 3

        kf_trigger2: &
        if ( tlcl+dtlcl+dtrh < tenv ) then
          !
          ! Parcel not buoyant, CYCLE back to start of trigger and
          ! evaluate next potential USL...
          !
          cycle usl
        else  ! Parcel is buoyant, determine updraft
          !
          ! Convective triggering criteria has been satisfied. Compute
          ! equivalent potential temperature
          ! (theteu) and vertical velocity of the rising parcel at the lcl.
          !
          theteu(k) = envirtht(pmix,tmix,qmix)
          !
          ! Modify calculation of initial parcel vertical velocity. jsk 11/26/97
          !
          dttot = dtlcl + dtrh
          if ( dttot > 1.0D-4 ) then
            gdt = d_two*egrav*dttot*500.0D0/tven ! Kain (2004) Eq. 3  (sort of)
            wlcl = d_one + d_half*sqrt(gdt)
            wlcl = min(wlcl,3.0D0)
          else
            wlcl = d_one
          end if
          plcl = p0(np,k) + (p0(np,klcl)-p0(np,k))*dlp
          wtw = wlcl*wlcl
          tvlcl = tlcl*(d_one + ep1*qmix)
          rholcl = plcl/(rdry*tvlcl)
          lcl = klcl
          let = lcl
          ! make RAD a function of background vertical velocity.
          ! (Kain (2004) Eq. 6)
          if ( wkl < d_zero ) then
            rad = d_1000
          else if ( wkl > 0.1D0 ) then
            rad = 2000.0D0
          else
            rad = d_1000 + d_1000*wkl/0.1D0
          end if
          !
          ! Compute updraft properties
          !
          ! Estimate initial updraft mass flux (umf(k)).
          !
          wu(k) = wlcl
          au0 = 0.01D0*dxsq
          umf(k) = rholcl*au0
          vmflcl = umf(k)
          upold = vmflcl
          upnew = upold
          !
          ! Ratio2 is the degree of glaciation in the cloud (0 to 1),
          ! uer is the envir entrainment rate, abe is available
          ! buoyant energy, trppt is the total rate of precipitation
          ! production.
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
          pptliq(np,k) = d_zero
          pptice(np,k) = d_zero
          iflag = 0
          !
          ! ttemp is used during calculation of the linear glaciation
          ! process; it is initially set to the temperature at which
          ! freezing is specified to begin.  Within the glaciation
          ! interval, it is set equal to the updraft temp at the
          ! previous model level.
          !
          ttemp = ttfrz
          !
          ! Enter the loop for updraft calculations. Calculate updraft temp,
          ! mixing ratio, vertical mass flux, lateral detrainment of mass and
          ! moisture, precipitation rates at each model level.
          !
          ! **1 variables indicate the bottom of a model layer
          ! **2 variables indicate the top of a model layer
          !
          ee1 = d_one
          ud1 = d_zero
          rei = d_zero
          dilbe = d_zero
          updraft: &
          do nk = k , kl-1
            nk1 = nk + 1
            ratio2(nk1) = ratio2(nk)
            frc1 = d_zero
            tu(nk1) = t0(np,nk1)
            theteu(nk1) = theteu(nk)
            qu(nk1) = qu(nk)
            qliq(nk1) = qliq(nk)
            qice(nk1) = qice(nk)
            call tpmix2(p0(np,nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1), &
                        qice(nk1),qnewlq,qnewic)
            !
            ! Check to see if updraft temp is above the temperature at which
            ! glaciation is assumed to initiate; if it is, calculate the
            ! fraction of remaining liquid water to freeze. ttfrz is the
            ! temp at which freezing begins, tbfrz the temp below which all
            ! liquid water is frozen at each level.
            !
            if ( tu(nk1) <= ttfrz ) then
              if ( tu(nk1) > tbfrz ) then
                if ( ttemp > ttfrz ) ttemp = ttfrz
                frc1 = (ttemp-tu(nk1))/(ttemp-tbfrz)
              else
                frc1 = d_one
                iflag = 1
              end if
              ttemp = tu(nk1)
              !
              ! Determine the effects of liquid water freezing when temperature
              ! is below ttfrz.
              !
              qfrz = (qliq(nk1)+qnewlq)*frc1
              qnewic = qnewic + qnewlq*frc1
              qnewlq = qnewlq - qnewlq*frc1
              qice(nk1) = qice(nk1) + qliq(nk1)*frc1
              qliq(nk1) = qliq(nk1) - qliq(nk1)*frc1
              call dtfrznew(tu(nk1),p0(np,nk1),theteu(nk1), &
                      qu(nk1),qfrz,qice(nk1))
            end if
            tvu(nk1) = tu(nk1)*(d_one + ep1*qu(nk1))
            !
            ! Calculate updraft vertical velocity and precipitation fallout.
            !
            if ( nk == k ) then
              be = (tvlcl+tvu(nk1))/(tven+tv0(nk1)) - d_one
              boterm = d_two*(z0(np,nk1)-zlcl)*egrav*be/1.5D0
              dzz = z0(np,nk1)-zlcl
            else
              be = (tvu(nk)+tvu(nk1))/(tv0(nk)+tv0(nk1)) - d_one
              boterm = d_two*dza(nk)*egrav*be/1.5D0
              dzz = dza(nk)
            end if
            enterm = d_two*rei*wtw/upold

            call condload(qliq(nk1),qice(nk1),wtw,dzz,boterm,enterm, &
                          qnewlq,qnewic,qlqout(nk1),qicout(nk1))
            !
            ! If vert velocity is less than zero, exit the updraft loop and,
            ! if cloud is tall enough, finalize updraft calculations.
            !
            if ( wtw < 1.0D-3 ) then
              exit updraft
            else
              wu(nk1) = sqrt(wtw)
            end if
            !
            ! Calculate value of theta-e in environment to entrain into updraft
            !
            thetee(nk1) = envirtht(p0(np,nk1),t0(np,nk1),q0(nk1))
            !
            ! rei is the rate of environmental inflow.
            !
            rei = vmflcl*dp(nk1)*0.03D0/rad ! KF (1990) Eq. 1; Kain (2004) Eq. 5
            tvqu(nk1) = tu(nk1)*(d_one + ep1*qu(nk1)-qliq(nk1)-qice(nk1))
            if ( nk == k ) then
              dilbe = ((tvlcl+tvqu(nk1))/(tven+tv0(nk1)) - d_one)*dzz
            else
              dilbe = ((tvqu(nk)+tvqu(nk1))/(tv0(nk)+tv0(nk1)) - d_one)*dzz
            end if
            if ( dilbe > d_zero ) abe = abe + dilbe*egrav
            !
            ! If cloud parcels are virtually colder than the environment,
            ! minimal entrainment (0.5*rei) is imposed.
            !
            if ( tvqu(nk1) <= tv0(nk1) ) then    ! Entrain/Detrain if block
              ee2 = d_half   ! Kain (2004)  Eq. 4
              ud2 = d_one
              eqfrc(nk1) = d_zero
            else
              let = nk1
              ttmp = tvqu(nk1)
              !
              ! Determine the critical mixed fraction of updraft and
              ! environmental air.
              !
              f1 = 0.95D0
              f2 = d_one - f1
              thttmp = f1*thetee(nk1) + f2*theteu(nk1)
              qtmp = f1*q0(nk1) + f2*qu(nk1)
              tmpliq = f2*qliq(nk1)
              tmpice = f2*qice(nk1)
              call tpmix2(p0(np,nk1),thttmp,ttmp,qtmp, &
                      tmpliq,tmpice,qnewlq,qnewic)
              tu95 = ttmp*(d_one + ep1*qtmp-tmpliq-tmpice)
              if ( tu95 > tv0(nk1) ) then
                ee2 = d_one
                ud2 = d_zero
                eqfrc(nk1) = d_one
              else
                f1 = 0.10D0
                f2 = d_one - f1
                thttmp = f1*thetee(nk1) + f2*theteu(nk1)
                qtmp = f1*q0(nk1) + f2*qu(nk1)
                tmpliq = f2*qliq(nk1)
                tmpice = f2*qice(nk1)
                call tpmix2(p0(np,nk1),thttmp,ttmp,qtmp, &
                        tmpliq,tmpice,qnewlq,qnewic)
                tu10 = ttmp*(d_one + ep1*qtmp-tmpliq-tmpice)
                tvdiff = abs(tu10-tvqu(nk1))
                if ( tvdiff < 1.0D-3 ) then
                  ee2 = d_one
                  ud2 = d_zero
                  eqfrc(nk1) = d_one
                else
                  eqfrc(nk1) = (tv0(nk1)-tvqu(nk1))*f1/(tu10-tvqu(nk1))
                  eqfrc(nk1) = max(d_zero,eqfrc(nk1))
                  eqfrc(nk1) = min(d_one,eqfrc(nk1))
                  if ( eqfrc(nk1)-d_one < dlowval ) then
                    ee2 = d_one
                    ud2 = d_zero
                  else if ( eqfrc(nk1) < dlowval ) then
                    ee2 = d_zero
                    ud2 = d_one
                  else
                    !
                    ! Subroutine prof5 integrates over the gaussian dist to
                    ! determine the fractional entrainment and detrainment rates
                    !
                    call prof5(eqfrc(nk1),ee2,ud2)
                  end if
                end if
              end if
            end if  ! End of Entrain/Detrain if BLOCK
            !
            ! Net entrainment and detrainment rates are given by the average
            ! fractional values in the layer.
            !
            ee2 = max(ee2,d_half)
            ud2 = 1.5D0*ud2
            uer(nk1) = d_half*rei*(ee1+ee2)
            udr(nk1) = d_half*rei*(ud1+ud2)
            !
            ! If the calculated updraft detrainment rate is greater than the
            ! total updraft mass flux, all cloud mass detrains, exit updraft
            ! calculations.
            !
            if ( umf(nk)-udr(nk1) < d_10 ) then
              !
              ! If the calculated detrained mass flux is greater than the
              ! total upd mass flux, impose total detrainment of updraft
              ! mass at the previous model lvl.
              ! First, correct abe calculation if needed.
              !
              if ( dilbe > d_zero ) then
                abe = abe - dilbe*egrav
              end if
              let=nk
              exit updraft
            else
              ee1 = ee2
              ud1 = ud2
              upold = umf(nk) - udr(nk1)
              upnew = upold + uer(nk1)
              umf(nk1) = upnew
              dilfrc(nk1) = upnew/upold
              !
              ! detlq and detic are the rates of detrainment of liquid and
              ! ice in the detraining updraft mass.
              !
              detlq(nk1) = qliq(nk1)*udr(nk1)
              detic(nk1) = qice(nk1)*udr(nk1)
              qdt(nk1) = qu(nk1)
              qu(nk1) = (upold*qu(nk1) + uer(nk1)*q0(nk1)) / upnew
              theteu(nk1) = (theteu(nk1)*upold + thetee(nk1)*uer(nk1)) / upnew
              qliq(nk1) = qliq(nk1)*upold/upnew
              qice(nk1) = qice(nk1)*upold/upnew
              !
              ! pptliq is the rate of generation (fallout) of
              ! liquid precip at a given model lvl, pptice the same for ice,
              ! trppt is the total rate of production of precip up to the
              ! current model level.
              !
              pptliq(np,nk1) = qlqout(nk1)*umf(nk)
              pptice(np,nk1) = qicout(nk1)*umf(nk)
              trppt = trppt + pptliq(np,nk1) + pptice(np,nk1)
              if ( nk1 <= kpbl ) uer(nk1) = uer(nk1) + vmflcl*dp(nk1)/dpthmx
            end if
          end do updraft
          !
          ! Check cloud depth. If cloud is tall enough, estimate the equilibrium
          ! temperature level (let) and adjust mass flux profile at cloud top so
          ! that mass flux decreases to zero as a linear function of pressure
          ! between the let and cloud top.
          !
          ! ltop is the model level just below the level at which vertical
          ! velocity first becomes negative.
          !
          ltop = nk
          cldhgt(lc) = z0(np,ltop) - zlcL
          !
          ! Instead of using the same minimum cloud height (for deep convection)
          ! everywhere, try specifying minimum cloud depth as a function of tlcl
          !
          ! Kain (2004)  Eq. 7
          !
          if ( tlcl > 293.0D0 ) then
            chmin = 4.0D3
          else if ( tlcl <= 293.0D0 .and. tlcl >= 273.0D0 ) then
            chmin = 2.0D3 + d_100 * (tlcl - 273.0D0)
          else if ( tlcl < 273.0D0 ) then
            chmin = 2.0D3
          end if
          do nk = k,ltop
            qc_kf(np,nk) = qliq(nk)
            qi_kf(np,nk) = qice(nk)
          end do
          !
          ! If cloud top height is less than the specified minimum for deep
          ! convection, save value to consider this level as source for
          ! shallow convection, go back up to check next level.
          !
          ! Try specifying minimum cloud depth as a function of TLCL
          !
          ! Do not allow any cloud from this layer if:
          !
          !  1.) if there is no CAPE, or
          !  2.) cloud top is at model level just above LCL, or
          !  3.) cloud top is within updraft source layer, or
          !  4.) cloud-top detrainment layer begins within
          !      updraft source layer.
          !
          if ( ltop <= klcl .or. &
               ltop <= kpbl .or. &
               let+1 <= kpbl ) then  ! No Convection Allowed
            cldhgt(lc) = d_zero
            do nk = k , ltop
              umf(nk) = d_zero
              udr(nk) = d_zero
              uer(nk) = d_zero
              detlq(nk) = d_zero
              detic(nk) = d_zero
              pptliq(np,nk) = d_zero
              pptice(np,nk) = d_zero
              cldfra_dp_kf(np,nk) = d_zero
              cldfra_sh_kf(np,nk) = d_zero
              qc_kf(np,nk) = d_zero
              qi_kf(np,nk) = d_zero
            end do
          else if ( cldhgt(lc) > chmin .and. abe > d_one ) then
            ! Deep Convection allowed
            ishall = 0
            do nk = k , ltop
              cldfra_sh_kf(np,nk) = d_zero
            end do
            exit usl
          else
            !
            ! To disallow shallow convection, comment out next line !!!!!!!!
            !
            ishall = 1
            do nk = k , ltop
              cldfra_dp_kf(np,nk) = d_zero
            end do
            if ( nu == nuchm ) then
              exit usl ! Shallow Convection from this layer
            else
              ! Remember this layer (by virtue of non-zero cldhgt) as potential
              ! shallow-cloud layer
              do nk = k , ltop
                umf(nk) = d_zero
                udr(nk) = d_zero
                uer(nk) = d_zero
                detlq(nk) = d_zero
                detic(nk) = d_zero
                pptliq(np,nk) = d_zero
                pptice(np,nk) = d_zero
                cldfra_dp_kf(np,nk) = d_zero
                cldfra_sh_kf(np,nk) = d_zero
                qc_kf(np,nk) = d_zero
                qi_kf(np,nk) = d_zero
              end do
            end if
          end if
        end if kf_trigger2
      end do usl

      if ( ishall == 1 ) then
        kstart = max(kpbl,klcl)
        let = kstart
      end if
      !
      ! If the let and ltop are the same, detrain all of the updraft mass fl
      ! this level.
      !
      if ( let == ltop ) then
        udr(ltop) = umf(ltop)+udr(ltop)-uer(ltop)
        detlq(ltop) = qliq(ltop)*udr(ltop)*upnew/upold
        detic(ltop) = qice(ltop)*udr(ltop)*upnew/upold
        uer(ltop) = d_zero
        umf(ltop) = d_zero
      else
        !
        ! Begin total detrainment at the level above the let.
        !
        dptt = d_zero
        do nj = let+1 , ltop
          dptt = dptt + dp(nj)
        end do
        dumfdp = umf(let)/dptt
        !
        ! Adjust mass flux profiles, detrainment rates, and precipitation fall
        ! rates to reflect the linear decrease in mass flx between the let and
        !
        do nk = let+1 , ltop
          !
          ! Entrainment is allowed at every level except for LTOP, so disallow
          ! entrainment at LTOP and adjust entrainment rates between
          ! LET and LTOP so the the dilution factor due to entrainment is
          ! not changed but the actual entrainment rate will change due to
          ! forced total detrainment in this layer.
          !
          if ( nk == ltop ) then
            udr(nk) = umf(nk-1)
            uer(nk) = d_zero
            detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
            detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
          else
            umf(nk) = umf(nk-1)-dp(nk)*dumfdp
            uer(nk) = umf(nk)*(d_one - d_one/dilfrc(nk))
            udr(nk) = umf(nk-1) - umf(nk) + uer(nk)
            detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
            detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
          end if
          if ( nk >= let+2 ) then
            trppt = trppt - pptliq(np,nk) - pptice(np,nk)
            pptliq(np,nk) = umf(nk-1)*qlqout(nk)
            pptice(np,nk) = umf(nk-1)*qicout(nk)
            trppt = trppt + pptliq(np,nk) + pptice(np,nk)
          end if
        end do
      end if
      !
      ! Initialize some arrays below cloud base and above cloud top.
      !
      do nk = 1 , ltop
        if ( t0(np,nk) > t00 ) ml = nk
      end do
      do nk = 1 , k
        if ( nk >= lc ) then
          if ( nk == lc ) then
            umf(nk) = vmflcl*dp(nk)/dpthmx
            uer(nk) = vmflcl*dp(nk)/dpthmx
          else if ( nk <= kpbl ) then
            uer(nk) = vmflcl*dp(nk)/dpthmx
            umf(nk) = umf(nk-1)+uer(nk)
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
          cldfra_dp_kf(np,nk) = d_zero
          cldfra_sh_kf(np,nK) = d_zero
          qc_kf(np,nk) = d_zero
          qi_kf(np,nk) = d_zero
        end if
        udr(nk) = d_zero
        qdt(nk) = d_zero
        qliq(nk) = d_zero
        qice(nk) = d_zero
        qlqout(nk) = d_zero
        qicout(nk) = d_zero
        pptliq(np,nk) = d_zero
        pptice(np,nk) = d_zero
        detlq(nk) = d_zero
        detic(nk) = d_zero
        ratio2(nk) = d_zero
        thetee(nk) = envirtht(p0(np,nk),t0(np,nk),q0(nk))
        eqfrc(nk) = d_one
      end do
      ltop1 = min(kx,ltop+1)
      ltopm1 = max(1,ltop-1)
      !
      ! Define variables above cloud top
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
        pptliq(np,nk) = d_zero
        pptice(np,nk) = d_zero
        if ( nk > ltop1 ) then
          tu(nk) = d_zero
          qu(nk) = d_zero
          wu(nk) = d_zero
          cldfra_dp_kf(np,nk) = d_zero
          cldfra_sh_kf(np,nk) = d_zero
          qc_kf(np,nk) = d_zero
          qi_kf(np,nk) = d_zero
        end if
        thta0(nk) = d_zero
        thtau(nk) = d_zero
        ems(nk) = d_zero
        emsd(nk) = d_zero
        tg(nk) = t0(np,nk)
        qg(nk) = q0(nk)
        qlg(nk) = d_zero
        qig(nk) = d_zero
        qrg(nk) = d_zero
        qsg(nk) = d_zero
        omg(nk) = d_zero
      end do
      omg(kx+1) = d_zero
      do nk = 1 , ltop
        ems(nk) = dp(nk)*dxsq*regrav
        emsd(nk) = d_one/ems(nk)
        !
        ! Initialize some variables to be used later
        ! in the vert advection scheme
        !
        exn(nk) = (p00/p0(np,nk))**(0.2854D0*(d_one-0.28D0*qdt(nk)))
        thtau(nk) = tu(nk)*exn(nk)
        exn(nk) = (p00/p0(np,nk))**(0.2854D0*(d_one-0.28D0*q0(nk)))
        thta0(nk) = t0(np,nk)*exn(nk)
        ddilfrc(nk) = d_one/dilfrc(nk)
        omg(nk) = d_zero
      end do
      !
      ! Compute convective time scale(timec). the mean wind at the lcl
      ! and midtroposphere is used.
      !
      wspd(klcl) = sqrt(u0(np,klcl)*u0(np,klcl) + v0(np,klcl)*v0(np,klcl))
      wspd(l5) = sqrt(u0(np,l5)*u0(np,l5) + v0(np,l5)*v0(np,l5))
      wspd(ltop) = sqrt(u0(np,ltop)*u0(np,ltop) + v0(np,ltop)*v0(np,ltop))
      vconv = d_half*(wspd(klcl)+wspd(l5))
      timec = dx/vconv
      tadvec = timec
      ! kf_min_dtcape >= TIMEC <= kf_max_dtcape
      timec = max(kf_min_dtcape,timec)
      timec = min(kf_max_dtcape,timec)
      ! shallow convection TIMEC
      if ( ishall == 1 ) then
        timec = max(d_half*(kf_max_dtcape+kf_min_dtcape)-300.0D0,300.0D0)
      end if
      nic = nint(timec/dt)
      timec = dble(nic)*dt
      !
      ! Compute wind shear and precipitation efficiency.
      !
      if ( wspd(ltop) > wspd(klcl) ) then
        shsign = d_one
      else
        shsign = -d_one
      end if
      vws = (u0(np,ltop)-u0(np,klcl))*(u0(np,ltop)-u0(np,klcl)) + &
            (v0(np,ltop)-v0(np,klcl))*(v0(np,ltop)-v0(np,klcl))
      vws = 1.0D3*shsign*sqrt(vws)/(z0(np,ltop)-z0(np,lcl))
      pef = 1.591D0 + vws*(-0.639D0 + vws*(9.53D-2 - vws*4.96D-3))
      pef = max(pef,kf_min_pef)
      pef = min(pef,kf_max_pef)
      !
      ! Precipitation efficiency is a function of the height of cloud base.
      !
      cbh = (zlcl-z0(np,1))*3.281D-3
      if ( cbh < 3.0D0 ) then
        rcbh = 0.02D0
      else
        rcbh = 0.96729352D0 + cbh*(-0.70034167D0 + cbh*(0.162179896D0 + &
                              cbh*(-1.2569798D-2 + cbh*(4.2772D-4 - &
                              cbh*5.44D-6))))
      end if
      if ( cbh > 25.0D0 ) rcbh = 2.4D0
      pefcbh = d_one/(d_one+rcbh)
      pefcbh = max(pefcbh,kf_min_pef)
      pefcbh = min(pefcbh,kf_max_pef)
      !
      ! mean pef. is used to compute rainfall.
      !
      peff = d_half*(pef+pefcbh)
      peff2 = peff   ! jsk mods
      if ( iprnt ) then
        write(stdout,1035) pef,pefcbh,lc,let,wkl,vwS
      end if
      !
      ! Compute downdraft properties
      !
      tder = d_zero
      devap: &
      if ( ishall == 1 ) then
        lfs = 1
      else
        !
        ! Start downdraft about above cloud base.
        !
        kstart = kpbl+1
        klfs = let-1
        findklfs: &
        do nk = kstart+1 , kl
          dppp = p0(np,kstart) - p0(np,nk)
          if ( dppp > kf_dpp * d_100 ) then
            klfs = nk
            exit findklfs
          end if
        end do findklfs
        klfs = min(klfs,let-1)
        lfs = klfs
        !
        ! If lfs is not at least 50 mb above cloud base (implying that the
        ! level of equil temp, let, is just above cloud base) do not allow a
        ! downdraft.
        !
        if ( (p0(np,kstart)-p0(np,lfs)) > 50.0D2 ) then
          theted(lfs) = thetee(lfs)
          qd(lfs) = q0(lfs)
          !
          ! Call tpmix2dd to find wet-bulb temp, qv
          !
          call tpmix2dd(p0(np,lfs),theted(lfs),tz(lfs),qss)
          thtad(lfs) = tz(lfs)*(p00/p0(np,lfs))**(0.2854D0*(d_one-0.28D0*qss))
          !
          ! Take a first guess at the initial downdraft mass flux
          !
          tvd(lfs) = tz(lfs)*(d_one + ep1*qss)
          rdd = p0(np,lfs)/(rdry*tvd(lfs))
          a1 = (d_one-peff)*au0
          dmf(lfs) = -a1*rdd
          der(lfs) = dmf(lfs)
          ddr(lfs) = d_zero
          rhbar = rh(lfs)*dp(lfs)
          dptt = dp(lfs)
          do nd = lfs-1 , kstart , -1
            nd1 = nd+1
            der(nd) = der(lfs)*ems(nd)/ems(lfs)
            ddr(nd) = d_zero
            dmf(nd) = dmf(nd1) + der(nd)
            theted(nd) = (theted(nd1)*dmf(nd1) + thetee(nd)*der(nd))/dmf(nd)
            qd(nd) = (qd(nd1)*dmf(nd1) + q0(nd)*der(nd))/dmf(nd)
            dptt = dptt + dp(nd)
            rhbar = rhbar + rh(nd)*dp(nd)
          end do
          rhbar = rhbar/dptt
          dmffrc = d_two*(d_one-rhbar) ! Kain (2004) eq. 11
          dpdd = d_zero
          ! Calculate melting effect
          ! first, compute total frozen precipitation generated.
          !
          pptmlt = 0.
          do nk = klcl , ltop
            pptmlt = pptmlt + pptice(np,nk)
          end do
          if ( lc < ml ) then
            ! For now, calculate melting effect as if dmf = -umf at klcl,
            ! i.e., as if dmffrc=1.
            ! Otherwise, for small dmffrc, dtmelt gets too large!
            ! 12/14/98 jsk...
            dtmelt = wlhf*pptmlt/(cpd*umf(klcl))
          else
            dtmelt = d_zero
          end if
          ldt = min(lfs-1,kstart-1)
          call tpmix2dd(p0(np,kstart),theted(kstart),tz(kstart),qss)
          tz(kstart) = tz(kstart) - dtmelt
          es = aliq*exp((bliq*tz(kstart)-cliq)/(tz(kstart)-dliq))
          qss = ep2*es/(p0(np,kstart)-es)
          theted(kstart) = tz(kstart) * &
               (p00/p0(np,kstart))**(0.2854D0*(d_one-0.28D0*qss))*    &
                exp((c1/tz(kstart)-c2)*qss*(d_one+c4*qss))
          ldt = min(lfs-1,kstart-1)
          findldb: &
          do nd = ldt , 1 , -1
            dpdd = dpdd + dp(nd)
            theted(nd) = theted(kstart)
            qd(nd) = qd(kstart)
            !
            ! Call tpmix2dd to find wet bulb temp, saturation mixing ratio
            !
            call tpmix2dd(p0(np,nd),theted(nd),tz(nd),qss)
            qsd(nd) = qss
            !
            ! Specify RH decrease of 20%/km in downdraft...
            !
            rhh = d_one-0.2D0/d_1000*(z0(np,kstart)-z0(np,nd))
            !
            ! Adjust downdraft TEMP, Q to specified RH:
            !
            if ( rhh < d_one ) then
              dssdt = (cliq-bliq*dliq)/((tz(nd)-dliq)*(tz(nd)-dliq))
              rl = xlv0 - xlv1 * tz(nd)
              dtmp = rl*qss*(d_one-rhh)/(cpd+rl*rhh*qss*dssdt)
              t1rh = tz(nd) + dtmp
              es = rhh*aliq*exp((bliq*t1rh-cliq)/(t1rh-dliq))
              qsrh = ep2*es/(p0(np,nd)-es)
              !
              ! Check to see if mixing ratio at specified rh is less than actual
              ! mixing ratio. If so, adjust to give zero evaporation.
              !
              if ( qsrh < qd(nd) ) then
                qsrh = qd(nd)
                t1rh = tz(nd) + (qss-qsrh)*rl*rcpd
              end if
              tz(nd) = t1rh
              qss = qsrh
              qsd(nd) = qss
            end if
            tvd(nd) = tz(nd)*(d_one + ep1*qsd(nd))
            if ( tvd(nd) > tv0(nd) .or. nd == 1 ) then
              ldb = nd
              exit findldb
            end if
          end do findldb
          if ( (p0(np,ldb)-p0(np,lfs)) > 50.0D2 ) then
            ! minimum Downdraft depth!
            do nd = ldt , ldb , -1
              nd1 = nd+1
              ddr(nd) = -dmf(kstart)*dp(nd)/dpdd
              der(nd) = d_zero
              dmf(nd) = dmf(nd1)+ddr(nd)
              tder = tder + (qsd(ND)-qd(nd))*ddr(nd)
              qd(nd) = qsd(nd)
              thtad(nd) = tz(nd) * &
                      (p00/p0(np,nd))**(0.2854D0*(d_one-0.28D0*qd(nd)))
            end do
          end if
        end if
      end if devap
      !
      ! If downdraft does not evaporate any water for specified relative
      ! humidity, no downdraft is allowed.
      !
      d_mf: &
      if ( tder < d_one ) then
        pptflx = trppt
        cpr = trppt
        tder = 0.
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
        aincm2 = d_100
      else
        ddinc = -dmffrc*umf(klcl)/dmf(kstart)
        updinc = d_one
        if ( tder*ddinc > trppt ) then
          ddinc = trppt/tder
        end if
        tder = tder*ddinc
        do nk = ldb , lfs
          dmf(nk) = dmf(nk)*ddinc
          der(nk) = der(nk)*ddinc
          ddr(nk) = ddr(nk)*ddinc
        end do
        cpr = trppt
        pptflx = trppt-tder
        peff = pptflx/trppt
        !
        ! Adjust updraft mass flux, mass detrainment rate, and liquid water an
        ! detrainment rates to be consistent with the transfer of the estimate
        ! from the updraft to the downdraft at the lfs.
        !
        ! Zero out the arrays for downdraft data at levels above and below the
        ! downdraft.
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
        do nk = lfs+1 , kx
          dmf(nk) = d_zero
          der(nk) = d_zero
          ddr(nk) = d_zero
          wd(nk) = d_zero
          tz(nk) = d_zero
          qd(nk) = d_zero
          thtad(nk) = d_zero
        end do
        do nk = ldt+1 , lfs-1
          tz(nk) = d_zero
          qd(nk) = d_zero
          thtad(nk) = d_zero
        end do
      end if d_mf
      !
      ! Set limits on the updraft and downdraft mass fluxes so that the inflow
      ! into convective drafts from a given layer is no more than is available
      ! in that layer initially
      !
      aincmx = d_1000
      lmax = max(klcl,lfs)
      do nk = lc , lmax
        if ( (uer(nk)-der(nk)) > 1.0D-3 ) then
          aincm1 = ems(nk)/((uer(nk)-der(nk))*timec)
          aincmx = min(aincmx,aincm1)
        end if
      end do
      ainc = d_one
      if ( aincmx < ainc ) ainc = aincmx
      !
      ! Save the relevent variables for a unit updraft and downdraft. They will
      ! be iteratively adjusted by the factor ainc to satisfy the stabilization
      ! closure.
      !
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
      stab = 0.95D0
      noitr = 0
      istop = 0
      if ( ishall == 1 ) then  ! First for shallow convection
        !
        ! No iteration for shallow convection; if turbulent kinetic energy
        ! (TKE) is available from a turbulence parameterization, scale
        ! cloud-base updraft mass flux as a function of TKE, but for now,
        ! just specify shallow-cloud mass flux using TKEMAX = 5...
        !
        ! find the maximum TKE value between LC and KLCL...
        evac = d_half*maxval(tke(np,lc:klcl))*0.1D0
        ainc = evac*dpthmx*dxsq/(vmflcl*egrav*timec)
        tder = tder2*ainc
        pptflx = pptfl2*ainc
        do nk = 1 , ltop
          umf(nk) = umf2(nk)*ainc
          dmf(nk) = dmf2(nk)*ainc
          detlq(nk) = detlq2(nk)*ainc
          detic(nk) = detic2(nk)*ainc
          udr(nk) = udr2(nk)*ainc
          uer(nk) = uer2(nk)*ainc
          der(nk) = der2(nk)*ainc
          ddr(nk) = ddr2(nk)*ainc
        end do
      end if  ! Otherwise for deep convection
              ! use iterative procedure to find mass fluxes...
      iter: &
      do ncount = 1 , maxiter
        !
        ! Compute properties for compensational subsidence
        !
        ! Determine omega value necessary at top and bottom of each layer to
        ! satisfy mass continuity.
        !
        dtt = timec
        do nk = 1 , ltop
          domgdp(nk) = -(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
          if ( nk > 1 ) then
            omg(nk) = omg(nk-1) - dp(nk-1)*domgdp(nk-1)
            absomg = abs(omg(nk))
            absomgtc = absomg*timec
            frdp = 0.75D0*dp(nk-1)
            if ( absomgtc > frdp ) then
              dtt1 = frdp/absomg
              dtt = min(dtt,dtt1)
            end if
          end if
        end do
        do nk = 1 , ltop
          thpa(nk) = thta0(nk)
          qpa(nk) = q0(nk)
          nstep = nint(timec/dtt+1)
          dtime = timec/dble(nstep)
          fxm(nk) = omg(nk)*dxsq*regrav
        end do
        !
        ! Do an upstream/forward-in-time advection of theta, qv
        !
        do ntc = 1 , nstep
          !
          ! Assign theta and q values at the top and bottom of each layer
          ! based on sign of omega
          !
          do nk = 1 , ltop
            thfxin(nk) = d_zero
            thfxout(nk) = d_zero
            qfxin(nk) = d_zero
            qfxout(nk) = d_zero
          end do
          do nk = 2 , ltop
            if ( omg(nk) <= d_zero ) then
              thfxin(nk) = -fxm(nk)*thpa(nk-1)
              qfxin(nk) = -fxm(nk)*qpa(nk-1)
              thfxout(nk-1) = thfxout(nk-1) + thfxin(nk)
              qfxout(nk-1) = qfxout(nk-1) + qfxin(nk)
            else
              thfxout(nk) = fxm(nk)*thpa(nk)
              qfxout(nk) = fxm(nk)*qpa(nk)
              thfxin(nk-1) = thfxin(nk-1) + thfxout(nk)
              qfxin(nk-1) = qfxin(nk-1) + qfxout(nk)
            end if
          end do
          !
          ! Update the theta and qv values at each level
          !
          do nk = 1 , ltop
            thpa(nk) = thpa(nk)+(thfxin(nk)+udr(nk)*thtau(nk)+ddr(nk)*      &
                       thtad(nk)-thfxout(nk)-(uer(nk)-der(nk))*thta0(nk))*  &
                       dtime*emsd(nk)
            qpa(nk) = qpa(nk)+(qfxin(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk)-    &
                      qfxout(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)
          end do
        end do
        do nk = 1 , ltop
          thtag(nk) = thpa(nk)
          qg(nk) = qpa(nk)
        end do
        !
        ! Check to see if mixing ratio dips below zero anywhere;  if so, borrow
        ! moisture from adjacent layers to bring it back up above zero
        !
        do nk = 1 , ltop
          if ( qg(nk) < d_zero ) then
            write(stderr,*) 'AT I = ',imap(np), ', J = ', jmap(np)
            write(stderr,*) 'KF: QG, QG(NK) < 0'
            nk1 = nk + 1
            if ( nk == ltop ) then
              nk1 = klcl
            end if
            tma = qg(nk1)*ems(nk1)
            tmb = qg(nk-1)*ems(nk-1)
            tmm = (qg(nk)-1.0D-9)*ems(nk)
            bcoeff = -tmm/((tma*tma)/tmb+tmb)
            acoeff = bcoeff*tma/tmb
            tmb = tmb*(d_one-bcoeff)
            tma = tma*(d_one-acoeff)
            qg(nk) = 1.0D-9
            qg(nk1) = tma*emsd(nk1)
            qg(nk-1) = tmb*emsd(nk-1)
          end if
        end do
        !
        ! Convert theta to t
        !
        do nk = 1 , ltop
          exn(nk) = (p00/p0(np,nk))**(0.2854D0*(d_one-0.28D0*qg(nk)))
          tg(nk) = thtag(nk)/exn(nk)
          tvg(nk) = tg(nk)*(d_one + ep1*qg(nk))
        end do
        if ( ishall == 1 ) then
          exit iter
        end if
        !
        ! Compute new cloud and change in available buoyant energy.
        !
        ! The following computations are similar to that for updraft
        !
        tmix = d_zero
        qmix = d_zero
        !
        ! Find the thermodynamic characteristics of the layer by
        ! mass-weighting the characteristics of the individual model
        ! layers.
        !
        do nk = lc , kpbl
          tmix = tmix + dp(nk)*tg(nk)
          qmix = qmix + dp(nk)*qg(nk)
        end do
        tmix = tmix/dpthmx
        qmix = qmix/dpthmx
        es = aliq*exp((tmix*bliq-cliq)/(tmix-dliq))
        qss = ep2*es/(pmix-es)
        !
        ! Remove supersaturation for diagnostic purposes, if necessary
        !
        if ( qmix > qss ) then
          rl = xlv0-xlv1*tmix
          cpm = cpd*(d_one+0.887D0*qmix)
          dssdt = qss*(cliq-bliq*dliq)/((tmix-dliq)*(tmix-dliq))
          dq = (qmix-qss)/(d_one+rl*dssdt/cpm)
          tmix = tmix + rl*rcpd*dq
          qmix = qmix - dq
          tlcl = tmix
        else
          qmix = max(qmix,d_zero)
          emix = qmix*pmix/(ep2+qmix)
          binc = aincb
          a1 = emix/aliq
          tp = (a1-astrt)/binc
          indlu = max(1, min(kfna-1,int(tp)+1))
          avalue = (indlu-1)*binc+astrt
          aintrp = (a1-avalue)/binc
          tlog = aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
          tdpt = (cliq-dliq*tlog)/(bliq-tlog)
          tlcl = tdpt - (0.212D0+1.571D-3*(tdpt-t00) - &
                         4.36D-4*(tmix-t00))*(tmix-tdpt)
          tlcl = min(tlcl,tmix)
        end if
        tvlcl = tlcl*(d_one + ep1*qmix)
        zlcl = zmix + (tlcl-tmix)/gdry
        findklcl2: &
        do nk = lc , kl
          klcl = nk
          if ( zlcl <= z0(np,nk) ) then
            exit findklcl2
          end if
        end do findklcl2
        klcl = max(2,klcl)
        k = klcl - 1
        dlp = (zlcl-z0(np,k))/(z0(np,klcl)-z0(np,k))
        !
        ! Estimate environmental temperature and mixing ratio at the lcl
        !
        tenv = tg(k) + (tg(klcl)-tg(k))*dlp
        qenv = qg(k) + (qg(klcl)-qg(k))*dlp
        tven = tenv*(d_one + ep1*qenv)
        plcl = p0(np,k) + (p0(np,klcl)-p0(np,k))*dlp
        theteu(k) = tmix*(p00/pmix)**(0.2854D0*(d_one-0.28D0*qmix))* &
               exp((c1/tlcl-c2)*qmix*(d_one+c4*qmix))
        !
        ! Compute adjusted abe(abeg).
        !
        abeg = d_zero
        do nk = k , ltopm1
          nk1 = nk + 1
          theteu(nk1) = theteu(nk)
          call tpmix2dd(p0(np,nk1),theteu(nk1),tgu(nk1),qgu(nk1))
          tvqu(nk1) = tgu(nk1)*(d_one + ep1*qgu(nk1)-qliq(nk1)-qice(nk1))
          if ( nk == k ) then
            dzz = z0(np,klcl) - zlcl
            dilbe = ((tvlcl+tvqu(nk1))/(tven+tvg(nk1))-d_one)*dzz
          else
            dzz = dza(nk)
            dilbe = ((tvqu(nk)+tvqu(nk1))/(tvg(nk)+tvg(nk1))-d_one)*dzz
          end if
          if ( dilbe > d_zero ) abeg = abeg + dilbe*egrav
          !
          ! Dilute by entrainment by the rate as original updraft
          !
          thteeg(nk1) = envirtht(p0(np,nk1),tg(nk1),qg(nk1))
          theteu(nk1) = theteu(nk1)*ddilfrc(nk1) + &
                        thteeg(nk1)*(d_one-ddilfrc(nk1))
        end do
        !
        ! Assume at least 90% of cape (abe) is removed by convection during
        ! the period timec
        !
        if ( noitr == 1 ) then
          exit iter
        end if
        dabe = max(abe-abeg, 0.1D0*abe)
        fabe = abeg/abe
        if ( fabe > d_one .and. ishall == 0 ) then
          cycle modelpoints
        end if
        if ( ncount /= 1 ) then
          if ( abs(ainc-aincold) < 0.00001D0 ) then
            noitr = 1
            ainc = aincold
            cycle iter
          end if
          dfda = (fabe-fabeold)/(ainc-aincold)
          if ( dfda > d_zero ) then
            noitr = 1
            ainc = aincold
            cycle iter
          end if
        end if
        aincold = ainc
        fabeold = fabe
        topomg = (udr(ltop)-uer(ltop))*dp(ltop)*emsd(ltop)
        if ( ainc/aincmx > 0.999D0 .and. fabe > 1.05D0-stab ) then
          exit iter
        end if
        if ( fabe <= 1.05D0-stab .and. fabe >= 0.95D0-stab ) then
          exit iter
        else
          !
          ! If more than 10% of the original cape remains, increase the
          ! convective mass flux by the factor ainc:
          !
          if ( abs(fabe) < dlowval ) then
            ainc = ainc*d_half
          else
            if ( dabe < 1.0D-4 ) then
              noitr = 1
              ainc = aincold
              cycle iter
            else
              ainc = ainc*stab*abe/dabe
            end if
          end if
          ainc = min(aincmx,ainc)
          ! If ainc becomes very small, effects of convection ! jsk mods
          ! will be minimal so just ignore it                 ! jsk mods
          if ( ainc < 0.05D0 ) then
            cycle modelpoints
          end if
          tder = tder2*ainc
          pptflx = pptfl2*ainc
          do nk = 1 , ltop
            umf(nk) = umf2(nk)*ainc
            dmf(nk) = dmf2(nk)*ainc
            detlq(nk) = detlq2(nk)*ainc
            detic(nk) = detic2(nk)*ainc
            udr(nk) = udr2(nk)*ainc
            uer(nk) = uer2(nk)*ainc
            der(nk) = der2(nk)*ainc
            ddr(nk) = ddr2(nk)*ainc
          end do
          !
          ! Go back up for another iteration.
          !
        end if
        if ( abs(topomg-omg(ltop)) > 1.0D-3 ) then
          iprnt = .true.
          write (stderr, *) 'POSSIBLE INSTABILITY IN KF CODE'
          write (stderr, *) 'MASS DOES NOT BALANCE IN KF SCHEME'
          write (stderr, *) 'NCOUNT = ', ncount
          if ( ncount == maxiter ) istop = 1
        else
          iprnt = .false.
        end if
      end do iter
      ! Get the cloud fraction for layer NK+1=NK1
      if ( ishall == 1 )  then
        do nk = klcl-1 , ltop1
          umf_new = umf(nk)/dxsq
          xcldfra = 0.07D0*log(d_one+(500.0D0*umf_new))
          xcldfra = max(0.01D0,xcldfra)
          cldfra_sh_kf(np,nk) = min(d_half*clfrcv,xcldfra)
        end do
      else
        do nk = klcl-1 , ltop1
          umf_new = umf(nk)/dxsq
          xcldfra = 0.14D0*log(d_one+(500.0D0*umf_new))
          xcldfra = max(0.01D0,xcldfra)
          cldfra_dp_kf(np,nk) = min(clfrcv,xcldfra)
        end do
      end if
      !
      ! Compute hydrometeor tendencies as is done for t, qv
      !
      ! frc2 is the fraction of total condensate      !  ppt fb mods
      ! generated that goes into precipitiation       !  ppt fb mods
      !
      ! Redistribute hydormeteors according to the final mass-flux values:
      !
      if ( cpr > d_zero ) then
        frc2 = pptflx/(cpr*ainc) !  ppt fb mods
      else
        frc2 = d_zero
      end if
      do nk = 1 , ltop
        qlpa(nk) = ql0(nk)
        qipa(nk) = qi0(nk)
        qrpa(nk) = qr0(nk)
        qspa(nk) = qs0(nk)
        rainfb(nk) = pptliq(np,nk)*ainc*fbfrc*frc2   !  ppt fb mods
        snowfb(nk) = pptice(np,nk)*ainc*fbfrc*frc2   !  ppt fb mods
      end do
      do ntc = 1 , nstep
        !
        ! Assign hydrometeors concentrations at the top and bottom of each layer
        ! based on the sign of omega
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
            qlfxout(nk-1) = qlfxout(nk-1) + qlfxin(nk)
            qifxout(nk-1) = qifxout(nk-1) + qifxin(nk)
            qrfxout(nk-1) = qrfxout(nk-1) + qrfxin(nk)
            qsfxout(nk-1) = qsfxout(nk-1) + qsfxin(nk)
          else
            qlfxout(nk) = fxm(nk)*qlpa(nk)
            qifxout(nk) = fxm(nk)*qipa(nk)
            qrfxout(nk) = fxm(nk)*qrpa(nk)
            qsfxout(nk) = fxm(nk)*qspa(nk)
            qlfxin(nk-1) = qlfxin(nk-1) + qlfxout(nk)
            qifxin(nk-1) = qifxin(nk-1) + qifxout(nk)
            qrfxin(nk-1) = qrfxin(nk-1) + qrfxout(nk)
            qsfxin(nk-1) = qsfxin(nk-1) + qsfxout(nk)
          end if
        end do
        !
        ! Update the hydrometeor concentration values at each level
        !
        do nk = 1 , ltop
          qlpa(nk) = qlpa(nk)+(qlfxin(nk)+detlq(nk)-qlfxout(nk))*dtime*emsd(nk)
          qipa(nk) = qipa(nk)+(qifxin(nk)+detic(nk)-qifxout(nk))*dtime*emsd(nk)
          ! ppt fb mods
          qrpa(nk) = qrpa(nk)+(qrfxin(nk)-qrfxout(nk)+rainfb(nk))*dtime*emsd(nk)
          qspa(nk) = qspa(nk)+(qsfxin(nk)-qsfxout(nk)+snowfb(nk))*dtime*emsd(nk)
          ! ppt fb mods
        end do
      end do
      do nk = 1 , ltop
        qlg(nk) = qlpa(nk)
        qig(nk) = qipa(nk)
        qrg(nk) = qrpa(nk)
        qsg(nk) = qspa(nk)
      end do
      !
      ! Clean things up, calculate convective feedback tendencies for this
      ! grid point
      !
      if ( iprnt ) then
        write(stderr,1080) LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,FABE,AINC
      end if
      !
      ! Send final parameterized values to output files
      !
      if ( iprnt .or. istop == 1 ) then
        write(stdout,*)
        write(stdout,*) 'P(LC), DTP, WKL, WKLCL =',p0(np,lc)/d_100, &
            tlcl+dtlcl+dtrh-tenv,wkl,wklcL
        write(stdout,*) 'TLCL, DTLCL, DTRH, TENV =',tlcl,dtlcl,    &
            dtrh,tenv
        write(stdout,1025) klcl,zlcl,dtlcl,ltop,p0(np,ltop),iflag,    &
            tmix-t00,pmix,qmix,abe
1025  format(5X,' KLCL=',I2,' ZLCL=',F7.1,'M',                            &
          ' DTLCL=',F5.2,' LTOP=',I2,' P0(LTOP)=',-2PF5.1,'MB FRZ LV=',   &
          I2,' TMIX=',0PF4.1,1X,'PMIX=',-2PF6.1,' QMIX=',3PF5.1,          &
          ' CAPE=',0PF7.1)
        write(stdout,1030) p0(np,let)/d_100,p0(np,ltop)/d_100, &
            vmflcl,plcl/d_100,wlcl,cldhgt(lc)
1030  format(' ',' P0(LET) = ',F6.1,' P0(LTOP) = ',F6.1,' VMFLCL =',    &
        E12.3,' PLCL =',F6.1,' WLCL =',F6.3,' CLDHGT =',F8.1)
        write(stdout,1035) pef,pefcbh,lc,let,wkl,vws
        write(stdout,*) 'PRECIP EFFICIENCY =' , peff
        write(stdout,1080) lfs,ldb,ldt,timec,tadvec,nstep,ncount,fabe,ainc
        write(stdout,'(16a8)') '  P  ','   DP ',' DT K/D ',' DR K/D ',   &
                '   OMG  ', ' DOMGDP ','   UMF  ','   UER  ','   UDR  ', &
                '   DMF  ','   DER  ','   DDR  ','   EMS  ','    W0  ',  &
                '  DETLQ ',' DETIC '
        do nk = 1 , ltop
          k = ltop - nk + 1
          dtt = (tg(k)-t0(np,k))*86400.0D0/timec
          rl = xlv0-xlv1*tg(k)
          dr = -(qg(k)-q0(k))*rl*86400.0D0/(timec*cpd)
          udfrc = udr(k)*timec*emsd(k)
          uefrc = uer(k)*timec*emsd(k)
          ddfrc = ddr(k)*timec*emsd(k)
          defrc = -der(k)*timec*emsd(k)
          write(stdout,1075) p0(np,k)/d_100,dp(k)/d_100,dtt,dr,omg(k),     &
                  domgdp(k)*1.0D4, umf(k)/1.0D6,uefrc,udfrc,dmf(k)/1.0D6,  &
                  defrc,ddfrc,ems(k)/1.0D11,w0avg(np,k)*1.0D2,             &
                  detlq(k)*timec*emsd(K)*1.0D3,detic(k)*timec*emsd(k)*1.0D3
   1075   format(F8.2,3(F8.2),2(F8.3),F8.2,2F8.3,F8.2,6F8.3)
        end do
        write(stdout,1085) 'K','P','Z','T0','TG','DT','TU','TD','Q0','QG', &
                    'DQ','QU','QD','QLG','QIG','QRG','QSG','RH0','RHG'
        do nk = 1 , kl
          k = kx - nk + 1
          dtt = tg(k) - t0(np,k)
          tuc = tu(k) - t00
          if ( k < lc .or. k > ltop ) tuc = d_zero
          tdc = tz(k)-t00
          if ( (k < ldb .or. k > ldt) .and. k /= lfs ) tdc = d_zero
          if ( t0(np,k) < t00 ) then
            es = aliq*exp((bliq*tg(k)-cliq)/(tg(k)-dliq))
          else
            es = aliq*exp((bliq*tg(k)-cliq)/(tg(k)-dliq))
          end if
          qgs = es*ep2/(p0(np,k)-es)
          rh0 = q0(k)/qes(k)
          rhg = qg(k)/qgs
          write(stdout,1090) k,p0(np,k)/d_100,z0(np,k),t0(np,k)-t00, &
                  tg(k)-t00,dtt,tuc,tdc,q0(k)*d_1000,qg(k)*d_1000,   &
                  (qg(k)-q0(k))*d_1000,qu(k)*d_1000,QD(K)*d_1000,    &
                  qlg(k)*d_1000,qig(k)*d_1000,qrg(k)*d_1000,         &
                  qsg(k)*d_1000,rh0,rhg
   1090   format(I3,F7.2,F7.0,10F7.2,4F7.3,2F8.3)
        end do
        !
        ! If calculations above show an error in the mass budget, print out a
        ! to be used later for diagnostic purposes, THEN abort run.
        !
        do nk = 1 , kl
          k = kl - nk + 1
          write(stdout,'(8f11.3)') p0(np,k)/d_100,t0(np,k)-t00, &
                  q0(k)*d_1000,u0(np,k),v0(np,k),w0avg(np,k)*d_100, &
                  dp(k),tke(np,k)
        end do
        if ( istop == 1 ) then
          write(stderr,*) 'AT I = ', imap(np), ', J = ', jmap(np)
          call fatal(__FILE__,__LINE__,'KAIN-FRITSCH, istop=1, diags')
        end if
      end if
      cndtnf = (d_one-eqfrc(lfs))*(qliq(lfs)+qice(lfs))*dmf(lfs)
      pratec(np) = pptflx*(d_one-fbfrc)/dxsq
      !
      ! Evaluate moisture budget
      !
      qinit = d_zero
      qfnl = d_zero
      do nk = 1 , ltop
        qinit = qinit + q0(nk)*ems(nk)/dxsq
        qfnl = qfnl + (qg(nk)+qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)/dxsq
      end do
      qfnl = qfnl + pptflx*timec*(d_one-fbfrc)/dxsq  !  ppt fb mods
      err2 = (qfnl-qinit)*d_100/qinit
      if ( abs(err2) > 0.05D0 .and. istop == 0 ) then
        istop = 1
        iprnt = .true.
        write(stdout,1110) qinit , qfnl , err2
        write(stdout,'(a,f12.3)') 'PPTFLX = ',(pptflx*(d_one-fbfrc)*timec)/dxsq
        write(stdout,'(a,i6,a,i2)') 'LEVELS = ', kl,', ISHALLOW = ',ishall
        write(stdout,'(a,i2)') 'NCOUNT : ',ncount
        do nk = 1 , kl
          k = kl - nk + 1
          write(stdout,'(8f12.3)') p0(np,k)/d_100,t0(np,k)-t00, &
                  q0(k)*d_1000,u0(np,k),v0(np,k),w0avg(np,k),   &
                  dp(k)/d_100,tke(np,k)
          write(stdout,'(5f12.3)') qg(nk)*d_1000,qlg(nk)*d_1000, &
                  qig(nk)*d_1000,qrg(nk)*d_1000,qsg(nk)*d_1000
        end do
      end if
      if ( pptflx > d_zero ) then
        relerr = err2*qinit/(pptflx*timec)
      else
        relerr = d_zero
      end if
      if ( iprnt .and. relerr > dlowval ) then
        write(stdout,'(a,f9.3,a)') &
                ' MOISTURE ERROR AS FUNCTION OF TOTAL PPT = ', relerr, '%'
        write(stdout,*) 'TDER, CPR, TRPPT = ', tder, cpr*ainc, trppt*ainc
      end if
      !
      ! Feedback to resolvable scale tendencies.
      !
      ! If the advective time period (tadvec) is less than specified minimum
      ! timec, allow feedback to occur only during tadvec.
      !
      if ( tadvec < timec ) nic = nint(tadvec/dt)
      if ( ishall == 1 ) then
        timec = max(d_half*(kf_max_dtcape+kf_min_dtcape)-300.0D0,300.0D0)
      end if
      do k = 1 , kx
        !
        ! If hydrometeors are not allowed, they must be evaporated or sublimated
        ! and fed back as vapor, along with associated changes in temperature.
        ! Note:  this will introduce changes in the convective temperature and
        ! water vapor feedback tendencies and may lead to supersaturated value
        ! of qg.
        !
        ! If ice phase is not allowed, melt all frozen hydrometeors.
        !
        if ( .not. f_qi .and. warm_rain ) then
          cpm = cpd*(d_one + 0.887D0*qg(k))
          tg(k) = tg(k) - (qig(k)+qsg(k))*wlhf/cpm
          dqcdt(np,k) = (qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
          dqidt(np,k) = d_zero
          dqrdt(np,k) = (qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
          dqsdt(np,k) = d_zero
        else if ( .not. f_qi .and. .not. warm_rain ) then
          !
          ! If ice phase is allowed, but mixed phase is not, melt frozen
          ! hydrometeors below the melting level, freeze liquid water above
          ! the melting level
          !
          cpm = cpd*(d_one + 0.887D0*qg(k))
          if ( k <= ml ) then
            tg(k) = tg(k) - (qig(k)+qsg(k))*wlhf/cpm
          else if ( k > ml ) then
            tg(k) = tg(k) + (qlg(k)+qrg(k))*wlhf/cpm
          end if
          dqcdt(np,k) = (qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
          dqidt(np,k) = d_zero
          dqrdt(np,k) = (qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
          dqsdt(np,k) = d_zero
        else if ( f_qi ) then
          !
          ! If mixed phase hydrometeors are allowed, feed back convective
          ! tendencies of hydrometeors directly.
          !
          dqcdt(np,k) = (qlg(k)-ql0(k))/timec
          dqidt(np,k) = (qig(k)-qi0(k))/timec
          dqrdt(np,k) = (qrg(k)-qr0(k))/timec
          if ( f_qs ) then
            dqsdt(np,k) = (qsg(k)-qs0(k))/timec
          else
            dqidt(np,k) = dqidt(np,k) + (qsg(k)-qs0(k))/timec
          end if
        else
          call fatal(__FILE__,__LINE__, &
              'KAIN-FRITSCH, THIS COMBINATION OF IMOIST,&
             & IEXICE, IICE NOT ALLOWED')
        end if
        dtdt(np,k) = (tg(k)-t0(np,k))/timec
        dqdt(np,k) = (qg(k)-q0(k))/timec
      end do
      pratec(np) = pptflx*(d_one-fbfrc)/dxsq
      ktop(np) = ltop
      kbot(np) = lcl
      if ( istop == 1 ) then
        call fatal(__FILE__,__LINE__,'MODEL STOPS')
      end if
    end do modelpoints

1035 format(1X,'PEF(WS)=',F5.2,'(CB)=',F5.2,'LC,LET=',2I3,'WKL=', &
            F6.3,'VWS=',F5.2)
1080 format(2X,'LFS,LDB,LDT =',3I3,' TIMEC, TADVEC, NSTEP=',      &
            2(1X,F5.0),I3,'NCOUNT, FABE, AINC=',I2,1X,F6.3,F6.2)
1085 format(A3,16A7,2A8)
1110 format(' ','INITIAL WATER =',E12.5,' FINAL WATER =',E12.5,   &
            ' TOTAL WATER CHANGE =',F8.2,'%')

    contains

    subroutine tpmix2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic)
      implicit none
      real(rk8) , intent(in) :: p , thes
      real(rk8) , intent(out) :: qnewlq , qnewic
      real(rk8) , intent(inout) :: tu , qu , qliq , qice
      real(rk8) :: tp , qq , bth , tth , pp
      real(rk8) :: t00 , t10 , t01 , t11
      real(rk8) :: q00 , q10 , q01 , q11
      real(rk8) :: temp , qs , qnew , dq , qtot , rll , cpp
      integer(ik4) :: iptb , ithtb
      !
      !*************************************
      ! scaling pressure and tt table index
      !*************************************
      !
      tp = (p-plutop) * rdpr
      qq = tp - dint(tp)
      iptb = max(1, min(kfnp-1,int(tp)+1))
      !
      !*********************************
      ! base and scaling factor for the
      ! scaling the and tt table index
      !*********************************
      !
      bth = (the0k(iptb+1) - the0k(iptb))*qq+the0k(iptb)
      tth = (thes-bth) * rdthk
      pp = tth - aint(tth)
      ithtb = max(1, min(kfnt-1,int(tth)+1))

      t00 = ttab(ithtb  ,iptb  )
      t10 = ttab(ithtb+1,iptb  )
      t01 = ttab(ithtb  ,iptb+1)
      t11 = ttab(ithtb+1,iptb+1)

      q00 = qstab(ithtb  ,iptb  )
      q10 = qstab(ithtb+1,iptb  )
      q01 = qstab(ithtb  ,iptb+1)
      q11 = qstab(ithtb+1,iptb+1)
      !
      !********************
      ! parcel temperature
      !********************
      !
      temp = (t00 + (t10-t00)*pp + (t01-t00)*qq + (t00-t10-t01+t11)*pp*qq)
      qs = (q00 + (q10-q00)*pp + (q01-q00)*qq + (q00-q10-q01+q11)*pp*qq)
      dq = qs-qu
      if ( dq <= d_zero ) then
        qnew = qu-qs
        qu = qs
      else
        !
        ! If the parcel is subsaturated, temperature and mixing ratio must be
        ! adjusted.  If liquid water is present, it is allowed to evaporate
        !
        qnew = d_zero
        qtot = qliq + qice
        !
        ! if there is enough liquid or ice to saturate the parcel, temp stays
        ! at its wet bulb value, vapor mixing ratio is at saturated level, and
        ! the mixing ratios of liquid and ice are adjusted to make up the
        ! original saturation deficit
        ! Otherwise, any available liq or ice vaporizes and appropriate
        ! adjustments to parcel temp; vapor, liquid, and ice mixing ratios
        ! are made.
        !
        ! Subsaturated values only occur in calculations involving various
        ! mixtures of updraft and environmental air for estimation of
        ! entrainment and detrainment.
        ! For these purposes, assume that reasonable estimates can be given
        ! using liquid water saturation calculations only - i.e., ignore the
        ! effect of the ice phase in this process only
        ! will not affect conservative properties
        !
        if ( qtot >= dq ) then
          qliq = qliq - dq*qliq/(qtot+1.0D-10)
          qice = qice - dq*qice/(qtot+1.0D-10)
          qu = qs
        else
          rll = xlv0-xlv1*temp
          cpp = cpd*(d_one+0.89D0*qu)
          if ( qtot < 1.0D-10 ) then
            ! If no liquid water or ice is available, temperature is given by:
            temp = temp + rll*(dq/(d_one+dq))/cpp
          else
            !
            ! If some liq water/ice is available, but not enough to achieve
            ! saturation, the temperature is given by:
            !
            temp = temp + rll*((dq-qtot)/(d_one+dq-qtot))/cpp
            qu = qu + qtot
            qtot = d_zero
            qliq = d_zero
            qice = d_zero
          end if
        end if
      end if
      tu = temp
      qnewlq = qnew
      qnewic = d_zero
    end subroutine tpmix2

    subroutine dtfrznew(tu,p,thteu,qu,qfrz,qice)
      implicit none
      real(rk8) , intent(in) :: p , qfrz
      real(rk8) , intent(inout) :: tu , thteu , qu , qice
      real(rk8) :: rlc , rls , rlf , cpp , a , dtfrz , es , qs , dqevap , pii
      !
      ! Allow the freezing of liquid water in the updraft to proceed as an
      ! approximately linear function of temperature in the temperature range
      ! ttfrz to tbfrz
      ! For colder temperatures, freeze all liquid water
      ! Thermodynamic properties are still calculated with respect to liquid
      ! water to allow the use of lookup table to extract tmp from thetae
      !
      rlc = 2.5D6 - 2369.276D0*(tu-t00)
      rls = 2833922.0D0 - 259.532D0*(tu-t00)
      rlf = rls - rlc
      cpp = cpd*(d_one + 0.89D0*qu)
      !
      ! a = d(es)/dt is that calculated from Buck (1981) emperical formulas
      ! for saturation vapor pressure
      !
      a = (cliq-bliq*dliq)/((tu-dliq)*(tu-dliq))
      dtfrz = rlf*qfrz/(cpp+rls*qu*a)
      tu = tu + dtfrz

      es = aliq*exp((bliq*tu-cliq)/(tu-dliq))
      qs = es*ep2/(p-es)
      !
      ! Freezing warms the air and it becomes unsaturated
      ! Assume that some of the liquid water that is available for freezing
      ! evaporates to maintain saturation
      ! Since this water has already been transferred to the ice category,
      ! subtract it from ice concentration, THEN set updraft mixing ratio
      ! at the new temperature to the saturation value
      !
      dqevap = qs - qu
      qice = qice - dqevap
      qu = qu + dqevap
      pii = (p00/p)**(0.2854D0*(d_one-0.28D0*qu))
      thteu = tu*pii*exp((c1/tu-c2)*qu*(d_one+c4*qu))
    end subroutine dtfrznew
    !
    ! 9/18/88...this precipitation fallout scheme is based on the scheme us
    ! by Ogura and Cho (1973).  Liquid water fallout from a parcel is cal-
    ! culated using the equation dq=-rate*q*dt, but to simulate a quasi-
    ! continuous process, and to eliminate a dependency on vertical
    ! resolution this is expressed as q=q*exp(-rate*dz).
    !
    subroutine condload(qliq,qice,wtw,dz,boterm,enterm,qnewlq, &
                        qnewic,qlqout,qicout)
      implicit none
      real(rk8) , intent(in) :: dz , boterm , enterm
      real(rk8) , intent(inout) :: qlqout , qicout , wtw
      real(rk8) , intent(inout) :: qliq , qice , qnewlq , qnewic
      real(rk8) :: qtot , qnew , qest , g1 , wavg , conv , ratio3
      real(rk8) :: oldq , ratio4 , dq , pptdrg

      qtot = qliq + qice
      qnew = qnewlq + qnewic
      !
      ! Estimate the vertical velocity so that an average vertical velocity
      ! be calculated to estimate the time required for ascent between model
      ! levels...
      !
      qest = d_half*(qtot+qnew)
      g1 = wtw + boterm - enterm - d_two*egrav*dz*qest/1.5D0
      if ( g1 < d_zero ) g1 = d_zero
      wavg = d_half * (sqrt(wtw) + sqrt(g1))
      conv = kf_entrate * dz/wavg ! KF90  Eq. 9
      !
      ! ratio3 is the fraction of liquid water in fresh condensate, ratio4 is
      ! the fraction of liquid water in the total amount of condensate involv
      ! in the precipitation process - note that only 60% of the fresh conden
      ! sate is is allowed to participate in the conversion process...
      !
      ratio3 = qnewlq / (qnew+1.0D-8)
      qtot = qtot + 0.6D0*qnew
      oldq = qtot
      ratio4 = (0.6D0*qnewlq + qliq) / (qtot+1.0D-8)
      if ( conv > 25.0D0 ) then
        qtot = d_zero
      else
        qtot = qtot * exp(-conv) ! KF90  Eq. 9
      end if
      !
      ! Determine the amount of precipitation that falls out of the updraft
      ! parcel at this level...
      !
      dq = oldq - qtot
      qlqout = ratio4 * dq
      qicout = (d_one-ratio4) * dq
      !
      ! Estimate the mean load of condensate on the updraft in the layer, cal
      ! late vertical velocity
      !
      pptdrg = d_half * (oldq+qtot-0.2D0*qnew)
      wtw = wtw + boterm - enterm - d_two*egrav*dz*pptdrg/1.5D0
      if ( abs(wtw) < 1.0D-4) wtw = 1.0D-4
      !
      ! Determine the new liquid water and ice concentrations including losses
      ! due to precipitation and gains from condensation...
      !
      qliq = ratio4*qtot + ratio3*0.4D0*qnew
      qice = (d_one-ratio4)*qtot + (d_one-ratio3)*0.4D0*qnew
      qnewlq = d_zero
      qnewic = d_zero
    end subroutine condload
    !
    ! GAUSSIAN TYPE MIXING PROFILE
    !
    ! THIS SUBROUTINE INTEGRATES THE AREA UNDER THE CURVE IN THE GAUSSIAN
    ! DISTRIBUTION...THE NUMERICAL APPROXIMATION TO THE INTEGRAL IS TAKEN FROM
    ! "HANDBOOK OF MATHEMATICAL FUNCTIONS WITH FORMULAS, GRAPHS AND
    !  MATHEMATICS TABLES"
    ! ED. BY ABRAMOWITZ AND STEGUN, NATL BUREAU OF STANDARDS APPLIED
    ! MATHEMATICS SERIES.  JUNE, 1964., MAY, 1968.
    !                        JACK KAIN
    !                          7/6/89
    ! Solves for KF90 Eq. 2
    !
    subroutine prof5(eq,ee,ud)
      implicit none
      real(rk8) , intent(in) :: eq
      real(rk8) , intent(inout) :: ee , ud
      real(rk8) :: x , y , ey , e45 , t1 , t2 , c1 , c2

      real(rk8) , parameter :: sqrt2p = 2.506628D0
      real(rk8) , parameter :: a1 = 0.4361836D0
      real(rk8) , parameter :: a2 = -0.1201676D0
      real(rk8) , parameter :: a3 = 0.9372980D0
      real(rk8) , parameter :: p = 0.33267D0
      real(rk8) , parameter :: sigma = 0.166666667D0
      real(rk8) , parameter :: fe = 0.202765151D0
      x = (eq-d_half)/sigma
      y = 6.0D0*eq - 3.0D0
      ey = exp(y*y/(-d_two))
      e45 = exp(-4.5D0)
      t2 = d_one / (d_one + p*abs(y))
      t1 = 0.500498D0
      c1 = a1*t1 + a2*t1*t1 + a3*t1*t1*t1
      c2 = a1*t2 + a2*t2*t2 + a3*t2*t2*t2
      if ( y >= d_zero ) then
        ee = sigma * (d_half*(sqrt2p-e45*c1-ey*c2) + &
                sigma*(e45-ey))-e45*eq*eq*d_half
        ud = sigma * (d_half*(ey*c2-e45*c1) + &
                sigma*(e45-ey)) - e45*(d_half+eq*eq*d_half - eq)
      else
        ee = sigma * (d_half*(ey*c2-e45*c1) + sigma*(e45-ey)) - e45*eq*eq*d_half
        ud = sigma * (d_half*(sqrt2p-e45*c1-ey*c2) + &
                sigma*(e45-ey)) - e45*(d_half+eq*eq*d_half - eq)
      end if
      ee = ee/fe
      ud = ud/fe
    end subroutine prof5

    subroutine tpmix2dd(p,thes,ts,qs)
      implicit none
      real(rk8) , intent(in) :: p , thes
      real(rk8) , intent(inout) :: ts , qs
      real(rk8) :: tp , qq , bth , tth , pp
      real(rk8) :: t00 , t10 , t01 , t11
      real(rk8) :: q00 , q10 , q01 , q11
      integer(ik4) :: iptb , ithtb
      !***************************************************************
      !
      !*************************************
      ! scaling pressure and tt table index
      !*************************************
      !
      tp = (p-plutop) * rdpr
      qq = tp - aint(tp)
      iptb = max(1, min(kfnp-1,int(tp)+1))
      !
      !**********************************
      ! base and scaling factor for the
      !**********************************
      !
      ! scaling the and tt table index
      !
      bth = (the0k(iptb+1)-the0k(iptb)) * qq + the0k(iptb)
      tth = (thes-bth) * rdthk
      pp = tth - aint(tth)
      ithtb = max(1, min(kfnt-1,int(tth)+1))

      t00 = ttab(ithtb  ,iptb  )
      t10 = ttab(ithtb+1,iptb  )
      t01 = ttab(ithtb  ,iptb+1)
      t11 = ttab(ithtb+1,iptb+1)

      q00 = qstab(ithtb  ,iptb  )
      q10 = qstab(ithtb+1,iptb  )
      q01 = qstab(ithtb  ,iptb+1)
      q11 = qstab(ithtb+1,iptb+1)
      !
      !************************************************
      ! parcel temperature and saturation mixing ratio
      !************************************************
      !
      ts = (t00 + (t10-t00)*pp + (t01-t00)*qq + (t00-t10-t01+t11)*pp*qq)
      qs = (q00 + (q10-q00)*pp + (q01-q00)*qq + (q00-q10-q01+q11)*pp*qq)
    end subroutine tpmix2dd

    pure real(rk8) function envirtht(p1,t1,q1) result(tht1)
      implicit none
      real(rk8) , intent(in) :: p1 , t1 , q1
      real(rk8) :: ee , tlog , a1 , tp , avalue , aintrp
      real(rk8) :: tdpt , tsat , tht
      integer(ik4) :: indlu
      !
      !  Calculate environmental equivalent potential temperature
      !
      ! NOTE: Calculations for mixed/ice phase no longer used. jsk 8/00
      !        For example, KF90 Eq. 10 no longer used
      !
      ee = q1*p1 / (ep2+q1)
      !
      ! Calculate log term using lookup table.
      !
      a1 = ee/aliq
      tp = (a1-astrt)/aincb
      indlu = max(1, min(kfna-1,int(tp)+1))
      avalue = (indlu-1)*aincb + astrt
      aintrp = (a1-avalue)/aincb
      tlog = aintrp*alu(indlu+1) + (1-aintrp)*alu(indlu)
      tdpt = (cliq-dliq*tlog)/(bliq-tlog)
      tsat = tdpt - (0.212D0+1.571D-3*(tdpt-t00)-4.36D-4*(t1-t00))*(t1-tdpt)
      tht = t1*(p00/p1)**(0.2854D0*(d_one-0.28D0*q1))
      tht1 = tht*exp((c1/tsat-c2)*q1*(d_one+c4*q1))
    end function envirtht

  end subroutine kfpara
  !
  ! This subroutine is a lookup table.
  ! Given a series of series of saturation equivalent potential
  ! temperatures, the temperature is calculated.
  !
  subroutine kf_lutab
    implicit none
    integer(ik4) :: kp , it , itcnt , i
    real(rk8) :: dpr , temp , p , es , qs , pi , thes , tgues , &
            thgues , f0 , t1 , t0 , f1 , dtx , a1 , thtgs

    ! minimum starting temp
    real(rk8) , parameter :: tmin = 150.0D0
    ! maximum bottom pressure (pascals)
    real(rk8) , parameter :: pbot = 1.1D5
    ! equivalent potential temperature increment
    real(rk8) , parameter :: dth = 1.0D0
    ! tolerance for accuracy of temperature
    real(rk8) , parameter :: toler = 0.001D0

    ! top pressure (pascals)
    plutop = max(ptop*d_1000,5000.0D0)
    ! pressure increment
    dpr = (pbot-plutop) / dble(kfnp-1)
    !
    ! compute parameters
    !
    ! 1._over_(sat. equiv. theta increment)
    rdthk = d_one / dth
    ! 1._over_(pressure increment)
    rdpr = d_one / dpr
    !
    ! calculate the starting sat. equiv. theta
    !
    temp = tmin
    p = plutop - dpr
    do kp = 1 , kfnp
      p = p + dpr
      es = aliq*exp((bliq*temp-cliq)/(temp-dliq))
      qs = ep2*es/(p-es)
      pi = (p00/p)**(0.2854D0*(d_one-0.28D0*qs))
      the0k(kp) = temp * pi * exp((c1/temp-c2)*qs*(d_one+c4*qs))
    end do
    !
    ! compute temperatures for each sat. equiv. potential temp.
    !
    p = plutop - dpr
    do kp = 1 , kfnp
      thes = the0k(kp) - dth
      p = p + dpr
      do it = 1 , kfnt
        ! define sat. equiv. pot. temp.
        thes = thes + dth
        ! iterate to find temperature
        ! find initial guess
        if ( it == 1 ) then
          tgues = tmin
        else
          tgues = ttab(it-1,kp)
        end if
        es = aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
        qs = ep2*es/(p-es)
        pi = (p00/p)**(0.2854D0*(d_one-0.28D0*qs))
        thgues = tgues * pi * exp((c1/tgues-c2)*qs*(d_one+c4*qs))
        f0 = thgues - thes
        t1 = tgues - d_half*f0
        t0 = tgues
        iter1: &
        do itcnt = 1 , 11
          es = aliq * exp((bliq*t1-cliq)/(t1-dliq))
          qs = ep2*es/(p-es)
          pi = (p00/p)**(0.2854D0*(d_one-0.28D0*qs))
          thtgs = t1 * pi * exp((c1/t1-c2)*qs*(d_one+c4*qs))
          f1 = thtgs - thes
          if ( abs(f1) < toler ) then
            exit iter1
          end if
          dtx = f1 * (t1-t0)/(f1-f0)
          t0 = t1
          f0 = f1
          t1 = t1 - dtx
        end do iter1
        ttab(it,kp) = t1
        qstab(it,kp) = qs
      end do
    end do
    !
    ! lookup table for tlog(emix/aliq)
    !
    ! set up intial values for lookup tables
    !
    a1 = astrt - aincb
    do i = 1 , kfna
      a1 = a1 + aincb
      alu(i) = log(a1)
    end do
  end subroutine kf_lutab

end module mod_cu_kf
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
