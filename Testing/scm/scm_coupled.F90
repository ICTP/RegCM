!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Single column test of the Holtslag boundary layer (ibltyp = 1) and the
! Tiedtke cumulus (icup = 5) schemes integrated together in time.
!
! Usage: scm_coupled [namelist file]       (default: scm_coupled.in)
!
! The columns of the profile file are laid out as a 1 x ncol grid. Every
! time step both schemes are called through their model entry points,
! holtbl and tiedtkedrv, on the same state, as the MOLOCH core does
! (idynamic = 3, no chemistry), and the state is advanced with the sum of
! their tendencies and of the large scale forcing:
!   - the input tendencies of the profile file (tten, qten: radiative
!     heating and horizontal advection), held fixed;
!   - the vertical advection by the omega of the profile file, computed
!     on the current state, which is also the omega seen by the cumulus
!     scheme;
!   - a relaxation of temperature and humidity with time scale taurelax
!     (default one day), above zrelax (linear ramp from zrelax/2; zero,
!     the default, relaxes all levels), to the initial profile or, if
!     tgtfile is given (default scm_amazon_target.txt), to the
!     profiles of that file (same format, for each column four profiles
!     at 0, 6, 12 and 18 UTC, linearly interpolated in time, repeated
!     every day).
! The surface fluxes follow a half sine diurnal cycle in local solar time
! (daylight from 6 to 18), with the daily means of the profile file:
! at night the sensible heat flux is -0.25 and the latent heat flux 0.1
! times the daily mean, the daytime half sine keeps the mean. There
! is no microphysics: the condensate detrained by the cumulus scheme
! evaporates in the layer where it is detrained, and supersaturation is
! removed by a saturation adjustment that rains out (large scale rain).
! The pressure of the levels is held fixed.
!
! The namelist has a &scmparam group and the model groups &tiedtkeparam,
! &cldparam (kfac_deep, kfac_shal, k2_const, iconvcf) and &holtslagparam,
! with the RegCM defaults. The output file is plain text, read by
! scm_coupled.py: a header with the settings used, then every nout
! steps per column a STEP line with the scalars and one line per level,
! top to bottom.

program scm_coupled
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_cu_common
  use mod_cu_tiedtke, only : allocate_mod_cu_tiedtke, tiedtkedrv, &
                             pmean, ktype, pmfu
  use mod_pbl_common, only : ricr, kmxpbl
  use mod_pbl_holtbl, only : allocate_mod_pbl_holtbl, holtbl
  implicit none

  character(len=256) :: nmlfile, infile, outfile, tgtfile
  real(rkx) :: dxkm, dtstep, ndays, taurelax, zrelax
  integer(ik4) :: nout
  integer(ik4) :: ncol, nk, j, k, n, iu, io, ios, istep, nsteps
  character(len=32), allocatable, dimension(:) :: col_name
  character(len=10), allocatable, dimension(:) :: col_date
  real(rkx), allocatable, dimension(:) :: col_lat, col_lon, col_ps, &
    col_shf, col_lhf, col_zpbl, col_cape, col_crr, col_topo, col_z0,  &
    col_hour, rainc, rainl, rainc_acc, rainl_acc, ustar
  integer(ik4), allocatable, dimension(:) :: col_land
  real(rkx), allocatable, dimension(:,:) :: t0, q0, u0, v0, omega, &
    ften, fqen, za, zf, pa, pf, t, q, u, v, mfu, ctt, cqt, ptt, pqt, &
    wrelax
  real(rkx), allocatable, dimension(:,:,:) :: ttgt, qtgt
  type(mod_2_cum) :: m2c
  type(mod_2_pbl) :: m2p
  type(pbl_2_mod) :: p2m
  real(rkx) :: dp, hour, utc, wt, tg, qg, sfac, rho, kbfs, hqfac, &
    zpblj, qs, dqs, cond, ws, zl, dtdp, dqdp, time
  character(len=512) :: line

  namelist /scmparam/ infile, outfile, dtstep, ndays, nout, dxkm, &
    taurelax, zrelax, tgtfile
  namelist /tiedtkeparam/ entrdd, entrpen_lnd, entrpen_ocn, &
    detrpen_lnd, detrpen_ocn, entshalp, rcuc_lnd, rcuc_ocn, &
    rcpec_lnd, rcpec_ocn, rhebc_lnd, rhebc_ocn, rprc_ocn, rprc_lnd, &
    revap_lnd, revap_ocn, lmfpen, lmfmid, lmfdd, lepcld, lmfdudv,   &
    lmfscv, lmfuvdis, lmftrac, lmfsmooth, lmfwstar, lmfcapebl, &
    lmfcapeforc, rtaucape
  namelist /cldparam/ kfac_deep, kfac_shal, k2_const, iconvcf
  namelist /holtslagparam/ ricr_ocn, ricr_lnd, zhnew_fac, &
    ifaholtth10, ifaholt, holtth10iter, lmoistn2, maxprandtl

  ! Defaults: as in Main/mod_params.F90
  infile = 'scm_amazon.txt'
  outfile = 'scm_coupled.out'
  tgtfile = 'scm_amazon_target.txt'
  dtstep = 300.0_rkx
  ndays = 5.0_rkx
  nout = 12
  dxkm = 12.0_rkx
  taurelax = 86400.0_rkx
  zrelax = 0.0_rkx
  entrpen_lnd = 1.75e-3_rkx
  entrpen_ocn = 1.75e-3_rkx
  detrpen_lnd = 0.75e-4_rkx
  detrpen_ocn = 0.75e-4_rkx
  entrdd = 3.0e-4_rkx
  entshalp = 2.0_rkx
  rcuc_lnd = 0.05_rkx
  rcuc_ocn = 0.05_rkx
  rcpec_lnd = 5.55e-5_rkx
  rcpec_ocn = 5.55e-5_rkx
  rhebc_lnd = 0.8_rkx
  rhebc_ocn = 0.8_rkx
  rprc_lnd = 1.4e-3_rkx
  rprc_ocn = 1.4e-3_rkx
  revap_lnd = 1.0e-5_rkx
  revap_ocn = 1.0e-5_rkx
  lmfpen = .true.
  lmfmid = .true.
  lmfdd = .true.
  lepcld = .true.
  lmfdudv = .true.
  lmfscv = .true.
  lmfuvdis = .true.
  lmftrac = .true.
  lmfsmooth = .false.
  lmfwstar = .false.
  lmfcapebl = .true.
  lmfcapeforc = .false.
  rtaucape = 43200.0_rkx
  k2_const = 500.0_rkx
  kfac_shal = 0.07_rkx
  kfac_deep = 0.14_rkx
  iconvcf = 0
  ricr_ocn = 0.275_rkx
  ricr_lnd = 0.250_rkx
  zhnew_fac = 0.25_rkx
  ifaholt = 0
  ifaholtth10 = 2
  holtth10iter = 1
  lmoistn2 = .false.
  maxprandtl = 1.0_rkx

  nmlfile = 'scm_coupled.in'
  if ( command_argument_count( ) > 0 ) call get_command_argument(1,nmlfile)
  open(newunit=iu,file=nmlfile,status='old',action='read',iostat=ios)
  if ( ios == 0 ) then
    read(iu,nml=scmparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading scmparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=tiedtkeparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading tiedtkeparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=cldparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading cldparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=holtslagparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading holtslagparam in '//trim(nmlfile))
    close(iu)
  else if ( command_argument_count( ) > 0 ) then
    call scm_error('cannot open '//trim(nmlfile))
  end if
  nout = max(nout,1)

  ! Model configuration seen by the schemes
  idynamic = 3
  icup_lnd = 5
  icup_ocn = 5
  icup(:) = 5
  ibltyp = 1
  ipptls = 2
  nqx = 5
  ichem = 0
  iaerosol = 0
  iindirect = 0
  icumcloud = 1
  ntr = 0
  dx = dxkm * 1000.0_rkx
  dt = dtstep
  rdt = d_one/dt
  dtcum = dtstep
  nsteps = nint(ndays*86400.0_rkx/dtstep)

  ! Read the column profiles
  open(newunit=iu,file=infile,status='old',action='read',iostat=ios)
  if ( ios /= 0 ) call scm_error('cannot open '//trim(infile))
  do
    read(iu,'(a)') line
    if ( line(1:1) /= '#' ) exit
  end do
  read(line,*) ncol, nk
  kz = nk
  kzp1 = nk + 1
  kzm1 = nk - 1
  call set_grid(ncol)
  call memory_init
  allocate(col_name(ncol), col_date(ncol), col_lat(ncol), col_lon(ncol), &
    col_ps(ncol), col_shf(ncol), col_lhf(ncol), col_zpbl(ncol),           &
    col_cape(ncol), col_crr(ncol), col_topo(ncol), col_land(ncol),        &
    col_z0(ncol), col_hour(ncol), rainc(ncol), rainl(ncol),               &
    rainc_acc(ncol), rainl_acc(ncol), ustar(ncol))
  allocate(t0(ncol,nk), q0(ncol,nk), u0(ncol,nk), v0(ncol,nk),       &
    omega(ncol,nk), ften(ncol,nk), fqen(ncol,nk), za(ncol,nk),       &
    zf(ncol,nk+1), pa(ncol,nk), pf(ncol,nk+1), t(ncol,nk),           &
    q(ncol,nk), u(ncol,nk), v(ncol,nk), mfu(ncol,nk), ctt(ncol,nk),  &
    cqt(ncol,nk), ptt(ncol,nk), pqt(ncol,nk), wrelax(ncol,nk))
  call getmem(cuscheme,jci1,jci2,ici1,ici2,'scm:cuscheme')
  cuscheme(:,:) = 5
  call init_mod_cumulus
  call allocate_mod_cu_tiedtke
  call getmem(ricr,jci1,jci2,ici1,ici2,'scm:ricr')
  call getmem(kmxpbl,jci1,jci2,ici1,ici2,'scm:kmxpbl')
  call allocate_mod_pbl_holtbl
  call allocate_fields
  do j = 1, ncol
    read(iu,*) col_name(j), col_lat(j), col_lon(j), col_land(j), &
      col_date(j), col_ps(j), col_shf(j), col_lhf(j), col_zpbl(j), &
      col_cape(j), col_crr(j), col_topo(j)
    do k = 1, nk
      read(iu,*) za(j,k), pa(j,k), t(j,k), q(j,k), u(j,k), v(j,k), &
        omega(j,k), ften(j,k), fqen(j,k)
    end do
    do k = 1, nk + 1
      read(iu,*) zf(j,k), pf(j,k)
    end do
    ! specific humidity to mixing ratio, as in the model
    do k = 1, nk
      fqen(j,k) = fqen(j,k)/(d_one-q(j,k))**2
      q(j,k) = q(j,k)/(d_one-q(j,k))
    end do
    read(col_date(j)(9:10),*) k
    ! local solar time of the start
    col_hour(j) = modulo(real(k,rkx) + col_lon(j)/15.0_rkx,24.0_rkx)
    col_z0(j) = merge(0.1_rkx,1.0e-3_rkx,col_land(j) == 1)
    ricr(j,1) = merge(ricr_lnd,ricr_ocn,col_land(j) == 1)
    kmxpbl(j,1) = kz
    do k = kzm1, 2, -1
      if ( za(j,k) >= 5500.0_rkx ) exit
      kmxpbl(j,1) = k
    end do
    do k = 1, nk
      if ( zrelax > d_zero ) then
        wrelax(j,k) = min(max((za(j,k)-d_half*zrelax)/(d_half*zrelax), &
          d_zero),d_one)/taurelax
      else
        wrelax(j,k) = d_one/taurelax
      end if
    end do
    m2p%ldmsk(j,1) = col_land(j)
    m2c%ldmsk(j,1) = col_land(j)
    m2p%coriol(j,1) = d_two*eomeg*sin(col_lat(j)*degrad)
    m2c%psb(j,1) = col_ps(j)
  end do
  close(iu)
  t0 = t
  q0 = q
  allocate(ttgt(ncol,nk,0:4), qtgt(ncol,nk,0:4))
  do n = 0, 4
    ttgt(:,:,n) = t0
    qtgt(:,:,n) = q0
  end do
  if ( len_trim(tgtfile) > 0 ) call read_targets
  u0 = u
  v0 = v
  do k = 1, nk
    pmean(k) = sum(pa(:,k))/real(ncol,rkx)
  end do
  do j = 1, ncol
    m2c%zas(j,1,:) = za(j,:)
    m2c%zfs(j,1,:) = zf(j,:)
    m2c%pas(j,1,:) = pa(j,:)
    m2c%pasf(j,1,:) = pf(j,:)
    m2p%za(j,1,:) = za(j,:)
    m2p%patm(j,1,:) = pa(j,:)
    m2p%patmf(j,1,:) = pf(j,:)
  end do
  p2m%zpbl(:,:) = 500.0_rkx
  rainc_acc = d_zero
  rainl_acc = d_zero

  open(newunit=io,file=outfile,status='replace',action='write')
  write(io,'(a)') '# RegCM Holtslag + Tiedtke coupled single column test'
  write(io,'(a,a)') '# profiles: ', trim(infile)
  write(io,'(a,a)') '# targets: ', trim(tgtfile)
  write(io,'(a,f8.2,a,f8.3,a,f10.1,a,f8.1,a,f8.3)') '# dt (s) ', dt, &
    '  ndays ', ndays, '  taurelax (s) ', taurelax, '  zrelax (m) ', &
    zrelax, '  dx (km) ', dxkm
  write(io,nml=tiedtkeparam)
  write(io,nml=cldparam)
  write(io,nml=holtslagparam)
  write(io,'(a)') '# STEP name time(h) local_hour shf lhf(W/m2) '// &
    'zpbl(m) kpbl hqfac ktype kctop kcbot rainc rainl(mm/day) '// &
    'rainc_acc rainl_acc(mm)'
  write(io,'(a)') '# z(m) p(hPa) T(K) q(g/kg) mfu(kg/m2/s) '// &
    'dTcu dTpbl(K/day) dqcu dqpbl(g/kg/day)'
  write(io,'(a,3i6)') 'NCOL ', ncol, nk, nsteps/nout

  do istep = 1, nsteps
    time = real(istep-1,rkx)*dt
    ! Surface fluxes
    do j = 1, ncol
      hour = modulo(col_hour(j) + time/3600.0_rkx,24.0_rkx)
      sfac = sin(mathpi*(hour-6.0_rkx)/12.0_rkx)
      if ( col_land(j) /= 1 ) then
        m2p%hfx(j,1) = col_shf(j)
        m2p%qfx(j,1) = col_lhf(j)/wlhv
      else if ( sfac > d_zero ) then
        m2p%hfx(j,1) = col_shf(j)*mathpi*1.125_rkx*sfac
        m2p%qfx(j,1) = col_lhf(j)*mathpi*0.95_rkx*sfac/wlhv
      else
        m2p%hfx(j,1) = -0.25_rkx*col_shf(j)
        m2p%qfx(j,1) = 0.1_rkx*col_lhf(j)/wlhv
      end if
    end do
    call update_derived

    ! Cumulus, on the state at the start of the step
    do j = 1, ncol
      m2c%tas(j,1,:) = t(j,:)
      m2c%qxas(j,1,:,iqv) = q(j,:)
      m2c%uas(j,1,:) = u(j,:)
      m2c%vas(j,1,:) = v(j,:)
      m2c%hfx(j,1) = m2p%hfx(j,1)
      m2c%qfx(j,1) = m2p%qfx(j,1)
      m2c%zpbl(j,1) = p2m%zpbl(j,1)
      m2c%ustar(j,1) = ustar(j)
    end do
    m2c%tten(:,:,:) = d_zero
    m2c%qxten(:,:,:,:) = d_zero
    avg_ww(:,1,:) = omega
    cu_prate(:,:) = d_zero
    cu_srate(:,:) = d_zero
    cu_ktop(:,:) = 0
    cu_kbot(:,:) = 0
    cu_tten(:,:,:) = d_zero
    cu_qten(:,:,:,:) = d_zero
    cu_uten(:,:,:) = d_zero
    cu_vten(:,:,:) = d_zero
    cu_cldfrc(:,:,:) = d_zero
    cu_cldlwc(:,:,:) = d_zero
    cu_qdetr(:,:,:) = d_zero
    cu_raincc(:,:,:) = d_zero
    total_precip_points = 0
    call tiedtkedrv(m2c)

    ! Boundary layer, on the same state
    p2m%tten(:,:,:) = d_zero
    p2m%uten(:,:,:) = d_zero
    p2m%vten(:,:,:) = d_zero
    p2m%qxten(:,:,:,:) = d_zero
    call holtbl(m2p,p2m)

    ! Advance the state
    do j = 1, ncol
      utc = modulo(col_hour(j) - col_lon(j)/15.0_rkx + &
        time/3600.0_rkx,24.0_rkx)
      n = min(int(utc/6.0_rkx),3)
      wt = (utc - 6.0_rkx*real(n,rkx))/6.0_rkx
      do k = 1, nk
        ! detrained condensate evaporates in the layer
        ctt(j,k) = cu_tten(j,1,k) - &
          (wlhv*cu_qten(j,1,k,iqc) + wlhs*cu_qten(j,1,k,iqi))*rcpd
        cqt(j,k) = cu_qten(j,1,k,iqv) + cu_qten(j,1,k,iqc) + &
          cu_qten(j,1,k,iqi)
        ptt(j,k) = p2m%tten(j,1,k)
        pqt(j,k) = p2m%qxten(j,1,k,iqv)
        mfu(j,k) = pmfu(j,k)
      end do
      do k = 1, nk
        ! vertical advection by omega, centred differences
        if ( k == 1 ) then
          dtdp = (t(j,2)-t(j,1))/(pa(j,2)-pa(j,1))
          dqdp = (q(j,2)-q(j,1))/(pa(j,2)-pa(j,1))
        else if ( k == nk ) then
          dtdp = (t(j,nk)-t(j,nk-1))/(pa(j,nk)-pa(j,nk-1))
          dqdp = (q(j,nk)-q(j,nk-1))/(pa(j,nk)-pa(j,nk-1))
        else
          dtdp = (t(j,k+1)-t(j,k-1))/(pa(j,k+1)-pa(j,k-1))
          dqdp = (q(j,k+1)-q(j,k-1))/(pa(j,k+1)-pa(j,k-1))
        end if
        tg = (d_one-wt)*ttgt(j,k,n) + wt*ttgt(j,k,n+1)
        qg = (d_one-wt)*qtgt(j,k,n) + wt*qtgt(j,k,n+1)
        t(j,k) = t(j,k) + dt*(ctt(j,k) + ptt(j,k) + ften(j,k) - &
          omega(j,k)*(dtdp - rovcp*t(j,k)/pa(j,k)) - &
          wrelax(j,k)*(t(j,k)-tg))
        q(j,k) = max(q(j,k) + dt*(cqt(j,k) + pqt(j,k) + fqen(j,k) - &
          omega(j,k)*dqdp - wrelax(j,k)*(q(j,k)-qg)),1.0e-8_rkx)
        u(j,k) = u(j,k) + dt*(cu_uten(j,1,k) + p2m%uten(j,1,k) + &
          m2p%coriol(j,1)*(v(j,k)-v0(j,k)))
        v(j,k) = v(j,k) + dt*(cu_vten(j,1,k) + p2m%vten(j,1,k) - &
          m2p%coriol(j,1)*(u(j,k)-u0(j,k)))
      end do
      ! saturation adjustment, the condensate rains out
      rainl(j) = d_zero
      do k = 1, nk
        qs = pfwsat(t(j,k),pa(j,k))
        if ( q(j,k) > qs ) then
          dqs = pfqsdt(t(j,k),pa(j,k))
          cond = (q(j,k)-qs)/(d_one + wlhv*rcpd*dqs)
          q(j,k) = q(j,k) - cond
          t(j,k) = t(j,k) + wlhv*rcpd*cond
          rainl(j) = rainl(j) + cond*(pf(j,k+1)-pf(j,k))*regrav
        end if
      end do
      rainl(j) = rainl(j)/dt
      rainc(j) = cu_prate(j,1)
      rainc_acc(j) = rainc_acc(j) + rainc(j)*dt
      rainl_acc(j) = rainl_acc(j) + rainl(j)*dt
    end do

    if ( mod(istep,nout) == 0 ) then
      time = real(istep,rkx)*dt
      do j = 1, ncol
        hour = modulo(col_hour(j) + time/3600.0_rkx,24.0_rkx)
        rho = pf(j,nk+1)/(rgas*t(j,nk))
        kbfs = m2p%hfx(j,1)/(rho*cpd) + ep1*t(j,nk)*m2p%qfx(j,1)/rho
        zpblj = p2m%zpbl(j,1)
        if ( kbfs > d_zero ) then
          hqfac = max(0.8_rkx - 8.5_rkx*zhnew_fac*kbfs,0.25_rkx)
        else
          hqfac = 0.8_rkx
        end if
        write(io,'(a,1x,a,2f9.3,2f9.2,f9.1,i4,f7.3,3i4,4es13.5)') 'STEP', &
          trim(col_name(j)), time/3600.0_rkx, hour, m2p%hfx(j,1), &
          m2p%qfx(j,1)*wlhv, zpblj, p2m%kpbl(j,1), hqfac, ktype(j), &
          cu_ktop(j,1), cu_kbot(j,1), rainc(j)*86400.0_rkx, &
          rainl(j)*86400.0_rkx, rainc_acc(j), rainl_acc(j)
        do k = 1, nk
          write(io,'(f10.2,f10.3,f9.3,7es13.5)') za(j,k), &
            pa(j,k)*0.01_rkx, t(j,k), q(j,k)/(d_one+q(j,k))*1.0e3_rkx, &
            mfu(j,k), ctt(j,k)*86400.0_rkx, ptt(j,k)*86400.0_rkx, &
            cqt(j,k)*8.64e7_rkx, pqt(j,k)*8.64e7_rkx
        end do
      end do
    end if
  end do
  close(io)
  write(*,'(a,i4,a,i8,a,a)') ' scm_coupled: ', ncol, ' columns, ', &
    nsteps, ' steps, columns written to ', trim(outfile)

  contains

  subroutine scm_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(a,a)') ' scm_coupled error: ', msg
    stop 1
  end subroutine scm_error

  ! Relaxation targets: per column the profiles at 0, 6, 12, 18 UTC
  subroutine read_targets
    implicit none
    integer(ik4) :: iut, nc, nkt, m, kk
    real(rkx) :: x1, x2, x3, x4
    character(len=32) :: cname
    open(newunit=iut,file=tgtfile,status='old',action='read',iostat=ios)
    if ( ios /= 0 ) call scm_error('cannot open '//trim(tgtfile))
    do
      read(iut,'(a)') line
      if ( line(1:1) /= '#' ) exit
    end do
    read(line,*) nc, nkt
    if ( nc /= 4*ncol .or. nkt /= nk ) &
      call scm_error('target file does not match the profiles')
    do j = 1, ncol
      do m = 0, 3
        read(iut,*) cname
        do kk = 1, nk
          read(iut,*) x1, x2, ttgt(j,kk,m), qtgt(j,kk,m), x3, x4
          qtgt(j,kk,m) = qtgt(j,kk,m)/(d_one-qtgt(j,kk,m))
        end do
        do kk = 1, nk + 1
          read(iut,*) x1, x2
        end do
      end do
      ttgt(j,:,4) = ttgt(j,:,0)
      qtgt(j,:,4) = qtgt(j,:,0)
    end do
    close(iut)
  end subroutine read_targets

  ! A 1 x n grid: every index range of the model domain is 1:n x 1:1
  subroutine set_grid(n)
    implicit none
    integer(ik4), intent(in) :: n
    jci1 = 1 ; jci2 = n ; ici1 = 1 ; ici2 = 1
    jce1 = 1 ; jce2 = n ; ice1 = 1 ; ice2 = 1
    jdi1 = 1 ; jdi2 = n ; idi1 = 1 ; idi2 = 1
    jci1ga = 1 ; jci2ga = n ; ici1ga = 1 ; ici2ga = 1
    jce1ga = 1 ; jce2ga = n ; ice1ga = 1 ; ice2ga = 1
    jce1gb = 1 ; jce2gb = n ; ice1gb = 1 ; ice2gb = 1
  end subroutine set_grid

  ! Boundary layer input derived from the state: potential temperature,
  ! surface air density, ground temperature (same potential temperature
  ! as the lowest level, used only with ifaholtth10 /= 2) and the surface
  ! drag from a neutral log law friction velocity.
  subroutine update_derived
    implicit none
    real(rkx) :: ps, tv, wspd
    do j = 1, ncol
      m2p%tatm(j,1,:) = t(j,:)
      m2p%qxatm(j,1,:,iqv) = q(j,:)
      m2p%uxatm(j,1,:) = u(j,:)
      m2p%vxatm(j,1,:) = v(j,:)
      do k = 1, nk
        m2p%thatm(j,1,k) = t(j,k)*(p00/pa(j,k))**rovcp
      end do
      ps = pf(j,nk+1)
      tv = t(j,nk)*(d_one+ep1*q(j,nk)/(d_one+q(j,nk)))
      m2p%rhox2d(j,1) = ps/(rgas*tv)
      m2p%tg(j,1) = m2p%thatm(j,1,nk)*(ps/p00)**rovcp
      wspd = max(sqrt(u(j,nk)**2+v(j,nk)**2),0.1_rkx)
      zl = max(za(j,nk),d_two*col_z0(j))
      ws = vonkar*wspd/log(zl/col_z0(j))
      ustar(j) = ws
      m2p%uvdrag(j,1) = m2p%rhox2d(j,1)*ws*ws/wspd
    end do
  end subroutine update_derived

  subroutine allocate_fields
    implicit none
    call getmem(m2c%pas,1,ncol,1,1,1,nk,'scm:pas')
    call getmem(m2c%pasf,1,ncol,1,1,1,nk+1,'scm:pasf')
    call getmem(m2c%zas,1,ncol,1,1,1,nk,'scm:zas')
    call getmem(m2c%zfs,1,ncol,1,1,1,nk+1,'scm:zfs')
    call getmem(m2c%tas,1,ncol,1,1,1,nk,'scm:tas')
    call getmem(m2c%uas,1,ncol,1,1,1,nk,'scm:uas')
    call getmem(m2c%vas,1,ncol,1,1,1,nk,'scm:vas')
    call getmem(m2c%qxas,1,ncol,1,1,1,nk,1,nqx,'scm:qxas')
    call getmem(m2c%tten,1,ncol,1,1,1,nk,'scm:ctten')
    call getmem(m2c%uten,1,ncol,1,1,1,nk,'scm:cuten')
    call getmem(m2c%vten,1,ncol,1,1,1,nk,'scm:cvten')
    call getmem(m2c%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:cqxten')
    call getmem(m2c%qfx,1,ncol,1,1,'scm:cqfx')
    call getmem(m2c%hfx,1,ncol,1,1,'scm:chfx')
    call getmem(m2c%zpbl,1,ncol,1,1,'scm:czpbl')
    call getmem(m2c%ustar,1,ncol,1,1,'scm:ustar')
    call getmem(m2c%psb,1,ncol,1,1,'scm:psb')
    call getmem(m2c%ldmsk,1,ncol,1,1,'scm:cldmsk')
    m2c%qxas(:,:,:,:) = d_zero
    m2c%uten(:,:,:) = d_zero
    m2c%vten(:,:,:) = d_zero
    call getmem(m2p%za,1,ncol,1,1,1,nk,'scm:za')
    call getmem(m2p%patm,1,ncol,1,1,1,nk,'scm:patm')
    call getmem(m2p%patmf,1,ncol,1,1,1,nk+1,'scm:patmf')
    call getmem(m2p%tatm,1,ncol,1,1,1,nk,'scm:tatm')
    call getmem(m2p%thatm,1,ncol,1,1,1,nk,'scm:thatm')
    call getmem(m2p%uxatm,1,ncol,1,1,1,nk,'scm:uxatm')
    call getmem(m2p%vxatm,1,ncol,1,1,1,nk,'scm:vxatm')
    call getmem(m2p%qxatm,1,ncol,1,1,1,nk,1,nqx,'scm:qxatm')
    call getmem(m2p%hfx,1,ncol,1,1,'scm:hfx')
    call getmem(m2p%qfx,1,ncol,1,1,'scm:qfx')
    call getmem(m2p%uvdrag,1,ncol,1,1,'scm:uvdrag')
    call getmem(m2p%rhox2d,1,ncol,1,1,'scm:rhox2d')
    call getmem(m2p%tg,1,ncol,1,1,'scm:tg')
    call getmem(m2p%coriol,1,ncol,1,1,'scm:coriol')
    call getmem(m2p%ldmsk,1,ncol,1,1,'scm:ldmsk')
    call getmem(p2m%tten,1,ncol,1,1,1,nk,'scm:tten')
    call getmem(p2m%uten,1,ncol,1,1,1,nk,'scm:uten')
    call getmem(p2m%vten,1,ncol,1,1,1,nk,'scm:vten')
    call getmem(p2m%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:qxten')
    call getmem(p2m%zpbl,1,ncol,1,1,'scm:zpbl')
    call getmem(p2m%kpbl,1,ncol,1,1,'scm:kpbl')
    m2p%qxatm(:,:,:,:) = d_zero
  end subroutine allocate_fields

#include <saturation.inc>

end program scm_coupled

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort
