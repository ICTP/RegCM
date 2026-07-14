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

module mod_bdycod
  !
  ! Subroutines for input of boundary values and tendencies
  ! Relaxation and Sponge Boundary Conditions routines
  !
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_regcm_types
  use mod_mppparam
  use mod_mpmessage
  use mod_memutil
  use mod_atm_interface
  use mod_pbl_interface, only : tkemin
  use mod_che_interface
  use mod_lm_interface
  use mod_mpmessage
  use mod_ncio
  use mod_service
  use mod_zita
  use mod_stdatm
  use mod_slabocean
  use mod_vertint, only : intz1
  use mod_humid, only : rh2mxr

  implicit none

  private

  public :: initideal
  public :: allocate_mod_bdycon, init_bdy, bdyin, bdyval
#ifdef ASYNC_NETCDF
  public :: bdyin_prefetch
#endif
  public :: sponge, nudge, morelax, setup_bdycon, raydamp
  public :: moupdate_norm, mospectral_nudge
  public :: is_present_qc, is_present_qi

  !
  ! West U External  = WUE
  ! West U Internal  = WUI
  ! East U External  = EUE
  ! .....
  ! North V Internal = NVI
  !
  public :: wue, wui, eue, eui
  public :: wve, wvi, eve, evi
  public :: sue, sui, nue, nui
  public :: sve, svi, nve, nvi
  ! fnudge : are the coefficients for the newtonian term.
  ! gnydge : are the coefficients for the diffusion term.
  public :: fnudge, gnudge
  !
  real(rkx), pointer, contiguous, dimension(:,:) :: sue => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: sui => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: nue => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: nui => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: sve => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: svi => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: nve => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: nvi => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: wue => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: wui => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: eue => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: eui => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: wve => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: wvi => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: eve => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: evi => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: psdot => null( )
  real(rkx), pointer, contiguous, dimension(:) :: fcx => null( )
  real(rkx), pointer, contiguous, dimension(:) :: gcx => null( )
  real(rkx), pointer, contiguous, dimension(:) :: fcd => null( )
  real(rkx), pointer, contiguous, dimension(:) :: gcd => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: hefc => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: hegc => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: hefd => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: hegd => null( )
  real(rkx), pointer, contiguous, dimension(:) :: wgtd => null( )
  real(rkx), pointer, contiguous, dimension(:) :: wgtx => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: fg1 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: fg2 => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: zn1 => null( )
  real(rkx), pointer, contiguous, dimension(:) :: cnudge => null( )
  integer(ik4) :: km, lm
  real(rkx), pointer, dimension(:,:), contiguous :: bvx, bvy
  real(rkx), pointer, dimension(:,:), contiguous :: sxg, syg
  real(rkx), pointer, dimension(:,:), contiguous :: sx, sy
  real(rkx) :: fnudge, gnudge, rtb, cn0
  integer(ik4) :: som_month
#ifdef ASYNC_NETCDF
  integer(ik4), save :: bdyin_prefetch_steps = 1
  logical, save :: bdyin_prefetch_steps_initialized = .false.
#endif

  real(rkx), parameter, dimension(10) :: qxbval = &
    [ 1.0e-8_rkx, 0.0_rkx, 0.0_rkx,       &  ! qv, qc, qi
      0.0_rkx, 0.0_rkx, 0.0_rkx, 0.0_rkx, &  ! qr, qs, qg, qh,
      1.0e10_rkx, 100.0_rkx, 0.01_rkx ]      ! ncc, nc, nr

  interface nudge
    module procedure nudge4d
    module procedure nudge4d3d
    module procedure nudge3d
    module procedure nudge2d
    module procedure nudgeuv
  end interface nudge

  interface sponge
    module procedure sponge4d
    module procedure sponge3d
    module procedure sponge2d
    module procedure spongeuv
  end interface sponge

  interface raydamp
    module procedure raydamp3
    module procedure raydamp3f
    module procedure raydampuv
    module procedure raydampqv
    module procedure raydampuv_c
  end interface raydamp

  logical, parameter :: bdyflow = .true.
  logical :: present_qc = .false.
  logical :: present_qi = .false.

  contains

  subroutine initideal
    implicit none
    integer(ik4) :: iunit, ierr
    integer(ik4) :: k
    real(rkx) :: ps, ts
    real(rkx), allocatable, dimension(:) :: zi, ti, qi
#ifndef RCEMIP
    real(rkx), allocatable, dimension(:) :: g, p, t, u, v, r
    real(rkx), allocatable, dimension(:) :: pi, ui, vi
    integer(ik4) :: i, nlev
    integer(ik4) :: nseed
    integer(ik8) :: sclock
    integer(ik4), allocatable, dimension(:) :: seed
    real(rkx) :: ht
    real(rkx), dimension(jce1ga:jce2ga,ice1ga:ice2ga) :: noise

    namelist /dimensions/ nlev
    namelist /surface/ ps, ts
    namelist /profile/ g, p, t, u, v, r

    call random_seed(size=nseed)
    call system_clock(sclock)
    seed = myid + int(sclock) + 37*[(i-1,i=1,nseed)]
    call random_seed(put = seed)
    call random_number(noise)

    open(newunit=iunit, file='profile.in', status='old', &
         action='read', iostat=ierr, err=100)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Open error for profile.in')
    end if

    read(iunit, nml=dimensions, iostat=ierr, err=200)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Read error for namelist dimensions in profile.in')
    end if
    rewind(iunit)

    allocate(g(nlev),p(nlev),t(nlev),u(nlev),v(nlev),r(nlev))
    allocate(zi(kz),pi(kz),ti(kz),ui(kz),vi(kz),qi(kz))

    read(iunit, nml=surface, iostat=ierr, err=300)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Read error for namelist surface in profile.in')
    end if
    rewind(iunit)

    read(iunit, nml=profile, iostat=ierr, err=400)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Read error for namelist profile in profile.in')
    end if

    close(iunit)

    if ( myid == italk ) then
      write(stdout,*) 'Successfully read in initial profile.'
    end if

    g = g * regrav
    r = r*d_r100
    p = p * d_100
    call rh2mxr(t,r,p,nlev)
    p = p * d_r100
    call invert_top_bottom(t)
    call invert_top_bottom(p)
    call invert_top_bottom(g)
    call invert_top_bottom(u)
    call invert_top_bottom(v)
    call invert_top_bottom(r)
    xtsb%b0 = ts
    if ( idynamic == 1 ) then
      xpsb%b0 = ps * 0.1_rkx
      ht = 0.0
      ! vertical interpolation for dub%b0,dvb%b0,xtb%b0,xqb%b0
      pi = sigma*(ps-ptop*10.0_rkx) + (ptop*10.0_rkx)
      call intz1(ui,u,pi,p,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,pi,p,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,pi,p,ht,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,pi,p,ht,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1, kz
        !$acc kernels
        dub%b0(:,:,k) = ui(k)
        dvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
        !$acc end kernels
      end do
    else if ( idynamic == 2 ) then
      ! vertical interpolation for dub%b0,dvb%b0,xtb%b0,xqb%b0
      xpsb%b0 = ps * 0.1_rkx
      ht = 0.0
      zi = atm0%z(jci1,ici1,1:kz)
      call intz1(ui,u,zi,g,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,zi,g,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,zi,g,ht,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,zi,g,ht,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1, kz
        !$acc kernels
        dub%b0(:,:,k) = ui(k)
        dvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
        !$acc end kernels
      end do
      xppb%b0 = 0.0_rkx
      xwwb%b0 = 0.0_rkx
    else if ( idynamic == 3 ) then
      ! vertical interpolation for dub%b0,dvb%b0,xtb%b0,xqb%b0
      xpsb%b0 = ps * 100.0_rkx
      ht = 0.0
      zi = mo_atm%zeta(jci1,ici1,1:kz)
      call intz1(ui,u,zi,g,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,zi,g,ht,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,zi,g,ht,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,zi,g,ht,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1, kz
        !$acc kernels
        dub%b0(:,:,k) = ui(k)
        dvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
        !$acc end kernels
      end do
      call paicompute(xpsb%b0,mo_atm%zeta,xtb%b0,xqb%b0,xpaib%b0)
    else
      call fatal(__FILE__,__LINE__, &
        'Should never get here....')
    end if

    deallocate(seed)
    deallocate(g,p,t,u,v,r)
    deallocate(zi,pi,ti,ui,vi,qi)
#else
    real(rkx) :: qs, tvs, tvt
    real(rkx), parameter, dimension(3) :: t0s = &
                [ 295.0_rkx, 300.0_rkx, 305.0_rkx ]
    real(rkx), parameter, dimension(3) :: q0s = &
                [ 0.012_rkx, 0.01865_rkx, 0.024_rkx ]
    real(rkx), parameter :: p0s = 101480

    namelist /surface/ ps, ts
    namelist /perturbation/ lrcemip_perturb, lrcemip_noise_level

    open(newunit=iunit, file='profile.in', status='old', &
         action='read', iostat=ierr, err=100)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Open error for profile.in')
    end if
    read(iunit, nml=surface, iostat=ierr, err=300)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__, &
                 'Read error for namelist surface in profile.in')
    end if
    rewind(iunit)
    lrcemip_perturb = .false.
    lrcemip_noise_level = 0.0001_rkx
    read(iunit, nml=perturbation, iostat=ierr)
    if ( ierr /= 0 ) then
      if ( myid == italk ) then
        write(stdout, *) 'No perturbation added.'
      end if
    end if
    close(iunit)

    if ( myid == italk ) then
      write(stdout,*) 'Successfully read in surface temperature.'
    end if

    ps = p0s
    if ( ts <= t0s(1) ) then
      ts = t0s(1)
      qs = q0s(1)
    else if ( ts <= t0s(2) ) then
      ts = t0s(2)
      qs = q0s(2)
    else
      ts = t0s(3)
      qs = q0s(3)
    end if
    tvs = ts * (1.0_rkx + ep1*qs)
    tvt = tvs - lrate * 15000.0_rkx

    allocate(zi(kz),ti(kz),qi(kz))

    if ( idynamic == 3 ) then
      zi = mo_atm%zeta(jci1,ici1,1:kz)
    else if ( idynamic == 2 ) then
      zi = atm0%z(jci1,ici1,1:kz)
    else
      call fatal(__FILE__,__LINE__, 'Not implemented...')
    end if

    do k = 1, kz
      if ( zi(k) < 15000.0_rkx ) then
        qi(k) = qs * exp(-zi(k)/4000.0_rkx) * exp((-zi(k)/7500.0_rkx)**2)
        ti(k) = (tvs - lrate * zi(k))/(1.0_rkx + ep1*qi(k))
      else
        qi(k) = minqq
        ti(k) = tvt
      end if
    end do

    xtsb%b0 = ts
    qi = qi/(1.0_rkx-qi)
    do k = 1, kz
      !$acc kernels
      dub%b0(:,:,k) = 0.0_rkx
      dvb%b0(:,:,k) = 0.0_rkx
      xtb%b0(:,:,k) = ti(k)
      xqb%b0(:,:,k) = qi(k)
      !$acc end kernels
    end do

    if ( idynamic == 1 ) then
      !$acc kernels
      xpsb%b0(:,:) = ps * 0.001_rkx
      !$acc end kernels
    else if ( idynamic == 2 ) then
      !$acc kernels
      xpsb%b0(:,:) = ps * 0.001_rkx
      xppb%b0(:,:,:) = 0.0_rkx
      xwwb%b0(:,:,:) = 0.0_rkx
      !$acc end kernels
    else if ( idynamic == 3 ) then
      !$acc kernels
      xpsb%b0(:,:) = ps
      !$acc end kernels
      call paicompute(xpsb%b0,mo_atm%zeta,xtb%b0,xqb%b0,xpaib%b0)
    else
      call fatal(__FILE__,__LINE__, &
        'Should never get here....')
    end if
    deallocate(zi,ti,qi)
#endif

    return

100 call fatal(__FILE__,__LINE__, 'Error opening namelist file  profile.in')
200 call fatal(__FILE__,__LINE__, 'Error reading namelist dimensions')
300 call fatal(__FILE__,__LINE__, 'Error reading namelist surface')
400 call fatal(__FILE__,__LINE__, 'Error reading namelist profile')

  end subroutine initideal

  logical function is_present_qc( )
    implicit none
    is_present_qc = present_qc
  end function is_present_qc

  logical function is_present_qi( )
    implicit none
    is_present_qi = present_qi
  end function is_present_qi

  subroutine allocate_mod_bdycon
    implicit none
    if ( idynamic == 3 ) then
      call getmem(fcx,1,nspgx,'bdycon:fcx')
      call getmem(zn1,jde1,jde2,ide1,ide2,'bdycon:zn1')
      call getmem(cnudge,1,kz,'bdycon:cnudge')
    else
      if ( iboudy == 1 .or. idynamic == 2 ) then
        call getmem(fcx,2,nspgx-1,'bdycon:fcx')
        call getmem(gcx,2,nspgx-1,'bdycon:gcx')
        call getmem(fcd,2,nspgd-1,'bdycon:fcd')
        call getmem(gcd,2,nspgd-1,'bdycon:gcd')
      end if
      if ( iboudy == 4 ) then
        call getmem(wgtd,2,nspgd-1,'bdycon:wgtd')
        call getmem(wgtx,2,nspgx-1,'bdycon:wgtx')
      end if
      if ( iboudy >= 5 ) then
        call getmem(hefc,2,nspgx-1,1,kz,'bdycon:hefc')
        call getmem(hegc,2,nspgx-1,1,kz,'bdycon:hegc')
        call getmem(hefd,2,nspgd-1,1,kz,'bdycon:hefd')
        call getmem(hegd,2,nspgd-1,1,kz,'bdycon:hegd')
      end if
      if ( ma%has_bdytop ) then
        call getmem(nue,jde1ga,jde2ga,1,kz,'bdycon:nue')
        call getmem(nui,jde1ga,jde2ga,1,kz,'bdycon:nui')
        call getmem(nve,jde1ga,jde2ga,1,kz,'bdycon:nve')
        call getmem(nvi,jde1ga,jde2ga,1,kz,'bdycon:nvi')
      end if
      if ( ma%has_bdybottom ) then
        call getmem(sue,jde1ga,jde2ga,1,kz,'bdycon:sue')
        call getmem(sui,jde1ga,jde2ga,1,kz,'bdycon:sui')
        call getmem(sve,jde1ga,jde2ga,1,kz,'bdycon:sve')
        call getmem(svi,jde1ga,jde2ga,1,kz,'bdycon:svi')
      end if
      if ( ma%has_bdyright ) then
        call getmem(eue,ide1ga,ide2ga,1,kz,'bdycon:eue')
        call getmem(eui,ide1ga,ide2ga,1,kz,'bdycon:eui')
        call getmem(eve,ide1ga,ide2ga,1,kz,'bdycon:eve')
        call getmem(evi,ide1ga,ide2ga,1,kz,'bdycon:evi')
      end if
      if ( ma%has_bdyleft ) then
        call getmem(wue,ide1ga,ide2ga,1,kz,'bdycon:wue')
        call getmem(wui,ide1ga,ide2ga,1,kz,'bdycon:wui')
        call getmem(wve,ide1ga,ide2ga,1,kz,'bdycon:wve')
        call getmem(wvi,ide1ga,ide2ga,1,kz,'bdycon:wvi')
      end if
      call getmem(psdot,jde1,jde2,ide1,ide2,'bdycon:psdot')
      call getmem(fg1,jde1ga,jde2ga,ide1ga,ide2ga,1,kzp1,'bdycon:fg1')
      call getmem(fg2,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'bdycon:fg2')
    end if
  end subroutine allocate_mod_bdycon

  subroutine setup_bdycon
    implicit none
    real(rkx), dimension(kz) :: anudge
    real(rkx) :: xfun, nb2, c1, c2
    integer(ik4) :: n, k
    ! Using Supposedly maximum Jet Stream velocity
    ! Note: The maximum at surface registered is 113 m/s
    !   (Tropical Cyclone Olivia on Barrow Island on April 10, 1996)
    real(rkx), parameter :: speedlimit = 130.0_rkx
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'setup_bdycon'
    integer(ik4), save :: idindx = 0
#endif
#ifdef DEBUG
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    ! Development of a Second-Generation Regional Climate Model (RegCM2).
    ! Part II: Convective Processes and Assimilation of Lateral Boundary
    ! Conditions
    !
    ! Article in Monthly Weather Review · October 1993
    ! DOI: 10.1175/1520-0493(1993)121<2814:DOASGR>2.0.CO;2
    !
    rtb = d_one / dtbdys
    if ( idynamic == 3 ) then
      fcx(1) = 1.0_rkx
      fcx(nspgx) = 0.0_rkx
      c1 = sqrt(d_two)*speedlimit*dtsec/dx/real(mo_nadv,rkx)
      c2 = sqrt(d_two)*sqrt(cpd/cvd*rgas*313.16_rkx)* &
              dtsec/dx/real(mo_nadv,rkx)/real(mo_nsound,rkx)
      call relax_coefficients(nspgx-2,0.0077_rkx,max(c1,c2),fcx(2:nspgx-1))
      if ( myid == 0 ) then
        call vprntv(fcx(2:nspgx-1),nspgx-2,'Boundary coefficients')
      end if
      if ( mo_spectral_nudging ) call lowpass_init( )
    else
      if ( iboudy == 1 .or. iboudy == 5 .or. iboudy == 6 ) then
        if ( bdy_nm > d_zero ) then
          fnudge = bdy_nm
        else
          fnudge = 0.1_rkx/dt2
        end if
        if ( bdy_dm > d_zero ) then
          gnudge = bdy_dm
        else
          ! The dxsq is simplified in below when dividing by dxsq
          gnudge = 0.02_rkx/dt2
        end if
        if ( myid == italk ) then
          write(stdout, '(a,f12.8,a,f12.8)') &
            ' Nudging coefficients F1=',fnudge,', F2=',gnudge
        end if
      end if
      if ( iboudy == 1 .or. idynamic == 2 ) then
        do n = 2, nspgx-1
          xfun = real(nspgx-n,rkx)/real(nspgx-2,rkx)
          fcx(n) = fnudge*xfun
          gcx(n) = gnudge*xfun
        end do
        do n = 2, nspgd-1
          xfun = real(nspgd-n,rkx)/real(nspgd-2,rkx)
          fcd(n) = fnudge*xfun
          gcd(n) = gnudge*xfun
        end do
      end if
      if ( iboudy == 4 ) then
        wgtd(2) = 0.20_rkx
        wgtd(3) = 0.55_rkx
        wgtd(4) = 0.80_rkx
        wgtd(5) = 0.95_rkx
        do k = 6, nspgd-1
          wgtd(k) = d_one
        end do
        wgtx(2) = 0.4_rkx
        wgtx(3) = 0.7_rkx
        wgtx(4) = 0.9_rkx
        do k = 5, nspgx-1
          wgtx(k) = 1.0_rkx
        end do
      end if
      if ( iboudy == 5 ) then
        call exponential_nudging(anudge)
        do k = 1, kz
          do n = 2, nspgx-1
            xfun = exp(-(real(n-2,rkx)/anudge(k)))
            hefc(n,k) = fnudge*xfun
            hegc(n,k) = gnudge*xfun
          end do
          do n = 2, nspgd-1
            xfun = exp(-(real(n-2,rkx)/anudge(k)))
            hefd(n,k) = fnudge*xfun
            hegd(n,k) = gnudge*xfun
          end do
        end do
      end if
      if ( iboudy == 6 ) then
        do k = 1, kz
          nb2 = d_two * nspgx
          do n = 2, nspgx-1
            xfun = d_half * (d_one - cos(mathpi/((n-nb2)/nb2)))
            hefc(n,k) = fnudge*xfun
            hegc(n,k) = gnudge*xfun
          end do
          nb2 = d_two * nspgd
          do n = 2, nspgd-1
            xfun = d_half * (d_one - cos(mathpi/((n-nb2)/nb2)))
            hefd(n,k) = fnudge*xfun
            hegd(n,k) = gnudge*xfun
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine setup_bdycon

  subroutine init_bdy
    implicit none
    integer(ik4) :: datefound, i, j, k, n
    character(len=32) :: appdat
    type (rcm_time_and_date) :: icbc_date
    type (rcm_time_interval) :: tdif
    integer(ik4), pointer, contiguous, dimension(:,:) :: dom_ldmsk
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: sub_ldmsk
    real(rkx), pointer, contiguous, dimension(:,:) :: dom_lndcat
    real(rkx), pointer, contiguous, dimension(:,:) :: ts0, ts1
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sfice
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sncv
    real(rkx), pointer, contiguous, dimension(:,:,:) :: snag
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'init_bdy'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    dom_ldmsk => mddom%ldmsk
    dom_lndcat => mddom%lndcat
    sub_ldmsk => mdsub%ldmsk
    ts0 => xtsb%b0
    ts1 => xtsb%b1
    sfice => lms%sfice
    sncv => lms%sncv
    snag => lms%snag

    bdydate1 = idate1
    bdydate2 = idate1
    xbctime = d_zero

    if ( bdydate1 == globidate1 ) then
      icbc_date = bdydate1
    else
      icbc_date = monfirst(bdydate1)
    end if

    call open_icbc(icbc_date)
    if ( islab_ocean == 1 .and. do_qflux_adj ) then
      call open_som
    end if
    call fixqcqi( )

    if ( we_have_qc( ) ) then
     call allocate_v3dbound(xlb,kz,cross)
     present_qc = .true.
    end if
    if ( we_have_qi( ) ) then
     call allocate_v3dbound(xib,kz,cross)
     present_qi = .true.
    end if

    datefound = icbc_search(bdydate1)
    if (datefound < 0) then
      !
      ! Cannot run without initial conditions
      !
      appdat = tochar(bdydate1)
      call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
    end if

    if ( idynamic == 2 ) then
      call read_icbc(nhbh0%ps,xtsb%b0,mddom%ldmsk,dub%b0,dvb%b0, &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0,&
                     xpaib%b0)

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh0%ps(j,i) = nhbh0%ps(j,i) * d_r10 - ptop
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh0%tvirt(j,i,k) = xtb%b0(j,i,k)*(d_one+ep1*xqb%b0(j,i,k))
        end do
      end if
    else if ( idynamic == 3 ) then
      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,dub%b0,dvb%b0,  &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0,&
                     xpaib%b0)
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b0,xqb%b0,dub%b0,dvb%b0,xpsb%b0,xtsb%b0)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b0,xqb%b0,dub%b0,dvb%b0,xpsb%b0,xtsb%b0)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh0%ps(j,i) = xpsb%b0(j,i)
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh0%tvirt(j,i,k) = xtb%b0(j,i,k)*(d_one+ep1*xqb%b0(j,i,k))
        end do
      end if
    else
      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,dub%b0,dvb%b0,  &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0,&
                     xpaib%b0)
    end if

    if ( islab_ocean == 1 .and. do_qflux_adj ) then
      som_month = rcmtimer%month
      datefound = som_search(som_month)
      if (datefound < 0) then
        appdat = tochar(bdydate1)
        call fatal(__FILE__,__LINE__,'SOM for '//appdat//' not found')
      end if
      call read_som(qflb0)
      where ( mddom%ldmsk == 1 ) qflb0 = d_zero
      tdif = bdydate1-monfirst(bdydate1)
      xslabtime = tohours(tdif)*secph
    end if

    if ( myid == italk ) then
      appdat = tochar(bdydate1)
      if ( rcmtimer%start( ) ) then
        write(stdout,*) 'READY IC DATA for ', appdat
      else
        write(stdout,*) 'READY BC DATA for ', appdat
      end if
    end if

    if ( idynamic == 1 ) then
      !$acc kernels
      xpsb%b0(:,:) = (xpsb%b0(:,:)*d_r10)-ptop
      !$acc end kernels
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b0,psdot)
    else if ( idynamic == 2 ) then
      !$acc kernels
      xpsb%b0(:,:) = atm0%ps(:,:) * d_r1000 ! Cb
      psdot(:,:) = atm0%psdot(jde1:jde2,ide1:ide2) * d_r1000
      xpsb%b1(:,:) = xpsb%b0(:,:)
      !$acc end kernels
    else
      !$acc kernels
      xpsb%b0(:,:) = xpsb%b0(:,:)*d_100
      !$acc end kernels
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Calculate P* on dot points
    ! Couple pressure u,v,t,q (pp,ww)
    !
    if ( idynamic /= 3 ) then
      call couple(dub%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(dvb%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      end if
      call exchange(dub%b0,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b0,1,jde1,jde2,ide1,ide2,1,kz)
    else
      call exchange(dub%b0,2,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b0,2,jde1,jde2,ide1,ide2,1,kz)
    end if
    call exchange(xtb%b0,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b0,1,jce1,jce2,ice1,ice2,1,kz)
    if ( present_qc ) then
      call exchange(xlb%b0,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( present_qi ) then
      call exchange(xib%b0,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( idynamic == 2 ) then
      call couple(xppb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      call couple(xwwb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%b0,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%b0,1,jce1,jce2,ice1,ice2,1,kzp1)
    else if ( idynamic == 3 ) then
      call exchange(xpaib%b0,1,jce1,jce2,ice1,ice2,1,kz)
    end if

    bdydate2 = bdydate2 + intbdy
    if ( myid == italk ) then
      write(stdout,'(a,a,a,i8)') ' SEARCH BC data for ', tochar10(bdydate2), &
                      ', step = ', rcmtimer%lcount
    end if
    datefound = icbc_search(bdydate2)
    if ( datefound < 0 ) then
      call open_icbc(monfirst(bdydate2))
      datefound = icbc_search(bdydate2)
      if ( datefound < 0 ) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
      end if
    end if

    if ( idynamic == 2 ) then
      call read_icbc(nhbh1%ps,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                     xpaib%b1)
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh1%ps(j,i) = nhbh1%ps(j,i) * d_r10 - ptop
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
        end do
      end if
    else if ( idynamic == 3 ) then
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1,  &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                     xpaib%b1)
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b1,xqb%b1,dub%b1,dvb%b1,xpsb%b1,xtsb%b1)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b1,xqb%b1,dub%b1,dvb%b1,xpsb%b1,xtsb%b1)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh1%ps(j,i) = xpsb%b1(j,i)
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
        end do
      end if
    else
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1,  &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                     xpaib%b1)
    end if

    if ( islab_ocean == 1 .and. do_qflux_adj ) then
      datefound = som_search(som_month+1)
      if (datefound < 0) then
        !
        ! Cannot run without initial conditions
        !
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'SOM for '//appdat//' not found')
      end if
      call read_som(qflb1)
      where ( mddom%ldmsk == 1 ) qflb1 = d_zero
      tdif = bdydate2-prevmon(bdydate2)
      qflbt = (qflb1-qflb0)/(real(tohours(tdif),rkx)*secph)
    end if

    if ( myid == italk ) then
      write (stdout,*) 'READY  BC from     ', &
            tochar10(bdydate1), ' to ', tochar10(bdydate2)
    end if

    bdydate1 = bdydate2
    !
    ! Repeat for T2
    !
    if ( idynamic == 1 ) then
      !$acc kernels
      xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
      !$acc end kernels
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    else if ( idynamic == 3 ) then
      !$acc kernels
      xpsb%b1(:,:) = xpsb%b1(:,:)*d_100
      !$acc end kernels
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(dub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(dvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      call exchange(dub%b1,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
    else
      call exchange(dub%b1,2,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b1,2,jde1,jde2,ide1,ide2,1,kz)
    end if
    call exchange(xtb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    if ( present_qc ) then
      call exchange(xlb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( present_qi ) then
      call exchange(xib%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( idynamic == 2 ) then
      call couple(xppb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xwwb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%b1,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%b1,1,jce1,jce2,ice1,ice2,1,kzp1)
    else if ( idynamic == 3 ) then
      call exchange(xpaib%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if

    if ( rcmtimer%start( ) ) then
      if ( iseaice == 1 ) then
        if ( islab_ocean == 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            if ( dom_ldmsk(j,i) == 1 ) cycle
            if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
            if ( iocncpl == 1 .or. iwavcpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            if ( ts0(j,i) <= icetriggert ) then
              ts0(j,i) = icetriggert
              dom_ldmsk(j,i) = 2
            end if
          end do
          do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
            if ( dom_ldmsk(j,i) == 1 ) cycle
            if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
            if ( iocncpl == 1 .or. iwavcpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            if ( ts0(j,i) <= icetriggert ) then
              if ( sub_ldmsk(n,j,i) == 0 ) then
                sub_ldmsk(n,j,i) = 2
                sfice(n,j,i) = 1.00_rkx
                sncv(n,j,i) = 1.0_rkx   ! 1 mm of snow over the ice
                snag(n,j,i) = 0.1_rkx
              end if
            end if
          end do
        end if
      end if
    end if

    if ( iseaice == 1 ) then
      if ( islab_ocean == 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          ! Update temperatures over ocean water and lakes
          if ( dom_ldmsk(j,i) == 1 ) cycle
          if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
          if ( iocncpl == 1 .or. iwavcpl == 1 ) then
            if ( cplmsk(j,i) /= 0 ) cycle
          end if
          if ( ts1(j,i) <= icetriggert ) then
            ts1(j,i) = icetriggert
            dom_ldmsk(j,i) = 2
          end if
        end do
        do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
          ! Update temperatures over ocean water and lakes
          if ( dom_ldmsk(j,i) == 1 ) cycle
          if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
          if ( iocncpl == 1 .or. iwavcpl == 1 ) then
            if ( cplmsk(j,i) /= 0 ) cycle
          end if
          if ( ts1(j,i) <= icetriggert ) then
            if ( sub_ldmsk(n,j,i) == 0 ) then
              sub_ldmsk(n,j,i) = 2
              sfice(n,j,i) = 1.00_rkx
            end if
          else
            if ( dom_ldmsk(j,i) == 2 ) then
              ! Decrease the surface ice to melt it
              if ( sub_ldmsk(n,j,i) == 2 ) then
                sfice(n,j,i) = sfice(n,j,i)*d_r10
              end if
            end if
          end if
        end do
      end if
    end if
#ifdef ASYNC_NETCDF
    call init_bdyin_prefetch()
#endif
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine init_bdy
  !
#ifdef ASYNC_NETCDF
  subroutine init_bdyin_prefetch()
    implicit none

    if ( .not. bdyin_prefetch_steps_initialized ) then
      call init_bdyin_prefetch_steps()
    end if
    if ( bdyin_prefetch_steps <= 0 ) return
    call warmup_icbc_prefetch()
  end subroutine init_bdyin_prefetch
  !
  subroutine bdyin_prefetch
    implicit none
    type(rcm_time_and_date) :: target_date
    real(rkx) :: lead_seconds
    logical :: should_prefetch

    if ( .not. bdyin_prefetch_steps_initialized ) then
      call init_bdyin_prefetch_steps()
    end if
    if ( bdyin_prefetch_steps <= 0 ) return

    lead_seconds = real(bdyin_prefetch_steps,rkx)*dtsec
    should_prefetch = alarm_in_bdy%will_act(lead_seconds)
    if ( .not. should_prefetch ) return

    target_date = bdydate2 + intbdy
    if ( target_date > idate2 ) return
    call prefetch_icbc(target_date)
  end subroutine bdyin_prefetch
  !
  subroutine init_bdyin_prefetch_steps()
    implicit none
    character(len=32) :: env_value
    integer(ik4) :: env_len, env_stat, read_stat

    bdyin_prefetch_steps = 1
    call get_environment_variable('RCM_BDYIN_PREFETCH_STEPS',env_value, &
      length=env_len,status=env_stat)
    if ( env_stat == 0 .and. env_len > 0 ) then
      read(env_value(1:env_len),*,iostat=read_stat) bdyin_prefetch_steps
      if ( read_stat /= 0 ) bdyin_prefetch_steps = 1
    end if
    bdyin_prefetch_steps_initialized = .true.
  end subroutine init_bdyin_prefetch_steps
#endif
  !
  ! this subroutine reads in the boundary conditions.
  !
  subroutine bdyin
    !@acc use nvtx
    implicit none
    integer(ik4) :: i, j, k, n, datefound
    character(len=32) :: appdat
    logical :: update_slabocn
#ifdef ASYNC_NETCDF
    logical :: prefetched
#endif
    type (rcm_time_interval) :: tdif
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: q0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: q1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: l0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: l1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: i0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: i1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ww0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ww1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pai0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pai1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ts0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ts1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ps0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ps1 => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: dom_ldmsk => null( )
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: sub_ldmsk => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dom_lndcat => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sfice => null( )
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyin'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !@acc call nvtxStartRange("bdyin")
    u0 => dub%b0
    u1 => dub%b1
    v0 => dvb%b0
    v1 => dvb%b1
    t0 => xtb%b0
    t1 => xtb%b1
    q0 => xqb%b0
    q1 => xqb%b1
    if ( present_qc ) then
      l0 => xlb%b0
      l1 => xlb%b1
    end if
    if ( present_qi ) then
      i0 => xib%b0
      i1 => xib%b1
    end if
    if ( idynamic == 2 ) then
      pp0 => xppb%b0
      pp1 => xppb%b1
      ww0 => xwwb%b0
      ww1 => xwwb%b1
    else if ( idynamic == 3 ) then
      pai0 => xpaib%b0
      pai1 => xpaib%b1
    end if
    ts0 => xtsb%b0
    ts1 => xtsb%b1
    ps0 => xpsb%b0
    ps1 => xpsb%b1
    dom_ldmsk => mddom%ldmsk
    dom_lndcat => mddom%lndcat
    sub_ldmsk => mdsub%ldmsk
    sfice => lms%sfice

    update_slabocn = ( islab_ocean == 1 .and. &
      do_qflux_adj .and. som_month /= rcmtimer%month )

    xbctime = d_zero
    !$acc kernels
    u0(:,:,:) = u1(:,:,:)
    v0(:,:,:) = v1(:,:,:)
    t0(:,:,:) = t1(:,:,:)
    q0(:,:,:) = q1(:,:,:)
    !$acc end kernels
    if ( present_qc ) then
      !$acc kernels
      l0(:,:,:) = l1(:,:,:)
      !$acc end kernels
    end if
    if ( present_qi ) then
      !$acc kernels
      i0(:,:,:) = i1(:,:,:)
      !$acc end kernels
    end if
    !$acc kernels
    ts0(:,:) = ts1(:,:)
    !$acc end kernels
    if ( idynamic == 2 ) then
      !$acc kernels
      pp0(:,:,:) = pp1(:,:,:)
      ww0(:,:,:) = ww1(:,:,:)
      !$acc end kernels
    else if ( idynamic == 3 ) then
      !$acc kernels
      pai0(:,:,:) = pai1(:,:,:)
      !$acc end kernels
    else
      !$acc kernels
      ps0(:,:) = ps1(:,:)
      !$acc end kernels
    end if

    if ( update_slabocn ) then
      ! Data are monthly
      som_month = rcmtimer%month
      qflb0 = qflb1
      tdif = bdydate1-monfirst(bdydate1)
      xslabtime = tohours(tdif)*secph
    end if

    bdydate2 = bdydate2 + intbdy
    if ( myid == italk ) then
      write(stdout,'(a,a,a,i8)') ' SEARCH BC data for ', tochar10(bdydate2), &
                      ', step = ', rcmtimer%lcount
    end if
    datefound = icbc_search(bdydate2)
    if ( datefound < 0 ) then
      call open_icbc(monfirst(bdydate2))
      datefound = icbc_search(bdydate2)
      if ( datefound < 0 ) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
      end if
    end if
    if ( idynamic == 2 ) then
#ifdef ASYNC_NETCDF
      call consume_icbc_prefetch(bdydate2,nhbh1%ps,xtsb%b1, &
                     mddom%ldmsk,dub%b1,dvb%b1,xtb%b1,xqb%b1, &
                     xlb%b1,xib%b1,xppb%b1,xwwb%b1,xpaib%b1,prefetched)
      if ( .not. prefetched ) then
#endif
        call read_icbc(nhbh1%ps,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1, &
                       xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                       xpaib%b1)
#ifdef ASYNC_NETCDF
      end if
#endif
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh1%ps(j,i) = nhbh1%ps(j,i) * d_r10 - ptop
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh1%tvirt(j,i,k) = t1(j,i,k)*(d_one+ep1*q1(j,i,k))
        end do
      end if
    else if ( idynamic == 3 ) then
#ifdef ASYNC_NETCDF
      call consume_icbc_prefetch(bdydate2,xpsb%b1,xtsb%b1, &
                     mddom%ldmsk,dub%b1,dvb%b1,xtb%b1,xqb%b1, &
                     xlb%b1,xib%b1,xppb%b1,xwwb%b1,xpaib%b1,prefetched)
      if ( .not. prefetched ) then
#endif
        call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1,  &
                       xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                       xpaib%b1)
#ifdef ASYNC_NETCDF
      end if
#endif
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b1,xqb%b1,dub%b1,dvb%b1,xpsb%b1,xtsb%b1)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b1,xqb%b1,dub%b1,dvb%b1,xpsb%b1,xtsb%b1)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          nhbh1%ps(j,i) = ps1(j,i)
        end do
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          nhbh1%tvirt(j,i,k) = t1(j,i,k)*(d_one+ep1*q1(j,i,k))
        end do
      end if
    else
#ifdef ASYNC_NETCDF
      call consume_icbc_prefetch(bdydate2,xpsb%b1,xtsb%b1, &
                     mddom%ldmsk,dub%b1,dvb%b1,xtb%b1,xqb%b1, &
                     xlb%b1,xib%b1,xppb%b1,xwwb%b1,xpaib%b1,prefetched)
      if ( .not. prefetched ) then
#endif
        call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,dub%b1,dvb%b1,  &
                       xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1,&
                       xpaib%b1)
#ifdef ASYNC_NETCDF
      end if
#endif
    end if

    if ( update_slabocn ) then
      datefound = som_search(som_month)
      if ( datefound < 0 ) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'SOM for '//appdat//' not found')
      end if
      call read_som(qflb1)
      where ( mddom%ldmsk == 1 ) qflb1 = d_zero
      tdif = bdydate2-prevmon(bdydate2)
      qflbt = (qflb1-qflb0)/(real(tohours(tdif),rkx)*secph)
    end if
    !
    ! Convert surface pressure to pstar
    !
    if ( idynamic == 1 ) then
      !$acc kernels
      ps1(:,:) = (ps1(:,:)*d_r10)-ptop
      !$acc end kernels
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    else if ( idynamic == 3 ) then
      !$acc kernels
      ps1(:,:) = ps1(:,:)*d_100
      !$acc end kernels
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(dub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(dvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      call exchange(dub%b1,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
    else
      call exchange(dub%b1,2,jde1,jde2,ide1,ide2,1,kz)
      call exchange(dvb%b1,2,jde1,jde2,ide1,ide2,1,kz)
    end if
    call exchange(xtb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    if ( present_qc ) then
      call exchange(xlb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( present_qi ) then
      call exchange(xib%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( idynamic == 2 ) then
      call couple(xppb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xwwb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%b1,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%b1,1,jce1,jce2,ice1,ice2,1,kzp1)
    else if ( idynamic == 3 ) then
      call exchange(xpaib%b1,1,jce1,jce2,ice1,ice2,1,kz)
    end if
    !
    ! Update ground temperature on Ocean/Lakes
    !
    if ( iseaice == 1 ) then
      if ( islab_ocean == 0 ) then
#ifdef OPENACC
        !$acc parallel loop collapse(2) gang vector
        do i = ici1, ici2
        do j = jci1, jci2
#else
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
#endif
          ! Update temperatures over ocean water only
          if ( dom_ldmsk(j,i) == 1 ) cycle
          ! Skip lake points if lake model active
          if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
          ! Do not update if coupling and ocean active here
          if ( iocncpl == 1 .or. iwavcpl == 1 ) then
            if ( cplmsk(j,i) /= 0 ) cycle
          end if
          ! Sea ice correction
          if ( ts1(j,i) <= icetriggert ) then
            ts1(j,i) = icetriggert
            dom_ldmsk(j,i) = 2
          end if
#ifdef OPENACC
        end do
#endif
        end do
#ifdef OPENACC
        !$acc parallel loop collapse(3) gang vector
        do i = ici1, ici2
        do j = jci1, jci2
        do n = 1, nnsg
#else
        do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
#endif
          ! Update temperatures over ocean water only
          if ( dom_ldmsk(j,i) == 1 ) cycle
          ! Skip lake points if lake model active
          if ( lakemod == 1 .and. islake(dom_lndcat(j,i)) ) cycle
          ! Do not update if coupling and ocean active here
          if ( iocncpl == 1 .or. iwavcpl == 1 ) then
            if ( cplmsk(j,i) /= 0 ) cycle
          end if
          ! Sea ice correction
          if ( ts1(j,i) <= icetriggert ) then
            if ( sub_ldmsk(n,j,i) == 0 ) then
              sub_ldmsk(n,j,i) = 2
              sfice(n,j,i) = 1.00_rkx
            end if
          else
            if ( dom_ldmsk(j,i) == 2 ) then
              ! Decrease the surface ice to melt it
              if ( sub_ldmsk(n,j,i) == 2 ) then
                sfice(n,j,i) = sfice(n,j,i)*d_r10
              end if
            end if
          end if
#ifdef OPENACC
        end do
        end do
#endif
        end do
      end if
    end if

    if ( myid == italk ) then
      write (stdout,*) 'READY  BC from     ', &
            tochar10(bdydate1), ' to ', tochar10(bdydate2)
    end if

    bdydate1 = bdydate2

    if ( ichem == 1 ) then
      call chem_bdyin
    end if
    !@acc call nvtxEndRange
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine bdyin
  !
  ! This subroutine sets the boundary values of u and v according
  ! to the boundary conditions specified.
  !
  !     xt : elapsed time from the initial boundary values.
  !
  subroutine bdyuv(x0,x1)
    implicit none
    real(rkx), intent(in) :: x0, x1
    integer(ik4) :: i, j, k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyuv'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Now compute last two points values in U and V
    ! Internal points
    !
    if ( ma%has_bdyleft ) then
      do concurrent ( i = idi1:idi2, k = 1:kz )
        wui(i,k) = atm1%u(jdi1,i,k)
        wvi(i,k) = atm1%v(jdi1,i,k)
      end do
    end if
    if ( ma%has_bdyright ) then
      do concurrent ( i = idi1:idi2, k = 1:kz )
        eui(i,k) = atm1%u(jdi2,i,k)
        evi(i,k) = atm1%v(jdi2,i,k)
      end do
    end if
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jdi1:jdi2, k = 1:kz )
        sui(j,k) = atm1%u(j,idi1,k)
        svi(j,k) = atm1%v(j,idi1,k)
      end do
    end if
    if ( ma%has_bdytop ) then
      do concurrent ( j = jdi1:jdi2, k = 1:kz )
        nui(j,k) = atm1%u(j,idi2,k)
        nvi(j,k) = atm1%v(j,idi2,k)
      end do
    end if
    !
    ! boundary slices:
    !
    if ( iboudy == 0 ) then
      !
      ! fixed boundary conditions:
      !
      ! west and east boundaries:
      !
      if ( ma%has_bdyleft ) then
        do concurrent ( i = ide1:ide2, k = 1:kz )
          wue(i,k) = dub%b0(jde1,i,k)
          wve(i,k) = dvb%b0(jde1,i,k)
        end do
      end if
      if ( ma%has_bdyright ) then
        do concurrent ( i = ide1:ide2, k = 1:kz )
          eue(i,k) = dub%b0(jde2,i,k)
          eve(i,k) = dvb%b0(jde2,i,k)
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%has_bdybottom ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          sue(j,k) = dub%b0(j,ide1,k)
          sve(j,k) = dvb%b0(j,ide1,k)
        end do
      end if
      if ( ma%has_bdytop ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          nue(j,k) = dub%b0(j,ide2,k)
          nve(j,k) = dvb%b0(j,ide2,k)
        end do
      end if
    else ! NOT Fixed
      !
      !     time-dependent boundary conditions:
      !
      ! west (j = 1) and east (j = jx) boundaries:
      !
      if ( ma%has_bdyleft ) then
        do concurrent ( i = idi1:idi2, k = 1:kz )
          wue(i,k) = (x0*dub%b0(jde1,i,k) + x1*dub%b1(jde1,i,k))
          wve(i,k) = (x0*dvb%b0(jde1,i,k) + x1*dvb%b1(jde1,i,k))
        end do
      end if
      if ( ma%has_bdyright ) then
        do concurrent ( i = idi1:idi2, k = 1:kz )
          eue(i,k) = (x0*dub%b0(jde2,i,k) + x1*dub%b1(jde2,i,k))
          eve(i,k) = (x0*dvb%b0(jde2,i,k) + x1*dvb%b1(jde2,i,k))
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%has_bdybottom ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          sue(j,k) = (x0*dub%b0(j,ide1,k) + x1*dub%b1(j,ide1,k))
          sve(j,k) = (x0*dvb%b0(j,ide1,k) + x1*dvb%b1(j,ide1,k))
        end do
      end if
      if ( ma%has_bdytop ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          nue(j,k) = (x0*dub%b0(j,ide2,k) + x1*dub%b1(j,ide2,k))
          nve(j,k) = (x0*dvb%b0(j,ide2,k) + x1*dvb%b1(j,ide2,k))
        end do
      end if
    end if
    !
    ! fill up the interior silces:
    !
    if ( ma%has_bdytopleft ) then
      do concurrent ( k = 1:kz )
        wui(ide2,k) = nue(jdi1,k)
        wvi(ide2,k) = nve(jdi1,k)
        nui(jde1,k) = wue(idi2,k)
        nvi(jde1,k) = wve(idi2,k)
      end do
    end if
    if ( ma%has_bdybottomleft ) then
      do concurrent ( k = 1:kz )
        wui(ide1,k) = sue(jdi1,k)
        wvi(ide1,k) = sve(jdi1,k)
        sui(jde1,k) = wue(idi1,k)
        svi(jde1,k) = wve(idi1,k)
      end do
    end if
    if ( ma%has_bdytopright ) then
      do concurrent ( k = 1:kz )
        eui(ide2,k) = nue(jdi2,k)
        evi(ide2,k) = nve(jdi2,k)
        nui(jde2,k) = eue(idi2,k)
        nvi(jde2,k) = eve(idi2,k)
      end do
    end if
    if ( ma%has_bdybottomright ) then
      do concurrent ( k = 1:kz )
        eui(ide1,k) = sue(jdi2,k)
        evi(ide1,k) = sve(jdi2,k)
        sui(jde2,k) = eue(idi1,k)
        svi(jde2,k) = eve(idi1,k)
      end do
    end if

    if ( ma%has_bdytop ) then
      call exchange_bdy_lr(nue,1,kz)
      call exchange_bdy_lr(nui,1,kz)
      call exchange_bdy_lr(nve,1,kz)
      call exchange_bdy_lr(nvi,1,kz)
    end if

    if ( ma%has_bdybottom ) then
      call exchange_bdy_lr(sue,1,kz)
      call exchange_bdy_lr(sui,1,kz)
      call exchange_bdy_lr(sve,1,kz)
      call exchange_bdy_lr(svi,1,kz)
    end if

    if ( ma%has_bdyleft ) then
      call exchange_bdy_bt(wue,1,kz)
      call exchange_bdy_bt(wui,1,kz)
      call exchange_bdy_bt(wve,1,kz)
      call exchange_bdy_bt(wvi,1,kz)
    end if
    if ( ma%has_bdyright ) then
      call exchange_bdy_bt(eue,1,kz)
      call exchange_bdy_bt(eui,1,kz)
      call exchange_bdy_bt(eve,1,kz)
      call exchange_bdy_bt(evi,1,kz)
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine bdyuv
  !
  ! This subroutine sets the boundary values for p*, p*u, p*v,
  ! p*t, p*qv, p*qc, and p*qr.
  !
  !     ---the boundary values of p*u and p*v are extrapolated from
  !        the interior points.
  !
  !     ---the boundary values of p* and p*t are specified.
  !
  !     ---the boundary values of p*qv, p*qc, and p*qr depend on
  !        inflow/outflow conditions, if iboudy = 3 or 4.
  !
  !     xt     : is the time in seconds the variables xxa represent.
  !
  subroutine bdyval
    implicit none
    integer(ik4) :: i, j, k, n
    real(rkx) :: qext, qint, qxint, qrat, windavg, tkeint
    real(rkx) :: x0, x1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyval'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Set time
    !
    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1
    !
    ! Set boundary conditions
    !
    ! MOLOCH
    !
    if ( idynamic == 3 ) then
      !
      ! west boundary: corners excluded
      !
      if ( ma%has_bdyleft ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          mo_atm%u(jde1,i,k) = x0*dub%b0(jde1,i,k) + x1*dub%b1(jde1,i,k)
        end do
        do concurrent ( i = idi1:idi2, k = 1:kz )
          mo_atm%v(jce1,i,k) = x0*dvb%b0(jce1,i,k) + x1*dvb%b1(jce1,i,k)
        end do
        do concurrent ( i = ici1:ici2, k = 1:kz )
          mo_atm%t(jce1,i,k) = x0*xtb%b0(jce1,i,k)+x1*xtb%b1(jce1,i,k)
          mo_atm%pai(jce1,i,k) = x0*xpaib%b0(jce1,i,k)+x1*xpaib%b1(jce1,i,k)
          mo_atm%qx(jce1,i,k,iqv) = x0*xqb%b0(jce1,i,k)+x1*xqb%b1(jce1,i,k)
        end do
        if ( present_qc ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            mo_atm%qx(jce1,i,k,iqc) = x0*xlb%b0(jce1,i,k)+x1*xlb%b1(jce1,i,k)
          end do
        end if
        if ( present_qi .and. ipptls > 1 ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            mo_atm%qx(jce1,i,k,iqi) = x0*xib%b0(jce1,i,k)+x1*xib%b1(jce1,i,k)
          end do
        end if
        do concurrent ( i = ici1:ici2, k = 1:kz, n = iqfrst:nqx )
          if ( present_qc .and. n == iqc ) cycle
          if ( present_qi .and. n == iqi ) cycle
          qxint = max(mo_atm%qx(jci1,i,k,n),d_zero)
          windavg = (mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k))
          if ( windavg > d_zero ) then
            mo_atm%qx(jce1,i,k,n) = 0.5_rkx*(qxbval(n)+qxint)
          else
            mo_atm%qx(jce1,i,k,n) = qxint
          end if
        end do
        if ( ibltyp == 2 ) then
          do concurrent ( i = ici1:ici2 )
            mo_atm%tke(jce1,i,1) = tkemin ! West boundary
          end do
          do concurrent ( i = ici1:ici2, k = 2:kz )
            tkeint = mo_atm%tke(jci1,i,k+1)
            windavg = mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k) + &
                      mo_atm%u(jde1,i,k-1) + mo_atm%u(jdi1,i,k-1)
            if ( windavg > d_zero ) then
              mo_atm%tke(jce1,i,k+1) = tkemin
            else
              mo_atm%tke(jce1,i,k+1) = tkeint
            end if
          end do
        end if
      end if
      !
      ! east boundary: corners excluded
      !
      if ( ma%has_bdyright ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          mo_atm%u(jde2,i,k) = x0*dub%b0(jde2,i,k) + x1*dub%b1(jde2,i,k)
        end do
        do concurrent ( i = idi1:idi2, k = 1:kz )
          mo_atm%v(jce2,i,k) = x0*dvb%b0(jce2,i,k) + x1*dvb%b1(jce2,i,k)
        end do
        do concurrent ( i = ici1:ici2, k = 1:kz )
          mo_atm%t(jce2,i,k) = x0*xtb%b0(jce2,i,k)+x1*xtb%b1(jce2,i,k)
          mo_atm%pai(jce2,i,k) = x0*xpaib%b0(jce2,i,k)+x1*xpaib%b1(jce2,i,k)
          mo_atm%qx(jce2,i,k,iqv) = x0*xqb%b0(jce2,i,k)+x1*xqb%b1(jce2,i,k)
        end do
        if ( present_qc ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            mo_atm%qx(jce2,i,k,iqc) = x0*xlb%b0(jce2,i,k)+x1*xlb%b1(jce2,i,k)
          end do
        end if
        if ( present_qi .and. ipptls > 1 ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            mo_atm%qx(jce2,i,k,iqi) = x0*xib%b0(jce2,i,k)+x1*xib%b1(jce2,i,k)
          end do
        end if
        do concurrent ( i = ici1:ici2, k = 1:kz, n = iqfrst:nqx )
          if ( present_qc .and. n == iqc ) cycle
          if ( present_qi .and. n == iqi ) cycle
          qxint = max(mo_atm%qx(jci2,i,k,n),d_zero)
          windavg = (mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k))
          if ( windavg < d_zero ) then
            mo_atm%qx(jce2,i,k,n) = 0.5_rkx*(qxbval(n)+qxint)
          else
            mo_atm%qx(jce2,i,k,n) = qxint
          end if
        end do
        if ( ibltyp == 2 ) then
          do concurrent ( i = ici1:ici2 )
            mo_atm%tke(jce2,i,1) = tkemin ! East boundary
          end do
          do concurrent ( i = ici1:ici2, k = 2:kz )
            tkeint = mo_atm%tke(jci2,i,k+1)
            windavg = mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k) + &
                      mo_atm%u(jde2,i,k-1) + mo_atm%u(jdi2,i,k-1)
            if ( windavg < d_zero ) then
              mo_atm%tke(jce2,i,k+1) = tkemin
            else
              mo_atm%tke(jce2,i,k+1) = tkeint
            end if
          end do
        end if
      end if
      !
      ! south boundary: corners included
      !
      if ( ma%has_bdybottom ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          mo_atm%u(j,ice1,k) = x0*dub%b0(j,ice1,k) + x1*dub%b1(j,ice1,k)
        end do
        do concurrent ( j = jce1:jce2, k = 1:kz )
          mo_atm%v(j,ide1,k) = x0*dvb%b0(j,ide1,k) + x1*dvb%b1(j,ide1,k)
        end do
        do concurrent ( j = jce1:jce2, k = 1:kz )
          mo_atm%t(j,ice1,k) = x0*xtb%b0(j,ice1,k)+x1*xtb%b1(j,ice1,k)
          mo_atm%pai(j,ice1,k) = x0*xpaib%b0(j,ice1,k)+x1*xpaib%b1(j,ice1,k)
          mo_atm%qx(j,ice1,k,iqv) = x0*xqb%b0(j,ice1,k)+x1*xqb%b1(j,ice1,k)
        end do
        if ( present_qc ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            mo_atm%qx(j,ice1,k,iqc) = x0*xlb%b0(j,ice1,k)+x1*xlb%b1(j,ice1,k)
          end do
        end if
        if ( present_qi .and. ipptls > 1 ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            mo_atm%qx(j,ice1,k,iqi) = x0*xib%b0(j,ice1,k)+x1*xib%b1(j,ice1,k)
          end do
        end if
        do concurrent ( j = jce1:jce2, k = 1:kz, n = iqfrst:nqx )
          if ( present_qc .and. n == iqc ) cycle
          if ( present_qi .and. n == iqi ) cycle
          qxint = max(mo_atm%qx(j,ici1,k,n),d_zero)
          windavg = (mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k))
          if ( windavg > d_zero ) then
            mo_atm%qx(j,ice1,k,n) = 0.5_rkx*(qxbval(n)+qxint)
          else
            mo_atm%qx(j,ice1,k,n) = qxint
          end if
        end do
        if ( ibltyp == 2 ) then
          do concurrent ( j = jce1:jce2 )
            mo_atm%tke(j,ice1,1) = tkemin  ! South boundary
          end do
          do concurrent ( j = jce1:jce2, k = 2:kz )
            tkeint = mo_atm%tke(j,ici1,k+1)
            windavg = mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k) + &
                      mo_atm%v(j,ide1,k-1) + mo_atm%v(j,idi1,k-1)
            if ( windavg > d_zero ) then
              mo_atm%tke(j,ice1,k+1) = tkemin
            else
              mo_atm%tke(j,ice1,k+1) = tkeint
            end if
          end do
        end if
      end if
      !
      ! north boundary: corners included
      !
      if ( ma%has_bdytop ) then
        do concurrent ( j = jde1:jde2, k = 1:kz )
          mo_atm%u(j,ice2,k) = x0*dub%b0(j,ice2,k) + x1*dub%b1(j,ice2,k)
        end do
        do concurrent ( j = jce1:jce2, k = 1:kz )
          mo_atm%v(j,ide2,k) = x0*dvb%b0(j,ide2,k) + x1*dvb%b1(j,ide2,k)
        end do
        do concurrent ( j = jce1:jce2, k = 1:kz )
          mo_atm%t(j,ice2,k) = x0*xtb%b0(j,ice2,k)+x1*xtb%b1(j,ice2,k)
          mo_atm%pai(j,ice2,k) = x0*xpaib%b0(j,ice2,k)+x1*xpaib%b1(j,ice2,k)
          mo_atm%qx(j,ice2,k,iqv) = x0*xqb%b0(j,ice2,k)+x1*xqb%b1(j,ice2,k)
        end do
        if ( present_qc ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            mo_atm%qx(j,ice2,k,iqc) = x0*xlb%b0(j,ice2,k)+x1*xlb%b1(j,ice2,k)
          end do
        end if
        if ( present_qi .and. ipptls > 1 ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            mo_atm%qx(j,ice2,k,iqi) = x0*xib%b0(j,ice2,k)+x1*xib%b1(j,ice2,k)
          end do
        end if
        do concurrent ( j = jce1:jce2, k = 1:kz, n = iqfrst:nqx )
          if ( present_qc .and. n == iqc ) cycle
          if ( present_qi .and. n == iqi ) cycle
          qxint = max(mo_atm%qx(j,ici2,k,n),d_zero)
          windavg = (mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k))
          if ( windavg < d_zero ) then
            mo_atm%qx(j,ice2,k,n) = 0.5_rkx*(qxbval(n)+qxint)
          else
            mo_atm%qx(j,ice2,k,n) = qxint
          end if
        end do
        if ( ibltyp == 2 ) then
          do concurrent ( j = jce1:jce2 )
            mo_atm%tke(j,ice2,1) = tkemin  ! North boundary
          end do
          do concurrent ( j = jce1:jce2, k = 2:kz )
            tkeint = mo_atm%tke(j,ici2,k+1)
            windavg = mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k) + &
                      mo_atm%v(j,ide2,k-1) + mo_atm%v(j,idi2,k-1)
            if ( windavg < d_zero ) then
              mo_atm%tke(j,ice2,k+1) = tkemin
            else
              mo_atm%tke(j,ice2,k+1) = tkeint
            end if
          end do
        end if
      end if

      if ( ichem == 1 ) then
        call chem_bdyval(mo_atm%u,mo_atm%v)
      end if

    else ! NO MOLOCH
      !
      ! Fill up the boundary value for xxb variables from xxa variables:
      ! if this subroutine is called for the first time, this part
      ! shall be skipped.
      !
      if ( rcmtimer%integrating( ) ) then
        !
        ! West boundary
        !
        if ( ma%has_bdyleft ) then
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm2%u(jde1,i,k) = atm1%u(jde1,i,k)
            atm2%v(jde1,i,k) = atm1%v(jde1,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm2%t(jce1,i,k) = atm1%t(jce1,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz, n = 1:nqx )
            atm2%qx(jce1,i,k,n) = atm1%qx(jce1,i,k,n)
          end do
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm2%pp(jce1,i,k) = atm1%pp(jce1,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm2%w(jce1,i,k) = atm1%w(jce1,i,k)
            end do
          else
            do concurrent ( i = ici1:ici2 )
              sfs%psb(jce1,i) = sfs%psa(jce1,i)
            end do
          end if
          if ( ibltyp == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm2%tke(jce1,i,k) = atm1%tke(jce1,i,k)
            end do
          end if
        end if
        !
        ! East boundary
        !
        if ( ma%has_bdyright ) then
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm2%u(jde2,i,k) = atm1%u(jde2,i,k)
            atm2%v(jde2,i,k) = atm1%v(jde2,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm2%t(jce2,i,k) = atm1%t(jce2,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz, n = 1:nqx )
            atm2%qx(jce2,i,k,n) = atm1%qx(jce2,i,k,n)
          end do
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm2%pp(jce2,i,k) = atm1%pp(jce2,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm2%w(jce2,i,k) = atm1%w(jce2,i,k)
            end do
          else
            do concurrent ( i = ici1:ici2 )
              sfs%psb(jce2,i) = sfs%psa(jce2,i)
            end do
          end if
          if ( ibltyp == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm2%tke(jce2,i,k) = atm1%tke(jce2,i,k)
            end do
          end if
        end if
        !
        ! North and South boundaries
        !
        if ( ma%has_bdybottom ) then
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm2%u(j,ide1,k) = atm1%u(j,ide1,k)
            atm2%v(j,ide1,k) = atm1%v(j,ide1,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm2%t(j,ice1,k) = atm1%t(j,ice1,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz, n = 1:nqx )
            atm2%qx(j,ice1,k,n) = atm1%qx(j,ice1,k,n)
          end do
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm2%pp(j,ice1,k) = atm1%pp(j,ice1,k)
            end do
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm2%w(j,ice1,k) = atm1%w(j,ice1,k)
            end do
          else
            do concurrent ( j = jce1:jce2 )
              sfs%psb(j,ice1) = sfs%psa(j,ice1)
            end do
          end if
          if ( ibltyp == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm2%tke(j,ice1,k) = atm1%tke(j,ice1,k)
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm2%u(j,ide2,k) = atm1%u(j,ide2,k)
            atm2%v(j,ide2,k) = atm1%v(j,ide2,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm2%t(j,ice2,k) = atm1%t(j,ice2,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz, n = 1:nqx )
            atm2%qx(j,ice2,k,n) = atm1%qx(j,ice2,k,n)
          end do
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm2%pp(j,ice2,k) = atm1%pp(j,ice2,k)
            end do
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm2%w(j,ice2,k) = atm1%w(j,ice2,k)
            end do
          else
            do concurrent ( j = jce1:jce2 )
              sfs%psb(j,ice2) = sfs%psa(j,ice2)
            end do
          end if
          if ( ibltyp == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm2%tke(j,ice2,k) = atm1%tke(j,ice2,k)
            end do
          end if
        end if
      end if

      if ( iboudy == 0 ) then
        !
        ! fixed boundary conditions:
        !
        if ( ma%has_bdyleft ) then
          if ( idynamic == 1 ) then
            do concurrent ( i = ici1:ici2 )
              sfs%psa(jce1,i) = xpsb%b0(jce1,i)
            end do
          end if
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm1%u(jde1,i,k) = dub%b0(jde1,i,k)
            atm1%v(jde1,i,k) = dvb%b0(jde1,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm1%t(jce1,i,k) = xtb%b0(jce1,i,k)
            atm1%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k)
          end do
          if ( present_qc ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce1,i,k,iqv) = xlb%b0(jce1,i,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce1,i,k,iqv) = xlb%b0(jce1,i,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%pp(jce1,i,k) = xppb%b0(jce1,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm1%w(jce1,i,k) = xwwb%b0(jce1,i,k)
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          if ( idynamic == 1 ) then
            do concurrent ( i = ici1:ici2 )
              sfs%psa(jce2,i) = xpsb%b0(jce2,i)
            end do
          end if
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm1%u(jde2,i,k) = dub%b0(jde2,i,k)
            atm1%v(jde2,i,k) = dvb%b0(jde2,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm1%t(jce2,i,k) = xtb%b0(jce2,i,k)
            atm1%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k)
          end do
          if ( present_qc ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce2,i,k,iqc) = xlb%b0(jce2,i,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce2,i,k,iqi) = xib%b0(jce2,i,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%pp(jce2,i,k) = xppb%b0(jce2,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm1%w(jce2,i,k) = xwwb%b0(jce2,i,k)
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          if ( idynamic == 1 ) then
            do concurrent ( j = jce1:jce2 )
              sfs%psa(j,ice1) = xpsb%b0(j,ice1)
            end do
          end if
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm1%u(j,ide1,k) = dub%b0(j,ide1,k)
            atm1%v(j,ide1,k) = dvb%b0(j,ide1,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm1%t(j,ice1,k) = xtb%b0(j,ice1,k)
            atm1%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k)
          end do
          if ( present_qc ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice1,k,iqc) = xlb%b0(j,ice1,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice1,k,iqi) = xib%b0(j,ice1,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%pp(j,ice1,k) = xppb%b0(j,ice1,k)
            end do
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm1%w(j,ice1,k) = xwwb%b0(j,ice1,k)
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          if (idynamic == 1 ) then
            do concurrent ( j = jce1:jce2 )
              sfs%psa(j,ice2) = xpsb%b0(j,ice2)
            end do
          end if
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm1%u(j,ide2,k) = dub%b0(j,ide2,k)
            atm1%v(j,ide2,k) = dvb%b0(j,ide2,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm1%t(j,ice2,k) = xtb%b0(j,ice2,k)
            atm1%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k)
          end do
          if ( present_qc ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice2,k,iqc) = xlb%b0(j,ice2,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice2,k,iqi) = xib%b0(j,ice2,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%pp(j,ice2,k) = xppb%b0(j,ice2,k)
            end do
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm1%w(j,ice2,k) = xwwb%b0(j,ice2,k)
            end do
          end if
        end if
      else
        !
        ! time-dependent boundary conditions:
        ! Set boundary values for p*:
        ! Set boundary conditions for p*u and p*v:
        !
        if ( ma%has_bdyleft ) then
          if ( idynamic == 1 ) then
            do concurrent ( i = ici1:ici2 )
              sfs%psa(jce1,i) = x0*xpsb%b0(jce1,i) + x1*xpsb%b1(jce1,i)
            end do
          end if
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm1%u(jde1,i,k) = x0*dub%b0(jde1,i,k) + x1*dub%b1(jde1,i,k)
            atm1%v(jde1,i,k) = x0*dvb%b0(jde1,i,k) + x1*dvb%b1(jde1,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm1%t(jce1,i,k)      = x0*xtb%b0(jce1,i,k)+x1*xtb%b1(jce1,i,k)
            atm1%qx(jce1,i,k,iqv) = x0*xqb%b0(jce1,i,k)+x1*xqb%b1(jce1,i,k)
          end do
          if ( present_qc ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce1,i,k,iqc) = x0*xlb%b0(jce1,i,k)+x1*xlb%b1(jce1,i,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce1,i,k,iqi) = x0*xib%b0(jce1,i,k)+x1*xib%b1(jce1,i,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%pp(jce1,i,k) = x0*xppb%b0(jce1,i,k)+x1*xppb%b1(jce1,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 2:kzp1 )
              atm1%w(jce1,i,k) = x0*xwwb%b0(jce1,i,k)+x1*xwwb%b1(jce1,i,k)
            end do
            do concurrent ( i = ici1:ici2 )
              atm1%w(jce1,i,1) = atm1%w(jci1,i,1)
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          if ( idynamic == 1 ) then
            do concurrent ( i = ici1:ici2 )
              sfs%psa(jce2,i) = x0*xpsb%b0(jce2,i) + x1*xpsb%b1(jce2,i)
            end do
          end if
          do concurrent ( i = idi1:idi2, k = 1:kz )
            atm1%u(jde2,i,k) = x0*dub%b0(jde2,i,k) + x1*dub%b1(jde2,i,k)
            atm1%v(jde2,i,k) = x0*dvb%b0(jde2,i,k) + x1*dvb%b1(jde2,i,k)
          end do
          do concurrent ( i = ici1:ici2, k = 1:kz )
            atm1%t(jce2,i,k)      = x0*xtb%b0(jce2,i,k)+x1*xtb%b1(jce2,i,k)
            atm1%qx(jce2,i,k,iqv) = x0*xqb%b0(jce2,i,k)+x1*xqb%b1(jce2,i,k)
          end do
          if ( present_qc ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce2,i,k,iqc) = x0*xlb%b0(jce2,i,k)+x1*xlb%b1(jce2,i,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%qx(jce2,i,k,iqi) = x0*xib%b0(jce2,i,k)+x1*xib%b1(jce2,i,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( i = ici1:ici2, k = 1:kz )
              atm1%pp(jce2,i,k) = x0*xppb%b0(jce2,i,k)+x1*xppb%b1(jce2,i,k)
            end do
            do concurrent ( i = ici1:ici2, k = 2:kzp1 )
              atm1%w(jce2,i,k) = x0*xwwb%b0(jce2,i,k)+x1*xwwb%b1(jce2,i,k)
            end do
            do concurrent ( i = ici1:ici2 )
              atm1%w(jce2,i,1) = atm1%w(jci2,i,1)
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          if ( idynamic == 1 ) then
            do concurrent ( j = jce1:jce2 )
              sfs%psa(j,ice1) = x0*xpsb%b0(j,ice1) + x1*xpsb%b1(j,ice1)
            end do
          end if
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm1%u(j,ide1,k) = x0*dub%b0(j,ide1,k) + x1*dub%b1(j,ide1,k)
            atm1%v(j,ide1,k) = x0*dvb%b0(j,ide1,k) + x1*dvb%b1(j,ide1,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm1%t(j,ice1,k)      = x0*xtb%b0(j,ice1,k)+x1*xtb%b1(j,ice1,k)
            atm1%qx(j,ice1,k,iqv) = x0*xqb%b0(j,ice1,k)+x1*xqb%b1(j,ice1,k)
          end do
          if ( present_qc ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice1,k,iqc) = x0*xlb%b0(j,ice1,k)+x1*xlb%b1(j,ice1,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice1,k,iqi) = x0*xib%b0(j,ice1,k)+x1*xib%b1(j,ice1,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%pp(j,ice1,k) = x0*xppb%b0(j,ice1,k)+x1*xppb%b1(j,ice1,k)
            end do
            do concurrent ( j = jce1:jce2, k = 2:kzp1 )
              atm1%w(j,ice1,k) = x0*xwwb%b0(j,ice1,k)+x1*xwwb%b1(j,ice1,k)
            end do
            do concurrent ( j = jce1:jce2 )
              atm1%w(j,ice1,1) = atm1%w(j,ici1,1)
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          if ( idynamic == 1 ) then
            do concurrent ( j = jce1:jce2 )
              sfs%psa(j,ice2) = x0*xpsb%b0(j,ice2) + x1*xpsb%b1(j,ice2)
            end do
          end if
          do concurrent ( j = jde1:jde2, k = 1:kz )
            atm1%u(j,ide2,k) = x0*dub%b0(j,ide2,k) + x1*dub%b1(j,ide2,k)
            atm1%v(j,ide2,k) = x0*dvb%b0(j,ide2,k) + x1*dvb%b1(j,ide2,k)
          end do
          do concurrent ( j = jce1:jce2, k = 1:kz )
            atm1%t(j,ice2,k)      = x0*xtb%b0(j,ice2,k)+x1*xtb%b1(j,ice2,k)
            atm1%qx(j,ice2,k,iqv) = x0*xqb%b0(j,ice2,k)+x1*xqb%b1(j,ice2,k)
          end do
          if ( present_qc ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice2,k,iqc) = x0*xlb%b0(j,ice2,k)+x1*xlb%b1(j,ice2,k)
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%qx(j,ice2,k,iqi) = x0*xib%b0(j,ice2,k)+x1*xib%b1(j,ice2,k)
            end do
          end if
          if ( idynamic == 2 ) then
            do concurrent ( j = jce1:jce2, k = 1:kz )
              atm1%pp(j,ice2,k) = x0*xppb%b0(j,ice2,k)+x1*xppb%b1(j,ice2,k)
            end do
            do concurrent ( j = jce1:jce2, k = 2:kzp1 )
               atm1%w(j,ice2,k) = x0*xwwb%b0(j,ice2,k)+x1*xwwb%b1(j,ice2,k)
            end do
            do concurrent ( j = jce1:jce2 )
              atm1%w(j,ice2,1) = atm1%w(j,ici2,1)
            end do
          end if
        end if
      end if

      call bdyuv(x0,x1)

      if ( iboudy == 3 .or. iboudy == 4 ) then
        !
        ! determine QV boundary values depends on inflow/outflow:
        !
        ! west boundary:
        !
        if ( ma%has_bdyleft ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            qext = atm1%qx(jce1,i,k,iqv)/sfs%psa(jce1,i)
            qint = atm1%qx(jci1,i,k,iqv)/sfs%psa(jci1,i)
            windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
            if ( windavg > d_zero ) then
              atm1%qx(jce1,i,k,iqv) = qext*sfs%psa(jce1,i)
            else
              atm1%qx(jce1,i,k,iqv) = qint*sfs%psa(jce1,i)
            end if
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            qext = atm1%qx(jce2,i,k,iqv)/sfs%psa(jce2,i)
            qint = atm1%qx(jci2,i,k,iqv)/sfs%psa(jci2,i)
            windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
            if ( windavg < d_zero ) then
              atm1%qx(jce2,i,k,iqv) = qext*sfs%psa(jce2,i)
            else
              atm1%qx(jce2,i,k,iqv) = qint*sfs%psa(jce2,i)
            end if
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            qext = atm1%qx(j,ice1,k,iqv)/sfs%psa(j,ice1)
            qint = atm1%qx(j,ici1,k,iqv)/sfs%psa(j,ici1)
            windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
            if ( windavg > d_zero ) then
              atm1%qx(j,ice1,k,iqv) = qext*sfs%psa(j,ice1)
            else
              atm1%qx(j,ice1,k,iqv) = qint*sfs%psa(j,ice1)
            end if
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do concurrent ( j = jce1:jce2, k = 1:kz )
            qext = atm1%qx(j,ice1,k,iqv)/sfs%psa(j,ice1)
            qext = atm1%qx(j,ice2,k,iqv)/sfs%psa(j,ice2)
            qint = atm1%qx(j,ici2,k,iqv)/sfs%psa(j,ici2)
            windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
            if ( windavg < d_zero ) then
              atm1%qx(j,ice2,k,iqv) = qext*sfs%psa(j,ice2)
            else
              atm1%qx(j,ice2,k,iqv) = qint*sfs%psa(j,ice2)
            end if
          end do
        end if
      end if

      if ( bdyflow ) then
        if ( ma%has_bdyleft ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( i = ici1:ici2, k = 1:kz )
              qxint = atm1%qx(jci1,i,k,n)
              windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
              if ( windavg > d_zero ) then
                atm1%qx(jce1,i,k,n) = qxbval(n)*sfs%psa(jce1,i)
              else
                atm1%qx(jce1,i,k,n) = qxint
              end if
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( i = ici1:ici2, k = 1:kz )
              qxint = atm1%qx(jci2,i,k,n)
              windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
              if ( windavg < d_zero ) then
                atm1%qx(jce2,i,k,n) = qxbval(n)**sfs%psa(jce2,i)
              else
                atm1%qx(jce2,i,k,n) = qxint
              end if
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( j = jce1:jce2, k = 1:kz )
              qxint = atm1%qx(j,ici1,k,n)
              windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
              if ( windavg > d_zero ) then
                atm1%qx(j,ice1,k,n) = qxbval(n)*sfs%psa(j,ice1)
              else
                atm1%qx(j,ice1,k,n) = qxint
              end if
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( j = jce1:jce2, k = 1:kz )
              qxint = atm1%qx(j,ici2,k,n)
              windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
              if ( windavg < d_zero ) then
                atm1%qx(j,ice2,k,n) = qxbval(n)*sfs%psa(j,ice2)
              else
                atm1%qx(j,ice2,k,n) = qxint
              end if
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( i = ici1:ici2, k = 1:kz )
              qxint = atm1%qx(jci1,i,k,n)/sfs%psa(jci1,i)
              qrat  = atm1%qx(jce1,i,k,iqv)/atm1%qx(jci1,i,k,iqv)
              atm1%qx(jce1,i,k,n) = qxint*sfs%psa(jce1,i)*qrat
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( i = ici1:ici2, k = 1:kz )
              qxint = atm1%qx(jci2,i,k,n)/sfs%psa(jci2,i)
              qrat  = atm1%qx(jce2,i,k,iqv)/atm1%qx(jci2,i,k,iqv)
              atm1%qx(jce2,i,k,n) = qxint*sfs%psa(jce2,i)*qrat
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( j = jce1:jce2, k = 1:kz )
              qxint = atm1%qx(j,ici1,k,n)/sfs%psa(j,ici1)
              qrat  = atm1%qx(j,ice1,k,iqv)/atm1%qx(j,ici1,k,iqv)
              atm1%qx(j,ice1,k,n) = qxint*sfs%psa(j,ice1)*qrat
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst, nqx
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do concurrent ( j = jce1:jce2, k = 1:kz )
              qxint = atm1%qx(j,ici2,k,n)/sfs%psa(j,ici2)
              qrat  = atm1%qx(j,ice2,k,iqv)/atm1%qx(j,ici2,k,iqv)
              atm1%qx(j,ice2,k,n) = qxint*sfs%psa(j,ice2)*qrat
            end do
          end do
        end if
      end if
      if ( ibltyp == 2 ) then
        if ( rcmtimer%start( ) ) then
          if ( ma%has_bdyleft ) then
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm1%tke(jce1,i,k) = tkemin ! East boundary
              atm2%tke(jce1,i,k) = tkemin ! East boundary
            end do
          end if
          if ( ma%has_bdyright ) then
            do concurrent ( i = ici1:ici2, k = 1:kzp1 )
              atm1%tke(jce2,i,k) = tkemin ! West boundary
              atm2%tke(jce2,i,k) = tkemin ! West boundary
            end do
          end if
          if ( ma%has_bdytop ) then
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm1%tke(j,ice2,k) = tkemin  ! North boundary
              atm2%tke(j,ice2,k) = tkemin  ! North boundary
            end do
          end if
          if ( ma%has_bdybottom ) then
            do concurrent ( j = jce1:jce2, k = 1:kzp1 )
              atm1%tke(j,ice1,k) = tkemin  ! South boundary
              atm2%tke(j,ice1,k) = tkemin  ! South boundary
            end do
          end if
        else
          ! if the boundary values and tendencies are not available,
          ! determine boundary values depends on inflow/outflow:
          ! inflow  : set it equal to zero.
          ! outflow : get from interior point.
          !
          ! west boundary:
          !
          if ( bdyflow ) then
            if ( ma%has_bdyleft ) then
              do concurrent ( i = ici1:ici2 )
                atm1%tke(jce1,i,1) = tkemin ! West boundary
                atm2%tke(jce1,i,1) = tkemin ! West boundary
              end do
              do concurrent ( i = ici1:ici2, k = 2:kz )
                tkeint = atm1%tke(jci1,i,k+1)
                windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k) + &
                  wue(i,k-1) + wue(i+1,k-1) + wui(i,k-1) + wui(i+1,k-1)
                if ( windavg > d_zero ) then
                  atm1%tke(jce1,i,k+1) = tkemin
                else
                  atm1%tke(jce1,i,k+1) = tkeint
                end if
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              do concurrent ( i = ici1:ici2 )
                atm1%tke(jce2,i,1) = tkemin ! East boundary
                atm2%tke(jce2,i,1) = tkemin ! East boundary
              end do
              do concurrent ( i = ici1:ici2, k = 2:kz )
                tkeint = atm1%tke(jci2,i,k+1)
                windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k) + &
                  eue(i,k-1) + eue(i+1,k-1) + eui(i,k-1) + eui(i+1,k-1)
                if ( windavg < d_zero ) then
                  atm1%tke(jce2,i,k+1) = tkemin
                else
                  atm1%tke(jce2,i,k+1) = tkeint
                end if
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              do concurrent ( j = jce1:jce2 )
                atm1%tke(j,ice1,1) = tkemin  ! South boundary
                atm2%tke(j,ice1,1) = tkemin  ! South boundary
              end do
              do concurrent ( j = jci1:jci2, k = 2:kz )
                tkeint = atm1%tke(j,ici1,k+1)
                windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k) + &
                  sve(j,k-1) + sve(j+1,k-1) + svi(j,k-1) + svi(j+1,k-1)
                if ( windavg > d_zero ) then
                  atm1%tke(j,ice1,k+1) = tkemin
                else
                  atm1%tke(j,ice1,k+1) = tkeint
                end if
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              do concurrent ( j = jce1:jce2 )
                atm1%tke(j,ice2,1) = tkemin  ! North boundary
                atm2%tke(j,ice2,1) = tkemin  ! North boundary
              end do
              do concurrent ( j = jce1:jce2, k = 2:kz )
                tkeint = atm1%tke(j,ici2,k+1)
                windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k) + &
                  nve(j,k-1) + nve(j+1,k-1) + nvi(j,k-1) + nvi(j+1,k-1)
                if ( windavg < d_zero ) then
                  atm1%tke(j,ice2,k+1) = tkemin
                else
                  atm1%tke(j,ice2,k+1) = tkeint
                end if
              end do
            end if
          else
            if ( ma%has_bdyleft ) then
              do concurrent ( i = ici1:ici2, k = 1:kzp1 )
                atm1%tke(jce1,i,k) = atm1%tke(jci1,i,k)
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              do concurrent ( i = ici1:ici2, k = 1:kzp1 )
                atm1%tke(jce2,i,k) = atm1%tke(jci2,i,k)
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              do concurrent ( j = jce1:jce2, k = 1:kzp1 )
                atm1%tke(j,ice1,k) = atm1%tke(j,ici1,k)
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              do concurrent ( j = jce1:jce2, k = 1:kzp1 )
                atm1%tke(j,ice2,k) = atm1%tke(j,ici2,k)
              end do
            end if
          end if
        end if
      end if

      if ( ichem == 1 ) then
        call chem_bdyval(sfs%psa,wue,wui,eue,eui,nve,nvi,sve,svi)
      end if

    end if ! MOLOCH/NO MOLOCH

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      ! Update temperatures over ocean water only
      if ( mddom%ldmsk(j,i) /= 0 ) cycle
      ! Skip lake points if lake model active
      if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
      ! Do not update if coupling and ocean active here
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        if ( cplmsk(j,i) /= 0 ) cycle
      end if
      ! FAB do not update if slaboc / adjust or restore run
      if ( islab_ocean == 1 ) cycle
      sfs%tg(j,i) = x0*xtsb%b0(j,i) + x1*xtsb%b1(j,i)
    end do

    xbctime = xbctime + dtsec

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine bdyval
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! this subroutine applies sponge boundary condition to the        c
  ! tendency term - ften.                                           c
  !                                                                 c
  ! nk    : is the number of vertical level to be adjusted.         c
  !                                                                 c
  ! ba    : is the boundary index structure                         c
  !                                                                 c
  ! wg    : are the weightings.                                     c
  !                                                                 c
  ! bnd   : Boundary condition data structure                       c
  !         2D or 3D (managed by interface declaration)             c
  !                                                                 c
  ! ften  : is the tendency calculated from the model.              c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine sponge4d(bnd,ften,m)
    implicit none
    integer(ik4), intent(in) :: m
    type(v3dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:,:) :: ften

    integer(ik4) :: i, j, k, ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge4d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      fg1(j,i,k) = (bnd%b1(j,i,k)-bnd%b0(j,i,k))/dtbdys
    end do

    if ( ba_cr%ns /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        if ( .not. ba_cr%bsouth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        if ( .not. ba_cr%bnorth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
       if ( .not. ba_cr%bwest(j,i) ) cycle
       ib = ba_cr%ibnd(j,i)
       ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        if ( .not. ba_cr%beast(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge4d

  subroutine spongeuv(bndu,bndv,ftenu,ftenv)
    implicit none
    type(v3dbound), intent(in) :: bndu, bndv
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: ftenu
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: ftenv
    integer(ik4) :: i, j, k, ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'spongeuv'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_dt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
      fg1(j,i,k) = (bndu%b1(j,i,k)-bndu%b0(j,i,k))/dtbdys
      fg2(j,i,k) = (bndv%b1(j,i,k)-bndv%b0(j,i,k))/dtbdys
    end do

    if ( ba_dt%ns /= 0 ) then
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
        if ( .not. ba_dt%bsouth(j,i) ) cycle
        ib = ba_dt%ibnd(j,i)
        ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k)+(d_one-wgtd(ib))*fg1(j,i,k)
        ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k)+(d_one-wgtd(ib))*fg2(j,i,k)
      end do
    end if
    if ( ba_dt%nn /= 0 ) then
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
        if ( .not. ba_dt%bnorth(j,i) ) cycle
        ib = ba_dt%ibnd(j,i)
        ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k)+(d_one-wgtd(ib))*fg1(j,i,k)
        ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k)+(d_one-wgtd(ib))*fg2(j,i,k)
      end do
    end if
    if ( ba_dt%nw /= 0 ) then
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
        if ( .not. ba_dt%bwest(j,i) ) cycle
        ib = ba_dt%ibnd(j,i)
        ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k)+(d_one-wgtd(ib))*fg1(j,i,k)
        ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k)+(d_one-wgtd(ib))*fg2(j,i,k)
      end do
    end if
    if ( ba_dt%ne /= 0 ) then
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
        if ( .not. ba_dt%beast(j,i) ) cycle
        ib = ba_dt%ibnd(j,i)
        ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k)+(d_one-wgtd(ib))*fg1(j,i,k)
        ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k)+(d_one-wgtd(ib))*fg2(j,i,k)
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine spongeuv

  subroutine sponge3d(bnd,ften)
    implicit none
    type(v3dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: ften
    integer(ik4) :: i, j, k, ib
    integer(ik4) :: nk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge3d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    nk = size(ften,3)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
      fg1(j,i,k) = (bnd%b1(j,i,k)-bnd%b0(j,i,k))/dtbdys
    end do

    if ( ba_cr%ns /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
        if ( .not. ba_cr%bsouth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k) = wgtx(ib)*ften(j,i,k)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
        if ( .not. ba_cr%bnorth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k) = wgtx(ib)*ften(j,i,k)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
        if ( .not. ba_cr%bwest(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k) = wgtx(ib)*ften(j,i,k)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:nk )
        if ( .not. ba_cr%beast(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i,k) = wgtx(ib)*ften(j,i,k)+(d_one-wgtx(ib))*fg1(j,i,k)
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge3d

  subroutine sponge2d(bnd,ften)
    implicit none
    type(v2dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:) :: ften
    integer(ik4) :: i, j, ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge2d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      fg1(j,i,1) = (bnd%b0(j,i)-bnd%b1(j,i))/dtbdys
    end do

    if ( ba_cr%ns /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( .not. ba_cr%bsouth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i) = wgtx(ib)*ften(j,i)+(d_one-wgtx(ib))*fg1(j,i,1)
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( .not. ba_cr%bnorth(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i) = wgtx(ib)*ften(j,i)+(d_one-wgtx(ib))*fg1(j,i,1)
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( .not. ba_cr%bwest(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i) = wgtx(ib)*ften(j,i)+(d_one-wgtx(ib))*fg1(j,i,1)
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        if ( .not. ba_cr%beast(j,i) ) cycle
        ib = ba_cr%ibnd(j,i)
        ften(j,i) = wgtx(ib)*ften(j,i)+(d_one-wgtx(ib))*fg1(j,i,1)
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge2d
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! These subroutines apply relaxation boundary conditions to the   c
  ! tendency term - ften - of variable f                            c
  !                                                                 c
  ! ldot  : logical dot (u,v) / cross (t,q,p) flag                  c
  !                                                                 c
  ! ip    : is the number of slices affected by nudging.            c
  !                                                                 c
  ! xt    : is the time in seconds for variable f                   c
  !                                                                 c
  ! ften  : is the tendency calculated from the model.              c
  !                                                                 c
  ! nk    : is the number of vertical level to be adjusted.         c
  !                                                                 c
  ! ibdy  : type of boundary condition relaxation, 1=linear         c
  !         5 = exponential                                         c
  !                                                                 c
  ! bnd   : Boundary condition data structure                       c
  !         2D or 3D (managed by interface declaration)             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine nudge4d3d(ibdy,f,bnd,ften,n)
    implicit none
    integer(ik4), intent(in) :: ibdy, n
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:,:) :: f
    type(v3dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:,:) :: ften
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, ib
    real(rkx) :: xf, xg, fls0, fls1, fls2, fls3, fls4
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge4d3d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kz )
      fg1(j,i,k) = x0*bnd%b0(j,i,k) + x1*bnd%b1(j,i,k) - f(j,i,k,n)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,k)
          xg = hegc(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,k)
          xg = hegc(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,k)
          xg = hegc(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,k)
          xg = hegc(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge4d3d

  subroutine nudgeuv(ibdy,fu,fv,bndu,bndv,ftenu,ftenv)
    implicit none
    integer(ik4), intent(in) :: ibdy
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:) :: fu, fv
    type(v3dbound), intent(in) :: bndu, bndv
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: ftenu, ftenv
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, ib
    real(rkx) :: xf, xg, fls0, fls1, fls2, fls3, fls4
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudgeuv'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_dt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    do concurrent ( j = jde1ga:jde2ga, i = ide1ga:ide2ga, k = 1:kz )
      fg1(j,i,k) = ((x0*bndu%b0(j,i,k) + x1*bndu%b1(j,i,k)) - fu(j,i,k))
      fg2(j,i,k) = ((x0*bndv%b0(j,i,k) + x1*bndv%b1(j,i,k)) - fv(j,i,k))
    end do

    if ( ibdy == 1 ) then
      if ( ba_dt%ns /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bsouth(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = fcd(ib)
          xg = gcd(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%nn /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bnorth(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = fcd(ib)
          xg = gcd(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%nw /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bwest(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = fcd(ib)
          xg = gcd(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%ne /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%beast(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = fcd(ib)
          xg = gcd(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 -  &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 -  &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    else
      if ( ba_dt%ns /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bsouth(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = hefc(ib,k)
          xg = hegc(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%nn /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bnorth(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = hefd(ib,k)
          xg = hegd(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%nw /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%bwest(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = hefd(ib,k)
          xg = hegd(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 - &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_dt%ne /= 0 ) then
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
          if ( .not. ba_dt%beast(j,i) ) cycle
          ib = ba_dt%ibnd(j,i)
          xf = hefd(ib,k)
          xg = hegd(ib,k)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ftenu(j,i,k) = ftenu(j,i,k) + xf*fls0 -  &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          fls0 = fg2(j,i,k)
          fls1 = fg2(j-1,i,k)
          fls2 = fg2(j+1,i,k)
          fls3 = fg2(j,i-1,k)
          fls4 = fg2(j,i+1,k)
          ftenv(j,i,k) = ftenv(j,i,k) + xf*fls0 -  &
                         xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudgeuv

  subroutine nudge4d(ibdy,f,bnd,ften,n1,n2)
    implicit none
    integer(ik4), intent(in) :: ibdy, n1, n2
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:,:) :: f
    type(v3dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:,:) :: ften
    integer(ik4) :: n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge4d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if
    do n = n1, n2
      call nudge4d3d(ibdy,f,bnd,ften,n)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge4d

  subroutine nudge3d(ibdy,f,bnd,ften)
    implicit none
    integer(ik4), intent(in) :: ibdy
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:) :: f
    type(v3dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: ften
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, ns, nk, ib
    real(rkx) :: xf, xg, fls0, fls1, fls2, fls3, fls4
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge3d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    ns = lbound(f,3)
    nk = ubound(f,3)
    !if ( nk == kzp1 ) ns = 2
    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = ns:nk )
      fg1(j,i,k) = (x0*bnd%b0(j,i,k) + x1*bnd%b1(j,i,k)) - f(j,i,k)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 -  &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,min(k,kz))
          xg = hegc(ib,min(k,kz))
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,min(k,kz))
          xg = hegc(ib,min(k,kz))
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,min(k,kz))
          xg = hegc(ib,min(k,kz))
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = ns:nk )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,min(k,kz))
          xg = hegc(ib,min(k,kz))
          fls0 = fg1(j,i,k)
          fls1 = fg1(j-1,i,k)
          fls2 = fg1(j+1,i,k)
          fls3 = fg1(j,i-1,k)
          fls4 = fg1(j,i+1,k)
          ften(j,i,k) = ften(j,i,k) + xf*fls0 -  &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge3d

  subroutine nudge2d(ibdy,f,bnd,ften)
    implicit none
    integer(ik4), intent(in) :: ibdy
    real(rkx), pointer, contiguous, intent(in), dimension(:,:) :: f
    type(v2dbound), intent(in) :: bnd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:) :: ften
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, ib
    real(rkx) :: xf, xg, fls0, fls1, fls2, fls3, fls4
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge2d'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga )
      fg1(j,i,1) = (x0*bnd%b0(j,i) + x1*bnd%b1(j,i)) - f(j,i)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = fcx(ib)
          xg = gcx(ib)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 -  &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,kz)
          xg = hegc(ib,kz)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,kz)
          xg = hegc(ib,kz)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                        xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,kz)
          xg = hegc(ib,kz)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 - &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          xf = hefc(ib,kz)
          xg = hegc(ib,kz)
          fls0 = fg1(j,i,1)
          fls1 = fg1(j-1,i,1)
          fls2 = fg1(j+1,i,1)
          fls3 = fg1(j,i-1,1)
          fls4 = fg1(j,i+1,1)
          ften(j,i) = ften(j,i) + xf*fls0 -  &
                      xg*(fls1+fls2+fls3+fls4-d_four*fls0)
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge2d

  subroutine couple(a,c,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: a
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: c
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    integer(ik4) :: i, j, k
    do concurrent ( j = j1:j2, i = i1:i2, k = k1:k2 )
      a(j,i,k) = a(j,i,k) * c(j,i)
    end do
  end subroutine couple

  subroutine raydampuv(z,u,v,uten,vten,ubnd,vbnd)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: u, v
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: uten
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vten
    type(v3dbound), intent(in) :: ubnd, vbnd
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, maxk
    real(rkx) :: bval

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1
    maxk = min(kzp1,rayndamp)
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:maxk )
      bval = x0*ubnd%b0(j,i,k) + x1*ubnd%b1(j,i,k)
      uten(j,i,k) = uten(j,i,k) + &
           tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd) * (bval-u(j,i,k))
    end do
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:maxk )
      bval = x0*vbnd%b0(j,i,k) + x1*vbnd%b1(j,i,k)
      vten(j,i,k) = vten(j,i,k) + &
        tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd) * (bval-v(j,i,k))
    end do
  end subroutine raydampuv

  subroutine raydampuv_c(z,u,v,uten,vten,sval)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: u, v
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: uten
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vten
    real(rkx), intent(in) :: sval
    integer(ik4) :: i, j, k, maxk
    maxk = min(kzp1,rayndamp)
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:maxk )
      uten(j,i,k) = uten(j,i,k) + &
        tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd) * (sval-u(j,i,k))
    end do
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:maxk )
      vten(j,i,k) = vten(j,i,k) + &
        tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd) * (sval-v(j,i,k))
    end do
  end subroutine raydampuv_c

  subroutine raydamp3f(z,var,vten,sval)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: var
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vten
    real(rkx), intent(in) :: sval
    integer(ik4) :: i, j, k, maxk
    maxk = min(kzp1,rayndamp)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:maxk )
        vten(j,i,k) = vten(j,i,k) + &
          tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd) * (sval-var(j,i,k))
    end do
  end subroutine raydamp3f

  subroutine raydamp3(z,var,vten,bnd)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: var
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: vten
    type(v3dbound), intent(in) :: bnd
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, maxk
    real(rkx) :: bval
    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1
    maxk = min(kz,rayndamp)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:maxk )
      bval = x0*bnd%b0(j,i,k) + x1*bnd%b1(j,i,k)
      vten(j,i,k) = vten(j,i,k) + &
            tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd)*(bval-var(j,i,k))
    end do
  end subroutine raydamp3

  subroutine raydampqv(z,var,vten,bnd)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var
    real(rkx), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: vten
    type(v3dbound), intent(in) :: bnd
    integer(ik4) :: i, j, k, maxk
    real(rkx) :: x0, x1
    real(rkx) :: bval
    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1
    maxk = min(kz,rayndamp)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:maxk )
      bval = x0*bnd%b0(j,i,k) + x1*bnd%b1(j,i,k)
      vten(j,i,k,iqv) = vten(j,i,k,iqv) + &
              tau(z(j,i,k),z(j,i,1),rayalpha0,rayhd)*(bval-var(j,i,k,iqv))
    end do
  end subroutine raydampqv

  pure real(rkx) function tau(z,zmax,r0,rhd)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: z, zmax, r0, rhd
    if ( z > zmax-rhd ) then
      tau = r0 * (sin(halfpi*(d_one-(zmax-z)/rhd)))**2
    else
      tau = d_zero
    end if
  end function tau

  subroutine paicompute(ps,z,t,q,pai)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: ps
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: z, t, q
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pai
    integer(ik4) :: i, j, k
    real(rkx) :: tv1, tv2, lrt, tv, zz, zb, p, zdelta, paikp1
    ! Hydrostatic initialization of pai
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      zdelta = z(j,i,kz)*egrav
      tv1 = t(j,i,kz) * (d_one + ep1*q(j,i,kz))
      tv2 = t(j,i,kz-1) * (d_one + ep1*q(j,i,kz-1))
      lrt = (tv2-tv1)/(z(j,i,kz-1)-z(j,i,kz))
      if ( lrt > govcp ) then
        lrt = govcp
      else if ( lrt < -0.005_rkx ) then
        lrt = 0.5_rkx*lrt - 0.5_rkx*lrate
      end if
      tv = tv1 - 0.5_rkx*z(j,i,kz)*lrt
      zz = d_one/(rgas*tv)
      p = ps(j,i) * exp(-zdelta*zz)
      paikp1 = (p/p00)**rovcp
      pai(j,i,kz) = paikp1
      !$acc loop seq
      do k = kzm1, 1, -1
        tv1 = t(j,i,k) * (d_one + ep1*q(j,i,k))
        tv2 = t(j,i,k+1) * (d_one + ep1*q(j,i,k+1))
        zb = d_two*egrav*mo_dzita/(mo_atm%fmzf(j,i,k+1)*cpd) + tv1 - tv2
        zdelta = sqrt(zb**2 + d_four * tv2 * tv1)
        paikp1 = -paikp1 / (d_two * tv2) * (zb - zdelta)
        pai(j,i,k) = paikp1
      end do
    end do
    call exchange(pai,1,jce1,jce2,ice1,ice2,1,kz)
  end subroutine paicompute

  subroutine moloch_static_test1(xt,xq,xu,xv,xps,xts)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xps, xts
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: xt, xq
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) ::  xu, xv
    integer(ik4) :: i, j, k
    xts = stdt
    xps = stdpmb
    xu = d_zero
    xv = d_zero
    xq = 1.0e-8_rkx
    do k = 1, kz
      do i = ice1, ice2
        do j = jce1, jce2
          xt(j,i,k) = max(xts(j,i) - lrate * mo_atm%zeta(j,i,k), 210.0_rkx)
        end do
      end do
    end do
  end subroutine moloch_static_test1

  subroutine moloch_static_test2(xt,xq,xu,xv,xps,xts)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xps, xts
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: xt, xq
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: xu, xv
    integer(ik4) :: i, j, k
    real(rkx) :: zlr
    xu = 10.0_rkx
    xv = d_zero
    xq = d_zero
    xts = stdt
    do k = 1, kz
      do i = ice1, ice2
        do j = jce1, jce2
          zlr = -lrate
          xt(j,i,k) = max(xts(j,i) + zlr * mo_atm%zeta(j,i,k), 210.0_rkx)
        end do
      end do
    end do
    do i = ice1, ice2
      do j = jce1, jce2
        xps(j,i) = stdpmb * exp(-govr*mo_atm%zeta(j,i,kz)/xt(j,i,kz))
      end do
    end do
  end subroutine moloch_static_test2

  subroutine lowpass_init
    implicit none
    real(rkx), dimension(:), allocatable :: px, py
    real(rkx) :: dx, dy
    integer(ik4) :: i, j, k, l
    real(rkx), parameter :: cutoff_wavelength_lon_km = 1500.0_rkx
    real(rkx), parameter :: cutoff_wavelength_lat_km =  750.0_rkx
    real(rk8) :: meanz, gmeanz, np

    km = max(nint((njcross*ds)/cutoff_wavelength_lon_km),1)
    lm = max(nint((nicross*ds)/cutoff_wavelength_lat_km),1)

    if ( myid == 0 ) then
      write(stdout, '(a)') ' Spectral nudging active.'
      dx = min((njcross*ds)/2,cutoff_wavelength_lon_km)
      dy = min((nicross*ds),cutoff_wavelength_lat_km)
      write(stdout, '(a,f8.2,a)') &
        '  Wavelenght cutoff for longitudinal waves : ',dx,' km'
      write(stdout, '(a,f8.2,a)') &
        '  Wavelenght cutoff for latitudinal waves  : ',dy,' km'
    end if

    dx = mathpi/real(njcross-1,rkx)
    dy = mathpi/real(nicross-1,rkx)

    allocate(px(2*km), py(2*lm))
    call getmem(bvx,jde1,jde2,1,2*km,'lowpass::bvx')
    call getmem(bvy,ide1,ide2,1,2*lm,'lowpass::bvy')
    call getmem(sx,ide1,ide2,1,2*km,'lowpass::sx')
    call getmem(sy,jde1,jde2,1,2*lm,'lowpass::sy')
    call getmem(sxg,ide1,ide2,1,2*km,'lowpass::sxg')
    call getmem(syg,jde1,jde2,1,2*lm,'lowpass::syg')
    do k = 1, 2*km
      px(k) = exp(-(real(k,rkx)/real(km,rkx))**2)
    end do
    do l = 1, 2*lm
      py(l) = exp(-(real(l,rkx)/real(lm,rkx))**2)
    end do
    do k = 1, 2*km
      do j = jde1, jde2
        bvx(j,k) = sqrt(2.0_rkx/real(jxm1,rkx)*px(k))*sin(k*(j-2)*dx)
      end do
    end do
    do l = 1, 2*lm
      do i = ide1, ide2
        bvy(i,l) = sqrt(2.0_rkx/real(iym1,rkx)*py(l))*sin(l*(i-2)*dy)
      end do
    end do
    deallocate(px,py)
    cn0 = dtrad/dtbdys
    np = real(jx*iy,rk8)
    do k = 1, kz
      meanz = 0.0_rkx
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
         meanz = meanz + real(mo_atm%zeta(j,i,k),rk8)/np
      end do
      call sumall(meanz,gmeanz)
      cnudge(k) = cn0 * min((real(gmeanz,rkx)/mo_h)**2, 1.0_rkx)
    end do
  end subroutine lowpass_init

  subroutine lowpass_filter(j1,j2,i1,i2,jj1,jj2,ii1,ii2,f)
    implicit none
    integer(ik4), intent(in) :: j1, j2, i1, i2, jj1, jj2, ii1, ii2
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:) :: f
    integer(ik4) :: i, j, k, l

    do k = 1, 2*km
      do i = i1, i2
        sx(i,k) = 0.0_rkx
        do j = jj1, jj2
          sx(i,k) = sx(i,k) + f(j,i)*bvx(j,k)
        end do
      end do
    end do
    call row_reduce(sx,sxg,i1,i2)
    f(:,:) = 0.0_rkx
    do k = 1, 2*km
      do i = i1, i2
        do j = j1, j2
          f(j,i) = f(j,i) + sxg(i,k)*bvx(j,k)
        end do
      end do
    end do
    do l = 1, 2*lm
      do j = j1, j2
        sy(j,l) = 0.
        do i = ii1, ii2
          sy(j,l) = sy(j,l) + f(j,i)*bvy(i,l)
        end do
      end do
    end do
    call column_reduce(sy,syg,j1,j2)
    f(:,:) = 0.0_rkx
    do l = 1, 2*lm
      do i = i1, i2
        do j = j1, j2
          f(j,i) = f(j,i) + syg(j,l)*bvy(i,l)
        end do
      end do
    end do
  end subroutine lowpass_filter

  subroutine mospectral_nudge(j1,j2,i1,i2,jj1,jj2,ii1,ii2,f,bnd)
    implicit none
    integer(ik4), intent(in) :: j1, j2, i1, i2, jj1, jj2, ii1, ii2
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: f
    type(v3dbound), intent(in) :: bnd
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    do k = 1, kz
      do concurrent ( j = j1:j2, i = i1:i2 )
        zn1(j,i) = f(j,i,k) - (x0*bnd%b0(j,i,k)+x1*bnd%b1(j,i,k))
      end do
      call lowpass_filter(j1,j2,i1,i2,jj1,jj2,ii1,ii2,zn1)
      do concurrent ( j = jj1:jj2, i = ii1:ii2 )
        f(j,i,k) = f(j,i,k) - cnudge(k)*zn1(j,i)
      end do
    end do
  end subroutine mospectral_nudge

  subroutine invert_top_bottom(v)
    implicit none
    real(rkx), dimension(:), intent(inout) :: v
    real(rkx), dimension(size(v)) :: swap
    integer(ik4) :: nk, k, kk
    swap = v
    nk = size(v)
    do k = 1, nk
      kk = nk-k+1
      v(k) = swap(kk)
    end do
  end subroutine invert_top_bottom

  subroutine moupdate_norm(ud,vd,unx)
    implicit none
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:) :: ud, vd
    real(rkx), pointer, contiguous, intent(inout), dimension(:,:,:) :: unx
    call outward_velocity(jci1,jci2,ici1,ici2,1,kz, &
                          jcross2,icross2,nspgx,ud,vd,unx)
  end subroutine moupdate_norm

  subroutine morelax(j1,j2,i1,i2,ba,f,bnd)
    implicit none
    integer(ik4) :: j1, j2, i1, i2
    real(rkx), pointer, contiguous, intent(in), dimension(:,:,:) :: f
    type(bound_area), intent(in) :: ba
    type(v3dbound), intent(in) :: bnd
    real(rkx) :: x0, x1
    integer(ik4) :: i, j, k, ib
    real(rkx) :: xf, fext
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'morelax'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    x1 = (xbctime + dt)/dtbdys
    x0 = 1.0_rkx - x1

    if ( ba%ns /= 0 ) then
      do concurrent ( j = j1:j2, i = i1:i2, k = 1:kz )
        if ( .not. ba%bsouth(j,i) ) cycle
        ib = ba%ibnd(j,i)
        xf = fcx(ib)
        fext = (x0*bnd%b0(j,i,k)+x1*bnd%b1(j,i,k))
        f(j,i,k) = (1.0_rkx-xf) * f(j,i,k) + xf*fext
      end do
    end if
    if ( ba%nn /= 0 ) then
      do concurrent ( j = j1:j2, i = i1:i2, k = 1:kz )
        if ( .not. ba%bnorth(j,i) ) cycle
        ib = ba%ibnd(j,i)
        xf = fcx(ib)
        fext = (x0*bnd%b0(j,i,k)+x1*bnd%b1(j,i,k))
        f(j,i,k) = (1.0_rkx-xf) * f(j,i,k) + xf*fext
      end do
    end if
    if ( ba%nw /= 0 ) then
      do concurrent ( j = j1:j2, i = i1:i2, k = 1:kz )
        if ( .not. ba%bwest(j,i) ) cycle
        ib = ba%ibnd(j,i)
        xf = fcx(ib)
        fext = (x0*bnd%b0(j,i,k)+x1*bnd%b1(j,i,k))
        f(j,i,k) = (1.0_rkx-xf) * f(j,i,k) + xf*fext
      end do
    end if
    if ( ba%ne /= 0 ) then
      do concurrent ( j = j1:j2, i = i1:i2, k = 1:kz )
        if ( .not. ba%beast(j,i) ) cycle
        ib = ba%ibnd(j,i)
        xf = fcx(ib)
        fext = (x0*bnd%b0(j,i,k)+x1*bnd%b1(j,i,k))
        f(j,i,k) = (1.0_rkx-xf) * f(j,i,k) + xf*fext
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine morelax

end module mod_bdycod

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
