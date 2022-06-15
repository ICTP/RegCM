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
  use mod_memutil
  use mod_atm_interface
  use mod_pbl_interface , only : tkemin
  use mod_che_interface
  use mod_lm_interface
  use mod_mpmessage
  use mod_ncio
  use mod_service
  use mod_zita
  use mod_stdatm
  use mod_slabocean
  use mod_vertint , only : intz1
  use mod_humid , only : rh2mxr

  implicit none

  private

  public :: initideal
  public :: allocate_mod_bdycon , init_bdy , bdyin , bdyval
  public :: sponge , nudge , setup_bdycon , raydamp
  public :: is_present_qc , is_present_qi

  !
  ! West U External  = WUE
  ! West U Internal  = WUI
  ! East U External  = EUE
  ! .....
  ! North V Internal = NVI
  !
  public :: wue , wui , eue , eui
  public :: wve , wvi , eve , evi
  public :: sue , sui , nue , nui
  public :: sve , svi , nve , nvi
  ! fnudge : are the coefficients for the newtonian term.
  ! gnydge : are the coefficients for the diffusion term.
  public :: fnudge , gnudge
  !
  real(rkx) , pointer , dimension(:,:) :: sue , sui , nue , nui , &
                                         sve , svi , nve , nvi
  real(rkx) , pointer , dimension(:,:) :: wue , wui , eue , eui , &
                                         wve , wvi , eve , evi
  real(rkx) , pointer , dimension(:,:) :: psdot
  real(rkx) , pointer , dimension(:) :: fcx , gcx
  real(rkx) , pointer , dimension(:) :: fcd , gcd
  real(rkx) , pointer , dimension(:,:) :: hefc , hegc , hefd , hegd
  real(rkx) , pointer , dimension(:) :: wgtd
  real(rkx) , pointer , dimension(:) :: wgtx
  real(rkx) , pointer , dimension(:,:,:) :: fg1 , fg2
  real(rkx) :: fnudge , gnudge , rdtbdy
  real(rk8) :: jday
  integer(ik4) :: som_month

  interface timeint
    module procedure timeint2 , timeint3
  end interface timeint

  interface nudge
    module procedure nudge4d
    module procedure nudge4d3d
    module procedure nudge3d
    module procedure nudge2d
    module procedure nudgeuv
    module procedure monudge4d
    module procedure monudge4d3d
    module procedure monudge3d
    module procedure monudge2d
    module procedure monudgeuv
  end interface nudge

  interface sponge
    module procedure sponge4d
    module procedure sponge3d
    module procedure sponge2d
    module procedure spongeuv
    module procedure mosponge4d
    module procedure mosponge3d
    module procedure mosponge2d
    module procedure mospongeuv
  end interface sponge

  interface raydamp
    module procedure raydamp3
    module procedure raydamp3f
    module procedure raydampuv
    module procedure raydampqv
    module procedure raydampuv_c
    module procedure moraydamp
  end interface raydamp

  logical , parameter :: bdyflow = .true.
  logical :: present_qc = .false.
  logical :: present_qi = .false.

  contains

  subroutine initideal
    implicit none
    integer(ik4) :: iunit , ierr
    integer(ik4) :: i , k , nlev , nseed
    real(rkx) :: ps , ts
    real(rkx) , allocatable , dimension(:) :: g , p , t , u , v , r
    real(rkx) , allocatable , dimension(:) :: zi , pi , ti , ui , vi , qi
    real(rkx) , dimension(jce1ga:jce2ga,ice1ga:ice2ga) :: noise
    integer(ik4) , allocatable , dimension(:) :: seed
    integer(ik8) :: sclock
    real(rkx) , dimension(1) :: ht

    namelist /dimensions/ nlev
    namelist /surface/ ps , ts
    namelist /profile/ g , p , t , u , v , r

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
      ! vertical interpolation for xub%b0,xvb%b0,xtb%b0,xqb%b0
      pi = sigma*(ps-ptop*10.0_rkx) + (ptop*10.0_rkx)
      call intz1(ui,u,pi,p,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,pi,p,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,pi,p,ht,1,1,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,pi,p,ht,1,1,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1 , kz
        xub%b0(:,:,k) = ui(k)
        xvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
      end do
    else if ( idynamic == 2 ) then
      ! vertical interpolation for xub%b0,xvb%b0,xtb%b0,xqb%b0
      xpsb%b0 = ps * 0.1_rkx
      ht = 0.0
      zi = atm0%z(jci1,ici1,1:kz)
      call intz1(ui,u,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,zi,g,ht,1,1,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1 , kz
        xub%b0(:,:,k) = ui(k)
        xvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
      end do
      xppb%b0 = 0.0_rkx
      xwwb%b0 = 0.0_rkx
    else if ( idynamic == 3 ) then
      ! vertical interpolation for xub%b0,xvb%b0,xtb%b0,xqb%b0
      xpsb%b0 = ps * 100.0_rkx
      ht = 0.0
      zi = mo_atm%zeta(jci1,ici1,1:kz)
      call intz1(ui,u,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(vi,v,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.2_rkx,0.2_rkx)
      call intz1(ti,t,zi,g,ht,1,1,kz,nlev,0.6_rkx,0.85_rkx,0.5_rkx)
      call intz1(qi,r,zi,g,ht,1,1,kz,nlev,0.7_rkx,0.7_rkx,0.4_rkx)
      do k = 1 , kz
        xub%b0(:,:,k) = ui(k)
        xvb%b0(:,:,k) = vi(k)
        xtb%b0(:,:,k) = ti(k)
        xqb%b0(:,:,k) = qi(k)
      end do
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
      call exchange(xub%b0,1,jde1,jde2,ice1,ice2,1,kz)
      call exchange(xvb%b0,1,jce1,jce2,ide1,ide2,1,kz)
      call exchange(xtb%b0,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xqb%b0,1,jce1,jce2,ice1,ice2,1,kz)
      call paicompute(mddom%xlat,xpsb%b0,mo_atm%zeta,xtb%b0,xqb%b0,xpaib%b0)
      call exchange(xpaib%b0,1,jce1,jce2,ice1,ice2,1,kz)
    else
      call fatal(__FILE__,__LINE__, &
        'Should never get here....')
    end if

    deallocate(seed)
    deallocate(g,p,t,u,v,r)
    deallocate(zi,pi,ti,ui,vi,qi)
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
    if ( iboudy == 1 .or. idynamic == 2 ) then
      call getmem1d(fcx,2,nspgx-1,'bdycon:fcx')
      call getmem1d(gcx,2,nspgx-1,'bdycon:gcx')
      call getmem1d(fcd,2,nspgd-1,'bdycon:fcd')
      call getmem1d(gcd,2,nspgd-1,'bdycon:gcd')
    end if
    if ( iboudy == 4 ) then
      call getmem1d(wgtd,2,nspgd-1,'bdycon:wgtd')
      call getmem1d(wgtx,2,nspgx-1,'bdycon:wgtx')
    end if
    if ( iboudy >= 5 ) then
      call getmem2d(hefc,2,nspgx-1,1,kz,'bdycon:hefc')
      call getmem2d(hegc,2,nspgx-1,1,kz,'bdycon:hegc')
      call getmem2d(hefd,2,nspgd-1,1,kz,'bdycon:hefd')
      call getmem2d(hegd,2,nspgd-1,1,kz,'bdycon:hegd')
    end if
    if ( idynamic /= 3 ) then
      if ( ma%has_bdytop ) then
        call getmem2d(nue,jde1ga,jde2ga,1,kz,'bdycon:nue')
        call getmem2d(nui,jde1ga,jde2ga,1,kz,'bdycon:nui')
        call getmem2d(nve,jde1ga,jde2ga,1,kz,'bdycon:nve')
        call getmem2d(nvi,jde1ga,jde2ga,1,kz,'bdycon:nvi')
      end if
      if ( ma%has_bdybottom ) then
        call getmem2d(sue,jde1ga,jde2ga,1,kz,'bdycon:sue')
        call getmem2d(sui,jde1ga,jde2ga,1,kz,'bdycon:sui')
        call getmem2d(sve,jde1ga,jde2ga,1,kz,'bdycon:sve')
        call getmem2d(svi,jde1ga,jde2ga,1,kz,'bdycon:svi')
      end if
      if ( ma%has_bdyright ) then
        call getmem2d(eue,ide1ga,ide2ga,1,kz,'bdycon:eue')
        call getmem2d(eui,ide1ga,ide2ga,1,kz,'bdycon:eui')
        call getmem2d(eve,ide1ga,ide2ga,1,kz,'bdycon:eve')
        call getmem2d(evi,ide1ga,ide2ga,1,kz,'bdycon:evi')
      end if
      if ( ma%has_bdyleft ) then
        call getmem2d(wue,ide1ga,ide2ga,1,kz,'bdycon:wue')
        call getmem2d(wui,ide1ga,ide2ga,1,kz,'bdycon:wui')
        call getmem2d(wve,ide1ga,ide2ga,1,kz,'bdycon:wve')
        call getmem2d(wvi,ide1ga,ide2ga,1,kz,'bdycon:wvi')
      end if
      call getmem2d(psdot,jde1,jde2,ide1,ide2,'bdycon:psdot')
    end if
    call getmem3d(fg1,jde1ga,jde2ga,ide1ga,ide2ga,1,kzp1,'bdycon:fg1')
    call getmem3d(fg2,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'bdycon:fg2')

  end subroutine allocate_mod_bdycon

  subroutine setup_bdycon
    implicit none
    real(rkx) , dimension(kz) :: anudge
    real(rkx) :: xfun , nb2
    integer(ik4) :: n , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'setup_bdycon'
    integer(ik4) , save :: idindx = 0
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
    rdtbdy = d_one / dtbdys
    if ( iboudy == 1 .or. iboudy >= 5 ) then
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
      if ( idynamic == 3 ) then
        fnudge = fnudge * mo_nadv * mo_nsound
        gnudge = 0.0_rkx
      end if
      if ( myid == italk ) then
        write(stdout, '(a,f12.6,a,f12.6)') &
          ' Nudging coefficients F1=',fnudge,', F2=',gnudge
      end if
    end if
    if ( iboudy == 1 .or. idynamic == 2 ) then
      do n = 2 , nspgx-1
        xfun = real(nspgx-n,rkx)/real(nspgx-2,rkx)
        fcx(n) = fnudge*xfun
        gcx(n) = gnudge*xfun
      end do
      do n = 2 , nspgd-1
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
      do k = 6 , nspgd-1
        wgtd(k) = d_one
      end do
      wgtx(2) = 0.4_rkx
      wgtx(3) = 0.7_rkx
      wgtx(4) = 0.9_rkx
      do k = 5 , nspgx-1
        wgtx(k) = 1.0_rkx
      end do
    end if
    if ( iboudy == 5 ) then
      call exponential_nudging(anudge)
      do k = 1 , kz
        do n = 2 , nspgx-1
          xfun = exp(-(real(n-2,rkx)/anudge(k)))
          hefc(n,k) = fnudge*xfun
          hegc(n,k) = gnudge*xfun
        end do
        do n = 2 , nspgd-1
          xfun = exp(-(real(n-2,rkx)/anudge(k)))
          hefd(n,k) = fnudge*xfun
          hegd(n,k) = gnudge*xfun
        end do
      end do
    end if
    if ( iboudy == 6 ) then
      do k = 1 , kz
        nb2 = d_two * nspgx
        do n = 2 , nspgx-1
          xfun = d_half * (d_one - cos(mathpi/((n-nb2)/nb2)))
          hefc(n,k) = fnudge*xfun
          hegc(n,k) = gnudge*xfun
        end do
        nb2 = d_two * nspgd
        do n = 2 , nspgd-1
          xfun = d_half * (d_one - cos(mathpi/((n-nb2)/nb2)))
          hefd(n,k) = fnudge*xfun
          hegd(n,k) = gnudge*xfun
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine setup_bdycon

  subroutine init_bdy
    implicit none
    integer(ik4) :: datefound , i , j , k , n
    character(len=32) :: appdat
    type (rcm_time_and_date) :: icbc_date
    type (rcm_time_interval) :: tdif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'init_bdy'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

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
      call read_icbc(nhbh0%ps,xtsb%b0,mddom%ldmsk,xub%b0,xvb%b0, &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0)

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh0%ps(j,i) = nhbh0%ps(j,i) * d_r10 - ptop
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh0%tvirt(j,i,k) = xtb%b0(j,i,k)*(d_one+ep1*xqb%b0(j,i,k))
            end do
          end do
        end do
      end if
    else if ( idynamic == 3 ) then
      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,xub%b0,xvb%b0, &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0)
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b0,xqb%b0,xub%b0,xvb%b0,xpsb%b0,xtsb%b0)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b0,xqb%b0,xub%b0,xvb%b0,xpsb%b0,xtsb%b0)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh0%ps(j,i) = xpsb%b0(j,i)
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh0%tvirt(j,i,k) = xtb%b0(j,i,k)*(d_one+ep1*xqb%b0(j,i,k))
            end do
          end do
        end do
      end if
    else
      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,xub%b0,xvb%b0, &
                     xtb%b0,xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0)
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
      xpsb%b0(:,:) = (xpsb%b0(:,:)*d_r10)-ptop
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b0,psdot)
    else if ( idynamic == 2 ) then
      xpsb%b0(:,:) = atm0%ps(:,:) * d_r1000 ! Cb
      psdot(:,:) = atm0%psdot(jde1:jde2,ide1:ide2) * d_r1000
      xpsb%b1(:,:) = xpsb%b0(:,:)
    else
      xpsb%b0(:,:) = xpsb%b0(:,:)*d_100
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Calculate P* on dot points
    ! Couple pressure u,v,t,q (pp,ww)
    !
    if ( idynamic /= 3 ) then
      call couple(xub%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xvb%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
      end if
    end if
    call exchange(xub%b0,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b0,1,jde1,jde2,ide1,ide2,1,kz)
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
      call read_icbc(nhbh1%ps,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh1%ps(j,i) = nhbh1%ps(j,i) * d_r10 - ptop
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
            end do
          end do
        end do
      end if
    else if ( idynamic == 3 ) then
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b1,xqb%b1,xub%b1,xvb%b1,xpsb%b1,xtsb%b1)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b1,xqb%b1,xub%b1,xvb%b1,xpsb%b1,xtsb%b1)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh1%ps(j,i) = xpsb%b1(j,i)
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
            end do
          end do
        end do
      end if
    else
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
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
      write (stdout,*) 'READY  BC from     ' , &
            tochar10(bdydate1) , ' to ' , tochar10(bdydate2)
    end if

    bdydate1 = bdydate2
    !
    ! Repeat for T2
    !
    if ( idynamic == 1 ) then
      xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    else if ( idynamic == 3 ) then
      xpsb%b1(:,:) = xpsb%b1(:,:)*d_100
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
    end if
    call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
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
    end if

    if ( rcmtimer%start( ) ) then
      if ( iseaice == 1 ) then
        if ( islab_ocean == 0 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( mddom%ldmsk(j,i) == 1 ) cycle
              if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
              if ( iocncpl == 1 .or. iwavcpl == 1 ) then
                if ( cplmsk(j,i) /= 0 ) cycle
              end if
              if ( xtsb%b0(j,i) <= icetriggert ) then
                xtsb%b0(j,i) = icetriggert
                mddom%ldmsk(j,i) = 2
                do n = 1 , nnsg
                  if ( mdsub%ldmsk(n,j,i) == 0 ) then
                    mdsub%ldmsk(n,j,i) = 2
                    lms%sfice(n,j,i) = 1.00_rkx
                    lms%sncv(n,j,i) = 1.0_rkx   ! 1 mm of snow over the ice
                    lms%snag(n,j,i) = 0.1_rkx
                  end if
                end do
              end if
            end do
          end do
        end if
      end if
    end if

    if ( iseaice == 1 ) then
      if ( islab_ocean == 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! Update temperatures over ocean water and lakes
            if ( mddom%ldmsk(j,i) == 1 ) cycle
            if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
            if ( iocncpl == 1 .or. iwavcpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            if ( xtsb%b1(j,i) <= icetriggert ) then
              xtsb%b1(j,i) = icetriggert
              mddom%ldmsk(j,i) = 2
              do n = 1 , nnsg
                if ( mdsub%ldmsk(n,j,i) == 0 ) then
                  mdsub%ldmsk(n,j,i) = 2
                  lms%sfice(n,j,i) = 1.00_rkx
                end if
              end do
            else
              if ( mddom%ldmsk(j,i) == 2 ) then
                ! Decrease the surface ice to melt it
                do n = 1 , nnsg
                  if ( mdsub%ldmsk(n,j,i) == 2 ) then
                    lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
                  end if
                end do
              end if
            end if
          end do
        end do
      end if
    end if
    !
    ! Calculate time varying component
    !
    call timeint(xub%b1,xub%b0,xub%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,rdtbdy)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,rdtbdy)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    if ( present_qc ) then
      call timeint(xlb%b1,xlb%b0,xlb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if
    if ( present_qi ) then
      call timeint(xib%b1,xib%b0,xib%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if
    call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2,rdtbdy)
    if ( idynamic == 1 ) then
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga,rdtbdy)
    else if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,rdtbdy)
    else if ( idynamic == 3 ) then
      jday = yeardayfrac(rcmtimer%idate)
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga,rdtbdy)
      call paicompute(mddom%xlat,xpsb%b0,mo_atm%zeta,xtb%b0,xqb%b0,xpaib%b0)
      call paicompute(mddom%xlat,xpsb%b1,mo_atm%zeta,xtb%b1,xqb%b1,xpaib%b1)
      call timeint(xpaib%b1,xpaib%b0,xpaib%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine init_bdy
  !
  ! this subroutine reads in the boundary conditions.
  !
  subroutine bdyin
    implicit none
    integer(ik4) :: i , j , k , n , datefound
    character(len=32) :: appdat
    logical :: update_slabocn
    type (rcm_time_interval) :: tdif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyin'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    update_slabocn = ( islab_ocean == 1 .and. &
      do_qflux_adj .and. som_month /= rcmtimer%month )

    xbctime = d_zero

    xub%b0(:,:,:) = xub%b1(:,:,:)
    xvb%b0(:,:,:) = xvb%b1(:,:,:)
    xtb%b0(:,:,:) = xtb%b1(:,:,:)
    xqb%b0(:,:,:) = xqb%b1(:,:,:)
    if ( present_qc ) then
      xlb%b0(:,:,:) = xlb%b1(:,:,:)
    end if
    if ( present_qi ) then
      xib%b0(:,:,:) = xib%b1(:,:,:)
    end if
    xtsb%b0(:,:) = xtsb%b1(:,:)
    if ( idynamic == 2 ) then
      xppb%b0(:,:,:) = xppb%b1(:,:,:)
      xwwb%b0(:,:,:) = xwwb%b1(:,:,:)
    else
      xpsb%b0(:,:) = xpsb%b1(:,:)
      if ( idynamic == 3 ) then
        xpaib%b0(:,:,:) = xpaib%b1(:,:,:)
      end if
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
      call read_icbc(nhbh1%ps,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh1%ps(j,i) = nhbh1%ps(j,i) * d_r10 - ptop
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
            end do
          end do
        end do
      end if
    else if ( idynamic == 3 ) then
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
      if ( moloch_do_test_1 ) then
        call moloch_static_test1(xtb%b1,xqb%b1,xub%b1,xvb%b1,xpsb%b1,xtsb%b1)
      end if
      if ( moloch_do_test_2 ) then
        call moloch_static_test2(xtb%b1,xqb%b1,xub%b1,xvb%b1,xpsb%b1,xtsb%b1)
      end if

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            nhbh1%ps(j,i) = xpsb%b1(j,i)
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              nhbh1%tvirt(j,i,k) = xtb%b1(j,i,k)*(d_one+ep1*xqb%b1(j,i,k))
            end do
          end do
        end do
      end if
    else
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1, &
                     xtb%b1,xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)
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
      xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    else if ( idynamic == 3 ) then
      xpsb%b1(:,:) = xpsb%b1(:,:)*d_100
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      if ( present_qc ) then
        call couple(xlb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( present_qi ) then
        call couple(xib%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      end if
    end if
    call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
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
    end if

    ! Linear time interpolation
    call timeint(xub%b1,xub%b0,xub%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,rdtbdy)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,rdtbdy)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    if ( present_qc ) then
      call timeint(xlb%b1,xlb%b0,xlb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if
    if ( present_qi ) then
      call timeint(xib%b1,xib%b0,xib%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if
    if ( idynamic == 1 ) then
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga,rdtbdy)
    else if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,rdtbdy)
    else if ( idynamic == 3 ) then
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga,rdtbdy)
      jday = yeardayfrac(rcmtimer%idate)
      call paicompute(mddom%xlat,xpsb%b1,mo_atm%zeta,xtb%b1,xqb%b1,xpaib%b1)
      call timeint(xpaib%b1,xpaib%b0,xpaib%bt, &
                   jce1ga,jce2ga,ice1ga,ice2ga,1,kz,rdtbdy)
    end if
    !
    ! Update ground temperature on Ocean/Lakes
    !
    if ( iseaice == 1 ) then
      if ( islab_ocean == 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! Update temperatures over ocean water only
            if ( mddom%ldmsk(j,i) == 1 ) cycle
            ! Skip lake points if lake model active
            if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
            ! Do not update if coupling and ocean active here
            if ( iocncpl == 1 .or. iwavcpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            ! Sea ice correction
            if ( xtsb%b1(j,i) <= icetriggert ) then
              xtsb%b1(j,i) = icetriggert
              mddom%ldmsk(j,i) = 2
              do n = 1 , nnsg
                if ( mdsub%ldmsk(n,j,i) == 0 ) then
                  mdsub%ldmsk(n,j,i) = 2
                  lms%sfice(n,j,i) = 1.00_rkx
                end if
              end do
            else
              if ( mddom%ldmsk(j,i) == 2 ) then
                ! Decrease the surface ice to melt it
                do n = 1 , nnsg
                  if ( mdsub%ldmsk(n,j,i) == 2 ) then
                    lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
                  end if
                end do
              end if
            end if
          end do
        end do
      end if
    end if

    call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2,rdtbdy)

    if ( myid == italk ) then
      write (stdout,*) 'READY  BC from     ' , &
            tochar10(bdydate1) , ' to ' , tochar10(bdydate2)
    end if

    bdydate1 = bdydate2

    if ( ichem == 1 ) then
      call chem_bdyin
    end if

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
  subroutine bdyuv(xt)
    implicit none
    real(rkx) , intent(in) :: xt
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyuv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Now compute last two points values in U and V
    ! Internal points
    !
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = idi1 , idi2
          wui(i,k) = atm1%u(jdi1,i,k)
          wvi(i,k) = atm1%v(jdi1,i,k)
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = idi1 , idi2
          eui(i,k) = atm1%u(jdi2,i,k)
          evi(i,k) = atm1%v(jdi2,i,k)
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          sui(j,k) = atm1%u(j,idi1,k)
          svi(j,k) = atm1%v(j,idi1,k)
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          nui(j,k) = atm1%u(j,idi2,k)
          nvi(j,k) = atm1%v(j,idi2,k)
        end do
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
        do k = 1 , kz
          do i = ide1 , ide2
            wue(i,k) = xub%b0(jde1,i,k)
            wve(i,k) = xvb%b0(jde1,i,k)
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ide1 , ide2
            eue(i,k) = xub%b0(jde2,i,k)
            eve(i,k) = xvb%b0(jde2,i,k)
          end do
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jde1 , jde2
            sue(j,k) = xub%b0(j,ide1,k)
            sve(j,k) = xvb%b0(j,ide1,k)
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jde1 , jde2
            nue(j,k) = xub%b0(j,ide2,k)
            nve(j,k) = xvb%b0(j,ide2,k)
          end do
        end do
      end if
    else ! NOT Fixed
      !
      !     time-dependent boundary conditions:
      !
      ! west (j = 1) and east (j = jx) boundaries:
      !
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = idi1 , idi2
            wue(i,k) = (xub%b0(jde1,i,k) + xt*xub%bt(jde1,i,k))
            wve(i,k) = (xvb%b0(jde1,i,k) + xt*xvb%bt(jde1,i,k))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = idi1 , idi2
            eue(i,k) = (xub%b0(jde2,i,k) + xt*xub%bt(jde2,i,k))
            eve(i,k) = (xvb%b0(jde2,i,k) + xt*xvb%bt(jde2,i,k))
          end do
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jde1 , jde2
            sue(j,k) = (xub%b0(j,ide1,k) + xt*xub%bt(j,ide1,k))
            sve(j,k) = (xvb%b0(j,ide1,k) + xt*xvb%bt(j,ide1,k))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jde1 , jde2
            nue(j,k) = (xub%b0(j,ide2,k) + xt*xub%bt(j,ide2,k))
            nve(j,k) = (xvb%b0(j,ide2,k) + xt*xvb%bt(j,ide2,k))
          end do
        end do
      end if
    end if
    !
    ! fill up the interior silces:
    !
    if ( ma%has_bdytopleft ) then
      do k = 1 , kz
        wui(ide2,k) = nue(jdi1,k)
        wvi(ide2,k) = nve(jdi1,k)
        nui(jde1,k) = wue(idi2,k)
        nvi(jde1,k) = wve(idi2,k)
      end do
    end if
    if ( ma%has_bdybottomleft ) then
      do k = 1 , kz
        wui(ide1,k) = sue(jdi1,k)
        wvi(ide1,k) = sve(jdi1,k)
        sui(jde1,k) = wue(idi1,k)
        svi(jde1,k) = wve(idi1,k)
      end do
    end if
    if ( ma%has_bdytopright ) then
      do k = 1 , kz
        eui(ide2,k) = nue(jdi2,k)
        evi(ide2,k) = nve(jdi2,k)
        nui(jde2,k) = eue(idi2,k)
        nvi(jde2,k) = eve(idi2,k)
      end do
    end if
    if ( ma%has_bdybottomright ) then
      do k = 1 , kz
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
    real(rkx) :: qxint , qrat , tkeint , qext , qint
    integer(ik4) :: i , j , k , n
    real(rkx) :: windavg
    real(rkx) :: xt
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyval'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Fill up the boundary value for xxb variables from xxa variables:
    ! if this subroutine is called for the first time, this part
    ! shall be skipped.
    !
    xt = xbctime + dt
    if ( rcmtimer%integrating( ) ) then
      if ( idynamic /= 3 ) then
        !
        ! West boundary
        !
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm2%u(jde1,i,k) = atm1%u(jde1,i,k)
              atm2%v(jde1,i,k) = atm1%v(jde1,i,k)
            end do
          end do
          do k = 1 , kz
            do i = ici1 , ici2
              atm2%t(jce1,i,k) = atm1%t(jce1,i,k)
            end do
          end do
          do n = 1 , nqx
            do k = 1 , kz
              do i = ici1 , ici2
                atm2%qx(jce1,i,k,n) = atm1%qx(jce1,i,k,n)
              end do
            end do
          end do
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm2%pp(jce1,i,k) = atm1%pp(jce1,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm2%w(jce1,i,k) = atm1%w(jce1,i,k)
              end do
            end do
          else
            do i = ici1 , ici2
              sfs%psb(jce1,i) = sfs%psa(jce1,i)
            end do
          end if
          if ( ibltyp == 2 ) then
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm2%tke(jce1,i,k) = atm1%tke(jce1,i,k)
              end do
            end do
          end if
        end if
        !
        ! East boundary
        !
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm2%u(jde2,i,k) = atm1%u(jde2,i,k)
              atm2%v(jde2,i,k) = atm1%v(jde2,i,k)
            end do
          end do
          do k = 1 , kz
            do i = ici1 , ici2
              atm2%t(jce2,i,k) = atm1%t(jce2,i,k)
            end do
          end do
          do n = 1 , nqx
            do k = 1 , kz
              do i = ici1 , ici2
                atm2%qx(jce2,i,k,n) = atm1%qx(jce2,i,k,n)
              end do
            end do
          end do
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm2%pp(jce2,i,k) = atm1%pp(jce2,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm2%w(jce2,i,k) = atm1%w(jce2,i,k)
              end do
            end do
          else
            do i = ici1 , ici2
              sfs%psb(jce2,i) = sfs%psa(jce2,i)
            end do
          end if
          if ( ibltyp == 2 ) then
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm2%tke(jce2,i,k) = atm1%tke(jce2,i,k)
              end do
            end do
          end if
        end if
        !
        ! North and South boundaries
        !
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm2%u(j,ide1,k) = atm1%u(j,ide1,k)
              atm2%v(j,ide1,k) = atm1%v(j,ide1,k)
            end do
          end do
          do k = 1 , kz
            do j = jce1 , jce2
              atm2%t(j,ice1,k) = atm1%t(j,ice1,k)
            end do
          end do
          do n = 1 , nqx
            do k = 1 , kz
              do j = jce1 , jce2
                atm2%qx(j,ice1,k,n) = atm1%qx(j,ice1,k,n)
              end do
            end do
          end do
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm2%pp(j,ice1,k) = atm1%pp(j,ice1,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm2%w(j,ice1,k) = atm1%w(j,ice1,k)
              end do
            end do
          else
            do j = jce1 , jce2
              sfs%psb(j,ice1) = sfs%psa(j,ice1)
            end do
          end if
          if ( ibltyp == 2 ) then
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm2%tke(j,ice1,k) = atm1%tke(j,ice1,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm2%u(j,ide2,k) = atm1%u(j,ide2,k)
              atm2%v(j,ide2,k) = atm1%v(j,ide2,k)
            end do
          end do
          do k = 1 , kz
            do j = jce1 , jce2
              atm2%t(j,ice2,k) = atm1%t(j,ice2,k)
            end do
          end do
          do n = 1 , nqx
            do k = 1 , kz
              do j = jce1 , jce2
                atm2%qx(j,ice2,k,n) = atm1%qx(j,ice2,k,n)
              end do
            end do
          end do
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm2%pp(j,ice2,k) = atm1%pp(j,ice2,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm2%w(j,ice2,k) = atm1%w(j,ice2,k)
              end do
            end do
          else
            do j = jce1 , jce2
              sfs%psb(j,ice2) = sfs%psa(j,ice2)
            end do
          end if
          if ( ibltyp == 2 ) then
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm2%tke(j,ice2,k) = atm1%tke(j,ice2,k)
              end do
            end do
          end if
        end if
      end if
    end if
    !
    ! Compute the boundary values for xxa variables:
    !
    ! Set boundary values for p*:
    ! Set boundary conditions for p*u and p*v:
    !
    if ( iboudy == 0 ) then
      !
      ! fixed boundary conditions:
      !
      if ( idynamic == 1 .or. idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do i = ici1 , ici2
            sfs%psa(jce1,i) = xpsb%b0(jce1,i)
          end do
        end if
        if ( ma%has_bdyright ) then
          do i = ici1 , ici2
            sfs%psa(jce2,i) = xpsb%b0(jce2,i)
          end do
        end if
        if ( ma%has_bdybottom ) then
          do j = jce1 , jce2
            sfs%psa(j,ice1) = xpsb%b0(j,ice1)
          end do
        end if
        if ( ma%has_bdytop ) then
          do j = jce1 , jce2
            sfs%psa(j,ice2) = xpsb%b0(j,ice2)
          end do
        end if
      end if
      if ( idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ice1 , ice2
              mo_atm%u(jde1,i,k) = xub%b0(jde1,i,k)
            end do
          end do
          do k = 1 , kz
            do i = ide1 , ide2
              mo_atm%v(jce1,i,k) = xvb%b0(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ice1 , ice2
              mo_atm%u(jde2,i,k) = xub%b0(jde2,i,k)
            end do
          end do
          do k = 1 , kz
            do i = ide1 , ide2
              mo_atm%v(jce2,i,k) = xvb%b0(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jde1 , jde2
              mo_atm%u(j,ice1,k) = xub%b0(j,ice1,k)
            end do
          end do
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%v(j,ide1,k) = xvb%b0(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jde1 , jde2
              mo_atm%u(j,ice2,k) = xub%b0(j,ice2,k)
            end do
          end do
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%v(j,ide2,k) = xvb%b0(j,ide2,k)
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm1%u(jde1,i,k) = xub%b0(jde1,i,k)
              atm1%v(jde1,i,k) = xvb%b0(jde1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm1%u(jde2,i,k) = xub%b0(jde2,i,k)
              atm1%v(jde2,i,k) = xvb%b0(jde2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm1%u(j,ide1,k) = xub%b0(j,ide1,k)
              atm1%v(j,ide1,k) = xvb%b0(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm1%u(j,ide2,k) = xub%b0(j,ide2,k)
              atm1%v(j,ide2,k) = xvb%b0(j,ide2,k)
            end do
          end do
        end if
      end if
    else
      !
      ! time-dependent boundary conditions:
      !
      if ( idynamic == 1 .or. idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do i = ici1 , ici2
            sfs%psa(jce1,i) = xpsb%b0(jce1,i) + xt*xpsb%bt(jce1,i)
          end do
        end if
        if ( ma%has_bdyright ) then
          do i = ici1 , ici2
            sfs%psa(jce2,i) = xpsb%b0(jce2,i) + xt*xpsb%bt(jce2,i)
          end do
        end if
        if ( ma%has_bdybottom ) then
          do j = jce1 , jce2
            sfs%psa(j,ice1) = xpsb%b0(j,ice1) + xt*xpsb%bt(j,ice1)
          end do
        end if
        if ( ma%has_bdytop ) then
          do j = jce1 , jce2
            sfs%psa(j,ice2) = xpsb%b0(j,ice2) + xt*xpsb%bt(j,ice2)
          end do
        end if
      end if
      if ( idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ice1 , ice2
              mo_atm%u(jde1,i,k) = xub%b0(jde1,i,k) + xt*xub%bt(jde1,i,k)
            end do
            do i = ide1 , ide2
              mo_atm%v(jce1,i,k) = xvb%b0(jce1,i,k) + xt*xvb%bt(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ice1 , ice2
              mo_atm%u(jde2,i,k) = xub%b0(jde2,i,k) + xt*xub%bt(jde2,i,k)
            end do
            do i = ide1 , ide2
              mo_atm%v(jce2,i,k) = xvb%b0(jce2,i,k) + xt*xvb%bt(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jde1 , jde2
              mo_atm%u(j,ice1,k) = xub%b0(j,ice1,k) + xt*xub%bt(j,ice1,k)
            end do
            do j = jce1 , jce2
              mo_atm%v(j,ide1,k) = xvb%b0(j,ide1,k) + xt*xvb%bt(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jde1 , jde2
              mo_atm%u(j,ice2,k) = xub%b0(j,ice2,k) + xt*xub%bt(j,ice2,k)
            end do
            do j = jce1 , jce2
              mo_atm%v(j,ide2,k) = xvb%b0(j,ide2,k) + xt*xvb%bt(j,ide2,k)
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm1%u(jde1,i,k) = xub%b0(jde1,i,k) + xt*xub%bt(jde1,i,k)
              atm1%v(jde1,i,k) = xvb%b0(jde1,i,k) + xt*xvb%bt(jde1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = idi1 , idi2
              atm1%u(jde2,i,k) = xub%b0(jde2,i,k) + xt*xub%bt(jde2,i,k)
              atm1%v(jde2,i,k) = xvb%b0(jde2,i,k) + xt*xvb%bt(jde2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm1%u(j,ide1,k) = xub%b0(j,ide1,k) + xt*xub%bt(j,ide1,k)
              atm1%v(j,ide1,k) = xvb%b0(j,ide1,k) + xt*xvb%bt(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jde1 , jde2
              atm1%u(j,ide2,k) = xub%b0(j,ide2,k) + xt*xub%bt(j,ide2,k)
              atm1%v(j,ide2,k) = xvb%b0(j,ide2,k) + xt*xvb%bt(j,ide2,k)
            end do
          end do
        end if
      end if
    end if

    if ( idynamic /= 3 ) call bdyuv(xt)

    !
    ! Set boundary values for p*t:
    ! Set boundary values for p*qv:
    !
    if ( iboudy == 0 ) then
      !
      ! fixed boundary conditions:
      !
      if ( idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce1,i,k) = xtb%b0(jce1,i,k)
              mo_atm%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k)
              mo_atm%pai(jce1,i,k) = xpaib%b0(jce1,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce1,i,k,iqc) = xlb%b0(jce1,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce1,i,k,iqi) = xib%b0(jce1,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce2,i,k) = xtb%b0(jce2,i,k)
              mo_atm%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k)
              mo_atm%pai(jce2,i,k) = xpaib%b0(jce2,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce2,i,k,iqc) = xlb%b0(jce2,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce2,i,k,iqi) = xib%b0(jce2,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice1,k) = xtb%b0(j,ice1,k)
              mo_atm%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k)
              mo_atm%pai(j,ice1,k) = xpaib%b0(j,ice1,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice1,k,iqc) = xlb%b0(j,ice1,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice1,k,iqi) = xib%b0(j,ice1,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice2,k) = xtb%b0(j,ice2,k)
              mo_atm%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k)
              mo_atm%pai(j,ice2,k) = xpaib%b0(j,ice2,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice2,k,iqc) = xlb%b0(j,ice2,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice2,k,iqi) = xib%b0(j,ice2,k)
              end do
            end do
          end if
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce1,i,k) = xtb%b0(jce1,i,k)
              atm1%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce1,i,k,iqv) = xlb%b0(jce1,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce1,i,k,iqv) = xlb%b0(jce1,i,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%pp(jce1,i,k) = xppb%b0(jce1,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm1%w(jce1,i,k) = xwwb%b0(jce1,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce2,i,k) = xtb%b0(jce2,i,k)
              atm1%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce2,i,k,iqc) = xlb%b0(jce2,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce2,i,k,iqi) = xib%b0(jce2,i,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%pp(jce2,i,k) = xppb%b0(jce2,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm1%w(jce2,i,k) = xwwb%b0(jce2,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              atm1%t(j,ice1,k) = xtb%b0(j,ice1,k)
              atm1%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice1,k,iqc) = xlb%b0(j,ice1,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice1,k,iqi) = xib%b0(j,ice1,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%pp(j,ice1,k) = xppb%b0(j,ice1,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm1%w(j,ice1,k) = xwwb%b0(j,ice1,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              atm1%t(j,ice2,k) = xtb%b0(j,ice2,k)
              atm1%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice2,k,iqc) = xlb%b0(j,ice2,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice2,k,iqi) = xib%b0(j,ice2,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%pp(j,ice2,k) = xppb%b0(j,ice2,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm1%w(j,ice2,k) = xwwb%b0(j,ice2,k)
              end do
            end do
          end if
        end if
      end if
    else
      !
      ! time-dependent boundary conditions:
      !
      if ( idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce1,i,k) = xtb%b0(jce1,i,k) + xt*xtb%bt(jce1,i,k)
              mo_atm%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k) + xt*xqb%bt(jce1,i,k)
              mo_atm%pai(jce1,i,k) = xpaib%b0(jce1,i,k) + xt*xpaib%bt(jce1,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce1,i,k,iqc) = xlb%b0(jce1,i,k) + xt*xlb%bt(jce1,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce1,i,k,iqi) = xib%b0(jce1,i,k) + xt*xib%bt(jce1,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce2,i,k) = xtb%b0(jce2,i,k) + xt*xtb%bt(jce2,i,k)
              mo_atm%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k) + xt*xqb%bt(jce2,i,k)
              mo_atm%pai(jce2,i,k) = xpaib%b0(jce2,i,k) + xt*xpaib%bt(jce2,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce2,i,k,iqc) = xlb%b0(jce2,i,k) + xt*xlb%bt(jce2,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                mo_atm%qx(jce2,i,k,iqi) = xib%b0(jce2,i,k) + xt*xib%bt(jce2,i,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice1,k) = xtb%b0(j,ice1,k) + xt*xtb%bt(j,ice1,k)
              mo_atm%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k) + xt*xqb%bt(j,ice1,k)
              mo_atm%pai(j,ice1,k) = xpaib%b0(j,ice1,k) + xt*xpaib%bt(j,ice1,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice1,k,iqc) = xlb%b0(j,ice1,k) + xt*xlb%bt(j,ice1,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice1,k,iqi) = xib%b0(j,ice1,k) + xt*xib%bt(j,ice1,k)
              end do
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice2,k) = xtb%b0(j,ice2,k) + xt*xtb%bt(j,ice2,k)
              mo_atm%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k) + xt*xqb%bt(j,ice2,k)
              mo_atm%pai(j,ice2,k) = xpaib%b0(j,ice2,k) + xt*xpaib%bt(j,ice2,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice2,k,iqc) = xlb%b0(j,ice2,k) + xt*xlb%bt(j,ice2,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                mo_atm%qx(j,ice2,k,iqi) = xib%b0(j,ice2,k) + xt*xib%bt(j,ice2,k)
              end do
            end do
          end if
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce1,i,k)      = xtb%b0(jce1,i,k) + xt*xtb%bt(jce1,i,k)
              atm1%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k) + xt*xqb%bt(jce1,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce1,i,k,iqc) = xlb%b0(jce1,i,k) + xt*xlb%bt(jce1,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce1,i,k,iqi) = xib%b0(jce1,i,k) + xt*xib%bt(jce1,i,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%pp(jce1,i,k) = xppb%b0(jce1,i,k) + xt*xppb%bt(jce1,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm1%w(jce1,i,k) = xwwb%b0(jce1,i,k) + xt*xwwb%bt(jce1,i,k)
              end do
            end do
            do i = ici1 , ici2
              atm1%w(jce1,i,1) = atm1%w(jci1,i,1)
            end do
          end if
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce2,i,k)      = xtb%b0(jce2,i,k) + xt*xtb%bt(jce2,i,k)
              atm1%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k) + xt*xqb%bt(jce2,i,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce2,i,k,iqc) = xlb%b0(jce2,i,k) + xt*xlb%bt(jce2,i,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%qx(jce2,i,k,iqi) = xib%b0(jce2,i,k) + xt*xib%bt(jce2,i,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                atm1%pp(jce2,i,k) = xppb%b0(jce2,i,k) + xt*xppb%bt(jce2,i,k)
              end do
            end do
            do k = 1 , kzp1
              do i = ici1 , ici2
                atm1%w(jce2,i,k) = xwwb%b0(jce2,i,k) + xt*xwwb%bt(jce2,i,k)
              end do
            end do
            do i = ici1 , ici2
              atm1%w(jce2,i,1) = atm1%w(jci2,i,1)
            end do
          end if
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              atm1%t(j,ice1,k)      = xtb%b0(j,ice1,k) + xt*xtb%bt(j,ice1,k)
              atm1%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k) + xt*xqb%bt(j,ice1,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice1,k,iqc) = xlb%b0(j,ice1,k) + xt*xlb%bt(j,ice1,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice1,k,iqi) = xib%b0(j,ice1,k) + xt*xib%bt(j,ice1,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%pp(j,ice1,k) = xppb%b0(j,ice1,k) + xt*xppb%bt(j,ice1,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm1%w(j,ice1,k) = xwwb%b0(j,ice1,k) + xt*xwwb%bt(j,ice1,k)
              end do
            end do
            do j = jce1 , jce2
              atm1%w(j,ice1,1) = atm1%w(j,ici1,1)
            end do
          end if
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              atm1%t(j,ice2,k)      = xtb%b0(j,ice2,k) + xt*xtb%bt(j,ice2,k)
              atm1%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k) + xt*xqb%bt(j,ice2,k)
            end do
          end do
          if ( present_qc ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice2,k,iqc) = xlb%b0(j,ice2,k) + xt*xlb%bt(j,ice2,k)
              end do
            end do
          end if
          if ( present_qi .and. ipptls > 1 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%qx(j,ice2,k,iqi) = xib%b0(j,ice2,k) + xt*xib%bt(j,ice2,k)
              end do
            end do
          end if
          if ( idynamic == 2 ) then
            do k = 1 , kz
              do j = jce1 , jce2
                atm1%pp(j,ice2,k) = xppb%b0(j,ice2,k) + xt*xppb%bt(j,ice2,k)
              end do
            end do
            do k = 1 , kzp1
              do j = jce1 , jce2
                atm1%w(j,ice2,k) = xwwb%b0(j,ice2,k) + xt*xwwb%bt(j,ice2,k)
              end do
            end do
            do j = jce1 , jce2
              atm1%w(j,ice2,1) = atm1%w(j,ici2,1)
            end do
          end if
        end if
      end if
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
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
        sfs%tg(j,i) = xtsb%b0(j,i) + xt*xtsb%bt(j,i)
      end do
    end do

    if ( iboudy == 3 .or. iboudy == 4 ) then
      !
      ! determine QV boundary values depends on inflow/outflow:
      !
      if ( idynamic == 3 ) then
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              qext = mo_atm%qx(jce1,i,k,iqv)
              qint = mo_atm%qx(jci1,i,k,iqv)
              windavg = mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k)
              if ( windavg > d_zero ) then
                mo_atm%qx(jce1,i,k,iqv) = qext
              else
                mo_atm%qx(jce1,i,k,iqv) = qint
              end if
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              qext = mo_atm%qx(jce2,i,k,iqv)
              qint = mo_atm%qx(jci2,i,k,iqv)
              windavg = mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k)
              if ( windavg < d_zero ) then
                mo_atm%qx(jce2,i,k,iqv) = qext
              else
                mo_atm%qx(jce2,i,k,iqv) = qint
              end if
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              qext = mo_atm%qx(j,ice1,k,iqv)
              qint = mo_atm%qx(j,ici1,k,iqv)
              windavg = mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k)
              if ( windavg > d_zero ) then
                mo_atm%qx(j,ice1,k,iqv) = qext
              else
                mo_atm%qx(j,ice1,k,iqv) = qint
              end if
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              qext = mo_atm%qx(j,ice2,k,iqv)
              qint = mo_atm%qx(j,ici2,k,iqv)
              windavg = mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k)
              if ( windavg < d_zero ) then
                mo_atm%qx(j,ice2,k,iqv) = qext
              else
                mo_atm%qx(j,ice2,k,iqv) = qint
              end if
            end do
          end do
        end if
      else
        !
        ! west boundary:
        !
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              qext = atm1%qx(jce1,i,k,iqv)/sfs%psa(jce1,i)
              qint = atm1%qx(jci1,i,k,iqv)/sfs%psa(jci1,i)
              windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
              if ( windavg > d_zero ) then
                atm1%qx(jce1,i,k,iqv) = qext*sfs%psa(jce1,i)
              else
                atm1%qx(jce1,i,k,iqv) = qint*sfs%psa(jce1,i)
              end if
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              qext = atm1%qx(jce2,i,k,iqv)/sfs%psa(jce2,i)
              qint = atm1%qx(jci2,i,k,iqv)/sfs%psa(jci2,i)
              windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
              if ( windavg < d_zero ) then
                atm1%qx(jce2,i,k,iqv) = qext*sfs%psa(jce2,i)
              else
                atm1%qx(jce2,i,k,iqv) = qint*sfs%psa(jce2,i)
              end if
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              qext = atm1%qx(j,ice1,k,iqv)/sfs%psa(j,ice1)
              qint = atm1%qx(j,ici1,k,iqv)/sfs%psa(j,ici1)
              windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
              if ( windavg > d_zero ) then
                atm1%qx(j,ice1,k,iqv) = qext*sfs%psa(j,ice1)
              else
                atm1%qx(j,ice1,k,iqv) = qint*sfs%psa(j,ice1)
              end if
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              qext = atm1%qx(j,ice2,k,iqv)/sfs%psa(j,ice2)
              qint = atm1%qx(j,ici2,k,iqv)/sfs%psa(j,ici2)
              windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
              if ( windavg < d_zero ) then
                atm1%qx(j,ice2,k,iqv) = qext*sfs%psa(j,ice2)
              else
                atm1%qx(j,ice2,k,iqv) = qint*sfs%psa(j,ice2)
              end if
            end do
          end do
        end if
      end if
    end if

    ! set boundary values for p*qx
    ! *** note ***
    ! for large domain, we assume the boundary tendencies are not available.
    !
    ! if the boundary values and tendencies are not available,
    ! determine boundary values depends on inflow/outflow:
    ! inflow  : set it equal to zero.
    ! outflow : get from interior point.
    ! west boundary:
    if ( idynamic == 3 ) then
      if ( bdyflow ) then
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = max(mo_atm%qx(jci1,i,k,n),d_zero)
                windavg = (mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k))
                if ( windavg > d_zero ) then
                  mo_atm%qx(jce1,i,k,n) = d_zero
                else
                  mo_atm%qx(jce1,i,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = max(mo_atm%qx(jci2,i,k,n),d_zero)
                windavg = (mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k))
                if ( windavg < d_zero ) then
                  mo_atm%qx(jce2,i,k,n) = d_zero
                else
                  mo_atm%qx(jce2,i,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = max(mo_atm%qx(j,ici1,k,n),d_zero)
                windavg = (mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k))
                if ( windavg > d_zero ) then
                  mo_atm%qx(j,ice1,k,n) = d_zero
                else
                  mo_atm%qx(j,ice1,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst , iqlst
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = max(mo_atm%qx(j,ici2,k,n),d_zero)
                windavg = (mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k))
                if ( windavg < d_zero ) then
                  mo_atm%qx(j,ice2,k,n) = d_zero
                else
                  mo_atm%qx(j,ice2,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = mo_atm%qx(jci1,i,k,n)
                qrat  = mo_atm%qx(jce1,i,k,iqv)/mo_atm%qx(jci1,i,k,iqv)
                mo_atm%qx(jce1,i,k,n) = qxint*qrat
              end do
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = mo_atm%qx(jci2,i,k,n)
                qrat  = mo_atm%qx(jce2,i,k,iqv)/mo_atm%qx(jci2,i,k,iqv)
                mo_atm%qx(jce2,i,k,n) = qxint*qrat
              end do
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = mo_atm%qx(j,ici1,k,n)
                qrat  = mo_atm%qx(j,ice1,k,iqv)/mo_atm%qx(j,ici1,k,iqv)
                mo_atm%qx(j,ice1,k,n) = qxint*qrat
              end do
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = mo_atm%qx(j,ici2,k,n)
                qrat  = mo_atm%qx(j,ice2,k,iqv)/mo_atm%qx(j,ici2,k,iqv)
                mo_atm%qx(j,ice2,k,n) = qxint*qrat
              end do
            end do
          end do
        end if
      end if
    else
      if ( bdyflow ) then
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = atm1%qx(jci1,i,k,n)
                windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
                if ( windavg > d_zero ) then
                  atm1%qx(jce1,i,k,n) = d_zero
                else
                  atm1%qx(jce1,i,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = atm1%qx(jci2,i,k,n)
                windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
                if ( windavg < d_zero ) then
                  atm1%qx(jce2,i,k,n) = d_zero
                else
                  atm1%qx(jce2,i,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = atm1%qx(j,ici1,k,n)
                windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
                if ( windavg > d_zero ) then
                  atm1%qx(j,ice1,k,n) = d_zero
                else
                  atm1%qx(j,ice1,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = atm1%qx(j,ici2,k,n)
                windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
                if ( windavg < d_zero ) then
                  atm1%qx(j,ice2,k,n) = d_zero
                else
                  atm1%qx(j,ice2,k,n) = qxint
                end if
              end do
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = atm1%qx(jci1,i,k,n)/sfs%psa(jci1,i)
                qrat  = atm1%qx(jce1,i,k,iqv)/atm1%qx(jci1,i,k,iqv)
                atm1%qx(jce1,i,k,n) = qxint*sfs%psa(jce1,i)*qrat
              end do
            end do
          end do
        end if
        !
        ! east boundary:
        !
        if ( ma%has_bdyright ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do i = ici1 , ici2
                qxint = atm1%qx(jci2,i,k,n)/sfs%psa(jci2,i)
                qrat  = atm1%qx(jce2,i,k,iqv)/atm1%qx(jci2,i,k,iqv)
                atm1%qx(jce2,i,k,n) = qxint*sfs%psa(jce2,i)*qrat
              end do
            end do
          end do
        end if
        !
        ! south boundary:
        !
        if ( ma%has_bdybottom ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = atm1%qx(j,ici1,k,n)/sfs%psa(j,ici1)
                qrat  = atm1%qx(j,ice1,k,iqv)/atm1%qx(j,ici1,k,iqv)
                atm1%qx(j,ice1,k,n) = qxint*sfs%psa(j,ice1)*qrat
              end do
            end do
          end do
        end if
        !
        ! north boundary:
        !
        if ( ma%has_bdytop ) then
          do n = iqfrst , iqlst
            if ( present_qc .and. n == iqc ) cycle
            if ( present_qi .and. n == iqi ) cycle
            do k = 1 , kz
              do j = jce1 , jce2
                qxint = atm1%qx(j,ici2,k,n)/sfs%psa(j,ici2)
                qrat  = atm1%qx(j,ice2,k,iqv)/atm1%qx(j,ici2,k,iqv)
                atm1%qx(j,ice2,k,n) = qxint*sfs%psa(j,ice2)*qrat
              end do
            end do
          end do
        end if
      end if
    end if

    if ( ibltyp == 2 ) then
      if ( idynamic == 3 ) then
        if ( rcmtimer%start( ) ) then
          if ( ma%has_bdyleft ) then
            mo_atm%tke(jce1,:,:) = tkemin ! East boundary
          end if
          if ( ma%has_bdyright ) then
            mo_atm%tke(jce2,:,:) = tkemin ! West boundary
          end if
          if ( ma%has_bdytop ) then
            mo_atm%tke(:,ice2,:) = tkemin  ! North boundary
          end if
          if ( ma%has_bdybottom ) then
            mo_atm%tke(:,ice1,:) = tkemin  ! South boundary
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
              mo_atm%tke(jce1,:,1) = tkemin ! West boundary
              do k = 2 , kz
                do i = ice1 , ice2
                  tkeint = mo_atm%tke(jci1,i,k+1)
                  windavg = mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k) + &
                            mo_atm%u(jde1,i,k-1) + mo_atm%u(jdi1,i,k-1)
                  if ( windavg > d_zero ) then
                    mo_atm%tke(jce1,i,k+1) = tkemin
                  else
                    mo_atm%tke(jce1,i,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              mo_atm%tke(jce2,:,1) = tkemin ! East boundary
              do k = 2 , kz
                do i = ice1 , ice2
                  tkeint = mo_atm%tke(jci2,i,k+1)
                  windavg = mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k) + &
                            mo_atm%u(jde2,i,k-1) + mo_atm%u(jdi2,i,k-1)
                  if ( windavg < d_zero ) then
                    mo_atm%tke(jce2,i,k+1) = tkemin
                  else
                    mo_atm%tke(jce2,i,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              mo_atm%tke(:,ice1,1) = tkemin  ! South boundary
              do k = 2 , kz
                do j = jci1 , jci2
                  tkeint = mo_atm%tke(j,ici1,k+1)
                  windavg = mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k) + &
                            mo_atm%v(j,ide1,k-1) + mo_atm%v(j,idi1,k-1)
                  if ( windavg > d_zero ) then
                    mo_atm%tke(j,ice1,k+1) = tkemin
                  else
                    mo_atm%tke(j,ice1,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              mo_atm%tke(:,ice2,1) = tkemin  ! North boundary
              do k = 2 , kz
                do j = jci1 , jci2
                  tkeint = mo_atm%tke(j,ici2,k+1)
                  windavg = mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k) + &
                            mo_atm%v(j,ide2,k-1) + mo_atm%v(j,idi2,k-1)
                  if ( windavg < d_zero ) then
                    mo_atm%tke(j,ice2,k+1) = tkemin
                  else
                    mo_atm%tke(j,ice2,k+1) = tkeint
                  end if
                end do
              end do
            end if
          else
            if ( ma%has_bdyleft ) then
              do k = 1 , kzp1
                do i = ice1 , ice2
                  mo_atm%tke(jce1,i,k) = mo_atm%tke(jci1,i,k)
                end do
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              do k = 1 , kzp1
                do i = ice1 , ice2
                  mo_atm%tke(jce2,i,k) = mo_atm%tke(jci2,i,k)
                end do
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              do k = 1 , kzp1
                do j = jci1 , jci2
                  mo_atm%tke(j,ice1,k) = mo_atm%tke(j,ici1,k)
                end do
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              do k = 1 , kzp1
                do j = jci1 , jci2
                  mo_atm%tke(j,ice2,k) = mo_atm%tke(j,ici2,k)
                end do
              end do
            end if
          end if
        end if
      else
        if ( rcmtimer%start( ) ) then
          if ( ma%has_bdyleft ) then
            atm1%tke(jce1,:,:) = tkemin ! East boundary
            atm2%tke(jce1,:,:) = tkemin ! East boundary
          end if
          if ( ma%has_bdyright ) then
            atm1%tke(jce2,:,:) = tkemin ! West boundary
            atm2%tke(jce2,:,:) = tkemin ! West boundary
          end if
          if ( ma%has_bdytop ) then
            atm1%tke(:,ice2,:) = tkemin  ! North boundary
            atm2%tke(:,ice2,:) = tkemin  ! North boundary
          end if
          if ( ma%has_bdybottom ) then
            atm1%tke(:,ice1,:) = tkemin  ! South boundary
            atm2%tke(:,ice1,:) = tkemin  ! South boundary
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
              atm1%tke(jce1,:,1) = tkemin ! West boundary
              atm2%tke(jce1,:,1) = tkemin ! West boundary
              do k = 2 , kz
                do i = ice1 , ice2
                  tkeint = atm1%tke(jci1,i,k+1)
                  windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k) + &
                    wue(i,k-1) + wue(i+1,k-1) + wui(i,k-1) + wui(i+1,k-1)
                  if ( windavg > d_zero ) then
                    atm1%tke(jce1,i,k+1) = tkemin
                  else
                    atm1%tke(jce1,i,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              atm1%tke(jce2,:,1) = tkemin ! East boundary
              atm2%tke(jce2,:,1) = tkemin ! East boundary
              do k = 2 , kz
                do i = ice1 , ice2
                  tkeint = atm1%tke(jci2,i,k+1)
                  windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k) + &
                    eue(i,k-1) + eue(i+1,k-1) + eui(i,k-1) + eui(i+1,k-1)
                  if ( windavg < d_zero ) then
                    atm1%tke(jce2,i,k+1) = tkemin
                  else
                    atm1%tke(jce2,i,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              atm1%tke(:,ice1,1) = tkemin  ! South boundary
              atm2%tke(:,ice1,1) = tkemin  ! South boundary
              do k = 2 , kz
                do j = jci1 , jci2
                  tkeint = atm1%tke(j,ici1,k+1)
                  windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k) + &
                    sve(j,k-1) + sve(j+1,k-1) + svi(j,k-1) + svi(j+1,k-1)
                  if ( windavg > d_zero ) then
                    atm1%tke(j,ice1,k+1) = tkemin
                  else
                    atm1%tke(j,ice1,k+1) = tkeint
                  end if
                end do
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              atm1%tke(:,ice2,1) = tkemin  ! North boundary
              atm2%tke(:,ice2,1) = tkemin  ! North boundary
              do k = 2 , kz
                do j = jci1 , jci2
                  tkeint = atm1%tke(j,ici2,k+1)
                  windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k) + &
                    nve(j,k-1) + nve(j+1,k-1) + nvi(j,k-1) + nvi(j+1,k-1)
                  if ( windavg < d_zero ) then
                    atm1%tke(j,ice2,k+1) = tkemin
                  else
                    atm1%tke(j,ice2,k+1) = tkeint
                  end if
                end do
              end do
            end if
          else
            if ( ma%has_bdyleft ) then
              do k = 1 , kzp1
                do i = ice1 , ice2
                  atm1%tke(jce1,i,k) = atm1%tke(jci1,i,k)
                end do
              end do
            end if
            !
            ! east boundary:
            !
            if ( ma%has_bdyright ) then
              do k = 1 , kzp1
                do i = ice1 , ice2
                  atm1%tke(jce2,i,k) = atm1%tke(jci2,i,k)
                end do
              end do
            end if
            !
            ! south boundary:
            !
            if ( ma%has_bdybottom ) then
              do k = 1 , kzp1
                do j = jci1 , jci2
                  atm1%tke(j,ice1,k) = atm1%tke(j,ici1,k)
                end do
              end do
            end if
            !
            ! north boundary:
            !
            if ( ma%has_bdytop ) then
              do k = 1 , kzp1
                do j = jci1 , jci2
                  atm1%tke(j,ice2,k) = atm1%tke(j,ici2,k)
                end do
              end do
            end if
          end if
        end if
      end if
    end if

    if ( ichem == 1 ) then
      if ( idynamic == 3 ) then
        call chem_bdyval(mo_atm%u,mo_atm%v)
      else
        call chem_bdyval(sfs%psa,wue,wui,eue,eui,nve,nvi,sve,svi)
      end if
    end if

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
    integer(ik4) , intent(in) :: m
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: ften

    integer(ik4) :: i , j , k
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_cr%ns /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m) + &
                            (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m) + &
                            (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m) + &
                            (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k,m) = wgtx(ib)*ften(j,i,k,m) + &
                            (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge4d

  subroutine mosponge4d(f,bnd,m)
    implicit none
    integer(ik4) , intent(in) :: m
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd

    integer(ik4) :: i , j , k
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mosponge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_cr%ns /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k,m) = f(j,i,k,m) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k,m) = f(j,i,k,m) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k,m) = f(j,i,k,m) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k,m) = f(j,i,k,m) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine mosponge4d

  subroutine spongeuv(bndu,bndv,ftenu,ftenv)
    implicit none
    type(v3dbound) , intent(in) :: bndu , bndv
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: ftenu , ftenv
    integer(ik4) :: i , j , k
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'spongeuv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_dt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_dt%ns /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bsouth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k) + &
                          (d_one-wgtd(ib))*bndu%bt(j,i,k)
            ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k) + &
                          (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%nn /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bnorth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k) + &
                           (d_one-wgtd(ib))*bndu%bt(j,i,k)
            ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k) + &
                           (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%nw /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bwest(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k) + &
                           (d_one-wgtd(ib))*bndu%bt(j,i,k)
            ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k) + &
                           (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%ne /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%beast(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            ftenu(j,i,k) = wgtd(ib)*ftenu(j,i,k) + &
                           (d_one-wgtd(ib))*bndu%bt(j,i,k)
            ftenv(j,i,k) = wgtd(ib)*ftenv(j,i,k) + &
                           (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine spongeuv

  subroutine mospongeuv(fu,fv,bndu,bndv)
    implicit none
    type(v3dbound) , intent(in) :: bndu , bndv
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: fu , fv
    integer(ik4) :: i , j , k
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mospongeuv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_ut%havebound .and. .not. ba_vt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_ut%ns /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_ut%bsouth(j,i) ) cycle
            ib = ba_ut%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtx(ib))*bndu%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_vt%ns /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_vt%bsouth(j,i) ) cycle
            ib = ba_vt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_ut%nn /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_ut%bnorth(j,i) ) cycle
            ib = ba_ut%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtx(ib))*bndu%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_vt%nn /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_vt%bnorth(j,i) ) cycle
            ib = ba_vt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_ut%nw /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_ut%bwest(j,i) ) cycle
            ib = ba_ut%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_vt%nw /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_vt%bwest(j,i) ) cycle
            ib = ba_vt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtx(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_ut%ne /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_ut%beast(j,i) ) cycle
            ib = ba_ut%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_vt%ne /= 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_vt%beast(j,i) ) cycle
            ib = ba_vt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtx(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine mospongeuv

  subroutine sponge3d(bnd,ften)
    implicit none
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: ften
    integer(ik4) :: i , j , k
    integer(ik4) :: ib , nk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    nk = size(ften,3)
    if ( ba_cr%ns /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k) = wgtx(ib)*ften(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k) = wgtx(ib)*ften(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k) = wgtx(ib)*ften(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            ften(j,i,k) = wgtx(ib)*ften(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge3d

  subroutine mosponge3d(f,bnd)
    implicit none
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: f
    integer(ik4) :: i , j , k
    integer(ik4) :: ib , nk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mosponge3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    nk = size(f,3)
    if ( ba_cr%ns /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k) = f(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k) = f(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k) = f(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do k = 1 , nk
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            f(j,i,k) = f(j,i,k) + (d_one-wgtx(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine mosponge3d

  subroutine sponge2d(bnd,ften)
    implicit none
    type(v2dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: ften
    integer(ik4) :: i , j
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_cr%ns /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          ften(j,i) = wgtx(ib)*ften(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          ften(j,i) = wgtx(ib)*ften(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          ften(j,i) = wgtx(ib)*ften(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          ften(j,i) = wgtx(ib)*ften(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge2d

  subroutine mosponge2d(f,bnd)
    implicit none
    type(v2dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: f
    integer(ik4) :: i , j
    integer(ik4) :: ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mosponge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_cr%ns /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bsouth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          f(j,i) = f(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%nn /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bnorth(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          f(j,i) = f(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%nw /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%bwest(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          f(j,i) = f(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba_cr%ne /= 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( .not. ba_cr%beast(j,i) ) cycle
          ib = ba_cr%ibnd(j,i)
          f(j,i) = f(j,i) + (d_one-wgtx(ib))*bnd%bt(j,i)
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine mosponge2d
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
    integer(ik4) , intent(in) :: ibdy , n
    real(rkx) , pointer , intent(in) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge4d3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = 1:kz )
      fg1(j,i,k) = bnd%b0(j,i,k) + xt*bnd%bt(j,i,k) - f(j,i,k,n)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge4d3d

  subroutine monudge4d3d(ibdy,f,bnd,n)
    implicit none
    integer(ik4) , intent(in) :: ibdy , n
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudge4d3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = 1:kz )
      fg1(j,i,k) = bnd%b0(j,i,k) + xt*bnd%bt(j,i,k) - f(j,i,k,n)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bsouth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bnorth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bwest(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%beast(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bsouth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bnorth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bwest(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%beast(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k,n) = f(j,i,k,n) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudge4d3d

  subroutine nudgeuv(ibdy,fu,fv,bndu,bndv,ftenu,ftenv)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: fu , fv
    type(v3dbound) , intent(in) :: bndu , bndv
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: ftenu , ftenv
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudgeuv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_dt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jde1ga:jde2ga , i = ide1ga:ide2ga , k = 1:kz )
      fg1(j,i,k) = ((bndu%b0(j,i,k) + xt*bndu%bt(j,i,k)) - fu(j,i,k))
      fg2(j,i,k) = ((bndv%b0(j,i,k) + xt*bndv%bt(j,i,k)) - fv(j,i,k))
    end do

    if ( ibdy == 1 ) then
      if ( ba_dt%ns /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%nn /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%nw /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%ne /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
    else
      if ( ba_dt%ns /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%nn /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%nw /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
      if ( ba_dt%ne /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
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
          end do
        end do
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudgeuv

  subroutine monudgeuv(ibdy,fu,fv,bndu,bndv)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: fu , fv
    type(v3dbound) , intent(in) :: bndu , bndv
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudgeuv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_ut%havebound .and. .not. ba_vt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jde1ga:jde2ga , i = ice1ga:ice2ga , k = 1:kz )
      fg1(j,i,k) = (bndu%b0(j,i,k) + xt*bndu%bt(j,i,k)) - fu(j,i,k)
    end do

    do concurrent ( j = jce1ga:jce2ga , i = ide1ga:ide2ga , k = 1:kz )
      fg2(j,i,k) = (bndv%b0(j,i,k) + xt*bndv%bt(j,i,k)) - fv(j,i,k)
    end do

    if ( ibdy == 1 ) then
      if ( ba_ut%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici1
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bsouth(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%ns /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%bsouth(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = fcd(ib)
              xg = gcd(ib)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bnorth(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%nn /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%bnorth(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = fcd(ib)
              xg = gcd(ib)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bwest(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = fcd(ib)
              xg = gcd(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%nw /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%bwest(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%beast(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = fcd(ib)
              xg = gcd(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%ne /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%beast(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    else
      if ( ba_ut%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bsouth(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%ns /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              if ( .not. ba_vt%bsouth(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = hefd(ib,k)
              xg = hegd(ib,k)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%nn /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bnorth(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%nn /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%bnorth(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = hefd(ib,k)
              xg = hegd(ib,k)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%bwest(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = hefd(ib,k)
              xg = hegd(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%nw /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%bwest(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_ut%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdi1 , jdi2
              if ( .not. ba_ut%beast(j,i) ) cycle
              ib = ba_ut%ibnd(j,i)
              xf = hefd(ib,k)
              xg = hegd(ib,k)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              fu(j,i,k) = fu(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_vt%ne /= 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_vt%beast(j,i) ) cycle
              ib = ba_vt%ibnd(j,i)
              xf = hefc(ib,k)
              xg = hegc(ib,k)
              fls0 = fg2(j,i,k)
              fls1 = fg2(j-1,i,k)
              fls2 = fg2(j+1,i,k)
              fls3 = fg2(j,i-1,k)
              fls4 = fg2(j,i+1,k)
              fv(j,i,k) = fv(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudgeuv

  subroutine monudge4d(ibdy,f,bnd,n1,n2)
    implicit none
    integer(ik4) , intent(in) :: ibdy , n1 , n2
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    integer(ik4) :: n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if
    do n = n1 , n2
      call monudge4d3d(ibdy,f,bnd,n)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudge4d

  subroutine nudge4d(ibdy,f,bnd,ften,n1,n2)
    implicit none
    integer(ik4) , intent(in) :: ibdy , n1 , n2
    real(rkx) , pointer , intent(in) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
    integer(ik4) :: n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if
    do n = n1 , n2
      call nudge4d3d(ibdy,f,bnd,ften,n)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge4d

  subroutine nudge3d(ibdy,f,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib , ns , nk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge3d'
    integer(ik4) , save :: idindx = 0
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
    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = ns:nk )
      fg1(j,i,k) = (bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)) - f(j,i,k)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
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
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge3d

  subroutine monudge3d(ibdy,f,bnd)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib , ns , nk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudge3d'
    integer(ik4) , save :: idindx = 0
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
    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = ns:nk )
      fg1(j,i,k) = (bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)) - f(j,i,k)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bsouth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bnorth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bwest(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%beast(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = fcx(ib)
              xg = gcx(ib)
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bsouth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,min(k,kz))
              xg = hegc(ib,min(k,kz))
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bnorth(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,min(k,kz))
              xg = hegc(ib,min(k,kz))
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%bwest(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,min(k,kz))
              xg = hegc(ib,min(k,kz))
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do k = ns , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( .not. ba_cr%beast(j,i) ) cycle
              ib = ba_cr%ibnd(j,i)
              xf = hefc(ib,min(k,kz))
              xg = hegc(ib,min(k,kz))
              fls0 = fg1(j,i,k)
              fls1 = fg1(j-1,i,k)
              fls2 = fg1(j+1,i,k)
              fls3 = fg1(j,i-1,k)
              fls4 = fg1(j,i+1,k)
              f(j,i,k) = f(j,i,k) + xf*fls0 -  &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudge3d

  subroutine nudge2d(ibdy,f,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(in) , dimension(:,:) :: f
    type(v2dbound) , intent(in) :: bnd
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: ften
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga )
      fg1(j,i,1) = (bnd%b0(j,i) + xt*bnd%bt(j,i)) - f(j,i)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
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
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge2d

  subroutine monudge2d(ibdy,f,bnd)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: f
    type(v2dbound) , intent(in) :: bnd
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba_cr%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga )
      fg1(j,i,1) = (bnd%b0(j,i) + xt*bnd%bt(j,i)) - f(j,i)
    end do

    if ( ibdy == 1 ) then
      if ( ba_cr%ns /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = fcx(ib)
            xg = gcx(ib)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = fcx(ib)
            xg = gcx(ib)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = fcx(ib)
            xg = gcx(ib)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = fcx(ib)
            xg = gcx(ib)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
    else
      if ( ba_cr%ns /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bsouth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%nn /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bnorth(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%nw /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%bwest(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 - &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba_cr%ne /= 0 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( .not. ba_cr%beast(j,i) ) cycle
            ib = ba_cr%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = fg1(j,i,1)
            fls1 = fg1(j-1,i,1)
            fls2 = fg1(j+1,i,1)
            fls3 = fg1(j,i-1,1)
            fls4 = fg1(j,i+1,1)
            f(j,i) = f(j,i) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine monudge2d

  subroutine couple(a,c,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: a
    real(rkx) , pointer , dimension(:,:) , intent(in) :: c
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: i , j , k
    do k = k1 , k2
      do i = i1 , i2
        do j = j1 , j2
          a(j,i,k) = a(j,i,k) * c(j,i)
        end do
      end do
    end do
  end subroutine couple

  subroutine raydampuv(z,u,v,uten,vten,ubnd,vbnd)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: uten , vten
    type(v3dbound) , intent(in) :: ubnd , vbnd
    real(rkx) :: xt , bval
    integer(ik4) :: i , j , k
    xt = xbctime + dt
    do k = 1 , min(kz,rayndamp)
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          bval = ubnd%b0(j,i,k) + xt*ubnd%bt(j,i,k)
          uten(j,i,k) = uten(j,i,k) + tau(z(j,i,k),z(j,i,1)) * (bval-u(j,i,k))
        end do
      end do
    end do
    do k = 1 , min(kz,rayndamp)
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          bval = vbnd%b0(j,i,k) + xt*vbnd%bt(j,i,k)
          vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k),z(j,i,1)) * (bval-v(j,i,k))
        end do
      end do
    end do
  end subroutine raydampuv

  subroutine raydampuv_c(z,u,v,uten,vten,sval)
  implicit none
  real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
  real(rkx) , pointer , dimension(:,:,:) , intent(in) :: u , v
  real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: uten , vten
  real(rkx) , intent(in) :: sval
  integer(ik4) :: i , j , k
  do k = 1 , min(kz,rayndamp)
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        uten(j,i,k) = uten(j,i,k) + tau(z(j,i,k),z(j,i,1)) * (sval-u(j,i,k))
      end do
    end do
  end do
  do k = 1 , min(kz,rayndamp)
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k),z(j,i,1)) * (sval-v(j,i,k))
      end do
    end do
  end do
  end subroutine raydampuv_c

  subroutine raydamp3f(z,var,vten,sval)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: var
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: vten
    real(rkx) , intent(in) :: sval
    integer(ik4) :: i , j , k
    do k = 1 , min(kzp1,rayndamp)
      do i = ici1 , ici2
        do j = jci1 , jci2
          vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k),z(j,i,1)) * (sval-var(j,i,k))
        end do
      end do
    end do
  end subroutine raydamp3f

  subroutine raydamp3(z,var,vten,bnd)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: var
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: vten
    type(v3dbound) , intent(in) :: bnd
    real(rkx) :: xt , bval
    integer(ik4) :: i , j , k
    xt = xbctime + dt
    do k = 1 , min(kz,rayndamp)
      do i = ici1 , ici2
        do j = jci1 , jci2
          bval = bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)
          vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k),z(j,i,1))*(bval-var(j,i,k))
        end do
      end do
    end do
  end subroutine raydamp3

  subroutine moraydamp(z,var,bnd,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: var
    type(v3dbound) , intent(in) :: bnd
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    real(rkx) :: xt , bval
    integer(ik4) :: i , j , k , k3
    xt = xbctime + dt
    k3 = max(k2,rayndamp)
    do k = k1 , k3
      do i = i1 , i2
        do j = j1 , j2
          bval = bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)
          var(j,i,k) = var(j,i,k) + dt*tau(z(j,i,k),z(j,i,1))*(bval-var(j,i,k))
        end do
      end do
    end do
  end subroutine moraydamp

  subroutine raydampqv(z,var,vten,bnd)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z
    real(rkx) , pointer , dimension(:,:,:,:) , intent(in) :: var
    real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: vten
    type(v3dbound) , intent(in) :: bnd
    integer(ik4) :: i , j , k
    real(rkx) :: xt , bval
    xt = xbctime + dt
    do k = 1 , min(kz,rayndamp)
      do i = ici1 , ici2
        do j = jci1 , jci2
          bval = bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)
          vten(j,i,k,iqv) = vten(j,i,k,iqv) + &
                  tau(z(j,i,k),z(j,i,1))*(bval-var(j,i,k,iqv))
        end do
      end do
    end do
  end subroutine raydampqv

  subroutine timeint2(a,b,c,j1,j2,i1,i2,rdtb)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: a , b
    real(rkx) , intent(in) :: rdtb
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: c
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: i , j
    do i = i1 , i2
      do j = j1 , j2
        c(j,i) = (a(j,i)-b(j,i))*rdtb
      end do
    end do
  end subroutine timeint2

  subroutine timeint3(a,b,c,j1,j2,i1,i2,k1,k2,rdtb)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: a , b
    real(rkx) , intent(in) :: rdtb
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: c
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: i , j , k
    do k = k1 , k2
      do i = i1 , i2
        do j = j1 , j2
          c(j,i,k) = (a(j,i,k)-b(j,i,k))*rdtb
        end do
      end do
    end do
  end subroutine timeint3

  pure real(rkx) function tau(z,zmax)
    implicit none
    real(rkx) , intent(in) :: z , zmax
    if ( z > zmax-rayhd ) then
      tau = rayalpha0 * (sin(halfpi*(d_one-(zmax-z)/rayhd)))**2
    else
      tau = d_zero
    end if
  end function tau

  subroutine paicompute(lat,ps,z,t,q,pai)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ps , lat
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: z , t , q
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: pai
    real(rkx) :: tv , tv1 , tv2 , p , zb , zdelta , zz , lrt
    integer(ik4) :: i , j , k
    ! Hydrostatic initialization of pai
    do i = ice1 , ice2
      do j = jce1 , jce2
        zdelta = z(j,i,kz)*egrav
        tv1 = t(j,i,kz) * (d_one + ep1*q(j,i,kz))
        tv2 = t(j,i,kz-1) * (d_one + ep1*q(j,i,kz-1))
        lrt = (tv2-tv1)/(z(j,i,kz-1)-z(j,i,kz))
        !lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,lat(j,i))
        lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
        tv = tv1 - 0.5_rkx*z(j,i,kz)*lrt
        zz = d_one/(rgas*tv)
        p = ps(j,i) * exp(-zdelta*zz)
        pai(j,i,kz) = (p/p00)**rovcp
      end do
    end do
    do k = kzm1 , 1 , -1
      do i = ice1 , ice2
        do j = jce1 , jce2
          tv1 = t(j,i,k) * (d_one + ep1*q(j,i,k))
          tv2 = t(j,i,k+1) * (d_one + ep1*q(j,i,k+1))
          zb = d_two*egrav*mo_dzita/(mo_atm%fmzf(j,i,k+1)*cpd) + tv1 - tv2
          zdelta = sqrt(zb**2 + d_four * tv2 * tv1)
          pai(j,i,k) = -pai(j,i,k+1) / (d_two * tv2) * (zb - zdelta)
        end do
      end do
    end do
    call exchange(pai,1,jce1,jce2,ice1,ice2,1,kz)
  end subroutine paicompute

  subroutine moloch_static_test1(xt,xq,xu,xv,xps,xts)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xps , xts
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: xt , xq , xu , xv
    integer(ik4) :: i , j , k
    xts = stdt
    xps = stdpmb
    xu = d_zero
    xv = d_zero
    xq = 1.0e-8_rkx
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xt(j,i,k) = max(xts(j,i) - lrate * mo_atm%zeta(j,i,k), 210.0_rkx)
        end do
      end do
    end do
  end subroutine moloch_static_test1

  subroutine moloch_static_test2(xt,xq,xu,xv,xps,xts)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xps , xts
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: xt , xq , xu , xv
    integer(ik4) :: i , j , k
    real(rkx) :: zlr
    xu = 10.0_rkx
    xv = d_zero
    xq = d_zero
    xts = stdt
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          zlr = -lrate
          xt(j,i,k) = max(xts(j,i) + zlr * mo_atm%zeta(j,i,k), 210.0_rkx)
        end do
      end do
    end do
    do i = ice1 , ice2
      do j = jce1 , jce2
        xps(j,i) = stdpmb * exp(-govr*mo_atm%zeta(j,i,kz)/xt(j,i,kz))
      end do
    end do
  end subroutine moloch_static_test2

  !  Computes optimal relaxation coefficients for lateral
  !  boundary conditions (Lehmann, MAP, 1993,1-14)
  !  See the paper for more comments
  !  NOTE : IS MUST BE POWER OF 2
  !  Input:  is       width of boundary relaxation zone (power of 2)
  !          gammin   minimal Courant number (c*dt/dx)
  !          gammax   maximal Courant number
  !  Output: alpha()  weight of externally specified values in the boundary
  !                   zone (corresponding to optimal relax. coefficients)
  subroutine relax (is, gammin, gammax, alpha)
    implicit none
    integer(ik4) , intent(in) :: is
    real(rkx) , intent(in) :: gammin , gammax
    real(rkx) , dimension(is) , intent(out) :: alpha
    real(rkx) , dimension(0:2*is) :: p , q , pp , qq
    real(rkx) :: my , kk , kdt2 , xxx
    integer(ik4) :: i , j , n

    n = 1
    p(0) = 0.0_rkx
    p(1) = 1.0_rkx
    q(0) = 1.0_rkx
    q(1) = 0.0_rkx
    my = sqrt(gammax/gammin)
    do
      my = sqrt((my+1.0_rkx/my)/2.0_rkx)
      do i = 0 , n+n
        pp(i) = 0.0_rkx
        qq(i) = 0.0_rkx
      end do
      do i = 0 , n
        do j = 0 , n
          pp(i+j) = pp(i+j) + p(i)*p(j) + q(i)*q(j)
          qq(i+j) = qq(i+j) + 2.0_rkx*my*p(i)*q(j)
        end do
      end do
      do i = 0 , n+n
        p(i) = pp(i)
        q(i) = qq(i)
      end do
      n = 2*n
      if ( n >= is ) exit
    end do
    do i = n , 1 , -1
      kk = p(i)/q(i-1)
      do j = i , 1 , -1
        xxx = q(j)
        q(j) = p(j) - kk*q(j-1)
        p(j) = xxx
      end do
      xxx = q(0)
      q(0) = p(0)
      p(0) = xxx
      kdt2 = kk*sqrt(gammin*gammax)
      alpha(i) = kdt2/(1.0_rkx+kdt2)
    end do
    !  Remark: this alpha corresponds to the leapfrog scheme,
    !  whereas kdt2 is independent of the integration scheme
  end subroutine relax

  subroutine invert_top_bottom(v)
    implicit none
    real(rkx) , dimension(:) , intent(inout) :: v
    real(rkx), dimension(size(v)) :: swap
    integer :: nk , k , kk
    swap = v
    nk = size(v)
    do k = 1 , nk
      kk = nk-k+1
      v(k) = swap(kk)
    end do
  end subroutine invert_top_bottom

end module mod_bdycod

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
