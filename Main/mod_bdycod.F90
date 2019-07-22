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
  use mod_slabocean

  implicit none

  private

  public :: allocate_mod_bdycon , init_bdy , bdyin , bdyval
  public :: sponge , nudge , setup_bdycon , raydamp

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
  end interface raydamp

  logical , parameter :: bdyflow = .true.

  contains

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
    if ( iboudy == 5 ) then
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
    real(rkx) :: xfun
    integer(ik4) :: n , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'setup_bdycon'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    rdtbdy = d_one / dtbdys
    if ( iboudy == 1 .or. iboudy == 5 ) then
      if ( bdy_nm > d_zero ) then
        fnudge = bdy_nm
      else
        fnudge = 0.1_rkx/(dtsec*2.0_rkx)
      end if
      if ( bdy_dm > d_zero ) then
        gnudge = bdy_dm
      else
        gnudge = d_one/(50.0_rkx*dtsec)
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
      do k = 1 , kz
        if ( hsigma(k) < 0.4_rkx ) then
          anudge(k) = high_nudge
        else if ( hsigma(k) < 0.8_rkx ) then
          anudge(k) = medium_nudge
        else
          anudge(k) = low_nudge
        end if
      end do
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
                     xtb%b0,xqb%b0,xppb%b0,xwwb%b0)
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
    else
      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,xub%b0,xvb%b0,xtb%b0,xqb%b0)
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
    end if
    call exchange(xub%b0,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b0,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xtb%b0,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b0,1,jce1,jce2,ice1,ice2,1,kz)
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
                     xtb%b1,xqb%b1,xppb%b1,xwwb%b1)
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
    else
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1,xtb%b1,xqb%b1)
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
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
    end if
    call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xtb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b1,1,jce1,jce2,ice1,ice2,1,kz)
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
    call timeint(xub%b1,xub%b0,xub%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
    call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2)
    if ( idynamic == 1 ) then
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga)
    else if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1)
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
    xtsb%b0(:,:) = xtsb%b1(:,:)
    if ( idynamic == 2 ) then
      xppb%b0(:,:,:) = xppb%b1(:,:,:)
      xwwb%b0(:,:,:) = xwwb%b1(:,:,:)
    else
      xpsb%b0(:,:) = xpsb%b1(:,:)
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
                     xtb%b1,xqb%b1,xppb%b1,xwwb%b1)
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
    else
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1,xtb%b1,xqb%b1)
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
    end if
    !
    ! Couple pressure u,v,t,q
    !
    if ( idynamic /= 3 ) then
      call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
      call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
    end if
    call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xtb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%b1,1,jce1,jce2,ice1,ice2,1,kz)
    if ( idynamic == 2 ) then
      call couple(xppb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
      call couple(xwwb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%b1,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%b1,1,jce1,jce2,ice1,ice2,1,kzp1)
    else
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1ga,jce2ga,ice1ga,ice2ga)
    end if

    ! Linear time interpolation
    call timeint(xub%b1,xub%b0,xub%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,kz)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
    if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kz)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1)
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

    call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2)

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
      if ( idynamic == 1 ) then
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
            do i = ici1 , ici2
              mo_atm%u(jde1,i,k) = xub%b0(jde1,i,k)
              mo_atm%v(jce1,i,k) = xub%b0(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%u(jde2,i,k) = xub%b0(jde2,i,k)
              mo_atm%u(jce2,i,k) = xub%b0(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%u(j,ice1,k) = xvb%b0(j,ice1,k)
              mo_atm%v(j,ide1,k) = xvb%b0(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%v(j,ice2,k) = xvb%b0(j,ice2,k)
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
      if ( idynamic == 1 ) then
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
            do i = ici1 , ici2
              mo_atm%u(jde1,i,k) = xub%b0(jde1,i,k) + xt*xub%bt(jde1,i,k)
              mo_atm%v(jce1,i,k) = xvb%b0(jce1,i,k) + xt*xvb%bt(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%u(jde2,i,k) = xub%b0(jde2,i,k) + xt*xub%bt(jde2,i,k)
              mo_atm%v(jce2,i,k) = xvb%b0(jce2,i,k) + xt*xvb%bt(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%u(j,ice1,k) = xub%b0(j,ice1,k) + xt*xub%bt(j,ice1,k)
              mo_atm%v(j,ide1,k) = xvb%b0(j,ide1,k) + xt*xvb%bt(j,ide1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%u(j,ice2,k) = xub%b0(j,ice2,k) + xt*xub%bt(j,ice2,k)
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
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce2,i,k) = xtb%b0(jce2,i,k)
              mo_atm%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice1,k) = xtb%b0(j,ice1,k)
              mo_atm%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice2,k) = xtb%b0(j,ice2,k)
              mo_atm%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k)
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce1,i,k) = xtb%b0(jce1,i,k)
              atm1%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k)
            end do
          end do
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
              mo_atm%t(jce1,i,k)      = xtb%b0(jce1,i,k) + xt*xtb%bt(jce1,i,k)
              mo_atm%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k) + xt*xqb%bt(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%t(jce2,i,k)      = xtb%b0(jce2,i,k) + xt*xtb%bt(jce2,i,k)
              mo_atm%qx(jce2,i,k,iqv) = xqb%b0(jce2,i,k) + xt*xqb%bt(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice1,k)      = xtb%b0(j,ice1,k) + xt*xtb%bt(j,ice1,k)
              mo_atm%qx(j,ice1,k,iqv) = xqb%b0(j,ice1,k) + xt*xqb%bt(j,ice1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              mo_atm%t(j,ice2,k)      = xtb%b0(j,ice2,k) + xt*xtb%bt(j,ice2,k)
              mo_atm%qx(j,ice2,k,iqv) = xqb%b0(j,ice2,k) + xt*xqb%bt(j,ice2,k)
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              atm1%t(jce1,i,k)      = xtb%b0(jce1,i,k) + xt*xtb%bt(jce1,i,k)
              atm1%qx(jce1,i,k,iqv) = xqb%b0(jce1,i,k) + xt*xqb%bt(jce1,i,k)
            end do
          end do
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
    !
    ! set boundary values for p*qx
    ! *** note ***
    ! for large domain, we assume the boundary tendencies are not available.
    !
    ! if the boundary values and tendencies are not available,
    ! determine boundary values depends on inflow/outflow:
    ! inflow  : set it equal to zero.
    ! outflow : get from interior point.
    !
    ! west boundary:
    !
    if ( idynamic == 3 ) then
      if ( bdyflow ) then
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            do k = 1 , kz
              do i = ice1 , ice2
                qxint = mo_atm%qx(jci1,i,k,n)
                windavg = mo_atm%u(jde1,i,k) + mo_atm%u(jdi1,i,k)
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
            do k = 1 , kz
              do i = ice1 , ice2
                qxint = mo_atm%qx(jci2,i,k,n)
                windavg = mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k)
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
            do k = 1 , kz
              do j = jci1 , jci2
                qxint = mo_atm%qx(j,ici1,k,n)
                windavg = mo_atm%v(j,ide1,k) + mo_atm%v(j,idi1,k)
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
              do j = jci1 , jci2
                qxint = mo_atm%qx(j,ici2,k,n)
                windavg = mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k)
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
            do k = 1 , kz
              do i = ice1 , ice2
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
            do k = 1 , kz
              do i = ice1 , ice2
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
            do k = 1 , kz
              do j = jci1 , jci2
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
            do k = 1 , kz
              do j = jci1 , jci2
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
            do k = 1 , kz
              do i = ice1 , ice2
                qxint = atm1%qx(jci1,i,k,n)/sfs%psa(jci1,i)
                windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
                if ( windavg > d_zero ) then
                  atm1%qx(jce1,i,k,n) = d_zero
                else
                  atm1%qx(jce1,i,k,n) = qxint*sfs%psa(jce1,i)
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
            do k = 1 , kz
              do i = ice1 , ice2
                qxint = atm1%qx(jci2,i,k,n)/sfs%psa(jci2,i)
                windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
                if ( windavg < d_zero ) then
                  atm1%qx(jce2,i,k,n) = d_zero
                else
                  atm1%qx(jce2,i,k,n) = qxint*sfs%psa(jce2,i)
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
            do k = 1 , kz
              do j = jci1 , jci2
                qxint = atm1%qx(j,ici1,k,n)/sfs%psa(j,ici1)
                windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
                if ( windavg > d_zero ) then
                  atm1%qx(j,ice1,k,n) = d_zero
                else
                  atm1%qx(j,ice1,k,n) = qxint*sfs%psa(j,ice1)
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
              do j = jci1 , jci2
                qxint = atm1%qx(j,ici2,k,n)/sfs%psa(j,ici2)
                windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
                if ( windavg < d_zero ) then
                  atm1%qx(j,ice2,k,n) = d_zero
                else
                  atm1%qx(j,ice2,k,n) = qxint*sfs%psa(j,ice2)
                end if
              end do
            end do
          end do
        end if
      else
        if ( ma%has_bdyleft ) then
          do n = iqfrst , iqlst
            do k = 1 , kz
              do i = ice1 , ice2
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
            do k = 1 , kz
              do i = ice1 , ice2
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
            do k = 1 , kz
              do j = jci1 , jci2
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
            do k = 1 , kz
              do j = jci1 , jci2
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

    if ( .not. ba_dt%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba_dt%ns /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bsouth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_dt%bsouth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%nn /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bnorth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_dt%bnorth(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%nw /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%bwest(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_dt%bwest(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba_dt%ne /= 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            if ( .not. ba_dt%beast(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fu(j,i,k) = fu(j,i,k) + (d_one-wgtd(ib))*bndu%bt(j,i,k)
          end do
        end do
        do i = idi1 , idi2
          do j = jci1 , jci2
            if ( .not. ba_dt%beast(j,i) ) cycle
            ib = ba_dt%ibnd(j,i)
            fv(j,i,k) = fv(j,i,k) + (d_one-wgtd(ib))*bndv%bt(j,i,k)
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
    real(rkx) , parameter :: nfac = 1.0e3_rkx
    real(rkx) , parameter :: rfac = d_one/nfac
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
      fg1(j,i,k) = nfac*(bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)) - nfac*f(j,i,k,n)
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              ften(j,i,k,n) = ften(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
    real(rkx) , pointer , intent(in) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    real(rkx) , parameter :: nfac = 1.0e3_rkx
    real(rkx) , parameter :: rfac = d_one/nfac
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
      fg1(j,i,k) = nfac*(bnd%b0(j,i,k) + xt*bnd%bt(j,i,k)) - nfac*f(j,i,k,n)
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
              f(j,i,k,n) = f(j,i,k,n) + rfac * (xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0))
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
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: fu , fv
    type(v3dbound) , intent(in) :: bndu , bndv
    real(rkx) :: xt , xf , xg , fls0 , fls1 , fls2 , fls3 , fls4
    integer(ik4) :: i , j , k , ib
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'monudgeuv'
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

    do concurrent ( j = jde1ga:jde2ga , i = ice1:ice2 , k = 1:kz )
      fg1(j,i,k) = ((bndu%b0(j,i,k) + xt*bndu%bt(j,i,k)) - fu(j,i,k))
    end do

    do concurrent ( j = jce1:jce2 , i = ide1ga:ide2ga , k = 1:kz )
      fg2(j,i,k) = ((bndv%b0(j,i,k) + xt*bndv%bt(j,i,k)) - fv(j,i,k))
    end do

    if ( ibdy == 1 ) then
      if ( ba_dt%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici1
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%bsouth(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
      if ( ba_dt%nn /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%bnorth(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
      if ( ba_dt%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%bwest(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
      if ( ba_dt%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%beast(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
              xf = fcd(ib)
              xg = gcd(ib)
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
      if ( ba_dt%ns /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              if ( .not. ba_dt%bsouth(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%bnorth(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
      if ( ba_dt%nw /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 - &
                            xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%bwest(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
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
      if ( ba_dt%ne /= 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
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
              fu(j,i,k) = fu(j,i,k) + xf*fls0 -  &
                          xg*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
          do i = idi1 , idi2
            do j = jci1 , jci2
              if ( .not. ba_dt%beast(j,i) ) cycle
              ib = ba_dt%ibnd(j,i)
              xf = hefd(ib,k)
              xg = hegd(ib,k)
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
    real(rkx) , pointer , intent(in) , dimension(:,:,:,:) :: f
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
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: f
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
      fg1(j,i,1) = ((bnd%b0(j,i) + xt*bnd%bt(j,i)) - f(j,i))
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
    real(rkx) , pointer , intent(in) , dimension(:,:) :: f
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
      fg1(j,i,1) = ((bnd%b0(j,i) + xt*bnd%bt(j,i)) - f(j,i))
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
    real(rkx) :: zz , xt , bval
    integer(ik4) :: i , j , k
    xt = xbctime + dt
    do k = 1 , min(kz,rayndamp)
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          bval = ubnd%b0(j,i,k) + xt*ubnd%bt(j,i,k)
          zz = d_rfour * (z(j,i,k) + z(j-1,i,k) + z(j,i-1,k) + z(j-1,i-1,k))
          uten(j,i,k) = uten(j,i,k) + tau(zz) * (bval-u(j,i,k))
        end do
      end do
    end do
    do k = 1 , min(kz,rayndamp)
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          bval = vbnd%b0(j,i,k) + xt*vbnd%bt(j,i,k)
          zz = d_rfour * (z(j,i,k) + z(j-1,i,k) + z(j,i-1,k) + z(j-1,i-1,k))
          vten(j,i,k) = vten(j,i,k) + tau(zz) * (bval-v(j,i,k))
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
  real(rkx) :: zz
  integer(ik4) :: i , j , k
  do k = 1 , min(kz,rayndamp)
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        zz = d_rfour * (z(j,i,k) + z(j-1,i,k) + z(j,i-1,k) + z(j-1,i-1,k))
        uten(j,i,k) = uten(j,i,k) + tau(zz) * (sval-u(j,i,k))
      end do
    end do
  end do
  do k = 1 , min(kz,rayndamp)
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        zz = d_rfour * (z(j,i,k) + z(j-1,i,k) + z(j,i-1,k) + z(j-1,i-1,k))
        vten(j,i,k) = vten(j,i,k) + tau(zz) * (sval-v(j,i,k))
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
          vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k)) * (sval-var(j,i,k))
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
          vten(j,i,k) = vten(j,i,k) + tau(z(j,i,k))*(bval-var(j,i,k))
        end do
      end do
    end do
  end subroutine raydamp3

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
                  tau(z(j,i,k))*(bval-var(j,i,k,iqv))
        end do
      end do
    end do
  end subroutine raydampqv

  subroutine timeint2(a,b,c,j1,j2,i1,i2)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: a , b
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: c
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: i , j
    do i = i1 , i2
      do j = j1 , j2
        c(j,i) = (a(j,i)-b(j,i))*rdtbdy
      end do
    end do
  end subroutine timeint2

  subroutine timeint3(a,b,c,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: a , b
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: c
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: i , j , k
    do k = k1 , k2
      do i = i1 , i2
        do j = j1 , j2
          c(j,i,k) = (a(j,i,k)-b(j,i,k))*rdtbdy
        end do
      end do
    end do
  end subroutine timeint3

  pure real(rkx) function tau(z)
    implicit none
    real(rkx) , intent(in) :: z
    if ( z > rayzd-rayhd ) then
      tau = rayalpha0 * (sin(halfpi*(d_one-(rayzd-z)/rayhd)))**2
    else
      tau = d_zero
    end if
  end function tau

end module mod_bdycod

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
