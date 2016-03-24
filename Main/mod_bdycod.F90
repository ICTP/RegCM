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
  real(rk8) , pointer , dimension(:,:) :: sue , sui , nue , nui , &
                                         sve , svi , nve , nvi
  real(rk8) , pointer , dimension(:,:) :: wue , wui , eue , eui , &
                                         wve , wvi , eve , evi
  real(rk8) , pointer , dimension(:,:) :: psdot
  real(rk8) , pointer , dimension(:) :: fcx , gcx
  real(rk8) , pointer , dimension(:) :: fcd , gcd
  real(rk8) , pointer , dimension(:) :: lfc , lgc
  real(rk8) , pointer , dimension(:,:) :: hefc , hegc
  real(rk8) , pointer , dimension(:,:) :: fefc , fegc
  real(rk8) , pointer , dimension(:,:) :: efc , egc
  real(rk8) , pointer , dimension(:) :: wgtd
  real(rk8) , pointer , dimension(:) :: wgtx
  real(rk8) :: fnudge , gnudge , rdtbdy
  integer(ik4) :: nbdm
  integer(ik4) :: som_month

  interface timeint
    module procedure timeint2 , timeint3
  end interface timeint

  interface nudge
    module procedure nudge4d , nudge3d , nudge2d
  end interface nudge

  interface sponge
    module procedure sponge4d , sponge3d , sponge2d
  end interface sponge

  public :: sponge , nudge , setup_bdycon

  contains

  subroutine allocate_mod_bdycon
    implicit none
    if ( iboudy == 1 ) then
      call getmem1d(fcx,1,nspgx,'bdycon:fcx')
      call getmem1d(gcx,1,nspgx,'bdycon:gcx')
      call getmem1d(fcd,1,nspgd,'bdycon:fcd')
      call getmem1d(gcd,1,nspgd,'bdycon:gcd')
    else if ( iboudy == 4 ) then
      call getmem1d(wgtd,1,nspgd,'bdycon:wgtd')
      call getmem1d(wgtx,1,nspgx,'bdycon:wgtx')
    else if ( iboudy == 5 ) then
      nbdm = max(nspgx,nspgd)
      call getmem2d(fefc,1,nbdm,1,kzp1,'bdycon:fefc')
      call getmem2d(fegc,1,nbdm,1,kzp1,'bdycon:fegc')
      call getmem2d(hefc,1,nbdm,1,kz,'bdycon:hefc')
      call getmem2d(hegc,1,nbdm,1,kz,'bdycon:hegc')
    end if

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
    if ( idynamic == 2 ) then
    end if
  end subroutine allocate_mod_bdycon

  subroutine setup_bdycon
    implicit none
    integer(ik4) :: i , j , n , k
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
      fnudge = 0.1D0/dt2
      gnudge = (dxsq/dt)/50.0D0
    end if
    if ( iboudy == 1 ) then
      do n = 2 , nspgx-1
        fcx(n) = fnudge*xfun(n,.false.)
        gcx(n) = gnudge*xfun(n,.false.)
      end do
      do n = 2 , nspgd-1
        fcd(n) = fnudge*xfun(n,.true.)
        gcd(n) = gnudge*xfun(n,.true.)
      end do
    else if ( iboudy == 4 ) then
      wgtd(1) = 0.00D0
      wgtd(2) = 0.20D0
      wgtd(3) = 0.55D0
      wgtd(4) = 0.80D0
      wgtd(5) = 0.95D0
      do k = 4 , nspgx
        wgtd(k) = d_one
      end do
      wgtx(1) = 0.0D0
      wgtx(2) = 0.4D0
      wgtx(3) = 0.7D0
      wgtx(4) = 0.9D0
      do k = 5 , nspgx
        wgtx(k) = 1.0D0
      end do
    else if ( iboudy == 5 ) then
      do k = 1 , kzp1
        if ( sigma(k) < 0.4D0 ) then
          anudgf(k) = high_nudge
        else if ( sigma(k) < 0.8D0 ) then
          anudgf(k) = medium_nudge
        else
          anudgf(k) = low_nudge
        end if
        do n = 2 , nbdm-1
          fefc(n,k) = fnudge*xfune(n,k,anudgf)
          fegc(n,k) = gnudge*xfune(n,k,anudgf)
        end do
      end do
      do k = 1 , kz
        if ( hsigma(k) < 0.4D0 ) then
          anudgh(k) = high_nudge
        else if ( hsigma(k) < 0.8D0 ) then
          anudgh(k) = medium_nudge
        else
          anudgh(k) = low_nudge
        end if
        do n = 2 , nbdm-1
          hefc(n,k) = fnudge*xfune(n,k,anudgh)
          hegc(n,k) = gnudge*xfune(n,k,anudgh)
        end do
      end do
    end if
    if ( idynamic == 2 .and. ifupr == 1 ) then
      wtbdy(jci1:jci2,ici1:ici2) = d_one
      if ( iboudy == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( .not. ba_cr%bwest(j,i)  .and. .not. ba_cr%beast(j,i) .and. &
                 .not. ba_cr%bsouth(j,i) .and. .not. ba_cr%bnorth(j,i) ) cycle
            wtbdy(j,i) = d_one-xfun(ba_cr%ibnd(j,i),.false.)
          end do
        end do
      else if ( iboudy == 4 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( .not. ba_cr%bwest(j,i)  .and. .not. ba_cr%beast(j,i) .and. &
                 .not. ba_cr%bsouth(j,i) .and. .not. ba_cr%bnorth(j,i) ) cycle
            wtbdy(j,i) = wgtx(ba_cr%ibnd(j,i))
          end do
        end do
      else if ( iboudy == 5 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( .not. ba_cr%bwest(j,i)  .and. .not. ba_cr%beast(j,i) .and. &
                 .not. ba_cr%bsouth(j,i) .and. .not. ba_cr%bnorth(j,i) ) cycle
            wtbdy(j,i) = d_one-xfune(ba_cr%ibnd(j,i),1,anudgf)
          end do
        end do
      end if
      call grid_fill(wtbdy,wtbdy_g,jce1,jce2,ice1,ice2)
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains

      pure real(rk8) function xfun(mm,ldot)
        implicit none
        integer(ik4) , intent(in) :: mm
        logical , intent(in) :: ldot
        if ( ldot ) then
          xfun = dble(nspgd-mm)/dble(nspgd-2)
        else
          xfun = dble(nspgx-mm)/dble(nspgx-2)
        end if
      end function xfun

      pure real(rk8) function xfune(mm,kk,an)
        implicit none
        integer(ik4) , intent(in) :: mm , kk
        real(rk8) , dimension(:) , intent(in) :: an
        xfune = dexp(-dble(mm-2)/an(kk))
      end function xfune

  end subroutine setup_bdycon

  subroutine init_bdy
    implicit none
    integer(ik4) :: datefound
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
    nbdytime = 0
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
      call read_icbc(hyps0,ts0,xub%b0,xvb%b0,xtb%b0,xqb%b0,xppb%b0,xwwb%b0)
      hyps0 = hyps0 * d_r10 - ptop
    else
      call read_icbc(xpsb%b0,ts0,xub%b0,xvb%b0,xtb%b0,xqb%b0)
    end if

    if ( islab_ocean == 1 .and. do_qflux_adj ) then
      som_month = xmonth
      datefound = som_search(som_month)
      if (datefound < 0) then
        appdat = tochar(bdydate1)
        call fatal(__FILE__,__LINE__,'SOM for '//appdat//' not found')
      end if
      call read_som(qflb0)
      where ( mddom%ldmsk > 0 ) qflb0 = d_zero
      tdif = bdydate1-monfirst(bdydate1)
      xslabtime = tohours(tdif)*secph
    end if

    if ( myid == italk ) then
      appdat = tochar(bdydate1)
      if ( .not. ifrest ) then
        write(stdout,*) 'READY IC DATA for ', appdat
      else
        write(stdout,*) 'READY BC DATA for ', appdat
      end if
    end if

    bdydate2 = bdydate2 + intbdy
    if ( myid == italk ) then
      write(stdout,'(a,i10,a,i8)') ' SEARCH BC data for ', toint10(bdydate2), &
                      ', ktau = ', ktau
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
      call read_icbc(hyps1,ts1,xub%b1,xvb%b1,xtb%b1,xqb%b1,xppb%b1,xwwb%b1)
      hyps1 = hyps1 * d_r10 - ptop
    else
      call read_icbc(xpsb%b1,ts1,xub%b1,xvb%b1,xtb%b1,xqb%b1)
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
      where ( mddom%ldmsk > 0 ) qflb1 = d_zero
      tdif = bdydate2-prevmon(bdydate2)
      qflbt = (qflb1-qflb0)/(tohours(tdif)*secph)
    end if

    if ( myid == italk ) then
      write (stdout,*) 'READY  BC from     ' , &
            toint10(bdydate1) , ' to ' , toint10(bdydate2)
    end if

    bdydate1 = bdydate2

    if ( idynamic == 2 ) then
      xpsb%b0(:,:) = atm0%ps(:,:) * d_r1000 ! Cb
      psdot(:,:) = atm0%psdot(jde1:jde2,ide1:ide2) * d_r1000
      xpsb%b1(:,:) = xpsb%b0(:,:)
    else
      xpsb%b0(:,:) = (xpsb%b0(:,:)*d_r10)-ptop
      xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
      call exchange(xpsb%b0,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b0,psdot)
    end if
    !
    ! Calculate P* on dot points
    !
    !
    ! Couple pressure u,v,t,q (pp,ww)
    !
    call couple(xub%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xvb%b0,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xtb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
    call couple(xqb%b0,xpsb%b0,jce1,jce2,ice1,ice2,1,kz)
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
    !
    ! Repeat for T2
    !
    if ( idynamic == 1 ) then
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
    call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
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
    !
    ! Calculate time varying component
    !
    call timeint(xub%b1,xub%b0,xub%bt,jde1,jde2,ide1,ide2,1,kz)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1,jde2,ide1,ide2,1,kz)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1,jce2,ice1,ice2,1,kz)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xtb%bt,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%bt,1,jce1,jce2,ice1,ice2,1,kz)
    if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%bt,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%bt,1,jce1,jce2,ice1,ice2,1,kzp1)
    else
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1,jce2,ice1,ice2)
      call exchange(xpsb%bt,1,jce1,jce2,ice1,ice2)
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
    integer(ik4) :: i , j , n , datefound
    character(len=32) :: appdat
    logical :: update_slabocn
    real(rk8) :: sfice_temp
    type (rcm_time_interval) :: tdif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bdyin'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    update_slabocn = ( islab_ocean == 1 .and. &
      do_qflux_adj .and. som_month /= xmonth )

    nbdytime = 0
    xbctime = d_zero

    xub%b0(:,:,:) = xub%b1(:,:,:)
    xvb%b0(:,:,:) = xvb%b1(:,:,:)
    xtb%b0(:,:,:) = xtb%b1(:,:,:)
    xqb%b0(:,:,:) = xqb%b1(:,:,:)
    ts0(:,:) = ts1(:,:)
    if ( idynamic == 2 ) then
      xppb%b0(:,:,:) = xppb%b1(:,:,:)
      xwwb%b0(:,:,:) = xwwb%b1(:,:,:)
    else
      xpsb%b0(:,:) = xpsb%b1(:,:)
    end if

    ! Data are monthly
    if ( update_slabocn ) then
      som_month = xmonth
      qflb0 = qflb1
      tdif = bdydate1-monfirst(bdydate1)
      xslabtime = tohours(tdif)*secph
    end if

    bdydate2 = bdydate2 + intbdy
    if ( myid == italk ) then
      write(stdout,'(a,i10,a,i8)') ' SEARCH BC data for ', toint10(bdydate2), &
                      ', ktau = ', ktau
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
      call read_icbc(hyps1,ts1,xub%b1,xvb%b1,xtb%b1,xqb%b1,xppb%b1,xwwb%b1)
      hyps1 = hyps1 * d_r10 - ptop
    else
      call read_icbc(xpsb%b1,ts1,xub%b1,xvb%b1,xtb%b1,xqb%b1)
    end if

    if ( update_slabocn ) then
      datefound = som_search(som_month)
      if ( datefound < 0 ) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'SOM for '//appdat//' not found')
      end if
      call read_som(qflb1)
      where ( mddom%ldmsk > 0 ) qflb1 = d_zero
      tdif = bdydate2-prevmon(bdydate2)
      qflbt = (qflb1-qflb0)/(tohours(tdif)*secph)
    end if
    !
    ! Convert surface pressure to pstar
    !
    if ( idynamic == 1 ) then
      xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
      call exchange(xpsb%b1,1,jce1,jce2,ice1,ice2)
      call psc2psd(xpsb%b1,psdot)
    else
      xpsb%b1(:,:) = xpsb%b0(:,:)
    end if
    !
    ! Couple pressure u,v,t,q
    !
    call couple(xub%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xvb%b1,psdot,jde1,jde2,ide1,ide2,1,kz)
    call couple(xtb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
    call couple(xqb%b1,xpsb%b1,jce1,jce2,ice1,ice2,1,kz)
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
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1,jce2,ice1,ice2)
      call exchange(xpsb%bt,1,jce1,jce2,ice1,ice2)
    end if
    call timeint(xub%b1,xub%b0,xub%bt,jde1,jde2,ide1,ide2,1,kz)
    call timeint(xvb%b1,xvb%b0,xvb%bt,jde1,jde2,ide1,ide2,1,kz)
    call timeint(xtb%b1,xtb%b0,xtb%bt,jce1,jce2,ice1,ice2,1,kz)
    call timeint(xqb%b1,xqb%b0,xqb%bt,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange(xtb%bt,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(xqb%bt,1,jce1,jce2,ice1,ice2,1,kz)
    if ( idynamic == 2 ) then
      call timeint(xppb%b1,xppb%b0,xppb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xwwb%b1,xwwb%b0,xwwb%bt,jce1,jce2,ice1,ice2,1,kzp1)
      call exchange(xppb%bt,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(xwwb%bt,1,jce1,jce2,ice1,ice2,1,kzp1)
    end if
    !
    ! Update ground temperature on Ocean/Lakes
    !
    if ( islab_ocean == 0 ) then
      sfice_temp = icetemp
      if ( idcsst == 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nnsg
              lms%sst(n,j,i) = ts1(j,i)
            end do
          end do
        end do
      end if
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Update temperatures over water
          if ( mddom%ldmsk(j,i) == 0 ) then
            if ( iocncpl == 1 .or. iwavcpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            sfs%tga(j,i) = ts1(j,i)
            sfs%tgb(j,i) = ts1(j,i)
          end if
          ! Sea ice correction
          if ( iseaice == 1 ) then
            if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
            if ( iocncpl == 1 ) then
              if ( cplmsk(j,i) /= 0 ) cycle
            end if
            if ( ts1(j,i) <= icetemp .and. mddom%ldmsk(j,i) == 0 ) then
              sfs%tga(j,i) = sfice_temp
              sfs%tgb(j,i) = sfice_temp
              ts1(j,i) = icetemp
              mddom%ldmsk(j,i) = 2
              do n = 1 , nnsg
                if ( mdsub%ldmsk(n,j,i) == 0 ) then
                  mdsub%ldmsk(n,j,i) = 2
                  lms%sfice(n,j,i) = 0.50D0 ! 10 cm
                end if
              end do
            else if ( ts1(j,i) > icetemp .and. mddom%ldmsk(j,i) == 2 ) then
              ! Decrease the surface ice to melt it
              sfs%tga(j,i) = ts1(j,i)
              sfs%tgb(j,i) = ts1(j,i)
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

    if ( myid == italk ) then
      write (stdout,*) 'READY  BC from     ' , &
            toint10(bdydate1) , ' to ' , toint10(bdydate2)
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
    real(rk8) , intent(in) :: xt
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
      call exchange_bdy_tb(wue,1,kz)
      call exchange_bdy_tb(wui,1,kz)
      call exchange_bdy_tb(wve,1,kz)
      call exchange_bdy_tb(wvi,1,kz)
    end if

    if ( ma%has_bdyright ) then
      call exchange_bdy_tb(eue,1,kz)
      call exchange_bdy_tb(eui,1,kz)
      call exchange_bdy_tb(eve,1,kz)
      call exchange_bdy_tb(evi,1,kz)
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
    real(rk8) :: qxint , tkeint , qext , qint
    integer(ik4) :: i , j , k , n
    real(rk8) :: windavg
    real(rk8) :: xt
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
    if ( ktau > 1 ) then
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

    call bdyuv(xt)

    !
    ! Set boundary values for p*t:
    ! Set boundary values for p*qv:
    !
    if ( iboudy == 0 ) then
      !
      ! fixed boundary conditions:
      !
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
    else
      !
      ! time-dependent boundary conditions:
      !
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
        end if
      end if
      !
      ! Top BC for QV
      !
      !do j = jci1 , jci2
      !  do i = ici1 , ici2
      !    atm1%qx(j,i,1,iqv) = (0.10D0 * atm1%qx(j,i,1,iqv) + &
      !              0.90D0 * (xqb%b0(j,i,1) + xt*xqb%bt(j,i,1)))
      !    atm1%qx(j,i,1,iqv) = max(atm1%qx(j,i,1,iqv), &
      !                xqb%b0(j,i,1) + xt*xqb%bt(j,i,1))
      !    atm1%qx(j,i,1,iqv) = xqb%b0(j,i,1) + xt*xqb%bt(j,i,1)
      !  end do
      !end do
      !
    end if

    if ( iboudy == 3 .or. iboudy == 4 ) then
      !
      ! determine QV boundary values depends on inflow/outflow:
      !
      ! west boundary:
      !
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            qext = atm1%qx(jce1,i,k,iqv)/sfs%psa(jce1,i)
            qint = atm1%qx(jci1,i,k,iqv)/sfs%psa(jci1,i)
            windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
            if ( windavg >= d_zero ) then
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
            if ( windavg <= d_zero ) then
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
            if ( windavg >= d_zero ) then
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
            if ( windavg <= d_zero ) then
              atm1%qx(j,ice2,k,iqv) = qext*sfs%psa(j,ice2)
            else
              atm1%qx(j,ice2,k,iqv) = qint*sfs%psa(j,ice2)
            end if
          end do
        end do
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
    if ( ma%has_bdyleft ) then
      do n = iqfrst , iqlst
        do k = 1 , kz
          do i = ice1 , ice2
            qxint = atm1%qx(jci1,i,k,n)/sfs%psa(jci1,i)
            windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
            if ( windavg >= d_zero ) then
              atm1%qx(jce1,i,k,n) = d_half*(minqx+qxint*sfs%psa(jce1,i))
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
            if ( windavg <= d_zero ) then
              atm1%qx(jce2,i,k,n) = d_half*(minqx+qxint*sfs%psa(jce2,i))
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
            if ( windavg >= d_zero ) then
              atm1%qx(j,ice1,k,n) = d_half*(minqx+qxint*sfs%psa(j,ice1))
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
            if ( windavg <= d_zero ) then
              atm1%qx(j,ice2,k,n) = d_half*(minqx+qxint*sfs%psa(j,ice2))
            else
              atm1%qx(j,ice2,k,n) = qxint*sfs%psa(j,ice2)
            end if
          end do
        end do
      end do
    end if

    if ( ibltyp == 2 ) then
      if ( ktau == 0 ) then
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
        if ( ma%has_bdyleft ) then
          atm1%tke(jce1,:,1) = tkemin ! East boundary
          atm2%tke(jce1,:,1) = tkemin ! East boundary
          do k = 1 , kz
            do i = ice1 , ice2
              tkeint = atm1%tke(jci1,i,k+1)
              windavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
              if ( windavg >= d_zero ) then
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
          atm1%tke(jce2,:,1) = tkemin ! West boundary
          atm2%tke(jce2,:,1) = tkemin ! West boundary
          do k = 1 , kz
            do i = ice1 , ice2
              tkeint = atm1%tke(jci2,i,k+1)
              windavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
              if ( windavg <= d_zero ) then
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
          do k = 1 , kz
            do j = jci1 , jci2
              tkeint = atm1%tke(j,ici1,k+1)
              windavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
              if ( windavg >= d_zero ) then
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
          atm1%tke(:,ice2,1) = tkemin  ! South boundary
          atm2%tke(:,ice2,1) = tkemin  ! South boundary
          do k = 1 , kz
            do j = jci1 , jci2
              tkeint = atm1%tke(j,ici2,k+1)
              windavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
              if ( windavg < d_zero ) then
                atm1%tke(j,ice2,k+1) = tkemin
              else
                atm1%tke(j,ice2,k+1) = tkeint
              end if
            end do
          end do
        end if
      end if
    end if

    if ( ichem == 1 ) then
      call chem_bdyval
    end if
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
  subroutine sponge4d(nk,m,ba,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: nk , m
    type(bound_area) , intent(in) :: ba
    type(v3dbound) , intent(in) :: bnd
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: ften

    integer(ik4) :: i , j , k
    integer(ik4) :: ib , i1 , i2 , j1 , j2
    real(rk8) , pointer , dimension(:) :: wg
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba%dotflag ) then
      wg => wgtd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      wg => wgtx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if

    if ( ba%ns /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bsouth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k,m) = wg(ib)*ften(j,i,k,m) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%nn /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bnorth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k,m) = wg(ib)*ften(j,i,k,m) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%nw /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bwest(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k,m) = wg(ib)*ften(j,i,k,m) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%ne /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%beast(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k,m) = wg(ib)*ften(j,i,k,m) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge4d

  subroutine sponge3d(nk,ba,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: nk
    type(bound_area) , intent(in) :: ba
    type(v3dbound) , intent(in) :: bnd
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    integer(ik4) :: i , j , k
    integer(ik4) :: ib , i1 , i2 , j1 , j2
    real(rk8) , pointer , dimension(:) :: wg
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba%dotflag ) then
      wg => wgtd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      wg => wgtx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if

    if ( ba%ns /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bsouth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k) = wg(ib)*ften(j,i,k) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%nn /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bnorth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k) = wg(ib)*ften(j,i,k) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%nw /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bwest(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k) = wg(ib)*ften(j,i,k) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
    if ( ba%ne /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%beast(j,i) ) cycle
            ib = ba%ibnd(j,i)
            ften(j,i,k) = wg(ib)*ften(j,i,k) + (d_one-wg(ib))*bnd%bt(j,i,k)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine sponge3d

  subroutine sponge2d(ba,bnd,ften)
    implicit none
    type(bound_area) , intent(in) :: ba
    type(v2dbound) , intent(in) :: bnd
    real(rk8) , pointer , intent(inout) , dimension(:,:) :: ften
    integer(ik4) :: i , j , i1 , i2 , j1 , j2
    integer(ik4) :: ib
    real(rk8) , pointer , dimension(:) :: wg
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'sponge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba%dotflag ) then
      wg => wgtd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      wg => wgtx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if

    if ( ba%ns /= 0 ) then
      do i = i1 , i2
        do j = j1 , j2
          if ( .not. ba%bsouth(j,i) ) cycle
          ib = ba%ibnd(j,i)
          ften(j,i) = wg(ib)*ften(j,i) + (d_one-wg(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba%nn /= 0 ) then
      do i = i1 , i2
        do j = j1 , j2
          if ( .not. ba%bnorth(j,i) ) cycle
          ib = ba%ibnd(j,i)
          ften(j,i) = wg(ib)*ften(j,i) + (d_one-wg(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba%nw /= 0 ) then
      do i = i1 , i2
        do j = j1 , j2
          if ( .not. ba%bwest(j,i) ) cycle
          ib = ba%ibnd(j,i)
          ften(j,i) = wg(ib)*ften(j,i) + (d_one-wg(ib))*bnd%bt(j,i)
        end do
      end do
    end if
    if ( ba%ne /= 0 ) then
      do i = i1 , i2
        do j = j1 , j2
          if ( .not. ba%beast(j,i) ) cycle
          ib = ba%ibnd(j,i)
          ften(j,i) = wg(ib)*ften(j,i) + (d_one-wg(ib))*bnd%bt(j,i)
        end do
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
  subroutine nudge4d(nk,n1,n2,ba,f,ibdy,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: ibdy , nk , n1 , n2
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    type(bound_area) , intent(in) :: ba
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
    real(rk8) :: xt , xf , fls0 , fls1 , fls2 , fls3 , fls4 , xg
    integer(ik4) :: i , j , k , ib , i1 , i2 , j1 , j2 , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    xt = xbctime + dt
    if ( ba%dotflag ) then
      lfc => fcd
      lgc => gcd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      lfc => fcx
      lgc => gcx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if
    if ( nk == size(hefc) ) then
      efc => hefc
      egc => hegc
    else
      efc => fefc
      egc => fegc
    end if

    if ( ba%ns /= 0 ) then
      do n = n1 , n2
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bsouth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              if ( ibdy == 1 ) then
                xf = lfc(ib)
                xg = lgc(ib)
              else
                xf = efc(ib,k)
                xg = egc(ib,k)
              end if
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k,n)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k,n)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k,n)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k,n)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( ba%nn /= 0 ) then
      do n = n1 , n2
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bnorth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              if ( ibdy == 1 ) then
                xf = lfc(ib)
                xg = lgc(ib)
              else
                xf = efc(ib,k)
                xg = egc(ib,k)
              end if
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k,n)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k,n)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k,n)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k,n)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( ba%nw /= 0 ) then
      do n = n1 , n2
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bwest(j,i) ) cycle
              ib = ba%ibnd(j,i)
              if ( ibdy == 1 ) then
                xf = lfc(ib)
                xg = lgc(ib)
              else
                xf = efc(ib,k)
                xg = egc(ib,k)
              end if
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k,n)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k,n)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k,n)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k,n)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
    if ( ba%ne /= 0 ) then
      do n = n1 , n2
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%beast(j,i) ) cycle
              ib = ba%ibnd(j,i)
              if ( ibdy == 1 ) then
                xf = lfc(ib)
                xg = lgc(ib)
              else
                xf = efc(ib,k)
                xg = egc(ib,k)
              end if
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k,n)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k,n)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k,n)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k,n)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k,n)
              ften(j,i,k,n) = ften(j,i,k,n) + xf*fls0 -  &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge4d

  subroutine nudge3d(nk,ba,f,ibdy,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: ibdy , nk
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    type(bound_area) , intent(in) :: ba
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rk8) :: xt , xf , fls0 , fls1 , fls2 , fls3 , fls4 , xg
    integer(ik4) :: i , j , k , ib , i1 , i2 , j1 , j2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    xt = xbctime + dt
    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba%dotflag ) then
      lfc => fcd
      lgc => gcd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      lfc => fcx
      lgc => gcx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if
    if ( nk == size(hefc) ) then
      efc => hefc
      egc => hegc
    else
      efc => fefc
      egc => fegc
    end if

    if ( ibdy == 1 ) then
      if ( ba%ns /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bsouth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = lfc(ib)
              xg = lgc(ib)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%nn /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bnorth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = lfc(ib)
              xg = lgc(ib)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%nw /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bwest(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = lfc(ib)
              xg = lgc(ib)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%ne /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%beast(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = lfc(ib)
              xg = lgc(ib)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 -  &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    else
      if ( ba%ns /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bsouth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = efc(ib,k)
              xg = egc(ib,k)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%nn /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bnorth(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = efc(ib,k)
              xg = egc(ib,k)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls2 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              fls3 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls4 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%nw /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%bwest(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = efc(ib,k)
              xg = egc(ib,k)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 - &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
      if ( ba%ne /= 0 ) then
        do k = 1 , nk
          do i = i1 , i2
            do j = j1 , j2
              if ( .not. ba%beast(j,i) ) cycle
              ib = ba%ibnd(j,i)
              xf = efc(ib,k)
              xg = egc(ib,k)
              fls0 = (bnd%b0(j,i,k)  +xt*bnd%bt(j,i,k))   - f(j,i,k)
              fls1 = (bnd%b0(j,i-1,k)+xt*bnd%bt(j,i-1,k)) - f(j,i-1,k)
              fls2 = (bnd%b0(j,i+1,k)+xt*bnd%bt(j,i+1,k)) - f(j,i+1,k)
              fls3 = (bnd%b0(j-1,i,k)+xt*bnd%bt(j-1,i,k)) - f(j-1,i,k)
              fls4 = (bnd%b0(j+1,i,k)+xt*bnd%bt(j+1,i,k)) - f(j+1,i,k)
              ften(j,i,k) = ften(j,i,k) + xf*fls0 -  &
                            xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
            end do
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge3d

  subroutine nudge2d(ba,f,ibdy,bnd,ften)
    implicit none
    integer(ik4) , intent(in) :: ibdy
    real(rk8) , pointer , intent(in) , dimension(:,:) :: f
    type(v2dbound) , intent(in) :: bnd
    type(bound_area) , intent(in) :: ba
    real(rk8) , pointer , intent(inout) , dimension(:,:) :: ften
    real(rk8) :: xt , xf , fls0 , fls1 , fls2 , fls3 , fls4 , xg
    integer(ik4) :: i , j , ib , i1 , i2 , j1 , j2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'nudge2d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    xt = xbctime + dt
    if ( .not. ba%havebound ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ba%dotflag ) then
      lfc => fcd
      lgc => gcd
      i1 = idi1
      i2 = idi2
      j1 = jdi1
      j2 = jdi2
    else
      lfc => fcx
      lgc => gcx
      i1 = ici1
      i2 = ici2
      j1 = jci1
      j2 = jci2
    end if

    if ( ibdy == 1 ) then
      if ( ba%ns /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bsouth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = lfc(ib)
            xg = lgc(ib)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls2 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            fls3 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls4 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%nn /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bnorth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = lfc(ib)
            xg = lgc(ib)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls2 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            fls3 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls4 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%nw /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bwest(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = lfc(ib)
            xg = lgc(ib)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls2 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            fls3 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls4 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%ne /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%beast(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = lfc(ib)
            xg = lgc(ib)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls2 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            fls3 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls4 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            ften(j,i) = ften(j,i) + xf*fls0 -  &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
    else
      if ( ba%ns /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bsouth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls2 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            fls3 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls4 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%nn /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bnorth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls2 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            fls3 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls4 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%nw /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bwest(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls2 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            fls3 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls4 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            ften(j,i) = ften(j,i) + xf*fls0 - &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
      if ( ba%ne /= 0 ) then
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%beast(j,i) ) cycle
            ib = ba%ibnd(j,i)
            xf = hefc(ib,kz)
            xg = hegc(ib,kz)
            fls0 = (bnd%b0(j,i)  +xt*bnd%bt(j,i))   - f(j,i)
            fls1 = (bnd%b0(j,i-1)+xt*bnd%bt(j,i-1)) - f(j,i-1)
            fls2 = (bnd%b0(j,i+1)+xt*bnd%bt(j,i+1)) - f(j,i+1)
            fls3 = (bnd%b0(j-1,i)+xt*bnd%bt(j-1,i)) - f(j-1,i)
            fls4 = (bnd%b0(j+1,i)+xt*bnd%bt(j+1,i)) - f(j+1,i)
            ften(j,i) = ften(j,i) + xf*fls0 -  &
                          xg*rdxsq*(fls1+fls2+fls3+fls4-d_four*fls0)
          end do
        end do
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine nudge2d

  subroutine couple(a,c,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: a
    real(rk8) , pointer , dimension(:,:) , intent(in) :: c
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

  subroutine timeint2(a,b,c,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: a , b
    real(rk8) , pointer , dimension(:,:) , intent(out) :: c
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
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: a , b
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: c
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

end module mod_bdycod

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
