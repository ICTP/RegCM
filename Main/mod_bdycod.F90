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
! Storage and subroutines for input of boundary values and tendencies
! of p*u, p*v, p*t, p*qv and p*, and outmost 2 slices of u and v for
! large domain.
! Relaxation and Sponge Boundary Conditions routines
!
  use mod_dynparam
  use mod_mppparam
  use mod_memutil
  use mod_atm_interface
  use mod_pbl_interface , only : set_tke_bc
  use mod_lm_interface
  use mod_mpmessage 
  use mod_ncio
  use mod_mppio
  use mod_service
!
  private
!
  public :: allocate_mod_bdycon , init_bdy , bdyin , bdyval
!
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
  public :: ts0
  public :: ts1 ! FOR DCSST
! fnudge : are the coefficients for the newtonian term.
! gnydge : are the coefficients for the diffusion term.
  public :: fnudge , gnudge
  public :: anudg
!
  real(dp) , pointer , dimension(:,:) :: sue , sui , nue , nui , &
                                         sve , svi , nve , nvi
  real(dp) , pointer , dimension(:,:) :: wue , wui , eue , eui , &
                                         wve , wvi , eve , evi
  real(dp) , pointer , dimension(:,:) :: ts0 , ts1
  real(dp) , pointer , dimension(:) :: anudg
  real(dp) , pointer , dimension(:) :: fcx , gcx
  real(dp) , pointer , dimension(:) :: fcd , gcd
  real(dp) , pointer , dimension(:) :: lfc , lgc
  real(dp) , pointer , dimension(:,:) :: efc , egc
  real(dp) , pointer , dimension(:) :: wgtd
  real(dp) , pointer , dimension(:) :: wgtx
  real(dp) :: fnudge , gnudge
  integer :: nbdm
!
  interface nudge
    module procedure nudge3d , nudge2d
  end interface nudge
!
  interface sponge
    module procedure sponge3d , sponge2d
  end interface sponge
!
  public :: sponge , nudge , setup_bdycon
!
  contains
!
  subroutine allocate_mod_bdycon(iboudy,lband)
    implicit none
    integer , intent(in) :: iboudy
    logical , intent(in) :: lband
    character (len=64) :: subroutine_name='allocate_mod_bdycon'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

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
      call getmem1d(anudg,1,kz,'bdycon:anudg')
      call getmem2d(efc,1,nbdm,1,kz,'bdycon:fcx')
      call getmem2d(egc,1,nbdm,1,kz,'bdycon:fcx')
    end if

    call getmem2d(ts0,1,jxp,icross1,icross2,'bdycon:ts0')
    call getmem2d(ts1,1,jxp,icross1,icross2,'bdycon:ts1')
!
    call getmem2d(sue,0,jxp+1,1,kz,'bdycon:sue')
    call getmem2d(sui,0,jxp+1,1,kz,'bdycon:sui')
    call getmem2d(nue,0,jxp+1,1,kz,'bdycon:nue')
    call getmem2d(nui,0,jxp+1,1,kz,'bdycon:nui')
    call getmem2d(nve,0,jxp+1,1,kz,'bdycon:nve')
    call getmem2d(nvi,0,jxp+1,1,kz,'bdycon:nvi')
    call getmem2d(sve,0,jxp+1,1,kz,'bdycon:sve')
    call getmem2d(svi,0,jxp+1,1,kz,'bdycon:svi')
!
    if ( .not. lband ) then
      call getmem2d(wue,idot1,idot2,1,kz,'bdycon:wue')
      call getmem2d(wui,idot1,idot2,1,kz,'bdycon:wui')
      call getmem2d(eue,idot1,idot2,1,kz,'bdycon:eue')
      call getmem2d(eui,idot1,idot2,1,kz,'bdycon:eui')
      call getmem2d(wve,idot1,idot2,1,kz,'bdycon:wve')
      call getmem2d(wvi,idot1,idot2,1,kz,'bdycon:wvi')
      call getmem2d(eve,idot1,idot2,1,kz,'bdycon:eve')
      call getmem2d(evi,idot1,idot2,1,kz,'bdycon:evi')
    end if
    call time_end(subroutine_name,idindx)

  end subroutine allocate_mod_bdycon
!
  subroutine setup_bdycon(hlev)
    implicit none
    real(dp) , pointer , dimension(:) , intent(in) :: hlev
    integer :: n , k
    character (len=64) :: subroutine_name='setup_bdycon'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! Specify the coefficients for nudging boundary conditions:
    !
    if ( iboudy == 1 .or. iboudy == 5 ) then
      fnudge = 0.1D0/dt2
      gnudge = (dxsq/dt)/50.0D0
    end if
    if ( iboudy == 1 ) then
      do n = 1 , nspgx
        fcx(n) = fnudge*xfun(n,.false.)
        gcx(n) = gnudge*xfun(n,.false.)
      end do
      do n = 1 , nspgd
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
      do k = 1 , kz
        if ( hlev(k) < 0.4D0 ) then
          anudg(k) = high_nudge
        else if ( hlev(k) < 0.8D0 ) then
          anudg(k) = medium_nudge
        else
          anudg(k) = low_nudge
        end if
        do n = 1 , nbdm
          efc(n,k) = fnudge*xfune(n,k)
          egc(n,k) = gnudge*xfune(n,k)
        end do
      end do
    end if
    call time_end(subroutine_name,idindx)
  end subroutine setup_bdycon

  subroutine init_bdy
    implicit none
    integer :: datefound , i , j , k , ierr
    character(len=32) :: appdat
    type (rcm_time_and_date) :: icbc_date
    character (len=64) :: subroutine_name='init_bdy'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    bdydate1 = idate1
    bdydate2 = idate1

    if ( myid == 0 ) then
      if ( bdydate1 == globidate1 ) then
        icbc_date = bdydate1
      else
        icbc_date = monfirst(bdydate1)
      end if

      call open_icbc(icbc_date)

      datefound = icbc_search(bdydate1)
      if (datefound < 0) then
        !
        ! Cannot run without initial conditions
        !
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
      end if

      call read_icbc(ps0_io,ts0_io,ub0_io,vb0_io,tb0_io,qb0_io)

      appdat = tochar(bdydate1)
      if ( .not. ifrest ) then
        write (6,*) 'READY IC DATA for ', appdat
      else
        write (6,*) 'READY BC DATA for ', appdat
      end if

      bdydate2 = bdydate2 + intbdy
      write (6,*) 'SEARCH BC data for ', toint10(bdydate2)
      datefound = icbc_search(bdydate2)
      if (datefound < 0) then
        call open_icbc(monfirst(bdydate2))
        datefound = icbc_search(bdydate2)
        if (datefound < 0) then
          appdat = tochar(bdydate2)
          call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
        end if
      end if
      call read_icbc(ps1_io,ts1_io,ub1_io,vb1_io,tb1_io,qb1_io)

      ps0_io(:,:) = (ps0_io(:,:)*d_r10)-ptop
      ps1_io(:,:) = (ps1_io(:,:)*d_r10)-ptop

      write (6,*) 'READY  BC from     ' , &
            toint10(bdydate1) , ' to ' , toint10(bdydate2)

    end if

    call date_bcast(bdydate2,0,mycomm,ierr)
    bdydate1 = bdydate2
    !
    ! Send each processor its computing slice
    !
    call deco1_scatter(ub0_io,xub%b0,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(vb0_io,xvb%b0,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(tb0_io,xtb%b0,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(qb0_io,xqb%b0,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(ps0_io,xpsb%b0,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ts0_io,ts0,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ub1_io,xub%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(vb1_io,xvb%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(tb1_io,xtb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(qb1_io,xqb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(ps1_io,xpsb%b1,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ts1_io,ts1,jcross1,jcross2,icross1,icross2)
    !
    ! Calculate P* on dot points
    !
    call deco1_exchange_left(xpsb%b0,1,ice1,ice2)
    call deco1_exchange_right(xpsb%b0,1,ice1,ice2)
    call psc2psd(xpsb%b0,psdot)
    !
    ! Couple pressure u,v,t,q
    !
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%b0(j,i,k) = xub%b0(j,i,k)*psdot(j,i)
          xvb%b0(j,i,k) = xvb%b0(j,i,k)*psdot(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xub%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%b0,1,ide1,ide2,1,kz)
!
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%b0(j,i,k) = xtb%b0(j,i,k)*xpsb%b0(j,i)
          xqb%b0(j,i,k) = xqb%b0(j,i,k)*xpsb%b0(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xtb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%b0,1,ice1,ice2,1,kz)
    !
    ! Repeat for T2
    !
    call deco1_exchange_left(xpsb%b1,1,ice1,ice2)
    call deco1_exchange_right(xpsb%b1,1,ice1,ice2)
    call psc2psd(xpsb%b1,psdot)
    !
    ! Couple pressure u,v,t,q
    !
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%b1(j,i,k) = xub%b1(j,i,k)*psdot(j,i)
          xvb%b1(j,i,k) = xvb%b1(j,i,k)*psdot(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xub%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%b1,1,ide1,ide2,1,kz)
!
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%b1(j,i,k) = xtb%b1(j,i,k)*xpsb%b1(j,i)
          xqb%b1(j,i,k) = xqb%b1(j,i,k)*xpsb%b1(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xtb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%b1,1,ice1,ice2,1,kz)
    !
    ! Calculate time varying component
    !
    do i = ice1 , ice2
      do j = jce1 , jce2
        xpsb%bt(j,i) = (xpsb%b1(j,i)-xpsb%b0(j,i))/dtbdys
      end do
    end do
    call deco1_exchange_left(xpsb%bt,1,ice1,ice2)
    call deco1_exchange_right(xpsb%bt,1,ice1,ice2)
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%bt(j,i,k) = (xub%b1(j,i,k)-xub%b0(j,i,k))/dtbdys
          xvb%bt(j,i,k) = (xvb%b1(j,i,k)-xvb%b0(j,i,k))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(xub%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%bt,1,ide1,ide2,1,kz)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%bt(j,i,k) = (xtb%b1(j,i,k)-xtb%b0(j,i,k))/dtbdys
          xqb%bt(j,i,k) = (xqb%b1(j,i,k)-xqb%b0(j,i,k))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(xtb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%bt,1,ice1,ice2,1,kz)
!
    call time_end(subroutine_name,idindx)

  end subroutine init_bdy

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! this subroutine reads in the boundary conditions.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bdyin
    implicit none
!
    integer :: i , j , k , n , mmrec
    integer :: ierr
    character(len=32) :: appdat
    character (len=64) :: subroutine_name='bdyin'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    xub%b0(:,:,:) = xub%b1(:,:,:)
    xvb%b0(:,:,:) = xvb%b1(:,:,:)
    xtb%b0(:,:,:) = xtb%b1(:,:,:)
    xqb%b0(:,:,:) = xqb%b1(:,:,:)
    xpsb%b0(:,:) = xpsb%b1(:,:)
    ts0(:,:) = ts1(:,:)
!
    if ( myid == 0 ) then
      bdydate2 = bdydate2 + intbdy
      write (6,'(a,i10)') 'SEARCH BC data for ', toint10(bdydate2)
      mmrec = icbc_search(bdydate2)
      if (mmrec < 0) then
        call open_icbc(monfirst(bdydate2))
        mmrec = icbc_search(bdydate2)
        if (mmrec < 0) then
          appdat = tochar(bdydate2)
          call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
        end if
      end if
      call read_icbc(ps1_io,ts1_io,ub1_io,vb1_io,tb1_io,qb1_io)
      !
      ! Convert surface pressure to pstar
      !
      ps1_io(:,:) = (ps1_io(:,:)*d_r10)-ptop
    end if
!
    call deco1_scatter(ub1_io,xub%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(vb1_io,xvb%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(tb1_io,xtb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(qb1_io,xqb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(ps1_io,xpsb%b1,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ts1_io,ts1,jcross1,jcross2,icross1,icross2)
!
    call deco1_exchange_left(xpsb%b1,1,ice1,ice2)
    call deco1_exchange_right(xpsb%b1,1,ice1,ice2)
    call psc2psd(xpsb%b1,psdot)
!
!   Couple pressure u,v,t,q
!
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%b1(j,i,k) = xub%b1(j,i,k)*psdot(j,i)
          xvb%b1(j,i,k) = xvb%b1(j,i,k)*psdot(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xub%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%b1,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%b1,1,ide1,ide2,1,kz)
!
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%b1(j,i,k) = xtb%b1(j,i,k)*xpsb%b1(j,i)
          xqb%b1(j,i,k) = xqb%b1(j,i,k)*xpsb%b1(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xtb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%b1,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%b1,1,ice1,ice2,1,kz)
!
    do i = ice1 , ice2
      do j = jce1 , jce2
        xpsb%bt(j,i) = (xpsb%b1(j,i)-xpsb%b0(j,i))/dtbdys
      end do
    end do
    call deco1_exchange_left(xpsb%bt,1,ice1,ice2)
    call deco1_exchange_right(xpsb%bt,1,ice1,ice2)
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%bt(j,i,k) = (xub%b1(j,i,k)-xub%b0(j,i,k))/dtbdys
          xvb%bt(j,i,k) = (xvb%b1(j,i,k)-xvb%b0(j,i,k))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(xub%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%bt,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%bt,1,ide1,ide2,1,kz)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%bt(j,i,k) = (xtb%b1(j,i,k)-xtb%b0(j,i,k))/dtbdys
          xqb%bt(j,i,k) = (xqb%b1(j,i,k)-xqb%b0(j,i,k))/dtbdys
        end do
      end do
    end do
    call deco1_exchange_left(xtb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%bt,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%bt,1,ice1,ice2,1,kz)
!
    !
    ! Update ground temperature on Ocean/Lakes
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( iswater(mddom%lndcat(j,i)) ) then
          if (idcsst == 1) then
            sfs%tga(j,i) = ts1(j,i) + dtskin(j,i)
            sfs%tgb(j,i) = ts1(j,i) + dtskin(j,i)
          else
            sfs%tga(j,i) = ts1(j,i)
            sfs%tgb(j,i) = ts1(j,i)
          end if
          if ( iseaice == 1 ) then
            if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
            if ( ts1(j,i) <= icetemp ) then
              sfs%tga(j,i) = icetemp
              sfs%tgb(j,i) = icetemp
              ts1(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ldmsk1(n,j,i) = 2
                sfice(n,j,i) = d_1000
                sncv(n,j,i) = d_zero
              end do
            else
              sfs%tga(j,i) = ts1(j,i)
              sfs%tgb(j,i) = ts1(j,i)
              ldmsk(j,i) = 0
              do n = 1, nnsg
                ldmsk1(n,j,i) = 0
                sfice(n,j,i) = d_zero
                sncv(n,j,i)  = d_zero
              end do
            end if
          end if
        end if
      end do
    end do

    if ( myid == 0 ) then
      write (6,'(a,i10,a,i10)') 'READY  BC from     ' , &
            toint10(bdydate1) , ' to ' , toint10(bdydate2)
    end if

    call date_bcast(bdydate2,0,mycomm,ierr)
    bdydate1 = bdydate2

    call time_end(subroutine_name,idindx)

  end subroutine bdyin
!
! This subroutine sets the boundary values of u and v according
! to the boundary conditions specified.
!
!     iboudy = 0 : fixed
!            = 1 : relaxation, linear technique
!            = 2 : time dependent
!            = 3 : time dependent and inflow/outflow dependent
!            = 4 : sponge
!            = 5 : relaxation, exponential technique
!
!     dtb        : elapsed time from the initial boundary values.
!
  subroutine bdyuv(iboudy,dtb)
!
    real(dp) , intent(in) :: dtb
    integer , intent(in) :: iboudy
!
    integer :: i , j , k
    character (len=64) :: subroutine_name='bdyuv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! First compute the p* at dot points to decouple U,V:
    !
    call deco1_exchange_left(sfs%psa,1,ice1,ice2)
    call deco1_exchange_right(sfs%psa,1,ice1,ice2)
    call psc2psd(sfs%psa,psdot)
    !
    ! Now compute last two points values in U and V
    ! Internal points
    !
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = idi1 , idi2
          wui(i,k) = atm1%u(jdi1,i,k)/psdot(jdi1,i)
          wvi(i,k) = atm1%v(jdi1,i,k)/psdot(jdi1,i)
        end do
      end do
    end if
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = idi1 , idi2
          eui(i,k) = atm1%u(jdi2,i,k)/psdot(jdi2,i)
          evi(i,k) = atm1%v(jdi2,i,k)/psdot(jdi2,i)
        end do
      end do
    end if
    if ( ma%hasbottom ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          sui(j,k) = atm1%u(j,idi1,k)/psdot(j,idi1)
          svi(j,k) = atm1%v(j,idi1,k)/psdot(j,idi1)
        end do
      end do
    end if
    if ( ma%hastop ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          nui(j,k) = atm1%u(j,idi2,k)/psdot(j,idi2)
          nvi(j,k) = atm1%v(j,idi2,k)/psdot(j,idi2)
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
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ide1 , ide2
            wue(i,k) = xub%b0(jde1,i,k)/psdot(jde1,i)
            wve(i,k) = xvb%b0(jde1,i,k)/psdot(jde1,i)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ide1 , ide2
            eue(i,k) = xub%b0(jde2,i,k)/psdot(jde2,i)
            eve(i,k) = xvb%b0(jde2,i,k)/psdot(jde2,i)
          end do
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jde1 , jde2
            sue(j,k) = xub%b0(j,ide1,k)/psdot(j,ide1)
            sve(j,k) = xvb%b0(j,ide1,k)/psdot(j,ide1)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jde1 , jde2
            nue(j,k) = xub%b0(j,ide2,k)/psdot(j,ide2)
            nve(j,k) = xvb%b0(j,ide2,k)/psdot(j,ide2)
          end do
        end do
      end if
!
    else ! NOT Fixed
!
      !
      !     time-dependent boundary conditions:
      !
      ! west (j = 1) and east (j = jx) boundaries:
      !
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ide1 , ide2
            wue(i,k) = (xub%b0(jde1,i,k)+dtb*xub%bt(jde1,i,k))/psdot(jde1,i)
            wve(i,k) = (xvb%b0(jde1,i,k)+dtb*xvb%bt(jde1,i,k))/psdot(jde1,i)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ide1 , ide2
            eue(i,k) = (xub%b0(jde2,i,k)+dtb*xub%bt(jde2,i,k))/psdot(jde2,i)
            eve(i,k) = (xvb%b0(jde2,i,k)+dtb*xvb%bt(jde2,i,k))/psdot(jde2,i)
          end do
        end do
      end if
      !
      ! south and north boundaries:
      !
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jde1 , jde2
            sue(j,k) = (xub%b0(j,ide1,k)+dtb*xub%bt(j,ide1,k))/psdot(j,ide1)
            sve(j,k) = (xvb%b0(j,ide1,k)+dtb*xvb%bt(j,ide1,k))/psdot(j,ide1)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jde1 , jde2
            nue(j,k) = (xub%b0(j,ide2,k)+dtb*xub%bt(j,ide2,k))/psdot(j,ide2)
            nve(j,k) = (xvb%b0(j,ide2,k)+dtb*xvb%bt(j,ide2,k))/psdot(j,ide2)
          end do
        end do
      end if
    end if
    !
    ! fill up the interior silces:
    !
    if ( ma%hasleft ) then
      do k = 1 , kz
        wui(ide1,k) = sue(jdi1,k)
        wui(ide2,k) = nue(jdi1,k)
        wvi(ide1,k) = sve(jdi1,k)
        wvi(ide2,k) = nve(jdi1,k)
        sui(jde1,k) = wue(idi1,k)
        nui(jde1,k) = wue(idi2,k)
        svi(jde1,k) = wve(idi1,k)
        nvi(jde1,k) = wve(idi2,k)
      end do
    end if
    if ( ma%hasright ) then
      do k = 1 , kz
        eui(ide1,k) = sue(jdi2,k)
        eui(ide2,k) = nue(jdi2,k)
        evi(ide1,k) = sve(jdi2,k)
        evi(ide2,k) = nve(jdi2,k)
        sui(jde2,k) = eue(idi1,k)
        nui(jde2,k) = eue(idi2,k)
        svi(jde2,k) = eve(idi1,k)
        nvi(jde2,k) = eve(idi2,k)
      end do
    end if

    call deco1_exchange_left(sue,1,1,kz)
    call deco1_exchange_right(sue,1,1,kz)
    call deco1_exchange_left(sui,1,1,kz)
    call deco1_exchange_right(sui,1,1,kz)
    call deco1_exchange_left(nue,1,1,kz)
    call deco1_exchange_right(nue,1,1,kz)
    call deco1_exchange_left(nui,1,1,kz)
    call deco1_exchange_right(nui,1,1,kz)
    call deco1_exchange_left(sve,1,1,kz)
    call deco1_exchange_right(sve,1,1,kz)
    call deco1_exchange_left(svi,1,1,kz)
    call deco1_exchange_right(svi,1,1,kz)
    call deco1_exchange_left(nve,1,1,kz)
    call deco1_exchange_right(nve,1,1,kz)
    call deco1_exchange_left(nvi,1,1,kz)
    call deco1_exchange_right(nvi,1,1,kz)

    call time_end(subroutine_name,idindx)

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
  subroutine bdyval(xt)
!
    implicit none
!
    real(dp) , intent(in) :: xt
!
    real(dp) :: qcx , qcint , qvx , qext , qint , vavg
    integer :: i , j , k
    real(dp) :: uavg
    character (len=64) :: subroutine_name='bdyval'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! Fill up the boundary value for xxb variables from xxa variables:
    ! if this subroutine is called for the first time, this part
    ! shall be skipped.
    !
    if ( ktau > 1 ) then
      !
      ! West boundary
      !
      if ( ma%hasleft ) then
        do i = ici1 , ici2
          sfs%psb(jce1,i) = sfs%psa(jce1,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm2%u(jde1,i,k) = atm1%u(jde1,i,k)/mddom%msfd(jde1,i)
            atm2%v(jde1,i,k) = atm1%v(jde1,i,k)/mddom%msfd(jde1,i)
          end do
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            atm2%t(jce1,i,k) = atm1%t(jce1,i,k)
            atm2%qv(jce1,i,k) = atm1%qv(jce1,i,k)
            atm2%qc(jce1,i,k) = atm1%qc(jce1,i,k)
          end do
        end do
      end if
      !
      ! East boundary
      !
      if ( ma%hasright ) then
        do i = ici1 , ici2
          sfs%psb(jce2,i) = sfs%psa(jce2,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm2%u(jde2,i,k) = atm1%u(jde2,i,k)/mddom%msfd(jde2,i)
            atm2%v(jde2,i,k) = atm1%v(jde2,i,k)/mddom%msfd(jde2,i)
          end do
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            atm2%t(jce2,i,k) = atm1%t(jce2,i,k)
            atm2%qv(jce2,i,k) = atm1%qv(jce2,i,k)
            atm2%qc(jce2,i,k) = atm1%qc(jce2,i,k)
          end do
        end do
      end if
      !
      ! North and South boundaries
      !
      if ( ma%hasbottom ) then
        do j = jce1 , jce2
          sfs%psb(j,ice1) = sfs%psa(j,ice1)
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atm2%u(j,ice1,k) = atm1%u(j,ice1,k)/mddom%msfd(j,ice1)
            atm2%v(j,ice1,k) = atm1%v(j,ice1,k)/mddom%msfd(j,ice1)
          end do
        end do
        do k = 1 , kz
          do j = jce1 , jce2
            atm2%t(j,ice1,k) = atm1%t(j,ice1,k)
            atm2%qv(j,ice1,k) = atm1%qv(j,ice1,k)
            atm2%qc(j,ice1,k) = atm1%qc(j,ice1,k)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do j = jce1 , jce2
          sfs%psb(j,ice2) = sfs%psa(j,ice2)
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atm2%u(j,ide2,k) = atm1%u(j,ide2,k)/mddom%msfd(j,ide2)
            atm2%v(j,ide2,k) = atm1%v(j,ide2,k)/mddom%msfd(j,ide2)
          end do
        end do
        do k = 1 , kz
          do j = jce1 , jce2
            atm2%t(j,ice2,k) = atm1%t(j,ice2,k)
            atm2%qv(j,ice2,k) = atm1%qv(j,ice2,k)
            atm2%qc(j,ice2,k) = atm1%qc(j,ice2,k)
          end do
        end do
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
      if ( ma%hasleft ) then
        do i = ici1 , ici2
          sfs%psa(jce1,i) = xpsb%b0(jce1,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm1%u(jde1,i,k) = xub%b0(jde1,i,k)
            atm1%v(jde1,i,k) = xvb%b0(jde1,i,k)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do i = ici1 , ici2
          sfs%psa(jce2,i) = xpsb%b0(jce2,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm1%u(jde2,i,k) = xub%b0(jde2,i,k)
            atm1%v(jde2,i,k) = xvb%b0(jde2,i,k)
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do j = jce1 , jce2
          sfs%psa(j,ice1) = xpsb%b0(j,ice1)
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atm1%u(j,ide1,k) = xub%b0(j,ide1,k)
            atm1%v(j,ide1,k) = xvb%b0(j,ide1,k)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do j = jce1 , jce2
          sfs%psa(j,ice2) = xpsb%b0(j,ice2)
        end do
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
      if ( ma%hasleft ) then
        do i = ici1 , ici2
          sfs%psa(jce1,i) = xpsb%b0(jce1,i) + xt*xpsb%bt(jce1,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm1%u(jde1,i,k) = xub%b0(jde1,i,k) + xt*xub%bt(jde1,i,k)
            atm1%v(jde1,i,k) = xvb%b0(jde1,i,k) + xt*xvb%bt(jde1,i,k)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do i = ici1 , ici2
          sfs%psa(jce2,i) = xpsb%b0(jce2,i) + xt*xpsb%bt(jce2,i)
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atm1%u(jde2,i,k) = xub%b0(jde2,i,k) + xt*xub%bt(jde2,i,k)
            atm1%v(jde2,i,k) = xvb%b0(jde2,i,k) + xt*xvb%bt(jde2,i,k)
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do j = jce1 , jce2
          sfs%psa(j,ice1) = xpsb%b0(j,ice1) + xt*xpsb%bt(j,ice1)
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atm1%u(j,ide1,k) = xub%b0(j,ide1,k) + xt*xub%bt(j,ide1,k)
            atm1%v(j,ide1,k) = xvb%b0(j,ide1,k) + xt*xvb%bt(j,ide1,k)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do j = jce1 , jce2
          sfs%psa(j,ice2) = xpsb%b0(j,ice2) + xt*xpsb%bt(j,ice2)
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atm1%u(j,ide2,k) = xub%b0(j,ide2,k) + xt*xub%bt(j,ide2,k)
            atm1%v(j,ide2,k) = xvb%b0(j,ide2,k) + xt*xvb%bt(j,ide2,k)
          end do
        end do
      end if
    end if
!
    call bdyuv(iboudy,xt)
!
    !
    ! Set boundary values for p*t:
    ! Set boundary values for p*qv:
    !
    if ( iboudy == 0 ) then
      !
      ! fixed boundary conditions:
      !
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            atm1%t(jce1,i,k) = xtb%b0(jce1,i,k)
            atm1%qv(jce1,i,k) = xqb%b0(jce1,i,k)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            atm1%t(jce2,i,k) = xtb%b0(jce2,i,k)
            atm1%qv(jce2,i,k) = xqb%b0(jce2,i,k)
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            atm1%t(j,ice1,k) = xtb%b0(j,ice1,k)
            atm1%qv(j,ice1,k) = xqb%b0(j,ice1,k)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            atm1%t(j,ice2,k) = xtb%b0(j,ice2,k)
            atm1%qv(j,ice2,k) = xqb%b0(j,ice2,k)
          end do
        end do
      end if
    else
      !
      ! time-dependent boundary conditions:
      !
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            atm1%t(jce1,i,k) = xtb%b0(jce1,i,k) + xt*xtb%bt(jce1,i,k)
            atm1%qv(jce1,i,k) = xqb%b0(jce1,i,k) + xt*xqb%bt(jce1,i,k)
          end do
        end do
      end if
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            atm1%t(jce2,i,k) = xtb%b0(jce2,i,k) + xt*xtb%bt(jce2,i,k)
            atm1%qv(jce2,i,k) = xqb%b0(jce2,i,k) + xt*xqb%bt(jce2,i,k)
          end do
        end do
      end if
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            atm1%t(j,ice1,k) = xtb%b0(j,ice1,k) + xt*xtb%bt(j,ice1,k)
            atm1%qv(j,ice1,k) = xqb%b0(j,ice1,k) + xt*xqb%bt(j,ice1,k)
          end do
        end do
      end if
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            atm1%t(j,ice2,k) = xtb%b0(j,ice2,k) + xt*xtb%bt(j,ice2,k)
            atm1%qv(j,ice2,k) = xqb%b0(j,ice2,k) + xt*xqb%bt(j,ice2,k)
          end do
        end do
      end if
    end if
!
    if ( iboudy == 3 .or. iboudy == 4 ) then
      !
      ! determine QV boundary values depends on inflow/outflow:
      !
      ! west boundary:
      !
      if ( ma%hasleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            qext = atm1%qv(jce1,i,k)/sfs%psa(jce1,i)
            qint = atm1%qv(jci1,i,k)/sfs%psa(jci1,i)
            uavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
            if ( uavg >= d_zero ) then
              qvx = qext
            else
              qvx = qint
            end if
            atm1%qv(jce1,i,k) = qvx*sfs%psa(jce1,i)
          end do
        end do
      end if
      !
      ! east boundary:
      !
      if ( ma%hasright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            qext = atm1%qv(jce2,i,k)/sfs%psa(jce2,i)
            qint = atm1%qv(jci2,i,k)/sfs%psa(jci2,i)
            uavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
            if ( uavg < d_zero ) then
              qvx = qext
            else
              qvx = qint
            end if
            atm1%qv(jce2,i,k) = qvx*sfs%psa(jce2,i)
          end do
        end do
      end if
      !
      ! south boundary:
      !
      if ( ma%hasbottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            qext = atm1%qv(j,ice1,k)/sfs%psa(j,ice1)
            qint = atm1%qv(j,ici1,k)/sfs%psa(j,ici1)
            vavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
            if ( vavg >= d_zero ) then
              qvx = qext
            else
              qvx = qint
            end if
            atm1%qv(j,ice1,k) = qvx*sfs%psa(j,ice1)
          end do
        end do
      end if
      !
      ! north boundary:
      !
      if ( ma%hastop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            qext = atm1%qv(j,ice2,k)/sfs%psa(j,ice2)
            qint = atm1%qv(j,ici2,k)/sfs%psa(j,ici2)
            vavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
            if ( vavg < d_zero ) then
              qvx = qext
            else
              qvx = qint
            end if
            atm1%qv(j,ice2,k) = qvx*sfs%psa(j,ice2)
          end do
        end do
      end if
    end if
!
    ! set boundary values for p*qc and p*qr:
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
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = ici1 , ici2
          qcint = atm1%qc(jci1,i,k)/sfs%psa(jci1,i)
          uavg = wue(i,k) + wue(i+1,k) + wui(i,k) + wui(i+1,k)
          if ( uavg >= d_zero ) then
            qcx = d_zero
          else
            qcx = qcint
          end if
          atm1%qc(jce1,i,k) = qcx*sfs%psa(jce1,i)
        end do
      end do
    end if
    !
    ! east boundary:
    !
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = ici1 , ici2
          qcint = atm1%qc(jci2,i,k)/sfs%psa(jci2,i)
          uavg = eue(i,k) + eue(i+1,k) + eui(i,k) + eui(i+1,k)
          if ( uavg < d_zero ) then
            qcx = d_zero
          else
            qcx = qcint
          end if
          atm1%qc(jce2,i,k) = qcx*sfs%psa(jce2,i)
        end do
      end do
    end if
    !
    ! south boundary:
    !
    if ( ma%hasbottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          qcint = atm1%qc(j,ici1,k)/sfs%psa(j,ici1)
          vavg = sve(j,k) + sve(j+1,k) + svi(j,k) + svi(j+1,k)
          if ( vavg >= d_zero ) then
            qcx = d_zero
          else
            qcx = qcint
          end if
          atm1%qc(j,ice1,k) = qcx*sfs%psa(j,ice1)
        end do
      end do
    end if
    !
    ! north boundary:
    !
    if ( ma%hastop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          qcint = atm1%qc(j,ici2,k)/sfs%psa(j,ici2)
          vavg = nve(j,k) + nve(j+1,k) + nvi(j,k) + nvi(j+1,k)
          if ( vavg < d_zero ) then
            qcx = d_zero
          else
            qcx = qcint
          end if
          atm1%qc(j,ice2,k) = qcx*sfs%psa(j,ice2)
        end do
      end do
    end if
   
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call set_tke_bc(atm1,atm2)
    end if
!
    call time_end(subroutine_name,idindx)
!
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
  subroutine sponge3d(nk,ba,bnd,ften)
!
    implicit none
!
    integer , intent(in) :: nk
    type(bound_area) , intent(in) :: ba
    type(v3dbound) , intent(in) :: bnd
    real(dp) , pointer , intent(inout) , dimension(:,:,:) :: ften
!
    integer :: i , j , k
    integer :: ib , i1 , i2 , j1 , j2
    real(dp) , pointer , dimension(:) :: wg
    character (len=64) :: subroutine_name='sponge3d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
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
!
!----------------------------------------------------------------------
!
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
    call time_end(subroutine_name,idindx)
  end subroutine sponge3d
!
  subroutine sponge2d(ba,bnd,ften)
!
    implicit none
!
    type(bound_area) , intent(in) :: ba
    type(v2dbound) , intent(in) :: bnd
    real(dp) , pointer , intent(inout) , dimension(:,:) :: ften
!
    integer :: i , j , i1 , i2 , j1 , j2
    integer :: ib
    real(dp) , pointer , dimension(:) :: wg
    character (len=64) :: subroutine_name='sponge2d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
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
!
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
    call time_end(subroutine_name,idindx)

  end subroutine sponge2d
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
  function xfun(mm,ldot)
    implicit none
    real(dp) :: xfun
    integer , intent(in) :: mm
    logical , intent(in) :: ldot
    if ( ldot ) then
      xfun = dble(nspgd-mm)/dble(nspgd-2)
    else
      xfun = dble(nspgx-mm)/dble(nspgx-2)
    end if
  end function xfun
!
  function xfune(mm,kk)
    implicit none
    real(dp) :: xfune
    integer , intent(in) :: mm , kk
    xfune = dexp(-dble(mm-2)/anudg(kk))
  end function xfune
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
! j     : is the j'th slice of the tendency to be adjusted.       c
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
  subroutine nudge3d(nk,ba,xt,f,ibdy,bnd,ften)
!
    implicit none
!
    integer , intent(in) :: ibdy , nk
    real(dp) , intent(in) :: xt
    real(dp) , pointer , intent(in) , dimension(:,:,:) :: f
    type(v3dbound) , intent(in) :: bnd
    type(bound_area) , intent(in) :: ba
    real(dp) , pointer , intent(inout) , dimension(:,:,:) :: ften
!
    real(dp) :: xf , fls0 , fls1 , fls2 , fls3 , fls4 , xg
    integer :: i , j , k , ib , i1 , i2 , j1 , j2
    character (len=64) :: subroutine_name='nudge3d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
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
!
    if ( ba%ns /= 0 ) then
      do k = 1 , nk
        do i = i1 , i2
          do j = j1 , j2
            if ( .not. ba%bsouth(j,i) ) cycle
            ib = ba%ibnd(j,i)
            if ( ibdy == 1 ) then
              xf = lfc(ib)
              xg = lgc(ib)
            else if ( ibdy == 5 ) then
              xf = efc(ib,k)
              xg = egc(ib,k)
            end if
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
            if ( ibdy == 1 ) then
              xf = lfc(ib)
              xg = lgc(ib)
            else if ( ibdy == 5 ) then
              xf = efc(ib,k)
              xg = egc(ib,k)
            end if
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
            if ( ibdy == 1 ) then
              xf = lfc(ib)
              xg = lgc(ib)
            else if ( ibdy == 5 ) then
              xf = efc(ib,k)
              xg = egc(ib,k)
            end if
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
            if ( ibdy == 1 ) then
              xf = lfc(ib)
              xg = lgc(ib)
            else if ( ibdy == 5 ) then
              xf = efc(ib,k)
              xg = egc(ib,k)
            end if
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

    call time_end(subroutine_name,idindx)
  end subroutine nudge3d
!
! ###################################################################
!
  subroutine nudge2d(ba,xt,f,ibdy,bnd,ften)
!
    implicit none
!
    integer , intent(in) :: ibdy
    real(dp) , intent(in) :: xt
    real(dp) , pointer , intent(in) , dimension(:,:) :: f
    type(v2dbound) , intent(in) :: bnd
    type(bound_area) , intent(in) :: ba
    real(dp) , pointer , intent(inout) , dimension(:,:) :: ften
!
    real(dp) :: xf , fls0 , fls1 , fls2 , fls3 , fls4 , xg
    integer :: i , j , ib , i1 , i2 , j1 , j2
    character (len=64) :: subroutine_name='nudge2d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
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

    if ( ba%ns /= 0 ) then
      do i = i1 , i2
        do j = j1 , j2
          if ( .not. ba%bsouth(j,i) ) cycle
          ib = ba%ibnd(j,i)
          if ( ibdy == 1 ) then
            xf = lfc(ib)
            xg = lgc(ib)
          else if ( ibdy == 5 ) then
            xf = efc(ib,kz)
            xg = egc(ib,kz)
          end if
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
          if ( ibdy == 1 ) then
            xf = lfc(ib)
            xg = lgc(ib)
          else if ( ibdy == 5 ) then
            xf = efc(ib,kz)
            xg = egc(ib,kz)
          end if
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
          if ( ibdy == 1 ) then
            xf = lfc(ib)
            xg = lgc(ib)
          else if ( ibdy == 5 ) then
            xf = efc(ib,kz)
            xg = egc(ib,kz)
          end if
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
          if ( ibdy == 1 ) then
            xf = lfc(ib)
            xg = lgc(ib)
          else if ( ibdy == 5 ) then
            xf = efc(ib,kz)
            xg = egc(ib,kz)
          end if
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

    call time_end(subroutine_name,idindx)

  end subroutine nudge2d
!
end module mod_bdycod
