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
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_moloch

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_stdio
  use mod_service
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_micro_interface
  use mod_bdycod
  use mod_slice
  use mod_sun
  use mod_slabocean
  use mod_massck
  use mod_stdatm
  use mod_zita

  implicit none

  private

  ! generalized vertical velocity
  real(rkx) , pointer , dimension(:,:,:) :: s
  ! nonhydrostatic term in pressure gradient force
  real(rkx) , pointer , dimension(:,:,:) :: deltaw
  ! tridiagonal inversion
  real(rkx) , pointer , dimension(:,:,:) :: wwkw
  real(rkx) , pointer , dimension(:,:,:) :: wx
  real(rkx) , pointer , dimension(:,:,:) :: tkex
  real(rkx) , pointer , dimension(:,:,:) :: wz
  real(rkx) , pointer , dimension(:,:) :: wfw
  real(rkx) , pointer , dimension(:,:) :: mx2
  real(rkx) , pointer , dimension(:,:) :: rmu
  real(rkx) , pointer , dimension(:,:) :: rmv
  real(rkx) , pointer , dimension(:,:,:) :: p0
  real(rkx) , pointer , dimension(:,:,:) :: zdiv2
  real(rkx) , pointer , dimension(:,:) :: zpby
  real(rkx) , pointer , dimension(:,:) :: zpbw

  real(rkx) , pointer , dimension(:,:,:) :: ten0
  real(rkx) , pointer , dimension(:,:,:) :: qen0
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0

  real(rkx) , dimension(:) , pointer :: gzitak
  real(rkx) , dimension(:) , pointer :: gzitakh
  real(rkx) , dimension(:,:) , pointer :: p2d
  real(rkx) , dimension(:,:) , pointer :: xlat , xlon , coru , corv
  real(rkx) , dimension(:,:) , pointer :: mu , hx , mx
  real(rkx) , dimension(:,:) , pointer :: mv , hy
  real(rkx) , dimension(:,:) , pointer :: ps , ht
  real(rkx) , dimension(:,:,:) , pointer :: fmz
  real(rkx) , dimension(:,:,:) , pointer :: fmzf
  real(rkx) , dimension(:,:,:) , pointer :: pai , pf
  real(rkx) , dimension(:,:,:) , pointer :: tetav , tf , tvirt
  real(rkx) , dimension(:,:,:) , pointer :: zeta , zetau , zetav
  real(rkx) , dimension(:,:,:) , pointer :: u , v , w
  real(rkx) , dimension(:,:,:) , pointer :: ux , vx
  real(rkx) , dimension(:,:,:) , pointer :: ud , vd
  real(rkx) , dimension(:,:,:) , pointer :: p , t , rho
  real(rkx) , dimension(:,:,:) , pointer :: qv , qc , qi , qr , qs , qsat
  real(rkx) , dimension(:,:,:) , pointer :: qwltot , qwitot
  real(rkx) , dimension(:,:,:) , pointer :: tke
  real(rkx) , dimension(:,:,:,:) , pointer :: qx , trac

  public :: allocate_moloch , init_moloch , moloch
  public :: uvstagtox , xtouvstag , wstagtox

  real(rkx) , parameter :: minden = 1.0e-30_rkx

  logical , parameter :: do_phys = .true.
  logical , parameter :: do_bdy = .true.
  logical , parameter :: do_fulleq = .true.
  logical , parameter :: do_filterpai = .false.
  logical , parameter :: do_vadvtwice = .true.
  logical , parameter :: do_filterqx = .false.
  logical , parameter :: do_filterdiv = .false.
  logical , parameter :: do_filtertheta = .false.
  logical :: moloch_realcase = (.not. moloch_do_test_1) .and. &
                               (.not. moloch_do_test_2)
  logical :: lrotllr

  integer :: nadv , nsound
  real(rkx) :: nupait , ddamp , dzita

  contains

#include <wlh.inc>
#include <pfesat.inc>
#include <pfwsat.inc>

  subroutine allocate_moloch
    implicit none
    call getmem1d(gzitak,1,kzp1,'moloch:gzitak')
    call getmem1d(gzitakh,1,kz,'moloch:gzitakh')
    call getmem2d(p2d,jdi1,jdi2,idi1,idi2,'moloch:p2d')
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:deltaw')
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(wx,jce1,jce2,ice1,ice2,1,kz,'moloch:wx')
    call getmem3d(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
    call getmem3d(wwkw,jci1,jci2,ici1,ici2,1,kzp1,'moloch:wwkw')
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem2d(wfw,jci1,jci2,1,kzp1,'moloch:wfw')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
    call getmem2d(zpby,jci1,jci2,ici1,ice2ga,'moloch:zpby')
    call getmem2d(zpbw,jci1,jce2ga,ici1,ici2,'moloch:zpbw')
    call getmem2d(mx2,jde1,jde2,ide1,ide2,'moloch:mx2')
    call getmem2d(rmu,jde1ga,jde2ga,ide1,ide2,'moloch:rmu')
    call getmem2d(rmv,jde1,jde2,ide1ga,ide2ga,'moloch:rmv')
    call getmem2d(coru,jde1,jde2,ice1,ice2,'moloch:coru')
    call getmem2d(corv,jce1,jce2,ide1,ide2,'moloch:corv')
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jce1,jce2,ice1,ice2,1,kz,'moloch:tkex')
    end if
    if ( idiag > 0 ) then
      call getmem3d(ten0,jci1,jci2,ici1,ici2,1,kz,'moloch:ten0')
      call getmem3d(qen0,jci1,jci2,ici1,ici2,1,kz,'moloch:qen0')
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        call getmem4d(chiten0,jci1,jci2,ici1,ici2,1,kz,1,ntr,'moloch:chiten0')
      end if
    end if
    call getmem3d(ud,jde1ga,jde2ga,ice1ga,ice2ga,1,kz,'moloch:ud')
    call getmem3d(vd,jce1ga,jce2ga,ide1ga,ide2ga,1,kz,'moloch:vd')
    if ( ifrayd == 1 ) then
      call getmem3d(zetau,jdi1,jdi2,ici1,ici2,1,kz,'moloch:zetau')
      call getmem3d(zetav,jci1,jci2,idi1,idi2,1,kz,'moloch:zetav')
    end if
    if ( do_fulleq ) then
      if ( ipptls > 0 ) then
        call getmem3d(qwltot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwltot')
        call getmem3d(qwitot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwitot')
      end if
    end if
    if ( do_filterpai ) then
      call getmem3d(pf,jce1,jce2,ice1,ice2,1,kz,'moloch:pf')
    end if
    if ( do_filtertheta ) then
      call getmem3d(tf,jce1,jce2,ice1,ice2,1,kz,'moloch:tf')
    end if
  end subroutine allocate_moloch

  subroutine init_moloch
    implicit none
    call assignpnt(mddom%msfu,mu)
    call assignpnt(mddom%msfv,mv)
    call assignpnt(mddom%msfx,mx)
    call assignpnt(mddom%hx,hx)
    call assignpnt(mddom%hy,hy)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(mddom%xlon,xlon)
    call assignpnt(mddom%ht,ht)
    call assignpnt(sfs%psa,ps)
    call assignpnt(mo_atm%fmz,fmz)
    call assignpnt(mo_atm%fmzf,fmzf)
    call assignpnt(mo_atm%pai,pai)
    call assignpnt(mo_atm%tetav,tetav)
    call assignpnt(mo_atm%u,u)
    call assignpnt(mo_atm%ux,ux)
    call assignpnt(mo_atm%v,v)
    call assignpnt(mo_atm%vx,vx)
    call assignpnt(mo_atm%w,w)
    call assignpnt(mo_atm%tvirt,tvirt)
    call assignpnt(mo_atm%zeta,zeta)
    call assignpnt(mo_atm%p,p)
    call assignpnt(mo_atm%t,t)
    call assignpnt(mo_atm%rho,rho)
    call assignpnt(mo_atm%qx,qx)
    call assignpnt(mo_atm%qs,qsat)
    call assignpnt(mo_atm%qx,qv,iqv)
    if ( ipptls > 0 ) then
      call assignpnt(mo_atm%qx,qc,iqc)
      if ( ipptls > 1 ) then
        call assignpnt(mo_atm%qx,qi,iqi)
        call assignpnt(mo_atm%qx,qr,iqr)
        call assignpnt(mo_atm%qx,qs,iqs)
      end if
    end if
    if ( ibltyp == 2 ) call assignpnt(mo_atm%tke,tke)
    if ( ichem == 1 ) call assignpnt(mo_atm%trac,trac)
    if ( ifrayd == 1 ) then
      call xtoustag(zeta,zetau)
      call xtovstag(zeta,zetav)
    end if
    coru = eomeg2*sin(mddom%ulat(jde1:jde2,ice1:ice2)*degrad)
    corv = eomeg2*sin(mddom%vlat(jce1:jce2,ide1:ide2)*degrad)
    mx2 = mx * mx
    rmu = d_one/mu
    rmv = d_one/mv
    gzitak = gzita(zita)
    gzitakh = gzita(zitah)
    nadv = mo_nadv
    nsound = mo_nsound
    dzita = mo_dzita
    nupait = 0.05_rkx
    wwkw(:,:,kzp1) = d_zero
    w(:,:,1) = d_zero
    lrotllr = (iproj == 'ROTLLR')
    ddamp = 0.25_rkx
  end subroutine init_moloch
  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch
    implicit none
    integer(ik4) :: jadv , jsound
    real(rkx) :: dtsound , dtstepa
    real(rkx) :: maxps , minps , pmax , pmin , zdgz
    real(rkx) :: tv , lrt , fice
    !real(rk8) :: jday
    integer(ik4) :: i , j , k
    integer(ik4) :: iconvec
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    dtstepa = dtsec / real(nadv,rkx)
    dtsound = dtstepa / real(nsound,rkx)
    iconvec = 0

    call reset_tendencies

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          p(j,i,k) = (pai(j,i,k)**cpovr) * p00
          qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
        end do
      end do
    end do

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
                                         qc(j,i,k) - qi(j,i,k) - &
                                         qr(j,i,k) - qs(j,i,k))
            end do
          end do
        end do
        if ( do_fulleq ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                qwltot(j,i,k) = qc(j,i,k) + qr(j,i,k)
                qwitot(j,i,k) = qi(j,i,k) + qs(j,i,k)
              end do
            end do
          end do
        end if
      else
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
        if ( do_fulleq ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( t(j,i,k) >= tzero ) then
                  qwltot(j,i,k) = qc(j,i,k)
                else if ( t(j,i,k) <= -20.0_rkx+tzero ) then
                  qwitot(j,i,k) = qc(j,i,k)
                else
                  fice = (tzero-t(j,i,k))/20.0_rkx
                  qwltot(j,i,k) = qc(j,i,k) * (1.0_rkx-fice)
                  qwitot(j,i,k) = qc(j,i,k) * fice
                end if
              end do
            end do
          end do
        end if
      end if
    else
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
      if ( do_fulleq ) then
        qwltot(:,:,:) = d_zero
        qwitot(:,:,:) = d_zero
      end if
    end if

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
        end do
      end do
    end do

    if ( idiag > 0 ) then
      ten0 = t(jci1:jci2,ici1:ici2,:)
      qen0 = qv(jci1:jci2,ici1:ici2,:)
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
      end if
    end if

    if ( do_filterpai ) then
      pf = pai(jce1:jce2,ice1:ice2,:)
    end if
    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        tf = tetav(jce1:jce2,ice1:ice2,:)
      end if
    end if

    do jadv = 1 , nadv

      call sound(dtsound)

      call advection(dtstepa)

    end do ! Advection loop

    if ( do_filterpai ) then
      pai(jce1:jce2,ice1:ice2,:) = pai(jce1:jce2,ice1:ice2,:) - pf
      call filtpai
      pai(jce1:jce2,ice1:ice2,:) = pai(jce1:jce2,ice1:ice2,:) + pf
    end if
    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        tetav(jce1:jce2,ice1:ice2,:) = tetav(jce1:jce2,ice1:ice2,:) - tf
        call filttheta
        tetav(jce1:jce2,ice1:ice2,:) = tetav(jce1:jce2,ice1:ice2,:) + tf
      end if
    end if

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
        end do
      end do
    end do

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
                             qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
            end do
          end do
        end do
      else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
      end if
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
    end if

    if ( idiag > 0 ) then
      tdiag%adh = (t(jci1:jci2,ici1:ici2,:) - ten0) * rdt
      qdiag%adh = (qv(jci1:jci2,ici1:ici2,:) - qen0) * rdt
    end if

    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        cadvhdiag = (trac(jci1:jci2,ici1:ici2,:,:) - chiten0) * rdt
      end if
    end if
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          p(j,i,k) = (pai(j,i,k)**cpovr) * p00
          rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
        end do
      end do
    end do

    !jday = yeardayfrac(rcmtimer%idate)
    do i = ice1 , ice2
      do j = jce1 , jce2
        zdgz = zeta(j,i,kz)*egrav
        lrt = (tvirt(j,i,kz-1)-tvirt(j,i,kz))/(zeta(j,i,kz-1)-zeta(j,i,kz))
        ! lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,xlat(j,i))
        lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
        tv = tvirt(j,i,kz) - d_half*zeta(j,i,kz)*lrt
        ps(j,i) = p(j,i,kz) * exp(zdgz/(rgas*tv))
      end do
    end do
    !
    ! Recompute saturation
    !
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
        end do
      end do
    end do
    !
    ! Lateral/damping boundary condition
    !
    if ( do_bdy .and. moloch_realcase .and. irceideal == 0 ) then
      call boundary
      if ( i_crm /= 1 ) then
        if ( ifrayd == 1 ) then
          call raydamp(zetau,u,xub,jdi1,jdi2,ici1,ici2,1,kz)
          call raydamp(zetav,v,xvb,jci1,jci2,idi1,idi2,1,kz)
          call raydamp(zeta,t,xtb,jci1,jci2,ici1,ici2,1,kz)
          call raydamp(zeta,pai,xpaib,jci1,jci2,ici1,ici2,1,kz)
        end if
      end if
    else
      if ( debug_level > 1 ) then
        if ( myid == italk .and. irceideal == 0 ) then
          write(stdout,*) 'WARNING: Physical boundary package disabled!!!'
        end if
      end if
    end if
    !
    ! Prepare fields to be used in physical parametrizations.
    !
    call mkslice
    !
    ! PHYSICS
    !
    if ( do_phys .and. moloch_realcase ) then
      call physical_parametrizations
    else
      if ( debug_level > 1 ) then
        if ( myid == italk ) then
          write(stdout,*) 'WARNING: Physical package disabled!!!'
        end if
      end if
    end if
    !
    ! Mass check
    !
    if ( debug_level > 0 ) call massck
    !
    ! Diagnostic and end timestep
    !
    if ( syncro_rep%act( ) .and. rcmtimer%integrating( ) ) then
      maxps = maxval(ps(jci1:jci2,ici1:ici2))
      minps = minval(ps(jci1:jci2,ici1:ici2))
      call maxall(maxps,pmax)
      call minall(minps,pmin)
      call sumall(total_precip_points,iconvec)
      if ( is_nan(pmax) .or. is_nan(pmin) ) then
        write (stderr,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
        write (stderr,*) 'No more atmosphere here....'
        write (stderr,*) 'CFL violation detected, so model STOP'
        write (stderr,*) '#####################################'
        write (stderr,*) '#            DECREASE DT !!!!       #'
        write (stderr,*) '#####################################'
        call fatal(__FILE__,__LINE__,'CFL VIOLATION')
      end if
      if ( myid == 0 ) then
        write(stdout,*) '$$$ ', rcmtimer%str( )
        write(stdout,'(a,2f8.2)') &
            ' $$$ max, min of ps (mb) = ', pmax*d_r100 , pmin*d_r100
        if ( any(icup > 0) ) then
          write(stdout,'(a,i7)') &
            ' $$$ no. of points with active convection = ', iconvec
        end if
      end if
    end if
    !
    ! Next timestep ready : increment elapsed forecast time
    !
    call rcmtimer%advance( )
    if ( islab_ocean == 1 ) xslabtime = xslabtime + dtsec
    !
    ! calculate new solar zenith angle
    !
    call zenitm(xlat,xlon,coszrs)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      subroutine boundary
        implicit none

        call exchange_lrbt(ps,1,jce1,jce2,ice1,ice2)
        call exchange_lrbt(u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange_lrbt(v,1,jce1,jce2,ide1,ide2,1,kz)
        call exchange_lrbt(t,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(qv,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
        if ( ichem == 1 ) then
          call exchange_lrbt(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if

        if ( idiag > 0 ) then
          ten0 = t(jci1:jci2,ici1:ici2,:)
          qen0 = qv(jci1:jci2,ici1:ici2,:)
        end if
        if ( ichem == 1 ) then
          if ( ichdiag > 0 ) then
            chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
          end if
        end if

        if ( iboudy == 1 .or. iboudy >= 5 ) then
          call nudge(iboudy,pai,xpaib)
          call nudge(iboudy,ps,xpsb)
          call nudge(iboudy,t,xtb)
          call nudge(iboudy,qv,xqb)
          call nudge(iboudy,u,v,xub,xvb)

          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
          if ( is_present_qc( ) ) then
            call nudge(iboudy,qc,xlb)
          end if
          if ( is_present_qi( ) ) then
            call nudge(iboudy,qi,xib)
          end if
        else if ( iboudy == 4 ) then
          call sponge(ps,xpsb)
          call sponge(pai,xpaib)
          call sponge(t,xtb)
          call sponge(qv,xqb)
          call sponge(u,v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
          if ( is_present_qc( ) ) then
            call sponge(qc,xlb)
          end if
          if ( is_present_qi( ) ) then
            call sponge(qi,xib)
          end if
        end if
        if ( ichem == 1 ) then
          if ( iboudy == 1 .or. iboudy >= 5 ) then
            call nudge_chi(trac)
          else if ( iboudy == 4 ) then
            ! Not implemented sponge_chi
          end if
          if ( ichdiag > 0 ) then
            cbdydiag = trac(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        end if
        call uvstagtox(u,v,ux,vx)
        if ( do_filterqx .and. ipptls > 0 ) then
          call exchange_lrbt(qx,1,jce1,jce2,ice1,ice2,1,kz,iqfrst,iqlst)
          call filt4d(qx,mo_anu2,iqfrst,iqlst)
        end if
      end subroutine boundary

      subroutine filt4d(p,nu,n1,n2)
        implicit none
        real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: p
        real(rkx) , intent(in) :: nu
        integer(ik4) , intent(in) :: n1 , n2
        integer(ik4) :: j , i , k , n

        do n = n1 , n2
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                p2d(j,i) = 0.125_rkx * (p(j-1,i,k,n) + p(j+1,i,k,n) + &
                                        p(j,i-1,k,n) + p(j,i+1,k,n)) - &
                           d_half   * p(j,i,k,n)
              end do
            end do
            do i = ici1 , ici2
              do j = jci1 , jci2
                p(j,i,k,n) = p(j,i,k,n) + nu * p2d(j,i)
              end do
            end do
          end do
        end do
      end subroutine filt4d

      subroutine filt3d(p,nu)
        implicit none
        real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: p
        real(rkx) , intent(in) :: nu
        integer(ik4) :: j , i , k

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (p(j-1,i,k) + p(j+1,i,k) + &
                                      p(j,i-1,k) + p(j,i+1,k)) - &
                         d_half   * p(j,i,k)
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              p(j,i,k) = p(j,i,k) + nu * p2d(j,i)
            end do
          end do
        end do
      end subroutine filt3d

      subroutine filtpai
        implicit none
        integer(ik4) :: j , i , k

        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (pai(j-1,i,k) + pai(j+1,i,k) + &
                                      pai(j,i-1,k) + pai(j,i+1,k)) - &
                         d_half   * pai(j,i,k)
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              pai(j,i,k) = pai(j,i,k) + nupait * p2d(j,i)
            end do
          end do
        end do
      end subroutine filtpai

      subroutine filttheta
        implicit none
        integer(ik4) :: j , i , k

        call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (tetav(j-1,i,k) + tetav(j+1,i,k) + &
                                      tetav(j,i-1,k) + tetav(j,i+1,k)) - &
                         d_half   * tetav(j,i,k)
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              tetav(j,i,k) = tetav(j,i,k) + nupait * p2d(j,i)
            end do
          end do
        end do
      end subroutine filttheta

      subroutine sound(dts)
        implicit none
        real(rkx) , intent(in) :: dts
        integer(ik4) :: i , j , k
        real(rkx) :: zuh , zvh , zcx , zcy
        real(rkx) :: zrfmzum , zrfmzup , zrfmzvm , zrfmzvp
        real(rkx) :: zup , zum , zvp , zvm , zqs , zdth
        real(rkx) :: zrom1w , zwexpl , zu , zd , zrapp
        real(rkx) :: zfz , zcor1u , zcor1v
        real(rkx) :: zrom1u , zrom1v
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz , zcs2

        zdtrdx = dts/dx
        zdtrdy = dts/dx
        zdtrdz = dts/dzita
        zcs2 = zdtrdz**2*rdrcv

        !  sound waves

        if ( .not. do_fulleq ) then
          call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
        end if

        do jsound = 1 , nsound

          ! partial definition of the generalized vertical velocity

          call exchange(u,1,jde1,jde2,ice1,ice2,1,kz)
          call exchange(v,1,jce1,jce2,ide1,ide2,1,kz)

          ud(:,:,1:kz) = u(jde1ga:jde2ga,ice1ga:ice2ga,1:kz)
          vd(:,:,1:kz) = v(jce1ga:jce2ga,ide1ga:ide2ga,1:kz)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
              zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
              w(j,i,kzp1) = d_half * (zuh+zvh)
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,kzp1) = -w(j,i,kzp1)
            end do
          end do

          ! Equation 10, generalized vertical velocity

          do k = kz , 2 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
                zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) + &
                      (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
                zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) + &
                      (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
                s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
              end do
            end do
          end do

          ! Part of divergence (except w contribution) put in zdiv2
          ! Equation 16

          if ( lrotllr ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
                  zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
                  zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
                  zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
                  zum = zdtrdx * u(j,i,k) * zrfmzum
                  zup = zdtrdx * u(j+1,i,k) * zrfmzup
                  zvm = zdtrdy * v(j,i,k) * zrfmzvm * rmv(j,i)
                  zvp = zdtrdy * v(j,i+1,k) * zrfmzvp * rmv(j,i+1)
                  zdiv2(j,i,k) = fmz(j,i,k) * mx(j,i) * ((zup-zum) + (zvp-zvm))
                end do
              end do
            end do
          else
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
                  zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
                  zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
                  zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
                  zum = zdtrdx * u(j,i,k)   * rmu(j,i)   * zrfmzum
                  zup = zdtrdx * u(j+1,i,k) * rmu(j+1,i) * zrfmzup
                  zvm = zdtrdy * v(j,i,k)   * rmv(j,i)   * zrfmzvm
                  zvp = zdtrdy * v(j,i+1,k) * rmv(j,i+1) * zrfmzvp
                  zdiv2(j,i,k) = mx2(j,i) * fmz(j,i,k) * ((zup-zum)+(zvp-zvm))
                end do
              end do
            end do
          end if
          call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
          call divdamp(dtsound)
          if ( do_filterdiv ) call filt3d(zdiv2,mo_anu2)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv2(j,i,k) = zdiv2(j,i,k) + fmz(j,i,k) * &
                       zdtrdz * (s(j,i,k) - s(j,i,k+1))
              end do
            end do
          end do

          ! new w (implicit scheme) from Equation 19

          do k = kz , 2 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
                deltaw(j,i,k) = -w(j,i,k)
                ! explicit w:
                !    it must be consistent with the initialization of pai
                zrom1w = d_half * cpd * fmzf(j,i,k) * &
                        (tetav(j,i,k-1) + tetav(j,i,k))
                zrom1w = zrom1w - cpd * w(j,i,k) * &
                         fmzf(j,i,k)*fmzf(j,i,k) * &
                         real(jsound,rkx) * zdtrdz * &
                         (tetav(j,i,k-1) - tetav(j,i,k)) !! GW
                if ( qv(j,i,k) > 0.96_rkx*qsat(j,i,k) .and. &
                     w(j,i,k) > 0.1_rkx ) then
                  zqs = d_half*(qsat(j,i,k)+qsat(j,i,k-1))
                  zdth = egrav*w(j,i,k)*real(jsound-1,rkx)*dts*wlhv*wlhv* &
                    zqs/(cpd*pai(j,i,k-1)*rwat*t(j,i,k-1)*t(j,i,k-1))
                  zrom1w = zrom1w + zdth*fmzf(j,i,k)
                end if
                zwexpl = w(j,i,k) - zrom1w * zdtrdz * &
                         (pai(j,i,k-1) - pai(j,i,k)) - egrav*dts
                zwexpl = zwexpl + rdrcv * zrom1w * zdtrdz * &
                         (pai(j,i,k-1) * zdiv2(j,i,k-1) - &
                          pai(j,i,k)   * zdiv2(j,i,k))
                ! computation of the tridiagonal matrix coefficients
                ! -zu*w(k+1) + (1+zu+zd)*w(k) - zd*w(k-1) = zwexpl
                zu = zcs2 * fmz(j,i,k-1) * zrom1w * pai(j,i,k-1) + ffilt(k)
                zd = zcs2 * fmz(j,i,k)   * zrom1w * pai(j,i,k)   + ffilt(k)
                ! 1st loop for the tridiagonal inversion
                ! a = -zd ; b = (1+zu+zd) ; c = -zu
                zrapp = d_one / (d_one + zd + zu - zd*wwkw(j,i,k+1))
                w(j,i,k) = zrapp * (zwexpl + zd * w(j,i,k+1))
                wwkw(j,i,k) = zrapp * zu
              end do
            end do
          end do

          ! 2nd loop for the tridiagonal inversion
          do k = 2 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
                deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
              end do
            end do
          end do

          ! new Exner function (Equation 19)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv2(j,i,k) = zdiv2(j,i,k) + zdtrdz * fmz(j,i,k) * &
                      (w(j,i,k) - w(j,i,k+1))
              end do
            end do
          end do

          if ( do_fulleq ) then
            if ( ipptls > 0 ) then
              do k = 1 , kz
                do i = ici1 , ici2
                  do j = jci1 , jci2
                    zdiv2(j,i,k) = zdiv2(j,i,k) * &
                       (d_one + 0.86_rkx * qv(j,i,k) + &
                                3.2_rkx * qc(j,i,k)) / &
                       (d_one + 0.96_rkx * qv(j,i,k) + &
                                4.8_rkx * qc(j,i,k))
                    tetav(j,i,k) = tetav(j,i,k) * &
                       (d_one + rdrcv*zdiv2(j,i,k) * &
                        (0.25_rkx * qv(j,i,k) +      &
                         4.2_rkx * qwltot(j,i,k) +   &
                         2.1_rkx * qwitot(j,i,k)))
                  end do
                end do
              end do
              call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
            end if
          end if

          ! horizontal momentum equations

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv2(j,i,k))
              end do
            end do
          end do

          call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
          call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1)

          if ( lrotllr ) then

            do k = 1 , kz
              do i = ici1 , ici2
                do j = jdi1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
            zcy = zdtrdy
            do k = 1 , kz
              do i = idi1 , idi2
                do j = jci1 , jci2
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do

          else

            do k = 1 , kz
              do i = ici1 , ici2
                do j = jdi1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
            do k = 1 , kz
              do i = idi1 , idi2
                do j = jci1 , jci2
                  zcy = zdtrdy * mv(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do
          end if

        end do ! sound loop

        ! complete computation of generalized vertical velocity
        ! Complete Equation 10

        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,1) = d_zero
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,kzp1) = d_zero
          end do
        end do

      end subroutine sound

      subroutine advection(dta)
        implicit none
        integer(ik4) :: n
        real(rkx) :: dta
        real(rkx) , pointer , dimension(:,:,:) :: ptr

        call uvstagtox(u,v,ux,vx)

        ! Compute W (and TKE if required) on zita levels

        call wstagtox(w,wx)

        if ( ibltyp == 2 ) then
          call wstagtox(tke,tkex)
        end if

        call wafone(tetav,dta)
        call wafone(pai,dta)
        call wafone(ux,dta)
        call wafone(vx,dta)
        call wafone(wx,dta)
        call wafone(qv,dta)
        if ( ipptls > 0 ) then
          do n = iqfrst , iqlst
            call assignpnt(qx,ptr,n)
            call wafone(ptr,dta)
          end do
        end if
        if ( ibltyp == 2 ) then
          call wafone(tkex,dta)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(trac,ptr,n)
            call wafone(ptr,dta)
          end do
        end if

        ! Interpolate on staggered points
        call xtouvstag(ux,vx,u,v)

        ! Back to half-levels
        call xtowstag(wx,w)
        if ( ibltyp == 2 ) then
          call xtowstag(tkex,tke)
        end if
      end subroutine advection

      pure real(rkx) function rdeno(t1,t2,t3,t4)
        implicit none
        real(rkx) , intent(in) :: t1 , t2 , t3 , t4
        real(rkx) :: zzden
        zzden = (t3-t4)
        rdeno = (t1-t2)/sign(max(abs(zzden),minden),zzden)
      end function rdeno

      subroutine wafone(pp,dta)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
        real(rkx) , intent(in) :: dta
        integer(ik4) :: j , i , k
        integer(ik4) :: k1 , k1p1 , ih , ihm1 , jh , jhm1
        real(rkx) :: zamu , r , b , zphi , is , zdv , zrfmu , zrfmd
        real(rkx) :: zrfmn , zrfmw , zrfme , zrfms
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz
        real(rkx) :: zhxvtn , zhxvts , zcostx

        real(rkx) , parameter :: wlow  = 0.0_rkx
        real(rkx) , parameter :: whigh = 2.0_rkx

        zdtrdx = dta/dx
        zdtrdy = dta/dx
        zdtrdz = dta/dzita
        if ( do_vadvtwice ) then
          zdtrdz = 0.5_rkx * zdtrdz
        end if

        ! Vertical advection

        do j = jci1 , jci2
          wfw(j,1) = d_zero
          wfw(j,kzp1) = d_zero
        end do

        do i = ici1 , ici2
          do k = 1 , kzm1
            do j = jci1 , jci2
              zamu = s(j,i,k+1) * zdtrdz
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k + 1
                k1p1 = k1 + 1
                if ( k1p1 > kz ) k1p1 = kz
              else
                is = -d_one
                k1 = k - 1
                k1p1 = k
                if ( k1 < 1 ) k1 = 1
              end if
              r = rdeno(pp(j,i,k1),pp(j,i,k1p1),pp(j,i,k),pp(j,i,k+1))
              b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
              zphi = is + zamu * b - is * b
              wfw(j,k+1) = d_half * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                                  (d_one-zphi)*pp(j,i,k))
              !wfw(j,k+1) = d_half * zamu * ((d_one+zphi)*pp(j,i,k+1) + &
              !                              (d_one-zphi)*pp(j,i,k))
            end do
          end do
          do k = 1 , kz
            do j = jci1 , jci2
              zrfmu = zdtrdz * fmz(j,i,k)/fmzf(j,i,k)
              zrfmd = zdtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
              zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
              wz(j,i,k) = pp(j,i,k) - wfw(j,k)*zrfmu + wfw(j,k+1)*zrfmd + zdv
            end do
          end do
          !do k = 1 , kz
          !  do j = jci1 , jci2
          !    zdv = (s(j,i,k) - s(j,i,k+1)) * zdtrdz * pp(j,i,k)
          !    wz(j,i,k) = pp(j,i,k) - wfw(j,k) + wfw(j,k+1) + zdv
          !  end do
          !end do
        end do

        if ( do_vadvtwice ) then
          do i = ici1 , ici2
            do k = 1 , kzm1
              do j = jci1 , jci2
                zamu = s(j,i,k+1) * zdtrdz
                if ( zamu >= d_zero ) then
                  is = d_one
                  k1 = k + 1
                  k1p1 = k1 + 1
                  if ( k1p1 > kz ) k1p1 = kz
                else
                  is = -d_one
                  k1 = k - 1
                  k1p1 = k
                  if ( k1 < 1 ) k1 = 1
                end if
                r = rdeno(wz(j,i,k1),wz(j,i,k1p1),wz(j,i,k),wz(j,i,k+1))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu * b - is * b
                wfw(j,k+1) = d_half*s(j,i,k+1) * ((d_one+zphi)*wz(j,i,k+1) + &
                                                  (d_one-zphi)*wz(j,i,k))
                !wfw(j,k+1) = d_half * zamu * ((d_one+zphi)*wz(j,i,k+1) + &
                !                              (d_one-zphi)*wz(j,i,k))
              end do
            end do
            do j = jci1 , jci2
              zrfmd = zdtrdz * fmz(j,i,1)/fmzf(j,i,2)
              zdv = -s(j,i,2) * zrfmd * wz(j,i,1)
              wz(j,i,1) = wz(j,i,1) + wfw(j,2) * zrfmd + zdv
            end do
            do k = 2 , kz
              do j = jci1 , jci2
                zrfmu = zdtrdz * fmz(j,i,k)/fmzf(j,i,k)
                zrfmd = zdtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
                zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * wz(j,i,k)
                wz(j,i,k) = wz(j,i,k) - wfw(j,k)*zrfmu + wfw(j,k+1)*zrfmd + zdv
              end do
            end do
            !do k = 1 , kz
            !  do j = jci1 , jci2
            !    zdv = (s(j,i,k) - s(j,i,k+1)) * zdtrdz * pp(j,i,k)
            !    wz(j,i,k) = pp(j,i,k) - wfw(j,k) + wfw(j,k+1) + zdv
            !  end do
            !end do
          end do
        end if

        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice1,k) = wz(j,ici1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice2,k) = wz(j,ici2,k)
            end do
          end do
        end if

        call exchange_bt(wz,2,jci1,jci2,ice1,ice2,1,kz)

        if ( lrotllr ) then

          ! Meridional advection

          do k = 1 , kz
            do i = ici1 , ice2ga
              do j = jci1 , jci2
                zamu = v(j,i,k) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i-1
                else
                  is = -d_one
                  ih = min(i+1,icross2)
                end if
                ihm1 = max(ih-1,icross1)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,i-1,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpby(j,i) = d_half * v(j,i,k) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
                !zpby(j,i) = d_half * zamu * &
                !  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
              end do
            end do
            do i = ici1 , ici2
              do j = jci1 , jci2
                zhxvtn = zdtrdy * rmv(j,i+1) * mx(j,i)
                zhxvts = zdtrdy * rmv(j,i) * mx(j,i)
                zrfmn = zhxvtn * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
                zrfms = zhxvts * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
                zdv = (v(j,i+1,k) * zrfmn - v(j,i,k) * zrfms) * pp(j,i,k)
                p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zrfms - zpby(j,i+1)*zrfmn + zdv
              end do
            end do
            !do i = ici1 , ici2
            !  do j = jci1 , jci2
            !    zhxvtn = zdtrdy * rmv(j,i+1) * mx(j,i)
            !    zhxvts = zdtrdy * rmv(j,i) * mx(j,i)
            !    zdv = (v(j,i+1,k) * zhxvtn - v(j,i,k) * zhxvts) * pp(j,i,k)
            !    p0(j,i,k) = wz(j,i,k) + zpby(j,i) - zpby(j,i+1) + zdv
            !  end do
            !end do
          end do

          if ( ma%has_bdyleft ) then
            do k = 1 , kz
              do i = ici1 , ici2
                p0(jce1,i,k) = p0(jci1,i,k)
              end do
            end do
          end if

          if ( ma%has_bdyright ) then
            do k = 1 , kz
              do i = ici1 , ici2
                p0(jce2,i,k) = p0(jci2,i,k)
              end do
            end do
          end if

          call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

          ! Zonal advection

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jce2ga
                zamu = u(j,i,k) * mu(j,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j-1
                else
                  is = -d_one
                  jh = min(j+1,jcross2)
                end if
                jhm1 = max(jh-1,jcross1)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(j-1,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbw(j,i) = d_half * u(j,i,k) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
                !zpbw(j,i) = d_half * zamu * &
                !     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
              end do
            end do
            do i = ici1 , ici2
              do j = jci1 , jci2
                zcostx = zdtrdx * mu(j,i)
                zrfme = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
                zrfmw = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
                zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
                pp(j,i,k) = p0(j,i,k) + &
                     zpbw(j,i)*zrfmw - zpbw(j+1,i)*zrfme + zdv
              end do
            end do
            !do i = ici1 , ici2
            !  do j = jci1 , jci2
            !    zrfme = mu(j+1,i) * zdtrdx
            !    zrfmw = mu(j,i) * zdtrdx
            !    zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
            !    pp(j,i,k) = p0(j,i,k) + zpbw(j,i) - zpbw(j+1,i) + zdv
            !  end do
            !end do
          end do

        else

          ! Meridional advection

          do k = 1 , kz
            do i = ici1 , ice2ga
              do j = jci1 , jci2
                zamu = v(j,i,k) * rmv(j,i) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i-1
                else
                  is = -d_one
                  ih = min(i+1,icross2)
                end if
                ihm1 = max(ih-1,icross1)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,i-1,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpby(j,i) = d_half * v(j,i,k) * rmv(j,i) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
                !zpby(j,i) = d_half * zamu * &
                !  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
              end do
            end do

            do i = ici1 , ici2
              do j = jci1 , jci2
                zrfmn = zdtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
                zrfms = zdtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
                zdv = (v(j,i+1,k) * rmv(j,i+1) * zrfmn - &
                       v(j,i,k)   * rmv(j,i)   * zrfms) * pp(j,i,k)
                p0(j,i,k) = wz(j,i,k) + &
                  mx2(j,i) * (zpby(j,i)*zrfms - zpby(j,i+1)*zrfmn + zdv)
              end do
            end do
            !do i = ici1 , ici2
            !  do j = jci1 , jci2
            !    zdv = (v(j,i+1,k) * rmv(j,i+1) - &
            !           v(j,i,k)   * rmv(j,i)   ) * zdtrdy * pp(j,i,k)
            !    p0(j,i,k) = wz(j,i,k) + &
            !      mx2(j,i) * (zpby(j,i) - zpby(j,i+1) + zdv)
            !  end do
            !end do
          end do

          if ( ma%has_bdyleft ) then
            do k = 1 , kz
              do i = ici1 , ici2
                p0(jce1,i,k) = p0(jci1,i,k)
              end do
            end do
          end if

          if ( ma%has_bdyright ) then
            do k = 1 , kz
              do i = ici1 , ici2
                p0(jce2,i,k) = p0(jci2,i,k)
              end do
            end do
          end if

          call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

          ! Zonal advection

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jce2ga
                zamu = u(j,i,k) * rmu(j,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j-1
                else
                  is = -d_one
                  jh = min(j+1,jcross2)
                end if
                jhm1 = max(jh-1,jcross1)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(j-1,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbw(j,i) = d_half * u(j,i,k) * rmu(j,i) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
                !zpbw(j,i) = d_half * zamu * &
                !     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
              end do
            end do

            do i = ici1 , ici2
              do j = jci1 , jci2
                zrfme = zdtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
                zrfmw = zdtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
                zdv = (u(j+1,i,k) * rmu(j+1,i) * zrfme - &
                       u(j,i,k)   * rmu(j,i)   * zrfmw) * pp(j,i,k)
                pp(j,i,k) = p0(j,i,k) + &
                  mx2(j,i) * (zpbw(j,i)*zrfmw - zpbw(j+1,i)*zrfme + zdv)
              end do
            end do
            !do i = ici1 , ici2
            !  do j = jci1 , jci2
            !    zdv = (u(j+1,i,k) * rmu(j+1,i) - &
            !           u(j,i,k)   * rmu(j,i)   ) * zdtrdx * pp(j,i,k)
            !    pp(j,i,k) = p0(j,i,k) + &
            !      mx2(j,i) * (zpbw(j,i) - zpbw(j+1,i) + zdv)
            !  end do
            !end do
          end do
        end if

      end subroutine wafone

      subroutine reset_tendencies
        implicit none

        s(:,:,:) = d_zero
        deltaw(:,:,:) = d_zero
        zdiv2(:,:,:) = d_zero
        mo_atm%tten = d_zero
        mo_atm%qxten = d_zero
        mo_atm%uten = d_zero
        mo_atm%vten = d_zero
        if ( ichem == 1 ) then
          mo_atm%chiten = d_zero
        end if
        if ( ibltyp == 2 ) then
          mo_atm%tketen = d_zero
        end if

        cldfra(:,:,:) = d_zero
        cldlwc(:,:,:) = d_zero

        if ( idiag > 0 ) then
          ten0 = t(jci1:jci2,ici1:ici2,:)
          qen0 = qv(jci1:jci2,ici1:ici2,:)
          if ( ichem == 1 ) then
            chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
          end if
        end if
      end subroutine reset_tendencies

      subroutine physical_parametrizations
        implicit none
        integer(ik4) :: i , j , k , n
        logical :: loutrad , labsem

#ifdef DEBUG
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              if ( (t(j,i,k) > 350.0_rkx) .or. t(j,i,k) < 170.0_rkx ) then
                write(100+myid,*) 'Before Phys On : ', myid
                write(100+myid,*) 'At : ', i,j,k
                write(100+myid,*) 'k pai u v w qv qc t tetav'
                do n = 1 , kz
                  write(100+myid,*) n, pai(j,i,n), u(j,i,n), v(j,i,n), &
                            w(j,i,n) , qv(j,i,n) , qc(j,i,n) , &
                            t(j,i,n) , tetav(j,i,n)
                end do
                flush(100+myid)
                call fatal(__FILE__,__LINE__, 'error')
              end if
            end do
          end do
        end do
#endif
        if ( any(icup > 0) ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%tten(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv)
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            chiten0 = mo_atm%chiten(jci1:jci2,ici1:ici2,:,:)
          end if
          call cumulus
          if ( ichem == 1 ) then
            if ( ichcumtra == 1 ) then
              if ( debug_level > 3 .and. myid == italk ) then
                write(stdout,*) 'Calling cumulus transport at ', &
                           trim(rcmtimer%str())
              end if
              call cumtran(trac)
            end if
          end if
          if ( idiag > 0 ) then
            tdiag%con = mo_atm%tten(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%con = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            cconvdiag = mo_atm%chiten(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        else
          if ( any(icup < 0) ) then
            call shallow_convection
            if ( idiag > 0 ) then
              tdiag%con = mo_atm%tten(jci1:jci2,ici1:ici2,:) - ten0
              qdiag%con = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv) - qen0
            end if
          end if
        end if
        !
        !------------------------------------------------
        ! Large scale precipitation microphysical schemes
        !------------------------------------------------
        !
        if ( ipptls > 0 ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%tten(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv)
          end if
          ! Cumulus clouds
          if ( icldfrac /= 2 ) then
            call cucloud
          end if
          ! Save cumulus cloud fraction for chemistry before it is
          ! overwritten in cldfrac
          if ( ichem == 1 ) then
            convcldfra(:,:,:) = cldfra(:,:,:)
          end if
          ! Clouds and large scale precipitation
          call cldfrac(cldlwc,cldfra)
          call microscheme
          if ( idiag > 0 ) then
            tdiag%lsc = mo_atm%tten(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%lsc = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        end if
        !
        !------------------------------------------------
        !       Call radiative transfer package
        !------------------------------------------------
        !
        if ( rcmtimer%start() .or. syncro_rad%will_act( ) ) then
          if ( debug_level > 3 .and. myid == italk ) then
            write(stdout,*) &
              'Calling radiative transfer at ',trim(rcmtimer%str())
          end if
          ! calculate albedo
          call surface_albedo
          ! Update / init Ozone profiles
          if ( iclimao3 == 1 ) then
            call updateo3(rcmtimer%idate,scenario)
          else
            if ( rcmtimer%start() ) call inito3
          end if
          if ( iclimaaer == 1 ) then
            call updateaerosol(rcmtimer%idate)
          else if ( iclimaaer == 2 ) then
            call updateaeropp(rcmtimer%idate)
          else if ( iclimaaer == 3 ) then
            call updateaeropp_cmip6(rcmtimer%idate)
          end if
          loutrad = ( rcmtimer%start() .or. alarm_out_rad%will_act(dtrad) )
          labsem = ( rcmtimer%start() .or. syncro_emi%will_act() )
          if ( debug_level > 3 .and. myid == italk ) then
            if ( labsem ) then
              write(stdout,*) 'Updating abs-emi at ',trim(rcmtimer%str())
            end if
            if ( loutrad ) then
              write(stdout,*) 'Collecting radiation at ',trim(rcmtimer%str())
            end if
          end if
          call radiation(rcmtimer%year,rcmtimer%month,loutrad,labsem)
        end if
        !
        ! Add radiative transfer package-calculated heating rates to
        ! temperature tendency (deg/sec)
        !
        do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
          mo_atm%tten(j,i,k) = mo_atm%tten(j,i,k) + heatrt(j,i,k)
        end do
        if ( idiag > 0 ) tdiag%rad = heatrt
        !
        !------------------------------------------------
        !            Call Surface model
        !------------------------------------------------
        !
        if ( rcmtimer%start() .or. syncro_srf%will_act( ) ) then
          if ( debug_level > 3 .and. myid == italk ) then
            write(stdout,*) 'Calling surface model at ',trim(rcmtimer%str())
          end if
          call surface_model
          !FAB now called in surface model
          ! if ( islab_ocean == 1 ) call update_slabocean(xslabtime)
        end if
        !
        !------------------------------------------------
        !             Call PBL scheme
        !------------------------------------------------
        !
        if ( ibltyp > 0 ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%tten(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv)
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            chiten0 = mo_atm%chiten(jci1:jci2,ici1:ici2,:,:)
          end if
          call pblscheme
          if ( idiag > 0 ) then
            tdiag%tbl = mo_atm%tten(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%tbl = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            ctbldiag = mo_atm%chiten(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        end if
        if ( ipptls == 1 ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%tten(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv)
          end if
          call condtq
          if ( idiag > 0 ) then
            tdiag%lsc = tdiag%lsc + mo_atm%tten(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%lsc = qdiag%lsc + &
                 mo_atm%qxten(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        end if
        !------------------------------------------------------------- 
        !call chemistry/aerosol schemes 
        !----------------------------------------------------------------------
        if ( ichem == 1 ) then
          call tractend2(rcmtimer%month,rcmtimer%day,declin)
        end if
        !
        ! Update status
        !
        do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
          t(j,i,k) = t(j,i,k) + dtsec * mo_atm%tten(j,i,k)
        end do
        do concurrent ( j = jdi1:jdi2 , i = ici1:ici2 , k = 1:kz )
          u(j,i,k) = u(j,i,k) + dtsec * mo_atm%uten(j,i,k)
        end do
        do concurrent ( j = jci1:jci2 , i = idi1:idi2 , k = 1:kz )
          v(j,i,k) = v(j,i,k) + dtsec * mo_atm%vten(j,i,k)
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              qx(j,i,k,iqv) = qx(j,i,k,iqv) + mo_atm%qxten(j,i,k,iqv)*dtsec
              qx(j,i,k,iqv) = max(qx(j,i,k,iqv),minqq)
            end do
          end do
        end do
        do n = iqfrst , iqlst
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                qx(j,i,k,n) = qx(j,i,k,n) + mo_atm%qxten(j,i,k,n)*dtsec
                qx(j,i,k,n) = max(qx(j,i,k,n),d_zero)
              end do
            end do
          end do
        end do
        if ( ibltyp == 2 ) then
          do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kzp1 )
            tke(j,i,k) = max(tke(j,i,k) + dtsec * mo_atm%tketen(j,i,k),tkemin)
          end do
        end if
        if ( ichem == 1 ) then
          do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz , n = 1:ntr )
            trac(j,i,k,n) = trac(j,i,k,n) + dtsec * mo_atm%chiten(j,i,k,n)
            trac(j,i,k,n) = max(trac(j,i,k,n),d_zero)
          end do
        end if
#ifdef DEBUG
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              if ( (t(j,i,k) > 350.0_rkx) .or. t(j,i,k) < 170.0_rkx ) then
                write(100+myid,*) 'On : ', myid
                write(100+myid,*) 'After Phys At : ', i,j,k
                write(100+myid,*) 'k pai u v w qv qc t tetav'
                do n = 1 , kz
                  write(100+myid,*) n, pai(j,i,n), u(j,i,n), v(j,i,n), &
                            w(j,i,n) , qv(j,i,n) , qc(j,i,n) , &
                            t(j,i,n) , tetav(j,i,n)
                end do
                flush(100+myid)
                call fatal(__FILE__,__LINE__, 'error')
              end if
            end do
          end do
        end do
#endif
      end subroutine physical_parametrizations

  end subroutine moloch

  subroutine wstagtox(w,wx)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: w
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    integer(ik4) :: i , j , k , i1 , i2 , j1 , j2
    i1 = lbound(wx,2)
    i2 = ubound(wx,2)
    j1 = lbound(wx,1)
    j2 = ubound(wx,1)
    do k = 2 , kzm1
      do i = i1 , i2
        do j = j1 , j2
          wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                      0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
        end do
      end do
    end do
    do i = i1 , i2
      do j = j1 , j2
        wx(j,i,1)  = d_half * (w(j,i,2)+w(j,i,1))
        wx(j,i,kz) = d_half * (w(j,i,kzp1)+w(j,i,kz))
      end do
    end do
  end subroutine wstagtox

  subroutine xtowstag(wx,w)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: w
    integer(ik4) :: i , j , k
    do k = 3 , kzm1
      do i = ice1 , ice2
        do j = jce1 , jce2
          w(j,i,k) = 0.5625_rkx * (wx(j,i,k)  +wx(j,i,k-1)) - &
                     0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
        end do
      end do
    end do
    do i = ice1 , ice2
      do j = jce1 , jce2
        w(j,i,2) = d_half * (wx(j,i,2)  +wx(j,i,1))
        w(j,i,kz) = d_half * (wx(j,i,kz)+wx(j,i,kzm1))
      end do
    end do
  end subroutine xtowstag

  subroutine xtoustag(ux,u)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u
    integer(ik4) :: i , j , k

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jdii1 , jdii2
          u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                     0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
        end do
      end do
    end do
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
        end do
      end do
    end if
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
        end do
      end do
    end if
  end subroutine xtoustag

  subroutine xtovstag(vx,v)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: v
    integer(ik4) :: i , j , k

    do k = 1 , kz
      do i = idii1 , idii2
        do j = jci1 , jci2
          v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                     0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
        end do
      end do
    end do
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
        end do
      end do
    end if
  end subroutine xtovstag

  subroutine xtouvstag(ux,vx,u,v)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    integer(ik4) :: i , j , k

    call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

    ! Back to wind points: U (fourth order)

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jdii1 , jdii2
          u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                     0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
        end do
      end do
    end do
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
        end do
      end do
    end if
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
        end do
      end do
    end if

    ! Back to wind points: V (fourth order)
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jci1 , jci2
          v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                     0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
        end do
      end do
    end do
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
        end do
      end do
    end if
  end subroutine xtouvstag

  subroutine uvstagtox(u,v,ux,vx)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    integer(ik4) :: i , j , k

    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)

    ! Compute U-wind on T points
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jci1 , jci2
          ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                      0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
        end do
      end do
    end do
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce1,i,k) = d_half * (u(jde1,i,k)+u(jdi1,i,k))
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce2,i,k) = d_half*(u(jde2,i,k) + u(jdi2,i,k))
        end do
      end do
    end if
    ! Compute V-wind on T points
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jce1 , jce2
          vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                      0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
        end do
      end do
    end do
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice1,k) = d_half * (v(j,ide1,k)+v(j,idi1,k))
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice2,k) = d_half*(v(j,ide2,k) + v(j,idi2,k))
        end do
      end do
    end if
  end subroutine uvstagtox

  subroutine divdamp(dts)
    implicit none
    real(rkx) , intent(in) :: dts
    integer(ik4) :: i , j , k
    real(rkx) :: ddamp1

    ddamp1 = ddamp*0.125_rkx*(dx**2)/dts
    if ( lrotllr ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                ddamp1/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                ddamp1/dx*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                ddamp1/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                ddamp1/(dx*rmv(j,i))*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
    end if
  end subroutine divdamp

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
