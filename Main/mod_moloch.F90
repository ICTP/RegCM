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
  real(rkx) , dimension(:) , pointer :: xknu , zprof
  real(rkx) , dimension(:,:) , pointer :: p2d
  real(rkx) , dimension(:,:) , pointer :: xlat , xlon , coru , corv
  real(rkx) , dimension(:,:) , pointer :: mu , hx , mx
  real(rkx) , dimension(:,:) , pointer :: mv , hy
  real(rkx) , dimension(:,:) , pointer :: ps , ts , ht
  real(rkx) , dimension(:,:,:) , pointer :: fmz
  real(rkx) , dimension(:,:,:) , pointer :: fmzf
  real(rkx) , dimension(:,:,:) , pointer :: pai , pf
  real(rkx) , dimension(:,:,:) , pointer :: tetav , tf , tvirt
  real(rkx) , dimension(:,:,:) , pointer :: zeta , zetau , zetav
  real(rkx) , dimension(:,:,:) , pointer :: u , v , w
  real(rkx) , dimension(:,:,:) , pointer :: ux , vx
  real(rkx) , dimension(:,:,:) , pointer :: ud , vd
  real(rkx) , dimension(:,:,:) , pointer :: p , t , rho
  real(rkx) , dimension(:,:,:) , pointer :: qv , qf , qc , qi , qr , qs , qsat
  real(rkx) , dimension(:,:,:) , pointer :: qwltot , qwitot
  real(rkx) , dimension(:,:,:) , pointer :: tke
  real(rkx) , dimension(:,:,:,:) , pointer :: qx , trac

  public :: allocate_moloch , init_moloch , moloch
  public :: uvstagtox , xtouvstag , wstagtox

  real(rkx) , parameter :: minden = 1.0e-30_rkx

  logical , parameter :: do_phys         = .true.
  logical , parameter :: do_fulleq       = .true.
  logical , parameter :: do_vadvtwice    = .true.
  logical , parameter :: do_bdy          = .true.
  logical , parameter :: do_divdamp      = .true.
  logical , parameter :: do_filterpai    = .false.
  logical , parameter :: do_filterqv     = .false.
  logical , parameter :: do_filterdiv    = .false.
  logical , parameter :: do_filtertheta  = .false.
#ifdef RCEMIP
  logical , parameter :: do_diffutend    = .false.
#endif

  logical :: moloch_realcase = (.not. moloch_do_test_1) .and. &
                               (.not. moloch_do_test_2)
  logical :: lrotllr

  real(rkx) , parameter :: nupaitq = 0.05_rkx
  real(rkx) , parameter :: ddamp = 0.25_rkx

  real(rkx) :: dzita
  integer(ik4) :: jmin , jmax , imin , imax

  contains

#include <pfesat_acc.inc>
#include <pfwsat_acc.inc>

  subroutine allocate_moloch
    implicit none
    integer(ik4) :: k
    call getmem1d(gzitak,1,kzp1,'moloch:gzitak')
!$acc enter data create(gzitak)
    call getmem1d(gzitakh,1,kz,'moloch:gzitakh')
!$acc enter data create(gzitakh)
    call getmem2d(p2d,jdi1,jdi2,idi1,idi2,'moloch:p2d')
!$acc enter data create(p2d)
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:deltaw')
!$acc enter data create(deltaw)
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
!$acc enter data create(s)
    call getmem3d(wx,jce1,jce2,ice1,ice2,1,kz,'moloch:wx')
!$acc enter data create(wx)
    call getmem3d(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
!$acc enter data create(zdiv2)
    call getmem3d(wwkw,jci1,jci2,ici1,ici2,1,kzp1,'moloch:wwkw')
!$acc enter data create(wwkw)
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
!$acc enter data create(wz)
    call getmem2d(wfw,jci1,jci2,1,kzp1,'moloch:wfw')
!$acc enter data create(wfw)
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
!$acc enter data create(p0)
    call getmem2d(zpby,jci1,jci2,ici1,ice2ga,'moloch:zpby')
!$acc enter data create(zpby)
    call getmem2d(zpbw,jci1,jce2ga,ici1,ici2,'moloch:zpbw')
!$acc enter data create(zpbw)
    call getmem2d(mx2,jde1,jde2,ide1,ide2,'moloch:mx2')
!$acc enter data create(mx2)
    call getmem2d(rmu,jde1ga,jde2ga,ide1,ide2,'moloch:rmu')
!$acc enter data create(rmu)
    call getmem2d(rmv,jde1,jde2,ide1ga,ide2ga,'moloch:rmv')
!$acc enter data create(rmv)
    call getmem2d(coru,jde1,jde2,ice1,ice2,'moloch:coru')
!$acc enter data create(coru)
    call getmem2d(corv,jce1,jce2,ide1,ide2,'moloch:corv')
!$acc enter data create(corv)
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jce1,jce2,ice1,ice2,1,kz,'moloch:tkex')
!$acc enter data create(tkex)
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
!$acc enter data create(ud, vd)
    if ( ifrayd == 1 ) then
      call getmem3d(zetau,jdi1,jdi2,ici1,ici2,1,kz,'moloch:zetau')
      call getmem3d(zetav,jci1,jci2,idi1,idi2,1,kz,'moloch:zetav')
!$acc enter data create(zetau, zetav)
    end if
    if ( do_fulleq ) then
      call getmem3d(qwltot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwltot')
      call getmem3d(qwitot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwitot')
!$acc enter data create(qwltot,qwitot)
    end if
    call getmem1d(xknu,1,kz,'moloch:xknu')
!$acc enter data create(xknu)
    call getmem1d(zprof,1,kz,'moloch:zprof')
!$acc enter data create(zprof)
    do k = 1 , kz
      xknu(k) = sin(d_half*mathpi*(1.0_rkx-real((k-1),rkx)/kz))*mo_anu2
    end do
    if ( do_filterpai ) then
      call getmem3d(pf,jce1,jce2,ice1,ice2,1,kz,'moloch:pf')
!$acc enter data create(pf)
    end if
    if ( do_filtertheta ) then
      call getmem3d(tf,jce1,jce2,ice1,ice2,1,kz,'moloch:tf')
!$acc enter data create(tf)
    end if
    if ( do_filterqv ) then
      call getmem3d(qf,jce1,jce2,ice1,ice2,1,kz,'moloch:qf')
!$acc enter data create(qf)
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
    call assignpnt(sfs%tg,ts)
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
    if ( ibltyp == 2 ) then
      call assignpnt(mo_atm%tke,tke)
    end if
    if ( ichem == 1 ) then
      call assignpnt(mo_atm%trac,trac)
    end if
    if ( ifrayd == 1 ) then
      call xtoustag(zeta,zetau)
      call xtovstag(zeta,zetav)
!$acc update self(zetau, zetav)
    end if
#ifdef RCEMIP
    coru = 0.0_rkx
    corv = 0.0_rkx
#else
    coru = eomeg2*sin(mddom%ulat(jde1:jde2,ice1:ice2)*degrad)
    corv = eomeg2*sin(mddom%vlat(jce1:jce2,ide1:ide2)*degrad)
#endif
    mx2 = mx * mx
    rmu = d_one/mu
    rmv = d_one/mv
    gzitak = gzita(zita)
    gzitakh = gzita(zitah)
    dzita = mo_dzita
    wwkw(:,:,kzp1) = d_zero
    w(:,:,1) = d_zero
    lrotllr = (iproj == 'ROTLLR')

    jmin = jcross1
    jmax = jcross2
    imin = icross1
    imax = icross2
    if ( ma%bandflag ) then
      jmin = jcross1 - 2
      jmax = jcross2 + 2
    end if
    if ( ma%crmflag ) then
      jmin = jcross1 - 2
      jmax = jcross2 + 2
      imin = icross1 - 2
      imax = icross2 + 2
    end if

! Update static arrays on device
!$acc update device(mu, mv, rmu, rmv, mx, mx2, fmz, fmzf, hx, hy,&
!$acc& gzitak, gzitakh, wwkw, w, coru, corv, mo_atm%zeta, mo_atm%zetaf)

! Update dynamic arrays to device
!$acc update device(u, v, pai, t, qx, ps, ts)
  end subroutine init_moloch
  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch
    implicit none
    real(rkx) :: dtsound , dtstepa
    real(rkx) :: maxps , minps , pmax , pmin , zdgz
    real(rkx) :: tv , lrt , fice
    !real(rk8) :: jday
    integer(ik4) :: i , j , k , nadv
    integer(ik4) :: iconvec
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    dtstepa = dtsec / real(mo_nadv,rkx)
    dtsound = dtstepa / real(mo_nsound,rkx)
    iconvec = 0
    !
    ! Start of accelerated section
    !
    on_device = .true.

    call reset_tendencies
!$acc parallel present(p, pai, qsat, t)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          p(j,i,k) = (pai(j,i,k)**cpovr) * p00
          qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
        end do
      end do
    end do
!$acc end parallel

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
!$acc parallel present(tvirt, t, qv, qc, qi, qr, qs)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
                                         qc(j,i,k) - qi(j,i,k) - &
                                         qr(j,i,k) - qs(j,i,k))
            end do
          end do
        end do
!$acc end parallel
        if ( do_fulleq ) then
!$acc parallel present(qwltot, qwitot, qc, qi, qr, qs)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                qwltot(j,i,k) = qc(j,i,k) + qr(j,i,k)
                qwitot(j,i,k) = qi(j,i,k)
              end do
            end do
          end do
!$acc end parallel
        end if
      else
!$acc parallel present(tvirt, t, qv, qc)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
!$acc end parallel
        if ( do_fulleq ) then
!$acc parallel present(t, qwltot, qwitot, qc) private(fice)
!$acc loop collapse(3)
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
!$acc end parallel
        end if
      end if
    else
!$acc parallel present(tvirt, t, qv)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
!$acc end parallel
      if ( do_fulleq ) then
!$acc kernels present(qwltot, qwitot)
        qwltot(:,:,:) = d_zero
        qwitot(:,:,:) = d_zero
!$acc end kernels
      end if
    end if

!$acc parallel present(tetav, tvirt, pai)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
        end do
      end do
    end do
!$acc end parallel

    if ( idiag > 0 ) then
!$acc kernels present(t, qv, ten0, qen0)
      ten0 = t(jci1:jci2,ici1:ici2,:)
      qen0 = qv(jci1:jci2,ici1:ici2,:)
!$acc end kernels
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
!$acc kernels present(trac, chiten0)
        chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
!$acc end kernels
      end if
    end if

    if ( do_filterpai ) then
!$acc kernels present(pf, pai)
      pf = pai(jce1:jce2,ice1:ice2,:)
!$acc end kernels
    end if
    if ( do_fulleq ) then
      if ( do_filtertheta ) then
!$acc kernels present(tf, tetav)
        tf = tetav(jce1:jce2,ice1:ice2,:)
!$acc end kernels
      end if
    end if
    if ( do_filterqv ) then
!$acc kernels present(qf, qv)
      qf = qv(jce1:jce2,ice1:ice2,:)
!$acc end kernels
    end if

    do nadv = 1 , mo_nadv

      call sound(dtsound)

      call advection(dtstepa)

    end do ! Advection loop
!$acc update self(u,v,qx,tke) async(2)

    if ( do_filterpai ) then
!$acc kernels present(pai, pf)
      pai(jce1:jce2,ice1:ice2,:) = pai(jce1:jce2,ice1:ice2,:) - pf
!$acc end kernels
      call filtpai
!$acc kernels present(pai, pf)
      pai(jce1:jce2,ice1:ice2,:) = pai(jce1:jce2,ice1:ice2,:) + pf
!$acc end kernels
    end if

!$acc update self(pai) async(2)

    if ( do_fulleq ) then
      if ( do_filtertheta ) then
!$acc kernels present(tetav, tf)
        tetav(jce1:jce2,ice1:ice2,:) = tetav(jce1:jce2,ice1:ice2,:) - tf
!$acc end kernels
        call filttheta
!$acc kernels present(tetav, tf)
        tetav(jce1:jce2,ice1:ice2,:) = tetav(jce1:jce2,ice1:ice2,:) + tf
!$acc end kernels
      end if
!$acc update self(tetav) async(2)
    end if

!$acc parallel present(tvirt, tetav, pai)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
        end do
      end do
    end do
!$acc end parallel

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
!$acc parallel present(t, tvirt, qv, qc, qi, qr, qs)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
                             qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
            end do
          end do
        end do
!$acc end parallel
      else
!$acc parallel present(t, tvirt, qv, qc)
!$acc loop collapse(3)
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
!$acc end parallel
      end if
    else
!$acc parallel present(t, tvirt, qv)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
!$acc end parallel
    end if

!$acc update self(t) async(2)

    if ( idiag > 0 ) then
!$acc wait(2)
      tdiag%adh = (t(jci1:jci2,ici1:ici2,:) - ten0) * rdt
      qdiag%adh = (qv(jci1:jci2,ici1:ici2,:) - qen0) * rdt
    end if

    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        cadvhdiag = (trac(jci1:jci2,ici1:ici2,:,:) - chiten0) * rdt
      end if
    end if

!$acc parallel present(p, pai, rho, t)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          p(j,i,k) = (pai(j,i,k)**cpovr) * p00
          rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
        end do
      end do
    end do
!$acc end parallel

    !jday = yeardayfrac(rcmtimer%idate)
!$acc parallel present(zeta, tvirt, ps, ts, p) private(zdgz, lrt, tv, invt)
!$acc loop collapse(2)
    do i = ice1 , ice2
      do j = jce1 , jce2
        zdgz = zeta(j,i,kz)*egrav
        lrt = (tvirt(j,i,kz)-tvirt(j,i,kz-1))/(zeta(j,i,kz-1)-zeta(j,i,kz))
        ! lrt = 0.65_rkx*lrt + 0.35_rkx*stdlrate(jday,xlat(j,i))
        lrt = 0.75_rkx*lrt + 0.25_rkx*lrate
        tv = tvirt(j,i,kz) + 0.5_rkx*zeta(j,i,kz)*lrt ! Mean temperature
        ps(j,i) = p(j,i,kz) * exp(zdgz/(rgas*tv))
      end do
    end do
!$acc end parallel
!$acc update self(ps) async(2)
    !
    ! Recompute saturation
    !
!$acc parallel present(qsat, t, p)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
        end do
      end do
    end do
!$acc end parallel
    !
    ! Lateral/damping boundary condition
    !
    if ( do_bdy .and. moloch_realcase .and. irceideal == 0 ) then
      call boundary
      if ( i_crm /= 1 ) then
        if ( ifrayd == 1 ) then
!$acc wait(2)
          call raydamp(zetau,u,xub,jdi1,jdi2,ici1,ici2,1,kz)
          call raydamp(zetav,v,xvb,jci1,jci2,idi1,idi2,1,kz)
          call raydamp(zeta,t,xtb,jci1,jci2,ici1,ici2,1,kz)
          call raydamp(zeta,pai,xpaib,jci1,jci2,ici1,ici2,1,kz)
!$acc update device(u, v, t, pai) async(2)
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
    call uvstagtox(u,v,ux,vx)
    call mkslice
    !
    ! PHYSICS
    !
    if ( do_phys .and. moloch_realcase ) then
!$acc update self(w, qsat, p, rho, ps, tvirt, tetav) async(2)
!$acc update self(tke) async(2) if(ibltyp == 2)
!$acc update self(trac) async(2) if(ichem == 1)
!$acc wait(2)
      on_device = .false.
      call physical_parametrizations
      on_device = .true.
!$acc update device(u, v, w, t, qx, qsat, p, rho, ps, tvirt, tetav) async(2)
!$acc update device(tke) async(2) if(ibltyp == 2)
!$acc update device(trac) async(2) if(ichem == 1)
!$acc wait(2)
    else
      if ( debug_level > 1 ) then
        if ( myid == italk ) then
          write(stdout,*) 'WARNING: Physical package disabled!!!'
        end if
      end if
    end if

    if ( do_filterqv ) then
!$acc kernels present(qv, qf)
      qv(jce1:jce2,ice1:ice2,:) = qv(jce1:jce2,ice1:ice2,:) - qf
!$acc end kernels
      call filtqv
!$acc kernels present(qv, qf)
      qv(jce1:jce2,ice1:ice2,:) = max(qv(jce1:jce2,ice1:ice2,:) + qf,minqq)
!$acc end kernels
    end if

!$acc update self(qv) async(2)

    !
    ! Mass check
    !
    if ( debug_level > 0 ) call massck
    !
    ! Diagnostic and end timestep
    !
    if ( syncro_rep%act( ) .and. rcmtimer%integrating( ) ) then
!$acc wait(2)
      maxps = maxval(ps(jci1:jci2,ici1:ici2))
      minps = minval(ps(jci1:jci2,ici1:ici2))
      call maxall(maxps,pmax)
      call minall(minps,pmin)
      call sumall(total_precip_points,iconvec)
      if ( is_nan(pmax) .or. is_nan(pmin) .or. &
           is_inf(pmax) .or. is_inf(pmin) ) then
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
    on_device = .false.

    contains

      subroutine boundary
        implicit none
#ifdef USE_MPI3
        type(commdata_real) :: comm1 , comm2 , comm3
#else
        logical :: do_nudge
        do_nudge = ( iboudy == 1 .or. iboudy >= 5 .or. iboudy == 4)
        call exchange_lrbt(ps,1,jce1,jce2,ice1,ice2)
!$acc update self(ps) async(2) if(do_nudge)
        call exchange_lrbt(u,1,jde1,jde2,ice1,ice2,1,kz)
!$acc update self(u) async(2) if(do_nudge)
        call exchange_lrbt(v,1,jce1,jce2,ide1,ide2,1,kz)
!$acc update self(v) async(2) if(do_nudge)
        call exchange_lrbt(t,1,jce1,jce2,ice1,ice2,1,kz)
!$acc update self(t) async(2) if(do_nudge)
        call exchange_lrbt(qv,1,jce1,jce2,ice1,ice2,1,kz)
!$acc update self(qv) async(2) if(do_nudge)
        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
!$acc update self(pai) async(2) if(do_nudge)
        if ( is_present_qc( ) ) then
!$acc update self(qc) async(2) if(do_nudge)
          call exchange_lrbt(qc,1,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( is_present_qi( ) ) then
!$acc update self(qi) async(2) if(do_nudge)
          call exchange_lrbt(qi,1,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( (iboudy == 1 .or. iboudy >= 5) .and. ichem == 1 ) then
          on_device = .false.
          call exchange_lrbt(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          on_device = .true.
        end if
#endif
        if ( idiag > 0 ) then
!$acc wait(2)
          ten0 = t(jci1:jci2,ici1:ici2,:)
          qen0 = qv(jci1:jci2,ici1:ici2,:)
        end if
        if ( ichem == 1 ) then
          if ( ichdiag > 0 ) then
            chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
          end if
        end if

        if ( iboudy == 1 .or. iboudy >= 5 ) then
#ifdef USE_MPI3
!$acc wait(2)
          call exchange_lrbt_pre(ps,1,jce1,jce2,ice1,ice2,comm1)
          call exchange_lrbt_pre(u,1,jde1,jde2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_pre(v,1,jce1,jce2,ide1,ide2,1,kz,comm3)
          call exchange_lrbt_post(ps,1,jce1,jce2,ice1,ice2,comm1)
          call nudge(iboudy,ps,xpsb)
!$acc update device(ps) async(2)
          call exchange_lrbt_pre(t,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call exchange_lrbt_post(u,1,jde1,jde2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_post(v,1,jce1,jce2,ide1,ide2,1,kz,comm3)
          call nudge(iboudy,u,v,xub,xvb)
!$acc update device(u,v) async(2)
          call exchange_lrbt_pre(qv,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_post(t,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call nudge(iboudy,t,xtb)
!$acc update device(t) async(2)
          call exchange_lrbt_pre(pai,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          call exchange_lrbt_post(qv,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          call nudge(iboudy,qv,xqb)
!$acc update device(qv) async(2)
          call exchange_lrbt_post(pai,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          call nudge(iboudy,pai,xpaib)
!$acc update device(pai) async(2)

          if ( ichem == 1 ) then
            call exchange_lrbt_pre(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr,comm1)
          end if
          if ( is_present_qc( ) ) then
            call exchange_lrbt_pre(qc,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          end if
          if ( is_present_qi( ) ) then
            call exchange_lrbt_pre(qi,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          end if
          if ( ichem == 1 ) then
            call exchange_lrbt_post(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr,comm1)
            call nudge_chi(trac)
            if ( ichdiag > 0 ) then
              cbdydiag = trac(jci1:jci2,ici1:ici2,:,:) - chiten0
            end if
          end if
          if ( is_present_qc( ) ) then
            call exchange_lrbt_post(qc,1,jce1,jce2,ice1,ice2,1,kz,comm2)
            call nudge(iboudy,qc,xlb)
!$acc update device(qc) async(2)
          end if
          if ( is_present_qi( ) ) then
            call exchange_lrbt_post(qi,1,jce1,jce2,ice1,ice2,1,kz,comm3)
            call nudge(iboudy,qi,xib)
!$acc update device(qi) async(2)
          end if
#else
!$acc wait(2)
          call nudge(iboudy,ps,xpsb)
!$acc update device(ps) async(2)
          call nudge(iboudy,u,v,xub,xvb)
!$acc update device(u,v) async(2)
          call nudge(iboudy,t,xtb)
!$acc update device(t) async(2)
          call nudge(iboudy,qv,xqb)
!$acc update device(qv) async(2)
          call nudge(iboudy,pai,xpaib)
!$acc update device(pai) async(2)
          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
          if ( is_present_qc( ) ) then
            call nudge(iboudy,qc,xlb)
!$acc update device(qc) async(2)
          end if
          if ( is_present_qi( ) ) then
            call nudge(iboudy,qi,xib)
!$acc update device(qi) async(2)
          end if
          if ( ichem == 1 ) then
            call nudge_chi(trac)
            if ( ichdiag > 0 ) then
              cbdydiag = trac(jci1:jci2,ici1:ici2,:,:) - chiten0
            end if
          end if
#endif
        else if ( iboudy == 4 ) then
#ifdef USE_MPI3
!$acc wait(2)
          call exchange_lrbt_pre(ps,1,jce1,jce2,ice1,ice2,comm1)
          call exchange_lrbt_pre(u,1,jde1,jde2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_pre(v,1,jce1,jce2,ide1,ide2,1,kz,comm3)
          call exchange_lrbt_post(ps,1,jce1,jce2,ice1,ice2,comm1)
          call sponge(ps,xpsb)
!$acc update device(ps) async(2)
          call exchange_lrbt_pre(t,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call exchange_lrbt_post(u,1,jde1,jde2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_post(v,1,jce1,jce2,ide1,ide2,1,kz,comm3)
          call sponge(u,v,xub,xvb)
!$acc update device(u,v) async(2)
          call exchange_lrbt_pre(qv,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          call exchange_lrbt_post(t,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call sponge(t,xtb)
!$acc update device(t) async(2)
          call exchange_lrbt_pre(pai,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          call exchange_lrbt_post(qv,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          call sponge(qv,xqb)
!$acc update device(qv) async(2)
          call exchange_lrbt_post(pai,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          call sponge(pai,xpaib)
!$acc update device(pai) async(2)

          if ( is_present_qc( ) ) then
            call exchange_lrbt_pre(qc,1,jce1,jce2,ice1,ice2,1,kz,comm2)
          end if
          if ( is_present_qi( ) ) then
            call exchange_lrbt_pre(qi,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          end if
          if ( is_present_qc( ) ) then
            call exchange_lrbt_post(qc,1,jce1,jce2,ice1,ice2,1,kz,comm2)
            call sponge(qc,xlb)
!$acc update device(qc) async(2)
          end if
          if ( is_present_qi( ) ) then
            call exchange_lrbt_post(qi,1,jce1,jce2,ice1,ice2,1,kz,comm3)
            call sponge(qi,xib)
!$acc update device(qi) async(2)
          end if
#else
!$acc wait(2)
          call sponge(ps,xpsb)
!$acc update device(ps) async(2)
          call sponge(u,v,xub,xvb)
!$acc update device(u,v) async(2)
          call sponge(t,xtb)
!$acc update device(t) async(2)
          call sponge(qv,xqb)
!$acc update device(qv) async(2)
          call sponge(pai,xpaib)
!$acc update device(pai) async(2)
          if ( is_present_qc( ) ) then
!$acc update self(qc)
            call sponge(qc,xlb)
!$acc update device(qc) async(2)
          end if
          if ( is_present_qi( ) ) then
!$acc update self(qi)
            call sponge(qi,xib)
!$acc update device(qi) async(2)
          end if
#endif
          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
        end if
!$acc wait(2)
      end subroutine boundary

#ifdef RCEMIP
      subroutine filt3d(p,nu)
        implicit none
        real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: p
        real(rkx) , intent(in) :: nu
        integer(ik4) :: j , i , k

!$acc parallel present(p2d, p) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = icii1 , icii2
            do j = jcii1 , jcii2
              p2d(j,i) = 0.125_rkx * (p(j-1,i,k) + p(j+1,i,k) + &
                                      p(j,i-1,k) + p(j,i+1,k)) - &
                         d_half   * p(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = icii1 , icii2
            do j = jcii1 , jcii2
              p(j,i,k) = p(j,i,k) + nu * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
      end subroutine filt3d

      subroutine filtuv(u,v,nu)
        implicit none
        real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
        real(rkx) , intent(in) :: nu
        integer(ik4) :: j , i , k

!$acc parallel present(p2d, u) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = icii1 , icii2
            do j = jdii1 , jdii2
              p2d(j,i) = 0.125_rkx * (u(j-1,i,k) + u(j+1,i,k) + &
                                      u(j,i-1,k) + u(j,i+1,k)) - &
                         d_half   * u(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = icii1 , icii2
            do j = jdii1 , jdii2
              u(j,i,k) = u(j,i,k) + nu * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel

!$acc parallel present(p2d, v) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = idii1 , idii2
            do j = jcii1 , jcii2
              p2d(j,i) = 0.125_rkx * (v(j-1,i,k) + v(j+1,i,k) + &
                                      v(j,i-1,k) + v(j,i+1,k)) - &
                         d_half   * v(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = idii1 , idii2
            do j = jcii1 , jcii2
              v(j,i,k) = v(j,i,k) + nu * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
      end subroutine filtuv

      subroutine filt4d(p,nu,n1,n2)
        implicit none
        real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: p
        real(rkx) , intent(in) :: nu
        integer(ik4) , intent(in) :: n1 , n2
        integer(ik4) :: j , i , k , n

!$acc parallel present(p2d, p) if(on_device)
!$acc loop collapse(2)
        do n = n1 , n2
          do k = 1 , kz
!$acc loop collapse(2)
            do i = icii1 , icii2
              do j = jcii1 , jcii2
                p2d(j,i) = 0.125_rkx * (p(j-1,i,k,n) + p(j+1,i,k,n) + &
                                        p(j,i-1,k,n) + p(j,i+1,k,n)) - &
                           d_half   * p(j,i,k,n)
              end do
            end do
!$acc loop collapse(2)
            do i = icii1 , icii2
              do j = jcii1 , jcii2
                p(j,i,k,n) = p(j,i,k,n) + nu * p2d(j,i)
              end do
            end do
          end do
        end do
!$acc end parallel
      end subroutine filt4d
#endif

      subroutine divergence_filter( )
        implicit none
        integer(ik4) :: j , i , k
#ifdef USE_MPI3
        type(commdata_real) :: comm

!!$acc update self(zdiv2)
        call exchange_lrbt_pre(zdiv2,1,jce1,jce2,ice1,ice2,1,kz,comm)
!!$acc update device(zdiv2)

!$acc parallel present(p2d, zdiv2, xknu)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = ici1+1 , ici2-1
            do j = jci1+1 , jci2-1
              p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                      zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                         d_half   * zdiv2(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = ici1+1 , ici2-1
            do j = jci1+1 , jci2-1
              zdiv2(j,i,k) = zdiv2(j,i,k) + xknu(k) * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel

!!$acc update self(zdiv2)
        call exchange_lrbt_post(zdiv2,1,jce1,jce2,ice1,ice2,1,kz,comm)
!!$acc update device(zdiv2)

!$acc parallel present(p2d, zdiv2, xknu)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector
          do i = ici1 , ici2
            p2d(jci1,i) = 0.125_rkx * (zdiv2(jci1-1,i,k) + zdiv2(jci1+1,i,k) + &
                                       zdiv2(jci1,i-1,k) + zdiv2(jci1,i+1,k)) -&
                         d_half   * zdiv2(jci1,i,k)
          end do
!$acc loop vector
          do i = ici1 , ici2
            zdiv2(jci1,i,k) = zdiv2(jci1,i,k) + xknu(k) * p2d(jci1,i)
          end do
!$acc loop vector
          do i = ici1 , ici2
            p2d(jci2,i) = 0.125_rkx * (zdiv2(jci2-1,i,k) + zdiv2(jci2+1,i,k) + &
                                       zdiv2(jci2,i-1,k) + zdiv2(jci2,i+1,k)) -&
                         d_half   * zdiv2(jci2,i,k)
          end do
!$acc loop vector
          do i = ici1 , ici2
            zdiv2(jci2,i,k) = zdiv2(jci2,i,k) + xknu(k) * p2d(jci2,i)
          end do
!$acc loop vector
          do j = jci1+1 , jci2-1
            p2d(j,ici1) = 0.125_rkx * (zdiv2(j-1,ici1,k) + zdiv2(j+1,ici1,k) + &
                                       zdiv2(j,ici1-1,k) + zdiv2(j,ici1+1,k)) -&
                       d_half   * zdiv2(j,ici1,k)
          end do
!$acc loop vector
          do j = jci1+1 , jci2-1
            zdiv2(j,ici1,k) = zdiv2(j,ici1,k) + xknu(k) * p2d(j,ici1)
          end do
!$acc loop vector
          do j = jci1+1 , jci2-1
            p2d(j,ici2) = 0.125_rkx * (zdiv2(j-1,ici2,k) + zdiv2(j+1,ici2,k) + &
                                       zdiv2(j,ici2-1,k) + zdiv2(j,ici2+1,k)) -&
                       d_half   * zdiv2(j,ici2,k)
          end do
!$acc loop vector
          do j = jci1+1 , jci2-1
            zdiv2(j,ici2,k) = zdiv2(j,ici2,k) + xknu(k) * p2d(j,ici2)
          end do
        end do
!$acc end parallel
#else
!!$acc update self(zdiv2)
        call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(zdiv2)

!$acc parallel present(p2d, zdiv2, xknu)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                      zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                         d_half   * zdiv2(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              zdiv2(j,i,k) = zdiv2(j,i,k) + xknu(k) * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
#endif
      end subroutine divergence_filter

      subroutine filtpai
        implicit none
        integer(ik4) :: j , i , k

!!$acc update self(pai)
        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(pai)

!$acc parallel present(p2d, pai) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (pai(j-1,i,k) + pai(j+1,i,k) + &
                                      pai(j,i-1,k) + pai(j,i+1,k)) - &
                         d_half   * pai(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              pai(j,i,k) = pai(j,i,k) + nupaitq * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
      end subroutine filtpai

      subroutine filttheta
        implicit none
        integer(ik4) :: j , i , k

!!$acc update self(tetav)
        call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(tetav)

!$acc parallel present(p2d, tetav) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (tetav(j-1,i,k) + tetav(j+1,i,k) + &
                                      tetav(j,i-1,k) + tetav(j,i+1,k)) - &
                         d_half   * tetav(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              tetav(j,i,k) = tetav(j,i,k) + nupaitq * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
      end subroutine filttheta

      subroutine filtqv
        implicit none
        integer(ik4) :: j , i , k

!!$acc update self(qv)
        call exchange_lrbt(qv,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(qv)

!$acc parallel present(p2d, qv) if(on_device)
!$acc loop gang
        do k = 1 , kz
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2d(j,i) = 0.125_rkx * (qv(j-1,i,k) + qv(j+1,i,k) + &
                                      qv(j,i-1,k) + qv(j,i+1,k)) - &
                         d_half   * qv(j,i,k)
            end do
          end do
!$acc loop vector collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              qv(j,i,k) = qv(j,i,k) + nupaitq * p2d(j,i)
            end do
          end do
        end do
!$acc end parallel
      end subroutine filtqv

      subroutine sound(dts)
        implicit none
        real(rkx) , intent(in) :: dts
        integer(ik4) :: i , j , k , nsound
        real(rkx) :: zuh , zvh , zcx , zcy
        real(rkx) :: zrfmzum , zrfmzup , zrfmzvm , zrfmzvp
        real(rkx) :: zup , zum , zvp , zvm , zqs , zdth
        real(rkx) :: zrom1w , zwexpl , zu , zd , zrapp
        real(rkx) :: zfz , zcor1u , zcor1v
        real(rkx) :: zrom1u , zrom1v
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz , zcs2
#ifdef USE_MPI3
        type(commdata_real) :: comm1, comm2 , comm3
#endif

        zdtrdx = dts/dx
        zdtrdy = dts/dx
        zdtrdz = dts/dzita
        zcs2 = zdtrdz**2*rdrcv

        !  sound waves

#ifndef USE_MPI3
        if ( .not. do_fulleq ) then
!!$acc update self(tetav)
          call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(tetav)
        end if
#endif

        do nsound = 1 , mo_nsound

#ifdef USE_MPI3
          if ( nsound == 1 .and. .not. do_fulleq ) then
!!$acc update self(tetav)
            call exchange_lrbt_pre(tetav,1,jce1,jce2,ice1,ice2,1,kz,comm3)
          end if
!!$acc update self(u)
          call exchange_lrbt_pre(u,1,jde1,jde2,ice1,ice2,1,kz,comm1)
!!$acc update self(v)
          call exchange_lrbt_pre(v,1,jce1,jce2,ide1,ide2,1,kz,comm2)

          ! partial definition of the generalized vertical velocity

!$acc parallel present(u, v, w, hx, hy) private(zuh, zvh)
!$acc loop collapse(2)
          do i = ici1 , ici2-1
            do j = jci1 , jci2-1
              zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
              zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
              w(j,i,kzp1) = d_half * (zuh+zvh)
            end do
          end do
!$acc end parallel
!$acc parallel present(w, s)
!$acc loop collapse(2)
          do i = ici1 , ici2-1
            do j = jci1 , jci2-1
              s(j,i,kzp1) = -w(j,i,kzp1)
            end do
          end do
!$acc end parallel

          ! Equation 10, generalized vertical velocity

!$acc parallel present(u, hx, v, hy, s, gzitak) private(zuh, zvh)
!$acc loop collapse(3)
!!$acc loop tile(4,8,8)
          do k = kz , 2 , -1
            do i = ici1 , ici2-1
              do j = jci1 , jci2-1
                zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) + &
                      (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
                zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) + &
                      (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
                s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
              end do
            end do
          end do
!$acc end parallel

          ! Part of divergence (except w contribution) put in zdiv2
          ! Equation 16

          if ( lrotllr ) then
!$acc parallel present(fmz, u, v, rmv, zdiv2, mx) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1 , ici2-1
                do j = jci1 , jci2-1
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
!$acc end parallel
          else
!$acc parallel present(fmz, u, v, rmv, rmu, zdiv2, mx2) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1 , ici2-1
                do j = jci1 , jci2-1
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
!$acc end parallel
          end if

          call exchange_lrbt_post(u,1,jde1,jde2,ice1,ice2,1,kz,comm1)
          call exchange_lrbt_post(v,1,jce1,jce2,ide1,ide2,1,kz,comm2)
!!$acc update device(u,v)

!$acc parallel present(u, v, w, hx, hy) private(zuh, zvh)
!$acc loop vector
          do i = ici1 , ici2-1
            zuh = u(jci2,i,kz) * hx(jci2,i) + u(jci2+1,i,kz) * hx(jci2+1,i)
            zvh = v(jci2,i,kz) * hy(jci2,i) + v(jci2,i+1,kz) * hy(jci2,i+1)
            w(jci2,i,kzp1) = d_half * (zuh+zvh)
          end do
!$acc end parallel
!$acc parallel present(w, s)
!$acc loop vector
          do i = ici1 , ici2-1
            s(jci2,i,kzp1) = -w(jci2,i,kzp1)
          end do
!$acc end parallel
!$acc parallel present(u, hx, v, hy, s, gzitak) private(zuh, zvh)
!$acc loop collapse(2)
          do k = kz , 2 , -1
            do i = ici1 , ici2-1
              zuh = (u(jci2,i,k)   + u(jci2,i,k-1))   * hx(jci2,i) + &
                    (u(jci2+1,i,k) + u(jci2+1,i,k-1)) * hx(jci2+1,i)
              zvh = (v(jci2,i,k)   + v(jci2,i,k-1))   * hy(jci2,i) + &
                    (v(jci2,i+1,k) + v(jci2,i+1,k-1)) * hy(jci2,i+1)
              s(jci2,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
            end do
          end do
!$acc end parallel
          if ( lrotllr ) then
!$acc parallel present(fmz, u, v, rmv, zdiv2, mx) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(2)
            do k = 1 , kz
              do i = ici1 , ici2-1
                zrfmzum = d_two / (fmz(jci2,i,k) + fmz(jci2-1,i,k))
                zrfmzvm = d_two / (fmz(jci2,i,k) + fmz(jci2,i-1,k))
                zrfmzup = d_two / (fmz(jci2,i,k) + fmz(jci2+1,i,k))
                zrfmzvp = d_two / (fmz(jci2,i,k) + fmz(jci2,i+1,k))
                zum = zdtrdx * u(jci2,i,k) * zrfmzum
                zup = zdtrdx * u(jci2+1,i,k) * zrfmzup
                zvm = zdtrdy * v(jci2,i,k) * zrfmzvm * rmv(jci2,i)
                zvp = zdtrdy * v(jci2,i+1,k) * zrfmzvp * rmv(jci2,i+1)
                zdiv2(jci2,i,k) = fmz(jci2,i,k) * &
                  mx(jci2,i) * ((zup-zum) + (zvp-zvm))
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(fmz, u, v, rmv, rmu, zdiv2, mx2) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(2)
            do k = 1 , kz
              do i = ici1 , ici2-1
                zrfmzum = d_two / (fmz(jci2,i,k) + fmz(jci2-1,i,k))
                zrfmzvm = d_two / (fmz(jci2,i,k) + fmz(jci2,i-1,k))
                zrfmzup = d_two / (fmz(jci2,i,k) + fmz(jci2+1,i,k))
                zrfmzvp = d_two / (fmz(jci2,i,k) + fmz(jci2,i+1,k))
                zum = zdtrdx * u(jci2,i,k)   * rmu(jci2,i)   * zrfmzum
                zup = zdtrdx * u(jci2+1,i,k) * rmu(jci2+1,i) * zrfmzup
                zvm = zdtrdy * v(jci2,i,k)   * rmv(jci2,i)   * zrfmzvm
                zvp = zdtrdy * v(jci2,i+1,k) * rmv(jci2,i+1) * zrfmzvp
                zdiv2(jci2,i,k) = mx2(jci2,i) * &
                  fmz(jci2,i,k) * ((zup-zum)+(zvp-zvm))
              end do
            end do
!$acc end parallel
          end if

!$acc parallel present(u, v, w, hx, hy) private(zuh, zvh)
!$acc loop vector
          do j = jci1 , jci2
            zuh = u(j,ici2,kz) * hx(j,ici2) + u(j+1,ici2,kz) * hx(j+1,ici2)
            zvh = v(j,ici2,kz) * hy(j,ici2) + v(j,ici2+1,kz) * hy(j,ici2+1)
            w(j,ici2,kzp1) = d_half * (zuh+zvh)
          end do
!$acc end parallel
!$acc parallel present(w, s)
!$acc loop vector
          do j = jci1 , jci2
            s(j,ici2,kzp1) = -w(j,ici2,kzp1)
          end do
!$acc end parallel
!$acc parallel present(u, hx, v, hy, s, gzitak) private(zuh, zvh)
!$acc loop collapse(2)
          do k = kz , 2 , -1
            do j = jci1 , jci2
              zuh = (u(j,ici2,k)   + u(j,ici2,k-1))   * hx(j,ici2) + &
                    (u(j+1,ici2,k) + u(j+1,ici2,k-1)) * hx(j+1,ici2)
              zvh = (v(j,ici2,k)   + v(j,ici2,k-1))   * hy(j,ici2) + &
                    (v(j,ici2+1,k) + v(j,ici2+1,k-1)) * hy(j,ici2+1)
              s(j,ici2,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
            end do
          end do
!$acc end parallel
          if ( lrotllr ) then
!$acc parallel present(fmz, u, v, rmv, zdiv2, mx) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jci1 , jci2
                zrfmzum = d_two / (fmz(j,ici2,k) + fmz(j-1,ici2,k))
                zrfmzvm = d_two / (fmz(j,ici2,k) + fmz(j,ici2-1,k))
                zrfmzup = d_two / (fmz(j,ici2,k) + fmz(j+1,ici2,k))
                zrfmzvp = d_two / (fmz(j,ici2,k) + fmz(j,ici2+1,k))
                zum = zdtrdx * u(j,ici2,k) * zrfmzum
                zup = zdtrdx * u(j+1,ici2,k) * zrfmzup
                zvm = zdtrdy * v(j,ici2,k) * zrfmzvm * rmv(j,ici2)
                zvp = zdtrdy * v(j,ici2+1,k) * zrfmzvp * rmv(j,ici2+1)
                zdiv2(j,ici2,k) = fmz(j,ici2,k) * &
                  mx(j,ici2) * ((zup-zum) + (zvp-zvm))
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(fmz, u, v, rmv, rmu, zdiv2, mx2) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jci1 , jci2
                zrfmzum = d_two / (fmz(j,ici2,k) + fmz(j-1,ici2,k))
                zrfmzvm = d_two / (fmz(j,ici2,k) + fmz(j,ici2-1,k))
                zrfmzup = d_two / (fmz(j,ici2,k) + fmz(j+1,ici2,k))
                zrfmzvp = d_two / (fmz(j,ici2,k) + fmz(j,ici2+1,k))
                zum = zdtrdx * u(j,ici2,k)   * rmu(j,ici2)   * zrfmzum
                zup = zdtrdx * u(j+1,ici2,k) * rmu(j+1,ici2) * zrfmzup
                zvm = zdtrdy * v(j,ici2,k)   * rmv(j,ici2)   * zrfmzvm
                zvp = zdtrdy * v(j,ici2+1,k) * rmv(j,ici2+1) * zrfmzvp
                zdiv2(j,ici2,k) = mx2(j,ici2) * &
                  fmz(j,ici2,k) * ((zup-zum)+(zvp-zvm))
              end do
            end do
!$acc end parallel
          end if
#else
!!$acc update self(u,v)
          call exchange_lrbt(u,1,jde1,jde2,ice1,ice2,1,kz)
          call exchange_lrbt(v,1,jce1,jce2,ide1,ide2,1,kz)
!!$acc update device(u,v)

          ! partial definition of the generalized vertical velocity

!$acc parallel present(u, v, w, hx, hy) private(zuh, zvh)
!$acc loop collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
              zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
              w(j,i,kzp1) = d_half * (zuh+zvh)
            end do
          end do
!$acc end parallel

!$acc parallel present(s, w)
!$acc loop collapse(2)
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,kzp1) = -w(j,i,kzp1)
            end do
          end do
!$acc end parallel

          ! Equation 10, generalized vertical velocity

!$acc parallel present(u, hx, v, hy, s, gzitak) private(zuh, zvh)
!$acc loop collapse(3)
!!$acc loop tile(4,8,8)
#ifdef OPENACC
          do k = 2 , kz
#else
          do k = kz , 2 , -1
#endif
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
!$acc end parallel

          ! Part of divergence (except w contribution) put in zdiv2
          ! Equation 16

          if ( lrotllr ) then
!$acc parallel present(fmz, u, v, rmv, zdiv2, mx) private(zrfmzum,&
!$acc& zrfmzvm, zrfmzup, zrfmzvp, zum, zup, zvm, zvp)
!$acc loop collapse(3)
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
!$acc end parallel
          else
!$acc parallel present(fmz, u, rmu, v, rmv, zdiv2, mx2) private(zum,&
!$acc& zup, zvm, zvp, zrfmzum, zrfmzvm, zrfmzup, zrfmzvp)
!$acc loop collapse(3)
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
!$acc end parallel
          end if

#endif

#ifdef USE_MPI3
          if ( nsound == 1 .and. .not. do_fulleq ) then
            call exchange_lrbt_post(tetav,1,jce1,jce2,ice1,ice2,1,kz,comm3)
!!$acc update device(tetav)
          end if
#endif
          if ( do_divdamp ) then
            do k = 1 , kz
              zprof(k) = (ddamp + (1.0_rkx-ddamp)*xknu(k)) * &
                      0.125_rkx*(dx**2)/dtsound
            end do
            call divdamp
          end if

!$acc parallel present(zdiv2, fmz, s)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv2(j,i,k) = zdiv2(j,i,k) + fmz(j,i,k) * &
                       zdtrdz * (s(j,i,k) - s(j,i,k+1))
              end do
            end do
          end do
!$acc end parallel

          ! new w (implicit scheme) from Equation 19

!$acc parallel present(deltaw, w, fmzf, tetav, qv, qsat, t, pai,&
!$acc& zdiv2, fmz, ffilt, wwkw) private(zrom1w, zqs, zdth, zwexpl,&
!$acc& zu, zd, zrapp)
!$acc loop collapse(2)
#ifdef OPENACC
          do i = ici1 , ici2
            do j = jci1 , jci2
!$acc loop seq
              do k = kz , 2 , -1
#else
          do k = kz , 2 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
#endif
                deltaw(j,i,k) = -w(j,i,k)
                ! explicit w:
                !    it must be consistent with the initialization of pai
                zrom1w = d_half * cpd * fmzf(j,i,k) * &
                        (tetav(j,i,k-1) + tetav(j,i,k))
                zrom1w = zrom1w - cpd * w(j,i,k) * &
                         fmzf(j,i,k)*fmzf(j,i,k) * &
                         real(nsound,rkx) * zdtrdz * &
                         (tetav(j,i,k-1) - tetav(j,i,k)) !! GW
                if ( qv(j,i,k) > 0.96_rkx*qsat(j,i,k) .and. &
                     w(j,i,k) > 0.1_rkx ) then
                  zqs = d_half*(qsat(j,i,k)+qsat(j,i,k-1))
                  zdth = egrav*w(j,i,k)*real(nsound-1,rkx)*dts*wlhv*wlhv* &
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
#ifdef OPENACC
!$acc loop seq
              do k = 2 , kz
                w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
                deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
              end do
#endif
            end do
          end do
!$acc end parallel

          ! 2nd loop for the tridiagonal inversion
#ifndef OPENACC
          do k = 2 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
                deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
              end do
            end do
          end do
#endif

          ! new Exner function (Equation 19)

!$acc parallel present(zdiv2, fmz, w)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv2(j,i,k) = zdiv2(j,i,k) + zdtrdz * fmz(j,i,k) * &
                      (w(j,i,k) - w(j,i,k+1))
              end do
            end do
          end do
!$acc end parallel

          if ( do_fulleq ) then
            if ( ipptls > 0 ) then
!$acc parallel present(zdiv2, qv, qc, tetav, qwltot, qwitot)
!$acc loop collapse(3)
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
!$acc end parallel
            end if
#ifdef USE_MPI3
!!$acc update self(tetav)
            call exchange_lrbt_pre(tetav,1,jce1,jce2,ice1,ice2,1,kz,comm3)
#else
!!$acc update self(tetav)
            call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(tetav)
#endif
          end if

          if ( do_filterdiv ) call divergence_filter( )

          ! horizontal momentum equations
!$acc parallel present(pai, zdiv2)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv2(j,i,k))
              end do
            end do
          end do
!$acc end parallel

!!$acc update self(u,v)
!$acc parallel present(u, ud)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ice1ga, ice2ga
              do j = jde1ga , jde2ga
                ud(j,i,k) = u(j,i,k)
              end do
            end do
          end do
!$acc end parallel
!$acc parallel present(v, vd)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ide1ga, ide2ga
              do j = jce1ga , jce2ga
                vd(j,i,k) = v(j,i,k)
              end do
            end do
          end do
!$acc end parallel

#ifdef USE_MPI3
!!$acc update self(pai, deltaw)
          call exchange_lrbt_pre(pai,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call exchange_lrbt_pre(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1,comm2)

          if ( do_fulleq ) then
            call exchange_lrbt_post(tetav,1,jce1,jce2,ice1,ice2,1,kz,comm3)
!!$acc update device(tetav)
          end if

          if ( lrotllr ) then
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, gzitakh,&
!$acc& pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1+1 , ici2
                do j = jdi1+1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + &
                        vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
!$acc end parallel
            zcy = zdtrdy
!$acc parallel present(deltaw, tetav, corv, ud, v, hy, gzitakh,&
!$acc& pai) private(zfz, zrom1v, zcor1v)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = idi1+1 , idi2
                do j = jci1+1 , jci2
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + &
                        ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, &
!$acc&   gzitakh, pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1+1 , ici2
                do j = jdi1+1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + &
                        vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
!$acc end parallel
!$acc parallel present(mv, deltaw, tetav, corv, ud, v, hy,&
!$acc&    gzitakh, pai) private(zcy, zfz, zrom1v, zcor1v)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = idi1+1 , idi2
                do j = jci1+1 , jci2
                  zcy = zdtrdy * mv(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + &
                        ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do
          end if
!$acc end parallel

          call exchange_lrbt_post(pai,1,jce1,jce2,ice1,ice2,1,kz,comm1)
          call exchange_lrbt_post(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1,comm2)
!!$acc update device(pai, deltaw)

          if ( lrotllr ) then
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, gzitakh,&
!$acc& pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(2)
            do k = 1 , kz
              do i = ici1+1 , ici2
                zcx = zdtrdx * mu(jdi1,i)
                zfz = 0.25_rkx * &
                  (deltaw(jdi1-1,i,k) + deltaw(jdi1-1,i,k+1) + &
                   deltaw(jdi1,i,k)   + deltaw(jdi1,i,k+1)) + egrav * dts
                zrom1u = d_half * cpd * (tetav(jdi1-1,i,k) + tetav(jdi1,i,k))
                zcor1u = coru(jdi1,i) * dts * 0.25_rkx * &
                     (vd(jdi1,i,k) + vd(jdi1-1,i,k) + &
                      vd(jdi1-1,i+1,k) + vd(jdi1,i+1,k))
                ! Equation 17
                u(jdi1,i,k) = u(jdi1,i,k) + zcor1u - &
                           zfz * hx(jdi1,i) * gzitakh(k) - &
                           zcx * zrom1u * (pai(jdi1,i,k) - pai(jdi1-1,i,k))
              end do
            end do
!$acc end parallel
!$acc parallel present(deltaw, tetav, corv, ud, v, hy, gzitakh,&
!$acc& pai) private(zfz, zrom1v, zcor1v)
!$acc loop collapse(2)
            zcy = zdtrdy
            do k = 1 , kz
              do i = idi1+1 , idi2
                zfz = 0.25_rkx * &
                  (deltaw(jci1,i-1,k) + deltaw(jci1,i-1,k+1) + &
                   deltaw(jci1,i,k)   + deltaw(jci1,i,k+1)) + egrav * dts
                zrom1v = d_half * cpd * (tetav(jci1,i-1,k) + tetav(jci1,i,k))
                zcor1v = corv(jci1,i) * dts * 0.25_rkx * &
                     (ud(jci1,i,k) + ud(jci1,i-1,k) + &
                      ud(jci1+1,i,k) + ud(jci1+1,i-1,k))
                ! Equation 18
                v(jci1,i,k) = v(jci1,i,k) - zcor1v - &
                           zfz * hy(jci1,i) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(jci1,i,k) - pai(jci1,i-1,k))
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, &
!$acc&   gzitakh, pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(2)
            do k = 1 , kz
              do i = ici1+1 , ici2
                zcx = zdtrdx * mu(jdi1,i)
                zfz = 0.25_rkx * &
                  (deltaw(jdi1-1,i,k) + deltaw(jdi1-1,i,k+1) + &
                   deltaw(jdi1,i,k)   + deltaw(jdi1,i,k+1)) + egrav * dts
                zrom1u = d_half * cpd * (tetav(jdi1-1,i,k) + tetav(jdi1,i,k))
                zcor1u = coru(jdi1,i) * dts * 0.25_rkx * &
                     (vd(jdi1,i,k) + vd(jdi1-1,i,k) + &
                      vd(jdi1-1,i+1,k) + vd(jdi1,i+1,k))
                ! Equation 17
                u(jdi1,i,k) = u(jdi1,i,k) + zcor1u - &
                           zfz * hx(jdi1,i) * gzitakh(k) - &
                           zcx * zrom1u * (pai(jdi1,i,k) - pai(jdi1-1,i,k))
              end do
            end do
!$acc end parallel
!$acc parallel present(mv, deltaw, tetav, corv, ud, v, hy,&
!$acc&    gzitakh, pai) private(zcy, zfz, zrom1v, zcor1v)
!$acc loop collapse(2)
            do k = 1 , kz
              do i = idi1+1 , idi2
                zcy = zdtrdy * mv(jci1,i)
                zfz = 0.25_rkx * &
                  (deltaw(jci1,i-1,k) + deltaw(jci1,i-1,k+1) + &
                   deltaw(jci1,i,k)   + deltaw(jci1,i,k+1)) + egrav * dts
                zrom1v = d_half * cpd * (tetav(jci1,i-1,k) + tetav(jci1,i,k))
                zcor1v = corv(jci1,i) * dts * 0.25_rkx * &
                     (ud(jci1,i,k) + ud(jci1,i-1,k) + &
                      ud(jci1+1,i,k) + ud(jci1+1,i-1,k))
                ! Equation 18
                v(jci1,i,k) = v(jci1,i,k) - zcor1v - &
                           zfz * hy(jci1,i) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(jci1,i,k) - pai(jci1,i-1,k))
              end do
            end do
!$acc end parallel
          end if

          if ( lrotllr ) then
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, gzitakh,&
!$acc& pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jdi1 , jdi2
                zcx = zdtrdx * mu(j,ici1)
                zfz = 0.25_rkx * &
                  (deltaw(j-1,ici1,k) + deltaw(j-1,ici1,k+1) + &
                   deltaw(j,ici1,k)   + deltaw(j,ici1,k+1)) + egrav * dts
                zrom1u = d_half * cpd * (tetav(j-1,ici1,k) + tetav(j,ici1,k))
                zcor1u = coru(j,ici1) * dts * 0.25_rkx * &
                     (vd(j,ici1,k) + vd(j-1,ici1,k) + &
                      vd(j-1,ici1+1,k) + vd(j,ici1+1,k))
                ! Equation 17
                u(j,ici1,k) = u(j,ici1,k) + zcor1u - &
                           zfz * hx(j,ici1) * gzitakh(k) - &
                           zcx * zrom1u * (pai(j,ici1,k) - pai(j-1,ici1,k))
              end do
            end do
!$acc end parallel
            zcy = zdtrdy
!$acc parallel present(deltaw, tetav, corv, ud, v, hy, gzitakh,&
!$acc& pai) private(zfz, zrom1v, zcor1v)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jci1 , jci2
                zfz = 0.25_rkx * &
                  (deltaw(j,idi1-1,k) + deltaw(j,idi1-1,k+1) + &
                   deltaw(j,idi1,k)   + deltaw(j,idi1,k+1)) + egrav * dts
                zrom1v = d_half * cpd * (tetav(j,idi1-1,k) + tetav(j,idi1,k))
                zcor1v = corv(j,idi1) * dts * 0.25_rkx * &
                     (ud(j,idi1,k) + ud(j,idi1-1,k) + &
                      ud(j+1,idi1,k) + ud(j+1,idi1-1,k))
                ! Equation 18
                v(j,idi1,k) = v(j,idi1,k) - zcor1v - &
                           zfz * hy(j,idi1) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(j,idi1,k) - pai(j,idi1-1,k))
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, &
!$acc&   gzitakh, pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jdi1 , jdi2
                zcx = zdtrdx * mu(j,ici1)
                zfz = 0.25_rkx * &
                  (deltaw(j-1,ici1,k) + deltaw(j-1,ici1,k+1) + &
                   deltaw(j,ici1,k)   + deltaw(j,ici1,k+1)) + egrav * dts
                zrom1u = d_half * cpd * (tetav(j-1,ici1,k) + tetav(j,ici1,k))
                zcor1u = coru(j,ici1) * dts * 0.25_rkx * &
                     (vd(j,ici1,k) + vd(j-1,ici1,k) + &
                      vd(j-1,ici1+1,k) + vd(j,ici1+1,k))
                ! Equation 17
                u(j,ici1,k) = u(j,ici1,k) + zcor1u - &
                           zfz * hx(j,ici1) * gzitakh(k) - &
                           zcx * zrom1u * (pai(j,ici1,k) - pai(j-1,ici1,k))
              end do
            end do
!$acc end parallel
!$acc parallel present(mv, deltaw, tetav, corv, ud, v, hy,&
!$acc&    gzitakh, pai) private(zcy, zfz, zrom1v, zcor1v)
!$acc loop collapse(2)
            do k = 1 , kz
              do j = jci1 , jci2
                zcy = zdtrdy * mv(j,idi1)
                zfz = 0.25_rkx * &
                  (deltaw(j,idi1-1,k) + deltaw(j,idi1-1,k+1) + &
                   deltaw(j,idi1,k)   + deltaw(j,idi1,k+1)) + egrav * dts
                zrom1v = d_half * cpd * (tetav(j,idi1-1,k) + tetav(j,idi1,k))
                zcor1v = corv(j,idi1) * dts * 0.25_rkx * &
                     (ud(j,idi1,k) + ud(j,idi1-1,k) + &
                      ud(j+1,idi1,k) + ud(j+1,idi1-1,k))
                ! Equation 18
                v(j,idi1,k) = v(j,idi1,k) - zcor1v - &
                           zfz * hy(j,idi1) * gzitakh(k) -  &
                           zcy * zrom1v * (pai(j,idi1,k) - pai(j,idi1-1,k))
              end do
            end do
!$acc end parallel
          end if
#else
!!$acc update self(pai, deltaw)
          call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
          call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1)
!!$acc update device(pai, deltaw)

          if ( lrotllr ) then
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, gzitakh,&
!$acc& pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jdi1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + &
                        vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
!$acc end parallel
            zcy = zdtrdy
!$acc parallel present(deltaw, tetav, corv, ud, v, hy, gzitakh,&
!$acc& pai) private(zfz, zrom1v, zcor1v)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = idi1 , idi2
                do j = jci1 , jci2
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + &
                        ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do
!$acc end parallel
          else
!$acc parallel present(mu, deltaw, tetav, coru, vd, u, hx, &
!$acc&   gzitakh, pai) private(zcx, zfz, zrom1u, zcor1u)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jdi1 , jdi2
                  zcx = zdtrdx * mu(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
                  zcor1u = coru(j,i) * dts * 0.25_rkx * &
                       (vd(j,i,k) + vd(j-1,i,k) + &
                        vd(j-1,i+1,k) + vd(j,i+1,k))
                  ! Equation 17
                  u(j,i,k) = u(j,i,k) + zcor1u - &
                             zfz * hx(j,i) * gzitakh(k) - &
                             zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
                end do
              end do
            end do
!$acc end parallel
!$acc parallel present(mv, deltaw, tetav, corv, ud, v, hy,&
!$acc&    gzitakh, pai) private(zcy, zfz, zrom1v, zcor1v)
!$acc loop collapse(3)
            do k = 1 , kz
              do i = idi1 , idi2
                do j = jci1 , jci2
                  zcy = zdtrdy * mv(j,i)
                  zfz = 0.25_rkx * &
                    (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                     deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav * dts
                  zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
                  zcor1v = corv(j,i) * dts * 0.25_rkx * &
                       (ud(j,i,k) + ud(j,i-1,k) + &
                        ud(j+1,i,k) + ud(j+1,i-1,k))
                  ! Equation 18
                  v(j,i,k) = v(j,i,k) - zcor1v - &
                             zfz * hy(j,i) * gzitakh(k) -  &
                             zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
                end do
              end do
            end do
!$acc end parallel
          end if
#endif

        end do ! sound loop

        ! complete computation of generalized vertical velocity
        ! Complete Equation 10
!$acc parallel present(s, w, fmzf)
!$acc loop collapse(3)
        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
            end do
          end do
        end do
!$acc end parallel
!$acc parallel present(s)
!$acc loop collapse(2)
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,1) = d_zero
          end do
        end do
!$acc end parallel
!$acc parallel present(s)
!$acc loop collapse(2)
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,kzp1) = d_zero
          end do
        end do
!$acc end parallel

      end subroutine sound

      subroutine advection(dta)
        implicit none
        integer(ik4) :: n
        real(rkx) :: dta
        real(rkx) , pointer , dimension(:,:,:) :: ptr

!!$acc update device(u,v)
        call uvstagtox(u,v,ux,vx)

        ! Compute W (and TKE if required) on zita levels

!!$acc update device(w)
        call wstagtox(w,wx)

        if ( ibltyp == 2 ) then
!!$acc update device(tke)
          call wstagtox(tke,tkex)
        end if

        call wafone(tetav,dta)
        call wafone(pai,dta)
        call wafone(ux,dta)
        call wafone(vx,dta)
        call wafone(wx,dta)
        call wafone(qv,dta,pmin=minqq)
        if ( ipptls > 0 ) then
          do n = iqfrst , iqlst
            call assignpnt(qx,ptr,n)
            call wafone(ptr,dta,pfac=1.0e4_rkx,pmin=epsilon(pmin))
          end do
        end if
        if ( ibltyp == 2 ) then
          call wafone(tkex,dta)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(trac,ptr,n)
            call wafone(ptr,dta,pfac=1.0e8_rkx,pmin=d_zero)
          end do
        end if

        ! Interpolate on staggered points
        call xtouvstag(ux,vx,u,v)

        ! Back to half-levels
!!$acc update device(wx)
        call xtowstag(wx,w)
        if ( ibltyp == 2 ) then
!!$acc update device(tkex)
          call xtowstag(tkex,tke)
        end if
      end subroutine advection

      pure real(rkx) function rdeno(t1,t2,t3,t4)
!$acc routine seq
        implicit none
        real(rkx) , intent(in) :: t1 , t2 , t3 , t4
        real(rkx) :: zzden
        zzden = (t3-t4)
        rdeno = (t1-t2)/sign(max(abs(zzden),minden),zzden)
      end function rdeno

      subroutine wafone(pp,dta,pfac,pmin)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
        real(rkx) , intent(in) :: dta
        real(rkx) , optional , intent(in) :: pfac , pmin
        integer(ik4) :: j , i , k
        integer(ik4) :: k1 , k1p1 , ih , ihm1 , jh , jhm1
        real(rkx) :: zamu , r , b , zphi , is , zdv , zrfmu , zrfmd
        real(rkx) :: zrfmn , zrfmw , zrfme , zrfms
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz
        real(rkx) :: zhxvtn , zhxvts , zcostx
        real(rkx) :: pfm
#ifdef OPENACC
        integer(ik4) :: k1m1 , k1p1m1
        real(rkx) :: rm1 , bm1 , ism1
        real(rkx) :: zpbys , zpbws
        real(rkx) :: wfwk , wfwkp1
        real(rkx) :: zpbysp1 , zpbwsp1
        real(rkx) :: zamum1 , zphim1
#endif
        real(rkx) , parameter :: wlow  = 0.0_rkx
        real(rkx) , parameter :: whigh = 2.0_rkx

        zdtrdx = dta/dx
        zdtrdy = dta/dx
        zdtrdz = dta/dzita
        if ( do_vadvtwice ) then
          zdtrdz = 0.5_rkx * zdtrdz
        end if

        pfm = 0.0_rkx
        if ( present(pmin) ) then
          pfm = pmin
        end if
        if ( present(pfac) ) then
!$acc parallel present(pp)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                pp(j,i,k) = pp(j,i,k) * pfac
              end do
            end do
          end do
!$acc end parallel
          if ( present(pmin) ) then
            pfm = pfm*pfac
          end if
        end if

        ! Vertical advection
!$acc parallel present(wfw) private(j)
        do j = jci1 , jci2
          wfw(j,1) = d_zero
          wfw(j,kzp1) = d_zero
        end do
!$acc end parallel

        if ( ma%has_bdybottom ) then
!$acc parallel present(pp,wz)
!$acc loop collapse(2)
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice1,k) = 0.5_rkx * (pp(j,ice1,k)+pp(j,ici1,k))
            end do
          end do
!$acc end parallel
        end if
        if ( ma%has_bdytop ) then
!$acc parallel present(pp,wz)
!$acc loop collapse(2)
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice2,k) = 0.5_rkx * (pp(j,ice2,k)+pp(j,ici2,k))
            end do
          end do
!$acc end parallel
        end if

!$acc parallel present(wfw, pp, s, fmzf, fmz, wz) private(zamu,&
!$acc&     is, k1, k1p1, r, b, zphi, wfwkp1, wfwk, zrfmu, zrfmd, zdv)
!$acc loop collapse(3)
#ifdef OPENACC
        do k = 1 , kz
          do i = ici1 , ici2
#else
        do i = ici1 , ici2
          do k = 1 , kzm1
#endif
            do j = jci1 , jci2
              if ( k < kz ) then
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
#ifdef OPENACC
                wfwkp1 = d_half * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                                (d_one-zphi)*pp(j,i,k))
              else
                wfwkp1 = d_zero
#else
                wfw(j,k+1) = d_half * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                                    (d_one-zphi)*pp(j,i,k))
#endif
              end if
              !wfw(j,k+1) = d_half * zamu * ((d_one+zphi)*pp(j,i,k+1) + &
              !                              (d_one-zphi)*pp(j,i,k))
#ifdef OPENACC
              if ( k > 1 ) then
                zamum1 = s(j,i,k) * zdtrdz
                if ( zamum1 >= d_zero ) then
                  is = d_one
                  k1 = k
                  k1p1 = k1 + 1
                  if ( k1p1 > kz ) k1p1 = kz
                else
                  is = -d_one
                  k1 = k - 2
                  k1p1 = k - 1
                  if ( k1 < 1 ) k1 = 1
                end if
                r = rdeno(pp(j,i,k1),pp(j,i,k1p1),pp(j,i,k-1),pp(j,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu * b - is * b
                wfwk = d_half * s(j,i,k) * ((d_one+zphi)*pp(j,i,k) + &
                                            (d_one-zphi)*pp(j,i,k-1))
              else
                wfwk = d_zero
              end if
              zrfmu = zdtrdz * fmz(j,i,k)/fmzf(j,i,k)
              zrfmd = zdtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
              zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
              wz(j,i,k) = pp(j,i,k) - wfwk * zrfmu + wfwkp1 * zrfmd + zdv
#endif
            end do
          end do
#ifndef OPENACC
          do k = 1 , kz
            do j = jci1 , jci2
              zrfmu = zdtrdz * fmz(j,i,k)/fmzf(j,i,k)
              zrfmd = zdtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
              zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
              wz(j,i,k) = pp(j,i,k) - wfw(j,k)*zrfmu + wfw(j,k+1)*zrfmd + zdv
            end do
          end do
#endif
          !do k = 1 , kz
          !  do j = jci1 , jci2
          !    zdv = (s(j,i,k) - s(j,i,k+1)) * zdtrdz * pp(j,i,k)
          !    wz(j,i,k) = pp(j,i,k) - wfw(j,k) + wfw(j,k+1) + zdv
          !  end do
          !end do
        end do
!$acc end parallel

        if ( do_vadvtwice ) then
!$acc parallel present(wz, wwkw, s) private(zamu, is, k1, k1p1, r, b, zphi)
!$acc loop collapse(3)
#ifdef OPENACC
          do k = 0 , kzm1
            do i = ici1 , ici2
#else
          do i = ici1 , ici2
            do k = 1 , kzm1
#endif
              do j = jci1 , jci2
                if ( k == 0 ) then
                  wwkw(j,i,1) = d_zero
                else
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
#ifdef OPENACC
                  wwkw(j,i,k+1) = d_half*s(j,i,k+1) * &
                    ((d_one+zphi)*wz(j,i,k+1) + (d_one-zphi)*wz(j,i,k))
#else
                  wfw(j,k+1) = d_half*s(j,i,k+1) * ((d_one+zphi)*wz(j,i,k+1) + &
                                                    (d_one-zphi)*wz(j,i,k))
#endif
                end if
                !wfw(j,k+1) = d_half * zamu * ((d_one+zphi)*wz(j,i,k+1) + &
                !                              (d_one-zphi)*wz(j,i,k))
              end do
            end do
#ifndef OPENACC
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
#endif
            !do k = 1 , kz
            !  do j = jci1 , jci2
            !    zdv = (s(j,i,k) - s(j,i,k+1)) * zdtrdz * pp(j,i,k)
            !    wz(j,i,k) = pp(j,i,k) - wfw(j,k) + wfw(j,k+1) + zdv
            !  end do
            !end do
          end do
#ifdef OPENACC
!$acc end parallel
!$acc parallel present(fmzf, wz, wwkw, s, fmz) private(zrfmd, zrfmu, zdv)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zrfmu = zdtrdz * fmz(j,i,k)/fmzf(j,i,k)
                zrfmd = zdtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
                zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * wz(j,i,k)
                wz(j,i,k) = wz(j,i,k) - wwkw(j,i,k)*zrfmu + &
                                        wwkw(j,i,k+1)*zrfmd + zdv
              end do
            end do
          end do
!$acc end parallel
#endif
        end if

        if ( present(pmin) ) then
!$acc kernels present(wz)
          wz = max(wz,pfm)
!$acc end kernels
        end if

!!$acc update self(wz)
        call exchange_bt(wz,2,jci1,jci2,ice1,ice2,1,kz)
!!$acc update device(wz)

        if ( ma%has_bdyleft ) then
!$acc parallel present(p0,pp)
!$acc loop collapse(2)
          do k = 1 , kz
            do i = ici1 , ici2
              p0(jce1,i,k) = 0.5_rkx * (pp(jce1,i,k)+pp(jci1,i,k))
            end do
          end do
!$acc end parallel
        end if

        if ( ma%has_bdyright ) then
!$acc parallel present(p0,pp)
!$acc loop collapse(2)
          do k = 1 , kz
            do i = ici1 , ici2
              p0(jce2,i,k) = 0.5_rkx * (pp(jce2,i,k)+pp(jci2,i,k))
            end do
          end do
!$acc end parallel
        end if

        if ( lrotllr ) then

          ! Meridional advection
!$acc parallel present(rmv, v, fmz, wz, mx2, p0, pp) private(zamu,&
!$acc&   is, ih, ihm1, r, b, zphi, zpbys, zpbysp1, zrfmn, zrfms, zdv)
!$acc loop collapse(3)
          do k = 1 , kz
#ifdef OPENACC
            do i = ici1 , ici2
#else
            do i = ici1 , ice2ga
#endif
              do j = jci1 , jci2
                zamu = v(j,i,k) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i-1
                else
                  is = -d_one
                  ih = min(i+1,imax)
                end if
                ihm1 = max(ih-1,imin)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,i-1,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
#ifdef OPENACC
                zpbys = d_half * v(j,i,k) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))

                zamu = v(j,i+1,k) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i
                else
                  is = -d_one
                  ih = min(i+2,imax)
                end if
                ihm1 = max(ih-1,imin)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i+1,k), wz(j,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbysp1 = d_half * v(j,i+1,k) * &
                  ((d_one+zphi)*wz(j,i,k) + (d_one-zphi)*wz(j,i+1,k))
                zhxvtn = zdtrdy * rmv(j,i+1) * mx(j,i)
                zhxvts = zdtrdy * rmv(j,i) * mx(j,i)
                zrfmn = zhxvtn * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
                zrfms = zhxvts * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
                zdv = (v(j,i+1,k) * zrfmn - v(j,i,k) * zrfms) * pp(j,i,k)
                p0(j,i,k) = wz(j,i,k) + zpbys*zrfms - zpbysp1*zrfmn + zdv
#else
                zpby(j,i) = d_half * v(j,i,k) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
#endif
                !zpby(j,i) = d_half * zamu * &
                !  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
              end do
            end do
#ifndef OPENACC
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
#endif
          end do
!$acc end parallel

          if ( present(pmin) ) then
!$acc kernels present(p0)
            p0 = max(p0,pfm)
!$acc end kernels
          end if

!!$acc update self(p0)
          call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)
!!$acc update device(p0)

          ! Zonal advection

!$acc parallel present(rmu, pp, mx2, u, p0, fmz) private(zamu, is,&
!$acc&      jh, jhm1, r, b, zphi, zpbws, zpbwsp1, zrfme, zrfmw, zdv)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
#ifdef OPENACC
              do j = jci1 , jci2
#else
              do j = jci1 , jce2ga
#endif
                zamu = u(j,i,k) * mu(j,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j-1
                else
                  is = -d_one
                  jh = min(j+1,jmax)
                end if
                jhm1 = max(jh-1,jmin)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(j-1,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
#ifdef OPENACC
                zpbws = d_half * u(j,i,k) * &
                ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))

                zamu = u(j+1,i,k) * rmu(j+1,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j
                else
                  is = -d_one
                  jh = min(j+2,jmax)
                end if
                jhm1 = max(jh-1,jmin)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j+1,i,k), p0(j,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbwsp1 = d_half * u(j+1,i,k) * &
                    ((d_one+zphi)*p0(j,i,k) + (d_one-zphi)*p0(j+1,i,k))
                zcostx = zdtrdx * mu(j,i)
                zrfme = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
                zrfmw = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
                zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
                pp(j,i,k) = p0(j,i,k) + zpbws*zrfmw - zpbwsp1*zrfme + zdv
#else
                zpbw(j,i) = d_half * u(j,i,k) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
#endif
                !zpbw(j,i) = d_half * zamu * &
                !     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
              end do
            end do
#ifndef OPENACC
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
#endif
          end do
!$acc end parallel

        else

          ! Meridional advection

!$acc parallel present(rmv, v, fmz, wz, mx2, p0, pp) private(zamu, is,&
!$acc&   ih, ihm1, r, b, zphi, zpbys, zpbysp1, zrfmn, zrfms, zdv)
!$acc loop collapse(3)
          do k = 1 , kz
#ifdef OPENACC
            do i = ici1 , ici2
#else
            do i = ici1 , ice2ga
#endif
              do j = jci1 , jci2
                zamu = v(j,i,k) * rmv(j,i) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i-1
                else
                  is = -d_one
                  ih = min(i+1,imax)
                end if
                ihm1 = max(ih-1,imin)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,i-1,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
#ifdef OPENACC
                zpbys = d_half * v(j,i,k) * rmv(j,i) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))

                zamu = v(j,i+1,k) * rmv(j,i+1) * zdtrdy
                if ( zamu > d_zero ) then
                  is = d_one
                  ih = i
                else
                  is = -d_one
                  ih = min(i+2,imax)
                end if
                ihm1 = max(ih-1,imin)
                r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i+1,k), wz(j,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbysp1 = d_half * v(j,i+1,k) * rmv(j,i+1) * &
                  ((d_one+zphi)*wz(j,i,k) + (d_one-zphi)*wz(j,i+1,k))

                zrfmn = zdtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
                zrfms = zdtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
                zdv = (v(j,i+1,k) * rmv(j,i+1) * zrfmn - &
                       v(j,i,k)   * rmv(j,i)   * zrfms) * pp(j,i,k)
                p0(j,i,k) = wz(j,i,k) + &
                  mx2(j,i) * (zpbys*zrfms - zpbysp1*zrfmn + zdv)
#else
                zpby(j,i) = d_half * v(j,i,k) * rmv(j,i) * &
                  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
                !zpby(j,i) = d_half * zamu * &
                !  ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
#endif
              end do
            end do
#ifndef OPENACC
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
#endif
          end do
!$acc end parallel

          if ( present(pmin) ) then
!$acc kernels present(p0)
            p0 = max(p0,pfm)
!$acc end kernels
          end if

!!$acc update self(p0)
          call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)
!!$acc update device(p0)

          ! Zonal advection
!$acc parallel present(rmu, pp, mx2, u, p0, fmz) private(zamu, is,&
!$acc&    jh, jhm1, r, b, zphi, zpbws, zpbwsp1, zrfme, zrfmw, zdv)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ici1 , ici2
#ifdef OPENACC
              do j = jci1 , jci2
#else
              do j = jci1 , jce2ga
#endif
                zamu = u(j,i,k) * rmu(j,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j-1
                else
                  is = -d_one
                  jh = min(j+1,jmax)
                end if
                jhm1 = max(jh-1,jmin)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(j-1,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
#ifdef OPENACC
                zpbws = d_half * u(j,i,k) * rmu(j,i) * &
                ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))

                zamu = u(j+1,i,k) * rmu(j+1,i) * zdtrdx
                if ( zamu > d_zero ) then
                  is = d_one
                  jh = j
                else
                  is = -d_one
                  jh = min(j+2,jmax)
                end if
                jhm1 = max(jh-1,jmin)
                r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j+1,i,k), p0(j,i,k))
                b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
                zphi = is + zamu*b - is*b
                zpbwsp1 = d_half * u(j+1,i,k) * rmu(j+1,i) * &
                    ((d_one+zphi)*p0(j,i,k) + (d_one-zphi)*p0(j+1,i,k))
                     !zpbw(j,i) = d_half * zamu * &
                !     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
                zrfme = zdtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
                zrfmw = zdtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
                zdv = (u(j+1,i,k) * rmu(j+1,i) * zrfme - &
                       u(j,i,k)   * rmu(j,i)   * zrfmw) * pp(j,i,k)
                pp(j,i,k) = p0(j,i,k) + &
                  mx2(j,i) * (zpbws*zrfmw - zpbwsp1*zrfme + zdv)
#else
                zpbw(j,i) = d_half * u(j,i,k) * rmu(j,i) * &
                     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
                !zpbw(j,i) = d_half * zamu * &
                !     ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
#endif
              end do
            end do
#ifndef OPENACC
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
#endif
          end do
!$acc end parallel
        end if

        if ( present(pmin) ) then
!$acc parallel present(pp)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                pp(j,i,k) = max(pp(j,i,k),pfm)
              end do
            end do
          end do
!$acc end parallel
        end if
        if ( present(pfac) ) then
!$acc parallel present(pp)
!$acc loop collapse(3)
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                pp(j,i,k) = pp(j,i,k) / pfac
              end do
            end do
          end do
!$acc end parallel
        end if
      end subroutine wafone

      subroutine reset_tendencies
        implicit none
!$acc kernels present(s, deltaw, zdiv2)
        s(:,:,:) = d_zero
        deltaw(:,:,:) = d_zero
        zdiv2(:,:,:) = d_zero
!$acc end kernels
        mo_atm%tten(:,:,:) = d_zero
        mo_atm%qxten(:,:,:,:) = d_zero
        mo_atm%uten(:,:,:) = d_zero
        mo_atm%vten(:,:,:) = d_zero
        if ( ichem == 1 ) then
          mo_atm%chiten(:,:,:,:) = d_zero
        end if
        if ( ibltyp == 2 ) then
          mo_atm%tketen(:,:,:) = d_zero
        end if

        cldfra(:,:,:) = d_zero
        cldlwc(:,:,:) = d_zero

        if ( idiag > 0 ) then
          ten0 = t(jci1:jci2,ici1:ici2,1:kz)
          qen0 = qv(jci1:jci2,ici1:ici2,1:kz)
          if ( ichem == 1 ) then
            chiten0 = trac(jci1:jci2,ici1:ici2,1:kz,1:ntr)
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
          call radiation(rcmtimer%year,rcmtimer%month, &
                         rcmtimer%day,loutrad,labsem)
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
#ifdef RCEMIP
        if ( do_diffutend ) then
          call exchange_lrbt(mo_atm%tten,1,jci1,jci2,ici1,ici2,1,kz)
          call filt3d(mo_atm%tten,mo_anu2)
          call exchange_lrbt(mo_atm%qxten,1,jci1,jci2,ici1,ici2,1,kz,1,nqx)
          call filt4d(mo_atm%qxten,mo_anu2,1,nqx)
          call exchange_lrbt(mo_atm%uten,1,jdi1,jdi2,ici1,ici2,1,kz)
          call exchange_lrbt(mo_atm%vten,1,jci1,jci2,idi1,idi2,1,kz)
          call filtuv(mo_atm%uten,mo_atm%vten,mo_anu2)
        end if
#endif
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

!$acc parallel present(wx, w) if(on_device)
!$acc loop collapse(3)
    do k = 2 , kzm1
      do i = i1 , i2
        do j = j1 , j2
          wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                      0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
        end do
      end do
    end do
!$acc end parallel
!$acc parallel present(wx, w) if(on_device)
!$acc loop collapse(2)
    do i = i1 , i2
      do j = j1 , j2
        wx(j,i,1)  = d_half * (w(j,i,2)+w(j,i,1))
        wx(j,i,kz) = d_half * (w(j,i,kzp1)+w(j,i,kz))
      end do
    end do
!$acc end parallel
  end subroutine wstagtox

  ! Fully device-resident: make sure all IO is on the GPU for this subroutine
  subroutine xtowstag(wx,w)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: w
    integer(ik4) :: i , j , k

!$acc parallel present(w, wx) if(on_device)
!$acc loop collapse(3)
    do k = 3 , kzm1
      do i = ice1 , ice2
        do j = jce1 , jce2
          w(j,i,k) = 0.5625_rkx * (wx(j,i,k)  +wx(j,i,k-1)) - &
                     0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
        end do
      end do
    end do
!$acc end parallel
!$acc parallel present(w, wx) if(on_device)
!$acc loop collapse(2)
    do i = ice1 , ice2
      do j = jce1 , jce2
        w(j,i,2) = d_half * (wx(j,i,2)  +wx(j,i,1))
        w(j,i,kz) = d_half * (wx(j,i,kz)+wx(j,i,kzm1))
      end do
    end do
!$acc end parallel
  end subroutine xtowstag

  ! Fully device-resident: make sure all IO is on the GPU for this subroutine
  subroutine xtoustag(ux,u)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u
    integer(ik4) :: i , j , k

!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jdii1 , jdii2
          u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                     0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdyright ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
        end do
      end do
!$acc end parallel
    end if

    if ( ma%has_bdyleft ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
        end do
      end do
!$acc end parallel
    end if
  end subroutine xtoustag

  ! Fully device-resident: make sure all IO is on the GPU for this subroutine
  subroutine xtovstag(vx,v)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: v
    integer(ik4) :: i , j , k

!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jci1 , jci2
          v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                     0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdytop ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdybottom ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
        end do
      end do
!$acc end parallel
    end if
  end subroutine xtovstag

  ! Fully device-resident: make sure all IO is on the GPU for this subroutine
  subroutine xtouvstag(ux,vx,u,v)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    integer(ik4) :: i , j , k
#ifdef USE_MPI3
    type(commdata_real) :: comm1, comm2

!!$acc update self(ux, vx)
    call exchange_lr_pre(ux,2,jce1,jce2,ice1,ice2,1,kz,comm1)
    call exchange_bt_pre(vx,2,jce1,jce2,ice1,ice2,1,kz,comm2)

    if ( ma%has_bdyright ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdyleft ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdytop ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdybottom ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
        end do
      end do
!$acc end parallel
    end if

    call exchange_lr_post(ux,2,jce1,jce2,ice1,ice2,1,kz,comm1)
    call exchange_bt_post(vx,2,jce1,jce2,ice1,ice2,1,kz,comm2)
!!$acc update device(ux, vx)

    ! Back to wind points: U/V (fourth order)

!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jdii1 , jdii2
          u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                     0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
        end do
      end do
    end do
!$acc end parallel
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jci1 , jci2
          v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                     0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
        end do
      end do
    end do
!$acc end parallel
#else

!!$acc update self(ux, vx)
    call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(ux, vx)

    ! Back to wind points: U (fourth order)

!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jdii1 , jdii2
          u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                     0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdyright ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdyleft ) then
!$acc parallel present(u, ux) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
        end do
      end do
!$acc end parallel
    end if

    ! Back to wind points: V (fourth order)

!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jci1 , jci2
          v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                     0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdytop ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdybottom ) then
!$acc parallel present(v, vx) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
        end do
      end do
!$acc end parallel
    end if
#endif
  end subroutine xtouvstag

  ! Fully device-resident: make sure all IO is on the GPU for this subroutine
  subroutine uvstagtox(u,v,ux,vx)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    integer(ik4) :: i , j , k
#ifdef USE_MPI3
    type(commdata_real) :: comm1, comm2

!!$acc update self(u, v)
    call exchange_lr_pre(u,2,jde1,jde2,ice1,ice2,1,kz,comm1)
    call exchange_bt_pre(v,2,jce1,jce2,ide1,ide2,1,kz,comm2)

    if ( ma%has_bdyleft ) then
!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce1,i,k) = d_half * (u(jde1,i,k)+u(jdi1,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdyright ) then
!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce2,i,k) = d_half*(u(jde2,i,k) + u(jdi2,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdybottom ) then
!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice1,k) = d_half * (v(j,ide1,k)+v(j,idi1,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdytop ) then
!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice2,k) = d_half*(v(j,ide2,k) + v(j,idi2,k))
        end do
      end do
!$acc end parallel
    end if

    call exchange_lr_post(u,2,jde1,jde2,ice1,ice2,1,kz,comm1)
    call exchange_bt_post(v,2,jce1,jce2,ide1,ide2,1,kz,comm2)
!!$acc update device(u, v)

    ! Compute U-wind on T points
!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jci1 , jci2
          ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                      0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
        end do
      end do
    end do
!$acc end parallel
    ! Compute V-wind on T points
!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jce1 , jce2
          vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                      0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
        end do
      end do
    end do
!$acc end parallel
#else

!!$acc update self(u, v)
    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)
!!$acc update device(u, v)

    ! Compute U-wind on T points

!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jci1 , jci2
          ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                      0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdyleft ) then
!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce1,i,k) = d_half * (u(jde1,i,k)+u(jdi1,i,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdyright ) then
!$acc parallel present(ux, u) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce2,i,k) = d_half*(u(jde2,i,k) + u(jdi2,i,k))
        end do
      end do
!$acc end parallel
    end if

    ! Compute V-wind on T points

!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(3)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jce1 , jce2
          vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                      0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
        end do
      end do
    end do
!$acc end parallel
    if ( ma%has_bdybottom ) then
!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice1,k) = d_half * (v(j,ide1,k)+v(j,idi1,k))
        end do
      end do
!$acc end parallel
    end if
    if ( ma%has_bdytop ) then
!$acc parallel present(vx, v) if(on_device)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice2,k) = d_half*(v(j,ide2,k) + v(j,idi2,k))
        end do
      end do
!$acc end parallel
    end if
#endif
  end subroutine uvstagtox

  subroutine divdamp
    implicit none
    integer(ik4) :: i , j , k
#ifdef USE_MPI3
    type(commdata_real) :: comm
#endif

#ifdef USE_MPI3
!!$acc update self(zdiv2)
    call exchange_lrbt_pre(zdiv2,1,jce1,jce2,ice1,ice2,1,kz,comm)

    if ( lrotllr ) then
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1+1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                zprof(k)/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
!$acc end parallel
!$acc parallel present(v, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = idi1+1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                zprof(k)/dx*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
!$acc end parallel
    else
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1+1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                zprof(k)/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
!$acc end parallel
!$acc parallel present(v, rmv, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = idi1+1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                zprof(k)/(dx*rmv(j,i))*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
!$acc end parallel
    end if

    call exchange_lrbt_post(zdiv2,1,jce1,jce2,ice1,ice2,1,kz,comm)
!!$acc update device(zdiv2)

    if ( lrotllr ) then
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = u(jdi1,i,k) + &
              zprof(k)/(dx*rmu(jdi1,i))*(zdiv2(jci1,i,k)-zdiv2(jci1-1,i,k))
        end do
      end do
!$acc end parallel
!$acc parallel present(v, zdiv2, zprof)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = v(j,idi1,k) + &
              zprof(k)/dx*(zdiv2(j,ici1,k)-zdiv2(j,ici1-1,k))
        end do
      end do
!$acc end parallel
    else
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(2)
      do k = 1 , kz
        do i = ici1 , ici2
          u(jdi1,i,k) = u(jdi1,i,k) + &
              zprof(k)/(dx*rmu(jdi1,i))*(zdiv2(jci1,i,k)-zdiv2(jci1-1,i,k))
        end do
      end do
!$acc end parallel
!$acc parallel present(v, rmv, zdiv2, zprof)
!$acc loop collapse(2)
      do k = 1 , kz
        do j = jci1 , jci2
          v(j,idi1,k) = v(j,idi1,k) + &
              zprof(k)/(dx*rmv(j,idi1))*(zdiv2(j,ici1,k)-zdiv2(j,ici1-1,k))
        end do
      end do
!$acc end parallel
    end if
#else
!!$acc update self(zdiv2)
    call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
!!$acc update device(zdiv2)

    if ( lrotllr ) then
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                zprof(k)/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
!$acc end parallel
!$acc parallel present(v, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                zprof(k)/dx*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
!$acc end parallel
    else
!$acc parallel present(u, rmu, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jdi1 , jdi2
            u(j,i,k) = u(j,i,k) + &
                zprof(k)/(dx*rmu(j,i))*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
          end do
        end do
      end do
!$acc end parallel
!$acc parallel present(v, rmv, zdiv2, zprof)
!$acc loop collapse(3)
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jci1 , jci2
            v(j,i,k) = v(j,i,k) + &
                zprof(k)/(dx*rmv(j,i))*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
          end do
        end do
      end do
!$acc end parallel
    end if
#endif
  end subroutine divdamp

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
