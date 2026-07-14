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
  real(rkx), pointer, contiguous, dimension(:,:,:) :: s => null( )
  ! nonhydrostatic term in pressure gradient force
  ! tridiagonal inversion
  real(rkx), pointer, contiguous, dimension(:,:,:) :: deltaw => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: wwkw => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: tkex => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: wz => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: mx2 => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: rmx => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: rmu => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: rmv => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: p0 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: wfw => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zpby => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zpbw => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zdiv2 => null( )

  real(rkx), pointer, contiguous, dimension(:,:,:) :: ten0 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qen0 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiten0 => null( )

  real(rkx), dimension(:), pointer, contiguous :: gzitak => null( )
  real(rkx), dimension(:), pointer, contiguous :: gzitakh => null( )
  real(rkx), dimension(:), pointer, contiguous :: xknu => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: laplacian => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: xlat => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: xlon => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: coru => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: corv => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: mu => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: hx => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: mx => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: mv => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: hy => null( )
  real(rkx), dimension(:), pointer, contiguous :: rlat => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: ps => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: ts => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: t2m => null( )
  real(rkx), dimension(:,:), pointer, contiguous :: ht => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: fmz => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: rfmzu => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: rfmzv => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: fmzf => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: pai => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: tetav => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: tvirt => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: z => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: u => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: v => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: w => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: ux => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: vx => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: wx => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: ud => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: vd => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: p => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: t => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: rho => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qv => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qc => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qi => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qr => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qs => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qsat => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qwltot => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qwitot => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: tke => null( )
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: qx => null( )
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: trac => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: uten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: vten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: tten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qvten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qcten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: qiten => null( )
  real(rkx), dimension(:,:,:), pointer, contiguous :: tketen => null( )
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: chiten => null( )
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: qxten => null( )

  public :: allocate_moloch, init_moloch, moloch

#ifdef SINGLE_PRECISION_REAL
  real(rk4), parameter :: minden = 1.0e-15_rkx
#else
  real(rk8), parameter :: minden = 1.0e-30_rkx
#endif

  logical, parameter :: do_bdy          = .true.
  logical, parameter :: do_vadvtwice    = .true.
  logical, parameter :: do_phys         = .true.
  logical, parameter :: do_convection   = .true.
  logical, parameter :: do_microphysics = .true.
  logical, parameter :: do_radiation    = .true.
  logical, parameter :: do_surface      = .true.
  logical, parameter :: do_pbl          = .true.
  logical, parameter :: do_fulleq       = .true.

  logical :: moloch_realcase = (.not. moloch_do_test_1) .and. &
                               (.not. moloch_do_test_2)
  logical :: lrotllr

  logical :: do_apply_bdy = .true.
  logical :: do_divdamp   = .true.
  logical :: do_divfilter = .true.

  ! Base damping coefficients
  real(rkx), parameter :: dcoff = 0.125_rkx
  real(rkx), parameter :: ddamp = 0.1_rkx
  real(rkx) :: numax

  real(rkx) :: rdzita
  integer(ik4) :: jmin, jmax, imin, imax

  real(rkx) :: dtsound, dtstepa

  contains

#include <pfwsat.inc>

  subroutine allocate_moloch
    implicit none
    integer(ik4) :: k
    call getmem(gzitak,1,kzp1,'moloch:gzitak')
    call getmem(gzitakh,1,kz,'moloch:gzitakh')
    call getmem(laplacian,jdi1,jdi2,idi1,idi2,1,kz,'moloch:laplacian')
    call getmem(wwkw,jce1,jce2,ice1,ice2,2,kzp1,'moloch:wwkw')
    call getmem(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:deltaw')
    call getmem(s,jce1,jce2,ice1,ice2,1,kzp1,'moloch:s')
    call getmem(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
    call getmem(wz,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem(p0,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'moloch:p0')
    call getmem(wfw,jce1,jce2,ice1,ice2,1,kzp1,'moloch:wfw')
    call getmem(wx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:wx')
    call getmem(zpby,jce1,jce2,ici1,ice2ga,1,kz,'moloch:zpby')
    call getmem(zpbw,jci1,jce2ga,ice1,ice2,1,kz,'moloch:zpbw')
    call getmem(mx2,jde1ga,jde2ga,ide1ga,ide2ga,'moloch:mx2')
    call getmem(rmx,jde1ga,jde2ga,ide1ga,ide2ga,'moloch:rmx')
    call getmem(rmu,jde1ga,jde2ga,ide1ga,ide2ga,'moloch:rmu')
    call getmem(rmv,jde1ga,jde2ga,ide1ga,ide2ga,'moloch:rmv')
    call getmem(coru,jde1,jde2,ice1,ice2,'moloch:coru')
    call getmem(corv,jce1,jce2,ide1,ide2,'moloch:corv')
    if ( ibltyp == 2 ) then
      call getmem(tkex,jce1,jce2,ice1,ice2,1,kz,'moloch:tkex')
    end if
    if ( idiag > 0 ) then
      call getmem(ten0,jci1,jci2,ici1,ici2,1,kz,'moloch:ten0')
      call getmem(qen0,jci1,jci2,ici1,ici2,1,kz,'moloch:qen0')
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        call getmem(chiten0,jci1,jci2,ici1,ici2,1,kz,1,ntr,'moloch:chiten0')
      end if
    end if
    call getmem(ud,jde1,jde2,ice1,ice2,1,kz,'moloch:ud')
    call getmem(vd,jce1,jce2,ide1,ide2,1,kz,'moloch:vd')
    if ( do_fulleq ) then
      call getmem(qwltot,jce1,jce2,ice1,ice2,1,kz,'moloch:qwltot')
      call getmem(qwitot,jce1,jce2,ice1,ice2,1,kz,'moloch:qwitot')
    end if
    call getmem(xknu,1,kz,'moloch:xknu')
    numax = 0.25_rkx*dx*dx*rdt
    do concurrent ( k = 1:kz )
      xknu(k) = (ddamp + (1.0_rkx-ddamp)/(k+2.0_rkx))
    end do
  end subroutine allocate_moloch

  subroutine init_moloch
    implicit none
    integer(ik4) :: i, j
    call assignpnt(mddom%msfu,mu)
    call assignpnt(mddom%msfv,mv)
    call assignpnt(mddom%msfx,mx)
    call assignpnt(mddom%rlat,rlat)
    call assignpnt(mddom%hx,hx)
    call assignpnt(mddom%hy,hy)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(mddom%xlon,xlon)
    call assignpnt(mddom%ht,ht)
    call assignpnt(sfs%psa,ps)
    call assignpnt(sfs%tg,ts)
    call assignpnt(sfs%t2m,t2m)
    call assignpnt(mo_atm%fmz,fmz)
    call assignpnt(mo_atm%rfmzu,rfmzu)
    call assignpnt(mo_atm%rfmzv,rfmzv)
    call assignpnt(mo_atm%fmzf,fmzf)
    call assignpnt(mo_atm%pai,pai)
    call assignpnt(mo_atm%tetav,tetav)
    call assignpnt(mo_atm%u,u)
    call assignpnt(mo_atm%ux,ux)
    call assignpnt(mo_atm%v,v)
    call assignpnt(mo_atm%vx,vx)
    call assignpnt(mo_atm%w,w)
    call assignpnt(mo_atm%tvirt,tvirt)
    call assignpnt(mo_atm%zeta,z)
    call assignpnt(mo_atm%p,p)
    call assignpnt(mo_atm%t,t)
    call assignpnt(mo_atm%rho,rho)
    call assignpnt(mo_atm%qx,qx)
    call assignpnt(mo_atm%qs,qsat)
    call assignpnt(mo_atm%qx,qv,iqv)
    call assignpnt(mo_atm%uten,uten)
    call assignpnt(mo_atm%vten,vten)
    call assignpnt(mo_atm%tten,tten)
    call assignpnt(mo_atm%qxten,qxten)
    call assignpnt(mo_atm%qxten,qvten,iqv)
    if ( ipptls > 0 ) then
      call assignpnt(mo_atm%qx,qc,iqc)
      call assignpnt(mo_atm%qxten,qcten,iqc)
      if ( ipptls > 1 ) then
        call assignpnt(mo_atm%qx,qi,iqi)
        call assignpnt(mo_atm%qxten,qiten,iqi)
        call assignpnt(mo_atm%qx,qr,iqr)
        call assignpnt(mo_atm%qx,qs,iqs)
      end if
    end if
    if ( ibltyp == 2 ) then
      call assignpnt(mo_atm%tke,tke)
      call assignpnt(mo_atm%tketen,tketen)
    end if
    if ( ichem == 1 ) then
      call assignpnt(mo_atm%trac,trac)
      call assignpnt(mo_atm%chiten,chiten)
    end if
#ifdef RCEMIP
    coru = 0.0_rkx
    corv = 0.0_rkx
#else
    coru = eomeg2*sin(mddom%ulat(jde1:jde2,ice1:ice2)*degrad)
    corv = eomeg2*sin(mddom%vlat(jce1:jce2,ide1:ide2)*degrad)
#endif
    do concurrent ( j = jde1:jde2, i = ide1:ide2 )
      mx2(j,i) = mx(j,i) * mx(j,i)
      rmx(j,i) = d_one/mx(j,i)
      rmu(j,i) = d_one/mu(j,i)
      rmv(j,i) = d_one/mv(j,i)
    end do
    call exchange_lrbt(mx2,1,jde1,jde2,ide1,ide2)
    call exchange_lrbt(rmx,1,jde1,jde2,ide1,ide2)
    call exchange_lrbt(rmu,1,jde1,jde2,ide1,ide2)
    call exchange_lrbt(rmv,1,jde1,jde2,ide1,ide2)
    gzitak = gzita(zita,mo_ztop,mo_a0)
    gzitakh = gzita(zitah,mo_ztop,mo_a0)
    rdzita = 1.0_rkx / mo_dzita
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      w(j,i,1) = d_zero
    end do
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
    do_divdamp = mo_divdamp
    do_divfilter = mo_divfilter
    do_apply_bdy = ( do_bdy .and. moloch_realcase .and. irceideal == 0 )
    dtstepa = dtsec / real(mo_nadv,rkx)
    dtsound = dtstepa / real(mo_nsound,rkx)
  end subroutine init_moloch
  !
  ! Moloch integration engine
  !
  subroutine moloch
    !@acc use nvtx
    implicit none
    real(rkx) :: maxps, minps, pmax, pmin
    integer(ik4) :: i, j, k, n
    integer(ik4) :: iconvec
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !@acc call nvtxStartRange("moloch")

    iconvec = 0

    call reset_tendencies

    !
    ! Dynamical core - update status variables to new timestep
    !
    call dynamical_core(dtstepa,dtsound)

    !@acc call nvtxStartRange("mkslice")
    call mkslice
    !@acc call nvtxEndRange
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
    ! Update status adding extra terms
    !
    call status_update(dtsec)
    !
    ! ===========================
    ! Diagnostic and end timestep
    ! ===========================
    !
    ! Mass check
    !
    if ( debug_level > 0 ) call massck
    !
#ifdef DEBUG
    do k = 1, kz
      do i = ice1, ice2
        do j = jce1, jce2
          if ( (t(j,i,k) > 350.0_rkx) .or. t(j,i,k) < 170.0_rkx ) then
            write(100+myid,*) 'On : ', myid
            write(100+myid,*) 'At : ', i,j,k
            write(100+myid,*) 'k pai u v w qv qc t tetav'
            do n = 1, kz
              write(100+myid,*) n, pai(j,i,n), u(j,i,n), v(j,i,n), &
                        w(j,i,n), qv(j,i,n), qc(j,i,n), &
                        t(j,i,n), tetav(j,i,n)
            end do
            flush(100+myid)
            call fatal(__FILE__,__LINE__, 'error')
          end if
        end do
      end do
    end do
#endif
    if ( syncro_rep%act( ) .and. rcmtimer%integrating( ) ) then
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
            ' $$$ max, min of ps (mb) = ', pmax*d_r100, pmin*d_r100
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
    !@acc call nvtxEndRange
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine moloch

  subroutine boundary
    !@acc use nvtx
    implicit none
    integer(ik4) :: i, j, k, n
    ! Newtonian factor
    real(rkx), save :: tspectral = 0.0_rkx
    !@acc call nvtxStartRange("boundary")

    if ( idiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ten0(j,i,k) = t(j,i,k)
        qen0(j,i,k) = qx(j,i,k,iqv)
      end do
    end if
    if ( ichem == 1 .and. ichdiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
        chiten0(j,i,k,n) = trac(j,i,k,n)
      end do
    end if

    ! Update external boundary point

    call bdyval

    ! Davies boundary condition on internal point

    call morelax(jdi1,jdi2,ici1,ici2,ba_ud,u,dub)
    call morelax(jci1,jci2,idi1,idi2,ba_vd,v,dvb)
    call morelax(jci1,jci2,ici1,ici2,ba_cr,t,xtb)
    call morelax(jci1,jci2,ici1,ici2,ba_cr,pai,xpaib)
    call morelax(jci1,jci2,ici1,ici2,ba_cr,qv,xqb)
    if ( is_present_qc( ) ) then
      call morelax(jci1,jci2,ici1,ici2,ba_cr,qc,xlb)
    end if
    if ( is_present_qi( ) ) then
      call morelax(jci1,jci2,ici1,ici2,ba_cr,qi,xib)
    end if
    if ( ichem == 1 ) then
      call morelax_chiten(trac)
    end if

    if ( mo_spectral_nudging ) then
      tspectral = tspectral + dtsec
      if ( int(mod(tspectral,dtrad)) == 0 ) then
        call mospectral_nudge(jce1,jce2,ice1,ice2,jci1,jci1,ici1,ici2,t,xtb)
        call mospectral_nudge(jde1,jde2,ice1,ice2,jdi1,jdi2,ici1,ici2,ud,dub)
        call mospectral_nudge(jce1,jce2,ide1,ide2,jci1,jci2,idi1,idi2,vd,dvb)
      end if
    end if

    if ( idiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        tdiag%bdy(j,i,k) = t(j,i,k) - ten0(j,i,k)
        qdiag%bdy(j,i,k) = qv(j,i,k) - qen0(j,i,k)
      end do
    end if
    if ( ichem == 1 .and. ichdiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
        cbdydiag(j,i,k,n) = trac(j,i,k,n) - chiten0(j,i,k,n)
      end do
    end if
    call uvstagtox(u,v,ux,vx)
    !@acc call nvtxEndRange
  end subroutine boundary

  subroutine divergence_diffusion(dts)
    implicit none
    integer(ik4) :: j, i, k
    real(rkx), intent(in) :: dts

    call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      laplacian(j,i,k) = &
       rdx**2 * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) - 2.0_rkx * zdiv2(j,i,k)) + &
       rdx**2 * (zdiv2(j,i-1,k) + zdiv2(j,i+1,k) - 2.0_rkx * zdiv2(j,i,k))
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      if ( j > nspgx+2 .or. j < jcross2+nspgx-1 ) cycle
      if ( i > nspgx+2 .or. i < icross2+nspgx-1 ) cycle
      zdiv2(j,i,k) = zdiv2(j,i,k) + dts * numax * xknu(k) * laplacian(j,i,k)
    end do
  end subroutine divergence_diffusion

  subroutine sound(dts)
    !@acc use nvtx
    implicit none
    real(rkx), intent(in) :: dts
    integer(ik4) :: i, j, k, nsound
    real(rkx) :: dtrdx, dtrdy, dtrdz, zcs2
    real(rkx) :: zum, zup, zvm, zvp, zuh, zvh
    real(rkx) :: zrom1w, zwexpl, zu, zd, zrapp
    real(rkx) :: zqs, zdth
    real(rkx) :: zcx, zcy, zfz
    real(rkx) :: zrom1u, zcor1u, zrom1v, zcor1v

    !@acc call nvtxStartRange("sound")
    dtrdx = dts*rdx
    dtrdy = dts*rdx
    dtrdz = dts*rdzita
    zcs2 = dtrdz**2*rdrcv

    if ( .not. do_fulleq ) then
      call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
    end if

    !  sound waves

    do nsound = 1, mo_nsound

      call exchange_lr(u,1,jde1,jde2,ice1,ice2,1,kz)
      call exchange_bt(v,1,jce1,jce2,ide1,ide2,1,kz)

      do concurrent ( j = jde1:jde2, i = ice1:ice2, k = 1:kz )
        ud(j,i,k) = u(j,i,k)
      end do
      do concurrent ( j = jce1:jce2, i = ide1:ide2, k = 1:kz )
        vd(j,i,k) = v(j,i,k)
      end do

      ! partial definition of the generalized vertical velocity

      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
        zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
        s(j,i,kzp1) = -0.5_rkx * (zuh+zvh)
        w(j,i,kzp1) = -s(j,i,kzp1)
      end do

      ! Equation 10, generalized vertical velocity

      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kz )
        zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) +    &
              (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
        zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) +    &
              (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
        s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
      end do

      ! Part of divergence (except w contribution) put in zdiv2
      ! Equation 16

      if ( lrotllr ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          zum = dtrdx * u(j,i,k)   * rfmzu(j,i,k)
          zup = dtrdx * u(j+1,i,k) * rfmzu(j+1,i,k)
          zvm = dtrdy * v(j,i,k)   * rfmzv(j,i,k)   * rmv(j,i)
          zvp = dtrdy * v(j,i+1,k) * rfmzv(j,i+1,k) * rmv(j,i+1)
          zdiv2(j,i,k) = fmz(j,i,k) * mx(j,i) * ((zup-zum) + (zvp-zvm))
        end do
      else
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          zum = dtrdx * u(j,i,k)   * rfmzu(j,i,k)   * rmu(j,i)
          zup = dtrdx * u(j+1,i,k) * rfmzu(j+1,i,k) * rmu(j+1,i)
          zvm = dtrdy * v(j,i,k)   * rfmzv(j,i,k)   * rmv(j,i)
          zvp = dtrdy * v(j,i+1,k) * rfmzv(j,i+1,k) * rmv(j,i+1)
          zdiv2(j,i,k) = fmz(j,i,k) * mx2(j,i) * ((zup-zum) + (zvp-zvm))
        end do
      end if

      if ( do_divdamp ) then
        call divergence_damping(dts)
      end if
      if ( do_divfilter ) then
        call divergence_diffusion(dts)
      end if

      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        zdiv2(j,i,k) = zdiv2(j,i,k) + dtrdz * fmz(j,i,k) * &
                  (s(j,i,k) - s(j,i,k+1))
      end do

      ! new w (implicit scheme) from Equation 19

      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        do k = kz, 2, -1
          deltaw(j,i,k) = -w(j,i,k)
          ! explicit w:
          !    it must be consistent with the initialization of pai
          zrom1w = d_half * cpd * fmzf(j,i,k) * &
                   (tetav(j,i,k-1) + tetav(j,i,k))
          zrom1w = zrom1w - cpd * w(j,i,k) * &
                   fmzf(j,i,k)*fmzf(j,i,k) * &
                   real(nsound,rkx) * dtrdz * &
                   (tetav(j,i,k-1) - tetav(j,i,k)) !! GW
          if ( qv(j,i,k) > 0.96_rkx*qsat(j,i,k) .and. &
               w(j,i,k) > 0.1_rkx ) then
            zqs = d_half*(qsat(j,i,k)+qsat(j,i,k-1))
            zdth = egrav*w(j,i,k)*real(nsound-1,rkx)*dts*wlhv*wlhv* &
                zqs/(cpd*pai(j,i,k-1)*rwat*t(j,i,k-1)*t(j,i,k-1))
            zrom1w = zrom1w + zdth*fmzf(j,i,k)
          end if
          zwexpl = w(j,i,k) - zrom1w * dtrdz * &
                   (pai(j,i,k-1) - pai(j,i,k)) - egrav*dts
          zwexpl = zwexpl + rdrcv * zrom1w * dtrdz * &
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

      ! 2nd loop for the tridiagonal inversion
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        do k = 2, kz
          w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
          deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
        end do
      end do

      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        zdiv2(j,i,k) = zdiv2(j,i,k) + dtrdz * fmz(j,i,k) * &
                  (w(j,i,k) - w(j,i,k+1))
      end do

      if ( do_fulleq ) then
        if ( ipptls > 0 ) then
          do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
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
          call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
        end if
      end if

      ! new Exner function (Equation 19)

      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv2(j,i,k))
      end do

      call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1)

      ! horizontal momentum equations

      if ( lrotllr ) then
        ! Equation 17
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          zcx = dtrdx * mu(j,i)
          zfz = egrav * dts + 0.25_rkx * &
              (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
               deltaw(j,i,k)   + deltaw(j,i,k+1))
          zrom1u = 0.5_rkx * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
          zcor1u = coru(j,i) * dts * vd(j,i,k)
          ! Equation 17
          u(j,i,k) = u(j,i,k) + zcor1u - &
                     zfz * hx(j,i) * gzitakh(k) - &
                     zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
        end do
        ! Equation 18
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          zcy = dtrdy
          zfz = egrav * dts + 0.25_rkx * &
              (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
               deltaw(j,i,k)   + deltaw(j,i,k+1))
          zrom1v = 0.5_rkx * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
          zcor1v = corv(j,i) * dts * ud(j,i,k)
          ! Equation 18
          v(j,i,k) = v(j,i,k) - zcor1v - &
                     zfz * hy(j,i) * gzitakh(k) -  &
                     zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
        end do
      else
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          zcx = dtrdx * mu(j,i)
          zfz = egrav * dts + 0.25_rkx * &
              (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
               deltaw(j,i,k)   + deltaw(j,i,k+1))
          zrom1u = 0.5_rkx * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
          zcor1u = coru(j,i) * dts * vd(j,i,k)
          ! Equation 17
          u(j,i,k) = u(j,i,k) + zcor1u - &
                     zfz * hx(j,i) * gzitakh(k) - &
                     zcx * zrom1u * (pai(j,i,k) - pai(j-1,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          zcy = dtrdy * mv(j,i)
          zfz = egrav * dts + 0.25_rkx * &
              (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
               deltaw(j,i,k)   + deltaw(j,i,k+1))
          zrom1v = 0.5_rkx * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
          zcor1v = corv(j,i) * dts * ud(j,i,k)
          ! Equation 18
          v(j,i,k) = v(j,i,k) - zcor1v - &
                     zfz * hy(j,i) * gzitakh(k) - &
                     zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
        end do
      end if

    end do ! sound loop

    ! complete computation of generalized vertical velocity
    ! Complete Equation 10
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kz )
      s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      s(j,i,1) = 0.0_rkx
      s(j,i,kzp1) = 0.0_rkx
    end do
    !@acc call nvtxEndRange
  end subroutine sound

  subroutine divergence_damping(dts)
    implicit none
    real(rkx), intent(in) :: dts
    integer(ik4) :: i, j, k
    real(rkx) :: xdam, dxrdt

    dxrdt = dcoff * dx/dts
    call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
    if ( lrotllr ) then
      do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
        xdam = dxrdt * xknu(k) * mu(j,i)*mu(j,i)*rmx(j,i)
        u(j,i,k) = u(j,i,k) + xdam * (zdiv2(j,i,k)-zdiv2(j-1,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
        xdam = dxrdt * xknu(k)
        v(j,i,k) = v(j,i,k) + xdam * (zdiv2(j,i,k)-zdiv2(j,i-1,k))
      end do
    else
      do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
        xdam = dxrdt * xknu(k) * mu(j,i)*mu(j,i)*rmx(j,i)
        u(j,i,k) = u(j,i,k) + xdam * (zdiv2(j,i,k)-zdiv2(j-1,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
        xdam = dxrdt * xknu(k) * mv(j,i)*mv(j,i)*rmx(j,i)
        v(j,i,k) = v(j,i,k) + xdam * (zdiv2(j,i,k)-zdiv2(j,i-1,k))
      end do
    end if
  end subroutine divergence_damping

  subroutine advection(dta)
    !@acc use nvtx
    implicit none
    real(rkx), intent(in) :: dta
    integer(ik4) :: i, j, k, n
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ptr => null( )
    real(rkx) :: tanx, tany, dlat

    !@acc call nvtxStartRange("advection")

    ! Compute U,V on cross points
    call uvstagtox(u,v,ux,vx)
    call zstagtoh(w,wx)

    ! Compute TKE if required on zita levels
    if ( ibltyp == 2 ) then
      call zstagtoh(tke,tkex)
    end if

    call wafone(tetav,dta)
    call wafone(pai,dta)
    call wafone(ux,dta)
    call wafone(vx,dta)
    call wafone(wx,dta)

    call wafone(qv,dta)
    if ( ipptls > 0 ) then
      do n = iqfrst, nqx
        call assignpnt(qx,ptr,n)
        call wafone(ptr,dta)
      end do
    end if
    if ( ibltyp == 2 ) then
      call wafone(tkex,dta)
    end if
    if ( ichem == 1 ) then
      do n = 1, ntr
        call assignpnt(trac,ptr,n)
        call wafone(ptr,dta)
      end do
    end if

    !!! Curvature terms !!!

    if ( lrotllr ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        dlat = degrad * 0.5_rkx * (rlat(i)+rlat(i+1))
        tanx = sin(dlat)*mx(j,i)*rearthrad
        ux(j,i,k) = ux(j,i,k) + ux(j,i,k) * vx(j,i,k) * tanx * dta
        vx(j,i,k) = vx(j,i,k) - ux(j,i,k) * ux(j,i,k) * tanx * dta
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        tanx = (mu(j-1,i)-mu(j,i))*rdx
        tany = (mv(j,i-1)-mv(j,i))*rdx
        ux(j,i,k) = ux(j,i,k) + ux(j,i,k) * vx(j,i,k) * tanx * dta
        vx(j,i,k) = vx(j,i,k) - ux(j,i,k) * ux(j,i,k) * tany * dta
      end do
    end if

    ! Interpolate on staggered points
    call xtouvstag(ux,vx,u,v)
    call htozstag(wx,w)

    if ( ibltyp == 2 ) then
      ! Back to half-levels
      call htozstag(tkex,tke)
    end if
    !@acc call nvtxEndRange
  end subroutine advection

  subroutine wafone(pp,dta)
    implicit none
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: pp
    real(rkx), intent(in) :: dta
    integer(ik4) :: j, i, k
    real(rkx) :: dtrdx, dtrdy, dtrdz
    real(rkx), parameter :: wlow  = 0.0_rkx
    real(rkx), parameter :: whigh = 2.0_rkx
    real(rkx) :: zamu, is, r, b, zphi, zzden, zdv
    real(rkx) :: zhxvtn, zhxvts, zcostx
    real(rkx) :: zrfmu, zrfmd
    real(rkx) :: zrfmn, zrfms
    real(rkx) :: zrfme, zrfmw
    integer(ik4) :: k1, k1p1
    integer(ik4) :: ih, ihm1
    integer(ik4) :: jh, jhm1

    dtrdx = dta*rdx
    dtrdy = dta*rdx
    dtrdz = dta*rdzita
    if ( do_vadvtwice ) then
      dtrdz = 0.5_rkx * dtrdz
    end if

    ! Vertical advection
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      wfw(j,i,1) = d_zero
      wfw(j,i,kzp1) = d_zero
    end do

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kzm1 )
      zamu = s(j,i,k+1) * dtrdz
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
      zzden = pp(j,i,k)-pp(j,i,k+1)
      zzden = sign(max(abs(zzden),minden),zzden)
      r = (pp(j,i,k1)-pp(j,i,k1p1))/zzden
      b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
      zphi = is + zamu * b - is * b
      wfw(j,i,k+1) = 0.5_rkx * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                             (d_one-zphi)*pp(j,i,k))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
      zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
      zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
      wz(j,i,k) = pp(j,i,k) - &
        wfw(j,i,k)*zrfmu + wfw(j,i,k+1)*zrfmd + zdv
    end do

    if ( do_vadvtwice ) then

      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kzm1 )
        zamu = s(j,i,k+1) * dtrdz
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
        zzden = wz(j,i,k)-wz(j,i,k+1)
        zzden = sign(max(abs(zzden),minden),zzden)
        r = (wz(j,i,k1)-wz(j,i,k1p1))/zzden
        b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
        zphi = is + zamu * b - is * b
        wfw(j,i,k+1) = 0.5_rkx * s(j,i,k+1) * &
          ((d_one+zphi)*wz(j,i,k+1) + (d_one-zphi)*wz(j,i,k))
      end do
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
        zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
        zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * wz(j,i,k)
        wz(j,i,k) = wz(j,i,k) - wfw(j,i,k)*zrfmu + &
                                wfw(j,i,k+1)*zrfmd + zdv
      end do

    end if

    call exchange_bt(wz,2,jce1,jce2,ice1,ice2,1,kz)

    if ( lrotllr ) then

      ! Meridional advection
      do concurrent ( j = jce1:jce2, i = ici1:ice2ga, k = 1:kz )
        zamu = v(j,i,k) * dtrdy
        if ( zamu > d_zero ) then
          is = d_one
          ih = i-1
        else
          is = -d_one
          ih = min(i+1,imax)
        end if
        ihm1 = max(ih-1,imin)
        zzden = wz(j,i,k)-wz(j,i-1,k)
        zzden = sign(max(abs(zzden),minden),zzden)
        r = (wz(j,ih,k)-wz(j,ihm1,k))/zzden
        b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
        zphi = is + zamu*b - is*b
        zpby(j,i,k) = 0.5_rkx * v(j,i,k) * &
            ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
      end do
      do concurrent ( j = jce1:jce2, i = ici1:ici2, k = 1:kz )
        zhxvtn = dtrdy * rmv(j,i+1) * mx(j,i)
        zhxvts = dtrdy * rmv(j,i)   * mx(j,i)
        zrfmn = zhxvtn * fmz(j,i,k) * rfmzv(j,i+1,k)
        zrfms = zhxvts * fmz(j,i,k) * rfmzv(j,i,k)
        zdv = (v(j,i+1,k) * zrfmn - v(j,i,k) * zrfms) * pp(j,i,k)
        p0(j,i,k) = wz(j,i,k) + &
              (zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv)
      end do

      call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

      ! Zonal advection

      do concurrent ( j = jci1:jce2ga, i = ici1:ici2, k = 1:kz )
        zamu = u(j,i,k) * mu(j,i) * dtrdx
        if ( zamu > d_zero ) then
          is = d_one
          jh = j-1
        else
          is = -d_one
          jh = min(j+1,jmax)
        end if
        jhm1 = max(jh-1,jmin)
        zzden = p0(j,i,k)-p0(j-1,i,k)
        zzden = sign(max(abs(zzden),minden),zzden)
        r = (p0(jh,i,k)-p0(jhm1,i,k))/zzden
        b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
        zphi = is + zamu*b - is*b
        zpbw(j,i,k) = 0.5_rkx * u(j,i,k) * &
               ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        zcostx = dtrdx * mx(j,i)
        zrfmw = zcostx * fmz(j,i,k) * rfmzu(j,i,k)
        zrfme = zcostx * fmz(j,i,k) * rfmzu(j+1,i,k)
        zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
        pp(j,i,k) = p0(j,i,k) + &
               zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv
      end do

    else ! Not the ROTLLR projection

      ! Meridional advection
      do concurrent ( j = jce1:jce2, i = ici1:ice2ga, k = 1:kz )
        zamu = v(j,i,k) * mv(j,i) * dtrdy
        if ( zamu > d_zero ) then
          is = d_one
          ih = i-1
        else
          is = -d_one
          ih = min(i+1,imax)
        end if
        ihm1 = max(ih-1,imin)
        zzden = wz(j,i,k)-wz(j,i-1,k)
        zzden = sign(max(abs(zzden),minden),zzden)
        r = (wz(j,ih,k)-wz(j,ihm1,k))/zzden
        b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
        zphi = is + zamu*b - is*b
        zpby(j,i,k) = 0.5_rkx * v(j,i,k) * &
            ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
      end do
      do concurrent ( j = jce1:jce2, i = ici1:ici2, k = 1:kz )
        zrfmn = dtrdy * fmz(j,i,k) * rfmzu(j,i+1,k)
        zrfms = dtrdy * fmz(j,i,k) * rfmzu(j,i,k)
        zdv = (v(j,i+1,k) * rmv(j,i+1) * zrfmn - &
               v(j,i,k)   * rmv(j,i)   * zrfms) * pp(j,i,k)
        p0(j,i,k) = wz(j,i,k) + &
          mx2(j,i) * (zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv)
      end do

      call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

      ! Zonal advection
      do concurrent ( j = jci1:jce2ga, i = ici1:ici2, k = 1:kz )
        zamu = u(j,i,k) * mu(j,i) * dtrdx
        if ( zamu > d_zero ) then
          is = d_one
           jh = j-1
        else
          is = -d_one
          jh = min(j+1,jmax)
        end if
        jhm1 = max(jh-1,jmin)
        zzden = p0(j,i,k)-p0(j-1,i,k)
        zzden = sign(max(abs(zzden),minden),zzden)
        r = (p0(jh,i,k)-p0(jhm1,i,k))/zzden
        b = max(wlow, min(whigh, max(r, min(d_two*r,d_one))))
        zphi = is + zamu*b - is*b
        zpbw(j,i,k) = 0.5_rkx * u(j,i,k) * &
               ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        zrfmw = dtrdx * fmz(j,i,k) * rfmzu(j,i,k)
        zrfme = dtrdx * fmz(j,i,k) * rfmzu(j+1,i,k)
        zdv = (u(j+1,i,k) * rmu(j+1,i) * zrfme - &
               u(j,i,k)   * rmu(j,i)   * zrfmw) * pp(j,i,k)
        pp(j,i,k) = p0(j,i,k) + &
            mx2(j,i) * (zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv)
      end do

    end if

  end subroutine wafone

  subroutine reset_tendencies
    !@acc use nvtx
    implicit none
    integer(ik4) :: i, j, k, n
    !@acc call nvtxStartRange("reset_tendencies")
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kzp1 )
      s(j,i,k) = d_zero
    end do
    do concurrent ( j = jci1ga:jci2ga, i = ici1ga:ici2ga, k = 1:kz )
      zdiv2(j,i,k) = d_zero
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kzp1 )
      wwkw(j,i,k) = d_zero
    end do
    do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kzp1 )
      deltaw(j,i,k) = d_zero
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      tten(j,i,k) = d_zero
      uten(j,i,k) = d_zero
      vten(j,i,k) = d_zero
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:nqx )
      qxten(j,i,k,n) = d_zero
    end do
    if ( ichem == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
        chiten(j,i,k,n) = d_zero
      end do
    end if
    if ( ibltyp == 2 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
        tketen(j,i,k) = d_zero
      end do
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      cldfra(j,i,k) = d_zero
      cldlwc(j,i,k) = d_zero
    end do

    if ( idiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ten0(j,i,k) = t(j,i,k)
        qen0(j,i,k) = qv(j,i,k)
      end do
      if ( ichem == 1 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
          chiten0(j,i,k,n) = trac(j,i,k,n)
        end do
      end if
    end if
    !@acc call nvtxEndRange
  end subroutine reset_tendencies

  subroutine dynamical_core(dta,dts)
    !@acc use nvtx
    implicit none
    real(rkx), intent(in) :: dta, dts
    real(rkx) :: fice
    integer(ik4) :: i, j, k, n, nadv
    !@acc call nvtxStartRange("dynamical_core")
    !
    ! Prepare fields
    !
    if ( do_fulleq ) then
      if ( ipptls > 0 ) then
        if ( ipptls > 1 ) then
          do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
            qwltot(j,i,k) = qc(j,i,k) + qr(j,i,k)
            qwitot(j,i,k) = qi(j,i,k) + qs(j,i,k)
          end do
        else
          do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
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
        end if
      else
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          qwltot(j,i,k) = d_zero
          qwitot(j,i,k) = d_zero
        end do
      end if
    end if

    if ( idiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ten0(j,i,k) = t(j,i,k)
        qen0(j,i,k) = qv(j,i,k)
      end do
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
          chiten0(j,i,k,n) = trac(j,i,k,n)
        end do
      end if
    end if
    !
    ! ############################################
    !
    do nadv = 1, mo_nadv

      call sound(dts)

      call advection(dta)

    end do ! Advection loop

    !
    ! Compute lateral boundary condition relaxation
    !
    if ( do_apply_bdy ) then
      call boundary
    end if

    !
    ! ############################################
    !
    ! Update "physical variables" to register change
    !
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
      p(j,i,k) = (pai(j,i,k)**cpovr) * p00
    end do

    call tvirt_to_temp( )

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
    end do

    call extrapolate_surface_pressure( )

    if ( idiag > 0 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        tdiag%adh(j,i,k) = (t(j,i,k) - ten0(j,i,k)) * rdt
        qdiag%adh(j,i,k) = (qv(j,i,k) - qen0(j,i,k)) * rdt
      end do
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
          cadvhdiag(j,i,k,n) = (trac(j,i,k,n) - chiten0(j,i,k,n)) * rdt
        end do
      end if
    end if
    !@acc call nvtxEndRange
  end subroutine dynamical_core

  subroutine physical_parametrizations
    !@acc use nvtx
    implicit none
    integer(ik4) :: i, j, k, n
    logical :: loutrad, labsem
    !@acc call nvtxStartRange("physical_parametrizations")
#ifdef DEBUG
    do k = 1, kz
      do i = ice1, ice2
        do j = jce1, jce2
          if ( (t(j,i,k) > 350.0_rkx) .or. t(j,i,k) < 170.0_rkx ) then
            write(100+myid,*) 'Before Phys On : ', myid
            write(100+myid,*) 'At : ', i,j,k
            write(100+myid,*) 'k pai u v w qv qc t tetav'
            do n = 1, kz
              write(100+myid,*) n, pai(j,i,n), u(j,i,n), v(j,i,n), &
                        w(j,i,n), qv(j,i,n), qc(j,i,n), &
                        t(j,i,n), tetav(j,i,n)
            end do
            flush(100+myid)
            call fatal(__FILE__,__LINE__, 'error')
          end if
        end do
      end do
    end do
#endif
    if ( do_convection ) then
      if ( any(icup > 0) ) then
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            ten0(j,i,k) = tten(j,i,k)
            qen0(j,i,k) = qxten(j,i,k,iqv)
          end do
        end if
        if ( ichem == 1 .and. ichdiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                          k = 1:kz, n = 1:ntr )
            chiten0(j,i,k,n) = chiten(j,i,k,n)
          end do
        end if
        !@acc call nvtxStartRange("cumulus")
        call cumulus
        !@acc call nvtxEndRange
        if ( ichem == 1 ) then
          if ( ichcumtra == 1 ) then
            if ( debug_level > 3 .and. myid == italk ) then
              write(stdout,*) 'Calling cumulus transport at ', &
                         trim(rcmtimer%str())
            end if
            !@acc call nvtxStartRange("cumtran")
            call cumtran(trac)
            !@acc call nvtxEndRange
          end if
        end if
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%con(j,i,k) = tten(j,i,k) - ten0(j,i,k)
            qdiag%con(j,i,k) = qxten(j,i,k,iqv) - qen0(j,i,k)
          end do
        end if
        if ( ichem == 1 .and. ichdiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                          k = 1:kz, n = 1:ntr )
            cconvdiag(j,i,k,n) = chiten(j,i,k,n) - chiten0(j,i,k,n)
          end do
        end if
      else
        if ( any(icup < 0) ) then
          !@acc call nvtxStartRange("shallow_convection")
          call shallow_convection
          !@acc call nvtxEndRange
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              tdiag%con(j,i,k) = tten(j,i,k) - ten0(j,i,k)
              qdiag%con(j,i,k) = qxten(j,i,k,iqv) - qen0(j,i,k)
            end do
          end if
        end if
      end if
    end if
    !
    !------------------------------------------------
    ! Large scale precipitation microphysical schemes
    !------------------------------------------------
    !
    if ( do_microphysics ) then
      if ( ipptls > 0 ) then
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            ten0(j,i,k) = tten(j,i,k)
            qen0(j,i,k) = qxten(j,i,k,iqv)
          end do
        end if
        ! Cumulus clouds
        if ( icldfrac /= 2 ) then
          !@acc call nvtxStartRange("cucloud")
          call cucloud
          !@acc call nvtxEndRange
        end if
        ! Save cumulus cloud fraction for chemistry before it is
        ! overwritten in cldfrac
        if ( ichem == 1 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            convcldfra(j,i,k) = cldfra(j,i,k)
          end do
        end if
        ! Clouds and large scale precipitation
        !@acc call nvtxStartRange("cldfrac")
        call cldfrac(cldlwc,cldfra)
        !@acc call nvtxEndRange
        !@acc call nvtxStartRange("microscheme")
        call microscheme
        !@acc call nvtxEndRange
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%lsc(j,i,k) = tten(j,i,k) - ten0(j,i,k)
            qdiag%lsc(j,i,k) = qxten(j,i,k,iqv) - qen0(j,i,k)
          end do
        end if
      end if
    end if
    !
    !------------------------------------------------
    !       Call radiative transfer package
    !------------------------------------------------
    !
    if ( do_radiation ) then
      if ( rcmtimer%start() .or. syncro_rad%will_act( ) ) then
        if ( debug_level > 3 .and. myid == italk ) then
          write(stdout,*) &
            'Calling radiative transfer at ',trim(rcmtimer%str())
        end if
        ! calculate albedo
        !@acc call nvtxStartRange("surface_albedo")
        call surface_albedo
        !@acc call nvtxEndRange
        if ( iclimao3 == 1 ) then
          !@acc call nvtxStartRange("updateo3")
          call updateo3(rcmtimer%idate,scenario)
          !@acc call nvtxEndRange
        end if
        if ( iclimaaer == 1 ) then
          !@acc call nvtxStartRange("updateaero_1")
          call updateaerosol(rcmtimer%idate)
          !@acc call nvtxEndRange
        else if ( iclimaaer == 2 ) then
          !@acc call nvtxStartRange("updateaero_2")
          call updateaeropp(rcmtimer%idate)
          !@acc call nvtxEndRange
        else if ( iclimaaer == 3 ) then
          !@acc call nvtxStartRange("updateaero_3")
          call updateaeropp_cmip6(rcmtimer%idate)
          !@acc call nvtxEndRange
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
        !@acc call nvtxStartRange("radiation")
        call radiation(rcmtimer%year,rcmtimer%month, &
                       rcmtimer%day,loutrad,labsem)
        !@acc call nvtxEndRange
      end if
      !
      ! Add radiative transfer package-calculated heating rates to
      ! temperature tendency (deg/sec)
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        tten(j,i,k) = tten(j,i,k) + heatrt(j,i,k)
      end do
      if ( idiag > 0 ) tdiag%rad = heatrt
    end if
    !
    !------------------------------------------------
    !            Call Surface model
    !------------------------------------------------
    !
    if ( do_surface ) then
      if ( rcmtimer%start() .or. syncro_srf%will_act( ) ) then
        if ( debug_level > 3 .and. myid == italk ) then
          write(stdout,*) 'Calling surface model at ',trim(rcmtimer%str())
        end if
        !@acc call nvtxStartRange("surface_model")
        call surface_model
        !@acc call nvtxEndRange
      end if
    end if
    !
    !------------------------------------------------
    !             Call PBL scheme
    !------------------------------------------------
    !
    if ( do_pbl ) then
      if ( ibltyp > 0 ) then
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            ten0(j,i,k) = tten(j,i,k)
            qen0(j,i,k) = qxten(j,i,k,iqv)
          end do
        end if
        if ( ichem == 1 .and. ichdiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                          k = 1:kz, n = 1:ntr )
            chiten0(j,i,k,n) = chiten(j,i,k,n)
          end do
        end if
        !@acc call nvtxStartRange("pblscheme")
        call pblscheme
        !@acc call nvtxEndRange
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%tbl(j,i,k) = tten(j,i,k) - ten0(j,i,k)
            qdiag%tbl(j,i,k) = qxten(j,i,k,iqv) - qen0(j,i,k)
          end do
        end if
        if ( ichem == 1 .and. ichdiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                          k = 1:kz, n = 1:ntr )
            ctbldiag(j,i,k,n) = chiten(j,i,k,n) - chiten0(j,i,k,n)
          end do
        end if
      end if
    end if

    if ( do_microphysics ) then
      if ( ipptls == 1 ) then
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            ten0(j,i,k) = tten(j,i,k)
            qen0(j,i,k) = qxten(j,i,k,iqv)
          end do
        end if
        call condtq
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%lsc(j,i,k) = tdiag%lsc(j,i,k)+tten(j,i,k)-ten0(j,i,k)
            qdiag%lsc(j,i,k) = qdiag%lsc(j,i,k)+qxten(j,i,k,iqv)-qen0(j,i,k)
          end do
        end if
      end if
    end if
    !-------------------------------------------------------------
    !call chemistry/aerosol schemes
    !----------------------------------------------------------------------
    if ( ichem == 1 ) then
      call tractend2(rcmtimer%month,rcmtimer%day,declin)
    end if
    !@acc call nvtxEndRange
  end subroutine physical_parametrizations

  subroutine status_update(dtinc)
    !@acc use nvtx
    implicit none
    real(rkx), intent(in) :: dtinc
    integer(ik4) :: i, j, k, n
    !@acc call nvtxStartRange("status_update")

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      t(j,i,k)   = t(j,i,k)   + dtinc * tten(j,i,k)
      ux(j,i,k)  = ux(j,i,k)  + dtinc * uten(j,i,k)
      vx(j,i,k)  = vx(j,i,k)  + dtinc * vten(j,i,k)
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      qx(j,i,k,iqv) = qx(j,i,k,iqv) + dtinc * qxten(j,i,k,iqv)
      if ( qx(j,i,k,iqv) < 1.0E-8_rkx ) qx(j,i,k,iqv) = 1.0E-8_rkx
    end do
    do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                    k = 1:kz, n = iqfrst:nqx)
      qx(j,i,k,n) = qx(j,i,k,n) + dtinc * qxten(j,i,k,n)
      if ( qx(j,i,k,n) < 1.0E-20_rkx ) qx(j,i,k,n) = 0.0_rkx
    end do
    if ( ibltyp == 2 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
        tke(j,i,k) = tke(j,i,k) + dtinc * tketen(j,i,k)
        if ( tke(j,i,k) < tkemin ) tke(j,i,k) = tkemin
      end do
    end if
    if ( ichem == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
        trac(j,i,k,n) = trac(j,i,k,n) + dtinc * chiten(j,i,k,n)
        if ( trac(j,i,k,n) < 0.0_rkx ) trac(j,i,k,n) = 0.0_rkx
      end do
    end if

    call temp_to_tvirt( )

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
      rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
    end do

    call xtouvstag(ux,vx,u,v)

    !@acc call nvtxEndRange
  end subroutine status_update

  subroutine zstagtoh(fl,hl)
    implicit none
    real(rkx), intent(in), dimension(:,:,:), pointer, contiguous :: fl
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: hl
    integer(ik4) :: i, j, k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kzm1 )
      hl(j,i,k) = 0.5625_rkx * (fl(j,i,k+1)+fl(j,i,k)) - &
                  0.0625_rkx * (fl(j,i,k+2)+fl(j,i,k-1))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      hl(j,i,1)  = 0.5_rkx * (fl(j,i,2)+fl(j,i,1))
      hl(j,i,kz) = 0.5_rkx * (fl(j,i,kzp1)+fl(j,i,kz))
    end do
  end subroutine zstagtoh

  subroutine htozstag(hl,fl)
    implicit none
    real(rkx), intent(in), dimension(:,:,:), pointer, contiguous :: hl
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: fl
    integer(ik4) :: i, j, k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 3:kzm1 )
      fl(j,i,k) = 0.5625_rkx * (hl(j,i,k)  +hl(j,i,k-1)) - &
                  0.0625_rkx * (hl(j,i,k+1)+hl(j,i,k-2))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      fl(j,i,2) = 0.5_rkx * (hl(j,i,2)  +hl(j,i,1))
      fl(j,i,kz) = 0.5_rkx * (hl(j,i,kz)+hl(j,i,kzm1))
    end do
  end subroutine htozstag

  subroutine xtouvstag(ux,vx,u,v)
    !@acc use nvtx
    implicit none
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: ux, vx
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: u, v
    integer(ik4) :: i, j, k
    !@acc call nvtxStartRange("xtouvstag")

    call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

    ! Back to wind points: U (fourth order)

    do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
      u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                 0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
    end do
    if ( ma%has_bdyright ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi2,i,k) = 0.5_rkx * (ux(jci2,i,k)+ux(jce2,i,k))
      end do
    end if
    if ( ma%has_bdyleft ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi1,i,k) = 0.5_rkx * (ux(jci1,i,k)+ux(jce1,i,k))
      end do
    end if

    ! Back to wind points: V (fourth order)

    do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
      v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                 0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
    end do
    if ( ma%has_bdytop ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi2,k) = 0.5_rkx * (vx(j,ici2,k)+vx(j,ice2,k))
      end do
    end if
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi1,k) = 0.5_rkx * (vx(j,ici1,k)+vx(j,ice1,k))
      end do
    end if
    !@acc call nvtxEndRange
  end subroutine xtouvstag

  subroutine uvstagtox(u,v,ux,vx)
    !@acc use nvtx
    implicit none
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: u, v
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: ux, vx
    integer(ik4) :: i, j, k
    !@acc call nvtxStartRange("uvstagtox")

    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)

    ! Compute U-wind on T points

    do concurrent ( j = jci1:jci2, i = ice1:ice2, k = 1:kz )
      ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                  0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
    end do
    if ( ma%has_bdyleft ) then
      do concurrent ( i = ice1:ice2, k = 1:kz )
        ux(jce1,i,k) = 0.5_rkx * (u(jde1,i,k)+u(jdi1,i,k))
      end do
    end if
    if ( ma%has_bdyright ) then
      do concurrent ( i = ice1:ice2, k = 1:kz )
        ux(jce2,i,k) = 0.5_rkx * (u(jde2,i,k)+u(jdi2,i,k))
      end do
    end if

    ! Compute V-wind on T points

    do concurrent ( j = jce1:jce2, i = ici1:ici2, k = 1:kz )
      vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                  0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
    end do
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jce1:jce2, k = 1:kz )
        vx(j,ice1,k) = 0.5_rkx * (v(j,ide1,k)+v(j,idi1,k))
      end do
    end if
    if ( ma%has_bdytop ) then
      do concurrent ( j = jce1:jce2, k = 1:kz )
        vx(j,ice2,k) = 0.5_rkx * (v(j,ide2,k)+v(j,idi2,k))
      end do
    end if
    !@acc call nvtxEndRange
  end subroutine uvstagtox

  subroutine extrapolate_surface_pressure( )
    implicit none
    real(rkx) :: zdgz, zh
    integer(ik4) :: i, j
    zh = 0.5_rkx*mo_dzita
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
     zdgz = egrav*md_zeta(zh,ht(j,i),mo_ztop,mo_h,mo_a0)
     ps(j,i) = p(j,i,kz) * exp(zdgz/(rgas*tvirt(j,i,kz)))
    end do
  end subroutine extrapolate_surface_pressure

  subroutine temp_to_tvirt( )
    implicit none
    integer(ik4) :: i, j, k
    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
          tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
             qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
        end do
      else
        do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
          tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - qc(j,i,k))
        end do
      end if
    else
      do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
        tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k))
      end do
    end if
  end subroutine temp_to_tvirt

  subroutine tvirt_to_temp( )
    implicit none
    integer(ik4) :: i, j, k
    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
          t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
             qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
        end do
      else
        do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
          t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - qc(j,i,k))
        end do
      end if
    else
      do concurrent( j=jce1:jce2, i = ice1:ice2, k = 1:kz )
        t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k))
      end do
    end if
  end subroutine tvirt_to_temp

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
