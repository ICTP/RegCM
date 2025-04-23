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
  real(rkx) , pointer , dimension(:,:,:) :: s => null( )
  ! nonhydrostatic term in pressure gradient force
  ! tridiagonal inversion
  real(rkx) , pointer , dimension(:,:,:) :: wx => null( )
  real(rkx) , pointer , dimension(:,:,:) :: deltaw => null( )
  real(rkx) , pointer , dimension(:,:,:) :: tkex => null( )
  real(rkx) , pointer , dimension(:,:,:) :: wz => null( )
  real(rkx) , pointer , dimension(:,:) :: mx2 => null( )
  real(rkx) , pointer , dimension(:,:) :: rmu => null( )
  real(rkx) , pointer , dimension(:,:) :: rmv => null( )
  real(rkx) , pointer , dimension(:,:,:) :: p0 => null( )
  real(rkx) , pointer , dimension(:,:,:) :: zdiv2 => null( )

  real(rkx) , pointer , dimension(:,:,:) :: ten0 => null( )
  real(rkx) , pointer , dimension(:,:,:) :: qen0 => null( )
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0 => null( )

  real(rkx) , dimension(:) , pointer :: gzitak => null( )
  real(rkx) , dimension(:) , pointer :: gzitakh => null( )
  real(rkx) , dimension(:) , pointer :: xknu => null( )
  real(rkx) , dimension(:,:) , pointer :: p2d => null( )
  real(rkx) , dimension(:,:) , pointer :: xlat => null( )
  real(rkx) , dimension(:,:) , pointer :: xlon => null( )
  real(rkx) , dimension(:,:) , pointer :: coru => null( )
  real(rkx) , dimension(:,:) , pointer :: corv => null( )
  real(rkx) , dimension(:,:) , pointer :: mu => null( )
  real(rkx) , dimension(:,:) , pointer :: hx => null( )
  real(rkx) , dimension(:,:) , pointer :: mx => null( )
  real(rkx) , dimension(:,:) , pointer :: mv => null( )
  real(rkx) , dimension(:,:) , pointer :: hy => null( )
  real(rkx) , dimension(:,:) , pointer :: ps => null( )
  real(rkx) , dimension(:,:) , pointer :: ts => null( )
  real(rkx) , dimension(:,:) , pointer :: ht => null( )
  real(rkx) , dimension(:,:,:) , pointer :: fmz => null( )
  real(rkx) , dimension(:,:,:) , pointer :: fmzf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: pai => null( )
  real(rkx) , dimension(:,:,:) , pointer :: pf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tetav => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tf => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tvirt => null( )
  real(rkx) , dimension(:,:,:) , pointer :: zeta => null( )
  real(rkx) , dimension(:,:,:) , pointer :: zetau => null( )
  real(rkx) , dimension(:,:,:) , pointer :: zetav => null( )
  real(rkx) , dimension(:,:,:) , pointer :: u => null( )
  real(rkx) , dimension(:,:,:) , pointer :: v => null( )
  real(rkx) , dimension(:,:,:) , pointer :: w => null( )
  real(rkx) , dimension(:,:,:) , pointer :: ux => null( )
  real(rkx) , dimension(:,:,:) , pointer :: vx => null( )
  real(rkx) , dimension(:,:,:) , pointer :: ud => null( )
  real(rkx) , dimension(:,:,:) , pointer :: vd => null( )
  real(rkx) , dimension(:,:,:) , pointer :: p => null( )
  real(rkx) , dimension(:,:,:) , pointer :: t => null( )
  real(rkx) , dimension(:,:,:) , pointer :: rho => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qv => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qc => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qi => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qr => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qs => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qsat => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qwltot => null( )
  real(rkx) , dimension(:,:,:) , pointer :: qwitot => null( )
  real(rkx) , dimension(:,:,:) , pointer :: tke => null( )
  real(rkx) , dimension(:,:,:,:) , pointer :: qx => null( )
  real(rkx) , dimension(:,:,:,:) , pointer :: trac => null( )

  public :: allocate_moloch , init_moloch , moloch

  real(rkx) , parameter :: minden = 1.0e-30_rkx
  real(rkx) , parameter :: xdamp = 0.0625_rkx

  logical , parameter :: do_fulleq       = .true.
  logical , parameter :: do_bdy          = .true.
  logical , parameter :: do_divdamp      = .true.
  logical , parameter :: do_vadvtwice    = .true.
  logical , parameter :: do_filterpai    = .false.
  logical , parameter :: do_filtertheta  = .false.
#ifdef RCEMIP
  logical , parameter :: do_gkfilter     = .false.
#endif
  logical , parameter :: do_phys         = .true.
  logical , parameter :: do_convection   = .true.
  logical , parameter :: do_microphysics = .true.
  logical , parameter :: do_radiation    = .true.
  logical , parameter :: do_surface      = .true.
  logical , parameter :: do_pbl          = .true.

  logical :: moloch_realcase = (.not. moloch_do_test_1) .and. &
                               (.not. moloch_do_test_2)
  logical :: lrotllr

  real(rkx) , parameter :: nupaitq = 0.05_rkx

  real(rkx) :: dzita
  integer(ik4) :: jmin , jmax , imin , imax

#ifdef RCEMIP
  ! real(rkx) , parameter :: sigma = 0.5_rkx*sqrt(2.0_rkx)
  ! do i = -2 , 2
  !   do j = -2 , 2
  !     kern(i+3,j+3) = (d_one/twopi)*exp(-(i*i+j*i)/(d_two*sigma**2))
  !   end do
  ! end do
  ! kern = kern/sum(kern)
  !real(rkx) , parameter , dimension(5,5) :: kern = reshape(&
  !       [ 1.0678874539336266E-004 , 2.1449092885792707E-003 , &
  !         5.8304679428380487E-003 , 2.1449092885792707E-003 , &
  !         1.0678874539336266E-004 , &
  !         2.1449092885792707E-003 , 4.3081654712647119E-002 , &
  !         0.1171080791453356E-000 , 4.3081654712647119E-002 , &
  !         2.1449092885792707E-003 , &
  !         5.8304679428380487E-003 , 0.1171080791453356E-000 , &
  !         0.3183327635065094e-000 , 0.1171080791453356E-000 , &
  !         5.8304679428380487E-003 , &
  !         2.1449092885792707E-003 , 4.3081654712647119E-002 , &
  !         0.1171080791453356e-000 , 4.3081654712647119E-002 , &
  !         2.1449092885792707E-003 , &
  !         1.0678874539336266E-004 , 2.1449092885792707E-003 , &
  !         5.8304679428380487E-003 , 2.1449092885792707E-003 , &
  !         1.0678874539336266E-004 ], [5,5])
  real(rkx) , parameter , dimension(5,5) :: kern = reshape(&
         [ 1.2477641543232609E-002 , 2.6415167354310421E-002 , &
           3.3917746268994867E-002 , 2.6415167354310421E-002 , &
           1.2477641543232609E-002 , &
           2.6415167354310421E-002 , 5.5920909727901751E-002 , &
           7.1803869414926613E-002 , 5.5920909727901751E-002 , &
           2.6415167354310421E-002 , &
           3.3917746268994867E-002 , 7.1803869414926613E-002 , &
           9.2197993345293314E-002 , 7.1803869414926613E-002 , &
           3.3917746268994867E-002 , &
           2.6415167354310421E-002 , 5.5920909727901751E-002 , &
           7.1803869414926613E-002 , 5.5920909727901751E-002 , &
           2.6415167354310421E-002 , &
           1.2477641543232609E-002 , 2.6415167354310421E-002 , &
           3.3917746268994867E-002 , 2.6415167354310421E-002 , &
           1.2477641543232609E-002 ], [5,5])
#endif

  contains

#include <pfesat.inc>
#include <pfwsat.inc>

  subroutine allocate_moloch
    implicit none
    integer(ik4) :: k
    call getmem1d(gzitak,1,kzp1,'moloch:gzitak')
    call getmem1d(gzitakh,1,kz,'moloch:gzitakh')
    call getmem2d(p2d,jdi1,jdi2,idi1,idi2,'moloch:p2d')
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:dw')
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
#ifdef RCEMIP
    call getmem3d(wx,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'moloch:wx')
#else
    call getmem3d(wx,jce1,jce2,ice1,ice2,1,kz,'moloch:wx')
#endif
    call getmem3d(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
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
      call getmem3d(qwltot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwltot')
      call getmem3d(qwitot,jci1,jci2,ici1,ici2,1,kz,'moloch:qwitot')
    end if
    call getmem1d(xknu,1,kz,'moloch:xknu')
    do concurrent ( k = 1:kz )
      xknu(k) = xdamp + &
        (1.0_rkx-xdamp) * sin(d_half*mathpi*(1.0_rkx-real(k-1,rkx)/kzm1))
    end do
    if ( do_filterpai ) then
      call getmem3d(pf,jce1,jce2,ice1,ice2,1,kz,'moloch:pf')
    end if
    if ( do_filtertheta ) then
      call getmem3d(tf,jce1,jce2,ice1,ice2,1,kz,'moloch:tf')
    end if
  end subroutine allocate_moloch

  subroutine init_moloch
    implicit none
    integer(ik4) :: i , j
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
    gzitak = gzita(zita,mo_ztop,mo_a0)
    gzitakh = gzita(zitah,mo_ztop,mo_a0)
    dzita = mo_dzita
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
  end subroutine init_moloch
  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch
    implicit none
    real(rkx) :: dtsound , dtstepa
    real(rkx) :: maxps , minps , pmax , pmin
    real(rkx) :: fice , zdgz , lrt , tv
    !real(rk8) :: jday
    integer(ik4) :: i , j , k , n , nadv
    integer(ik4) :: iconvec
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    dtstepa = dtsec / real(mo_nadv,rkx)
    dtsound = dtstepa / real(mo_nsound,rkx)
    iconvec = 0

    call reset_tendencies

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      p(j,i,k) = (pai(j,i,k)**cpovr) * p00
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
    end do

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
                                     qc(j,i,k) - qi(j,i,k) - &
                                     qr(j,i,k) - qs(j,i,k))
        end do
        if ( do_fulleq ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            qwltot(j,i,k) = qc(j,i,k) + qr(j,i,k)
            qwitot(j,i,k) = qi(j,i,k) + qs(j,i,k)
          end do
        end if
      else
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - qc(j,i,k))
        end do
        if ( do_fulleq ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
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
      end if
    else
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k))
      end do
      if ( do_fulleq ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          qwltot(j,i,k) = d_zero
          qwitot(j,i,k) = d_zero
        end do
      end if
    end if

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
    end do

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

    if ( do_filterpai ) then
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pf(j,i,k) = pai(j,i,k)
      end do
    end if
    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tf(j,i,k) = tetav(j,i,k)
        end do
      end if
    end if

    do nadv = 1 , mo_nadv

      call sound(dtsound)

      call advection(dtstepa)

    end do ! Advection loop

    if ( do_filterpai ) then
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pai(j,i,k) = pai(j,i,k) - pf(j,i,k)
      end do
      call filtpai
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        pai(j,i,k) = pai(j,i,k) + pf(j,i,k)
      end do
    end if

    if ( do_fulleq ) then
      if ( do_filtertheta ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tetav(j,i,k) = tetav(j,i,k) - tf(j,i,k)
        end do
        call filttheta
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          tetav(j,i,k) = tetav(j,i,k) + tf(j,i,k)
        end do
      end if
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
    end do

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
                         qc(j,i,k) - qi(j,i,k) - qr(j,i,k) - qs(j,i,k))
        end do
      else
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - qc(j,i,k))
        end do
      end if
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k))
      end do
    end if

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

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      p(j,i,k) = (pai(j,i,k)**cpovr) * p00
      rho(j,i,k) = p(j,i,k)/(rgas*t(j,i,k))
    end do

    !jday = yeardayfrac(rcmtimer%idate)

    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      zdgz = zeta(j,i,kz)*egrav
      lrt = (tvirt(j,i,kz-1)-tvirt(j,i,kz))/(zeta(j,i,kz-1)-zeta(j,i,kz))
      lrt = 0.65_rkx*lrt - 0.35_rkx*lrate
      tv = tvirt(j,i,kz) - 0.5_rkx*zeta(j,i,kz)*lrt ! Mean temperature
      ps(j,i) = p(j,i,kz) * exp(zdgz/(rgas*tv))
    end do
    !
    ! Recompute saturation
    !
    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
      qsat(j,i,k) = pfwsat(t(j,i,k),p(j,i,k))
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
    call uvstagtox(u,v,ux,vx)
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
  end subroutine moloch

  subroutine boundary
    implicit none
    logical :: do_nudge
    integer(ik4) :: i , j , k , n
    do_nudge = ( iboudy == 1 .or. iboudy >= 5 .or. iboudy == 4)
    call exchange_lrbt(ps,1,jce1,jce2,ice1,ice2)
      call exchange_lrbt(u,1,jde1,jde2,ice1,ice2,1,kz)
      call exchange_lrbt(v,1,jce1,jce2,ide1,ide2,1,kz)
      call exchange_lrbt(t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange_lrbt(qv,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
      if ( is_present_qc( ) ) then
        call exchange_lrbt(qc,1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( is_present_qi( ) ) then
        call exchange_lrbt(qi,1,jce1,jce2,ice1,ice2,1,kz)
      end if
      if ( (iboudy == 1 .or. iboudy >= 5) .and. ichem == 1 ) then
        call exchange_lrbt(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
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

      if ( iboudy == 1 .or. iboudy >= 5 ) then
        call nudge(iboudy,ps,xpsb)
        call nudge(iboudy,u,v,xub,xvb)
        call nudge(iboudy,t,xtb)
        call nudge(iboudy,qv,xqb)
        call nudge(iboudy,pai,xpaib)
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%bdy(j,i,k) = t(j,i,k) - ten0(j,i,k)
            qdiag%bdy(j,i,k) = qv(j,i,k) - qen0(j,i,k)
          end do
        end if
        if ( is_present_qc( ) ) then
          call nudge(iboudy,qc,xlb)
        end if
        if ( is_present_qi( ) ) then
          call nudge(iboudy,qi,xib)
        end if
        if ( ichem == 1 ) then
          call nudge_chi(trac)
          if ( ichdiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                            k = 1:kz, n = 1:ntr )
              cbdydiag(j,i,k,n) = trac(j,i,k,n) - chiten0(j,i,k,n)
            end do
          end if
        end if
      else if ( iboudy == 4 ) then
        call sponge(ps,xpsb)
        call sponge(u,v,xub,xvb)
        call sponge(t,xtb)
        call sponge(qv,xqb)
        call sponge(pai,xpaib)
        if ( is_present_qc( ) ) then
          call sponge(qc,xlb)
        end if
        if ( is_present_qi( ) ) then
          call sponge(qi,xib)
        end if
        if ( idiag > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            tdiag%bdy(j,i,k) = t(j,i,k) - ten0(j,i,k)
            qdiag%bdy(j,i,k) = qv(j,i,k) - qen0(j,i,k)
          end do
        end if
      end if
    end subroutine boundary

#ifdef RCEMIP
    subroutine gkfilter2d(f,irep)
      implicit none
      integer(ik4) , intent(in) :: irep
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: f
      real(rkx) , allocatable , dimension(:,:) :: temp
      integer(ik4) :: i , j , n

      call exchange(f,2,jce1,jce2,ice1,ice2)
      allocate(temp, mold=f)
      temp(:,:) = f(:,:)
      do n = 1 , irep
        do concurrent( j = jcii1:jcii2, i = icii1:icii2 )
          temp(j,i) = kern(1,1)*f(j-2,i-2) + &
            kern(2,1)*f(j-1,i-2) + kern(3,1)*f(j,  i-2) + &
            kern(4,1)*f(j+1,i-2) + kern(5,1)*f(j+2,i-2) + &
            kern(1,2)*f(j-2,i-1) + kern(2,2)*f(j-1,i-1) + &
            kern(3,2)*f(j,  i-1) + kern(4,2)*f(j+1,i-1) + &
            kern(5,2)*f(j+2,i-1) + kern(1,3)*f(j-2,i  ) + &
            kern(2,3)*f(j-1,i  ) + kern(3,3)*f(j,  i  ) + &
            kern(4,3)*f(j+1,i  ) + kern(5,3)*f(j+2,i  ) + &
            kern(1,4)*f(j-2,i+1) + kern(2,4)*f(j-1,i+1) + &
            kern(3,4)*f(j,  i+1) + kern(4,4)*f(j+1,i+1) + &
            kern(5,4)*f(j+2,i+1) + kern(1,5)*f(j-2,i+2) + &
            kern(2,5)*f(j-1,i+2) + kern(3,5)*f(j,  i+2) + &
            kern(4,5)*f(j+1,i+2) + kern(5,5)*f(j+2,i+2)
        end do
      end do
      do concurrent( j = jcii1:jcii2, i = icii1:icii2 )
        f(j,i) = temp(j,i)
      end do
      deallocate(temp)
    end subroutine gkfilter2d

    subroutine gkfilter3d(f,irep)
      implicit none
      integer(ik4) , intent(in) :: irep
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: f
      real(rkx) , allocatable , dimension(:,:,:) :: temp
      integer(ik4) :: i , j , k , n

      call exchange(f,2,jce1,jce2,ice1,ice2,1,kz)
      allocate(temp, mold=f)
      temp(:,:,:) = f(:,:,:)
      do n = 1 , irep
        do concurrent( j = jcii1:jcii2, i = icii1:icii2, k = 1:kz )
          temp(j,i,k) = kern(1,1)*f(j-2,i-2,k) + &
            kern(2,1)*f(j-1,i-2,k) + kern(3,1)*f(j,  i-2,k) + &
            kern(4,1)*f(j+1,i-2,k) + kern(5,1)*f(j+2,i-2,k) + &
            kern(1,2)*f(j-2,i-1,k) + kern(2,2)*f(j-1,i-1,k) + &
            kern(3,2)*f(j,  i-1,k) + kern(4,2)*f(j+1,i-1,k) + &
            kern(5,2)*f(j+2,i-1,k) + kern(1,3)*f(j-2,i  ,k) + &
            kern(2,3)*f(j-1,i  ,k) + kern(3,3)*f(j,  i  ,k) + &
            kern(4,3)*f(j+1,i  ,k) + kern(5,3)*f(j+2,i  ,k) + &
            kern(1,4)*f(j-2,i+1,k) + kern(2,4)*f(j-1,i+1,k) + &
            kern(3,4)*f(j,  i+1,k) + kern(4,4)*f(j+1,i+1,k) + &
            kern(5,4)*f(j+2,i+1,k) + kern(1,5)*f(j-2,i+2,k) + &
            kern(2,5)*f(j-1,i+2,k) + kern(3,5)*f(j,  i+2,k) + &
            kern(4,5)*f(j+1,i+2,k) + kern(5,5)*f(j+2,i+2,k)
        end do
      end do
      do concurrent( j = jcii1:jcii2, i = icii1:icii2, k = 1:kz )
        f(j,i,k) = temp(j,i,k)
      end do
      deallocate(temp)
    end subroutine gkfilter3d

    subroutine gkfilter4d(f,irep)
      implicit none
      integer(ik4) , intent(in) :: irep
      real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: f
      real(rkx) , allocatable , dimension(:,:,:,:) :: temp
      integer(ik4) :: i , j , k , n , m
      call exchange(f,2,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      allocate(temp, mold=f)
      temp(:,:,:,:) = f(:,:,:,:)
      do m = 1 , irep
        do concurrent(j=jcii1:jcii2, i=icii1:icii2, k=1:kz , n=1:nqx)
          temp(j,i,k,n) = kern(1,1)*f(j-2,i-2,k,n) + &
            kern(2,1)*f(j-1,i-2,k,n) + kern(3,1)*f(j,  i-2,k,n) + &
            kern(4,1)*f(j+1,i-2,k,n) + kern(5,1)*f(j+2,i-2,k,n) + &
            kern(1,2)*f(j-2,i-1,k,n) + kern(2,2)*f(j-1,i-1,k,n) + &
            kern(3,2)*f(j,  i-1,k,n) + kern(4,2)*f(j+1,i-1,k,n) + &
            kern(5,2)*f(j+2,i-1,k,n) + kern(1,3)*f(j-2,i  ,k,n) + &
            kern(2,3)*f(j-1,i  ,k,n) + kern(3,3)*f(j,  i  ,k,n) + &
            kern(4,3)*f(j+1,i  ,k,n) + kern(5,3)*f(j+2,i  ,k,n) + &
            kern(1,4)*f(j-2,i+1,k,n) + kern(2,4)*f(j-1,i+1,k,n) + &
            kern(3,4)*f(j,  i+1,k,n) + kern(4,4)*f(j+1,i+1,k,n) + &
            kern(5,4)*f(j+2,i+1,k,n) + kern(1,5)*f(j-2,i+2,k,n) + &
            kern(2,5)*f(j-1,i+2,k,n) + kern(3,5)*f(j,  i+2,k,n) + &
            kern(4,5)*f(j+1,i+2,k,n) + kern(5,5)*f(j+2,i+2,k,n)
        end do
      end do
      do concurrent(j=jcii1:jcii2, i=icii1:icii2, k=1:kz, n=1:nqx)
        f(j,i,k,n) = temp(j,i,k,n)
      end do
      deallocate(temp)
    end subroutine gkfilter4d

    subroutine gkfilteruv(u,v,irep)
      implicit none
      integer(ik4) , intent(in) :: irep
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
      real(rkx) , allocatable , dimension(:,:,:) :: temp
      integer(ik4) :: i , j , k , n

      call exchange(u,2,jce1,jce2,ice1,ice2,1,kz)
      allocate(temp, mold=u)
      temp(:,:,:) = u(:,:,:)
      do n = 1 , irep
        do concurrent( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          temp(j,i,k) = kern(1,1)*u(j-2,i-2,k) + &
            kern(2,1)*u(j-1,i-2,k) + kern(3,1)*u(j,  i-2,k) + &
            kern(4,1)*u(j+1,i-2,k) + kern(5,1)*u(j+2,i-2,k) + &
            kern(1,2)*u(j-2,i-1,k) + kern(2,2)*u(j-1,i-1,k) + &
            kern(3,2)*u(j,  i-1,k) + kern(4,2)*u(j+1,i-1,k) + &
            kern(5,2)*u(j+2,i-1,k) + kern(1,3)*u(j-2,i  ,k) + &
            kern(2,3)*u(j-1,i  ,k) + kern(3,3)*u(j,  i  ,k) + &
            kern(4,3)*u(j+1,i  ,k) + kern(5,3)*u(j+2,i  ,k) + &
            kern(1,4)*u(j-2,i+1,k) + kern(2,4)*u(j-1,i+1,k) + &
            kern(3,4)*u(j,  i+1,k) + kern(4,4)*u(j+1,i+1,k) + &
            kern(5,4)*u(j+2,i+1,k) + kern(1,5)*u(j-2,i+2,k) + &
            kern(2,5)*u(j-1,i+2,k) + kern(3,5)*u(j,  i+2,k) + &
            kern(4,5)*u(j+1,i+2,k) + kern(5,5)*u(j+2,i+2,k)
        end do
      end do
      do concurrent( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        u(j,i,k) = temp(j,i,k)
      end do
      deallocate(temp)
      call exchange(v,2,jce1,jce2,ice1,ice2,1,kz)
      allocate(temp, mold=v)
      temp(:,:,:) = v(:,:,:)
      do n = 1 , irep
        do concurrent( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          temp(j,i,k) = kern(1,1)*v(j-2,i-2,k) + &
            kern(2,1)*v(j-1,i-2,k) + kern(3,1)*v(j,  i-2,k) + &
            kern(4,1)*v(j+1,i-2,k) + kern(5,1)*v(j+2,i-2,k) + &
            kern(1,2)*v(j-2,i-1,k) + kern(2,2)*v(j-1,i-1,k) + &
            kern(3,2)*v(j,  i-1,k) + kern(4,2)*v(j+1,i-1,k) + &
            kern(5,2)*v(j+2,i-1,k) + kern(1,3)*v(j-2,i  ,k) + &
            kern(2,3)*v(j-1,i  ,k) + kern(3,3)*v(j,  i  ,k) + &
            kern(4,3)*v(j+1,i  ,k) + kern(5,3)*v(j+2,i  ,k) + &
            kern(1,4)*v(j-2,i+1,k) + kern(2,4)*v(j-1,i+1,k) + &
            kern(3,4)*v(j,  i+1,k) + kern(4,4)*v(j+1,i+1,k) + &
            kern(5,4)*v(j+2,i+1,k) + kern(1,5)*v(j-2,i+2,k) + &
            kern(2,5)*v(j-1,i+2,k) + kern(3,5)*v(j,  i+2,k) + &
            kern(4,5)*v(j+1,i+2,k) + kern(5,5)*v(j+2,i+2,k)
        end do
      end do
      do concurrent( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        v(j,i,k) = temp(j,i,k)
      end do
      deallocate(temp)
    end subroutine gkfilteruv

    subroutine filt3d(p,nu)
      implicit none
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: p
      real(rkx) , intent(in) :: nu
      integer(ik4) :: j , i , k

      do k = 1 , kz
        do concurrent ( j = jcii1:jcii2, i = icii1:icii2 )
          p2d(j,i) = 0.125_rkx * (p(j-1,i,k) + p(j+1,i,k) + &
                                  p(j,i-1,k) + p(j,i+1,k)) - &
                       d_half   * p(j,i,k)
        end do
        do concurrent ( j = jcii1:jcii2, i = icii1:icii2 )
          p(j,i,k) = p(j,i,k) + nu * p2d(j,i)
        end do
      end do
    end subroutine filt3d

    subroutine filtuv(u,v,nu)
      implicit none
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
      real(rkx) , intent(in) :: nu
      integer(ik4) :: j , i , k

      do k = 1 , kz
        do concurrent ( j = jdii1:jdii2, i = icii1:icii2 )
          p2d(j,i) = 0.125_rkx * (u(j-1,i,k) + u(j+1,i,k) + &
                                  u(j,i-1,k) + u(j,i+1,k)) - &
                       d_half   * u(j,i,k)
        end do
        do concurrent ( j = jdii1:jdii2, i = icii1:icii2 )
          u(j,i,k) = u(j,i,k) + nu * p2d(j,i)
        end do
      end do
      do k = 1 , kz
        do concurrent ( j = jcii1:jcii2, i = idii1:idii2 )
          p2d(j,i) = 0.125_rkx * (v(j-1,i,k) + v(j+1,i,k) + &
                                  v(j,i-1,k) + v(j,i+1,k)) - &
                       d_half   * v(j,i,k)
        end do
        do concurrent ( j = jcii1:jcii2, i = idii1:idii2 )
          v(j,i,k) = v(j,i,k) + nu * p2d(j,i)
        end do
      end do
    end subroutine filtuv

    subroutine filt4d(p,nu,n1,n2)
      implicit none
      real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: p
      real(rkx) , intent(in) :: nu
      integer(ik4) , intent(in) :: n1 , n2
      integer(ik4) :: j , i , k , n

      do n = n1 , n2
        do k = 1 , kz
          do concurrent ( j = jcii1:jcii2, i = icii1:icii2 )
            p2d(j,i) = 0.125_rkx * (p(j-1,i,k,n) + p(j+1,i,k,n) + &
                                    p(j,i-1,k,n) + p(j,i+1,k,n)) - &
                           d_half * p(j,i,k,n)
          end do
          do concurrent ( j = jcii1:jcii2, i = icii1:icii2 )
            p(j,i,k,n) = p(j,i,k,n) + nu * p2d(j,i)
          end do
        end do
      end do
    end subroutine filt4d
#endif

    subroutine divergence_filter( )
      implicit none
      integer(ik4) :: j , i , k
      if ( ma%has_bdybottom ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          zdiv2(j,ice1,k) = zdiv2(j,ici1,k)
        end do
      end if
      if ( ma%has_bdytop ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          zdiv2(j,ice2,k) = zdiv2(j,ici2,k)
        end do
      end if
      if ( ma%has_bdyleft ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          zdiv2(jce1,i,k) = zdiv2(jci1,i,k)
        end do
      end if
      if ( ma%has_bdyright ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          zdiv2(jce2,i,k) = zdiv2(jci2,i,k)
        end do
      end if
      call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)
      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                  zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                       d_half   * zdiv2(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          zdiv2(j,i,k) = zdiv2(j,i,k) + mo_anu2 * xknu(k) * p2d(j,i)
        end do
      end do
    end subroutine divergence_filter

    subroutine filtpai
      implicit none
      integer(ik4) :: j , i , k

      call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)

      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (pai(j-1,i,k) + pai(j+1,i,k) + &
                                  pai(j,i-1,k) + pai(j,i+1,k)) - &
                       d_half   * pai(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          pai(j,i,k) = pai(j,i,k) + nupaitq * p2d(j,i)
        end do
      end do
    end subroutine filtpai

    subroutine filttheta
      implicit none
      integer(ik4) :: j , i , k

      call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)

      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (tetav(j-1,i,k) + tetav(j+1,i,k) + &
                                  tetav(j,i-1,k) + tetav(j,i+1,k)) - &
                       d_half   * tetav(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          tetav(j,i,k) = tetav(j,i,k) + nupaitq * p2d(j,i)
        end do
      end do
    end subroutine filttheta

    subroutine sound(dts)
      implicit none
      real(rkx) , intent(in) :: dts
      integer(ik4) :: i , j , k , nsound
      real(rkx) :: dtrdx , dtrdy , dtrdz , zcs2
      real(rkx) , dimension(jci1:jci2,ici1:ici2,2:kzp1) :: wwkw
      real(rkx) :: zrfmzum , zrfmzvm , zrfmzup , zrfmzvp
      real(rkx) :: zum , zup , zvm , zvp , zuh , zvh
      real(rkx) :: zrom1w , zwexpl , zqs , zdth , zu , zd , zrapp
      real(rkx) :: zcx , zcy , zfz
      real(rkx) :: zrom1u , zcor1u , zrom1v , zcor1v

      dtrdx = dts/dx
      dtrdy = dts/dx
      dtrdz = dts/dzita
      zcs2 = dtrdz**2*rdrcv

      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        wwkw(j,i,kzp1) = d_zero
      end do

      !  sound waves

      if ( .not. do_fulleq ) then
        call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
      end if

      do nsound = 1 , mo_nsound

        call exchange(u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange(v,1,jce1,jce2,ide1,ide2,1,kz)

        ! partial definition of the generalized vertical velocity

        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
          zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
          w(j,i,kzp1) = d_half * (zuh+zvh)
        end do

        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          s(j,i,kzp1) = -w(j,i,kzp1)
        end do

        ! Equation 10, generalized vertical velocity

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
          zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) +    &
                (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
          zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) +    &
                (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
          s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak(k)
        end do

        ! Part of divergence (except w contribution) put in zdiv2
        ! Equation 16

        if ( lrotllr ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
            zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
            zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
            zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
            zum = dtrdx * u(j,i,k) * zrfmzum
            zup = dtrdx * u(j+1,i,k) * zrfmzup
            zvm = dtrdy * v(j,i,k) * zrfmzvm * rmv(j,i)
            zvp = dtrdy * v(j,i+1,k) * zrfmzvp * rmv(j,i+1)
            zdiv2(j,i,k) = fmz(j,i,k) * mx(j,i) * ((zup-zum) + (zvp-zvm))
          end do
        else
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            zrfmzum = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
            zrfmzvm = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
            zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
            zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
            zum = dtrdx * u(j,i,k)   * rmu(j,i)   * zrfmzum
            zup = dtrdx * u(j+1,i,k) * rmu(j+1,i) * zrfmzup
            zvm = dtrdy * v(j,i,k)   * rmv(j,i)   * zrfmzvm
            zvp = dtrdy * v(j,i+1,k) * rmv(j,i+1) * zrfmzvp
            zdiv2(j,i,k) = mx2(j,i) * fmz(j,i,k) * ((zup-zum)+(zvp-zvm))
          end do
        end if

        if ( do_divdamp ) then
          call divdamp(dts)
        end if

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zdiv2(j,i,k) = zdiv2(j,i,k) + fmz(j,i,k) * &
                     dtrdz * (s(j,i,k) - s(j,i,k+1))
        end do

        ! new w (implicit scheme) from Equation 19

        do k = kz , 2 , -1
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            deltaw(j,i,k) = -w(j,i,k)
            ! explicit w:
            !    it must be consistent with the initialization of pai
            zrom1w = d_half * cpd * fmzf(j,i,k) * &
                     (tetav(j,i,k-1)+tetav(j,i,k))
            zrom1w = zrom1w - cpd * w(j,i,k) * &
                     fmzf(j,i,k)*fmzf(j,i,k) * &
                     real(nsound,rkx) * dtrdz * &
                     (tetav(j,i,k-1)-tetav(j,i,k)) !! GW
            if ( qv(j,i,k) > 0.96_rkx*qsat(j,i,k) .and. &
                  w(j,i,k) > 0.1_rkx ) then
              zqs = d_half*(qsat(j,i,k)+qsat(j,i,k-1))
              zdth = egrav*w(j,i,k)*real(nsound-1,rkx)*dts*wlhv*wlhv* &
                zqs/(cpd*pai(j,i,k-1)*rwat*t(j,i,k-1)*t(j,i,k-1))
              zrom1w = zrom1w + zdth*fmzf(j,i,k)
            end if
            ! explicit w:
            !    it must be consistent with the initialization of pai
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
        do k = 2 , kz
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
            deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
          end do
        end do

        ! new Exner function (Equation 19)

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zdiv2(j,i,k) = zdiv2(j,i,k) + dtrdz * fmz(j,i,k) * &
                    (w(j,i,k) - w(j,i,k+1))
        end do

        if ( do_fulleq ) then
          if ( ipptls > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
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
          end if
          call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)
        end if

        if ( mo_divfilter ) call divergence_filter( )

        ! horizontal momentum equations
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv2(j,i,k))
        end do

        call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kzp1)

        do concurrent ( j = jde1ga:jde2ga , i = ice1ga:ice2ga , k = 1:kz )
          ud(j,i,k) = u(j,i,k)
        end do
        do concurrent ( j = jce1ga:jce2ga , i = ide1ga:ide2ga , k = 1:kz )
          vd(j,i,k) = v(j,i,k)
        end do

        if ( lrotllr ) then
          ! Equation 17
          do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
            zcx = dtrdx * mu(j,i)
            zfz = egrav * dts + 0.25_rkx * &
                (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                 deltaw(j,i,k)   + deltaw(j,i,k+1))
            zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
            zcor1u = coru(j,i) * dts * 0.25_rkx * &
                 (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
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
            zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
            zcor1v = corv(j,i) * dts * 0.25_rkx * &
                   (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
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
            zrom1u = d_half * cpd * (tetav(j-1,i,k) + tetav(j,i,k))
            zcor1u = coru(j,i) * dts * 0.25_rkx * &
                 (vd(j,i,k) + vd(j-1,i,k) + vd(j-1,i+1,k) + vd(j,i+1,k))
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
            zrom1v = d_half * cpd * (tetav(j,i-1,k) + tetav(j,i,k))
            zcor1v = corv(j,i) * dts * 0.25_rkx * &
                   (ud(j,i,k) + ud(j,i-1,k) + ud(j+1,i,k) + ud(j+1,i-1,k))
            ! Equation 18
            v(j,i,k) = v(j,i,k) - zcor1v - &
                       zfz * hy(j,i) * gzitakh(k) -  &
                       zcy * zrom1v * (pai(j,i,k) - pai(j,i-1,k))
          end do
        end if

      end do ! sound loop

      ! complete computation of generalized vertical velocity
      ! Complete Equation 10
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        s(j,i,1) = d_zero
        s(j,i,kzp1) = d_zero
      end do

    end subroutine sound

    subroutine divdamp(dts)
      implicit none
      real(rkx) , intent(in) :: dts
      integer(ik4) :: i , j , k
      real(rkx) :: ddamp

      if ( ma%has_bdybottom ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          zdiv2(j,ice1,k) = zdiv2(j,ici1,k)
        end do
      end if
      if ( ma%has_bdytop ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          zdiv2(j,ice2,k) = zdiv2(j,ici2,k)
        end do
      end if
      if ( ma%has_bdyleft ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          zdiv2(jce1,i,k) = zdiv2(jci1,i,k)
        end do
      end if
      if ( ma%has_bdyright ) then
        do concurrent ( i = ici1:ici2, k = 1:kz )
          zdiv2(jce2,i,k) = zdiv2(jci2,i,k)
        end do
      end if
      call exchange_lrbt(zdiv2,1,jce1,jce2,ice1,ice2,1,kz)

      ddamp = 0.125_rkx * (dx/dts)
      if ( lrotllr ) then
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          u(j,i,k) = u(j,i,k) + &
                  xknu(k)*ddamp*mu(j,i)*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          v(j,i,k) = v(j,i,k) + &
                  xknu(k)*ddamp*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
        end do
      else
        do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
          u(j,i,k) = u(j,i,k) + &
                  xknu(k)*ddamp*mu(j,i)*(zdiv2(j,i,k)-zdiv2(j-1,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
          v(j,i,k) = v(j,i,k) + &
                  xknu(k)*ddamp*mv(j,i)*(zdiv2(j,i,k)-zdiv2(j,i-1,k))
        end do
      end if
      ! Horizontal diffusion
      ddamp = xdamp * 0.015625_rkx/dts
      do k = 1 , kz
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          p2d(j,i) = 0.125_rkx * (zdiv2(j-1,i,k) + zdiv2(j+1,i,k) + &
                                  zdiv2(j,i-1,k) + zdiv2(j,i+1,k)) - &
                       d_half   * zdiv2(j,i,k)
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          zdiv2(j,i,k) = zdiv2(j,i,k) + ddamp * p2d(j,i)
        end do
      end do
    end subroutine divdamp

    subroutine advection(dta)
      implicit none
      real(rkx) , intent(in) :: dta
      integer(ik4) :: n
      real(rkx) , pointer , dimension(:,:,:) :: ptr => null( )

      ! Compute U,V on cross points

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
        do n = iqfrst , nqx
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

#ifdef RCEMIP
      if ( do_gkfilter ) then
        call gkfilteruv(ux,vx,1)
        call gkfilter3d(wx,1)
      end if
#endif
      ! Interpolate on staggered points
      call xtouvstag(ux,vx,u,v)

      ! Back to half-levels
      call xtowstag(wx,w)
      if ( ibltyp == 2 ) then
        call xtowstag(tkex,tke)
      end if
    end subroutine advection

    subroutine wafone(pp,dta,pfac,pmin)
      implicit none
      real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
      real(rkx) , intent(in) :: dta
      real(rkx) , optional , intent(in) :: pfac , pmin
      integer(ik4) :: j , i , k
      real(rkx) :: dtrdx , dtrdy , dtrdz
      real(rkx) , parameter :: wlow  = 0.0_rkx
      real(rkx) , parameter :: whigh = 2.0_rkx
      real(rkx) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: wfw
      real(rkx) , dimension(jci1:jci2,ici1:ice2ga,1:kz) :: zpby
      real(rkx) , dimension(jci1:jce2ga,ici1:ici2,1:kz) :: zpbw
      real(rkx) :: zamu , is , r , b , zphi , zzden , zdv
      real(rkx) :: zhxvtn , zhxvts , zcostx
      real(rkx) :: zrfmu , zrfmd
      real(rkx) :: zrfmn , zrfms
      real(rkx) :: zrfme , zrfmw
      real(rkx) :: zfac , zmin
      integer(ik4) :: k1 , k1p1
      integer(ik4) :: ih , ihm1
      integer(ik4) :: jh , jhm1

      dtrdx = dta/dx
      dtrdy = dta/dx
      dtrdz = dta/dzita
      if ( do_vadvtwice ) then
        dtrdz = 0.5_rkx * dtrdz
      end if

      if ( present(pfac) ) then
        zfac = pfac
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = pp(j,i,k) * zfac
        end do
      end if

      ! Vertical advection
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 )
        wfw(j,i,1) = d_zero
        wfw(j,i,kzp1) = d_zero
      end do

      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 )
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
        wfw(j,i,k+1) = d_half * s(j,i,k+1) * ((d_one+zphi)*pp(j,i,k+1) + &
                                              (d_one-zphi)*pp(j,i,k))
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
        zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
        zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * pp(j,i,k)
        wz(j,i,k) = pp(j,i,k) - &
          wfw(j,i,k)*zrfmu + wfw(j,i,k+1)*zrfmd + zdv
      end do

      if ( do_vadvtwice ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzm1 )
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
          wfw(j,i,k+1) = d_half * s(j,i,k+1) * &
            ((d_one+zphi)*wz(j,i,k+1) + (d_one-zphi)*wz(j,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zrfmu = dtrdz * fmz(j,i,k)/fmzf(j,i,k)
          zrfmd = dtrdz * fmz(j,i,k)/fmzf(j,i,k+1)
          zdv = (s(j,i,k)*zrfmu - s(j,i,k+1)*zrfmd) * wz(j,i,k)
          wz(j,i,k) = wz(j,i,k) - wfw(j,i,k)*zrfmu + &
                                  wfw(j,i,k+1)*zrfmd + zdv
        end do
      end if

      if ( ma%has_bdybottom ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          wz(j,ice1,k) = wz(j,ici1,k)
        end do
      end if
      if ( ma%has_bdytop ) then
        do concurrent ( j = jci1:jci2, k = 1:kz )
          wz(j,ice2,k) = wz(j,ici2,k)
        end do
      end if

      call exchange_bt(wz,2,jci1,jci2,ice1,ice2,1,kz)

      if ( lrotllr ) then

        ! Meridional advection
        do concurrent ( j = jci1:jci2, i = ici1:ice2ga, k = 1:kz )
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
          zpby(j,i,k) = d_half * v(j,i,k) * &
              ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zhxvtn = dtrdy * rmv(j,i+1) * mx(j,i)
          zhxvts = dtrdy * rmv(j,i) * mx(j,i)
          zrfmn = zhxvtn * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
          zrfms = zhxvts * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
          zdv = (v(j,i+1,k) * zrfmn - v(j,i,k) * zrfms) * pp(j,i,k)
          p0(j,i,k) = wz(j,i,k) + &
                zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv
        end do

        if ( ma%has_bdyleft ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            p0(jce1,i,k) = p0(jci1,i,k)
          end do
        end if

        if ( ma%has_bdyright ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            p0(jce2,i,k) = p0(jci2,i,k)
          end do
        end if

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
          zpbw(j,i,k) = d_half * u(j,i,k) * &
                 ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zcostx = dtrdx * mu(j,i)
          zrfme = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
          zrfmw = zcostx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
          zdv = (u(j+1,i,k) * zrfme - u(j,i,k) * zrfmw) * pp(j,i,k)
          pp(j,i,k) = p0(j,i,k) + &
                 zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv
        end do

      else

        ! Meridional advection

        do concurrent ( j = jci1:jci2, i = ici1:ice2ga, k = 1:kz )
          zamu = v(j,i,k) * rmv(j,i) * dtrdy
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
          zpby(j,i,k) = d_half * v(j,i,k) * rmv(j,i) * &
              ((d_one+zphi)*wz(j,i-1,k) + (d_one-zphi)*wz(j,i,k))
        end do
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zrfmn = dtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i+1,k))
          zrfms = dtrdy * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j,i-1,k))
          zdv = (v(j,i+1,k) * rmv(j,i+1) * zrfmn - &
                 v(j,i,k)   * rmv(j,i)   * zrfms) * pp(j,i,k)
          p0(j,i,k) = wz(j,i,k) + &
            mx2(j,i) * (zpby(j,i,k)*zrfms - zpby(j,i+1,k)*zrfmn + zdv)
        end do

        if ( ma%has_bdyleft ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            p0(jce1,i,k) = p0(jci1,i,k)
          end do
        end if

        if ( ma%has_bdyright ) then
          do concurrent ( i = ici1:ici2, k = 1:kz )
            p0(jce2,i,k) = p0(jci2,i,k)
          end do
        end if

        call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

        ! Zonal advection

        do concurrent ( j = jci1:jce2ga, i = ici1:ici2, k = 1:kz )
          zamu = u(j,i,k) * rmu(j,i) * dtrdx
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
          zpbw(j,i,k) = d_half * u(j,i,k) * rmu(j,i) * &
                 ((d_one+zphi)*p0(j-1,i,k) + (d_one-zphi)*p0(j,i,k))
        end do

        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          zrfme = dtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j+1,i,k))
          zrfmw = dtrdx * d_two * fmz(j,i,k)/(fmz(j,i,k)+fmz(j-1,i,k))
          zdv = (u(j+1,i,k) * rmu(j+1,i) * zrfme - &
                 u(j,i,k)   * rmu(j,i)   * zrfmw) * pp(j,i,k)
          pp(j,i,k) = p0(j,i,k) + &
              mx2(j,i) * (zpbw(j,i,k)*zrfmw - zpbw(j+1,i,k)*zrfme + zdv)
        end do
      end if

      if ( present(pfac) ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = pp(j,i,k) / zfac
        end do
      end if
      if ( present(pmin) ) then
        zmin = pmin
        do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
          pp(j,i,k) = max(pp(j,i,k),zmin)
        end do
      end if
    end subroutine wafone

    subroutine reset_tendencies
      implicit none
      integer(ik4) :: i , j , k , n
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
        s(j,i,k) = d_zero
      end do
      do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kz )
        zdiv2(j,i,k) = d_zero
      end do
      do concurrent ( j = jce1ga:jce2ga, i = ice1ga:ice2ga, k = 1:kzp1 )
        deltaw(j,i,k) = d_zero
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        mo_atm%tten(j,i,k) = d_zero
        mo_atm%uten(j,i,k) = d_zero
        mo_atm%vten(j,i,k) = d_zero
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:nqx )
        mo_atm%qxten(j,i,k,n) = d_zero
      end do
      if ( ichem == 1 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
          mo_atm%chiten(j,i,k,n) = d_zero
        end do
      end if
      if ( ibltyp == 2 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
          mo_atm%tketen(j,i,k) = d_zero
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
      if ( do_convection ) then
        if ( any(icup > 0) ) then
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              ten0(j,i,k) = mo_atm%tten(j,i,k)
              qen0(j,i,k) = mo_atm%qxten(j,i,k,iqv)
            end do
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                            k = 1:kz, n = 1:ntr )
              chiten0(j,i,k,n) = mo_atm%chiten(j,i,k,n)
            end do
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
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              tdiag%con(j,i,k) = mo_atm%tten(j,i,k) - ten0(j,i,k)
              qdiag%con(j,i,k) = mo_atm%qxten(j,i,k,iqv) - qen0(j,i,k)
            end do
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                            k = 1:kz, n = 1:ntr )
              cconvdiag(j,i,k,n) = mo_atm%chiten(j,i,k,n) - chiten0(j,i,k,n)
            end do
          end if
        else
          if ( any(icup < 0) ) then
            call shallow_convection
            if ( idiag > 0 ) then
              do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
                tdiag%con(j,i,k) = mo_atm%tten(j,i,k) - ten0(j,i,k)
                qdiag%con(j,i,k) = mo_atm%qxten(j,i,k,iqv) - qen0(j,i,k)
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
              ten0(j,i,k) = mo_atm%tten(j,i,k)
              qen0(j,i,k) = mo_atm%qxten(j,i,k,iqv)
            end do
          end if
          ! Cumulus clouds
          if ( icldfrac /= 2 ) then
            call cucloud
          end if
          ! Save cumulus cloud fraction for chemistry before it is
          ! overwritten in cldfrac
          if ( ichem == 1 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              convcldfra(j,i,k) = cldfra(j,i,k)
            end do
          end if
          ! Clouds and large scale precipitation
          call cldfrac(cldlwc,cldfra)
          call microscheme
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              tdiag%lsc(j,i,k) = mo_atm%tten(j,i,k) - ten0(j,i,k)
              qdiag%lsc(j,i,k) = mo_atm%qxten(j,i,k,iqv) - qen0(j,i,k)
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
          call surface_albedo
          if ( iclimao3 == 1 ) then
            call updateo3(rcmtimer%idate,scenario)
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
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          mo_atm%tten(j,i,k) = mo_atm%tten(j,i,k) + heatrt(j,i,k)
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
          call surface_model
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
              ten0(j,i,k) = mo_atm%tten(j,i,k)
              qen0(j,i,k) = mo_atm%qxten(j,i,k,iqv)
            end do
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                            k = 1:kz, n = 1:ntr )
              chiten0(j,i,k,n) = mo_atm%chiten(j,i,k,n)
            end do
          end if
          call pblscheme
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              tdiag%tbl(j,i,k) = mo_atm%tten(j,i,k) - ten0(j,i,k)
              qdiag%tbl(j,i,k) = mo_atm%qxten(j,i,k,iqv) - qen0(j,i,k)
            end do
          end if
          if ( ichem == 1 .and. ichdiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                            k = 1:kz, n = 1:ntr )
              ctbldiag(j,i,k,n) = mo_atm%chiten(j,i,k,n) - chiten0(j,i,k,n)
            end do
          end if
        end if
      end if

      if ( do_microphysics ) then
        if ( ipptls == 1 ) then
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              ten0(j,i,k) = mo_atm%tten(j,i,k)
              qen0(j,i,k) = mo_atm%qxten(j,i,k,iqv)
            end do
          end if
          call condtq
          if ( idiag > 0 ) then
            do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
              tdiag%lsc(j,i,k) = tdiag%lsc(j,i,k) + &
                mo_atm%tten(j,i,k) - ten0(j,i,k)
              qdiag%lsc(j,i,k) = qdiag%lsc(j,i,k) + &
                 mo_atm%qxten(j,i,k,iqv) - qen0(j,i,k)
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
      !
      ! Update status
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        t(j,i,k) = t(j,i,k) + dtsec * mo_atm%tten(j,i,k)
      end do
      do concurrent ( j = jdi1:jdi2, i = ici1:ici2, k = 1:kz )
        u(j,i,k) = u(j,i,k) + dtsec * mo_atm%uten(j,i,k)
      end do
      do concurrent ( j = jci1:jci2, i = idi1:idi2, k = 1:kz )
        v(j,i,k) = v(j,i,k) + dtsec * mo_atm%vten(j,i,k)
      end do
#ifdef RCEMIP
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        if ( k > ktrop(j,i) ) then
          qx(j,i,k,iqv) = qx(j,i,k,iqv) + mo_atm%qxten(j,i,k,iqv)*dtsec
          qx(j,i,k,iqv) = max(qx(j,i,k,iqv),minqq)
        end if
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                      k = 1:kz, n = iqfrst:nqx)
        if ( k > ktrop(j,i) ) then
          qx(j,i,k,n) = qx(j,i,k,n) + mo_atm%qxten(j,i,k,n)*dtsec
          qx(j,i,k,n) = max(qx(j,i,k,n),d_zero)
        end if
      end do
#else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        qx(j,i,k,iqv) = qx(j,i,k,iqv) + mo_atm%qxten(j,i,k,iqv)*dtsec
        qx(j,i,k,iqv) = max(qx(j,i,k,iqv),minqq)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, &
                      k = 1:kz, n = iqfrst:nqx)
        qx(j,i,k,n) = qx(j,i,k,n) + mo_atm%qxten(j,i,k,n)*dtsec
        qx(j,i,k,n) = max(qx(j,i,k,n),d_zero)
      end do
#endif
      if ( ibltyp == 2 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
          tke(j,i,k) = max(tke(j,i,k) + dtsec * mo_atm%tketen(j,i,k),tkemin)
        end do
      end if
      if ( ichem == 1 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:ntr )
          trac(j,i,k,n) = trac(j,i,k,n) + dtsec * mo_atm%chiten(j,i,k,n)
          trac(j,i,k,n) = max(trac(j,i,k,n), d_zero)
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

  subroutine wstagtox(w,wx)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: w
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    integer(ik4) :: i , j , k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 2:kzm1 )
      wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                  0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      wx(j,i,1)  = d_half * (w(j,i,2)+w(j,i,1))
      wx(j,i,kz) = d_half * (w(j,i,kzp1)+w(j,i,kz))
    end do
  end subroutine wstagtox

  subroutine xtowstag(wx,w)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: wx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: w
    integer(ik4) :: i , j , k

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 3:kzm1 )
      w(j,i,k) = 0.5625_rkx * (wx(j,i,k)  +wx(j,i,k-1)) - &
                 0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
    end do
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      w(j,i,2) = d_half * (wx(j,i,2)  +wx(j,i,1))
      w(j,i,kz) = d_half * (wx(j,i,kz)+wx(j,i,kzm1))
    end do
  end subroutine xtowstag

  subroutine xtoustag(ux,u)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: ux
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u
    integer(ik4) :: i , j , k

    do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
      u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                 0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
    end do
    if ( ma%has_bdyright ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
      end do
    end if

    if ( ma%has_bdyleft ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
      end do
    end if
  end subroutine xtoustag

  subroutine xtovstag(vx,v)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: vx
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: v
    integer(ik4) :: i , j , k

    do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
      v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                 0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
    end do
    if ( ma%has_bdytop ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
      end do
    end if
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
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

    do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
      u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                 0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
    end do
    if ( ma%has_bdyright ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi2,i,k) = d_half * (ux(jci2,i,k)+ux(jce2,i,k))
      end do
    end if
    if ( ma%has_bdyleft ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        u(jdi1,i,k) = d_half * (ux(jci1,i,k)+ux(jce1,i,k))
      end do
    end if

    ! Back to wind points: V (fourth order)

    do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
      v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                 0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
    end do
    if ( ma%has_bdytop ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi2,k) = d_half * (vx(j,ici2,k)+vx(j,ice2,k))
      end do
    end if
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jci1:jci2, k = 1:kz )
        v(j,idi1,k) = d_half * (vx(j,ici1,k)+vx(j,ice1,k))
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

    do concurrent ( j = jci1:jci2, i = ice1:ice2, k = 1:kz )
      ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                  0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
    end do
    if ( ma%has_bdyleft ) then
      do concurrent ( i = ice1:ice2, k = 1:kz )
        ux(jce1,i,k) = d_half * (u(jde1,i,k)+u(jdi1,i,k))
      end do
    end if
    if ( ma%has_bdyright ) then
      do concurrent ( i = ice1:ice2, k = 1:kz )
        ux(jce2,i,k) = d_half*(u(jde2,i,k) + u(jdi2,i,k))
      end do
    end if

    ! Compute V-wind on T points

    do concurrent ( j = jce1:jce2, i = ici1:ici2, k = 1:kz )
      vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                  0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
    end do
    if ( ma%has_bdybottom ) then
      do concurrent ( j = jce1:jce2, k = 1:kz )
        vx(j,ice1,k) = d_half * (v(j,ide1,k)+v(j,idi1,k))
      end do
    end if
    if ( ma%has_bdytop ) then
      do concurrent ( j = jce1:jce2, k = 1:kz )
        vx(j,ice2,k) = d_half*(v(j,ide2,k) + v(j,idi2,k))
      end do
    end if
  end subroutine uvstagtox

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
