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

module mod_cu_tiedtke
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_stdio
  use mod_constants
  use mod_cu_common
  use mod_cu_tables
  use mod_service
  use mod_runparams, only : sigma, iqc, iqv, iqi, dx, ichem, iaerosol, &
         iindirect, ipptls, ichcumtra, rcmtimer, icup, dtcum, io3
  use mod_runparams, only : k2_const, kfac_deep, kfac_shal
  use mod_mpmessage
  use mod_runparams, only : rcrit
  use mod_runparams, only : entrdd, entshalp, entrpen_lnd, entrpen_ocn
  use mod_runparams, only : detrpen_lnd, detrpen_ocn
  use mod_runparams, only : rprc_ocn, rprc_lnd, revap_lnd, revap_ocn
  use mod_runparams, only : rhebc_ocn, rhebc_lnd, rcuc_lnd, rcuc_ocn
  use mod_runparams, only : rcpec_ocn, rcpec_lnd
  use mod_runparams, only : lmfpen, lmfmid, lmfdd, lepcld, lmfdudv, &
          lmfscv, lmfuvdis, lmftrac, lmfsmooth, lmfwstar
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_cu_tiedtke, tiedtkedrv

  real(rkx), parameter :: rhow = 1000.0_rkx
  real(rkx), parameter :: rkap = 0.4_rkx
  real(rkx), parameter :: qsmax = 0.5_rkx
  real(rkx), parameter :: cwdrag = (3.0_rkx/8.0_rkx)*0.506_rkx/0.200_rkx

  real(rkx) :: rtau, dtc
  real(rkx) :: rmfcfl ! Massflux multiple of cfl stability criterium
  integer(ik4) :: nk350, nk060, nk950

  ! Mass flux solver for momentum (1) or no (0)
  real(rkx), parameter :: rmfsoluv = 1.0_rkx
  ! Mass flux solver for T and q (1) or no (0)
  real(rkx), parameter :: rmfsoltq = 1.0_rkx
  ! Mass flux solver for tracer (1) or no (0)
  real(rkx), parameter :: rmfsolct = 1.0_rkx

  ! Use CFL mass flux limit (1) or absolut limit (0)
  real(rkx), parameter :: rmflic = 1.0_rkx
  ! Value of absolut mass flux limit
  real(rkx), parameter :: rmflia = 0.0_rkx

  ! Relaxation time for melting of snow
  real(rkx), parameter :: rtaumel = 5.0_rkx*3600.0_rkx*1.5_rkx

  ! Updraught velocity perturbation for implicit (m/s)
  real(rkx), parameter :: ruvper = 0.3_rkx

  ! Maximum allowed cloud thickness for shallow cloud (Pa)
  real(rkx), parameter :: rdepths = 4.0e4_rkx

  ! Fractional massflux for downdrafts at lfs
  real(rkx), parameter :: rmfdeps = 0.3_rkx

  ! evaporation coefficient for kuo0
  integer(ik4), pointer, contiguous, dimension(:,:) :: ilab
  integer(ik4), pointer, contiguous, dimension(:) :: ktype
  logical, pointer, contiguous, dimension(:) :: ldland
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pxtm1, pxtte
  real(rkx), pointer, contiguous, dimension(:,:) :: ptm1, pqm1, pum1, pvm1,    &
        pxlm1, pxim1, pxite, papp1, paphp1, pxtec, pqtec, zlude, &
        pmflxr, pccn, zcvrout
  real(rkx), pointer, contiguous, dimension(:) :: prsfc, pssfc, &
        ptopmax, xphfx, xpqfx
  real(rkx), pointer, contiguous, dimension(:,:) :: ptte, pvom, pvol, pqte, &
        pxlte, pverv, pmfu, xpg, xpgh
  integer(ik4), pointer, contiguous, dimension(:) :: kctop, kcbot
  integer(ik4) :: nipoi
  integer(ik4), pointer, contiguous, dimension(:) :: imap, jmap
  integer(ik4) :: nskmax

  integer(ik4), public :: nmctop
  real(rkx), pointer, contiguous, public, dimension(:) :: pmean

  ! Max massflux value
  real(rkx), parameter :: cmfcmax = 1.0_rkx
  ! Min massflux value (for safety)
  real(rkx), parameter :: cmfcmin = 1.0e-10_rkx
  ! Rainfall max elevation
  real(rkx), parameter :: zdlev   = 1.5e4_rkx
  ! Relat. cloud massflux at level above nonbuoyancy
  real(rkx), parameter :: cmfctop = 0.35_rkx

  ! Fractional massflux for downdrafts at lfs
  real(rkx), parameter :: cmfdeps = 0.3_rkx

  ! Trigger parameter for the convection 0.0 <-> 1.0
  real(rkx), parameter :: ctrigger = -1.1_rkx

  real(rkx), parameter :: almostzero = dlowval

  contains
  !
  ! This subroutines allocates space
  !
  subroutine allocate_mod_cu_tiedtke
    implicit none
    integer(ik4) :: i, j, ii
    call getmem(pmean,1,kz,'mod_cu_tiedtke:pmean')
    nipoi = 0
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cuscheme(j,i) == 5 ) then
          nipoi = nipoi + 1
        end if
      end do
    end do
    if ( nipoi == 0 ) then
      return
    end if
    call getmem(imap,1,nipoi,'mod_cu_tiedtke:imap')
    call getmem(jmap,1,nipoi,'mod_cu_tiedtke:jmap')
    ii = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cuscheme(j,i) == 5 ) then
          imap(ii) = i
          jmap(ii) = j
          ii = ii + 1
        end if
      end do
    end do
    call getmem(ptte,1,nipoi,1,kz,'mod_cu_tiedtke:ptte')
    call getmem(pvom,1,nipoi,1,kz,'mod_cu_tiedtke:pvom')
    call getmem(pvol,1,nipoi,1,kz,'mod_cu_tiedtke:pvol')
    call getmem(pqte,1,nipoi,1,kz,'mod_cu_tiedtke:pqte')
    call getmem(pxlte,1,nipoi,1,kz,'mod_cu_tiedtke:pxlte')
    call getmem(pverv,1,nipoi,1,kz,'mod_cu_tiedtke:pverv')
    call getmem(pmfu,1,nipoi,1,kz,'mod_cu_tiedtke:pmfu')
    call getmem(xpg,1,nipoi,1,kz,'mod_cu_tiedtke:xpg')
    call getmem(xpgh,1,nipoi,1,kz+1,'mod_cu_tiedtke:xpgh')
    call getmem(pxtm1,1,nipoi,1,kz,1,ntr,'mod_cu_tiedtke:pxtm1')
    call getmem(pxtte,1,nipoi,1,kz,1,ntr,'mod_cu_tiedtke:pxtte')
    call getmem(ilab,1,nipoi,1,kz,'mod_cu_tiedtke:ilab')
    call getmem(ktype,1,nipoi,'mod_cu_tiedtke:ktype')
    call getmem(ptm1,1,nipoi,1,kz,'mod_cu_tiedtke:ptm1')
    call getmem(pqm1,1,nipoi,1,kz,'mod_cu_tiedtke:pqm1')
    call getmem(pum1,1,nipoi,1,kz,'mod_cu_tiedtke:pum1')
    call getmem(pvm1,1,nipoi,1,kz,'mod_cu_tiedtke:pvm1')
    call getmem(pxlm1,1,nipoi,1,kz,'mod_cu_tiedtke:pxlm1')
    call getmem(pxim1,1,nipoi,1,kz,'mod_cu_tiedtke:pxim1')
    call getmem(pxite,1,nipoi,1,kz,'mod_cu_tiedtke:pxite')
    call getmem(papp1,1,nipoi,1,kz,'mod_cu_tiedtke:papp1')
    call getmem(zlude,1,nipoi,1,kz,'mod_cu_tiedtke:zlude')
    call getmem(pxtec,1,nipoi,1,kz,'mod_cu_tiedtke:pxtec')
    call getmem(pqtec,1,nipoi,1,kz,'mod_cu_tiedtke:pqtec')
    call getmem(pccn,1,nipoi,1,kz,'mod_cu_tiedtke:pccn')
    call getmem(zcvrout,1,nipoi,1,kz,'mod_cu_tiedtke:zcvrout')
    call getmem(pmflxr,1,nipoi,1,kz+1,'mod_cu_tiedtke:pmflxr')
    call getmem(kctop,1,nipoi,'mod_cu_tiedtke:kctop')
    call getmem(kcbot,1,nipoi,'mod_cu_tiedtke:kcbot')
    call getmem(paphp1,1,nipoi,1,kzp1,'mod_cu_tiedtke:paphp1')
    call getmem(prsfc,1,nipoi,'mod_cu_tiedtke:prsfc')
    call getmem(pssfc,1,nipoi,'mod_cu_tiedtke:pssfc')
    call getmem(ptopmax,1,nipoi,'mod_cu_tiedtke:ptopmax')
    call getmem(xphfx,1,nipoi,'mod_cu_tiedtke:xphfx')
    call getmem(xpqfx,1,nipoi,'mod_cu_tiedtke:xpqfx')
    call getmem(ldland,1,nipoi,'mod_cu_tiedtke:ldland')

  end subroutine allocate_mod_cu_tiedtke
  !
  ! This subroutines calls cucall
  !
  subroutine tiedtkedrv(m2c)
    implicit none
    type(mod_2_cum), intent(in) :: m2c
    integer(ik4) :: i, j, k, n, ii
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tiedtkedrv'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nipoi == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    dtc = dtcum
    do concurrent( ii = 1:nipoi, k = 1:kz )
      ilab(ii,k) = 2
    end do

    if ( ichem == 1 ) then
      do n = 1, ntr
        ! prevent excessive convectivedownward transport for o3
        if ( n == io3 ) cycle
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          i = imap(ii)
          j = jmap(ii)
          ! tracers input profile : implicit loop on tracer
          pxtm1(ii,k,n) = m2c%chias(j,i,k,n)
        end do
      end do
      if ( idynamic == 3 ) then
        do concurrent ( ii = 1:nipoi, k = 1:kz, n = 1:ntr )
          i = imap(ii)
          j = jmap(ii)
          pxtte(ii,k,n) = m2c%chiten(j,i,k,n)
        end do
      else
        do concurrent ( ii = 1:nipoi, k = 1:kz, n = 1:ntr )
          i = imap(ii)
          j = jmap(ii)
          pxtte(ii,k,n) = m2c%chiten(j,i,k,n)/m2c%psb(j,i)
        end do
      end if
      if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          i = imap(ii)
          j = jmap(ii)
          pccn(ii,k) = m2c%ccn(j,i,k)
        end do
      end if
    else
      do concurrent ( ii = 1:nipoi, k = 1:kz, n = 1:ntr )
        pxtm1(ii,k,n) = d_zero ! tracers input profiles
        pxtte(ii,k,n) = d_zero ! tracer tendencies
      end do
    end if
    do concurrent ( ii = 1:nipoi )
      ! AMT NOTE: This is used in the switch between deep and shallow
      ! convection. The simpler switch on pressure difference still
      ! used in ECMWF is commented out - possibly tests should be made
      ! to reinstate the ECMWF version of this
      ! this array will then be obsolete
      i = imap(ii)
      j = jmap(ii)
      xpqfx(ii) = m2c%qfx(j,i)
      xphfx(ii) = m2c%hfx(j,i)
      ldland(ii) = (m2c%ldmsk(j,i) == 1)
    end do

    do concurrent ( ii = 1:nipoi, k = 1:kz )
      i = imap(ii)
      j = jmap(ii)
      papp1(ii,k) = m2c%pas(j,i,k)       ! Pressure in Pa
      xpg(ii,k)   = m2c%zas(j,i,k)*egrav ! geopotential
      ptm1(ii,k)  = m2c%tas(j,i,k)  ! temperature
      pum1(ii,k)  = m2c%uas(j,i,k)  ! u (guessing!)
      pvm1(ii,k)  = m2c%vas(j,i,k)  ! v     "
      pverv(ii,k) = avg_ww(j,i,k)
      pqm1(ii,k)  = m2c%qxas(j,i,k,iqv) ! humidity
      pxlm1(ii,k) = m2c%qxas(j,i,k,iqc) ! cloud liquid water
      pvom(ii,k)  = m2c%uten(j,i,k)
      pvol(ii,k)  = m2c%vten(j,i,k)
      ! scheme diagnostic output - tendencies due to convection
      pxtec(ii,k) = d_zero ! detrained cloud water tendency
      pqtec(ii,k) = d_zero ! detrained humidity tendency
      ! pressure top limit for convection: the level above tropopause
    end do
    do concurrent ( ii = 1:nipoi, k = 1:kz )
      i = imap(ii)
      j = jmap(ii)
      if ( idynamic == 3 ) then
        ptte(ii,k)  = m2c%tten(j,i,k)
        pqte(ii,k)  = m2c%qxten(j,i,k,iqv)
        pxlte(ii,k) = m2c%qxten(j,i,k,iqc)
        if ( ipptls > 1 ) then
          pxim1(ii,k) = m2c%qxas(j,i,k,iqi)      ! cloud ice water
          pxite(ii,k) = m2c%qxten(j,i,k,iqi)
        else
          pxim1(ii,k) = d_zero
          pxite(ii,k) = d_zero
        end if
      else
        ptte(ii,k)  = m2c%tten(j,i,k)/m2c%psb(j,i)
        pqte(ii,k)  = m2c%qxten(j,i,k,iqv)/m2c%psb(j,i)
        pxlte(ii,k) = m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
        if ( ipptls > 1 ) then
          pxim1(ii,k) = m2c%qxas(j,i,k,iqi)      ! cloud ice water
          pxite(ii,k) = m2c%qxten(j,i,k,iqi)/m2c%psb(j,i)
        else
          pxim1(ii,k) = d_zero
          pxite(ii,k) = d_zero
        end if
      end if
    end do

    ! Transform to specific humidity
    do concurrent ( ii = 1:nipoi, k = 1:kz )
      pqte(ii,k) = pqte(ii,k)/(d_one + pqm1(ii,k))**2
      pqm1(ii,k) = pqm1(ii,k)/(d_one + pqm1(ii,k))
    end do

    nskmax = nint(17747.5/(dx*d_r1000))
    call setup(nskmax,pmean)

    do concurrent ( ii = 1:nipoi, k = 1:kzp1 )
      i = imap(ii)
      j = jmap(ii)
      ! 1st guess pressure at full levels
      paphp1(ii,k) = m2c%pasf(j,i,k)
      xpgh(ii,k) = m2c%zfs(j,i,k)*egrav ! geopotential
    end do

    ! Output variables (1d)
    do concurrent ( ii = 1:nipoi )
      prsfc(ii) = d_zero ! CHECK - surface rain flux
      pssfc(ii) = d_zero ! CHECK - surface snow flux
      ktype(ii) = 0
      kctop(ii) = 0
      kcbot(ii) = 0
    end do
    do concurrent ( ii = 1:nipoi, k = 1:kz )
      pmfu(ii,k) = d_zero
    end do

    call cucall(nipoi,nipoi,kz,kzp1,kzm1,ilab,ntr,pxtm1,pxtte,ptm1,   &
                pqm1,pum1,pvm1,pxlm1,pxim1,ptte,pqte,pvom,pvol,pxlte, &
                pxite,pverv,pxtec,pqtec,xphfx,xpqfx,papp1,paphp1,xpg, &
                xpgh,prsfc,pssfc,pmflxr,pmfu,zlude,ktype,ldland,kctop,&
                kcbot,ptopmax,pccn,zcvrout)
    !
    ! postprocess some fields including precipitation fluxes
    !
    !   KTYPE = 1 => DEEP CONVECTION
    !   KTYPE = 2 => SHALLOW CONVECTION
    !   KTYPE = 3 => MIDLEVEL CONVECTION
    !
    if ( ipptls > 1 ) then
      do concurrent ( ii = 1:nipoi )
        if (ktype(ii) > 0) then
          i = imap(ii)
          j = jmap(ii)
          total_precip_points = total_precip_points + 1
          ! rainfall for surface
          cu_prate(j,i) = cu_prate(j,i) + prsfc(ii) + pssfc(ii)
          cu_srate(j,i) = cu_srate(j,i) + pssfc(ii)
        end if
      end do
    else
      do concurrent ( ii = 1:nipoi )
        if (ktype(ii) > 0) then
          i = imap(ii)
          j = jmap(ii)
          total_precip_points = total_precip_points + 1
          ! rainfall for surface
          cu_prate(j,i)= cu_prate(j,i) + prsfc(ii) + pssfc(ii)
        end if
      end do
    end if
    !
    ! update tendencies - note that rate were ADDED in cudtdq
    !                     thus here we must reset the rates.
    if ( idynamic == 3 ) then
      do concurrent ( ii = 1:nipoi, k = 1:kz )
        if ( ktype(ii) > 0 ) then
          i = imap(ii)
          j = jmap(ii)
          cu_tten(j,i,k) = ptte(ii,k) - m2c%tten(j,i,k)
          ! Tendency in specific humidity to mixing ratio tendency.
          cu_qten(j,i,k,iqv) = pqte(ii,k)/(d_one-pqm1(ii,k))**2 - &
                               m2c%qxten(j,i,k,iqv)
        end if
      end do
    else
      do concurrent ( ii = 1:nipoi, k = 1:kz )
        if ( ktype(ii) > 0 ) then
          i = imap(ii)
          j = jmap(ii)
          cu_tten(j,i,k) = ptte(ii,k) - m2c%tten(j,i,k)/m2c%psb(j,i)
          ! Tendency in specific humidity to mixing ratio tendency.
          cu_qten(j,i,k,iqv) = pqte(ii,k)/(d_one-pqm1(ii,k))**2 - &
                               m2c%qxten(j,i,k,iqv)/m2c%psb(j,i)
        end if
      end do
    end if
    do concurrent ( ii = 1:nipoi, k = 1:kz )
      if ( ktype(ii) > 0 ) then
        i = imap(ii)
        j = jmap(ii)
        cu_uten(j,i,k) = pvom(ii,k) - m2c%uten(j,i,k)
        cu_vten(j,i,k) = pvol(ii,k) - m2c%vten(j,i,k)
        cu_qdetr(j,i,k) = zlude(ii,k)*egrav/(paphp1(ii,k+1)-paphp1(ii,k))
        cu_raincc(j,i,k) = pmflxr(ii,k)
      end if
    end do
    if ( ipptls > 1 ) then
      if ( idynamic == 3 ) then
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_qten(j,i,k,iqc) = pxlte(ii,k) - m2c%qxten(j,i,k,iqc)
            cu_qten(j,i,k,iqi) = pxite(ii,k) - m2c%qxten(j,i,k,iqi)
          end if
        end do
      else
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_qten(j,i,k,iqc) = pxlte(ii,k) - &
                m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
            cu_qten(j,i,k,iqi) = pxite(ii,k) - &
                m2c%qxten(j,i,k,iqi)/m2c%psb(j,i)
          end if
        end do
      end if
    else
      if ( idynamic == 3 ) then
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_qten(j,i,k,iqc) = (pxlte(ii,k)+pxite(ii,k)) - &
                m2c%qxten(j,i,k,iqc)
            cu_tten(j,i,k) = cu_tten(j,i,k) + pxite(ii,k)*wlhf/cpd
          end if
        end do
      else
        do concurrent ( ii = 1:nipoi, k = 1:kz )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_qten(j,i,k,iqc) = (pxlte(ii,k) + pxite(ii,k)) - &
                m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
            cu_tten(j,i,k) = cu_tten(j,i,k) + pxite(ii,k)*wlhf/cpd
          end if
        end do
      end if
    end if

    if ( ichem == 1 .and. ichcumtra == 1 .and. &
         .not. any(icup == 2 .or. icup == 6) ) then
      if ( idynamic == 3 ) then
        do concurrent ( ii = 1:nipoi, k = 1:kz, n = 1:ntr )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_chiten(j,i,k,n) = pxtte(ii,k,n) - m2c%chiten(j,i,k,n)
          end if
        end do
      else
        do concurrent ( ii = 1:nipoi, k = 1:kz, n = 1:ntr )
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_chiten(j,i,k,n) = pxtte(ii,k,n) - &
                        m2c%chiten(j,i,k,n)/m2c%psb(j,i)
          end if
        end do
      end if
      ! build for chemistry 3d table of conv. precipitation
      ! tend. Kg/kg/s
      do concurrent ( ii = 1:nipoi, k = 1:kz )
        if (ktype(ii) > 0) then
          i = imap(ii)
          j = jmap(ii)
          cu_convpr(j,i,k) = zcvrout(ii,k)
        end if
      end do
    end if

    !
    ! Xu, K.-M., and S. K. Krueger:
    !   Evaluation of cloudiness parameterizations using a cumulus
    !   ensemble model, Mon. Wea. Rev., 119, 342-367, 1991.
    !
    do concurrent ( ii = 1:nipoi )
      if (ktype(ii) > 0) then
        i = imap(ii)
        j = jmap(ii)
        cu_ktop(j,i) = kctop(ii)
        cu_kbot(j,i) = kcbot(ii)
        !$acc loop seq
        do k = kctop(ii), kcbot(ii)
          if ( ktype(ii) == 1 ) then
            cu_cldfrc(j,i,k) = kfac_deep*log(d_one+k2_const*pmfu(ii,k))
            cu_cldfrc(j,i,k) = min(max(0.0_rkx,cu_cldfrc(j,i,k)),0.6_rkx)
          else if ( ktype(ii) == 2 ) then
            cu_cldfrc(j,i,k) = kfac_shal*log(d_one+k2_const*pmfu(ii,k))
            cu_cldfrc(j,i,k) = min(max(0.0_rkx,cu_cldfrc(j,i,k)),0.2_rkx)
          else
            cu_cldfrc(j,i,k) = d_half*(kfac_deep+kfac_shal) * &
                           log(d_one+k2_const*pmfu(ii,k))
            cu_cldfrc(j,i,k) = min(max(0.0_rkx,cu_cldfrc(j,i,k)),0.4_rkx)
          end if
        end do
      end if
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains
    !
    ! This routine defines parameters for massflux scheme
    !
    subroutine setup(smax,pmean)
      implicit none
      integer(ik4), intent(in) :: smax
      real(rkx), dimension(kz), intent(in) :: pmean
      integer(ik4) :: k
      rtau = d_one+264.0_rkx/real(smax,rkx)
      !rtau = d_one+528.0_rkx/real(smax,rkx)
      rtau = min(3.0_rkx,rtau)
      if ( smax >= 511 ) then
        rmfcfl = 3.0_rkx
      else
        rmfcfl = 5.0_rkx
      end if
      nk350 = kz
      nk060 = 1
      nk950 = kz
      !$acc update host(pmean)
      do k = kz, 2, -1
        if ( pmean(k)/pmean(kz)*stdp >= 350.e2_rkx ) nk350 = k
        if ( pmean(k)/pmean(kz)*stdp >=  60.e2_rkx ) nk060 = k
        if ( pmean(k)/pmean(kz)*stdp >= 950.e2_rkx ) nk950 = k
      end do
      !$acc update device(pmean)
    end subroutine setup

  end subroutine tiedtkedrv
  !
  ! - MASTER ROUTINE - PROVIDES INTERFACE FOR:
  !        *CUMASTRX* (CUMULUS PARAMETERIZATION)
  !
  !  PROVIDES INPUT FOR CUMASTR
  !  RECEIVES UPDATED TENDENCIES, PRECIPITATION.
  !
  subroutine cucall(kproma,kbdim,klev,klevp1,klevm1,ilab,ktrac,  &
                    pxtm1,pxtte,ptm1,pqm1,pum1,pvm1,pxlm1,pxim1, &
                    ptte,pqte,pvom,pvol,pxlte,pxite,pverv,pxtec, &
                    pqtec,pshfla,pqhfla,papp1,paphp1,pgeo,pgeoh, &
                    prsfc,pssfc,pmflxr,pmfu,zlude,ktype,ldland,  &
                    kctop,kcbot,ptopmax,pccn,zcvrout)
    implicit none
    integer(ik4), intent(in) :: kbdim, klev, klevm1, klevp1, kproma, ktrac
    integer(ik4), dimension(kbdim,klev), intent(inout) :: ilab
    integer(ik4), dimension(kbdim), intent(inout) :: ktype
    logical, dimension(kbdim), intent(inout) :: ldland
    real(rkx), dimension(kbdim,klevp1), intent(in) :: paphp1, pgeoh
    real(rkx), dimension(kbdim,klevp1), intent(out) :: pmflxr
    real(rkx), dimension(kbdim,klev), intent(inout) :: ptte, pxite, &
           pxlte, pmfu, pqte, pvom, pvol
    real(rkx), dimension(kbdim,klev), intent(in) :: papp1, pgeo,   &
           pqm1, pqtec, ptm1, pum1, pverv, pvm1, pxim1, pxlm1, pxtec, pccn
    real(rkx), dimension(kbdim,klev), intent(out) :: zlude
    real(rkx), dimension(kbdim), intent(in) :: pqhfla, pshfla
    real(rkx), dimension(kbdim), intent(inout) :: prsfc, pssfc
    real(rkx), dimension(kbdim), intent(inout) :: ptopmax
    integer(ik4), dimension(kbdim), intent(out) :: kctop, kcbot
    real(rkx), dimension(kbdim,klev,ktrac), intent(in) :: pxtm1
    real(rkx), dimension(kbdim,klev,ktrac), intent(inout) :: pxtte
    integer(ik4), dimension(kbdim) :: itopec2
    integer(ik4) :: ilevmin, jk, jl, jt
    logical, dimension(kbdim) :: locum
    real(rkx), dimension(kbdim,klev) :: zlu, zmfd, zqp1, zqsat, zqu, &
      zqude, ztp1, ztu, ztvp1, zup1, zvp1, zxp1, zcvrout
    real(rkx), dimension(kbdim) :: zrain, ztopmax
    integer(ik4), dimension(kbdim) :: kbotsc
    logical, dimension(kbdim) :: ldsc
    real(rkx), dimension(kbdim,klev,ktrac) :: zxtp1, zxtu
    real(rkx), dimension(kbdim,klev+1) :: pmflxs
    real(rkx), dimension(kbdim,klev) :: pmfude_rate, pmfdde_rate
    real(rkx), dimension(kbdim) :: pcape
    real(rkx), dimension(kbdim,klev+1) :: pqhfl, pahfs
    integer(ik4) :: it
    real(rkx) :: zxlp, zxip

    lookupoverflow = .false.
    !
    !---------------------------------------
    ! 1. CALCULATE T,Q AND QS AT MAIN LEVELS
    ! --------------------------------------
    !
    do concurrent ( jl = 1:kproma, jk = 1:klev )
      ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*dtc
      zqp1(jl,jk) = max(1.0e-8_rkx,pqm1(jl,jk) + pqte(jl,jk)*dtc)
      zxlp = max(d_zero,pxlm1(jl,jk) + pxlte(jl,jk)*dtc)
      zxip = max(d_zero,pxim1(jl,jk) + pxite(jl,jk)*dtc)
      zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*dtc
      zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*dtc
      zxp1(jl,jk) = max(d_zero,zxlp+zxip)
      it = int(ztp1(jl,jk)*d_1000)
#ifdef DEBUG
      if ( it < jptlucu1 .or. it > jptlucu2 ) then
        lookupoverflow = .true.
      end if
#endif
      it = max(min(it,jptlucu2),jptlucu1)
      zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
      zqsat(jl,jk) = min(qsmax,zqsat(jl,jk))
      zqsat(jl,jk) = zqsat(jl,jk)/(d_one-ep1*zqsat(jl,jk))
      ztvp1(jl,jk) = ztp1(jl,jk)*d_one+ep1*(zqp1(jl,jk)-zxp1(jl,jk))
    end do
    if ( lookupoverflow ) then
      do jk = 1, klev
        do jl = 1, kproma
          it = int(ztp1(jl,jk)*d_1000)
          if ( it < jptlucu1 .or. it > jptlucu2 ) then
            write(stderr,'(a,f12.4,a,i8)') &
              '! LOOKUP PROBLEM FOR T = ',real(ztp1(jl,jk)), &
              ' at ', rcmtimer%str( )
          end if
        end do
      end do
      call fatal(__FILE__,__LINE__, &
                 'Cumulus Tables lookup error: OVERFLOW')
    end if

    do concurrent ( jl = 1:kproma, jk = 1:klev, jt = 1:ktrac)
      zxtp1(jl,jk,jt) = max(d_zero,pxtm1(jl,jk,jt) + pxtte(jl,jk,jt)*dtc)
    end do

    do concurrent ( jl = 1:kproma )
      zrain(jl) = d_zero
      locum(jl) = .false.
    end do
    !
    !---------------------------------------------------------------
    ! 2. CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
    ! --------------------------------------------------------------
    !
    ! =====================
    ! ABOUT pqhfl and pahfs
    ! =====================
    ! The variables in question are defined in the turbulence scheme and
    ! are the turbulence flux of water vapour and sensible heat as expected...
    ! The naming changes in the turbulence scheme; they are PDIFTQ and PDIFTS,
    ! respectively...
    ! The flux values are just backed out from the new-old tendency values
    ! differences over the scheme (because of course the scheme employs an
    ! implicit solver.)
    !
    ! These are the definitions at the start of vdfmain.f90
    !
    ! *PDIFTQ*       TURBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)
    ! *PDIFTS*       TURBULENT FLUX OF HEAT                        J/(M2*S)
    !
    ! And the code where they are defined is here:
    !
    !DO JK=KLEV-1,1,-1
    !  DO JL=KIDIA,KFDIA
    !   ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
    !   !...changes in q,l,i tendencies are converted to fluxes
    !   PDIFTQ(JL,JK) = (PQE (JL,JK+1) - ZQEA(JL,JK+1)) * ZGDPH +
    !         PDIFTQ(JL,JK+1)
    !   PDIFTL(JL,JK) = (PLE (JL,JK+1) - ZLEA(JL,JK+1)) * ZGDPH +
    !         PDIFTL(JL,JK+1)
    !   PDIFTI(JL,JK) = (PIE (JL,JK+1) - ZIEA(JL,JK+1)) * ZGDPH +
    !         PDIFTI(JL,JK+1)
    !   !...slg=s-Lc*ql-Ld*qi (same for fluxes)
    !   PDIFTS(JL,JK) = ZDIFTSLG(JL,JK) &
    !         & + RLVTT * PDIFTL(JL,JK) + RLSTT * PDIFTI(JL,JK)
    !  ENDDO
    !ENDDO
    ! ZGDPH - is just the mass of the layers (g/dp) defined at half level
    ! PQE and ZQEA are the new and old tendencies
    ! *PQE*          MOISTURE TENDENCY KG/(KG S)
    ! The sensible heat calculation is slightly more complicated as it is
    ! backed out of the conserved variables used in the turbulence scheme.
    ! I presume REGCM's turbulence scheme already or can easily provide these
    ! quantities...
    !
    do concurrent ( jl = 1:kproma )
      do jk = kz, 1, -1
        pqhfl(jl,jk) = xpqfx(jl) * exp(-egrav*(d_one - sigma(jk)))
        pahfs(jl,jk) = xphfx(jl) * exp(-egrav*(d_one - sigma(jk)))
      end do
      pqhfl(jl,kzp1) = xpqfx(jl)
      pahfs(jl,kzp1) = xphfx(jl)
    end do
    !
    ! Reset to zero output precipitation fluxes
    !
    !$acc kernels
    pmflxr = d_zero
    pmflxs = d_zero
    !$acc end kernels
    !
    ! Call scheme
    !
    call ntiedtke(1,kproma,kbdim,klev,ldland,ztp1,zqp1,          &
                  zup1,zvp1,zxp1,pverv,pqhfl,pahfs,papp1,paphp1, &
                  pgeo,pgeoh,ptte,pqte,pvom,pvol,pxlte,pxite,    &
                  locum,ktype,kcbot,kctop,kbotsc,ldsc,ztu,zqu,   &
                  zlu,pmflxr,pmflxs,zrain,pmfu,zmfd,zlude,       &
                  pmfude_rate,pmfdde_rate,pcape,ktrac,pxtm1,     &
                  pxtte,pccn,zcvrout)
    !$acc kernels
    prsfc = pmflxr(:,klev+1)*1.e3_rkx
    pssfc = pmflxs(:,klev+1)*1.e3_rkx
    !$acc end kernels
    !
    ! ----------------------------------------------
    ! 3. PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
    ! ----------------------------------------------
    !
    ilevmin = klev - 4
    do concurrent ( jl = 1:kproma )
      itopec2(jl) = klevp1
    end do
    do concurrent ( jl = 1:kproma )
      !$acc loop seq
      do jk = 1, ilevmin
        if ( ilab(jl,jk) == 2 .and. itopec2(jl) == klevp1 ) itopec2(jl) = jk
      end do
    end do
    !$acc kernels
    ztopmax(1:kproma) = ptopmax(1:kproma)
    !$acc end kernels
    do concurrent ( jl = 1:kproma )
      if ( itopec2(jl) == 1 ) then
        ptopmax(jl) = papp1(jl,1)
      else if ( itopec2(jl) /= klevp1 ) then
        ptopmax(jl) = paphp1(jl,itopec2(jl))
      else
        ptopmax(jl) = 99999.0_rkx
      end if
      ptopmax(jl) = min(ptopmax(jl),ztopmax(jl))
    end do
  end subroutine cucall
  !
  ! This routine computes the physical tendencies of the prognostic
  ! variables t,q,u,v and tracers due to convective processes.
  ! Processes considered are: convective fluxes, formation of precipitation,
  ! evaporation of falling rain below cloud base, saturated cumulus downdrafts.
  !
  ! Parameterization is done using a massflux-scheme.
  !
  !  (1) Define constants and parameters
  !  (2) Specify values (t,q,qs...) at half levels and
  !      initialize updraft- and downdraft-values
  !  (3) Calculate cloud base and specify cloud base massflux from
  !      PBL moisture budget
  !  (4) Do cloud ascent in absence of downdrafts
  !  (5) Do downdraft calculations:
  !      (A) Determine values at lfs
  !      (B) Determine moist descent
  !      (C) Recalculate cloud base massflux considering the
  !          effect of downdrafts
  !  (6) Do final cloud ascent
  !  (7) Do final adjusments to convective fluxes
  !      Do evaporation in subcloud layer
  !  (8) Calculate increments of t and q
  !  (9) Calculate increments of u and v
  !
  subroutine ntiedtke(n1,n2,np,nk,ldland,                           &
                      t,q,u,v,qctot,omega,qhfl,ahfs,ph,pf,geo,geof, &
                      tent,tenq,tenu,tenv,tenl,teni,                &
                      ldcum,ktype,kcbot,kctop,kbotsc,ldsc,          &
                      tu,qu,lu,mflxr,mflxs,rain,mfu,mfd,            &
                      lude,mfude_rate,mfdde_rate,cape,              &
                      ntrac,qtrac,tenc,ccn,cvrainout)
    implicit none
    integer(ik4), intent(in) :: np    ! number of points
    integer(ik4), intent(in) :: nk    ! number of levels
    integer(ik4), intent(in) :: n1    ! start point
    integer(ik4), intent(in) :: n2    ! end point
    integer(ik4), intent(in) :: ntrac ! number of tracers
    ! Land/sea mask flag
    logical, dimension(np), intent(in) :: ldland
    !-----------------------------------------------------
    ! ::::::::::: Input atmospheric condition ::::::::::::
    !-----------------------------------------------------
    ! Temperature K
    real(rkx), dimension(np,nk), intent(in) :: t
    ! Water Vapor Specific Humidity kg/kg
    real(rkx), dimension(np,nk), intent(in) :: q
    ! U wind m/s
    real(rkx), dimension(np,nk), intent(in) :: u
    ! V wind m/s
    real(rkx), dimension(np,nk), intent(in) :: v
    ! Total (Liquid + Ice) Mixing ratio kg/kg
    real(rkx), dimension(np,nk), intent(in) :: qctot
    ! Vertical Pressure Pa/s
    real(rkx), dimension(np,nk), intent(in) :: omega
    ! Turbolent moisture flux kg/(m^2*s)
    real(rkx), dimension(np,nk+1), intent(in) :: qhfl
    ! Turbolent sensible heat flux W/(m^2*s)
    real(rkx), dimension(np,nk+1), intent(in) :: ahfs
    ! Pressure Pa
    real(rkx), dimension(np,nk), intent(in) :: ph
    ! Pressure on sigma levels Pa
    real(rkx), dimension(np,nk+1), intent(in) :: pf
    ! Geopotential m^2/s^2
    real(rkx), dimension(np,nk), intent(in) :: geo
    ! Geopotential m^2/s^2 on sigma levels
    real(rkx), dimension(np,nk+1), intent(in) :: geof
    ! Tracer Mixing ratio kg/kg
    real(rkx), dimension(np,nk,ntrac), intent(in) :: qtrac
    ! Cloud condensation nuclei #
    real(rkx), dimension(np,nk), intent(in) :: ccn
    !-----------------------------------------------------
    ! ::::::::::::::: Output Tendencies ::::::::::::::::::
    !-----------------------------------------------------
    ! Temperature  tendency K/s
    real(rkx), dimension(np,nk), intent(inout) :: tent
    ! specific humidity tendency kg/(kg*s)
    real(rkx), dimension(np,nk), intent(inout) :: tenq
    ! U wind m/s^2
    real(rkx), dimension(np,nk), intent(inout) :: tenu
    ! V wind m/s^2
    real(rkx), dimension(np,nk), intent(inout) :: tenv
    ! Tracer mixing ratio tendency kg/(kg*s)
    real(rkx), dimension(np,nk,ntrac), intent(inout) :: tenc
    ! Cloud liquid water mixing ratio tendency kg/(kg*s)
    real(rkx), dimension(np,nk), intent(out) :: tenl
    ! Ice water mixing ratio tendency kg/(kg*s)
    real(rkx), dimension(np,nk), intent(out) :: teni
    ! Detrained liquid water kg/(m^2*s)
    real(rkx), dimension(np,nk), intent(out) :: lude
    ! Cumulus activation flag
    logical, dimension(np), intent(out) :: ldcum
    ! Type of cumulus convection
    integer(ik4), dimension(np), intent(out) :: ktype
    ! Bottom level of cumulus cloud
    integer(ik4), dimension(np), intent(out) :: kcbot
    ! Top level of cumulus cloud
    integer(ik4), dimension(np), intent(out) :: kctop
    ! Bottom level of shallow cumulus cloud
    integer(ik4), dimension(np), intent(out) :: kbotsc
    ! Shallow cumulus flag
    logical, dimension(np), intent(out) :: ldsc
    ! Updraft temperature K
    real(rkx), dimension(np,nk), intent(out) :: tu
    ! Updraft Water Vapor specific humidity kg/kg
    real(rkx), dimension(np,nk), intent(out) :: qu
    ! Updraft Cloud Liquid Water mixing ratio kg/kg
    real(rkx), dimension(np,nk), intent(out) :: lu
    ! rainfall Production in updraft
    real(rkx), dimension(np,nk), intent(out) :: cvrainout
    ! Rain Mass Flux kg/(m^2*s)
    real(rkx), dimension(np,nk+1), intent(out) :: mflxr
    ! Snow Mass Flux kg/(m^2*s)
    real(rkx), dimension(np,nk+1), intent(out) :: mflxs
    ! Total precipitation produced in convective updrafts
    real(rkx), dimension(np), intent(out) :: rain
    ! Massflux updrafts kg/(m^2*s)
    real(rkx), dimension(np,nk), intent(out) :: mfu
    ! Massflux downdrafts kg/(m^2*s)
    real(rkx), dimension(np,nk), intent(out) :: mfd
    ! Massflux updrafts rate kg/(m^3*s)
    real(rkx), dimension(np,nk), intent(out) :: mfude_rate
    ! Massflux downdrafts rate kg/(m^3*s)
    real(rkx), dimension(np,nk), intent(out) :: mfdde_rate
    ! Convectve available potential energy J/kg
    real(rkx), dimension(np), intent(out) :: cape
    !-----------------------------------------------------
    ! ::::::::::::::::: Local Variables ::::::::::::::::::
    !-----------------------------------------------------
    ! Vertically averaged updraught velocity m/s
    real(rkx), dimension(np) :: wmean
    ! Saturation mixing ratio kg/kg
    real(rkx), dimension(np,nk) :: qs
    ! Temperature on full sigma levels (where massfluxes are computed) K
    real(rkx), dimension(np,nk) :: tf
    ! Mixing ratio on full sigma levels kg/kg
    real(rkx), dimension(np,nk) :: qf
    ! Saturation mixing ratio on full sigma levels kg/kg
    real(rkx), dimension(np,nk) :: qsf
    ! Temperature in Downdrafts K
    real(rkx), dimension(np,nk) :: td
    ! Mixing ratio in downdrafts kg/kg
    real(rkx), dimension(np,nk) :: qd
    ! Flux of dry static energy in updrafts J/(m^2*s)
    real(rkx), dimension(np,nk) :: mfus
    ! Flux of dry static energy in downdrafts J/(m^2*s)
    real(rkx), dimension(np,nk) :: mfds
    ! Flux of humidity in updrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: mfuq
    ! Flux of humidity in downdrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: mfdq
    ! Flux difference of precipitation in updrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: dmfup
    ! Flux difference of precipitation in downdrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: dmfdp
    ! Flux of liquid water in updrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: mful
    ! Updraft winds m/s
    real(rkx), dimension(np,nk) :: uu, vu
    ! Downdraft winds m/s
    real(rkx), dimension(np,nk) :: ud, vd
    ! Updraft kinetic energy m^2/s^2
    real(rkx), dimension(np,nk) :: kineu
    ! Downdraft kinetic energy m^2/s^2
    real(rkx), dimension(np,nk) :: kined
    ! Precipitation Flux kg/(m^2*s)
    real(rkx), dimension(np) :: rfl
    ! Mass Flux Updraft base of the cloud
    real(rkx), dimension(np) :: mfub, mfub1
    ! Change in precipitation fluxes due to melting kg/(m^2*s)
    real(rkx), dimension(np,nk) :: dpmel
    ! Precipitation flux
    real(rkx), dimension(np,nk) :: xrain
    ! Flux of frozen cloudwater in updrafts kg/(m^2*s)
    real(rkx), dimension(np,nk) :: lglac
    real(rkx), dimension(np) :: dqcv
    real(rkx), dimension(np) :: dhpbl
    real(rkx), dimension(np) :: wubase
    ! Entrainment, detrainment
    real(rkx), dimension(np,nk) :: dmfen, dmfde
    ! Flag
    integer(ik4), dimension(np,nk) :: ilab
    ! Level of max wind
    integer(ik4), dimension(np) :: ilwmin
    integer(ik4), dimension(np) :: idtop, ictop0
    ! Departure level for convection
    integer(ik4), dimension(np) :: idpl
    real(rkx), dimension(np) :: xcape, rheat
    logical, dimension(np) :: lddraf, llddraf3, lldcum, llo2
    logical :: llo1
    integer(ik4) :: ikb, itopm2, k, ik, n
    real(rkx) :: cons2, cons, dh, dqmin, dz, eps, fac, &
                 mfmax, pbmpt, qumqe, xro, mfa, erate,  &
                 derate, duten, dvten, tdis, xalv
    real(rkx), dimension(np) :: sfl
    real(rkx), dimension(np) :: tau ! adjustment time
    ! scaling factor for momentum and tracer massflux
    real(rkx), dimension(np) :: mfs, mfuub, mfuvb, xsum12, xsum22
    real(rkx), dimension(np) :: mf_shal
    real(rkx), dimension(np,nk) :: mfuus, mfdus, mfudr, xtent, &
                                    mfddr, xtenu, xtenv, xtenq, &
                                    uv2
    logical :: llconscheck = .false.
    integer(ik4) :: nt
    real(rkx), dimension(:,:), allocatable :: xsumc
    real(rkx), dimension(:,:,:), allocatable :: xtenc
    real(rkx) :: tfl

    ! ---------------------------------------
    ! 0. Compute Saturation specific humidity
    ! ---------------------------------------
    !$acc kernels
    ldcum(:) = .false.
    qs(:,:) = q(:,:)
    !$acc end kernels
    call satur(n1,n2,nk,nk060,ph,t,qs)
    !------------------------------------
    ! 1. Specify constants and parameters
    ! -----------------------------------
    cons2 = rmfcfl/(egrav*dtc)
    cons = d_one/(egrav*dtc)
    !---------------------------------------------
    ! 2. Initialize values at vertical grid points
    ! --------------------------------------------
    do concurrent ( n = n1:n2, k = 1:nk )
      xtent(n,k) = tent(n,k)
      xtenq(n,k) = tenq(n,k)
      xtenu(n,k) = tenu(n,k)
      xtenv(n,k) = tenv(n,k)
    end do
    !
    ! This routine interpolates large-scale fields of t,q etc.
    ! to sigma levels (i.e. grid for massflux scheme), determines level
    ! of maximum vertical velocity and initializes values for updrafts
    ! and downdrafts interface
    !
    initcum: block
      real(rkx), dimension(np) :: wmax
      real(rkx), dimension(np) :: xph
      logical, dimension(np) :: llflag
      integer(ik4) :: icall, ik, k, n, km1
      real(rkx) :: zs
      !--------------------------------------------------
      ! 1. Specify large scale parameters at half levels
      !    adjust temperature fields if staticly unstable
      !    find level of maximum vertical velocity
      ! -------------------------------------------------
      do k = 2, nk
        do concurrent ( n = n1:n2 )
          tf(n,k) = (max(cpd*t(n,k-1) + &
            geo(n,k-1),cpd*t(n,k)+geo(n,k))-geof(n,k))*rcpd
          qf(n,k) = q(n,k-1)
          qsf(n,k) = qs(n,k-1)
          xph(n) = pf(n,k)
          llflag(n) = .true.
        end do
        if ( k >= nk-1 .or. k < nk060 ) cycle
        ik = k
        icall = 3
        call moistadj(ik,xph,tf,qsf,llflag,icall)
        do concurrent ( n = n1:n2 )
          qf(n,k) = min(q(n,k-1),qs(n,k-1))+(qsf(n,k)-qs(n,k-1))
          qf(n,k) = max(qf(n,k),d_zero)
        end do
      end do
      do concurrent ( n = n1:n2 )
        tf(n,nk) = (cpd*t(n,nk)+geo(n,nk)-geof(n,nk))*rcpd
        qf(n,nk) = q(n,nk)
        tf(n,1) = t(n,1)
        qf(n,1) = q(n,1)
        ilwmin(n) = nk
        wmax(n) = d_zero
      end do
      do concurrent ( n = n1:n2, k = nk-1:2:-1 )
        zs = max(cpd*tf(n,k)+geof(n,k),cpd*tf(n,k+1)+geof(n,k+1))
        tf(n,k) = (zs-geof(n,k))*rcpd
      end do
      do k = nk, 3, -1
        do concurrent ( n = n1:n2 )
          if ( omega(n,k) < wmax(n) ) then
            wmax(n) = omega(n,k)
            ilwmin(n) = k
          end if
        end do
      end do
      !-------------------------------------------------
      ! 2. Initialize values for updrafts and downdrafts
      ! ------------------------------------------------
      do concurrent ( n = n1:n2, k = 1:nk )
        km1 = max(k-1,1)
        tu(n,k) = tf(n,k)
        td(n,k) = tf(n,k)
        qu(n,k) = qf(n,k)
        qd(n,k) = qf(n,k)
        lu(n,k) = d_zero
        uu(n,k) = u(n,km1)
        ud(n,k) = u(n,km1)
        vu(n,k) = v(n,km1)
        vd(n,k) = v(n,km1)
        ilab(n,k) = 0
        cvrainout(n,k) = d_zero
      end do
    end block initcum
    !---------------------------
    ! 3. Cloud base calculations
    ! --------------------------
    ! ---------------------------------
    !   (A) Determine cloud base values
    ! ---------------------------------
    !
    ! This routine calculates cloud base fields, cloud base height and cloud
    ! top height. Produce cloud base and cloud top values for parametrization
    ! It returns cloud fields values and flags as follows:
    !    ilab = 0 for stable layers
    !    ilab = 1 for subcloud levels
    !    ilab = 2 for cloud levels level
    ! Include cycling over levels to find unstable departure/base level+
    ! mixed layer properties +w Trigger
    ! Lift surface air dry-adiabatically to cloud top
    ! (entraining plume, with entrainment proportional to (1/z))
    !
    cloudbase: block
      integer(ik4), dimension(np) :: ictop, icbot, ibotsc, iidpl
      integer(ik4), dimension(np,nk) :: iilab
      logical, dimension(np) :: ll_ldbase, llgo_on, lldeep, lldcum, &
                lldsc, llfirst, llresetn
      logical :: llreset
      integer(ik4) :: icall, ik, is, k, n, kk, kt1, kt2, kt, kb
      real(rkx), dimension(np,nk) :: suh, wu2h, buoh
      real(rkx), dimension(np,nk+1) :: xsenh, xqenh
      real(rkx), dimension(np) :: qold, xph
      real(rkx), dimension(np) :: zmix
      real(rkx), dimension(np) :: dz, cbase
      real(rkx), dimension(np,nk) :: xlu, xqu, xtu, xuu, xvu
      ! local for CAPE at every departure level
      real(rkx), dimension(np,nk) :: xcape
      real(rkx) :: buof, c2, epsadd
      real(rkx) :: rho    ! Density at surface (kg/m^3)
      real(rkx) :: khvfl  ! Surface buoyancy flux (K m/s)
      real(rkx) :: ws     ! Sigma_w at lowest model level (m/s)
      real(rkx) :: qexc   ! Humidity excess at lowest model level (kg/kg)
      real(rkx) :: texc   ! Temperature excess at lowest model level (K)
      real(rkx) :: eps    ! Fractional entrainment rate [m^-1]
      real(rkx) :: tvenh  ! Environment virtual temperature (K)
      real(rkx) :: tvuh   ! Updraft virtual temperature (K)
      real(rkx) :: xlglac ! Updraft liquid water frozen in one layer
      real(rkx) :: qsu, cor, dq, alfaw, facw, faci, fac,     &
                   esdp, dqsdt, dtdp, dp, pdifftop, pdiffbot, &
                   sf, xqf, xaw, xbw
      real(rkx), dimension(np) :: tven1, tvu1
      real(rkx) :: tven2, tvu2 ! pseudoadiabatic T_v
      real(rkx) :: work1, work2 ! work for T and w perturbations
      real(rkx) :: xtmp, entrpen, mymax
      !-----------------------------------
      ! 0. Initialize constants and fields
      ! ----------------------------------
      c2 = 0.55_rkx
      xaw = d_one
      xbw = d_one
      epsadd = 1.e-4_rkx
      do concurrent ( n = n1:n2 )
        wubase(n) = d_zero
        llgo_on(n) = .true.
        llfirst(n) = .true.
        idpl(n) = nk
      end do
      kt1 = nk350
      kt2 = nk060
      do concurrent ( n = n1:n2, k = 1:nk )
        xtu(n,k) = tu(n,k)
        xqu(n,k) = qu(n,k)
        xlu(n,k) = lu(n,k)
        xuu(n,k) = uu(n,k)
        xvu(n,k) = vu(n,k)
        iilab(n,k) = ilab(n,k)
        xcape(n,k) = d_zero
      end do
      !-------------------------------------------
      ! 1. Prepare fields by linear interpolation
      !    of specific humidity and static energy
      ! ------------------------------------------
      do concurrent ( n = n1:n2, k = 1:nk )
        wu2h(n,k) = d_zero
        xqenh(n,k) = qf(n,k)
        xsenh(n,k) = cpd*tf(n,k) + geof(n,k)
      end do
      do kk = nk, kt1, -1
        !
        ! Find first departure level that produces deepest cloud top
        ! or take surface level for shallow convection and Sc
        !
        ! ----------------------------------------------
        ! 1.2 Initialise fields at departure model level
        ! ----------------------------------------------
        !
        is = 0
#ifdef STDPAR
        do concurrent ( n = n1:n2 ) reduce(+:is)
#else
        do n = n1, n2
#endif
          if ( llgo_on(n) ) then
            is = is + 1
            iidpl(n) = kk   ! departure level
            icbot(n) = kk   ! cloud base level for convection,
                            ! (-1 if not found)
            ibotsc(n) = nk - 1
                            ! sc    base level for sc-clouds,
                            ! (-1 if not found)
            ictop(n) = nk - 1
                            ! cloud top for convection (-1 if not found)
            lldcum(n) = .false.
                            ! on exit: true if cloudbase=found
            lldsc(n) = .false.
                            ! on exit: true if cloudbase=found
            ll_ldbase(n) = .false.
                            ! on exit: true if cloudbase=found
            xuu(n,kk) = u(n,kk)*(pf(n,kk+1)-pf(n,kk))
            xvu(n,kk) = v(n,kk)*(pf(n,kk+1)-pf(n,kk))
          end if
        end do
        if ( is /= 0 ) then
          if ( kk == nk ) then
            do concurrent ( n = n1:n2 )
              if ( llgo_on(n) ) then
                rho = pf(n,kk+1)/(rgas*(t(n,kk) * (d_one+ep1*q(n,kk))))
                khvfl = (ahfs(n,kk+1)*rcpd + ep1*t(n,kk)*qhfl(n,kk+1))/rho
                ws = 0.001_rkx - 1.5_rkx*rkap*khvfl * &
                         (geof(n,nk)-geof(n,nk+1))/t(n,nk)
                if ( khvfl < d_zero ) then
                  ws = 1.2_rkx*ws**.3333_rkx
                  iilab(n,kk) = 1
                  texc = max(-1.5_rkx*ahfs(n,kk+1)/(rho*ws*cpd),d_zero)
                  qexc = max(-1.5_rkx*qhfl(n,kk+1)/(rho*ws),d_zero)
                  xqu(n,kk) = xqenh(n,kk) + qexc
                  suh(n,kk) = xsenh(n,kk) + cpd*texc
                  xtu(n,kk) = (xsenh(n,kk)-geof(n,kk))*rcpd + texc
                  xlu(n,kk) = d_zero
                  wu2h(n,kk) = ws**2
                  !
                  ! Determine buoyancy at lowest half level
                  !
                  tvenh = (d_one+ep1*xqenh(n,kk)) * &
                           (xsenh(n,kk)-geof(n,kk))*rcpd
                  tvuh = (d_one+ep1*xqu(n,kk))*xtu(n,kk)
                  buoh(n,kk) = (tvuh-tvenh)*egrav/tvenh
                else
                  llgo_on(n) = .false. ! Non-convective point
                end if
              end if
            end do
          else
            do concurrent ( n = n1:n2 )
              if ( llgo_on(n) ) then
                rho = pf(n,kk+1) / (rgas*(t(n,kk)*(d_one+ep1*q(n,kk))))
                iilab(n,kk) = 1
                texc = 0.2_rkx
                qexc = 1.e-4_rkx
                xqu(n,kk) = xqenh(n,kk) + qexc
                suh(n,kk) = xsenh(n,kk) + cpd*texc
                xtu(n,kk) = (xsenh(n,kk)-geof(n,kk))*rcpd + texc
                xlu(n,kk) = d_zero
                ! Construct mixed layer for parcels emanating in lowest 60 hPa
                if ( pf(n,nk+1)-pf(n,kk-1) < 60.e2_rkx ) then
                  xqu(n,kk) = d_zero
                  suh(n,kk) = d_zero
                  work1 = d_zero
                  do k = kk + 1, kk - 1, -1
                    if ( work1 < 50.e2_rkx ) then
                      work2 = pf(n,k) - pf(n,k-1)
                      work1 = work1 + work2
                      xqu(n,kk) = xqu(n,kk) + xqenh(n,k)*work2
                      suh(n,kk) = suh(n,kk) + xsenh(n,k)*work2
                    end if
                  end do
                  xqu(n,kk) = xqu(n,kk)/work1 + qexc
                  suh(n,kk) = suh(n,kk)/work1 + cpd*texc
                  xtu(n,kk) = (suh(n,kk)-geof(n,kk))*rcpd + texc
                end if
                wu2h(n,kk) = d_one
                !
                ! Determine buoyancy at lowest level
                !
                tvenh = (d_one+ep1*xqenh(n,kk)) * &
                         (xsenh(n,kk)-geof(n,kk))*rcpd
                tvuh = (d_one+ep1*xqu(n,kk))*xtu(n,kk)
                buoh(n,kk) = (tvuh-tvenh)*egrav/tvenh
              end if
            end do
          end if
        end if
        !----------------------------------------------
        ! 2. Do ascent in subcloud and layer,
        !    check for existence of condensation level,
        !    adjust t,q and l accordingly
        !    Check for buoyancy and set flags
        ! ---------------------------------------------
        ! ----------------------------------------------------------
        ! 2.1 Do the vertical ascent until velocity becomes negative
        ! ----------------------------------------------------------
        do k = kk - 1, kt2, -1
          is = 0
          if ( kk == nk ) then
            ! 1/z mixing for shallow
#ifdef STDPAR
            do concurrent ( n = n1:n2 ) reduce(+:is)
#else
            do n = n1, n2
#endif
              if ( llgo_on(n) ) then
                is = is + 1
                dz(n) = (geof(n,k)-geof(n,k+1))*regrav
                eps = c2/((geof(n,k)-geof(n,nk+1))*regrav) + epsadd
                zmix(n) = d_half*dz(n)*eps
                xqf = (qf(n,k+1)+qf(n,k))*d_half
                sf = (xsenh(n,k+1)+xsenh(n,k))*d_half
                xtmp = d_one/(d_one+zmix(n))
                xqu(n,k) = (xqu(n,k+1)*(d_one-zmix(n)) + &
                  d_two*zmix(n)*xqf)*xtmp
                suh(n,k) = (suh(n,k+1)*(d_one-zmix(n)) + &
                  d_two*zmix(n)*sf)*xtmp
                qold(n) = xqu(n,k)
                xtu(n,k) = (suh(n,k)-geof(n,k))*rcpd
                xph(n) = pf(n,k)
              end if
            end do
          else
#ifdef STDPAR
            do concurrent ( n = n1:n2 ) reduce(+:is)
#else
            do n = n1, n2
#endif
              if ( llgo_on(n) ) then
                is = is + 1
                dz(n) = (geof(n,k)-geof(n,k+1))*regrav
                xqf = (qf(n,k+1)+qf(n,k))*d_half
                sf = (xsenh(n,k+1)+xsenh(n,k))*d_half
                entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(n))
                zmix(n) = 0.4_rkx*entrpen*dz(n) * &
                           min(d_one,(qs(n,k)/qs(n,nk))**3)
                xqu(n,k) = xqu(n,k+1)*(d_one-zmix(n)) + xqf*zmix(n)
                suh(n,k) = suh(n,k+1)*(d_one-zmix(n)) + sf*zmix(n)
                qold(n) = xqu(n,k)
                xtu(n,k) = (suh(n,k)-geof(n,k))*rcpd
                xph(n) = pf(n,k)
              end if
            end do
          end if
          if ( is == 0 ) exit
          ik = k
          icall = 1
          call moistadj(ik,xph,xtu,xqu,llgo_on,icall)
          do concurrent ( n = n1:n2 )
            if ( llgo_on(n) ) then
              ! Add condensation to water
              dq = max(qold(n)-xqu(n,k),d_zero)
              xlu(n,k) = xlu(n,k+1) + dq
              ! Freezing
              xlglac = dq*((d_one-xalpha(xtu(n,k))) - &
                       (d_one-xalpha(xtu(n,k+1))))
              ! Pseudo-Microphysics
              if ( kk == nk ) then
                ! No precip for shallow
                xlu(n,k) = min(xlu(n,k),5.e-3_rkx)
                ! Chose a more pseudo-adiabatic formulation as
                ! original overestimates water loading effect and
                ! therefore strongly underestimates cloud thickness
              else
                xlu(n,k) = d_half*xlu(n,k)
              end if
              ! Update dry static energy after condensation + freezing
              suh(n,k) = cpd*(xtu(n,k)+wlhfocp*xlglac) + geof(n,k)
              ! Buoyancy on half and full levels
              tvuh = (d_one+ep1*xqu(n,k)-xlu(n,k))*xtu(n,k) + &
                      wlhfocp*xlglac
              tvenh = (d_one+ep1*xqenh(n,k)) * (xsenh(n,k)-geof(n,k))*rcpd
              buoh(n,k) = (tvuh-tvenh)*egrav/tvenh
              buof = (buoh(n,k)+buoh(n,k+1))*d_half
              ! Solve kinetic energy equation
              xtmp = d_one/(d_one+d_two*xbw*zmix(n))
              wu2h(n,k) = (wu2h(n,k+1)*(d_one-d_two*xbw*zmix(n)) + &
                          d_two*xaw*buof*dz(n))*xtmp
              ! Compute pseudoadiabatique CAPE for diagnostics
              tvu2 = xtu(n,k)*(d_one+ep1*xqu(n,k))
              tven2 = tf(n,k)*(d_one+ep1*qf(n,k))
              if ( k == kk-1 ) then
                tvu1(n) = tvu2
                tven1(n) = tven2
              end if
              buof = (tvu2+tvu1(n)-tven1(n)-tven2)/tven2
              buof = buof*dz(n)*egrav
              xcape(n,kk) = xcape(n,kk) + max(d_zero,buof)
              tvu1(n) = tvu2
              tven1(n) = tven2
              ! First layer with liquid water - find exact cloud base
              if ( xlu(n,k) > d_zero .and. iilab(n,k+1) == 1 ) then
                ik = k + 1
                qsu = fesat(xtu(n,ik))/pf(n,ik)
                qsu = min(qsmax,qsu)
                cor = d_one/(d_one-ep1*qsu)
                qsu = qsu*cor
                dq = min(d_zero,xqu(n,ik)-qsu)
                alfaw = xalpha(xtu(n,ik))
                facw = c5les/((xtu(n,ik)-c4les)**2)
                faci = c5ies/((xtu(n,ik)-c4ies)**2)
                fac = alfaw*facw + (d_one-alfaw)*faci
                esdp = fesat(xtu(n,ik))/pf(n,ik)
                cor = d_one/(d_one-ep1*esdp)
                dqsdt = fac*cor*qsu
                dtdp = rgas*xtu(n,ik)/(cpd*pf(n,ik))
                dp = dq/(dqsdt*dtdp)
                cbase(n) = pf(n,ik) + dp
                ! Chose nearest half level as cloud base
                pdifftop = cbase(n) - pf(n,k)
                pdiffbot = pf(n,k+1) - cbase(n)
                if ( pdifftop > pdiffbot .and. wu2h(n,k+1) > d_zero ) then
                  kb = min(nk-1,k+1)
                  iilab(n,kb) = 2
                  iilab(n,k) = 2
                  ll_ldbase(n) = .true.
                  lldsc(n) = .true.
                  ibotsc(n) = kb
                  icbot(n) = kb
                  xlu(n,k+1) = 0.0_rkx
                else if ( pdifftop <= pdiffbot .and. &
                          wu2h(n,k) > d_zero ) then
                  iilab(n,k) = 2
                  ll_ldbase(n) = .true.
                  lldsc(n) = .true.
                  ibotsc(n) = k
                  icbot(n) = k
                end if
                kb = icbot(n)
              end if
              ! Decide on presence of convection, cloud base and cloud
              ! top based on kinetic energy
              if ( wu2h(n,k) < d_zero ) then
                llgo_on(n) = .false.
                if ( xlu(n,k+1) > d_zero ) then
                  ictop(n) = k
                  lldcum(n) = .true.
                else
                  lldcum(n) = .false.
                end if
              else if ( xlu(n,k) > d_zero ) then
                iilab(n,k) = 2
              else
                iilab(n,k) = 1
              end if
            end if
          end do
          if ( lmfdudv .and. kk == nk ) then
            do concurrent ( n = n1:n2 )
              if ( .not. ll_ldbase(n) .and. llgo_on(n) ) then
                xuu(n,kk) = xuu(n,kk) + &
                  u(n,k)*(pf(n,k+1)-pf(n,k))
                xvu(n,kk) = xvu(n,kk) + &
                  v(n,k)*(pf(n,k+1)-pf(n,k))
              end if
            end do
          end if
        end do
        if ( kk == nk ) then
          ! Set values for departure level for PBL clouds = first model level
          do concurrent ( n = n1:n2 )
            ldsc(n) = lldsc(n)
            if ( ldsc(n) ) then
              kbotsc(n) = ibotsc(n)
            else
              kbotsc(n) = -1
            end if
            llgo_on(n) = .false.
            kt = ictop(n)
            kb = icbot(n)
            lldeep(n) = pf(n,kb) - pf(n,kt)>rdepths
            if ( lldeep(n) ) lldcum(n) = .false. ! No deep allowed for KLEV
            lldeep(n) = .false.
            ! For deep convection start only at level NK-1
            ! and form mixed layer, so go on
            ! Test further for deep convective columns as not yet found
            if ( lldeep(n) ) llfirst(n) = .false.
            llgo_on(n) = .not. lldeep(n)
            if ( lldcum(n) ) then
              kcbot(n) = icbot(n)
              ictop0(n) = ictop(n)
              idpl(n) = iidpl(n)
              ldcum(n) = lldcum(n)
              wubase(n) = sqrt(max(wu2h(n,kb),d_zero))
            else
              ictop0(n) = -1
              kcbot(n) = -1
              idpl(n) = nk - 1
              ldcum(n) = .false.
              wubase(n) = d_zero
            end if
          end do
          do k = nk, 1, -1
            do concurrent ( n = n1:n2 )
              kt = ictop(n)
              if ( k >= kt ) then
                ilab(n,k) = iilab(n,k)
                tu(n,k) = xtu(n,k)
                qu(n,k) = xqu(n,k)
                lu(n,k) = xlu(n,k)
              end if
            end do
          end do
        end if
        if ( kk < nk ) then
          llreset = .false.
          do concurrent ( n = n1:n2 )
            if ( .not.lldeep(n) ) then
              kt = ictop(n)
              kb = icbot(n)
              ! Test on cloud thickness and buoyancy
              lldeep(n) = pf(n,kb) - pf(n,kt)>=rdepths
            end if
            llresetn(n) = lldeep(n) .and. llfirst(n)
            llreset = llreset .or. llresetn(n)
          end do
          if ( llreset ) then
            do concurrent ( n = n1:n2, k = nk:1:-1 )
              ! Keep first departure level that produces deep cloud
              if ( llresetn(n) ) then
                kt = ictop(n)
                kb = iidpl(n)
                if ( k <= kb .and. k >= kt ) then
                  ilab(n,k) = iilab(n,k)
                  tu(n,k) = xtu(n,k)
                  qu(n,k) = xqu(n,k)
                  lu(n,k) = xlu(n,k)
                else
                  ilab(n,k) = 1
                  tu(n,k) = tf(n,k)
                  qu(n,k) = qf(n,k)
                  lu(n,k) = d_zero
                end if
                if ( k < kt ) ilab(n,k) = 0
              end if
            end do
          end if
          do concurrent ( n = n1:n2 )
            if ( lldeep(n) .and. llfirst(n) ) then
              idpl(n) = iidpl(n)
              ictop0(n) = ictop(n)
              kcbot(n) = icbot(n)
              ldcum(n) = lldcum(n)
              ldsc(n) = .false.
              kbotsc(n) = -1
              kb = kcbot(n)
              wubase(n) = sqrt(max(wu2h(n,kb),d_zero))
              ! No initialization of wind for deep here
              llfirst(n) = .false.
            end if
            llgo_on(n) = .not. lldeep(n)
          end do
        end if
      end do ! End of big loop for search of departure level
      ! Chose maximum CAPE value
      do concurrent ( n = n1:n2 )
        mymax = xcape(n,1)
        !$acc loop seq
        do k = 2, nk
          mymax = max(mymax,xcape(n,k))
        end do
        cape(n) = mymax
      end do
    end block cloudbase
    !-----------------------------------------------------------
    !   (B) Determine total moisture convergence and
    !       decide on type of cumulus convection on the basis of
    !       the depth of the convection
    !           deep    if cloud depth > 200mb
    !           shallow if cloud depth < 200mb
    ! ----------------------------------------------------------
    ! Calculate column and sub cloud layer moisture convergence
    ! and sub cloud layer moist static energy convergence
    do concurrent ( n = n1:n2 )
      dqcv(n) = d_zero
      dhpbl(n) = d_zero
      idtop(n) = 0
    end do
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = nk060, nk
        dqcv(n) = dqcv(n) + max(d_zero,tenq(n,k)) * (pf(n,k+1)-pf(n,k))
        if ( ldcum(n) .and. k >= kcbot(n) ) then
          dhpbl(n) = dhpbl(n) + &
            (wlhv*tenq(n,k)+cpd*tent(n,k))*(pf(n,k+1)-pf(n,k))
        end if
      end do
    end do
    !-----------------------------------------------------
    ! Estimate cloud height for entrainment/detrainment
    ! Calculations and initial determination of cloud type
    ! ----------------------------------------------------
    ! Specify initial cloud type
    !
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) ) then
        ikb = kcbot(n)
        itopm2 = ictop0(n)
        pbmpt = pf(n,ikb) - pf(n,itopm2)
        if ( pbmpt >= rdepths ) then
          ktype(n) = 1
        else
          ktype(n) = 2
        end if
      else
        ktype(n) = 0
      end if
    end do
    !-----------------------------------------------------------------------
    ! (C) calculate initial updraught mass flux and set lateral mixing rates
    !     for deep convection assume it is 10% of maximum value
    !     which is determined by the thickness of the layer and timestep
    !     for shallow convection calculated assuming a balance of moist
    !     static energy in the sub-cloud layer (ignores present of
    !     downdraughts)
    ! ----------------------------------------------------------------------
    if ( lmfwstar ) then
      do concurrent ( n = n1:n2 )
        if ( ldcum(n) ) then
          ikb = kcbot(n)
          dz = max(d_zero,min(1.5e3_rkx,(geof(n,ikb)-geof(n,nk+1))*regrav))
          mf_shal(n) = 0.07_rkx*(egrav/t(n,nk)*dz *      &
                       max(d_zero,-ahfs(n,nk+1)*rcpd - &
                       ep1*t(n,nk)*qhfl(n,nk+1)))**.3333_rkx
          mfmax = (pf(n,ikb)-pf(n,ikb-1))*cons2*rmflic + rmflia
          mf_shal(n) = min(mf_shal(n),mfmax)
        end if
      end do
    end if
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) ) then
        ikb = kcbot(n)
        mfmax = (pf(n,ikb)-pf(n,ikb-1))*cons2*rmflic + rmflia
        !     deep convection
        if ( ktype(n) == 1 ) then
          mfub(n) = mfmax*0.1_rkx
        else if ( ktype(n) == 2 ) then
          !       shallow convection
          qumqe = qu(n,ikb) + lu(n,ikb) - qf(n,ikb)
          dqmin = max(0.01_rkx*qf(n,ikb),1.0e-10_rkx)
          dh = cpd*(tu(n,ikb)-tf(n,ikb)) + wlhv*qumqe
          dh = egrav*max(dh,1.e5_rkx*dqmin)
          if ( dhpbl(n) > d_zero ) then
            mfub(n) = dhpbl(n)/dh
            mfub(n) = min(mfub(n),mfmax)
          else
            mfub(n) = mfmax*0.1_rkx
            ldcum(n) = .false.
          end if
          if ( lmfwstar ) mfub(n) = mf_shal(n)
        end if
      else
        ! No buoyancy cloud base from surface: set cloud base mass flux
        ! and mixing rate to default value for safety
        mfub(n) = d_zero
      end if
    end do
    !-----------------------------------------------
    ! 4. Determine cloud ascent for entraining plume
    ! ----------------------------------------------
    !  (B) Do ascent in absence of downdrafts
    ! --------------------------------------------
    !
    ! This routine does the calculations for cloud ascents, i.e. produce
    ! cloud ascents vertical profiles of T,Q,L,U AND V and corresponding
    ! fluxes as well as precipitation rates.
    ! Lift surface air dry-adiabatically to cloud base and then calculate
    ! moist ascent for entraining/detraining plume.
    ! Entrainment and detrainment rates differ for shallow and deep cumulus
    ! convection. In case there is no penetrative or shallow convection
    ! check for possibility of mid level convection.
    !
    ascent: block
      real(rkx), dimension(np,nk) :: buo
      real(rkx), dimension(np) :: ndmfen, ndmfde, qold, luold, precip
      real(rkx), dimension(np) :: dpmean
      real(rkx), dimension(np) :: zoentr, xph
      logical, dimension(np) :: llflag, llo1
      logical :: llo3
      integer(ik4) :: icall, ik, is, k, n, ikb
      integer(ik4) :: nll, nlm, mynlm
      integer(ik4), dimension(np) :: nlx
      real(rkx) :: cldmax, cprc2, cwifrac, alfaw, bc,   &
                   be, buoc, zc, cbf, cons2, zd, zdfi, dkbuo, &
                   dken, dnoprc, zdt, fac, facbuo, zint,       &
                   kedke, lcrit, leen, lnew, mfmax, mftest,     &
                   mfulk, mfun, mfuqk, mfusk, prcdgw, prcon,    &
                   qeen, qude, rnew, rold, scde, seen, &
                   zvi, zvv, zvw, zwu, zco, arg, entrpen
      real(rkx) :: change, zxs, zxe
      logical, dimension(np) :: llklab
      !----------------------
      ! 1. Specify parameters
      ! ---------------------
      cons2 = rmfcfl/(egrav*dtc)
      facbuo = d_half/(d_one+d_half)
      cldmax = 5.e-3_rkx
      cwifrac = d_half
      cprc2 = d_half
      !----------------------
      ! 2. Set default values
      ! ---------------------
      llo3 = .false.
      do concurrent ( n = n1:n2 )
        luold(n) = d_zero
        if ( .not. ldcum(n) ) then
          kcbot(n) = -1
          mfub(n) = d_zero
          qu(n,nk) = d_zero
          ktype(n) = 0
        end if
        wmean(n) = d_zero
        dpmean(n) = d_zero
        zoentr(n) = d_zero
      end do
      ! Initalize various quantities
      ! Note that liquid water and kinetic energy at cloud base is preserved
      do concurrent ( n = n1:n2 )
        llklab(n) = .false.
        if ( .not. ldcum(n) .or. ktype(n) == 3 ) llklab(n) = .true.
      end do
      do concurrent ( n = n1:n2, k = 1:nk )
        if ( k /= kcbot(n) ) lu(n,k) = d_zero
        kineu(n,k) = d_zero
        mfu(n,k) = d_zero
        mfus(n,k) = d_zero
        mfuq(n,k) = d_zero
        mful(n,k) = d_zero
        lude(n,k) = d_zero
        lglac(n,k) = d_zero
        dmfup(n,k) = d_zero
        xrain(n,k) = d_zero
        buo(n,k) = d_zero
        if ( llklab(n) ) ilab(n,k) = 0
        dmfen(n,k) = d_zero
        mfude_rate(n,k) = d_zero
      end do
      do concurrent ( n = n1:n2 )
        !$acc loop seq
        do k = 1, nk
          if ( .not. ldcum(n) .and. ph(n,k) < 4.0e4_rkx ) ictop0(n) = k
        end do
      end do
      do concurrent ( n = n1:n2 )
        if ( ktype(n) == 3 ) ldcum(n) = .false.
      end do
      !-----------------------------------------
      ! 3. Initialize values at cloud base level
      ! ----------------------------------------
      do concurrent ( n = n1:n2 )
        kctop(n) = kcbot(n)
        if ( ldcum(n) ) then
          ikb = kcbot(n)
          kineu(n,ikb) = d_half*wubase(n)**2
          mfu(n,ikb) = mfub(n)
          mfus(n,ikb) = mfub(n)*(cpd*tu(n,ikb)+geof(n,ikb))
          mfuq(n,ikb) = mfub(n)*qu(n,ikb)
          mful(n,ikb) = mfub(n)*lu(n,ikb)
        end if
      end do
      !-------------------------------------------------------
      ! 4. Do ascent: subcloud layer (ilab=1) ,clouds (ilab=2)
      !    By doing first dry-adiabatic ascent and then
      !    by adjusting t,q and l accordingly, then check for
      !    buoyancy and set flags accordingly
      ! ------------------------------------------------------
      do k = nk - 1, 3, -1
        ! Specify cloud base values for midlevel convection
        ! in case there is not already convection
        ! ----------------------------------------------------
        ik = k
        call mcbase(ik)
        is = 0
        nlm = 0
        do concurrent ( n = n1:n2 )
          llflag(n) = .false.
          precip(n) = d_zero
          llo1(n) = .false.
          is = is + ilab(n,k+1)
          if ( ilab(n,k+1) == 0 ) ilab(n,k) = 0
          if ( (ldcum(n) .and. ilab(n,k+1) == 2) .or.   &
               (ktype(n) == 3 .and. ilab(n,k+1) == 1) ) then
            llflag(n) = .true.
            !$acc atomic capture
            nlm = nlm + 1
            mynlm = nlm
            !$acc end atomic
            nlx(mynlm) = n
          end if
          xph(n) = ph(n,k)
          if ( ktype(n) == 3 .and. k == kcbot(n) ) then
            mfmax = (pf(n,k)-pf(n,k-1))*cons2*rmflic + rmflia
            if ( mfub(n) > mfmax ) then
              fac = mfmax/mfub(n)
              mfu(n,k+1) = mfu(n,k+1)*fac
              mfus(n,k+1) = mfus(n,k+1)*fac
              mfuq(n,k+1) = mfuq(n,k+1)*fac
              mfub(n) = mfmax
            end if
          end if
        end do
        if ( is > 0 ) llo3 = .true.
        !--------------------------
        ! Specify entrainment rates
        ! -------------------------
        ik = k
        call entrainm(ik,llo3,ndmfen,ndmfde)
        !----------------------------------------------------
        ! Do adiabatic ascent for entraining/detraining plume
        ! ---------------------------------------------------
        if ( llo3 ) then
          do concurrent ( n = n1:n2 )
            qold(n) = d_zero
          end do
          do concurrent ( nll = 1:nlm )
            n = nlx(nll)
            ndmfde(n) = min(ndmfde(n),0.75_rkx*mfu(n,k+1))
            if ( k == kcbot(n) ) then
              entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(n))
              zoentr(n) = -entrpen*(min(d_one,q(n,k)/qs(n,k)) - &
                           d_one)*(geof(n,k)-geof(n,k+1))*regrav
              zoentr(n) = min(0.4_rkx,zoentr(n))*mfu(n,k+1)
            end if
            if ( k < kcbot(n) ) then
              mfmax = (pf(n,k)-pf(n,k-1))*cons2*rmflic + rmflia
              zxs = max(mfu(n,k+1)-mfmax,d_zero)
              wmean(n) = wmean(n) + kineu(n,k+1)*(ph(n,k+1)-ph(n,k))
              dpmean(n) = dpmean(n) + ph(n,k+1) - ph(n,k)
              ndmfen(n) = zoentr(n)
              if ( ktype(n) >= 2 ) then
                ndmfen(n) = entshalp*ndmfen(n)
                ndmfde(n) = ndmfen(n)
              end if
              ndmfde(n) = ndmfde(n) * (1.6_rkx-min(d_one,q(n,k)/qs(n,k)))
              mftest = mfu(n,k+1) + ndmfen(n) - ndmfde(n)
              change = max(mftest-mfmax,d_zero)
              zxe = max(change-zxs,d_zero)
              ndmfen(n) = ndmfen(n) - zxe
              change = change - zxe
              ndmfde(n) = ndmfde(n) + change
            end if
            dmfen(n,k) = ndmfen(n) - ndmfde(n)
            mfu(n,k) = mfu(n,k+1) + ndmfen(n) - ndmfde(n)
            qeen = qf(n,k+1)*ndmfen(n)
            seen = (cpd*tf(n,k+1)+geof(n,k+1))*ndmfen(n)
            if ( qctot(n,k) > 1.0e-10_rkx ) then
              leen = qctot(n,k)*ndmfen(n)
            else
              leen = d_zero
            end if
            scde = (cpd*tu(n,k+1)+geof(n,k+1))*ndmfde(n)
            qude = qu(n,k+1)*ndmfde(n)
            lude(n,k) = lu(n,k+1)*ndmfde(n)
            mfusk = mfus(n,k+1) + seen - scde
            mfuqk = mfuq(n,k+1) + qeen - qude
            mfulk = mful(n,k+1) + leen - lude(n,k)
            lu(n,k) = mfulk*(d_one/max(cmfcmin,mfu(n,k)))
            qu(n,k) = mfuqk*(d_one/max(cmfcmin,mfu(n,k)))
            tu(n,k) = (mfusk * &
              (d_one/max(cmfcmin,mfu(n,k)))-geof(n,k))*rcpd
            tu(n,k) = max(100.0_rkx,tu(n,k))
            tu(n,k) = min(400.0_rkx,tu(n,k))
            qold(n) = qu(n,k)
            xrain(n,k) = xrain(n,k+1)*(mfu(n,k+1)-ndmfde(n)) * &
                            (d_one/max(cmfcmin,mfu(n,k)))
            luold(n) = lu(n,k)
          end do
          ! reset to environmental values if below departure level
          do concurrent ( n = n1:n2 )
            if ( k > idpl(n) ) then
              tu(n,k) = tf(n,k)
              qu(n,k) = qf(n,k)
              lu(n,k) = d_zero
              luold(n) = lu(n,k)
            end if
          end do
          ! Do corrections for moist ascent by adjusting t,q and l
          ! ------------------------------------------------------
          ik = k
          icall = 1
          if ( nlm > 0 ) then
            call moistadj(ik,xph,tu,qu,llflag,icall)
          end if
          do concurrent ( nll = 1:nlm )
            n = nlx(nll)
            if ( qu(n,k) /= qold(n) ) then
              lglac(n,k) = lu(n,k) * ((d_one-xalpha(tu(n,k)))- &
                                     (d_one-xalpha(tu(n,k+1))))
              tu(n,k) = tu(n,k) + wlhfocp*lglac(n,k)
            end if
          end do
          do concurrent ( nll = 1:nlm )
            n = nlx(nll)
            if ( qu(n,k) /= qold(n) ) then
              ilab(n,k) = 2
              lu(n,k) = lu(n,k) + qold(n) - qu(n,k)
              bc = tu(n,k)*(d_one+ep1*qu(n,k)-lu(n,k+1) - xrain(n,k+1))
              be = tf(n,k)*(d_one+ep1*qf(n,k))
              buo(n,k) = bc - be
              ! set flags in case of midlevel convection
              if ( ktype(n) == 3 .and. ilab(n,k+1) == 1 ) then
                if ( buo(n,k) > -d_half ) then
                  ldcum(n) = .true.
                  kctop(n) = k
                  kineu(n,k) = d_half
                else
                  ilab(n,k) = 0
                  mfu(n,k) = d_zero
                  lude(n,k) = d_zero
                  lu(n,k) = d_zero
                end if
              end if
              if ( ilab(n,k+1) == 2 ) then
                if ( buo(n,k) < d_zero .and. ilab(n,k+1) == 2 ) then
                  tf(n,k) = d_half*(t(n,k)+t(n,k-1))
                  qf(n,k) = d_half*(q(n,k)+q(n,k-1))
                  buo(n,k) = bc - tf(n,k)*(d_one+ep1*qf(n,k))
                end if
                buoc = (buo(n,k) / (tf(n,k)*(d_one+ep1*qf(n,k)))+buo(n,k+1) / &
                                   (tf(n,k+1)*(d_one+ep1*qf(n,k+1))))*d_half
                dkbuo = (geof(n,k)-geof(n,k+1))*facbuo*buoc
                ! mixing and "pressure" gradient term in upper
                ! troposphere
                if ( ndmfen(n) > d_zero ) then
                  dken = min(d_one,(d_one+cwdrag)*ndmfen(n) / &
                         max(cmfcmin,mfu(n,k+1)))
                else
                  dken = min(d_one,(d_one+cwdrag)*ndmfde(n) / &
                         max(cmfcmin,mfu(n,k+1)))
                end if
                kineu(n,k) = (kineu(n,k+1)*(d_one-dken)+dkbuo) / (d_one+dken)
                if ( buo(n,k) < d_zero .and. ilab(n,k+1) == 2 ) then
                  kedke = kineu(n,k)/max(1.0e-10_rkx,kineu(n,k+1))
                  kedke = max(d_zero,min(d_one,kedke))
                  mfun = sqrt(kedke)*mfu(n,k+1)
                  ndmfde(n) = max(ndmfde(n),mfu(n,k+1)-mfun)
                  lude(n,k) = lu(n,k+1)*ndmfde(n)
                  mfu(n,k) = mfu(n,k+1) + ndmfen(n) - ndmfde(n)
                end if
                if ( buo(n,k) >- 0.2_rkx .and. ilab(n,k+1) == 2 ) then
                  ikb = kcbot(n)
                  entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(n))
                  zoentr(n) = entrpen*(0.3_rkx-(min(d_one,q(n,k-1) /    &
                         qs(n,k-1))-d_one))*(geof(n,k-1)-geof(n,k)) * &
                         regrav*min(d_one,qs(n,k)/qs(n,ikb))**3
                  zoentr(n) = min(0.4_rkx,zoentr(n))*mfu(n,k)
                else
                  zoentr(n) = d_zero
                end if
               ! Erase values if below departure level
                if ( k > idpl(n) ) then
                  mfu(n,k) = mfu(n,k+1)
                  kineu(n,k) = d_half
                end if
                if ( kineu(n,k) > d_zero .and. mfu(n,k) > d_zero ) then
                  kctop(n) = k
                  llo1(n) = .true.
                else
                  ilab(n,k) = 0
                  mfu(n,k) = d_zero
                  kineu(n,k) = d_zero
                  ndmfde(n) = mfu(n,k+1)
                  lude(n,k) = lu(n,k+1)*ndmfde(n)
                end if
                ! store detrainment rates for updraught
                if ( mfu(n,k+1) > d_zero ) mfude_rate(n,k) = ndmfde(n)
              end if
            else if ( ktype(n) == 2 .and. qu(n,k) == qold(n) ) then
              ilab(n,k) = 0
              mfu(n,k) = d_zero
              kineu(n,k) = d_zero
              ndmfde(n) = mfu(n,k+1)
              lude(n,k) = lu(n,k+1)*ndmfde(n)
              mfude_rate(n,k) = ndmfde(n)
            end if
          end do
          ! Calculate precipitation rate by analytic integration
          ! of equation for l
          do concurrent ( n = n1:n2 )
            if ( llo1(n) ) then
              if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
                dnoprc = ccn(n,k)*(4.0_rkx/3.0_rkx)*mathpi * &
                         ((rcrit*1e-6_rkx)**3)*rhow
                if ( ldland(n) ) then
                  prcdgw = rprc_lnd*regrav
                else
                  prcdgw = rprc_ocn*regrav
                end if
              else
                if ( ldland(n) ) then
                  dnoprc = 5.e-4_rkx
                  prcdgw = rprc_lnd*regrav
                else
                  dnoprc = 3.e-4_rkx
                  prcdgw = rprc_ocn*regrav
                end if
              end if
              if ( lu(n,k) > dnoprc ) then
                zwu = min(15.0_rkx,sqrt(d_two*max(0.1_rkx,kineu(n,k+1))))
                prcon = prcdgw/(0.75_rkx*zwu)
                ! Parameters for bergeron-findeisen process (T < -5C)
                zdt = min(rtber-rtice,max(rtber-tu(n,k),d_zero))
                cbf = d_one + cprc2*sqrt(zdt)
                zco = prcon*cbf
                lcrit = dnoprc/cbf
                zdfi = geof(n,k) - geof(n,k+1)
                zc = (lu(n,k)-luold(n))
                arg = (lu(n,k)/lcrit)**2
                if ( arg < 25.0_rkx ) then
                  zd = zco*(d_one-exp(-arg))*zdfi
                else
                  zd = zco*zdfi
                end if
                zint = exp(-zd)
                lnew = luold(n)*zint + zc/zd*(d_one-zint)
                lnew = max(d_zero,min(lu(n,k),lnew))
                lnew = min(cldmax,lnew)
                precip(n) = max(d_zero,luold(n)+zc-lnew)
                dmfup(n,k) = precip(n)*mfu(n,k)
                xrain(n,k) = xrain(n,k) + precip(n)
                lu(n,k) = lnew
                !FAB save precip water (kg/kg) created from
                ! autoconverson in updraft
                if ( ichem == 1 ) cvrainout(n,k) = precip(n)
              end if
            end if
          end do
          do concurrent ( n = n1:n2 )
            if ( llo1(n) ) then
              if ( xrain(n,k) > d_zero ) then
                zvw = 21.18_rkx*xrain(n,k)**0.2_rkx
                zvi = cwifrac*zvw
                alfaw = xalpha(tu(n,k))
                zvv = alfaw*zvw + (d_one-alfaw)*zvi
                rold = xrain(n,k) - precip(n)
                zc = precip(n)
                zwu = min(15.0_rkx,sqrt(d_two*max(0.1_rkx,kineu(n,k))))
                zd = zvv/zwu
                zint = exp(-zd)
                rnew = rold*zint + zc/zd*(d_one-zint)
                rnew = max(d_zero,min(xrain(n,k),rnew))
                xrain(n,k) = rnew
              end if
            end if
          end do
          do concurrent ( nll = 1:nlm )
            n = nlx(nll)
            mful(n,k) = lu(n,k)*mfu(n,k)
            mfus(n,k) = (cpd*tu(n,k)+geof(n,k))*mfu(n,k)
            mfuq(n,k) = qu(n,k)*mfu(n,k)
          end do
        end if
      end do
      !----------------------
      ! 5. Final calculations
      ! ---------------------
      do concurrent ( n = n1:n2 )
        if ( kctop(n) == -1 ) ldcum(n) = .false.
        kcbot(n) = max(kcbot(n),kctop(n))
        if ( ldcum(n) ) then
          wmean(n) = max(1.e-2_rkx,wmean(n)/max(d_one,dpmean(n)))
          wmean(n) = sqrt(d_two*wmean(n))
        end if
      end do
    end block ascent
    !
    !  (C) Check cloud depth and change entrainment rate accordingly
    !      Calculate precipitation rate (for downdraft calculation)
    ! -----------------------------------------------------------------------
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) ) then
        ikb = kcbot(n)
        itopm2 = kctop(n)
        pbmpt = pf(n,ikb) - pf(n,itopm2)
        if ( ktype(n) == 1 .and. pbmpt <  rdepths ) ktype(n) = 2
        if ( ktype(n) == 2 .and. pbmpt >= rdepths ) ktype(n) = 1
        ictop0(n) = kctop(n)
      end if
      rfl(n) = dmfup(n,1)
    end do
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 2, nk
        rfl(n) = rfl(n) + dmfup(n,k)
      end do
    end do
    do concurrent ( n = n1:n2, k = 1:nk )
      mfd(n,k) = d_zero
      mfds(n,k) = d_zero
      mfdq(n,k) = d_zero
      dmfdp(n,k) = d_zero
      dpmel(n,k) = d_zero
    end do
    !----------------------------------
    ! 5. Cumulus downdraft calculations
    ! ---------------------------------
    if ( lmfdd ) then
      !-------------------
      ! (A) Determine lfs
      ! --=---------------
      !
      ! This routine calculates level of free sinking for cumulus downdrafts
      ! and specifies t,q,u and v values to produce lfs-values for cumulus
      ! downdrafts for massflux cumulus parameterization interface
      ! Check for negative buoyancy of air of equal parts of moist
      ! environmental air and cloud air.
      !
      lfs: block
        integer(ik4), dimension(np) :: ikhsmin
        real(rkx), dimension(np,nk) :: tenwb, qenwb
        real(rkx), dimension(np) :: cond, xph, hsmin
        logical, dimension(np) :: llo2
        integer(ik4) :: icall, ik, ike, is, k, n
        real(rkx) :: buo, hsk, mftop, qtest, ttest
        !-------------------------------------
        ! 1. Set default values for downdrafts
        ! ------------------------------------
        do concurrent ( n = n1:n2 )
          lddraf(n) = .false.
          idtop(n) = nk + 1
          ikhsmin(n) = nk + 1
          hsmin(n) = 1.e8_rkx
        end do
        if ( lmfdd ) then
          !------------------------------------------------------------------
          ! 2. Determine level of free sinking:
          !    downdrafts shall start at model level of minimum of saturation
          !    moist static energy or below respectively for every point and
          !    proceed as follows:
          !      (1) Determine level of minimum of hs
          !      (2) Determine wet bulb environmental t and q
          !      (3) Do mixing with cumulus cloud air
          !      (4) Check for negative buoyancy
          !      (5) If buoyancy>0 repeat (2) to (4) for next
          !          level below
          ! The assumption is that air of downdrafts is mixture of 50% cloud
          ! air + 50% environmental air at wet bulb temperature (i.e. which
          ! became saturated due to evaporation of rain and cloud water)
          ! ----------------------------------------------------------------
          do concurrent ( n = n1:n2 )
            !$acc loop seq
            do k = 3, nk - 2
              hsk = cpd*t(n,k) + geo(n,k) + mlwt(t(n,k))*qs(n,k)
              if ( hsk < hsmin(n) ) then
                hsmin(n) = hsk
                ikhsmin(n) = k
              end if
            end do
          end do
          ike = nk - 3
          do k = 3, ike
            ! -----------------------------------------------
            ! 2.1 Calculate wet-bulb temperature and moisture
            !     for environmental air
            ! -----------------------------------------------
            is = 0
#ifdef STDPAR
            do concurrent ( n = n1:n2 ) reduce(+:is)
#else
            do n = n1, n2
#endif
              tenwb(n,k) = tf(n,k)
              qenwb(n,k) = qf(n,k)
              xph(n) = pf(n,k)
              llo2(n) = ldcum(n) .and.                     &
                rfl(n) > d_zero .and. .not.lddraf(n) .and. &
                (k < kcbot(n) .and. k > kctop(n)) .and.    &
                k >= ikhsmin(n)
              if ( llo2(n) ) is = is + 1
            end do
            if ( is == 0 ) cycle
            ik = k
            icall = 2
            call moistadj(ik,xph,tenwb,qenwb,llo2,icall)
            !-----------------------------------------------
            ! 2.2 Do mixing of cumulus and environmental air
            !     and check for negative buoyancy.
            !     Then set values for downdraft at lfs.
            ! ----------------------------------------------
            do concurrent ( n = n1:n2 )
              if ( llo2(n) ) then
                ttest = d_half*(tu(n,k)+tenwb(n,k))
                qtest = d_half*(qu(n,k)+qenwb(n,k))
                buo = ttest*(d_one+ep1*qtest) - tf(n,k)*(d_one+ep1*qf(n,k))
                cond(n) = qf(n,k) - qenwb(n,k)
                mftop = -rmfdeps*mfub(n)
                if ( buo < d_zero .and. rfl(n) > 10.0_rkx*mftop*cond(n) ) then
                  idtop(n) = k
                  lddraf(n) = .true.
                  td(n,k) = ttest
                  qd(n,k) = qtest
                  mfd(n,k) = mftop
                  mfds(n,k) = mfd(n,k)*(cpd*td(n,k)+geof(n,k))
                  mfdq(n,k) = mfd(n,k)*qd(n,k)
                  dmfdp(n,k-1) = -d_half*mfd(n,k)*cond(n)
                  rfl(n) = rfl(n) + dmfdp(n,k-1)
                end if
              end if
            end do
          end do
        end if
      end block lfs
      !---------------------------------------
      ! (B) Determine downdraft t,q and fluxes
      ! --------------------------------------
      !
      ! This routine calculates cumulus downdraft descent to produce the
      ! vertical profiles for cumulus downdrafts (i.e. t,q,u and v and fluxes)
      ! It calculates moist descent for entraining/detraining plume by
      !   A) Moving air dry-adiabatically to next level below and
      !   B) Correcting for evaporation to obtain saturated state.
      !
      ddrafdsc: block
        real(rkx), dimension(np) :: dmfen, xdmfde, cond, oentr, buoy
        real(rkx), dimension(np) :: xph
        logical, dimension(np) :: llo2
        integer(ik4) :: icall, ik, is, itopde, k, n
        real(rkx) :: buo, buoyz, buoyv, xdmfdp, dz, entr, mfdqk,  &
                     mfdsk, qdde, qeen, rain, sdde, seen, zentr, &
                     facbuo, dkbuo, dken
        itopde = nk950
        facbuo = d_half/(d_one+d_half)
        !-------------------------------------------------------
        ! 1. Calculate moist descent for cumulus downdraft by
        !      (A) Calculating entrainment/detrainment rates,
        !          including organized entrainment dependent on
        !          negative buoyancy and assuming
        !          linear decrease of massflux in PBL
        !      (B) Doing moist descent - evaporative cooling
        !          and moistening
        !      (C) Checking for negative buoyancy and
        !          specifying final t,q,u,v and downward fluxes
        ! -----------------------------------------------------
        do concurrent ( n = n1:n2 )
          oentr(n) = d_zero
          buoy(n) = d_zero
          dmfen(n) = d_zero
          xdmfde(n) = d_zero
        end do
        do concurrent ( n = n1:n2, k = 1:nk )
          dmfde(n,k) = d_zero
          mfdde_rate(n,k) = d_zero
          kined(n,k) = d_zero
        end do
        do k = 3, nk
          is = 0
#ifdef STDPAR
          do concurrent ( n = n1:n2 ) reduce(+:is)
#else
          do n = n1, n2
#endif
            xph(n) = pf(n,k)
            llo2(n) = lddraf(n) .and. mfd(n,k-1) < d_zero
            if ( llo2(n) ) is = is + 1
          end do
          if ( is == 0 ) cycle
          do concurrent ( n = n1:n2 )
            if ( llo2(n) ) then
              entr = entrdd*mfd(n,k-1)*(geof(n,k-1)-geof(n,k))*regrav
              dmfen(n) = entr
              xdmfde(n) = entr
            end if
          end do
          if ( k > itopde ) then
            do concurrent ( n = n1:n2 )
              if ( llo2(n) ) then
                dmfen(n) = d_zero
                xdmfde(n) = mfd(n,itopde)*(pf(n,k)-pf(n,k-1)) / &
                  (pf(n,nk+1)-pf(n,itopde))
              end if
            end do
          end if
          if ( k <= itopde ) then
            do concurrent ( n = n1:n2 )
              if ( llo2(n) ) then
                dz = -(geof(n,k-1)-geof(n,k))*regrav
                zentr = oentr(n)*dz*mfd(n,k-1)
                dmfen(n) = dmfen(n) + zentr
                dmfen(n) = max(dmfen(n),0.3_rkx*mfd(n,k-1))
                dmfen(n) = max(dmfen(n), &
                  -0.75_rkx*mfu(n,k)-(mfd(n,k-1)-xdmfde(n)))
                dmfen(n) = min(dmfen(n),d_zero)
              end if
              dmfde(n,k) = dmfen(n) - xdmfde(n)
            end do
          end if
          do concurrent ( n = n1:n2 )
            if ( llo2(n) ) then
              mfd(n,k) = mfd(n,k-1) + dmfen(n) - xdmfde(n)
              seen = (cpd*tf(n,k-1)+geof(n,k-1))*dmfen(n)
              qeen = qf(n,k-1)*dmfen(n)
              sdde = (cpd*td(n,k-1)+geof(n,k-1))*xdmfde(n)
              qdde = qd(n,k-1)*xdmfde(n)
              mfdsk = mfds(n,k-1) + seen - sdde
              mfdqk = mfdq(n,k-1) + qeen - qdde
              qd(n,k) = mfdqk*(d_one/min(-cmfcmin,mfd(n,k)))
              td(n,k) = (mfdsk*(d_one/min(-cmfcmin,mfd(n,k)))-geof(n,k))*rcpd
              td(n,k) = min(400.0_rkx,td(n,k))
              td(n,k) = max(100.0_rkx,td(n,k))
              cond(n) = qd(n,k)
            end if
          end do
          ik = k
          icall = 2
          call moistadj(ik,xph,td,qd,llo2,icall)
          do concurrent ( n = n1:n2 )
            if ( llo2(n) ) then
              cond(n) = cond(n) - qd(n,k)
              buo = td(n,k)*(d_one+ep1*qd(n,k)) - tf(n,k)*(d_one+ep1*qf(n,k))
              if ( rfl(n) > d_zero .and. mfu(n,k) > d_zero ) then
                rain = rfl(n)/mfu(n,k)
                buo = buo - td(n,k)*rain
              end if
              if ( buo >= d_zero .or. rfl(n)<=(mfd(n,k)*cond(n)) ) then
                mfd(n,k) = d_zero
                buo = d_zero
              end if
              mfds(n,k) = (cpd*td(n,k)+geof(n,k))*mfd(n,k)
              mfdq(n,k) = qd(n,k)*mfd(n,k)
              xdmfdp = -mfd(n,k)*cond(n)
              dmfdp(n,k-1) = xdmfdp
              rfl(n) = rfl(n) + xdmfdp
              ! Compute organized entrainment for use at next level
              buoyz = buo/tf(n,k)
              buoyv = buoyz
              buoyz = min(buoyz,d_zero)
              dz = -(geo(n,k-1)-geo(n,k))
              buoy(n) = buoy(n) + buoyz*dz
              oentr(n) = egrav*buoyz*d_half/(d_one+buoy(n))
              ! Store downdraught detrainment rates
              mfdde_rate(n,k) = -xdmfde(n)
              ! Compute kinetic energy
              dkbuo = dz*buoyv*facbuo
              if ( dmfen(n) < d_zero ) then
                dken = min(d_one,(d_one+cwdrag)*dmfen(n) / &
                       min(-cmfcmin,mfd(n,k-1)))
              else
                dken = min(d_one,(d_one+cwdrag)*xdmfde(n) / &
                       min(-cmfcmin,mfd(n,k-1)))
              end if
              kined(n,k) = max(d_zero, &
                (kined(n,k-1)*(d_one-dken)+dkbuo)/(d_one+dken))
            end if
          end do
        end do
      end block ddrafdsc
    end if
    !-----------
    ! 6. Closure
    ! ----------
    !   Recalculate cloud base massflux from a cape closure for deep
    !   convection (ktype=1) and by PBL equilibrum taking downdrafts
    !   into account for shallow convection (ktype=2)
    ! --------------------------------------------------------------
    ! Deep convection
    do concurrent ( n = n1:n2 )
      rheat(n) = d_zero
      xcape(n) = d_zero
      mfub1(n) = mfub(n)
    end do
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 1, nk
        llo1 = ldcum(n) .and. ktype(n) == 1
        if ( llo1 .and. k <= kcbot(n) .and. k > kctop(n) ) then
          ikb = kcbot(n)
          xro = pf(n,k)/(rgas*tf(n,k)*(d_one+ep1*qf(n,k)))
          dz = (geof(n,k-1)-geof(n,k))
          rheat(n) = rheat(n) + ((t(n,k-1)-t(n,k)+dz*rcpd) / &
                      tf(n,k)+ep1*(q(n,k-1)-q(n,k))) *       &
                      (egrav*(mfu(n,k)+mfd(n,k)))/xro
          xcape(n) = xcape(n) + ((tu(n,k)-tf(n,k))/tf(n,k) + &
                      ep1*(qu(n,k)-qf(n,k))-lu(n,k))*dz
        end if
      end do
    end do
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) .and. ktype(n) == 1 ) then
        ikb = kcbot(n)
        ik = kctop(n)
        xcape(n) = max(d_zero,min(xcape(n),5000.0_rkx))
        rheat(n) = max(1.e-4_rkx,rheat(n))
        tau(n) = (geof(n,ik)-geof(n,ikb)) / &
                   ((d_two+min(15.0_rkx,wmean(n)))*egrav)*rtau
        tau(n) = max(dtc,min(10800.0_rkx,tau(n)))
        tau(n) = max(720.0_rkx,tau(n))
        mfub1(n) = (xcape(n)*mfub(n))/(rheat(n)*tau(n))
        mfub1(n) = max(mfub1(n),0.001_rkx)
        mfmax = (pf(n,ikb)-pf(n,ikb-1))*cons2*rmflic + rmflia
        mfub1(n) = min(mfub1(n),mfmax)
      end if
    end do
    ! Shallow convection and mid_level
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) .and. (ktype(n) == 2 .or. ktype(n) == 3) ) then
        ikb = kcbot(n)
        if ( mfd(n,ikb) < d_zero ) then
          eps = -mfd(n,ikb)/max(mfub(n),1.0e-10_rkx)
        else
          eps = 0.0_rkx
        end if
        qumqe = qu(n,ikb) + lu(n,ikb) - eps*qd(n,ikb) - (d_one-eps)*qf(n,ikb)
        dqmin = max(0.01_rkx*qf(n,ikb),1.0e-10_rkx)
        ! Maximum permisable value of ud base mass flux
        mfmax = (pf(n,ikb)-pf(n,ikb-1))*cons2*rmflic + rmflia
        ! Shallow convection
        if ( ktype(n) == 2 ) then
          dh = cpd*(tu(n,ikb)-eps*td(n,ikb) - &
                   (d_one-eps)*tf(n,ikb)) + wlhv*qumqe
          dh = egrav*max(dh,p00*dqmin)
          if ( dhpbl(n) > d_zero ) then
            mfub1(n) = dhpbl(n)/dh
          else
            mfub1(n) = mfub(n)
          end if
          mfub1(n) = min(mfub1(n),mfmax)
          if ( lmfwstar ) mfub1(n) = mf_shal(n)
        end if
        ! Mid-level convection
        if ( ktype(n) == 3 ) then
          mfub1(n) = mfub(n)*(d_one+eps)
          mfub1(n) = min(mfub1(n),mfmax)
        end if
      end if
    end do
    ! Rescale DD fluxes if deep and shallow convection
    do concurrent ( n = n1:n2, k = 1:nk )
      if ( lddraf(n) .and. (ktype(n) == 1 .or. ktype(n) == 2) ) then
        fac = mfub1(n)/max(mfub(n),1.0e-10_rkx)
        mfd(n,k) = mfd(n,k)*fac
        mfds(n,k) = mfds(n,k)*fac
        mfdq(n,k) = mfdq(n,k)*fac
        dmfdp(n,k) = dmfdp(n,k)*fac
        ! also rescale detrainment flux for ERA pp
        mfdde_rate(n,k) = mfdde_rate(n,k)*fac
      end if
    end do
    !--------------------------
    ! 6.1 Final closure=scaling
    ! -------------------------
    do concurrent ( n = n1:n2 )
      if ( ldcum(n) ) mfs(n) = mfub1(n)/max(cmfcmin,mfub(n))
    end do
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 2, nk
        if ( ldcum(n) .and. k >= kctop(n)-1 ) then
          ikb = kcbot(n)
          if ( k>ikb ) then
            dz = ((pf(n,nk+1)-pf(n,k))/(pf(n,nk+1)-pf(n,ikb)))
            mfu(n,k) = mfu(n,ikb)*dz
          end if
          mfmax = (pf(n,k)-pf(n,k-1))*cons2*rmflic + rmflia
          if ( mfu(n,k)*mfs(n) > mfmax ) then
            mfs(n) = min(mfs(n),mfmax/mfu(n,k))
          end if
        end if
      end do
    end do
    do concurrent ( n = n1:n2, k = 2:nk )
      if ( ldcum(n) .and. k <= kcbot(n) .and. k >= kctop(n)-1 ) then
        mfu(n,k) = mfu(n,k)*mfs(n)
        mfus(n,k) = mfus(n,k)*mfs(n)
        mfuq(n,k) = mfuq(n,k)*mfs(n)
        mful(n,k) = mful(n,k)*mfs(n)
        dmfup(n,k) = dmfup(n,k)*mfs(n)
        dmfen(n,k) = dmfen(n,k)*mfs(n)
        lude(n,k) = lude(n,k)*mfs(n)
        mfude_rate(n,k) = mfude_rate(n,k)*mfs(n)
      end if
    end do
    !--------------------------------------------------------
    ! 6.2 In case that either deep or shallow is switched off
    !      reset ldcum to false-> fluxes set to zero
    ! -------------------------------------------------------
    ! exclude pathological ktype=2 kcbot=kctop=nk-1
    do concurrent ( n = n1:n2 )
      if ( ktype(n) == 2 .and.        &
           kcbot(n) == kctop(n) .and. &
           kcbot(n) >= nk-1 ) then
        ldcum(n) = .false.
        ktype(n) = 0
      end if
    end do
    if ( .not. lmfscv .or. .not. lmfpen ) then
      do concurrent ( n = n1:n2 )
        llo2(n) = .false.
        if ( (.not. lmfscv .and. ktype(n) == 2) .or. &
             (.not. lmfpen .and. ktype(n) == 1) ) then
          llo2(n) = .true.
          ldcum(n) = .false.
        end if
      end do
    end if
    !-------------------------------------
    ! 7. Determine final convective fluxes
    ! ------------------------------------
    !- set DD mass fluxes to zero above cloud top
    ! (because of inconsistency with second updraught)
    do concurrent ( n = n1:n2 )
      if ( lddraf(n) .and. idtop(n) <= kctop(n) ) then
        idtop(n) = kctop(n) + 1
      end if
    end do
    do concurrent ( n = n1:n2, k = 2:nk )
      if ( lddraf(n) ) then
        if ( k < idtop(n) ) then
          mfd(n,k) = d_zero
          mfds(n,k) = d_zero
          mfdq(n,k) = d_zero
          mfdde_rate(n,k) = d_zero
          dmfdp(n,k) = d_zero
        else if ( k == idtop(n) ) then
          mfdde_rate(n,k) = d_zero
        end if
      end if
    end do
    !
    ! This routine does the final calculation of convective fluxes in the
    ! cloud layer and in the subcloud layer
    !
    cfluxes: block
      real(rkx), dimension(np) :: rhebc
      integer(ik4) :: ik, ikb, k, n
      integer(ik4), dimension(np) :: idbas
      logical :: llddraf
      real(rkx) :: alfaw, cons1, cons1a, cons2, denom, drfl, &
                   drfl1, fac, pdr, pds, xrfl, xrfln, rmin,  &
                   rnew, snmlt, zp, rcpecons, rcucov
      !-------------------
      ! 0. Setup constants
      ! ------------------
      cons1a = cpd/(wlhf*egrav*rtaumel)
      cons2 = rmfcfl/(egrav*dtc)
      !-------------------------------------
      ! 1. Determine final convective fluxes
      ! ------------------------------------
      do concurrent ( n = n1:n2 )
        rain(n) = d_zero
        if ( .not.ldcum(n) .or. idtop(n) < kctop(n) ) lddraf(n) = .false.
        if ( .not.ldcum(n) ) ktype(n) = 0
        idbas(n) = nk
        if ( ldland(n) ) then
          rhebc(n) = rhebc_lnd
        else
          rhebc(n) = rhebc_ocn
        end if
      end do
      ! To get identical results force itopm2 to 2
      itopm2 = 2
      do k = itopm2, nk
        ikb = min(k+1,nk)
        do concurrent ( n = n1:n2 )
          mflxr(n,k) = d_zero
          mflxs(n,k) = d_zero
          dpmel(n,k) = d_zero
          if ( ldcum(n) .and. k >= kctop(n) ) then
            mfus(n,k) = mfus(n,k) - mfu(n,k)*(cpd*tf(n,k)+geof(n,k))
            mfuq(n,k) = mfuq(n,k) - mfu(n,k)*qf(n,k)
            lglac(n,k) = mfu(n,k)*lglac(n,k)
            llddraf = lddraf(n) .and. k >= idtop(n)
            if ( llddraf ) then
              mfds(n,k) = mfds(n,k) - mfd(n,k)*(cpd*tf(n,k)+geof(n,k))
              mfdq(n,k) = mfdq(n,k) - mfd(n,k)*qf(n,k)
            else
              mfd(n,k) = d_zero
              mfds(n,k) = d_zero
              mfdq(n,k) = d_zero
              dmfdp(n,k-1) = d_zero
            end if
            if ( llddraf .and.           &
                 mfd(n,k) < d_zero .and. &
                 abs(mfd(n,ikb)) < almostzero ) then
              idbas(n) = k
            end if
          else
            mfu(n,k) = d_zero
            mfd(n,k) = d_zero
            mfus(n,k) = d_zero
            mfds(n,k) = d_zero
            mfuq(n,k) = d_zero
            mfdq(n,k) = d_zero
            mful(n,k) = d_zero
            lglac(n,k) = d_zero
            dmfup(n,k-1) = d_zero
            dmfdp(n,k-1) = d_zero
            lude(n,k-1) = d_zero
          end if
        end do
      end do
      do concurrent( n = n1:n2 )
        mflxr(n,nk+1) = d_zero
        mflxs(n,nk+1) = d_zero
      end do
      !----------------------------------
      ! 1.1 Scale fluxes below cloud base
      !      linear dcrease
      ! ---------------------------------
      do n = n1, n2
        if ( ldcum(n) ) then
          ikb = kcbot(n)
          ik = ikb + 1
          zp = ((pf(n,nk+1)-pf(n,ik))/(pf(n,nk+1)-pf(n,ikb)))
          if ( ktype(n) == 3 ) zp = zp*zp
          mfu(n,ik) = mfu(n,ikb)*zp
          ! passing a scalar to avoid a compiler bug with OpenACC/Stdpar
          tfl=tf(n,ikb)
          mfus(n,ik) = (mfus(n,ikb) - mlwt(tfl) * mful(n,ikb))*zp
          !mfus(n,ik) = (mfus(n,ikb) - mlwt(tf(n,ikb)) * mful(n,ikb))*zp
          mfuq(n,ik) = (mfuq(n,ikb)+mful(n,ikb))*zp
          mful(n,ik) = d_zero
        end if
      end do
      do k = itopm2, nk
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) .and. k > kcbot(n)+1 ) then
            ikb = kcbot(n) + 1
            zp = ((pf(n,nk+1)-pf(n,k))/(pf(n,nk+1)-pf(n,ikb)))
            if ( ktype(n) == 3 ) zp = zp*zp
            mfu(n,k) = mfu(n,ikb)*zp
            mfus(n,k) = mfus(n,ikb)*zp
            mfuq(n,k) = mfuq(n,ikb)*zp
            mful(n,k) = d_zero
          end if
          ik = idbas(n)
          llddraf = lddraf(n) .and. k > ik .and. ik < nk
          if ( llddraf .and. ik == kcbot(n)+1 ) then
            zp = ((pf(n,nk+1)-pf(n,k))/(pf(n,nk+1)-pf(n,ik)))
            if ( ktype(n) == 3 ) zp = zp*zp
            mfd(n,k) = mfd(n,ik)*zp
            mfds(n,k) = mfds(n,ik)*zp
            mfdq(n,k) = mfdq(n,ik)*zp
            mfdde_rate(n,k) = -(mfd(n,k-1)-mfd(n,k))
          else if ( llddraf .and. ik /= kcbot(n)+1 .and. k == ik+1 ) then
            mfdde_rate(n,k) = -(mfd(n,k-1)-mfd(n,k))
          end if
        end do
      end do
      !-----------------------------------
      ! 2. Calculate rain/snow fall rates
      !    Calculate melting of snow
      !    Calculate evaporation of precip
      ! ----------------------------------
      do k = itopm2, nk
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) .and. k >= kctop(n)-1 ) then
            rain(n) = rain(n) + dmfup(n,k)
            if ( mflxs(n,k) > d_zero .and. t(n,k) > tzero ) then
              cons1 = cons1a*(d_one+d_half*(t(n,k)-tzero))
              fac = cons1*(pf(n,k+1)-pf(n,k))
              snmlt = min(mflxs(n,k),fac*(t(n,k)-tzero))
              dpmel(n,k) = snmlt
              qs(n,k) = fesat(t(n,k)-snmlt/fac)/ph(n,k)
            end if
            alfaw = xalpha(t(n,k))
            ! no liquid precipitation above melting level
            if ( t(n,k) < tzero .and. alfaw > d_zero ) then
              lglac(n,k) = lglac(n,k)+alfaw*(dmfup(n,k)+dmfdp(n,k))
              alfaw = d_zero
            end if
            mflxr(n,k+1) = mflxr(n,k) + alfaw*(dmfup(n,k)+dmfdp(n,k))+dpmel(n,k)
            mflxs(n,k+1) = mflxs(n,k) + &
                      (d_one-alfaw)*(dmfup(n,k)+dmfdp(n,k)) - dpmel(n,k)
            if ( mflxr(n,k+1)+mflxs(n,k+1) < d_zero ) then
              dmfdp(n,k) = -(mflxr(n,k)+mflxs(n,k)+dmfup(n,k))
              mflxr(n,k+1) = d_zero
              mflxs(n,k+1) = d_zero
              dpmel(n,k) = d_zero
            else if ( mflxr(n,k+1) < d_zero ) then
              mflxs(n,k+1) = mflxs(n,k+1) + mflxr(n,k+1)
              mflxr(n,k+1) = d_zero
            else if ( mflxs(n,k+1) < d_zero ) then
              mflxr(n,k+1) = mflxr(n,k+1) + mflxs(n,k+1)
              mflxs(n,k+1) = d_zero
            end if
          end if
        end do
      end do
      do k = itopm2, nk
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) .and. k >= kcbot(n) ) then
            xrfl = mflxr(n,k) + mflxs(n,k)
            if ( xrfl > almostzero ) then
              rcpecons = merge(rcpec_lnd,rcpec_ocn,ldland(n))
              rcucov = merge(rcuc_lnd,rcuc_ocn,ldland(n))
              drfl1 = rcpecons * max(d_zero,qs(n,k)-q(n,k))*rcucov * &
                      (sqrt(pf(n,k)/pf(n,nk+1)) / &
                      5.09e-3_rkx*xrfl/rcucov)**0.5777_rkx*(pf(n,k+1)-pf(n,k))
              rnew = xrfl - drfl1
              rmin = xrfl - rcucov*max(d_zero,rhebc(n)*qs(n,k) - &
                      q(n,k))*cons2*(pf(n,k+1)-pf(n,k))
              rnew = max(rnew,rmin)
              xrfln = max(rnew,d_zero)
              drfl = min(d_zero,xrfln-xrfl)
              alfaw = xalpha(t(n,k))
              if ( t(n,k) < tzero ) alfaw = d_zero
              pdr = alfaw*dmfdp(n,k)
              pds = (d_one-alfaw)*dmfdp(n,k)
              denom = d_one/max(almostzero,mflxr(n,k)+mflxs(n,k))
              mflxr(n,k+1) = mflxr(n,k) + pdr + dpmel(n,k) + &
                             drfl*mflxr(n,k)*denom
              mflxs(n,k+1) = mflxs(n,k) + pds - dpmel(n,k) + &
                             drfl*mflxs(n,k)*denom
              dmfup(n,k) = dmfup(n,k) + drfl
              if ( mflxr(n,k+1)+mflxs(n,k+1) < d_zero ) then
                dmfup(n,k) = dmfup(n,k)-(mflxr(n,k+1)+mflxs(n,k+1))
                mflxr(n,k+1) = d_zero
                mflxs(n,k+1) = d_zero
                dpmel(n,k) = d_zero
              else if ( mflxr(n,k+1) < d_zero ) then
                mflxs(n,k+1) = mflxs(n,k+1) + mflxr(n,k+1)
                mflxr(n,k+1) = d_zero
              else if ( mflxs(n,k+1) < d_zero ) then
                mflxr(n,k+1) = mflxr(n,k+1) + mflxs(n,k+1)
                mflxs(n,k+1) = d_zero
              end if
            else
              mflxr(n,k+1) = d_zero
              mflxs(n,k+1) = d_zero
              dmfdp(n,k) = d_zero
              dpmel(n,k) = d_zero
            end if
          end if
        end do
      end do
    end block cfluxes
    !- rescale DD fluxes if total mass flux becomes negative
    !- correct DD detrainment rates if entrainment becomes negative
    !- correct UD detrainment rates if entrainment becomes negative
    !- conservation correction for precip
    !$acc kernels
    mfs(:) = d_one
    !$acc end kernels
    ! change for stability
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 2, nk
        if ( lddraf(n) .and. k >= idtop(n)-1 ) then
          mfmax = mfu(n,k)*0.98_rkx
          if ( mfd(n,k)+mfmax+1.e-15_rkx < d_zero ) then
            mfs(n) = min(mfs(n),-mfmax/mfd(n,k))
          end if
        end if
      end do
    end do
    !$acc kernels
    mfuub(:) = d_zero
    !$acc end kernels
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 2, nk
        if ( mfs(n) < d_one .and. k >= idtop(n)-1 ) then
          mfd(n,k) = mfd(n,k)*mfs(n)
          mfds(n,k) = mfds(n,k)*mfs(n)
          mfdq(n,k) = mfdq(n,k)*mfs(n)
          mfdde_rate(n,k) = mfdde_rate(n,k)*mfs(n)
          mfuub(n) = mfuub(n) - (d_one-mfs(n))*dmfdp(n,k)
          mflxr(n,k+1) = mflxr(n,k+1) + mfuub(n)
          dmfdp(n,k) = dmfdp(n,k)*mfs(n)
        end if
      end do
    end do
    do concurrent ( n = n1:n2, k = 2:nk-1 )
      if ( lddraf(n) .and. k >= idtop(n)-1 ) then
        erate = -mfd(n,k) + mfd(n,k-1) + mfdde_rate(n,k)
        if ( erate < d_zero ) then
          mfdde_rate(n,k) = mfdde_rate(n,k) - erate
        end if
      end if
      if ( ldcum(n) .and. k >= kctop(n)-1 ) then
        erate = mfu(n,k) - mfu(n,k+1) + mfude_rate(n,k)
        if ( erate < d_zero ) then
          mfude_rate(n,k) = mfude_rate(n,k) - erate
        end if
        dmfup(n,k) = mflxr(n,k+1) + mflxs(n,k+1) - mflxr(n,k) - mflxs(n,k)
        dmfdp(n,k) = d_zero
      end if
    end do
    ! Avoid negative humidities at ddraught top
    do concurrent ( n = n1:n2 )
      if ( lddraf(n) ) then
        k = idtop(n)
        ik = min(k+1,nk)
        if ( mfdq(n,k) < 0.3_rkx*mfdq(n,ik) ) then
          if ( abs(rmfsoltq) < almostzero ) then
            mfdq(n,k) = 0.3_rkx*mfdq(n,ik)
          else
            mfd(n,k) = 0.3_rkx*mfd(n,ik)
          end if
        end if
      end if
    end do
    ! Avoid negative humidities near cloud top because gradient of precip flux
    ! and detrainment / liquid water flux too large
    do concurrent ( n = n1:n2, k = 2:nk )
      if ( ldcum(n) .and. k >= kctop(n)-1 .and. k < kcbot(n) ) then
        dz = dtc*egrav/(pf(n,k+1)-pf(n,k))
        mfa = mfuq(n,k+1) + mfdq(n,k+1) - mfuq(n,k) - mfdq(n,k) + &
              mful(n,k+1) - mful(n,k) + dmfup(n,k)
        mfa = (mfa-lude(n,k))*dz
        if ( q(n,k)+mfa < d_zero ) then
          lude(n,k) = lude(n,k) + d_two*(q(n,k)+mfa)/dz
        end if
        if ( lude(n,k) < d_zero ) lude(n,k) = d_zero
      end if
      if ( .not. ldcum(n) ) mfude_rate(n,k) = d_zero
      if ( abs(mfd(n,k-1)) < almostzero ) mfdde_rate(n,k) = d_zero
    end do
    if ( llconscheck ) then
      if ( lmftrac .and. ntrac > 0 ) then
        allocate (xtenc(np,nk,ntrac))
        allocate (xsumc(np,4+ntrac))
        do concurrent ( n = n1:n2, k = 2:nk, nt = 1:ntrac )
          if ( ldcum(n) ) xtenc(n,k,nt) = tenc(n,k,nt)
        end do
      else
        allocate (xsumc(np,4))
      end if
    end if
    !---------------------------------
    ! 8. Update tendencies for t and q
    ! --------------------------------
    if ( rmfsoltq > d_zero ) then
      !   derive draught properties for implicit
      do concurrent ( n = n1:n2, k = nk:2:-1 )
        if ( ldcum(n) ) then
          if ( k > kcbot(n) ) then
            mfa = d_one/max(1.e-15_rkx,mfu(n,k))
            qu(n,k) = qf(n,k) + mfuq(n,k)*mfa
            tu(n,k) = tf(n,k) + mfus(n,k)*mfa*rcpd
            mfus(n,k) = mfu(n,k)*(cpd*tu(n,k)+geof(n,k))
            mfuq(n,k) = mfu(n,k)*qu(n,k)
            if ( lddraf(n) ) then
              mfa = d_one/min(-1.e-15_rkx,mfd(n,k))
              qd(n,k) = qf(n,k) + mfdq(n,k)*mfa
              td(n,k) = tf(n,k) + mfds(n,k)*mfa*rcpd
              mfdq(n,k) = mfd(n,k)*qd(n,k)
              mfds(n,k) = mfd(n,k)*(cpd*td(n,k)+geof(n,k))
            end if
          else if ( k <= kcbot(n) .and. k >= kctop(n) ) then
            mfus(n,k) = mfu(n,k)*(cpd*tu(n,k)+geof(n,k))
            mfuq(n,k) = mfu(n,k)*qu(n,k)
            mfds(n,k) = mfd(n,k)*(cpd*td(n,k)+geof(n,k))
            mfdq(n,k) = mfd(n,k)*qd(n,k)
          end if
        end if
      end do
    end if
    !
    ! Updates t and q tendencies, precipitation rates, global diagnostics
    !
    dtdqc: block
      logical :: lltest
      integer(ik4) :: k, ik, n
      real(rkx) :: ximp, xalv, zp, gq, gs, gh, zs, zq
      real(rkx), dimension(np,nk) :: xmfus, xmfuq, xmfds, xmfdq
      real(rkx), dimension(:,:), allocatable :: dtdt, dqdt, dp
      real(rkx), dimension(:,:), allocatable :: bb, r1, r2
      logical, dimension(:,:), allocatable :: llcumbas
      !-----------------------------
      ! 1. Setup and initializations
      ! ----------------------------
      ximp = d_one - rmfsoltq
      allocate (dtdt(np,nk))
      allocate (dqdt(np,nk))
      allocate (dp(np,nk))
      ! Zero detrained liquid water if diagnostic cloud scheme to be used
      ! This means that detrained liquid water will be evaporated in the
      ! cloud environment and not fed directly into a cloud liquid water
      ! variable
      lltest = .not. lepcld
      if ( lltest ) then
        do concurrent ( n = n1:n2, k = 1:nk )
          lude(n,k) = d_zero
        end do
      end if
      do concurrent ( n = n1:n2, k = 1:nk )
        if ( ldcum(n) ) then
          dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
          xmfus(n,k) = mfus(n,k)
          xmfds(n,k) = mfds(n,k)
          xmfuq(n,k) = mfuq(n,k)
          xmfdq(n,k) = mfdq(n,k)
        end if
      end do
      if ( rmfsoltq > d_zero ) then
        !-------------------------------------------
        ! 2. Recompute convective fluxes if implicit
        ! ------------------------------------------
        do concurrent ( n = n1:n2, k = itopm2:nk )
          ik = k - 1
          if ( ldcum(n) .and. k >= kctop(n)-1 ) then
            ! Compute interpolating coefficients GS and GQ
            ! for half-level values
            gq = (qf(n,k)-q(n,ik))/qs(n,k)
            gh = cpd*t(n,k) + geo(n,k)
            gs = (cpd*(tf(n,k)-t(n,ik)) + geof(n,k)-geo(n,ik))/gh
            ! Half-level environmental values for ZS and ZQ
            zs = cpd*(ximp*t(n,ik)+gs*t(n,k)) + geo(n,ik) + gs*geo(n,k)
            zq = ximp*q(n,ik) + gq*qs(n,k)
            xmfus(n,k) = mfus(n,k) - mfu(n,k)*zs
            xmfuq(n,k) = mfuq(n,k) - mfu(n,k)*zq
            if ( lddraf(n) .and. k >= idtop(n) ) then
              xmfds(n,k) = mfds(n,k) - mfd(n,k)*zs
              xmfdq(n,k) = mfdq(n,k) - mfd(n,k)*zq
            end if
          end if
        end do
      end if
      !----------------------
      ! 3. Compute tendencies
      ! ---------------------
      do k = itopm2, nk
        if ( k < nk ) then
          do concurrent ( n = n1:n2 )
            xalv = mlwt(t(n,k))
            if ( ldcum(n) ) then
              dtdt(n,k) = dp(n,k)*rcpd *                      &
                (xmfus(n,k+1)-xmfus(n,k)+xmfds(n,k+1) -       &
                 xmfds(n,k)+wlhf*lglac(n,k)-wlhf*dpmel(n,k) - &
                 xalv*(mful(n,k+1)-mful(n,k)-lude(n,k)-dmfup(n,k)))
              dqdt(n,k) = dp(n,k)*(xmfuq(n,k+1) -                &
                xmfuq(n,k)+xmfdq(n,k+1)-xmfdq(n,k)+mful(n,k+1) - &
                mful(n,k)-lude(n,k)-dmfup(n,k))
            end if
          end do
        else
          do concurrent ( n = n1:n2 )
            xalv = mlwt(t(n,k))
            if ( ldcum(n) ) then
              dtdt(n,k) = -dp(n,k)*rcpd *                &
                (xmfus(n,k)+xmfds(n,k)+wlhf*dpmel(n,k) - &
                 xalv*(mful(n,k)+dmfup(n,k)))
              dqdt(n,k) = -dp(n,k)*(xmfuq(n,k) + &
                xmfdq(n,k)+(mful(n,k)+dmfup(n,k)))
            end if
          end do
        end if
      end do
      if ( abs(rmfsoltq) < almostzero ) then
        !----------------------
        ! 3.1 Update tendencies
        ! ---------------------
        do concurrent ( n = n1:n2, k = itopm2:nk )
          if ( ldcum(n) ) then
            tent(n,k) = tent(n,k) + dtdt(n,k)
            tenq(n,k) = tenq(n,k) + dqdt(n,k)
          end if
        end do
      else
        !----------------------
        ! 3.2 Implicit solution
        ! ---------------------
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
        ! DTDT and DQDT correspond to the RHS ("constants") of the equation
        ! The solution is in R1 and R2
        allocate (bb(np,nk), source=1.0_rkx)
        allocate (r1(np,nk), source=0.0_rkx)
        allocate (r2(np,nk), source=0.0_rkx)
        allocate (llcumbas(np,nk), source=.false.)
        do concurrent ( n = n1:n2, k = 1:nk )
          xmfus(n,k) = d_zero
        end do
        ! Fill vectors A, B and RHS
        do concurrent ( n = n1:n2, k = itopm2:nk )
          ik = k + 1
          llcumbas(n,k) = ldcum(n) .and. k >= kctop(n) - 1
          if ( llcumbas(n,k) ) then
            zp = rmfsoltq*dp(n,k)*dtc
            xmfus(n,k) = -zp*(mfu(n,k)+mfd(n,k))
            dtdt(n,k) = dtdt(n,k)*dtc + t(n,k)
            dqdt(n,k) = dqdt(n,k)*dtc + q(n,k)
            if ( k < nk ) then
              bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
            else
              bb(n,k) = d_one
            end if
          end if
        end do
        call solver(n1,n2,nk,kctop,llcumbas,xmfus,bb,dtdt,r1)
        call solver(n1,n2,nk,kctop,llcumbas,xmfus,bb,dqdt,r2)
        ! Compute tendencies
        do concurrent ( n = n1:n2, k = itopm2:nk )
          if ( llcumbas(n,k) ) then
            tent(n,k) = tent(n,k) + (r1(n,k)-t(n,k))/dtc
            tenq(n,k) = tenq(n,k) + (r2(n,k)-q(n,k))/dtc
          end if
        end do
        deallocate (llcumbas)
        deallocate (r2)
        deallocate (r1)
        deallocate (bb)
      end if
      deallocate (dp)
      deallocate (dqdt)
      deallocate (dtdt)
    end block dtdqc
    !-------------------------------------------------
    ! 9. Compute momentum in updraught and downdraught
    ! ------------------------------------------------
    if ( lmfdudv ) then
      do concurrent ( n = n1:n2 )
        !$acc loop seq
        do k = nk - 1, 2, -1
          ik = k + 1
          if ( ldcum(n) ) then
            if ( k == kcbot(n) .and. ktype(n) < 3 ) then
              ikb = idpl(n)
              uu(n,k) = u(n,ikb-1)
              vu(n,k) = v(n,ikb-1)
            else if ( k == kcbot(n) .and. ktype(n) == 3 ) then
              uu(n,k) = u(n,k-1)
              vu(n,k) = v(n,k-1)
            end if
            if ( k < kcbot(n) .and. k >= kctop(n) ) then
              fac = d_zero
              if ( ktype(n) == 1 .or. ktype(n) == 3 ) fac = d_two
              if ( ktype(n) == 1 .and. k <= kctop(n)+2 ) fac = d_three
              erate = mfu(n,k) - mfu(n,ik) + (d_one+fac)*mfude_rate(n,k)
              derate = (d_one+fac)*mfude_rate(n,k)
              mfa = d_one/max(cmfcmin,mfu(n,k))
              uu(n,k) = (uu(n,ik)*mfu(n,ik) + erate*u(n,k)-derate*uu(n,ik))*mfa
              vu(n,k) = (vu(n,ik)*mfu(n,ik) + erate*v(n,k)-derate*vu(n,ik))*mfa
            end if
          end if
        end do
      end do
      do concurrent ( n = n1:n2 )
        !$acc loop seq
        do k = 3, nk
          ik = k - 1
          if ( ldcum(n) ) then
            if ( k == idtop(n) ) then
              ud(n,k) = d_half*(uu(n,k)+u(n,ik))
              vd(n,k) = d_half*(vu(n,k)+v(n,ik))
            else if ( k > idtop(n) ) then
              erate = -mfd(n,k) + mfd(n,ik) + mfdde_rate(n,k)
              mfa = d_one/min(-cmfcmin,mfd(n,k))
              ud(n,k) = (ud(n,ik)*mfd(n,ik) - &
                        erate*u(n,ik)+mfdde_rate(n,k)*ud(n,ik))*mfa
              vd(n,k) = (vd(n,ik)*mfd(n,ik) - &
                        erate*v(n,ik)+mfdde_rate(n,k)*vd(n,ik))*mfa
            end if
          end if
        end do
      end do
      !----------------------------------
      ! 9.1 Update tendencies for u and v
      ! ---------------------------------
      !---------------------------------------------------------------
      ! For explicit/semi-implicit rescale massfluxes for stability in
      ! momentum
      !---------------------------------------------------------------
      !$acc kernels
      mfs(:) = d_one
      !$acc end kernels
      if ( rmfsoluv <= d_one ) then
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = 2, nk
            if ( ldcum(n) .and. k >= kctop(n)-1 ) then
              mfmax = (pf(n,k)-pf(n,k-1))*cons
              if ( mfu(n,k) > mfmax .and. k >= kctop(n) ) then
                mfs(n) = min(mfs(n),mfmax/mfu(n,k))
              end if
            end if
          end do
        end do
      end if
      do concurrent ( n = n1:n2, k = 1:nk )
        mfuus(n,k) = mfu(n,k)
        mfdus(n,k) = mfd(n,k)
        if ( ldcum(n) .and. k >= kctop(n)-1 ) then
          mfuus(n,k) = mfu(n,k)*mfs(n)
          mfdus(n,k) = mfd(n,k)*mfs(n)
        end if
      end do
      ! Recompute Draught properties below for Implicit
      ! based on linear flux profiles
      if ( rmfsoluv > d_zero ) then
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) ) then
            k = kcbot(n)
            ik = k - 1
            mfuub(n) = mfuus(n,k)*(uu(n,k)-u(n,ik))
            mfuvb(n) = mfuus(n,k)*(vu(n,k)-v(n,ik))
          end if
        end do
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = 2, nk
            ik = k - 1
            if ( ldcum(n) .and. k > kcbot(n) ) then
              ikb = kcbot(n)
              dz = ((pf(n,nk+1)-pf(n,k)) / (pf(n,nk+1)-pf(n,ikb)))
              if ( ktype(n) == 3 ) dz = dz*dz
              mfa = d_one/max(cmfcmin,mfuus(n,k))
              uu(n,k) = u(n,ik) + mfuub(n)*dz*mfa
              vu(n,k) = v(n,ik) + mfuvb(n)*dz*mfa
              mfdus(n,k) = mfdus(n,ikb)*dz
              ud(n,k) = u(n,ik) + ud(n,ikb) - u(n,ikb-1)
              vd(n,k) = v(n,ik) + vd(n,ikb) - v(n,ikb-1)
            end if
            ! add UV perturb to correct wind bias
            if ( ldcum(n) .and. k >= kctop(n) ) then
              uu(n,k) = uu(n,k) - ruvper*sign(d_one,uu(n,k))
              vu(n,k) = vu(n,k) - ruvper*sign(d_one,vu(n,k))
            end if
          end do
        end do
      end if
      !------------------------------------
      ! End
      ! Intermediate Solution for stability
      !------------------------------------
      call dudvx(mfuus,mfdus)
      if ( lmfuvdis ) then
        ! add KE dissipation
        do concurrent ( n = n1:n2 )
          xsum12(n) = d_zero
          xsum22(n) = d_zero
        end do
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = 1, nk
            uv2(n,k) = d_zero
            if ( ldcum(n) .and. k >= kctop(n)-1 ) then
              dz = (pf(n,k+1)-pf(n,k))
              duten = tenu(n,k) - xtenu(n,k)
              dvten = tenv(n,k) - xtenv(n,k)
              uv2(n,k) = sqrt(duten**2+dvten**2)
              xsum22(n) = xsum22(n) + uv2(n,k)*dz
              xsum12(n) = xsum12(n) - (u(n,k)*duten+v(n,k)*dvten)*dz
            end if
          end do
        end do
        do concurrent ( n = n1:n2, k = 1:nk )
          if ( ldcum(n) .and. k>=kctop(n)-1 ) then
            dz = (pf(n,k+1)-pf(n,k))
            tdis = rcpd*xsum12(n)*uv2(n,k)/max(1.e-15_rkx,xsum22(n))
            tent(n,k) = tent(n,k) + tdis
          end if
        end do
      end if
    end if
    !--------------------------------------------------------
    ! 10. In case that either deep or shallow is switched off
    !     need to set some variables a posteriori to zero
    !--------------------------------------------------------
    if ( .not. lmfscv .or. .not. lmfpen ) then
      do concurrent ( n = n1:n2, k = 2:nk )
        if ( llo2(n) .and. k >= kctop(n)-1 ) then
          tu(n,k) = t(n,k)
          qu(n,k) = q(n,k)
          lu(n,k) = d_zero
          mfude_rate(n,k) = d_zero
          mfdde_rate(n,k) = d_zero
        end if
      end do
      do concurrent ( n = n1:n2 )
        if ( llo2(n) ) then
          kctop(n) = nk - 1
          kcbot(n) = nk - 1
        end if
      end do
    end if
    !------------------------------
    ! 11. Chemical tracer transport
    ! -----------------------------
    if ( lmftrac .and. ntrac > 0 ) then
      ! transport switched off for mid-level convection
      do concurrent ( n = n1:n2 )
        if ( ldcum(n) .and. ktype(n) /= 3 .and. kcbot(n)-kctop(n) >= 1 ) then
          lldcum(n) = .true.
          llddraf3(n) = lddraf(n)
        else
          lldcum(n) = .false.
          llddraf3(n) = .false.
        end if
      end do
      ! check and correct mass fluxes for CFL criterium
      !$acc kernels
      mfs(:) = d_one
      !$acc end kernels
      if ( rmfsolct <= d_three ) then
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = 2, nk
            if ( lldcum(n) .and. k >= kctop(n) ) then
              mfmax = (pf(n,k)-pf(n,k-1))*0.8_rkx*cons
              if ( mfu(n,k) > mfmax ) then
                mfs(n) = min(mfs(n),mfmax/mfu(n,k))
              end if
            end if
          end do
        end do
      end if
      do concurrent ( n = n1:n2, k = 1:nk )
        if ( lldcum(n) .and. k >= kctop(n)-1 ) then
          mfuus(n,k) = mfu(n,k)*mfs(n)
          mfudr(n,k) = mfude_rate(n,k)*mfs(n)
        else
          mfuus(n,k) = d_zero
          mfudr(n,k) = d_zero
        end if
        if ( llddraf3(n) .and. k >= idtop(n)-1 ) then
          mfdus(n,k) = mfd(n,k)*mfs(n)
          mfddr(n,k) = mfdde_rate(n,k)*mfs(n)
        else
          mfdus(n,k) = d_zero
          mfddr(n,k) = d_zero
        end if
      end do
      if ( lmfsmooth ) then
        ! smoothing of mass fluxes (gradients) at top and bottom of draughts
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = 2, nk - 1
            if ( llddraf3(n) .and. mfdus(n,k) < d_zero .and. &
                 abs(mfdus(n,k+1)) < almostzero ) then
              erate = min(d_zero,mfdus(n,k)-d_half*mfdus(n,k-1))
              mfdus(n,k) = mfdus(n,k) - erate
              mfddr(n,k) = mfddr(n,k) - erate
              mfddr(n,k+1) = -mfdus(n,k)
            end if
            if ( lldcum(n) .and. k==kctop(n) ) then
              erate = max(d_zero,mfuus(n,k)-d_half*mfuus(n,k+1))
              mfuus(n,k) = mfuus(n,k) - erate
              mfudr(n,k) = mfudr(n,k) + erate
              mfudr(n,k-1) = mfuus(n,k)
            end if
          end do
        end do
        do concurrent ( n = n1:n2 )
          !$acc loop seq
          do k = nk - 1, 2, -1
            if ( lldcum(n) ) then
              if ( abs(mfudr(n,k)) < almostzero .and. &
                       mfudr(n,k-1) > d_zero ) then
                mfudr(n,k) = d_half*mfudr(n,k-1)
              end if
            end if
          end do
        end do
      end if
      call ctracer(lldcum,llddraf3,mfuus,mfdus,mfudr,mfddr)
    end if
    !----------------------------------------------------------
    ! 12. Put detrainment rates from mflx units in units mflx/m
    ! ---------------------------------------------------------
    do concurrent ( n = n1:n2, k = 2:nk )
      if ( ldcum(n) ) then
        xro = egrav/(geof(n,k)-geof(n,k+1)) ! 1/dz
        mfude_rate(n,k) = mfude_rate(n,k)*xro
        mfdde_rate(n,k) = mfdde_rate(n,k)*xro
        if ( k < kctop(n) ) then
          lu(n,k) = d_zero
          tu(n,k) = t(n,k)
          qu(n,k) = q(n,k)
        end if
      end if
    end do

    if ( llconscheck ) then
      !--------------------------------------
      ! 13. Conservation check AND correction
      ! -------------------------------------
      do n = n1, n2
        do k = 1, 4
          xsumc(n,k) = d_zero
        end do
      end do
      do n = n1, n2
        do k = nk, 2, -1
          xalv = mlwt(t(n,k))
          if ( ldcum(n) .and. k >= kctop(n)-1 ) then
            dz = (pf(n,k+1)-pf(n,k))*regrav
            xsumc(n,1) = xsumc(n,1) + &
                         (tenq(n,k)-xtenq(n,k))*dz + lude(n,k)
            xsumc(n,2) = xsumc(n,2) + &
                         cpd*(tent(n,k)-xtent(n,k))*dz - xalv*lude(n,k)
            xsumc(n,3) = xsumc(n,3) + (tenu(n,k)-xtenu(n,k))*dz
            xsumc(n,4) = xsumc(n,4) + (tenv(n,k)-xtenv(n,k))*dz
          end if
        end do
      end do
      if ( lmftrac .and. ntrac > 0 ) then
        do nt = 1, ntrac
          do n = n1, n2
            do k = nk, 2, -1
              if ( ldcum(n) .and. k >= kctop(n)-1 ) then
                dz = (pf(n,k+1)-pf(n,k))*regrav
                xsumc(n,4+nt) = xsumc(n,4+nt) + (tenc(n,k,nt)-xtenc(n,k,nt))*dz
              end if
            end do
          end do
        end do
      end if
      do n = n1, n2
        xalv = mlwt(t(n,nk))
        if ( ldcum(n) ) then
          sfl(n) = mflxr(n,nk+1) + mflxs(n,nk+1)
          write (61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')     &
            n, ' CONS q: ', -xsumc(n,1)*xalv, sfl(n)*xalv,   &
            ktype(n), ' CONS h: ', xsumc(n,2), ' CONS uv: ', &
            xsumc(n,3), xsumc(n,4)
          if ( lmftrac .and. ntrac > 0 ) then
            write (61,*) ' Conserv Error Tracers 1-', ntrac, ' :'
            do nt = 1, ntrac
              write (61,'(i4,e12.4)') nt, xsumc(n,4+nt)
            end do
#ifndef TESTME
            call fatal(__FILE__,__LINE__, &
                       'ERROR IN TRACER CONSERVATION')
#endif
          end if
          ikb = kctop(n)
          dz = (pf(n,nk+1)-pf(n,ikb-1))*regrav
          xsumc(n,1) = (xsumc(n,1)+sfl(n))/dz
          xsumc(n,2) = (xsumc(n,2)-xalv*sfl(n))/(dz*cpd)
        end if
      end do
      deallocate (xsumc)
      if ( lmftrac .and. ntrac > 0 ) deallocate (xtenc)
    end if
    !-------------------------------------------------------
    ! 14. Compute convective tendencies for liquid and solid
    !     cloud condensate
    ! ------------------------------------------------------
    do concurrent ( n = n1:n2, k = 1:nk )
      tenl(n,k) = lude(n,k)*egrav/(pf(n,k+1)-pf(n,k))
      teni(n,k) = (d_one-xalpha(t(n,k)))*tenl(n,k)
      tenl(n,k) = tenl(n,k) - teni(n,k)
      mflxr(n,k) = mflxr(n,k)*1.e-3_rkx
      mflxs(n,k) = mflxs(n,k)*1.e-3_rkx
    end do
    do concurrent ( n = n1:n2 )
      mflxr(n,nk+1) = mflxr(n,nk+1)*1.e-3_rkx
      mflxs(n,nk+1) = mflxs(n,nk+1)*1.e-3_rkx
    end do

    contains
    !
    ! Produce adjusted t,q and l values
    !
    subroutine moistadj(kk,sp,t,q,ldflag,jcall)
      implicit none
      integer(ik4), intent(in) :: kk ! actual level
      real(rkx), dimension(np), intent(in) :: sp ! Surface pressure Pa
      real(rkx), dimension(np,nk), intent(inout) :: t
      real(rkx), dimension(np,nk), intent(inout) :: q
      logical, dimension(np), intent(in) :: ldflag ! Active convection flag
      integer(ik4), intent(in) :: jcall ! Calling method
      integer(ik4) :: jl
      real(rkx) :: cond, cond1, cor, qs, rp
      real(rkx) :: zl, zi, zf, za
      !---------------------------------------------------------
      ! 1. Calculate condensation and adjust t and q accordingly
      ! --------------------------------------------------------
      if ( jcall == 1 ) then
        do concurrent ( jl = n1:n2 )
          if ( ldflag(jl) ) then
            rp = d_one/sp(jl)
            zl = d_one/(t(jl,kk)-c4les)
            zi = d_one/(t(jl,kk)-c4ies)
            za = xalpha(t(jl,kk))
            qs = c2es*(za*exp(c3les*(t(jl,kk)-tzero)*zl) + &
                  (d_one-za)*exp(c3ies*(t(jl,kk)-tzero)*zi))
            qs = qs*rp
            qs = min(qsmax,qs)
            cor = d_one - ep1*qs
            zf = za*c5alvcp*zl**2 + (d_one-za)*c5alscp*zi**2
            cond = (q(jl,kk)*cor**2-qs*cor)/(cor**2+qs*zf)
            if ( cond > d_zero ) then
              t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond
              q(jl,kk) = q(jl,kk) - cond
              zl = d_one/(t(jl,kk)-c4les)
              zi = d_one/(t(jl,kk)-c4ies)
              za = xalpha(t(jl,kk))
              qs = c2es*(za*exp(c3les*(t(jl,kk)-tzero)*zl) + &
                    (d_one-za)*exp(c3ies*(t(jl,kk)-tzero)*zi))
              qs = qs*rp
              qs = xmin(qsmax,qs)
              cor = d_one - ep1*qs
              zf = za*c5alvcp*zl**2 + (d_one-za)*c5alscp*zi**2
              cond1 = (q(jl,kk)*cor**2-qs*cor)/(cor**2+qs*zf)
              if ( abs(cond) < almostzero ) cond1 = d_zero
              t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
              q(jl,kk) = q(jl,kk) - cond1
            end if
          end if
        end do
      else if ( jcall == 2 ) then
        do concurrent ( jl = n1:n2 )
          if ( ldflag(jl) ) then
            rp = d_one/sp(jl)
            qs = fesat(t(jl,kk))*rp
            qs = min(qsmax,qs)
            cor = d_one/(d_one-ep1*qs)
            qs = qs*cor
            cond = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
            cond = min(cond,d_zero)
            t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond
            q(jl,kk) = q(jl,kk) - cond
            qs = fesat(t(jl,kk))*rp
            qs = min(qsmax,qs)
            cor = d_one/(d_one-ep1*qs)
            qs = qs*cor
            cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
            if ( abs(cond) < almostzero ) cond1 = min(cond1,d_zero)
            t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
            q(jl,kk) = q(jl,kk) - cond1
          end if
        end do
      else if ( jcall == 0 ) then
        do concurrent ( jl = n1:n2 )
          rp = d_one/sp(jl)
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
          q(jl,kk) = q(jl,kk) - cond1
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
          q(jl,kk) = q(jl,kk) - cond1
        end do
      else if ( jcall == 4 ) then
        do concurrent ( jl = n1:n2 )
          if ( ldflag(jl) ) then
            rp = d_one/sp(jl)
            qs = fesat(t(jl,kk))*rp
            qs = min(qsmax,qs)
            cor = d_one/(d_one-ep1*qs)
            qs = qs*cor
            cond = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
            t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond
            q(jl,kk) = q(jl,kk) - cond
            qs = fesat(t(jl,kk))*rp
            qs = min(qsmax,qs)
            cor = d_one/(d_one-ep1*qs)
            qs = qs*cor
            cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
            t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
            q(jl,kk) = q(jl,kk) - cond1
          end if
        end do
      else if ( jcall == 5 ) then ! Same as 4 but with LDFLAG all true
        do concurrent ( jl = n1:n2 )
          rp = d_one/sp(jl)
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond
          q(jl,kk) = q(jl,kk) - cond
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
          q(jl,kk) = q(jl,kk) - cond1
        end do
      else if ( jcall == 3 ) then
        do concurrent ( jl = n1:n2 )
          rp = d_one/sp(jl)
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
          q(jl,kk) = q(jl,kk) - cond1
          qs = fesat(t(jl,kk))*rp
          qs = min(qsmax,qs)
          cor = d_one/(d_one-ep1*qs)
          qs = qs*cor
          cond1 = (q(jl,kk)-qs)/(d_one+qs*cor*fdqsat(t(jl,kk)))
          t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
          q(jl,kk) = q(jl,kk) - cond1
        end do
      else
#ifndef TESTME
        call fatal(__FILE__,__LINE__, 'Unknown method jcall in moistadj')
#endif
      end if

    end subroutine moistadj
    !
    ! This routine calculates entrainment/detrainment rates for updrafts
    ! in cumulus parameterization. Input are environmental values t,q etc
    ! and updraft values t,q etc. It returns entrainment/detrainment rates
    ! Turbulent entrainment is simulated by a constant multiplied by a
    ! vertical scaling function
    !
    subroutine entrainm(kk,ldwork,dmfen,dmfde)
      implicit none
      integer(ik4), intent(in) :: kk ! input level
      logical, intent(in) :: ldwork                      ! enable/disable
      real(rkx), dimension(np), intent(out) :: dmfen    ! Entrainment
      real(rkx), dimension(np), intent(out) :: dmfde    ! Detrainment
      integer(ik4) :: n
      real(rkx) :: mf, entrpen, detrpen
      !-----------------------------------------------
      ! 1. Calculate entrainment and detrainment rates
      ! ----------------------------------------------
      if ( ldwork ) then
        do concurrent ( n = n1:n2 )
          dmfen(n) = d_zero
          dmfde(n) = d_zero
        end do
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) ) then
            if ( kk < kcbot(n) ) then
              mf = mfu(n,kk+1)*(geof(n,kk)-geof(n,kk+1))*regrav
              entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(n))
              detrpen = merge(detrpen_lnd,detrpen_ocn,ldland(n))
              dmfen(n) = entrpen*mf
              dmfde(n) = detrpen*mf
            end if
          end if
        end do
      end if
    end subroutine entrainm
    !
    ! This routine calculates cloud base values for midlevel convection
    ! Input are environmental values t,q etc. It returns cloudbase values
    ! for midlevel convection
    !
    subroutine mcbase(kk)
      implicit none
      integer(ik4), intent(in) :: kk                     ! current level
      integer(ik4) :: n
      !-----------------------------------------------
      ! 1. Calculate entrainment and detrainment rates
      ! ----------------------------------------------
      do concurrent ( n = n1:n2 )
        if ( .not.ldcum(n) .and. ilab(n,kk+1) == 0 ) then
          if ( lmfmid .and. geo(n,kk) >  5000.0_rkx .and. &
                            geo(n,kk) < 10000.0_rkx .and. &
                            q(n,kk) > 0.80_rkx*qs(n,kk) ) then
            tu(n,kk+1) = (cpd*t(n,kk)+geo(n,kk)-geof(n,kk+1))*rcpd
            qu(n,kk+1) = q(n,kk)
            lu(n,kk+1) = d_zero
            mfub(n) = min(max(cmfcmin,-omega(n,kk)*regrav),cmfcmax)
            mfu(n,kk+1) = mfub(n)
            mfus(n,kk+1) = mfub(n)*(cpd*tu(n,kk+1)+geof(n,kk+1))
            mfuq(n,kk+1) = mfub(n)*qu(n,kk+1)
            mful(n,kk+1) = d_zero
            dmfup(n,kk+1) = d_zero
            xrain(n,kk+1) = d_zero
            kcbot(n) = kk
            ilab(n,kk+1) = 1
            ktype(n) = 3
          end if
        end if
      end do
    end subroutine mcbase
    !
    ! Updates u and v tendencies, with a global diagnostic of dissipation
    ! Explicit upstream and implicit solution of vertical advection depending
    ! on value of rmfsoluv:
    !       0 = Explicit
    !     0-1 = Semi-Implicit
    !     >=1 = Implicit
    !
    ! For explicit solution: only one single iteration
    ! For implicit solution: first implicit solver, then explicit solver
    ! to correct tendencies below cloud base
    !
    subroutine dudvx(mfu,mfd)
      implicit none
      real(rkx), dimension(np,nk), intent(in) :: mfu   ! Mass Flux Updraft
      real(rkx), dimension(np,nk), intent(in) :: mfd   ! Mass Flux Downdraft
      real(rkx), dimension(np,nk) :: uen, ven, mfuu, mfdu, mfuv, mfdv
      integer(ik4) :: ik, ikb, k, n
      real(rkx) :: zp, ximp
      real(rkx), dimension(:,:), allocatable :: dudt, dvdt, dp
      real(rkx), dimension(:,:), allocatable :: bb, r1, r2
      logical, dimension(:,:), allocatable :: llcumbas
      ximp = d_one - rmfsoluv
      allocate (dudt(np,nk))
      allocate (dvdt(np,nk))
      allocate (dp(np,nk))
      do concurrent ( n = n1:n2, k = 1:nk )
        if ( ldcum(n) ) then
          uen(n,k) = u(n,k)
          ven(n,k) = v(n,k)
          dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
        end if
      end do
      !--------------------------------------------------
      ! 1. Calculate fluxes and update u and v tendencies
      ! -------------------------------------------------
      do concurrent ( n = n1:n2, k = itopm2:nk )
        ik = k - 1
        if ( ldcum(n) ) then
          mfuu(n,k) = mfu(n,k) * (uu(n,k) - ximp*uen(n,ik))
          mfuv(n,k) = mfu(n,k) * (vu(n,k) - ximp*ven(n,ik))
          mfdu(n,k) = mfd(n,k) * (ud(n,k) - ximp*uen(n,ik))
          mfdv(n,k) = mfd(n,k) * (vd(n,k) - ximp*ven(n,ik))
        end if
      end do
      ! linear fluxes below cloud
      if ( abs(rmfsoluv) < almostzero ) then
        do concurrent ( n = n1:n2, k = itopm2:nk )
          if ( ldcum(n) .and. k > kcbot(n) ) then
            ikb = kcbot(n)
            zp = ((pf(n,nk+1)-pf(n,k))/(pf(n,nk+1)-pf(n,ikb)))
            if ( ktype(n) == 3 ) zp = zp*zp
            mfuu(n,k) = mfuu(n,ikb)*zp
            mfuv(n,k) = mfuv(n,ikb)*zp
            mfdu(n,k) = mfdu(n,ikb)*zp
            mfdv(n,k) = mfdv(n,ikb)*zp
          end if
        end do
      end if
      ! ----------------------
      ! 1.2 Compute tendencies
      ! ----------------------
      do concurrent ( n = n1:n2, k = itopm2:nk )
        if ( k < nk ) then
          ik = k + 1
          if ( ldcum(n) ) then
            dudt(n,k) = dp(n,k)*(mfuu(n,ik)-mfuu(n,k)+mfdu(n,ik)-mfdu(n,k))
            dvdt(n,k) = dp(n,k)*(mfuv(n,ik)-mfuv(n,k)+mfdv(n,ik)-mfdv(n,k))
          end if
        else
          if ( ldcum(n) ) then
            dudt(n,k) = -dp(n,k)*(mfuu(n,k)+mfdu(n,k))
            dvdt(n,k) = -dp(n,k)*(mfuv(n,k)+mfdv(n,k))
          end if
        end if
      end do
      if ( abs(rmfsoluv) < almostzero ) then
        !----------------------
        ! 1.3 Update tendencies
        ! ---------------------
        do concurrent ( n = n1:n2, k = itopm2:nk )
          if ( ldcum(n) ) then
            tenu(n,k) = tenu(n,k) + dudt(n,k)
            tenv(n,k) = tenv(n,k) + dvdt(n,k)
          end if
        end do
      else
        !------------------------------------------------------------------
        ! 1.4 Implicit solution
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k; reuse A and B
        ! DUDT and DVDT correspond to the RHS ("constants") of the equation
        ! The solution is in R1 and R2
        ! -----------------------------------------------------------------
        allocate (bb(np,nk), source=1.0_rkx)
        allocate (r1(np,nk), source=0.0_rkx)
        allocate (r2(np,nk), source=0.0_rkx)
        allocate (llcumbas(np,nk), source=.false.)
        do concurrent ( n = n1:n2, k = 1:nk )
          mfuu(n,k) = d_zero
        end do
        ! Fill vectors A, B and RHS
        do concurrent ( n = n1:n2, k = itopm2:nk )
          ik = k + 1
          llcumbas(n,k) = ldcum(n) .and. k >= kctop(n) - 1
          if ( llcumbas(n,k) ) then
            zp = rmfsoluv*dp(n,k)*dtc
            mfuu(n,k) = -zp*(mfu(n,k)+mfd(n,k))
            dudt(n,k) = dudt(n,k)*dtc + uen(n,k)
            dvdt(n,k) = dvdt(n,k)*dtc + ven(n,k)
            if ( k < nk ) then
              bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
            else
              bb(n,k) = d_one
            end if
          end if
        end do
        call solver(n1,n2,nk,kctop,llcumbas,mfuu,bb,dudt,r1)
        call solver(n1,n2,nk,kctop,llcumbas,mfuu,bb,dvdt,r2)
        do concurrent ( n = n1:n2, k = itopm2:nk )
          if ( llcumbas(n,k) ) then
            tenu(n,k) = tenu(n,k) + (r1(n,k)-uen(n,k))/dtc
            tenv(n,k) = tenv(n,k) + (r2(n,k)-ven(n,k))/dtc
          end if
        end do
        deallocate (llcumbas)
        deallocate (r2)
        deallocate (r1)
        deallocate (bb)
      end if
      deallocate (dp)
      deallocate (dvdt)
      deallocate (dudt)
    end subroutine dudvx
    !
    ! Compute convective transport of chemical tracers
    ! Explicit upstream and implicit solution of vertical advection
    ! depending on value of rmfsolct:
    !         0 = explicit
    !       0-1 = semi-implicit
    !       >=1 = implicit
    !
    ! For explicit solution: only one single iteration
    ! For implicit solution: first implicit solver, then explicit solver
    !                        To correct tendencies below cloud base
    !
    subroutine ctracer(ldcum,lddraf,mfu,mfd,pudrate,pddrate)
      implicit none
      logical, dimension(np), intent(in) :: ldcum
      logical, dimension(np), intent(in) :: lddraf
      real(rkx), dimension(np,nk), intent(in) :: mfu
      real(rkx), dimension(np,nk), intent(in) :: mfd
      real(rkx), dimension(np,nk), intent(in) :: pudrate
      real(rkx), dimension(np,nk), intent(in) :: pddrate
      integer(ik4) :: ik, k, n, nt
      real(rkx) :: zp, mfa, ximp, erate, posi
      real(rkx), dimension(:,:,:), allocatable :: cen, cu, cd, xtenc, mfc
      real(rkx), dimension(:,:), allocatable :: dp, bb, r1
      logical, dimension(:,:), allocatable :: llcumask, llcumbas

      ximp = d_one - rmfsolct
      allocate (cen(np,nk,ntrac))   ! Half-level environmental values
      allocate (cu(np,nk,ntrac))    ! Updraft values
      allocate (cd(np,nk,ntrac))    ! Downdraft values
      allocate (xtenc(np,nk,ntrac)) ! Tendency
      allocate (mfc(np,nk,ntrac))   ! Fluxes
      allocate (dp(np,nk))          ! Pressure difference
      allocate (llcumask(np,nk))    ! Mask for convection
      !-----------------------------------------
      ! 0. Initialize Cumulus mask + some setups
      ! ----------------------------------------
      do concurrent ( n = n1:n2, k = 2:nk )
        llcumask(n,k) = .false.
        if ( ldcum(n) ) then
          dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
          if ( k >= kctop(n)-1 ) llcumask(n,k) = .true.
        end if
      end do
      do nt = 1, ntrac
        !---------------------------------------
        ! 1. Define tracers at full sigma levels
        ! --------------------------------------
        do concurrent ( n = n1:n2, k = 2:nk )
          ik = k - 1
          cen(n,k,nt) = qtrac(n,k,nt)
          cd(n,k,nt) = qtrac(n,ik,nt)
          cu(n,k,nt) = qtrac(n,ik,nt)
          mfc(n,k,nt) = d_zero
          xtenc(n,k,nt) = d_zero
        end do
        do concurrent ( n = n1:n2 )
          cu(n,nk,nt) = qtrac(n,nk,nt)
        end do
        !--------------------------
        ! 2. Compute updraft values
        ! -------------------------
        do concurrent ( n = n1:n2, k = nk-1:3:-1 )
          ik = k + 1
          if ( llcumask(n,k) ) then
            erate = mfu(n,k) - mfu(n,ik) + pudrate(n,k)
            mfa = d_one/max(cmfcmin,mfu(n,k))
            if ( k >= kctop(n) ) then
              cu(n,k,nt) = (mfu(n,ik)*cu(n,ik,nt) + &
                           erate*qtrac(n,k,nt) - &
                           pudrate(n,k)*cu(n,ik,nt))*mfa
            end if
          end if
        end do
        !----------------------------
        ! 3. Compute downdraft values
        ! ---------------------------
        do concurrent ( n = n1:n2, k = 3:nk )
          ik = k - 1
          if ( lddraf(n) .and. k == idtop(n) ) then
            ! Note: in order to avoid final negative tracer values
            ! the allowed value of cd depends on the jump in mass flux
            ! at the LFS
            cd(n,k,nt) = 0.1_rkx*cu(n,k,nt) + 0.9_rkx*qtrac(n,ik,nt)
          else if ( lddraf(n) .and. k > idtop(n) ) then
            erate = -mfd(n,k) + mfd(n,ik) + pddrate(n,k)
            mfa = d_one/min(-cmfcmin,mfd(n,k))
            cd(n,k,nt) = (mfd(n,ik)*cd(n,ik,nt) - &
                         erate*qtrac(n,ik,nt) + &
                         pddrate(n,k)*cd(n,ik,nt))*mfa
          end if
        end do
        ! In order to avoid negative Tracer at nk adjust cd
        k = nk
        ik = k - 1
        do concurrent ( n = n1:n2 )
          if ( lddraf(n) ) then
            posi = -dp(n,k) *(mfu(n,k)*cu(n,k,nt) + &
              mfd(n,k)*cd(n,k,nt)-(mfu(n,k)+mfd(n,k))*qtrac(n,ik,nt))
            if ( qtrac(n,k,nt)+posi*dtc < d_zero ) then
              mfa = d_one/min(-cmfcmin,mfd(n,k))
              cd(n,k,nt) = ((mfu(n,k)+mfd(n,k))*qtrac(n,ik,nt) - &
                mfu(n,k)*cu(n,k,nt)+qtrac(n,k,nt) / &
                (dtc*dp(n,k)))*mfa
            end if
          end if
        end do
      end do

      do nt = 1, ntrac
        !------------------
        ! 4. Compute fluxes
        ! -----------------
        do concurrent ( n = n1:n2, k = 2:nk )
          ik = k - 1
          if ( llcumask(n,k) ) then
            mfa = mfu(n,k) + mfd(n,k)
            mfc(n,k,nt) = mfu(n,k)*cu(n,k,nt) + &
              mfd(n,k)*cd(n,k,nt) - ximp*mfa*cen(n,ik,nt)
          end if
        end do
        !----------------------------
        ! 5. Compute tendencies = rhs
        ! ---------------------------
        do concurrent ( n = n1:n2, k = 2:nk-1 )
          ik = k + 1
          if ( llcumask(n,k) ) then
            xtenc(n,k,nt) = dp(n,k)*(mfc(n,ik,nt)-mfc(n,k,nt))
          end if
        end do
        k = nk
        do concurrent ( n = n1:n2 )
          if ( ldcum(n) ) xtenc(n,k,nt) = -dp(n,k)*mfc(n,k,nt)
        end do
      end do
      if ( abs(rmfsolct) < almostzero ) then
        !---------------------
        ! 6. Update tendencies
        ! --------------------
        do concurrent ( n = n1:n2, k = 2:nk, nt = 1:ntrac )
          if ( llcumask(n,k) ) then
            tenc(n,k,nt) = tenc(n,k,nt)+xtenc(n,k,nt)
          end if
        end do
      else
        !---------------------
        ! 7. Implicit solution
        ! ---------------------
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
        ! tenc corresponds to the RHS ("constants") of the equation
        ! The solution is in R1
        allocate (bb(np,nk), source=1.0_rkx)
        allocate (r1(np,nk), source=0.0_rkx)
        allocate (llcumbas(np,nk), source=.false.)
        do nt = 1, ntrac
          ! Fill vectors A, B and RHS
          do concurrent ( n = n1:n2, k = 2:nk )
            ik = k + 1
            llcumbas(n,k) = llcumask(n,k)
            if ( llcumbas(n,k) ) then
              zp = rmfsolct*dp(n,k)*dtc
              mfc(n,k,nt) = -zp*(mfu(n,k)+mfd(n,k))
              xtenc(n,k,nt) = xtenc(n,k,nt)*dtc + qtrac(n,k,nt)
              ! for implicit solution including tendency source term
              if ( k < nk ) then
                bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
              else
                bb(n,k) = d_one
              end if
            end if
          end do
          call solver(n1,n2,nk,kctop,llcumbas,mfc(:,:,nt),bb,xtenc(:,:,nt),r1)
          ! Compute tendencies
          do concurrent ( n = n1:n2, k = 2:nk )
            !  for implicit solution including tendency source term
            !  tenc(n,k,nt) = (r1(n,k)-qtrac(n,k,nt))/dtc
            if ( llcumbas(n,k) ) then
              tenc(n,k,nt) = tenc(n,k,nt) + (r1(n,k)-qtrac(n,k,nt))/dtc
            end if
          end do
        end do
        deallocate (llcumbas)
        deallocate (bb)
        deallocate (r1)
      end if

      deallocate (llcumask)
      deallocate (dp)
      deallocate (mfc)
      deallocate (xtenc)
      deallocate (cd)
      deallocate (cu)
      deallocate (cen)
    end subroutine ctracer

  end subroutine ntiedtke
  !
  ! Computes specific humidity at saturation used by the diagnostic
  ! cloud scheme to compute relative humidity and liquid water content
  !
  subroutine satur(n1,n2,nk,kt,pr,t,qsat)
    implicit none
    integer(ik4), intent(in) :: n1, n2, nk              ! dimensions
    integer(ik4), intent(in) :: kt                      ! top level
    real(rkx), dimension(n1:n2,nk), intent(in) :: pr    ! Pressure Pa
    real(rkx), dimension(n1:n2,nk), intent(in) :: t     ! Temperature K
    real(rkx), dimension(n1:n2,nk), intent(out) :: qsat ! Satur. MXR kg/kg
    integer(ik4) :: k, n
    real(rkx) :: qs
    do concurrent ( n = n1:n2, k = kt:nk )
      qs = min(qsmax,fesat(t(n,k))/pr(n,k))
      qsat(n,k) = qs/(d_one-ep1*qs)
    end do
  end subroutine satur
  !
  ! Solves bidiagonal system for implicit solution of advection equation
  ! interface. It returns updated value of quantity.
  !
  ! Numerical Recipes (Cambridge Press)
  ! Derived from tridiagonal algorihm with c=0, only one forward
  ! substitution necessary.
  !
  !          M  x  U  = R
  !          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C
  !          (  B(ktop-1)     C(ktop-1)     0          0        )
  !          (  A(ktop)       B(ktop)       C(ktop)    0        )
  !          (  0             A(k)          B(k)       C(k)     )
  !          (  0             0             A(nk)      B(nk)    )
  !
  subroutine solver(n1,n2,nk,ktop,wmask,a,b,r,u)
    implicit none
    integer(ik4), intent(in) :: n1, n2, nk          ! dimensions
    integer(ik4), dimension(n1:n2), intent(in) :: ktop ! cloud top
    logical, dimension(n1:n2,nk), intent(in) :: wmask  ! mask
    real(rkx), dimension(n1:n2,nk), intent(in) :: a    ! lower diagonal
    real(rkx), dimension(n1:n2,nk), intent(in) :: b    ! diagonal
    real(rkx), dimension(n1:n2,nk), intent(in) :: r    ! right hand side
    real(rkx), dimension(n1:n2,nk), intent(out) :: u   ! solution
    integer(ik4) :: k, n
    !$acc kernels
    u(:,:) = d_zero
    !$acc end kernels
    do concurrent ( n = n1:n2 )
      !$acc loop seq
      do k = 2, nk
        if ( wmask(n,k) ) then
          if ( abs(b(n,k)) > almostzero ) then
            if ( k == ktop(n)-1 ) then
              u(n,k) = r(n,k)/b(n,k)
            else if ( k > ktop(n)-1 ) then
              u(n,k) = (r(n,k)-a(n,k)*u(n,k-1))/b(n,k)
            end if
          end if
        end if
      end do
    end do
  end subroutine solver

  pure real(rkx) function xmin(x,y)
    implicit none
    real(rkx), intent(in) :: x, y
    xmin = y - d_half*(abs(x-y)-(x-y))
  end function xmin

  pure real(rkx) function lwocp(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rkx) :: gtzero
    gtzero = max(d_zero,sign(d_one,t-tzero))
    lwocp = gtzero*wlhvocp + (d_one-gtzero)*wlhsocp
  end function lwocp

  pure real(rkx) function xalpha(t)
    implicit none
    real(rkx), intent(in) :: t
    xalpha = min(d_one,((max(rtice,min(rtwat,t))-rtice)*rtwat_rtice_r)**2)
  end function xalpha
  !
  ! Magnus Tetens formula
  !
  pure real(rkx) function fesat(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rkx) :: xa
    xa = xalpha(t)
    fesat = c2es*(xa*exp((c3les*((t-tzero)/(t-c4les)))) + &
          (d_one-xa)*exp((c3ies*((t-tzero)/(t-c4ies)))))
  end function fesat

  pure real(rkx) function fdqsat(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rkx) :: xa
    xa = xalpha(t)
    fdqsat = xa*c5alvcp*(d_one/(t-c4les)**2) + &
          (d_one-xa)*c5alscp*(d_one/(t-c4ies)**2)
  end function fdqsat

  pure real(rkx) function mlwocp(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rkx) :: xa
    xa = xalpha(t)
    mlwocp = xa*wlhvocp+(d_one-xa)*wlhsocp
  end function mlwocp

  pure real(rkx) function esw(t)
    implicit none
    real(rkx), intent(in) :: t
    esw = c3les*(t-tzero)/(t-c4les)
  end function esw

  pure real(rkx) function esi(t)
    implicit none
    real(rkx), intent(in) :: t
    esi = c3ies*(t-tzero)/(t-c4ies)
  end function esi

  pure real(rkx) function mlwt(t)
    implicit none
    real(rkx), intent(in) :: t
    real(rkx) :: xa
    xa = xalpha(t)
    mlwt = xa*wlhv+(d_one-xa)*wlhs
  end function mlwt

end module mod_cu_tiedtke
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
