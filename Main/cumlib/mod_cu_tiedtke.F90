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
  use mod_runparams , only : iqc , dt , iqv , iqi , entrmax , dx , &
         entrdd , entrmid , cprcon , entrpen_lnd , entrpen_ocn ,   &
         entrscv , iconv , ichem , iaerosol , iindirect , ipptls , &
         hsigma , sigma , ichcumtra , rcmtimer , icup , dtcum
  use mod_runparams , only : k2_const , kfac_deep , kfac_shal
  use mod_mpmessage
  use mod_runparams , only : rcrit
  use mod_runparams , only : rprc_ocn , rprc_lnd , revap_lnd , revap_ocn
  use mod_runparams , only : detrpen_lnd , detrpen_ocn , entshalp , entrdd
  use mod_runparams , only : rhebc_ocn , rhebc_lnd , rcuc_lnd , rcuc_ocn
  use mod_runparams , only : rcpec_ocn , rcpec_lnd , cmtcape
  use mod_runparams , only : lmfpen , lmfmid , lmfdd , lepcld , lmfdudv , &
          lmfscv , lmfuvdis , lmftrac , lmfsmooth , lmfwstar
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_cu_tiedtke , tiedtkedrv

  real(rkx) , parameter :: rhow = 1000.0_rkx
  real(rkx) , parameter :: rkap = 0.4_rkx
  real(rkx) , parameter :: qsmax = 0.5_rkx
  real(rkx) , parameter :: cwdrag = (3.0_rkx/8.0_rkx)*0.506_rkx/0.200_rkx

  real(rkx) :: rtau
  real(rkx) :: rmfcfl ! Massflux multiple of cfl stability criterium
  integer(ik4) :: nk350 , nk060 , nk950

  ! Mass flux solver for momentum (1) or no (0)
  real(rkx) , parameter :: rmfsoluv = 1.0_rkx
  ! Mass flux solver for T and q (1) or no (0)
  real(rkx) , parameter :: rmfsoltq = 1.0_rkx
  ! Mass flux solver for tracer (1) or no (0)
  real(rkx) , parameter :: rmfsolct = 1.0_rkx

  ! Use CFL mass flux limit (1) or absolut limit (0)
  real(rkx) , parameter :: rmflic = 1.0_rkx
  ! Value of absolut mass flux limit
  real(rkx) , parameter :: rmflia = 0.0_rkx

  ! Relaxation time for melting of snow
  real(rkx) , parameter :: rtaumel = 5.0_rkx*3600.0_rkx*1.5_rkx

  ! Updraught velocity perturbation for implicit (m/s)
  real(rkx) , parameter :: ruvper = 0.3_rkx

  ! Maximum allowed cloud thickness for shallow cloud (Pa)
  real(rkx) , parameter :: rdepths = 4.0e4_rkx

  ! Fractional massflux for downdrafts at lfs
  real(rkx) , parameter :: rmfdeps = 0.3_rkx

  ! evaporation coefficient for kuo0
  integer(ik4) , pointer , dimension(:,:) :: ilab
  integer(ik4) , pointer , dimension(:) :: ktype
  logical , pointer , dimension(:) :: ldland
  real(rkx) , pointer , dimension(:,:,:) :: pxtm1 , pxtte
  real(rkx) , pointer , dimension(:,:) :: ptm1 , pqm1 , pum1 , pvm1 ,    &
        pxlm1 , pxim1 , pxite , papp1 , paphp1 , pxtec , pqtec , zlude , &
        pmflxr , pccn
  real(rkx) , pointer , dimension(:) :: prsfc , pssfc , &
        ptopmax , xphfx , xpqfx
  real(rkx) , pointer , dimension(:,:) :: ptte , pvom , pvol , pqte , &
        pxlte , pverv , pmfu , xpg , xpgh
  integer(ik4) , pointer , dimension(:) :: kctop , kcbot
  integer(ik4) :: nipoi
  integer(ik4) , pointer , dimension(:) :: imap , jmap
  integer(ik4) :: nskmax

  integer(ik4) , public :: nmctop
  real(rkx) , pointer , public , dimension(:) :: pmean

  ! Max massflux value
  real(rkx) , parameter :: cmfcmax = 1.0_rkx
  ! Min massflux value (for safety)
  real(rkx) , parameter :: cmfcmin = 1.0e-10_rkx
  ! Rainfall max elevation
  real(rkx) , parameter :: zdlev   = 1.5e4_rkx
  ! Relat. cloud massflux at level above nonbuoyancy
  real(rkx) , parameter :: cmfctop = 0.35_rkx

  ! Fractional massflux for downdrafts at lfs
  real(rkx) , parameter :: cmfdeps = 0.3_rkx

  ! Trigger parameter for the convection 0.0 <-> 1.0
  real(rkx) , parameter :: ctrigger = -1.1_rkx

  real(rkx) , parameter :: almostzero = dlowval

  contains
  !
  ! This subroutines allocates space
  !
  subroutine allocate_mod_cu_tiedtke
    implicit none
    integer(ik4) :: i , j , ii
    call getmem1d(pmean,1,kz,'mod_cu_tiedtke:pmean')
    nipoi = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( cuscheme(j,i) == 5 ) then
          nipoi = nipoi + 1
        end if
      end do
    end do
    if ( nipoi == 0 ) then
      return
    end if
    call getmem1d(imap,1,nipoi,'mod_cu_tiedtke:imap')
    call getmem1d(jmap,1,nipoi,'mod_cu_tiedtke:jmap')
    ii = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( cuscheme(j,i) == 5 ) then
          imap(ii) = i
          jmap(ii) = j
          ii = ii + 1
        end if
      end do
    end do
    call getmem2d(ptte,1,nipoi,1,kz,'mod_cu_tiedtke:ptte')
    call getmem2d(pvom,1,nipoi,1,kz,'mod_cu_tiedtke:pvom')
    call getmem2d(pvol,1,nipoi,1,kz,'mod_cu_tiedtke:pvol')
    call getmem2d(pqte,1,nipoi,1,kz,'mod_cu_tiedtke:pqte')
    call getmem2d(pxlte,1,nipoi,1,kz,'mod_cu_tiedtke:pxlte')
    call getmem2d(pverv,1,nipoi,1,kz,'mod_cu_tiedtke:pverv')
    call getmem2d(pmfu,1,nipoi,1,kz,'mod_cu_tiedtke:pmfu')
    call getmem2d(xpg,1,nipoi,1,kz,'mod_cu_tiedtke:xpg')
    call getmem2d(xpgh,1,nipoi,1,kz+1,'mod_cu_tiedtke:xpgh')
    call getmem3d(pxtm1,1,nipoi,1,kz,1,ntr,'mod_cu_tiedtke:pxtm1')
    call getmem3d(pxtte,1,nipoi,1,kz,1,ntr,'mod_cu_tiedtke:pxtte')
    call getmem2d(ilab,1,nipoi,1,kz,'mod_cu_tiedtke:ilab')
    call getmem1d(ktype,1,nipoi,'mod_cu_tiedtke:ktype')
    call getmem2d(ptm1,1,nipoi,1,kz,'mod_cu_tiedtke:ptm1')
    call getmem2d(pqm1,1,nipoi,1,kz,'mod_cu_tiedtke:pqm1')
    call getmem2d(pum1,1,nipoi,1,kz,'mod_cu_tiedtke:pum1')
    call getmem2d(pvm1,1,nipoi,1,kz,'mod_cu_tiedtke:pvm1')
    call getmem2d(pxlm1,1,nipoi,1,kz,'mod_cu_tiedtke:pxlm1')
    call getmem2d(pxim1,1,nipoi,1,kz,'mod_cu_tiedtke:pxim1')
    call getmem2d(pxite,1,nipoi,1,kz,'mod_cu_tiedtke:pxite')
    call getmem2d(papp1,1,nipoi,1,kz,'mod_cu_tiedtke:papp1')
    call getmem2d(zlude,1,nipoi,1,kz,'mod_cu_tiedtke:zlude')
    call getmem2d(pxtec,1,nipoi,1,kz,'mod_cu_tiedtke:pxtec')
    call getmem2d(pqtec,1,nipoi,1,kz,'mod_cu_tiedtke:pqtec')
    if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
      call getmem2d(pccn,1,nipoi,1,kz,'mod_cu_tiedtke:pccn')
    end if
    call getmem2d(pmflxr,1,nipoi,1,kz+1,'mod_cu_tiedtke:pmflxr')
    call getmem1d(kctop,1,nipoi,'mod_cu_tiedtke:kctop')
    call getmem1d(kcbot,1,nipoi,'mod_cu_tiedtke:kcbot')
    call getmem2d(paphp1,1,nipoi,1,kzp1,'mod_cu_tiedtke:paphp1')
    call getmem1d(prsfc,1,nipoi,'mod_cu_tiedtke:prsfc')
    call getmem1d(pssfc,1,nipoi,'mod_cu_tiedtke:pssfc')
    call getmem1d(ptopmax,1,nipoi,'mod_cu_tiedtke:ptopmax')
    call getmem1d(xphfx,1,nipoi,'mod_cu_tiedtke:xphfx')
    call getmem1d(xpqfx,1,nipoi,'mod_cu_tiedtke:xpqfx')
    call getmem1d(ldland,1,nipoi,'mod_cu_tiedtke:ldland')

  end subroutine allocate_mod_cu_tiedtke
  !
  ! This subroutines calls cucall
  !
  subroutine tiedtkedrv(m2c,uxten,vxten)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: uxten
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: vxten
    integer(ik4) :: i , j , k , n , ii
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tiedtkedrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nipoi == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    ilab(:,:) = 2

    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do ii = 1 , nipoi
            i = imap(ii)
            j = jmap(ii)
            ! tracers input profile : implicit loop on tracer
            pxtm1(ii,k,n) = m2c%chias(j,i,k,n)
          end do
        end do
      end do
      if ( idynamic == 3 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do ii = 1 , nipoi
              i = imap(ii)
              j = jmap(ii)
              pxtte(ii,k,n) = m2c%chiten(j,i,k,n)
            end do
          end do
        end do
      else
        do n = 1 , ntr
          do k = 1 , kz
            do ii = 1 , nipoi
              i = imap(ii)
              j = jmap(ii)
              pxtte(ii,k,n) = m2c%chiten(j,i,k,n)/m2c%psb(j,i)
            end do
          end do
        end do
      end if
      if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
        do k = 1 , kz
          do ii = 1 , nipoi
            i = imap(ii)
            j = jmap(ii)
            pccn(ii,k) = m2c%ccn(j,i,k)
          end do
        end do
      end if
    else
      pxtm1(:,:,:) = d_zero ! tracers input profiles
      pxtte(:,:,:) = d_zero ! tracer tendencies
    end if

    do ii = 1 , nipoi
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

    do k = 1 , kz
      do ii = 1 , nipoi
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
        pvom(ii,k)  = uxten(j,i,k)
        pvol(ii,k)  = vxten(j,i,k)
        ! scheme diagnostic output - tendencies due to convection
        pxtec(ii,k) = d_zero ! detrained cloud water tendency
        pqtec(ii,k) = d_zero ! detrained humidity tendency
        ! pressure top limit for convection: the level above tropopause
      end do
    end do
    if ( idynamic == 3 ) then
      do k = 1 , kz
        do ii = 1 , nipoi
          i = imap(ii)
          j = jmap(ii)
          ptte(ii,k)  = m2c%tten(j,i,k)
          pqte(ii,k)  = m2c%qxten(j,i,k,iqv)
          pxlte(ii,k) = m2c%qxten(j,i,k,iqc)
        end do
      end do
    else
      do k = 1 , kz
        do ii = 1 , nipoi
          i = imap(ii)
          j = jmap(ii)
          ptte(ii,k)  = m2c%tten(j,i,k)/m2c%psb(j,i)
          pqte(ii,k)  = m2c%qxten(j,i,k,iqv)/m2c%psb(j,i)
          pxlte(ii,k) = m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
        end do
      end do
    end if
    if ( ipptls > 1 ) then
      if ( idynamic == 3 ) then
        do k = 1 , kz
          do ii = 1 , nipoi
            i = imap(ii)
            j = jmap(ii)
            pxim1(ii,k) = m2c%qxas(j,i,k,iqi)      ! cloud ice water
            pxite(ii,k) = m2c%qxten(j,i,k,iqi)
          end do
        end do
      else
        do k = 1 , kz
          do ii = 1 , nipoi
            i = imap(ii)
            j = jmap(ii)
            pxim1(ii,k) = m2c%qxas(j,i,k,iqi)      ! cloud ice water
            pxite(ii,k) = m2c%qxten(j,i,k,iqi)/m2c%psb(j,i)
          end do
        end do
      end if
    else
      pxim1(:,:) = d_zero
      pxite(:,:) = d_zero
    end if

    ! Transform to specific humidity
    pqm1(:,:) = pqm1(:,:)/(d_one + pqm1(:,:))
    pqte(:,:) = pqte(:,:)/(d_one + pqm1(:,:))**2

    nskmax = nint(17747.5/(dx*d_r1000))
    if ( iconv == 4 ) then
      call setup(nskmax,pmean)
    else
      do ii = 1 , nipoi
        i = imap(ii)
        j = jmap(ii)
        ptopmax(ii) = papp1(ii,max(1,m2c%ktrop(j,i)-1))
      end do
    end if

    do k = 1 , kzp1
      do ii = 1 , nipoi
        i = imap(ii)
        j = jmap(ii)
        ! 1st guess pressure at full levels
        paphp1(ii,k) = m2c%pasf(j,i,k)
        xpgh(ii,k) = m2c%zfs(j,i,k)*egrav ! geopotential
      end do
    end do

    ! Output variables (1d)
    prsfc(:) = d_zero ! CHECK - surface rain flux
    pssfc(:) = d_zero ! CHECK - surface snow flux

    ktype(:) = 0
    kctop(:) = 0
    kcbot(:) = 0
    pmfu(:,:) = d_zero

    call cucall(nipoi,nipoi,kz,kzp1,kzm1,ilab,ntr,pxtm1,pxtte,ptm1,   &
                pqm1,pum1,pvm1,pxlm1,pxim1,ptte,pqte,pvom,pvol,pxlte, &
                pxite,pverv,pxtec,pqtec,xphfx,xpqfx,papp1,paphp1,xpg, &
                xpgh,prsfc,pssfc,pmflxr,pmfu,zlude,ktype,ldland,kctop,&
                kcbot,ptopmax,pccn)
    !
    ! postprocess some fields including precipitation fluxes
    !
    !   KTYPE = 1 => DEEP CONVECTION
    !   KTYPE = 2 => SHALLOW CONVECTION
    !   KTYPE = 3 => MIDLEVEL CONVECTION
    !
    do ii = 1 , nipoi
      if (ktype(ii) > 0) then
        i = imap(ii)
        j = jmap(ii)
        total_precip_points = total_precip_points + 1
        ! rainfall for surface
        cu_prate(j,i)= cu_prate(j,i) + prsfc(ii) + pssfc(ii)
      end if
    end do
    !
    ! update tendencies - note that rate were ADDED in cudtdq
    !                     thus here we must reset the rates.
    if ( idynamic == 3 ) then
      do k = 1 , kz
        do ii = 1 , nipoi
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_tten(j,i,k) = ptte(ii,k) - m2c%tten(j,i,k)
            ! Tendency in specific humidity to mixing ratio tendency.
            cu_qten(j,i,k,iqv) = pqte(ii,k)/(d_one-pqm1(ii,k))**2 - &
                                 m2c%qxten(j,i,k,iqv)
          end if
        end do
      end do
    else
      do k = 1 , kz
        do ii = 1 , nipoi
          if ( ktype(ii) > 0 ) then
            i = imap(ii)
            j = jmap(ii)
            cu_tten(j,i,k) = ptte(ii,k) - m2c%tten(j,i,k)/m2c%psb(j,i)
            ! Tendency in specific humidity to mixing ratio tendency.
            cu_qten(j,i,k,iqv) = pqte(ii,k)/(d_one-pqm1(ii,k))**2 - &
                                 m2c%qxten(j,i,k,iqv)/m2c%psb(j,i)
          end if
        end do
      end do
    end if
    do k = 1 , kz
      do ii = 1 , nipoi
        if ( ktype(ii) > 0 ) then
          i = imap(ii)
          j = jmap(ii)
          cu_uten(j,i,k) = pvom(ii,k) - uxten(j,i,k)
          cu_vten(j,i,k) = pvol(ii,k) - vxten(j,i,k)
          cu_qdetr(j,i,k) = zlude(ii,k)*egrav/(paphp1(ii,k+1)-paphp1(ii,k))
          cu_raincc(j,i,k) = pmflxr(ii,k)
        end if
      end do
    end do

    if ( ipptls > 1 ) then
      if ( idynamic == 3 ) then
        do k = 1 , kz
          do ii = 1 , nipoi
            if ( ktype(ii) > 0 ) then
              i = imap(ii)
              j = jmap(ii)
              cu_qten(j,i,k,iqc) = pxlte(ii,k) - m2c%qxten(j,i,k,iqc)
              cu_qten(j,i,k,iqi) = pxite(ii,k) - m2c%qxten(j,i,k,iqi)
            end if
          end do
        end do
      else
        do k = 1 , kz
          do ii = 1 , nipoi
            if ( ktype(ii) > 0 ) then
              i = imap(ii)
              j = jmap(ii)
              cu_qten(j,i,k,iqc) = pxlte(ii,k) - &
                  m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
              cu_qten(j,i,k,iqi) = pxite(ii,k) - &
                  m2c%qxten(j,i,k,iqi)/m2c%psb(j,i)
            end if
          end do
        end do
      end if
    else
      if ( idynamic == 3 ) then
        do k = 1 , kz
          do ii = 1 , nipoi
            if ( ktype(ii) > 0 ) then
              i = imap(ii)
              j = jmap(ii)
              cu_qten(j,i,k,iqc) = (pxlte(ii,k)+pxite(ii,k)) - &
                  m2c%qxten(j,i,k,iqc)
              cu_tten(j,i,k) = cu_tten(j,i,k) + pxite(ii,k)*wlhf/cpd
            end if
          end do
        end do
      else
        do k = 1 , kz
          do ii = 1 , nipoi
            if ( ktype(ii) > 0 ) then
              i = imap(ii)
              j = jmap(ii)
              cu_qten(j,i,k,iqc) = (pxlte(ii,k) + pxite(ii,k)) - &
                  m2c%qxten(j,i,k,iqc)/m2c%psb(j,i)
              cu_tten(j,i,k) = cu_tten(j,i,k) + pxite(ii,k)*wlhf/cpd
            end if
          end do
        end do
      end if
    end if

    if ( ichem == 1 .and. ichcumtra == 1 .and. &
         .not. any(icup == 2 .or. icup == 6) ) then
      if ( idynamic == 3 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do ii = 1 , nipoi
              if ( ktype(ii) > 0 ) then
                i = imap(ii)
                j = jmap(ii)
                cu_chiten(j,i,k,n) = pxtte(ii,k,n) - m2c%chiten(j,i,k,n)
              end if
            end do
          end do
        end do
      else
        do n = 1 , ntr
          do k = 1 , kz
            do ii = 1 , nipoi
              if ( ktype(ii) > 0 ) then
                i = imap(ii)
                j = jmap(ii)
                cu_chiten(j,i,k,n) = pxtte(ii,k,n) - &
                            m2c%chiten(j,i,k,n)/m2c%psb(j,i)
              end if
            end do
          end do
        end do
      end if
      ! build for chemistry 3d table of precipitation rate
      ! from the surface to the top of the convection
      do k = 1 , kz
        do ii = 1 , nipoi
          if (ktype(ii) > 0) then
            i = imap(ii)
            j = jmap(ii)
            if ( k > kctop(ii) ) then
              cu_convpr(j,i,k) = (prsfc(ii)+pssfc(ii))
            end if
          end if
        end do
      end do
    end if

    !
    ! Xu, K.-M., and S. K. Krueger:
    !   Evaluation of cloudiness parameterizations using a cumulus
    !   ensemble model, Mon. Wea. Rev., 119, 342-367, 1991.
    !
    do ii = 1 , nipoi
      if (ktype(ii) > 0) then
        i = imap(ii)
        j = jmap(ii)
        cu_ktop(j,i) = kctop(ii)
        cu_kbot(j,i) = kcbot(ii)
        do k = kctop(ii) , kcbot(ii)
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
      integer(ik4) , intent(in) :: smax
      real(rkx) , dimension(kz) , intent(in) :: pmean
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
      do k = kz , 2 , -1
        if ( pmean(k)/pmean(kz)*stdp >= 350.e2_rkx ) nk350 = k
        if ( pmean(k)/pmean(kz)*stdp >=  60.e2_rkx ) nk060 = k
        if ( pmean(k)/pmean(kz)*stdp >= 950.e2_rkx ) nk950 = k
      end do
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
                    kctop,kcbot,ptopmax,pccn)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim,klev) :: ilab
    integer(ik4) , dimension(kbdim) :: ktype
    logical , dimension(kbdim) :: ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1 , pgeoh
    real(rkx) , dimension(kbdim,klevp1) , intent(out) :: pmflxr
    real(rkx) , dimension(kbdim,klev) :: papp1 , pgeo , pqm1 , pqte ,   &
           pqtec , ptm1 , ptte , pum1 , pverv , pvm1 , pvol , pvom ,  &
           pxim1 , pxite , pxlm1 , pxlte , pxtec , pmfu , zlude , pccn
    real(rkx) , dimension(kbdim) :: pqhfla , pshfla , prsfc , &
           pssfc , ptopmax
    integer(ik4) , dimension(kbdim) , intent(out) :: kctop , kcbot
    real(rkx) , dimension(kbdim,klev,ktrac) :: pxtm1 , pxtte
    intent (in) papp1 , pqm1 , ptm1 , pum1 , pvm1 , pxim1 , &
           pxlm1 , pxtm1 , pccn
    intent(out) :: zlude
    intent (inout) ptopmax
    integer(ik4) , dimension(kbdim) :: itopec2
    integer(ik4) :: ilevmin , it , jk , jl , jt
    logical , dimension(kbdim) :: locum
    real(rkx) , dimension(kbdim,klev) :: zlu , zmfd , zqp1 , zqsat , zqu , &
      zqude , ztp1 , ztu , ztvp1 , zup1 , zvp1 , zxp1
    real(rkx) , dimension(kbdim) :: zrain , ztopmax
    real(rkx) :: zxip1 , zxlp1
    integer(ik4) , dimension(kbdim) :: kbotsc
    logical , dimension(kbdim) :: ldsc
    real(rkx) , dimension(kbdim,klev,ktrac) :: zxtp1 , zxtu
    real(rkx) , dimension(kbdim,klev+1) :: pmflxs
    real(rkx) , dimension(kbdim,klev) :: pmfude_rate , pmfdde_rate
    real(rkx) , dimension(kbdim) :: pcape
    real(rkx) , dimension(kbdim,klev+1) :: pqhfl , pahfs

    lookupoverflow = .false.
    !
    !---------------------------------------
    ! 1. CALCULATE T,Q AND QS AT MAIN LEVELS
    ! --------------------------------------
    !
    if ( iconv /= 4 ) then
      do jk = 1 , klev
        do jl = 1 , kproma
          ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*dt
          zqp1(jl,jk) = max(1.0e-8_rkx,pqm1(jl,jk) + pqte(jl,jk)*dt)
          zxlp1 = max(d_zero,pxlm1(jl,jk) + pxlte(jl,jk)*dt)
          zxip1 = max(d_zero,pxim1(jl,jk) + pxite(jl,jk)*dt)
          zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*dt
          zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*dt
          zxp1(jl,jk) = max(d_zero,zxlp1+zxip1)
          it = int(ztp1(jl,jk)*d_1000)
          if ( it < jptlucu1 .or. it > jptlucu2 ) then
            write(stderr,'(a,f12.4,a,i8)') &
              '! LOOKUP PROBLEM FOR T = ',real(ztp1(jl,jk)), &
              ' at ', rcmtimer%str( )
            lookupoverflow = .true.
          end if
          it = max(min(it,jptlucu2),jptlucu1)
          zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
          zqsat(jl,jk) = min(qsmax,zqsat(jl,jk))
          zqsat(jl,jk) = zqsat(jl,jk)/(d_one-ep1*zqsat(jl,jk))
          ztvp1(jl,jk) = ztp1(jl,jk)*d_one+ep1*(zqp1(jl,jk)-zxp1(jl,jk))
        end do
      end do
      if ( lookupoverflow ) then
        call fatal(__FILE__,__LINE__, &
                   'Cumulus Tables lookup error: OVERFLOW')
      end if
    else
      do jk = 1 , klev
        do jl = 1 , kproma
          ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*dt
          zqp1(jl,jk) = max(1.0e-8_rkx,pqm1(jl,jk) + pqte(jl,jk)*dt)
          zxlp1 = max(d_zero,pxlm1(jl,jk) + pxlte(jl,jk)*dt)
          zxip1 = max(d_zero,pxim1(jl,jk) + pxite(jl,jk)*dt)
          zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*dt
          zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*dt
          zxp1(jl,jk) = max(d_zero,zxlp1+zxip1)
        end do
      end do
    end if

    do jt = 1 , ktrac
      do jk = 1 , klev
        do jl = 1 , kproma
          zxtp1(jl,jk,jt) = max(d_zero,pxtm1(jl,jk,jt) + pxtte(jl,jk,jt)*dt)
        end do
      end do
    end do

    do jl = 1 , kproma
      zrain(jl) = d_zero
      locum(jl) = .false.
    end do
    !
    !---------------------------------------------------------------
    ! 2. CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
    ! --------------------------------------------------------------
    !
    select case (iconv)
    case (1)
      call cumastr(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,    &
                   zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,pxtte,&
                   pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,pvom,pvol,&
                   prsfc,pssfc,pxtec,pqtec,zqude,locum,ktype,kcbot,   &
                   kctop,ztu,zqu,zlu,zlude,pmfu,zmfd,zrain)
    case (2)
      call cumastrt(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                    zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                    pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                    pvom,pvol,prsfc,pssfc,pxtec,pqtec,zqude,locum,    &
                    ktype,kcbot,kctop,ztu,zqu,zlu,zlude,pmfu,zmfd,zrain)
    case (3)
      call cumastrh(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                    zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                    pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                    pvom,pvol,prsfc,pssfc,pxtec,pqtec,zqude,locum,    &
                    ktype,kcbot,kctop,ztu,zqu,zlu,zlude,pmfu,zmfd,zrain)
    case (4)
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
      ! ### for this test case we just specify an arbitrary decrease of
      ! turb. fluxes with height ###
      !do jk = klev+1 , 1 , -1
      !  pqhfl(:,jk) = pqhfla(:) * exp(-egrav*(d_one-sigma(jk)))
      !  pahfs(:,jk) = pshfla(:) * exp(-egrav*(d_one-sigma(jk)))
      !end do
      !
      ! Set to zero BL fluxes: we call BL after CU
      !
      pqhfl = d_zero
      pahfs = d_zero
      !
      ! Reset to zero output precipitation fluxes
      !
      pmflxr = d_zero
      pmflxs = d_zero
      !
      ! Call scheme
      !
      call ntiedtke(1,kproma,kbdim,klev,ldland,ztp1,zqp1,          &
                    zup1,zvp1,zxp1,pverv,pqhfl,pahfs,papp1,paphp1, &
                    pgeo,pgeoh,ptte,pqte,pvom,pvol,pxlte,pxite,    &
                    locum,ktype,kcbot,kctop,kbotsc,ldsc,ztu,zqu,   &
                    zlu,pmflxr,pmflxs,zrain,pmfu,zmfd,zlude,       &
                    pmfude_rate,pmfdde_rate,pcape,ktrac,pxtm1,     &
                    pxtte,pccn)
      prsfc = pmflxr(:,klev+1)*1.e3_rkx
      pssfc = pmflxs(:,klev+1)*1.e3_rkx
    case default
      call fatal(__FILE__,__LINE__, &
                 'ICONV must be in the range 1-4')
    end select
    !
    ! ----------------------------------------------
    ! 3. PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
    ! ----------------------------------------------
    !
    if ( iconv /= 4 ) then
      ilevmin = klev - 4
      do jl = 1 , kproma
        itopec2(jl) = klevp1
      end do
      do jk = 1 , ilevmin
        do jl = 1 , kproma
          if ( ilab(jl,jk) == 2 .and. itopec2(jl) == klevp1 ) itopec2(jl) = jk
        end do
      end do
      ztopmax(1:kproma) = ptopmax(1:kproma)
      do jl = 1 , kproma
        if ( itopec2(jl) == 1 ) then
          ptopmax(jl) = papp1(jl,1)
        else if ( itopec2(jl) /= klevp1 ) then
          ptopmax(jl) = paphp1(jl,itopec2(jl))
        else
          ptopmax(jl) = 99999.0_rkx
        end if
        ptopmax(jl) = min(ptopmax(jl),ztopmax(jl))
      end do
    end if
  end subroutine cucall
  !
  ! CUMASTR  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
  !
  ! THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
  ! PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
  ! PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
  ! PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
  ! SATURATED CUMULUS DOWNDRAFTS.
  !
  ! THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
  ! T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
  ! IT RETURNS ITS OUTPUT TO THE SAME SPACE
  !   1.MODIFIED TENDENCIES OF MODEL VARIABLES
  !   2.RATES OF CONVECTIVE PRECIPITATION (USED IN SURF SCHEME)
  !
  ! PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
  !   (1) DEFINE CONSTANTS AND PARAMETERS
  !   (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
  !       INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
  !   (3) CALCULATE CLOUD BASE IN 'CUBASE'
  !       AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
  !   (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
  !   (5) DO DOWNDRAFT CALCULATIONS:
  !       (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
  !       (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
  !       (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
  !           EFFECT OF CU-DOWNDRAFTS
  !   (6) DO FINAL CLOUD ASCENT IN 'CUASC'
  !   (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
  !       DO EVAPORATION IN SUBCLOUD LAYER
  !   (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
  !   (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
  !
  !   CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
  !   CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
  !   CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
  !   CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
  !   CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
  !   CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
  !   CUDQDT: UPDATES TENDENCIES FOR T AND Q
  !   CUDUDV: UPDATES TENDENCIES FOR U AND V
  !
  !   LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
  !   LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
  !   LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
  !   LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
  !   LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
  !
  !   MODEL PARAMETERS
  !   ------------------------------------------------
  !   ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
  !   ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
  !   ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
  !   ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
  !   CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
  !   CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
  !   CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
  !   CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  !   CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
  !
  !   REFERENCE.
  !   ----------
  !
  !   PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
  !
  subroutine cumastr(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,pqen,&
                     pxen,puen,pven,ptven,ktrac,ldland,pxten,pxtu,  &
                     pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,ptte,pqte,&
                     pvom,pvol,prsfc,pssfc,pxtec,pqtec,pqude,ldcum, &
                     ktype,kcbot,kctop,ptu,pqu,plu,plude,pmfu,pmfd,prain)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , &
                                 klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim,klev) :: ilab
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , ktype
    logical , dimension(kbdim) :: ldcum , ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: pqhfla , prain , prsfc , pssfc
    real(rkx) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
           pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
           ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
           pxen , pxtec
    real(rkx) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
    intent (in) pqhfla
    intent (inout) ktype , ldcum , pmfd
    integer(ik4) , dimension(kbdim) :: ictop0 , idtop , ihmin , ilwmin
    integer(ik4) :: ikb , it , it1 , itopm2 , jk , jl , jt
    logical :: llo1 , lo
    logical , dimension(kbdim) :: loddraf
    real(rkx) :: zalvdcp , zalvs , zb , zbi , zcons2 , zcor , zdepth ,  &
               zdhdz , zdqmin , zdqsdt , zdz , zeps , zes , zfac ,    &
               zgam , zhhat , zhsat , zmfmax , zpbmpt , zqalv ,       &
               zqsat , zqst1 , qumqe , zrh , zro , ztau , zzz , entrpen
    real(rkx) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                  zhcbase , zheat , zhmin , zmfub ,   &
                                  zmfub1 , zrfl , zsfl
    real(rkx) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
           zdmfup , zdpmel , zgeoh , zhhatt , zmfdq , zmfds , zmful , &
           zmfuq , zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , &
           zuu , zvd , zvu , zxenh
    real(rkx) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
           zxtenh
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cumastr'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lookupoverflow = .false.
    !
    !------------------------------------
    ! 1. SPECIFY CONSTANTS AND PARAMETERS
    ! -----------------------------------
    !
    zcons2 = d_one/(egrav*dtcum)

    ! *AMT* NOTE!
    ! this paramter is the CAPE adjustment timescale which in the global model
    ! was a function of horizontal resolution
    ! (nn wavenumber of a spectral model)
    ! this is translated roughly into horizontal resolution in meters
    !
    ztau = min(cmtcape,453600._rkx/nskmax)
    !
    !--------------------------------------------------------
    ! 2. INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
    ! -------------------------------------------------------
    !
    call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
               puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
               zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
               zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
               zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
               zdpmel,plu,plude,pqude,ilab)
    !
    !----------------------------
    ! 3.0 CLOUD BASE CALCULATIONS
    ! ---------------------------
    !
    ! (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
    ! -------------------------------------------
    !
    call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
                paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
                kcbot,ilab)
    !
    ! (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
    !     THEN DECIDE ON TYPE OF CUMULUS CONVECTION
    ! ---------------------------------------------
    !
    jk = 1
    do jl = 1 , kproma
      zdqpbl(jl) = 0.00_rkx
      zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
      idtop(jl) = 0
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        if ( jk >= kcbot(jl) ) then
          zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
        end if
      end do
    end do
    !
    ! (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
    ! AND DETERMINE CLOUD BASE MASSFLUX IGNORING
    ! THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
    ! ------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      qumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
      zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
      llo1 = zdqpbl(jl) > d_zero .and. qumqe > zdqmin .and. ldcum(jl)
      zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(qumqe,zdqmin)),0.010_rkx,llo1)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub(jl) = min(zmfub(jl),zmfmax)
      if ( .not.llo1 ) ldcum(jl) = .false.
      ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,ctrigger*pqhfla(jl)*egrav))
      entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(jl))
      zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    end do
    !
    !------------------------------------------------
    ! 4.0 DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
    ! -----------------------------------------------
    !
    ! (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
    ! CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
    ! FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
    ! -------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      zalvs = merge(wlhv,wlhs,ptu(jl,ikb) > tzero)
      zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
      ictop0(jl) = kcbot(jl) - 1
    end do
    do jk = klevm1 , 3 , -1
      do jl = 1 , kproma
        zalvs = merge(wlhv,wlhs,ztenh(jl,jk) > tzero)
        zalvdcp = zalvs/zcpcu(jl,jk)
        zqalv = d_one/zalvs
        zhsat = zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = nint(ztenh(jl,jk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/paphp1(jl,jk)
        zes = min(qsmax,zes)
        lo = zes < 0.40_rkx
        zcor = d_one/(d_one-ep1*zes)
        zqsat = zes*zcor
        it1 = it + 1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/paphp1(jl,jk)
        zqst1 = min(qsmax,zqst1)
        zqst1 = zqst1/(d_one-ep1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
        zzz = zcpcu(jl,jk)*ztenh(jl,jk)*ep1
        zhhat = zhsat - (zzz+zgam*zzz)/(d_one+zgam*zzz*zqalv) * &
                max(zqsenh(jl,jk)-zqenh(jl,jk),d_zero)
        zhhatt(jl,jk) = zhhat
        if ( jk < ictop0(jl) .and. zhcbase(jl) > zhhat ) ictop0(jl) = jk
      end do
    end do
    !
    ! DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
    ! (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
    !
    do jl = 1 , kproma
      ktype(jl) = merge(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl))>2.0e4_rkx)
      entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(jl))
      zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    end do

    if ( lookupoverflow ) then
      call fatal(__FILE__,__LINE__, &
                 'Cumulus Tables lookup error: OVERFLOW')
    end if
    !
    ! FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
    ! -------------------------------------------
    !
    do jl = 1 , kproma
      zhmin(jl) = d_zero
      ihmin(jl) = 0
      llo1 = ldcum(jl) .and. ktype(jl) == 1
      if ( llo1 ) then
        ikb = kcbot(jl)
        ihmin(jl) = ikb
      end if
    end do

    zb = 25.0_rkx
    zbi = d_one/(zb*egrav)
    do jk = klev , 1 , -1
      do jl = 1 , kproma
        llo1 = ldcum(jl) .and. ktype(jl) == 1 .and. ihmin(jl) == kcbot(jl)
        if ( llo1 .and. jk < kcbot(jl) .and. jk >= ictop0(jl) ) then
          zalvs = merge(wlhv,wlhs,ztenh(jl,jk) > tzero)
          ikb = kcbot(jl)
          zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+ep1*zqenh(jl,jk)))
          zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
          zdhdz = (zcpen(jl,jk-1)*pten(jl,jk-1)-zcpen(jl,jk)      &
                  *pten(jl,jk)+zalvs*(pqen(jl,jk-1)-pqen(jl,jk))  &
                  +(pgeo(jl,jk-1)-pgeo(jl,jk)))                   &
                  *egrav/(pgeo(jl,jk-1)-pgeo(jl,jk))
          zdepth = zgeoh(jl,jk) - zgeoh(jl,ikb)
          zfac = sqrt(d_one+zdepth*zbi)
          zhmin(jl) = zhmin(jl) + zdhdz*zfac*zdz
          zrh = -zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
          if ( zhmin(jl) > zrh ) ihmin(jl) = jk
        end if
      end do
    end do

    do jl = 1 , kproma
      if ( ldcum(jl) .and. ktype(jl) == 1 ) then
        if ( ihmin(jl) < ictop0(jl) ) ihmin(jl) = ictop0(jl)
      end if
    end do
    !
    ! (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
    ! ------------------------------------------------
    !
    call cuasc(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven, &
               ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,   &
               zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,ilab,&
               ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,zmful,&
               plude,pqude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,zcpen,  &
               zcpcu,kcbot,kctop,ictop0)
    !
    ! (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
    ! CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
    ! -------------------------------------------------------------
    !
    do jl = 1 , kproma
      zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
      if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.e4_rkx ) then
        ktype(jl) = 2
      end if
      if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
      if ( ktype(jl) == 2 ) zentr(jl) = entrscv
      zrfl(jl) = zdmfup(jl,1)
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zrfl(jl) = zrfl(jl) + zdmfup(jl,jk)
      end do
    end do
    !
    !--------------------------------------------
    ! 5.0          CUMULUS DOWNDRAFT CALCULATIONS
    ! -------------------------------------------
    !
    if ( lmfdd ) then
      !
      ! (A) DETERMINE LFS IN 'CUDLFS'
      ! -----------------------------
      !
      call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                  ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                  zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                  zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
      !
      ! (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
      ! ----------------------------------------------------
      !
      call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                   ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                   zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)

    end if
    !
    ! 5.5 RECALCULATE CLOUD BASE MASSFLUX FROM A
    !     CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
    !     AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
    !     ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
    ! ------------------------------------------------
    !
    do jl = 1 , kproma
      zheat(jl) = d_zero
      zcape(jl) = d_zero
      zmfub1(jl) = zmfub(jl)
    end do

    do jk = 1 , klev
      do jl = 1 , kproma
        llo1 = ldcum(jl) .and. ktype(jl) == 1
        if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
          ikb = kcbot(jl)
          zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+ep1*zqenh(jl,jk)))
          zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
          zheat(jl) = zheat(jl)                                    &
                      + ((pten(jl,jk-1)-pten(jl,jk)                &
                      +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)        &
                      +ep1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                      *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
          zcape(jl) = zcape(jl)                                       &
                      + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) &
                      +egrav*ep1*(pqu(jl,jk)-zqenh(jl,jk))            &
                      -egrav*plu(jl,jk))*zdz
        end if
      end do
    end do

    do jl = 1 , kproma
      if ( ldcum(jl) .and. ktype(jl) == 1 ) then
        ikb = kcbot(jl)
        zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
        zmfub1(jl) = max(zmfub1(jl),0.0010_rkx)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl) = min(zmfub1(jl),zmfmax)
      end if
    end do
    !
    ! RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
    ! DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
    ! FOR SHALLOW CONVECTION (KTYPE=2)
    ! ----------------------------------------------
    !
    do jl = 1 , kproma
      if ( ktype(jl) == 2 ) then
        ikb = kcbot(jl)
        llo1 = pmfd(jl,ikb) < d_zero .and. loddraf(jl)
        zeps = merge(cmfdeps,d_zero,llo1)
        qumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb) &
                 - (d_one-zeps)*zqenh(jl,ikb)
        zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1 = zdqpbl(jl) > d_zero .and. qumqe > zdqmin .and. &
               ldcum(jl) .and. zmfub(jl) < zmfmax
        zmfub1(jl) = merge(zdqpbl(jl)/(egrav*            &
                     max(qumqe,zdqmin)),zmfub(jl),llo1)
        zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl)  &
                     -zmfub(jl)) < 0.20_rkx*zmfub(jl))
      end if
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
          pmfd(jl,jk) = pmfd(jl,jk)*zfac
          zmfds(jl,jk) = zmfds(jl,jk)*zfac
          zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
            zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
          end if
        end do
      end do

    end do
    !
    ! NEW VALUES OF CLOUD BASE MASS FLUX
    ! ----------------------------------
    !
    do jl = 1 , kproma
      if ( ldcum(jl) ) zmfub(jl) = zmfub1(jl)
    end do
    !
    !------------------------------------------------------
    ! 6.0 DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
    !     FOR PENETRATIVE CONVECTION (TYPE=1),
    !     FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
    !     AND FOR MID-LEVEL CONVECTION (TYPE=3).
    ! -----------------------------------------------------
    !
    call cuasc(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven, &
               ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,   &
               zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,ilab,&
               ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,zmful,&
               plude,pqude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,zcpen,  &
               zcpcu,kcbot,kctop,ictop0)
    !
    !-------------------------------------------------
    ! 7.0 DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
    ! ------------------------------------------------
    !
    call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
               zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
               ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
               zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
               itopm2)
    !
    !-------------------------------------------------------
    ! 8.0 UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
    ! ------------------------------------------------------
    !
    call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
                pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
                zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
                zsfl,zcpen,pqtec,pqude,prsfc,pssfc)

    !
    !-------------------------------------------------------
    ! 9.0 UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
    ! ------------------------------------------------------
    !
    if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                               kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                               zuu,zud,zvu,zvd,pmfu,pmfd)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cumastr
  !
  ! THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
  ! PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
  ! PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
  ! PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
  ! SATURATED CUMULUS DOWNDRAFTS.
  !
  ! *CUMASTRH* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 3
  ! THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
  ! T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
  ! IT RETURNS ITS OUTPUT TO THE SAME SPACE
  !   1.MODIFIED TENDENCIES OF MODEL VARIABLES
  !   2.RATES OF CONVECTIVE PRECIPITATION (USED IN SURF SCHEME)
  !
  ! METHOD
  ! -------
  !
  ! PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
  !   (1) DEFINE CONSTANTS AND PARAMETERS
  !   (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
  !       INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
  !   (3) CALCULATE CLOUD BASE IN 'CUBASE'
  !       AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
  !   (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
  !   (5) DO DOWNDRAFT CALCULATIONS:
  !       (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
  !       (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
  !       (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
  !           EFFECT OF CU-DOWNDRAFTS
  !   (6) DO FINAL CLOUD ASCENT IN 'CUASC'
  !   (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
  !       DO EVAPORATION IN SUBCLOUD LAYER
  !   (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
  !   (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
  !
  !   CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
  !   CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
  !   CUASCT:  CLOUD ASCENT FOR ENTRAINING PLUME
  !   CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
  !   CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
  !   CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
  !   CUDQDT: UPDATES TENDENCIES FOR T AND Q
  !   CUDUDV: UPDATES TENDENCIES FOR U AND V
  !
  !   LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
  !   LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
  !   LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
  !   LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
  !   LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
  !
  !   MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
  !   ----------------------------------------------------
  !   ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
  !   ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
  !   ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
  !   ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
  !   CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
  !   CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
  !   CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
  !   CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  !   CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
  !
  !   REFERENCE.
  !   ----------
  !
  !   PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
  !
  subroutine cumastrh(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,    &
                      pqen,pxen,puen,pven,ptven,ktrac,ldland,pxten, &
                      pxtu,pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,    &
                      ptte,pqte,pvom,pvol,prsfc,pssfc,pxtec,pqtec,  &
                      pqude,ldcum,ktype,kcbot,kctop,ptu,pqu,plu,    &
                      plude,pmfu,pmfd,prain)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , &
                                 klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim,klev) :: ilab
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , ktype
    logical , dimension(kbdim) :: ldcum , ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: pqhfla , prain , prsfc , pssfc
    real(rkx) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
           pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
           ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
           pxen , pxtec
    real(rkx) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
    intent (in) pqhfla
    intent (inout) ktype , ldcum , pmfd
    integer(ik4) , dimension(kbdim) :: ictop0 , idtop , ilwmin
    integer(ik4) :: ikb , it , it1 , itopm2 , jk , jl , jt
    logical :: llo1 , lo
    logical , dimension(kbdim) :: loddraf
    real(rkx) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
               zdz , zeps , zes , zfac , zgam , zhhat , zhsat ,       &
               zmfmax , zpbmpt , zqalv , zqsat , zqst1 , zqumqe ,     &
               zro , ztau , zzz , entrpen
    real(rkx) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                  zhcbase , zheat , zmfub , zmfub1 ,  &
                                  zrfl , zsfl
    real(rkx) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
           zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
           zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
           zvd , zvu , zxenh
    real(rkx) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
           zxtenh
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cumastrh'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lookupoverflow = .false.
    !
    !------------------------------------
    ! 1. SPECIFY CONSTANTS AND PARAMETERS
    ! -----------------------------------
    !
    zcons2 = d_one/(egrav*dtcum)

    ! *AMT* NOTE!
    ! this paramter is the CAPE adjustment timescale which in the global model
    ! was a function of horizontal resolution
    ! (nn wavenumber of a spectral model)
    ! this is translated roughly into horizontal resolution in meters
    !
    ztau = min(cmtcape,453600._rkx/nskmax)
    !
    !--------------------------------------------------------
    ! 2. INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
    ! -------------------------------------------------------
    !
    call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
               puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
               zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
               zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
               zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
               zdpmel,plu,plude,pqude,ilab)
    !
    !----------------------------
    ! 3.0 CLOUD BASE CALCULATIONS
    ! ---------------------------
    !
    ! (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
    ! -------------------------------------------
    !
    call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
                paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
                kcbot,ilab)
    !
    ! (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
    !     THEN DECIDE ON TYPE OF CUMULUS CONVECTION
    ! ---------------------------------------------
    !
    jk = 1
    do jl = 1 , kproma
      zdqpbl(jl) = 0.00_rkx
      zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
      idtop(jl) = 0
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)                           &
                    *(paphp1(jl,jk+1)-paphp1(jl,jk))
        if ( jk >= kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)  &
             *(paphp1(jl,jk+1)-paphp1(jl,jk))
      end do
    end do
    !
    ! (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
    !     AND DETERMINE CLOUD BASE MASSFLUX IGNORING
    !     THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
    ! ------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
      zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
      llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl)
      zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010_rkx,llo1)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub(jl) = min(zmfub(jl),zmfmax)
      if ( .not.llo1 ) ldcum(jl) = .false.
      ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,ctrigger*pqhfla(jl)*egrav))
      entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(jl))
      zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    end do
    !
    !------------------------------------------------
    ! 4.0 DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
    ! -----------------------------------------------
    !
    ! (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
    !     CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
    !     FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
    ! -----------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      zalvs = merge(wlhv,wlhs,ptu(jl,ikb) > tzero)
      zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb)+zgeoh(jl,ikb)+zalvs*pqu(jl,ikb)
      ictop0(jl) = kcbot(jl) - 1
    end do
    do jk = klevm1 , 3 , -1
      do jl = 1 , kproma
        zalvs = merge(wlhv,wlhs,ztenh(jl,jk) > tzero)
        zalvdcp = zalvs/zcpcu(jl,jk)
        zqalv = d_one/zalvs
        zhsat = zcpcu(jl,jk)*ztenh(jl,jk)+zgeoh(jl,jk)+zalvs*zqsenh(jl,jk)
        it = nint(ztenh(jl,jk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/paphp1(jl,jk)
        zes = min(qsmax,zes)
        lo = zes < 0.40_rkx
        zcor = d_one/(d_one-ep1*zes)
        zqsat = zes*zcor
        it1 = it + 1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/paphp1(jl,jk)
        zqst1 = min(qsmax,zqst1)
        zqst1 = zqst1/(d_one-ep1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
        zzz = zcpcu(jl,jk)*ztenh(jl,jk)*ep1
        zhhat = zhsat - (zzz+zgam*zzz)/(d_one+zgam*zzz*zqalv)* &
                max(zqsenh(jl,jk)-zqenh(jl,jk),d_zero)
        if ( jk < ictop0(jl) .and. zhcbase(jl) > zhhat ) ictop0(jl) = jk
      end do
    end do

    if ( lookupoverflow ) then
      call fatal(__FILE__,__LINE__, &
                 'Cumulus Tables lookup error: OVERFLOW')
    end if
    !!
    !! DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
    !! (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
    !!
    !  do jl = 1 , kproma
    !    ktype(jl) = merge(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl))>2.e4_rkx)
    !    entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(jl))
    !    zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    !  end do
    !!
    ! (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
    ! --------------------------------------------------
    !
    call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
                ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
                zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
                ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
                zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
                ictop0)
    !
    ! (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
    !     CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
    ! -------------------------------------------------------------
    !
    do jl = 1 , kproma
      zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
      if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.e4_rkx ) then
        ktype(jl) = 2
      end if
      if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
      if ( ktype(jl) == 2 ) zentr(jl) = entrscv
      zrfl(jl) = zdmfup(jl,1)
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zrfl(jl) = zrfl(jl) + zdmfup(jl,jk)
      end do
    end do
    !
    !-----------------------------------
    ! 5.0 CUMULUS DOWNDRAFT CALCULATIONS
    ! ----------------------------------
    !
    if ( lmfdd ) then
      !
      ! (A) DETERMINE LFS IN 'CUDLFS'
      ! -----------------------------
      !
      call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                  ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                  zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                  zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
      !
      ! (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
      ! ----------------------------------------------------
      !
      call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                   ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                   zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)

    end if
    !
    ! 5.5 RECALCULATE CLOUD BASE MASSFLUX FROM A
    !     CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
    !     AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
    !     ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
    ! ------------------------------------------------
    !
    do jl = 1 , kproma
      zheat(jl) = d_zero
      zcape(jl) = d_zero
      zmfub1(jl) = zmfub(jl)
    end do

    do jk = 1 , klev
      do jl = 1 , kproma
        llo1 = ldcum(jl) .and. ktype(jl) == 1
        if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
          ikb = kcbot(jl)
          zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+ep1*zqenh(jl,jk)))
          zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
          zheat(jl) = zheat(jl)                                    &
                      + ((pten(jl,jk-1)-pten(jl,jk)                &
                      +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)        &
                      +ep1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                      *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
          zcape(jl) = zcape(jl)                                       &
                      + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) &
                      +egrav*ep1*(pqu(jl,jk)-zqenh(jl,jk))            &
                      -egrav*plu(jl,jk))*zdz
        end if
      end do
    end do

    do jl = 1 , kproma
      if ( ldcum(jl) .and. ktype(jl) == 1 ) then
        ikb = kcbot(jl)
        zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
        zmfub1(jl) = max(zmfub1(jl),0.0010_rkx)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        zmfub1(jl) = min(zmfub1(jl),zmfmax)
      end if
    end do
    !
    ! RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
    ! DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
    ! FOR SHALLOW CONVECTION (KTYPE=2)
    ! ----------------------------------------------
    !
    do jl = 1 , kproma
      if ( ktype(jl) == 2 ) then
        ikb = kcbot(jl)
        llo1 = pmfd(jl,ikb) < d_zero .and. loddraf(jl)
        zeps = merge(cmfdeps,d_zero,llo1)
        zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)  &
                 - (d_one-zeps)*zqenh(jl,ikb)
        zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. &
               ldcum(jl) .and. zmfub(jl) < zmfmax
        zmfub1(jl) = merge(zdqpbl(jl)/(egrav*                  &
               max(zqumqe,zdqmin)),zmfub(jl),llo1)
        zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl) &
                     -zmfub(jl)) < 0.20_rkx*zmfub(jl))
      end if
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
          pmfd(jl,jk) = pmfd(jl,jk)*zfac
          zmfds(jl,jk) = zmfds(jl,jk)*zfac
          zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
            zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
          end if
        end do
      end do

    end do
    !
    ! NEW VALUES OF CLOUD BASE MASS FLUX
    ! ----------------------------------
    !
    do jl = 1 , kproma
      if ( ldcum(jl) ) zmfub(jl) = zmfub1(jl)
    end do
    !
    !------------------------------------------------------
    ! 6.0 DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
    !     FOR PENETRATIVE CONVECTION (TYPE=1),
    !     FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
    !     AND FOR MID-LEVEL CONVECTION (TYPE=3).
    ! -----------------------------------------------------
    !
    call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
                ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
                zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
                ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
                zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
                ictop0)
    !
    !-------------------------------------------------
    ! 7.0 DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
    ! ------------------------------------------------
    !
    call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
               zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
               ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
               zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
               itopm2)
    !
    !-------------------------------------------------------
    ! 8.0 UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
    ! ------------------------------------------------------
    !
    call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
                pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
                zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
                zsfl,zcpen,pqtec,pqude,prsfc,pssfc)
    !
    !-------------------------------------------------------
    ! 9.0 UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
    ! ------------------------------------------------------
    !
    if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                               kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                               zuu,zud,zvu,zvd,pmfu,pmfd)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cumastrh
  !
  !  THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
  !  PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
  !  PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
  !  PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
  !  SATURATED CUMULUS DOWNDRAFTS.
  !
  !  *CUMASTRT* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 2
  !  THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
  !  T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
  !  IT RETURNS ITS OUTPUT TO THE SAME SPACE
  !    1.MODIFIED TENDENCIES OF MODEL VARIABLES
  !    2.RATES OF CONVECTIVE PRECIPITATION (USED IN SURF SCHEME)
  !
  !  PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
  !    (1) DEFINE CONSTANTS AND PARAMETERS
  !    (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
  !        INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
  !    (3) CALCULATE CLOUD BASE IN 'CUBASE'
  !        AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
  !    (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
  !    (5) DO DOWNDRAFT CALCULATIONS:
  !        (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
  !        (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
  !        (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
  !            EFFECT OF CU-DOWNDRAFTS
  !    (6) DO FINAL CLOUD ASCENT IN 'CUASC'
  !    (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
  !        DO EVAPORATION IN SUBCLOUD LAYER
  !    (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
  !    (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
  !
  ! CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
  ! CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
  ! CUASCT: CLOUD ASCENT FOR ENTRAINING PLUME
  ! CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
  ! CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
  ! CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
  ! CUDQDT: UPDATES TENDENCIES FOR T AND Q
  ! CUDUDV: UPDATES TENDENCIES FOR U AND V
  !
  ! LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
  ! LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
  ! LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
  ! LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
  ! LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
  !
  ! ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
  ! ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
  ! ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
  ! ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
  ! CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
  ! CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
  ! CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
  ! CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
  ! CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
  !
  ! REFERENCE.
  ! ----------
  !
  ! PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
  !
  subroutine cumastrt(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,    &
                      pqen,pxen,puen,pven,ptven,ktrac,ldland,pxten, &
                      pxtu,pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,    &
                      ptte,pqte,pvom,pvol,prsfc,pssfc,pxtec,pqtec,  &
                      pqude,ldcum,ktype,kcbot,kctop,ptu,pqu,plu,    &
                      plude,pmfu,pmfd,prain)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim,klev) :: ilab
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , ktype
    logical , dimension(kbdim) :: ldcum , ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: pqhfla , prain , prsfc , pssfc
    real(rkx) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
           pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
           ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
           pxen , pxtec
    real(rkx) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
    intent (in) pqhfla
    intent (inout) ktype , ldcum , pmfd

    integer(ik4) , dimension(kbdim) :: ictop0 , idtop , ilwmin
    integer(ik4) :: ikb , it , it1 , itopm2 , jk , jl , jt
    logical :: llo1 , lo
    logical , dimension(kbdim) :: loddraf
    real(rkx) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
               zeps , zes , zfac , zgam , zhhat , zhsat , zmfmax ,    &
               zpbmpt , zqalv , zqsat , zqst1 , zqumqe , zzz , entrpen
    real(rkx) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
           zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
           zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
           zvd , zvu , zxenh
    real(rkx) , dimension(kbdim) :: zdqcv , zdqpbl , zentr , zhcbase ,  &
                                  zmfub , zmfub1 , zrfl , zsfl
    real(rkx) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
           zxtenh
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cumastrt'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lookupoverflow = .false.
    !
    !------------------------------------
    ! 1. SPECIFY CONSTANTS AND PARAMETERS
    ! -----------------------------------
    !
    !
    zcons2 = d_one/(egrav*dtcum)
    !
    !--------------------------------------------------------
    ! 2. INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
    ! -------------------------------------------------------
    !
    call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
               puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
               zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
               zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
               zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
               zdpmel,plu,plude,pqude,ilab)
    !
    !----------------------------
    ! 3.0 CLOUD BASE CALCULATIONS
    ! ---------------------------
    !
    ! (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
    ! -------------------------------------------
    !
    call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
                paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
                kcbot,ilab)
    !
    ! (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
    !     THEN DECIDE ON TYPE OF CUMULUS CONVECTION
    ! ---------------------------------------------
    !
    jk = 1
    do jl = 1 , kproma
      zdqpbl(jl) = 0.00_rkx
      zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
      idtop(jl) = 0
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)                           &
                    *(paphp1(jl,jk+1)-paphp1(jl,jk))
        if ( jk >= kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)  &
             *(paphp1(jl,jk+1)-paphp1(jl,jk))
      end do
    end do
    !
    ! (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
    !     AND DETERMINE CLOUD BASE MASSFLUX IGNORING
    !     THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
    ! ------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
      zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
      llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl)
      zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010_rkx,llo1)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub(jl) = min(zmfub(jl),zmfmax)
      if ( .not.llo1 ) ldcum(jl) = .false.
      ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,ctrigger*pqhfla(jl)*egrav))
      entrpen = merge(entrpen_lnd,entrpen_ocn,ldland(jl))
      zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    end do
    !
    !------------------------------------------------
    ! 4.0 DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
    ! -----------------------------------------------
    !
    ! (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
    !     CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
    !     FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
    ! -----------------------------------------------------
    !
    do jl = 1 , kproma
      ikb = kcbot(jl)
      zalvs = merge(wlhv,wlhs,ptu(jl,ikb) > tzero)
      zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb)         &
                    + zalvs*pqu(jl,ikb)
      ictop0(jl) = kcbot(jl) - 1
    end do
    do jk = klevm1 , 3 , -1
      do jl = 1 , kproma
        zalvs = merge(wlhv,wlhs,ztenh(jl,jk) > tzero)
        zalvdcp = zalvs/zcpcu(jl,jk)
        zqalv = d_one/zalvs
        zhsat = zcpcu(jl,jk)*ztenh(jl,jk) + zgeoh(jl,jk)              &
                + zalvs*zqsenh(jl,jk)
        it = nint(ztenh(jl,jk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/paphp1(jl,jk)
        zes = min(qsmax,zes)
        lo = zes < 0.40_rkx
        zcor = d_one/(d_one-ep1*zes)
        zqsat = zes*zcor
        it1 = it + 1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/paphp1(jl,jk)
        zqst1 = min(qsmax,zqst1)
        zqst1 = zqst1/(d_one-ep1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
        zzz = zcpcu(jl,jk)*ztenh(jl,jk)*ep1
        zhhat = zhsat - (zzz+zgam*zzz)/(d_one+zgam*zzz*zqalv)         &
                *max(zqsenh(jl,jk)-zqenh(jl,jk),d_zero)
        if ( jk < ictop0(jl) .and. zhcbase(jl) > zhhat ) ictop0(jl) = jk
      end do
    end do

    if ( lookupoverflow ) then
      call fatal(__FILE__,__LINE__, &
                 'Cumulus Tables lookup error: OVERFLOW')
    end if
    !!
    !! DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
    !! (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
    !!
    ! do jl = 1 , kproma
    !   ktype(jl) = merge(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl))>2.e4_rkx)
    !   zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
    ! end do
    !!
    ! (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
    ! --------------------------------------------------
    !
    call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
                ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
                zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
                ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
                zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
                ictop0)
    !
    ! (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
    !     CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
    ! -------------------------------------------------------------
    !
    do jl = 1 , kproma
      zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
      if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.e4_rkx ) then
        ktype(jl) = 2
      end if
      if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
      if ( ktype(jl) == 2 ) zentr(jl) = entrscv
      zrfl(jl) = zdmfup(jl,1)
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zrfl(jl) = zrfl(jl) + zdmfup(jl,jk)
      end do
    end do
    !
    !-----------------------------------
    ! 5.0 CUMULUS DOWNDRAFT CALCULATIONS
    ! ----------------------------------
    !
    if ( lmfdd ) then
      !
      ! (A) DETERMINE LFS IN 'CUDLFS'
      ! -----------------------------
      !
      call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                  ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                  zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                  zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
      !
      ! (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
      ! ----------------------------------------------------
      !
      call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                   ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                   zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)
      !
      ! (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
      !      DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
      ! ---------------------------------------------------
      !
      do jl = 1 , kproma
        if ( loddraf(jl) ) then
          ikb = kcbot(jl)
          llo1 = pmfd(jl,ikb) < d_zero
          zeps = merge(cmfdeps,d_zero,llo1)
          zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)       &
                   - (d_one-zeps)*zqenh(jl,ikb)
          zdqmin = max(0.010_rkx*zqenh(jl,ikb),1.0e-10_rkx)
          zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
          llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl) &
                 .and. zmfub(jl) < zmfmax
          zmfub1(jl) = merge(zdqpbl(jl) / &
                        (egrav*max(zqumqe,zdqmin)),zmfub(jl),llo1)
          zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),(ktype(jl) == 1 .or.  &
                       ktype(jl) == 2) .and. &
                       abs(zmfub1(jl)-zmfub(jl)) < 0.20_rkx*zmfub(jl))
        end if
      end do
      do jk = 1 , klev
        do jl = 1 , kproma
          if ( loddraf(jl) ) then
            zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
            pmfd(jl,jk) = pmfd(jl,jk)*zfac
            zmfds(jl,jk) = zmfds(jl,jk)*zfac
            zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
            zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
          end if
        end do
        do jt = 1 , ktrac
          do jl = 1 , kproma
            if ( loddraf(jl) ) then
              zfac = zmfub1(jl)/max(zmfub(jl),1.0e-10_rkx)
              zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
            end if
          end do
        end do
      end do
      !
      ! NEW VALUES OF CLOUD BASE MASS FLUX
      ! ----------------------------------
      do jl = 1 , kproma
        if ( loddraf(jl) ) zmfub(jl) = zmfub1(jl)
      end do

    end if
    !
    !------------------------------------------------------
    ! 6.0 DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
    !     FOR PENETRATIVE CONVECTION (TYPE=1),
    !     FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
    !     AND FOR MID-LEVEL CONVECTION (TYPE=3).
    ! -----------------------------------------------------
    !
    call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
                ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
                zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
                ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
                zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
                ictop0)
    !
    !-------------------------------------------------
    ! 7.0 DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
    ! ------------------------------------------------
    !
    call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
               zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
               ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
               zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
               itopm2)
    !
    !-------------------------------------------------------
    ! 8.0 UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
    ! ------------------------------------------------------
    !
    call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
                pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
                zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
                zsfl,zcpen,pqtec,pqude,prsfc,pssfc)
    !
    !-------------------------------------------------------
    ! 9.0 UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
    ! ------------------------------------------------------
    !
    if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                               kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                               zuu,zud,zvu,zvd,pmfu,pmfd)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cumastrt
  !
  ! THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
  ! TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
  ! DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
  ! AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
  ! FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
  !
  subroutine cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen, &
                   pxen,puen,pven,ptven,ktrac,pxten,pxtenh,pxtu,    &
                   pxtd,pmfuxt,pmfdxt,pverv,pgeo,paphp1,pgeoh,ptenh,&
                   pqenh,pqsenh,pxenh,klwmin,ptu,pqu,ptd,pqd,puu,   &
                   pvu,pud,pvd,pmfu,pmfd,pmfus,pmfds,pmfuq,pmfdq,   &
                   pdmfup,pdmfdp,pcpen,pcpcu,pdpmel,plu,plude,pqude,&
                   klab)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , kproma , ktrac , klevm1 , klevp1
    integer(ik4) , dimension(kbdim,klev) :: klab
    integer(ik4) , dimension(kbdim) :: klwmin
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfdp ,       &
           pdmfup , pdpmel , pgeo , pgeoh , plu , plude , pmfd ,      &
           pmfdq , pmfds , pmfu , pmfuq , pmfus , pqd , pqen , pqenh ,&
           pqsen , pqsenh , pqu , pqude , ptd , pten , ptenh , ptu ,  &
           ptven , pud , puen , puu , pvd , pven , pverv , pvu ,      &
           pxen , pxenh
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtd , &
           pxten , pxtenh , pxtu
    intent (in) paphp1 , pgeo , pqen , pqsen , pten , ptven , puen , &
                pven , pverv , pxen , pxten
    intent (out) klab , klwmin , pdmfdp , pdmfup , pdpmel , plu ,     &
                 plude , pmfd , pmfdq , pmfds , pmfdxt , pmfu ,       &
                 pmfuq , pmfus , pmfuxt , pqd , pqu , pqude , ptd ,   &
                 ptu , pud , puu , pvd , pvu , pxenh , pxtd , pxtu
    intent (inout) pcpcu , pcpen , pgeoh , pqenh , pqsenh , ptenh ,   &
                   pxtenh
    integer(ik4) :: icall , ik , jk , jl , jt
    logical , dimension(kbdim) :: loflag
    real(rkx) :: zarg , zcpm , zzs
    real(rkx) , dimension(kbdim) :: zph , zwmax
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cumini'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !--------------------------------------------------
    ! 1. SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
    !    ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
    !    FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
    ! -------------------------------------------------
    !
    do jk = 1 , klev
      do jl = 1 , kproma
        ! pcpen(jl,jk) = cpd*(d_one+hcrm1*max(pqen(jl,jk),0.00_rkx))
        pcpen(jl,jk) = cpd
      end do
    end do
    do jl = 1 , kproma
      zarg = paphp1(jl,klevp1)/paphp1(jl,klev)
      pgeoh(jl,klev) = rgas*ptven(jl,klev)*log(zarg)
    end do
    do jk = klevm1 , 2 , -1
      do jl = 1 , kproma
        zarg = paphp1(jl,jk+1)/paphp1(jl,jk)
        pgeoh(jl,jk) = pgeoh(jl,jk+1) + rgas*ptven(jl,jk)*log(zarg)
      end do
    end do
    do jk = 2 , klev
      do jl = 1 , kproma
        zcpm = (pcpen(jl,jk)+pcpen(jl,jk-1))*d_half
        ptenh(jl,jk) = (max(pcpen(jl,jk-1)*pten(jl,jk-1)+ &
               pgeo(jl,jk-1),pcpen(jl,jk)*pten(jl,jk)+    &
               pgeo(jl,jk))-pgeoh(jl,jk))/zcpm
        pqsenh(jl,jk) = pqsen(jl,jk-1)
        zph(jl) = paphp1(jl,jk)
        loflag(jl) = .true.
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          pxtenh(jl,jk,jt) = (pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*d_half
        end do
      end do

      ik = jk
      icall = 0
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptenh,pqsenh,loflag,icall)

      do jl = 1 , kproma
        pxenh(jl,jk) = (pxen(jl,jk)+pxen(jl,jk-1))*d_half
        pqenh(jl,jk) = min(pqen(jl,jk-1),pqsen(jl,jk-1))              &
                       + (pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk) = max(pqenh(jl,jk),d_zero)
        ! pcpcu(jl,jk) = cpd*(d_one+hcrm1*pqenh(jl,jk))
        pcpcu(jl,jk) = cpd
      end do
    end do

    do jl = 1 , kproma
      ptenh(jl,klev) = (pcpen(jl,klev)*pten(jl,klev)+pgeo(jl,klev)    &
                       -pgeoh(jl,klev))/pcpen(jl,klev)
      pxenh(jl,klev) = pxen(jl,klev)
      pqenh(jl,klev) = pqen(jl,klev)
      pcpcu(jl,1) = pcpen(jl,1)
      ptenh(jl,1) = pten(jl,1)
      pxenh(jl,1) = pxen(jl,1)
      pqenh(jl,1) = pqen(jl,1)
      pgeoh(jl,1) = pgeo(jl,1)
      klwmin(jl) = klev
      zwmax(jl) = d_zero
    end do

    do jt = 1 , ktrac
      do jl = 1 , kproma
        pxtenh(jl,klev,jt) = pxten(jl,klev,jt)
        pxtenh(jl,1,jt) = pxten(jl,1,jt)
      end do
    end do

    do jk = klevm1 , 2 , -1
      do jl = 1 , kproma
        zzs = max(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),pcpcu(jl,jk+1)* &
               ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk) = (zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
      end do
    end do

    do jk = klev , 3 , -1
      do jl = 1 , kproma
        if ( pverv(jl,jk) < zwmax(jl) ) then
          zwmax(jl) = pverv(jl,jk)
          klwmin(jl) = jk
        end if
      end do
    end do
    !
    !--------------------------------------------------
    ! 2.0 INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
    ! -------------------------------------------------
    !
    do jk = 1 , klev
      ik = jk - 1
      if ( jk == 1 ) ik = 1
      do jl = 1 , kproma
        ptu(jl,jk) = ptenh(jl,jk)
        ptd(jl,jk) = ptenh(jl,jk)
        pqu(jl,jk) = pqenh(jl,jk)
        pqd(jl,jk) = pqenh(jl,jk)
        plu(jl,jk) = d_zero
        puu(jl,jk) = puen(jl,ik)
        pud(jl,jk) = puen(jl,ik)
        pvu(jl,jk) = pven(jl,ik)
        pvd(jl,jk) = pven(jl,ik)
        pmfu(jl,jk) = d_zero
        pmfd(jl,jk) = d_zero
        pmfus(jl,jk) = d_zero
        pmfds(jl,jk) = d_zero
        pmfuq(jl,jk) = d_zero
        pmfdq(jl,jk) = d_zero
        pdmfup(jl,jk) = d_zero
        pdmfdp(jl,jk) = d_zero
        pdpmel(jl,jk) = d_zero
        plude(jl,jk) = d_zero
        pqude(jl,jk) = d_zero
        klab(jl,jk) = 0
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          pxtu(jl,jk,jt) = pxtenh(jl,jk,jt)
          pxtd(jl,jk,jt) = pxtenh(jl,jk,jt)
          pmfuxt(jl,jk,jt) = d_zero
          pmfdxt(jl,jk,jt) = d_zero
        end do
      end do

    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuini
  !
  ! THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
  ! FOR CUMULUS PARAMETERIZATION
  ! TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
  ! (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
  ! FLUXES AS WELL AS PRECIPITATION RATES)
  !
  ! LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
  ! AND THEN CALCULATE MOIST ASCENT FOR
  ! ENTRAINING/DETRAINING PLUME.
  ! ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
  ! SHALLOW AND DEEP CUMULUS CONVECTION.
  ! IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
  ! CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
  ! (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
  !
  !  *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
  !  *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
  !  *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
  !
  !  REFERENCE
  !  ---------
  !          (TIEDTKE,1989)
  !
  subroutine cuasc(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,puen,&
                   pven,ktrac,pxtenh,pxten,pxtu,pmfuxt,pten,pqen,   &
                   pqsen,pgeo,pgeoh,paphp1,pqte,pverv,klwmin,ldcum, &
                   ldland,ktype,klab,ptu,pqu,plu,puu,pvu,pmfu,pmfub,&
                   pentr,pmfus,pmfuq,pmful,plude,pqude,pdmfup,khmin,&
                   phhatt,phcbase,pqsenh,pcpen,pcpcu,kcbot,kctop,   &
                   kctop0)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac , &
               klevm1
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , kctop0 , khmin ,    &
                                  klwmin , ktype
    integer(ik4) , dimension(kbdim,klev) :: klab
    logical , dimension(kbdim) :: ldcum , ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
           pgeoh , phhatt , plu , plude , pmfu , pmful , pmfuq ,      &
           pmfus , pqen , pqenh , pqsen , pqsenh , pqte , pqu ,       &
           pqude , pten , ptenh , ptu , puen , puu , pven , pverv ,   &
           pvu
    real(rkx) , dimension(kbdim) :: pentr , phcbase , pmfub
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
           pxtu
    intent (in) ldland , pcpcu , phcbase , phhatt , pqsenh , pxtenh
    intent (inout) kcbot , kctop , kctop0 , klab , ktype , ldcum ,    &
                   plu , plude , pmfu , pmfub , pmful , pmfuq ,       &
                   pmfus , pmfuxt , pqu , pqude , ptu , puu , pvu ,   &
                   pxtu
    integer(ik4) :: icall , ik , ikb , ikt , jk , jl , jt
    logical , dimension(kbdim) :: loflag
    real(rkx) :: zalvs , zbuo , zbuoyz , zcons2 , zdmfdu ,      &
               zdmfeu , zdnoprc , zdprho , zdrodz , zdt , zdz , zfac ,&
               zga , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,     &
               zmfusk , zmfuxtk , zmse , znevn , zodmax , zprcon ,    &
               zqcod , zqeen , zqude , zscde , zscod , zseen ,        &
               ztglace , zxteen , zxtude , zz , zzdmf
    real(rkx) , dimension(kbdim) :: zbuoy , zdmfde , zdmfen , zmfuu ,   &
                                  zmfuv , zpbase , zph , zqold
    real(rkx) , dimension(kbdim,klev) :: zodetr , zoentr
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuasc'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !----------------------
    ! 1. SPECIFY PARAMETERS
    ! ---------------------
    !
    zcons2 = d_one/(egrav*dtcum)
    ztglace = tzero - 13.0_rkx
    zqold(1:kproma) = d_zero
    !
    ! AMT NOTE!!! in the original scheme, this level which restricts rainfall
    ! below a certain pressure (from the surface) is hard wired according
    ! to the vertical resolution of the model
    ! if ( klev /= 11 ) then
    !   zdlev = 3.0e4_rkx
    ! else if ( nn == 21 ) then
    !   zdlev = 1.5e4_rkx
    ! else if ( nn == 31 ) then
    !   zdlev = 2.0e4_rkx
    ! else
    !   zdlev = 3.0e4_rkx
    ! end if
    !
    !----------------------
    ! 2. SET DEFAULT VALUES
    ! ---------------------
    !
    do jl = 1 , kproma
      zmfuu(jl) = d_zero
      zmfuv(jl) = d_zero
      if ( .not.ldcum(jl) ) ktype(jl) = 0
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        plu(jl,jk) = d_zero
        pmfu(jl,jk) = d_zero
        pmfus(jl,jk) = d_zero
        pmfuq(jl,jk) = d_zero
        pmful(jl,jk) = d_zero
        plude(jl,jk) = d_zero
        pqude(jl,jk) = d_zero
        pdmfup(jl,jk) = d_zero
        if ( .not.ldcum(jl) .or. ktype(jl) == 3 ) klab(jl,jk) = 0
        if ( .not.ldcum(jl) .and. paphp1(jl,jk) < 4.e4_rkx ) kctop0(jl) = jk
        if ( jk < kcbot(jl) ) klab(jl,jk) = 0
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          pmfuxt(jl,jk,jt) = d_zero
        end do
      end do
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        zoentr(jl,jk) = d_zero
        zodetr(jl,jk) = d_zero
      end do
    end do
    !
    !---------------------------------------
    ! 3.0 INITIALIZE VALUES AT LIFTING LEVEL
    ! --------------------------------------
    !
    do jl = 1 , kproma
      kctop(jl) = klevm1
      if ( .not.ldcum(jl) ) then
        kcbot(jl) = klevm1
        pmfub(jl) = d_zero
        pqu(jl,klev) = d_zero
      end if
      pmfu(jl,klev) = pmfub(jl)
      pmfus(jl,klev) = pmfub(jl)                                      &
                       *(pcpcu(jl,klev)*ptu(jl,klev)+pgeoh(jl,klev))
      pmfuq(jl,klev) = pmfub(jl)*pqu(jl,klev)
      if ( lmfdudv ) then
        zmfuu(jl) = pmfub(jl)*puu(jl,klev)
        zmfuv(jl) = pmfub(jl)*pvu(jl,klev)
      end if
    end do

    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = d_zero
        pmfuxt(jl,klev,jt) = pmfub(jl)*pxtu(jl,klev,jt)
      end do
    end do

    do jl = 1 , kproma
      ldcum(jl) = .false.
    end do
    !
    !---------------------------------------------
    ! 3.5 FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
    ! --------------------------------------------
    !
    do jl = 1 , kproma
      if ( ktype(jl) == 1 ) then
        ikb = kcbot(jl)
        zbuoy(jl) = egrav*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb)     &
                    + egrav*ep1*(pqu(jl,ikb)-pqenh(jl,ikb))
        if ( zbuoy(jl) > d_zero ) then
          zdz = (pgeo(jl,ikb-1)-pgeo(jl,ikb))*regrav
          zdrodz = -log(pten(jl,ikb-1)/pten(jl,ikb))                  &
                   /zdz - egrav/(rgas*ptenh(jl,ikb)                   &
                   *(d_one+ep1*pqenh(jl,ikb)))
          ! nb zoentr is here a fractional value
          zoentr(jl,ikb-1) = zbuoy(jl)*d_half/(d_one+zbuoy(jl)*zdz)   &
                             + zdrodz
          zoentr(jl,ikb-1) = min(zoentr(jl,ikb-1),entrmax)
          zoentr(jl,ikb-1) = max(zoentr(jl,ikb-1),d_zero)
        end if
      end if
    end do
    !
    !-------------------------------------------------------
    ! 4. DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
    !    BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
    !    BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
    !    THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
    ! ------------------------------------------------------
    !
    do jk = klevm1 , 2 , -1
      !
      ! SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
      ! IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
      ! ----------------------------------------------------
      !
      ik = jk
      if ( lmfmid .and. ik < klevm1 .and. ik > nmctop )               &
        call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,       &
           puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
           ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
           pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)

      do jl = 1 , kproma
        if ( klab(jl,jk+1) == 0 ) klab(jl,jk) = 0
        loflag(jl) = klab(jl,jk+1) > 0
        zph(jl) = paphp1(jl,jk)
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          if ( pmfub(jl) > zmfmax ) then
            zfac = zmfmax/pmfub(jl)
            pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
            pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
            pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
            zmfuu(jl) = zmfuu(jl)*zfac
            zmfuv(jl) = zmfuv(jl)*zfac
          end if
        end if
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
            zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
            if ( pmfub(jl) > zmfmax ) then
              zfac = zmfmax/pmfub(jl)
              pmfuxt(jl,jk+1,jt) = pmfuxt(jl,jk+1,jt)*zfac
            end if
          end if
        end do
      end do
      !
      ! RESET PMFUB IF NECESSARY
      !
      do jl = 1 , kproma
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          pmfub(jl) = min(pmfub(jl),zmfmax)
        end if
      end do
      !
      !
      ! SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
      ! RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
      ! --------------------------------------------------
      !
      ik = jk
      call cuentr(kproma,kbdim,klev,klevp1,ik,ptenh,pqenh,pqte,paphp1,&
                  klwmin,ldcum,ktype,kcbot,kctop0,zpbase,pmfu,pentr,  &
                  zodetr,khmin,pgeoh,zdmfen,zdmfde)
      !
      ! DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
      ! THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
      ! IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
      ! ARE DETRAINED
      ! IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
      ! MOISTURE THAT ARE NEUTRAL COMPARED TO THE
      ! ENVIRONMENTAL AIR ARE DETRAINED
      ! ---------------------------------------------------
      !
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( jk < kcbot(jl) ) then
            zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
            zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))      &
                     *zcons2)
            zdmfen(jl) = max(zdmfen(jl)-max(zmftest-zmfmax,d_zero),   &
                         d_zero)
          end if
          zdmfde(jl) = min(zdmfde(jl),0.750_rkx*pmfu(jl,jk+1))
          pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          if ( ktype(jl) == 1 .and. jk < kcbot(jl) ) then
            zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
            zoentr(jl,jk) = zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
            zmftest = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
            zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))      &
                     *zcons2)
            zoentr(jl,jk) = max(zoentr(jl,jk)-max(zmftest-zmfmax,d_zero),d_zero)
          else
            zoentr(jl,jk) = d_zero
          end if
          if ( ktype(jl) == 1 .and. jk < kcbot(jl) .and. jk <= khmin(jl) ) then
            ! limit organized detrainment to not allowing for too
            ! deep clouds
            zalvs = merge(wlhv,wlhs,ptu(jl,jk+1) > tzero)
            zmse = pcpcu(jl,jk+1)*ptu(jl,jk+1) + zalvs*pqu(jl,jk+1)   &
                   + pgeoh(jl,jk+1)
            ikt = kctop0(jl)
            znevn = (pgeoh(jl,ikt)-pgeoh(jl,jk+1))                    &
                    *(zmse-phhatt(jl,jk+1))*regrav
            if ( znevn <= d_zero ) znevn = d_one
            zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
            zodmax = ((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
            zodmax = max(zodmax,d_zero)
            zodetr(jl,jk) = min(zodetr(jl,jk),zodmax)
          end if
          zodetr(jl,jk) = min(zodetr(jl,jk),0.750_rkx*pmfu(jl,jk))
          pmfu(jl,jk) = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
          zqeen = pqenh(jl,jk+1)*zdmfen(jl)
          zqeen = zqeen + pqenh(jl,jk+1)*zoentr(jl,jk)
          zseen = (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
          zseen = zseen +                                             &
                  (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zoentr(jl,jk)
          zscde = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
          ! find moist static energy that give nonbuoyant air
          zalvs = merge(wlhv,wlhs,ptenh(jl,jk+1) > tzero)
          zga = zalvs*pqsenh(jl,jk+1)/(rwat*(ptenh(jl,jk+1)**2))
          zdt = (plu(jl,jk+1)-ep1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))&
                /(d_one/ptenh(jl,jk+1)+ep1*zga)
          zscod = pcpcu(jl,jk+1)*ptenh(jl,jk+1) + pgeoh(jl,jk+1)      &
                  + pcpcu(jl,jk+1)*zdt
          zscod = max(zscod,d_zero)
          zscde = zscde + zodetr(jl,jk)*zscod
          zqude = pqu(jl,jk+1)*zdmfde(jl)
          zqcod = pqsenh(jl,jk+1) + zga*zdt
          zqcod = max(zqcod,d_zero)
          zqude = zqude + zodetr(jl,jk)*zqcod
          pqude(jl,jk) = zqude
          plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
          plude(jl,jk) = plude(jl,jk) + plu(jl,jk+1)*zodetr(jl,jk)
          zmfusk = pmfus(jl,jk+1) + zseen - zscde
          zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
          zmfulk = pmful(jl,jk+1) - plude(jl,jk)
          plu(jl,jk) = zmfulk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          pqu(jl,jk) = zmfuqk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          ptu(jl,jk) = (zmfusk*(d_one/max(cmfcmin,pmfu(jl,jk)))       &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          ptu(jl,jk) = max(100.0_rkx,ptu(jl,jk))
          ptu(jl,jk) = min(400.0_rkx,ptu(jl,jk))
          zqold(jl) = pqu(jl,jk)
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( loflag(jl) ) then
            zxteen = pxtenh(jl,jk+1,jt)*(zdmfen(jl)+zoentr(jl,jk))
            zxtude = pxtu(jl,jk+1,jt)*(zdmfde(jl)+zodetr(jl,jk))
            zmfuxtk = pmfuxt(jl,jk+1,jt) + zxteen - zxtude
            pxtu(jl,jk,jt) = zmfuxtk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          end if
        end do
      end do
      !
      ! DO CORRECTIONS FOR MOIST ASCENT
      ! BY ADJUSTING T,Q AND L IN *CUADJTQ*
      ! -----------------------------------
      !
      ik = jk
      icall = 1
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)

      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( pqu(jl,jk) < zqold(jl) ) then
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbuo = ptu(jl,jk)*(d_one+ep1*pqu(jl,jk)-plu(jl,jk))    &
                   - ptenh(jl,jk)*(d_one+ep1*pqenh(jl,jk))
            if ( klab(jl,jk+1) == 1 ) zbuo = zbuo + d_half
            if ( zbuo > d_zero .and. pmfu(jl,jk) >= 0.010_rkx*pmfub(jl) .and.&
                 jk >= kctop0(jl) ) then
              kctop(jl) = jk
              ldcum(jl) = .true.
              zdnoprc = merge(zdlev,1.5e4_rkx,ldland(jl))
              zprcon = merge(d_zero,cprcon,zpbase(jl)-paphp1(jl,jk) < zdnoprc)
              zlnew = plu(jl,jk)/(d_one+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
              pdmfup(jl,jk) = max(d_zero,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
              plu(jl,jk) = zlnew
            else
              klab(jl,jk) = 0
              pmfu(jl,jk) = d_zero
            end if
          end if
        end if
      end do
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
          pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
          pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        end if
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( loflag(jl) ) pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
        end do
      end do

      if ( lmfdudv ) then
        do jl = 1 , kproma
          zdmfen(jl) = zdmfen(jl) + zoentr(jl,jk)
          zdmfde(jl) = zdmfde(jl) + zodetr(jl,jk)
        end do
        do jl = 1 , kproma
          if ( loflag(jl) ) then
            if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) then
              zz = merge(3.0_rkx,2.0_rkx,zdmfen(jl) == d_zero)
            else
              zz = merge(d_one,d_zero,zdmfen(jl) == d_zero)
            end if
            zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
            zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
            zdmfdu = min(zdmfdu,0.750_rkx*pmfu(jl,jk+1))
            zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk)                &
                        - zdmfdu*puu(jl,jk+1)
            zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk)                &
                        - zdmfdu*pvu(jl,jk+1)
            if ( pmfu(jl,jk) > d_zero ) then
              puu(jl,jk) = zmfuu(jl)*(d_one/pmfu(jl,jk))
              pvu(jl,jk) = zmfuv(jl)*(d_one/pmfu(jl,jk))
            end if
          end if
        end do
      end if
      !
      ! COMPUTE ORGANIZED ENTRAINMENT
      ! FOR USE AT NEXT LEVEL
      ! ------------------------------
      !
      do jl = 1 , kproma
        if ( loflag(jl) .and. ktype(jl) == 1 ) then
          zbuoyz = egrav*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk)      &
                   + egrav*ep1*(pqu(jl,jk)-pqenh(jl,jk))             &
                   - egrav*plu(jl,jk)
          zbuoyz = max(zbuoyz,0.00_rkx)
          zdz = (pgeo(jl,jk-1)-pgeo(jl,jk))*regrav
          zdrodz = -log(pten(jl,jk-1)/pten(jl,jk))                   &
                   /zdz - egrav/(rgas*ptenh(jl,jk)                   &
                   *(d_one+ep1*pqenh(jl,jk)))
          zbuoy(jl) = zbuoy(jl) + zbuoyz*zdz
          zoentr(jl,jk-1) = zbuoyz*d_half/(d_one+zbuoy(jl)) + zdrodz
          zoentr(jl,jk-1) = min(zoentr(jl,jk-1),entrmax)
          zoentr(jl,jk-1) = max(zoentr(jl,jk-1),d_zero)

        end if
      end do

    end do
    !
    !--------------------------------------------------------
    ! 5. DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
    ! -------------------------------------------------------
    !
    ! (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
    ! AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
    ! FROM PREVIOUS CALCULATIONS ABOVE)
    !
    do jl = 1 , kproma
      if ( kctop(jl) == klevm1 ) ldcum(jl) = .false.
      kcbot(jl) = max(kcbot(jl),kctop(jl))
    end do
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        jk = kctop(jl) - 1
        zzdmf = cmfctop
        zdmfde(jl) = (d_one-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk) = zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk) = zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk) = pmfu(jl,jk+1) - zdmfde(jl)
        pdmfup(jl,jk) = d_zero
        pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
        if ( jk >= 2 ) then
          plude(jl,jk-1) = pmful(jl,jk)
          pqude(jl,jk-1) = pmfuq(jl,jk)
        else
          plude(jl,jk) = plude(jl,jk) + pmful(jl,jk)
          pqude(jl,jk) = pqude(jl,jk) + pmfuq(jl,jk)
        end if
      end if
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          jk = kctop(jl) - 1
          pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
        end if
      end do
    end do

    if ( lmfdudv ) then
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          jk = kctop(jl) - 1
          puu(jl,jk) = puu(jl,jk+1)
          pvu(jl,jk) = pvu(jl,jk+1)
        end if
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuasc
  !
  ! THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
  ! FOR CUMULUS PARAMETERIZATION
  ! TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
  ! (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
  ! FLUXES AS WELL AS PRECIPITATION RATES)
  !
  ! LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
  ! AND THEN CALCULATE MOIST ASCENT FOR
  ! ENTRAINING/DETRAINING PLUME.
  ! ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
  ! SHALLOW AND DEEP CUMULUS CONVECTION.
  ! IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
  ! CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
  ! (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
  !
  !  *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
  !  *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
  !  *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
  !
  !  REFERENCE
  !  ---------
  !          (TIEDTKE,1989)
  !
  subroutine cuasct(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,    &
                    puen,pven,ktrac,pxtenh,pxten,pxtu,pmfuxt,pten,  &
                    pqen,pqsen,pgeo,pgeoh,paphp1,pqte,pverv,klwmin, &
                    ldcum,ldland,ktype,klab,ptu,pqu,plu,puu,pvu,    &
                    pmfu,pmfub,pentr,pmfus,pmfuq,pmful,plude,pqude, &
                    pdmfup,pcpen,pcpcu,kcbot,kctop,kctop0)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , &
                                 klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , kctop0 , klwmin ,   &
                                  ktype
    integer(ik4) , dimension(kbdim,klev) :: klab
    logical , dimension(kbdim) :: ldcum , ldland
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
           pgeoh , plu , plude , pmfu , pmful , pmfuq , pmfus , pqen ,&
           pqenh , pqsen , pqte , pqu , pqude , pten , ptenh , ptu ,  &
           puen , puu , pven , pverv , pvu
    real(rkx) , dimension(kbdim) :: pentr , pmfub
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
           pxtu
    intent (in) ldland , pcpcu , pxtenh
    intent (out) pqude
    intent (inout) kcbot , kctop , klab , ktype , ldcum , plu ,       &
                   plude , pmfu , pmfub , pmful , pmfuq , pmfus ,     &
                   pmfuxt , pqu , ptu , puu , pvu , pxtu
    integer(ik4) :: icall , ik , jk , jl , jt
    logical , dimension(kbdim) :: loflag
    real(rkx) :: zbuo , zcons2 , zdmfdu , zdmfeu , zdnoprc ,    &
               zfac , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,    &
               zmfusk , zmfuxtk , zprcon , zqeen , zqude , zscde ,    &
               zseen , ztglace , zxteen , zxtude , zz , zzdmf
    real(rkx) , dimension(kbdim) :: zdmfde , zdmfen , zmfuu , zmfuv ,   &
                                  zpbase , zph , zqold
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuasct'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !----------------------
    ! 1. SPECIFY PARAMETERS
    ! ---------------------
    !
    zcons2 = d_one/(egrav*dtcum)
    ztglace = tzero - 13.0_rkx
    !
    ! AMT NOTE!!! in the original scheme, this level which restricts rainfall
    ! below a certain pressure (from the surface) is hard wired according
    ! to the vertical resolution of the model
    ! if ( klev /= 11 ) then
    !   zdlev = 3.0e4_rkx
    ! else if ( nn == 21 ) then
    !   zdlev = 1.5e4_rkx
    ! else if ( nn == 31 ) then
    !   zdlev = 2.0e4_rkx
    ! else
    !   zdlev = 3.0e4_rkx
    ! end if
    !
    !----------------------
    ! 2. SET DEFAULT VALUES
    ! ---------------------
    !
    do jl = 1 , kproma
      zmfuu(jl) = d_zero
      zmfuv(jl) = d_zero
      if ( .not.ldcum(jl) ) ktype(jl) = 0
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        plu(jl,jk) = d_zero
        pmfu(jl,jk) = d_zero
        pmfus(jl,jk) = d_zero
        pmfuq(jl,jk) = d_zero
        pmful(jl,jk) = d_zero
        plude(jl,jk) = d_zero
        pqude(jl,jk) = d_zero
        pdmfup(jl,jk) = d_zero
        if ( .not.ldcum(jl) .or. ktype(jl) == 3 ) klab(jl,jk) = 0
        if ( .not.ldcum(jl) .and. paphp1(jl,jk) < 4.e4_rkx ) kctop0(jl) = jk
        if ( jk < kcbot(jl) ) klab(jl,jk) = 0
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          pmfuxt(jl,jk,jt) = d_zero
        end do
      end do
    end do
    !
    !---------------------------------------
    ! 3.0 INITIALIZE VALUES AT LIFTING LEVEL
    ! --------------------------------------
    !
    do jl = 1 , kproma
      kctop(jl) = klevm1
      if ( .not.ldcum(jl) ) then
        kcbot(jl) = klevm1
        pmfub(jl) = d_zero
        pqu(jl,klev) = d_zero
      end if
      pmfu(jl,klev) = pmfub(jl)
      pmfus(jl,klev) = pmfub(jl)                                      &
                       *(pcpcu(jl,klev)*ptu(jl,klev)+pgeoh(jl,klev))
      pmfuq(jl,klev) = pmfub(jl)*pqu(jl,klev)
      if ( lmfdudv ) then
        zmfuu(jl) = pmfub(jl)*puu(jl,klev)
        zmfuv(jl) = pmfub(jl)*pvu(jl,klev)
      end if
    end do

    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = d_zero
        pmfuxt(jl,klev,jt) = pmfub(jl)*pxtu(jl,klev,jt)
      end do
    end do

    do jl = 1 , kproma
      ldcum(jl) = .false.
    end do
    !
    !-------------------------------------------------------
    ! 4. DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
    !    BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
    !    BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
    !    THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
    ! ------------------------------------------------------
    !
    do jk = klevm1 , 2 , -1
      !
      ! SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
      ! IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
      ! ----------------------------------------------------
      !
      ik = jk
      if ( lmfmid .and. ik < klevm1 .and. ik > nmctop )               &
        call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,       &
           puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
           ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
           pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)

      do jl = 1 , kproma
        if ( klab(jl,jk+1) == 0 ) klab(jl,jk) = 0
        loflag(jl) = klab(jl,jk+1) > 0
        zph(jl) = paphp1(jl,jk)
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          if ( pmfub(jl) > zmfmax ) then
            zfac = zmfmax/pmfub(jl)
            pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
            pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
            pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
            zmfuu(jl) = zmfuu(jl)*zfac
            zmfuv(jl) = zmfuv(jl)*zfac
          end if
        end if
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
            zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
            if ( pmfub(jl) > zmfmax ) then
              zfac = zmfmax/pmfub(jl)
              pmfuxt(jl,jk+1,jt) = pmfuxt(jl,jk+1,jt)*zfac
            end if
          end if
        end do
      end do
      !
      ! RESET PMFUB IF NECESSARY
      !
      do jl = 1 , kproma
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          pmfub(jl) = min(pmfub(jl),zmfmax)
        end if
      end do
      !
      ! SPECIFY ENTRAINMENT RATES IN *CUENTRT*
      ! --------------------------------------
      !
      ik = jk
      call cuentrt(kproma,kbdim,klev,klevp1,ik,ptenh,pqenh,pqte,      &
                   paphp1,klwmin,ldcum,ktype,kcbot,kctop0,zpbase,pmfu,&
                   pentr,zdmfen,zdmfde)
      !
      ! DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
      ! ---------------------------------------------------
      !
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( jk < kcbot(jl) ) then
            zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
            zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
            zdmfen(jl) = max(zdmfen(jl)-max(zmftest-zmfmax,d_zero),d_zero)
          end if
          zdmfde(jl) = min(zdmfde(jl),0.750_rkx*pmfu(jl,jk+1))
          pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zqeen = pqenh(jl,jk+1)*zdmfen(jl)
          zseen = (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
          zscde = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
          zqude = pqu(jl,jk+1)*zdmfde(jl)
          pqude(jl,jk) = zqude
          plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
          zmfusk = pmfus(jl,jk+1) + zseen - zscde
          zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
          zmfulk = pmful(jl,jk+1) - plude(jl,jk)
          plu(jl,jk) = zmfulk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          pqu(jl,jk) = zmfuqk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          ptu(jl,jk) = (zmfusk*(d_one/max(cmfcmin,pmfu(jl,jk)))       &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          ptu(jl,jk) = max(100.0_rkx,ptu(jl,jk))
          ptu(jl,jk) = min(400.0_rkx,ptu(jl,jk))
          zqold(jl) = pqu(jl,jk)
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( loflag(jl) ) then
            zxteen = pxtenh(jl,jk+1,jt)*zdmfen(jl)
            zxtude = pxtu(jl,jk+1,jt)*zdmfde(jl)
            zmfuxtk = pmfuxt(jl,jk+1,jt) + zxteen - zxtude
            pxtu(jl,jk,jt) = zmfuxtk*(d_one/max(cmfcmin,pmfu(jl,jk)))
          end if
        end do
      end do
      !
      ! DO CORRECTIONS FOR MOIST ASCENT
      ! BY ADJUSTING T,Q AND L IN *CUADJTQ*
      ! -----------------------------------
      !
      ik = jk
      icall = 1
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)

      do jl = 1 , kproma
        if ( loflag(jl) .and. pqu(jl,jk) < zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+ep1*pqu(jl,jk)-plu(jl,jk))      &
                 - ptenh(jl,jk)*(d_one+ep1*pqenh(jl,jk))
          if ( klab(jl,jk+1) == 1 ) zbuo = zbuo + d_half
          if ( zbuo > d_zero .and. pmfu(jl,jk) >= 0.10_rkx*pmfub(jl) ) then
            kctop(jl) = jk
            ldcum(jl) = .true.
            zdnoprc = merge(zdlev,1.5e4_rkx,ldland(jl))
            zprcon = merge(d_zero,cprcon,zpbase(jl)-paphp1(jl,jk) < zdnoprc)
            zlnew = plu(jl,jk)                                        &
                    /(d_one+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
            pdmfup(jl,jk) = max(d_zero,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
            plu(jl,jk) = zlnew
          else
            klab(jl,jk) = 0
            pmfu(jl,jk) = d_zero
          end if
        end if
      end do
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
          pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
          pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        end if
      end do
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( loflag(jl) ) pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
        end do
      end do

      if ( lmfdudv ) then
        do jl = 1 , kproma
          if ( loflag(jl) ) then
            if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) then
              zz = merge(3.0_rkx,2.0_rkx,zdmfen(jl) == d_zero)
            else
              zz = merge(d_one,d_zero,zdmfen(jl) == d_zero)
            end if
            zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
            zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
            zdmfdu = min(zdmfdu,0.750_rkx*pmfu(jl,jk+1))
            zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk) - zdmfdu*puu(jl,jk+1)
            zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk) - zdmfdu*pvu(jl,jk+1)
            if ( pmfu(jl,jk) > d_zero ) then
              puu(jl,jk) = zmfuu(jl)*(d_one/pmfu(jl,jk))
              pvu(jl,jk) = zmfuv(jl)*(d_one/pmfu(jl,jk))
            end if
          end if
        end do
      end if

    end do
    !
    !
    !--------------------------------------------------------
    ! 5. DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
    ! -------------------------------------------------------
    !  (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
    !  AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
    !  FROM PREVIOUS CALCULATIONS ABOVE)
    !
    do jl = 1 , kproma
      if ( kctop(jl) == klevm1 ) ldcum(jl) = .false.
      kcbot(jl) = max(kcbot(jl),kctop(jl))
    end do
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        jk = kctop(jl) - 1
        zzdmf = cmfctop
        zdmfde(jl) = (d_one-zzdmf)*pmfu(jl,jk+1)
        plude(jl,jk) = zdmfde(jl)*plu(jl,jk+1)
        pqude(jl,jk) = zdmfde(jl)*pqu(jl,jk+1)
        pmfu(jl,jk) = pmfu(jl,jk+1) - zdmfde(jl)
        pdmfup(jl,jk) = d_zero
        pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
        pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
        plude(jl,jk-1) = pmful(jl,jk)
        pqude(jl,jk-1) = pmfuq(jl,jk)
      end if
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          jk = kctop(jl) - 1
          pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)*pmfu(jl,jk)
        end if
      end do
    end do

    if ( lmfdudv ) then
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          jk = kctop(jl) - 1
          puu(jl,jk) = puu(jl,jk+1)
          pvu(jl,jk) = pvu(jl,jk+1)
        end if
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuasct
  !
  ! THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
  ! FOR CUMULUS PARAMETERIZATION
  ! TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
  !
  ! INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
  ! IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
  !   KLAB=1 FOR SUBCLOUD LEVELS
  !   KLAB=2 FOR CONDENSATION LEVEL
  !
  ! LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
  ! (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
  !
  !  *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
  !
  subroutine cubase(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,    &
                    pgeoh,paph,ptu,pqu,plu,puen,pven,puu,pvu,pcpcu, &
                    ldcum,kcbot,klab)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma
    integer(ik4) , dimension(kbdim) :: kcbot
    integer(ik4) , dimension(kbdim,klev) :: klab
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klevp1) :: paph
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pgeoh , plu , pqenh ,  &
           pqu , ptenh , ptu , puen , puu , pven , pvu
    intent (in) paph , pcpcu , pgeoh , pqenh , ptenh , puen , pven
    intent (inout) kcbot , klab , ldcum , plu , pqu , ptu , puu , pvu
    integer(ik4) :: icall , ik , ikb , is , jk , jl
    logical , dimension(kbdim) :: loflag
    real(rkx) :: zbuo , zz
    real(rkx) , dimension(kbdim) :: zph , zqold
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cubase'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !--------------------------------------
    ! 1. INITIALIZE VALUES AT LIFTING LEVEL
    ! -------------------------------------
    !
    do jl = 1 , kproma
      klab(jl,klev) = 1
      kcbot(jl) = klevm1
      ldcum(jl) = .false.
      puu(jl,klev) = puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
      pvu(jl,klev) = pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
    end do
    !
    !-----------------------------------------------
    ! 2.0 DO ASCENT IN SUBCLOUD LAYER,
    !     CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
    !     ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
    !     CHECK FOR BUOYANCY AND SET FLAGS
    ! ----------------------------------------------
    !
    do jk = klevm1 , 2 , -1
      is = 0
      do jl = 1 , kproma
        is = is + merge(1,0,klab(jl,jk+1) == 1)
        loflag(jl) = klab(jl,jk+1) == 1
        zph(jl) = paph(jl,jk)
      end do
      if ( is == 0 ) cycle
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          pqu(jl,jk) = pqu(jl,jk+1)
          ptu(jl,jk) = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)  &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+ep1*pqu(jl,jk)) - ptenh(jl,jk)&
                 *(d_one+ep1*pqenh(jl,jk)) + d_half
          if ( zbuo > d_zero ) klab(jl,jk) = 1
          zqold(jl) = pqu(jl,jk)
        end if
      end do

      ik = jk
      icall = 1
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)

      do jl = 1 , kproma
        if ( loflag(jl) .and. pqu(jl,jk) < zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+ep1*pqu(jl,jk)-plu(jl,jk))    &
                 - ptenh(jl,jk)*(d_one+ep1*pqenh(jl,jk)) + d_half
          if ( zbuo > d_zero ) then
            kcbot(jl) = jk
            ldcum(jl) = .true.
          end if
        end if
      end do
      !
      ! CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
      ! THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
      !
      if ( lmfdudv ) then
        do jl = 1 , kproma
          if ( jk >= kcbot(jl) ) then
            puu(jl,klev) = puu(jl,klev) + puen(jl,jk)               &
                           *(paph(jl,jk+1)-paph(jl,jk))
            pvu(jl,klev) = pvu(jl,klev) + pven(jl,jk)               &
                           *(paph(jl,jk+1)-paph(jl,jk))
          end if
        end do
      end if
    end do

    if ( lmfdudv ) then
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          ikb = kcbot(jl)
          zz = d_one/(paph(jl,klevp1)-paph(jl,ikb))
          puu(jl,klev) = puu(jl,klev)*zz
          pvu(jl,klev) = pvu(jl,klev)*zz
        else
          puu(jl,klev) = puen(jl,klevm1)
          pvu(jl,klev) = pven(jl,klevm1)
        end if
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cubase
  !
  ! THIS ROUTINE CALCULATES CLOUD BASE VALUES
  ! FOR MIDLEVEL CONVECTION
  ! INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
  ! IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
  !
  subroutine cubasmc(kproma,kbdim,klev,kk,klab,pten,pqen,pqsen,puen,&
                     pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh, &
                     ldcum,ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,&
                     puu,pvu,pmfus,pmfuq,pmful,pdmfup,pmfuu,pcpen,  &
                     pmfuv)
    implicit none
    integer(ik4) , intent(in) :: kbdim , kk , klev , kproma , ktrac
    integer(ik4) , dimension(kbdim) :: kcbot , ktype
    integer(ik4) , dimension(kbdim,klev) :: klab
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klev) :: pcpen , pdmfup , pgeo , pgeoh ,&
           plu , pmfu , pmful , pmfuq , pmfus , pqen , pqsen , pqu ,  &
           pten , ptu , puen , puu , pven , pverv , pvu
    real(rkx) , dimension(kbdim) :: pentr , pmfub , pmfuu , pmfuv
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtu
    intent (in) ldcum , pcpen , pgeo , pgeoh , pqen , pqsen , pten , &
                puen , pven , pverv , pxten
    intent (out) kcbot , ktype , pdmfup , pentr , plu , pmfu , pmful ,&
                 pmfuq , pmfus , pmfuu , pmfuv , pmfuxt
    intent (inout) klab , pmfub , pqu , ptu , puu , pvu , pxtu
    integer(ik4) :: jl , jt
    logical , dimension(kbdim) :: llo3
    real(rkx) :: zzzmb
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cubasmc'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-----------------------------------------------
    ! 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! ----------------------------------------------
    !
    do jl = 1 , kproma
      llo3(jl) = .false.
      if ( .not.ldcum(jl) .and. klab(jl,kk+1) == 0 .and. &
            pqen(jl,kk) > 0.900_rkx*pqsen(jl,kk) ) then
        llo3(jl) = .true.
        ptu(jl,kk+1) = (pcpen(jl,kk)*pten(jl,kk)+pgeo(jl,kk)- &
                        pgeoh(jl,kk+1))/pcpen(jl,kk)
        pqu(jl,kk+1) = pqen(jl,kk)
        plu(jl,kk+1) = d_zero
        zzzmb = max(cmfcmin,-pverv(jl,kk)*regrav)
        zzzmb = min(zzzmb,cmfcmax)
        pmfub(jl) = zzzmb
        pmfu(jl,kk+1) = pmfub(jl)
        pmfus(jl,kk+1) = pmfub(jl)*(pcpen(jl,kk)*ptu(jl,kk+1)+pgeoh(jl,kk+1))
        pmfuq(jl,kk+1) = pmfub(jl)*pqu(jl,kk+1)
        pmful(jl,kk+1) = d_zero
        pdmfup(jl,kk+1) = d_zero
        kcbot(jl) = kk
        klab(jl,kk+1) = 1
        ktype(jl) = 3
        pentr(jl) = entrmid
        if ( lmfdudv ) then
          puu(jl,kk+1) = puen(jl,kk)
          pvu(jl,kk+1) = pven(jl,kk)
          pmfuu(jl) = pmfub(jl)*puu(jl,kk+1)
          pmfuv(jl) = pmfub(jl)*pvu(jl,kk+1)
        end if
      end if
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( llo3(jl) ) then
          pxtu(jl,kk+1,jt) = pxten(jl,kk,jt)
          pmfuxt(jl,kk+1,jt) = pmfub(jl)*pxtu(jl,kk+1,jt)
        end if
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cubasmc
  !
  ! THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
  ! TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
  ! (I.E. T,Q,U AND V AND FLUXES)
  ! INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
  ! IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
  ! AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
  !
  ! CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
  !   A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
  !   B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
  !
  ! *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
  ! SATURATED DESCENT
  !
  ! REFERENCE
  ! ---------
  !          (TIEDTKE,1989)
  !
  subroutine cuddraf(kproma,kbdim,klev,klevp1,ptenh,pqenh,puen,pven,&
                     ktrac,pxtenh,pxtd,pmfdxt,pgeoh,paphp1,prfl,ptd,&
                     pqd,pud,pvd,pmfd,pmfds,pmfdq,pdmfdp,pcpcu,     &
                     lddraf)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
    logical , dimension(kbdim) :: lddraf
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
           pmfdq , pmfds , pqd , pqenh , ptd , ptenh , pud , puen ,   &
           pvd , pven
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh
    real(rkx) , dimension(kbdim) :: prfl
    intent (in) lddraf , paphp1 , pcpcu , pgeoh , pqenh , ptenh , &
                puen , pven , pxtenh
    intent (out) pdmfdp
    intent (inout) pmfd , pmfdq , pmfds , pmfdxt , pqd , prfl , ptd , &
                   pud , pvd , pxtd
    integer(ik4) :: icall , ik , is , itopde , jk , jl , jt
    logical :: llo1
    logical , dimension(kbdim) :: llo2
    real(rkx) :: zbuo , zdmfdp , zentr , zmfdqk , zmfdsk , zmfduk ,     &
               zmfdvk , zmfdxtk , zqdde , zqeen , zsdde , zseen ,     &
               zxtdde , zxteen
    real(rkx) , dimension(kbdim) :: zcond , zdmfde , zdmfen , zph
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuddraf'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-----------------------------------------------------
    ! 1. CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
    !     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
    !         LINEAR DECREASE OF MASSFLUX IN PBL
    !     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
    !         AND MOISTENING IS CALCULATED IN *CUADJTQ*
    !     (C) CHECKING FOR NEGATIVE BUOYANCY AND
    !         SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
    ! ----------------------------------------------------
    !
    do jk = 3 , klev
      is = 0
      do jl = 1 , kproma
        zph(jl) = paphp1(jl,jk)
        llo2(jl) = lddraf(jl) .and. pmfd(jl,jk-1) < d_zero
        is = is + merge(1,0,llo2(jl))
      end do
      if ( is == 0 ) cycle
      do jl = 1 , kproma
        if ( llo2(jl) ) then
          zentr = entrdd*pmfd(jl,jk-1)*rgas*ptenh(jl,jk-1)          &
                  /(egrav*paphp1(jl,jk-1))*(paphp1(jl,jk)-paphp1(jl,jk-1))
          zdmfen(jl) = zentr
          zdmfde(jl) = zentr
        end if
      end do
      itopde = klev - 2
      if ( jk > itopde ) then
        do jl = 1 , kproma
          if ( llo2(jl) ) then
            zdmfen(jl) = d_zero
            zdmfde(jl) = pmfd(jl,itopde)*(paphp1(jl,jk)-paphp1(jl,jk-1)) &
                         /(paphp1(jl,klevp1)-paphp1(jl,itopde))
          end if
        end do
      end if

      do jl = 1 , kproma
        if ( llo2(jl) ) then
          pmfd(jl,jk) = pmfd(jl,jk-1) + zdmfen(jl) - zdmfde(jl)
          zseen = (pcpcu(jl,jk-1)*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
          zqeen = pqenh(jl,jk-1)*zdmfen(jl)
          zsdde = (pcpcu(jl,jk-1)*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
          zqdde = pqd(jl,jk-1)*zdmfde(jl)
          zmfdsk = pmfds(jl,jk-1) + zseen - zsdde
          zmfdqk = pmfdq(jl,jk-1) + zqeen - zqdde
          pqd(jl,jk) = zmfdqk*(d_one/min(-cmfcmin,pmfd(jl,jk)))
          ptd(jl,jk) = (zmfdsk*(d_one/min(-cmfcmin,pmfd(jl,jk)))    &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          ptd(jl,jk) = min(400.0_rkx,ptd(jl,jk))
          ptd(jl,jk) = max(100.0_rkx,ptd(jl,jk))
          zcond(jl) = pqd(jl,jk)
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( llo2(jl) ) then
            zxteen = pxtenh(jl,jk-1,jt)*zdmfen(jl)
            zxtdde = pxtd(jl,jk-1,jt)*zdmfde(jl)
            zmfdxtk = pmfdxt(jl,jk-1,jt) + zxteen - zxtdde
            pxtd(jl,jk,jt) = zmfdxtk*(d_one/min(-cmfcmin,pmfd(jl,jk)))
          end if
        end do
      end do

      ik = jk
      icall = 2
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptd,pqd,llo2,icall)

      do jl = 1 , kproma
        if ( llo2(jl) ) then
          zcond(jl) = zcond(jl) - pqd(jl,jk)
          zbuo = ptd(jl,jk)*(d_one+ep1*pqd(jl,jk)) - &
                 ptenh(jl,jk)*(d_one+ep1*pqenh(jl,jk))
          llo1 = zbuo < d_zero .and.                                   &
                 (prfl(jl)-pmfd(jl,jk)*zcond(jl) > d_zero)
          pmfd(jl,jk) = merge(pmfd(jl,jk),d_zero,llo1)
          pmfds(jl,jk) = (pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
          pmfdq(jl,jk) = pqd(jl,jk)*pmfd(jl,jk)
          zdmfdp = -pmfd(jl,jk)*zcond(jl)
          pdmfdp(jl,jk-1) = zdmfdp
          prfl(jl) = prfl(jl) + zdmfdp
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( llo2(jl) ) pmfdxt(jl,jk,jt) = pxtd(jl,jk,jt)*pmfd(jl,jk)
        end do
      end do

      if ( lmfdudv ) then
        do jl = 1 , kproma
          if ( llo2(jl) .and. pmfd(jl,jk) < d_zero ) then
            zmfduk = pmfd(jl,jk-1)*pud(jl,jk-1) + zdmfen(jl)        &
                     *puen(jl,jk-1) - zdmfde(jl)*pud(jl,jk-1)
            zmfdvk = pmfd(jl,jk-1)*pvd(jl,jk-1) + zdmfen(jl)        &
                     *pven(jl,jk-1) - zdmfde(jl)*pvd(jl,jk-1)
            pud(jl,jk) = zmfduk*(d_one/min(-cmfcmin,pmfd(jl,jk)))
            pvd(jl,jk) = zmfdvk*(d_one/min(-cmfcmin,pmfd(jl,jk)))
          end if
        end do
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuddraf
  !
  ! THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
  ! CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
  ! TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
  ! FOR MASSFLUX CUMULUS PARAMETERIZATION
  !
  ! INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
  ! AND UPDRAFT VALUES T,Q,U AND V AND ALSO
  ! CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
  ! IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
  !
  ! CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
  ! MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
  !
  ! *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
  !
  subroutine cudlfs(kproma,kbdim,klev,klevp1,ptenh,pqenh,puen,pven, &
                    ktrac,pxtenh,pxtu,pxtd,pmfdxt,pgeoh,paphp1,ptu, &
                    pqu,puu,pvu,ldcum,kcbot,kctop,pmfub,prfl,ptd,   &
                    pqd,pud,pvd,pmfd,pmfds,pmfdq,pdmfdp,pcpcu,kdtop,&
                    lddraf)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , kdtop
    logical , dimension(kbdim) :: ldcum , lddraf
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
           pmfdq , pmfds , pqd , pqenh , pqu , ptd , ptenh , ptu ,    &
           pud , puen , puu , pvd , pven , pvu
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh , &
           pxtu
    real(rkx) , dimension(kbdim) :: pmfub , prfl
    intent (in) kcbot , kctop , ldcum , paphp1 , pcpcu , pgeoh , &
                pmfub , pqenh , pqu , ptenh , ptu , puen , puu , &
                pven , pvu , pxtenh , pxtu
    intent (out) kdtop , pmfdq , pmfds , pmfdxt , pud , pvd
    intent (inout) lddraf , pdmfdp , pmfd , pqd , prfl , ptd , pxtd
    integer(ik4) :: icall , ik , is , jk , jl , jt , ke
    logical , dimension(kbdim) :: llo2 , llo3
    real(rkx) :: zbuo , zmftop , zqtest , zttest
    real(rkx) , dimension(kbdim) :: zcond , zph
    real(rkx) , dimension(kbdim,klev) :: zqenwb , ztenwb
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cudlf'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-------------------------------------
    ! 1. SET DEFAULT VALUES FOR DOWNDRAFTS
    ! ------------------------------------
    !
    do jl = 1 , kproma
      lddraf(jl) = .false.
      kdtop(jl) = klevp1
    end do

    if ( .not. lmfdd ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if
    !
    !---------------------------------------------------------
    ! 2. DETERMINE LEVEL OF FREE SINKING BY
    !    DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
    !
    !    FOR EVERY POINT AND PROCEED AS FOLLOWS:
    !
    !    (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
    !    (2) DO MIXING WITH CUMULUS CLOUD AIR
    !    (3) CHECK FOR NEGATIVE BUOYANCY
    !
    !    THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
    !    OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
    !    TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
    !    EVAPORATION OF RAIN AND CLOUD WATER)
    ! -------------------------------------------------------
    !
    ke = klev - 3
    do jk = 3 , ke
      !
      ! 2.1 CALCULATE WET-BULB TEMPERATURE AND MOISTURE
      !     FOR ENVIRONMENTAL AIR IN *CUADJTQ*
      ! -----------------------------------------------
      !
      is = 0
      do jl = 1 , kproma
        ztenwb(jl,jk) = ptenh(jl,jk)
        zqenwb(jl,jk) = pqenh(jl,jk)
        zph(jl) = paphp1(jl,jk)
        llo2(jl) = ldcum(jl) .and. prfl(jl) > d_zero .and.  &
                   .not.lddraf(jl) .and.                    &
                   (jk < kcbot(jl) .and. jk > kctop(jl))
        is = is + merge(1,0,llo2(jl))
      end do
      if ( is == 0 ) cycle
      ik = jk
      icall = 2
      call cuadjtq(kproma,kbdim,klev,ik,zph,ztenwb,zqenwb,llo2,icall)
      !
      ! 2.2 DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
      !     AND CHECK FOR NEGATIVE BUOYANCY.
      !     THEN SET VALUES FOR DOWNDRAFT AT LFS.
      ! ----------------------------------------------
      !
      do jl = 1 , kproma
        llo3(jl) = .false.
        if ( llo2(jl) ) then
          zttest = d_half*(ptu(jl,jk)+ztenwb(jl,jk))
          zqtest = d_half*(pqu(jl,jk)+zqenwb(jl,jk))
          zbuo = zttest*(d_one+ep1*zqtest) - ptenh(jl,jk)      &
                 *(d_one+ep1*pqenh(jl,jk))
          zcond(jl) = pqenh(jl,jk) - zqenwb(jl,jk)
          zmftop = -cmfdeps*pmfub(jl)
          if ( zbuo < d_zero .and. prfl(jl) > 10.0_rkx*zmftop*zcond(jl) ) then
            llo3(jl) = .true.
            kdtop(jl) = jk
            lddraf(jl) = .true.
            ptd(jl,jk) = zttest
            pqd(jl,jk) = zqtest
            pmfd(jl,jk) = zmftop
            pmfds(jl,jk) = pmfd(jl,jk)*(pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk) = pmfd(jl,jk)*pqd(jl,jk)
            pdmfdp(jl,jk-1) = -d_half*pmfd(jl,jk)*zcond(jl)
            prfl(jl) = prfl(jl) + pdmfdp(jl,jk-1)
          end if
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( llo3(jl) ) then
            pxtd(jl,jk,jt) = d_half*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
            pmfdxt(jl,jk,jt) = pmfd(jl,jk)*pxtd(jl,jk,jt)
          end if
        end do
      end do

      if ( lmfdudv ) then
        do jl = 1 , kproma
          if ( pmfd(jl,jk) < d_zero ) then
            pud(jl,jk) = d_half*(puu(jl,jk)+puen(jl,jk-1))
            pvd(jl,jk) = d_half*(pvu(jl,jk)+pven(jl,jk-1))
          end if
        end do
      end if

    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cudlfs
  !
  ! UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
  ! DOES GLOBAL DIAGNOSTICS
  !
  subroutine cudtdq(kproma,kbdim,klev,klevp1,ktopm2,ldcum,ktrac,    &
                    paphp1,pten,ptte,pqte,pxtte,pxtec,pmfuxt,pmfdxt,&
                    pmfus,pmfds,pmfuq,pmfdq,pmful,pdmfup,pdmfdp,    &
                    plude,pdpmel,prfl,psfl,pcpen,pqtec,pqude,prsfc, &
                    pssfc)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktopm2 , ktrac
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: prfl , prsfc , psfl , pssfc
    real(rkx) , dimension(kbdim,klev) :: pcpen , pdmfdp , pdmfup ,      &
           pdpmel , plude , pmfdq , pmfds , pmful , pmfuq , pmfus ,   &
           pqte , pqtec , pqude , pten , ptte , pxtec
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtte
    intent (in) ldcum , paphp1 , pcpen , pdmfdp , pdmfup , pdpmel ,   &
                plude , pmfdq , pmfds , pmfdxt , pmful , pmfuq ,      &
                pmfus , pmfuxt , pqude , prfl , psfl , pten
    intent (out) pqtec , prsfc , pssfc , pxtec
    intent (inout) pqte , ptte , pxtte
    integer(ik4) :: jk , jl , jt
    logical :: llo1
    real(rkx) :: zalv , zdqdt , zdtdt , zdxtdt , zrcpm
    real(rkx) , dimension(kbdim) :: zmelt , zsheat
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cudtdq'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-----------------------------------------
    ! 2.0 INCREMENTATION OF T AND Q TENDENCIES
    ! ----------------------------------------
    !
    do jl = 1 , kproma
      zmelt(jl) = d_zero
      zsheat(jl) = d_zero
    end do

    do jk = ktopm2 , klev
      if ( jk < klev ) then
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            llo1 = (pten(jl,jk)-tzero) > d_zero
            zalv = merge(wlhv,wlhs,llo1)
            zrcpm = d_one/pcpen(jl,jk)
            zdtdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))            &
                    *zrcpm*(pmfus(jl,jk+1)-pmfus(jl,jk)+pmfds(jl,jk+1) &
                    -pmfds(jl,jk)-wlhf*pdpmel(jl,jk)                   &
                    -zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl,jk)    &
                    -(pdmfup(jl,jk)+pdmfdp(jl,jk))))
            ptte(jl,jk) = ptte(jl,jk) + zdtdt
            zdqdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))             &
                    *(pmfuq(jl,jk+1)-pmfuq(jl,jk)+pmfdq(jl,jk+1)        &
                    -pmfdq(jl,jk)+pmful(jl,jk+1)-pmful(jl,jk)           &
                    -plude(jl,jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk)))
            pqte(jl,jk) = pqte(jl,jk) + zdqdt
            pxtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                           *plude(jl,jk)
            pqtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                           *pqude(jl,jk)
            zsheat(jl) = zsheat(jl)+zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zmelt(jl) = zmelt(jl) + pdpmel(jl,jk)
          end if
        end do
        do jt = 1 , ktrac
          do jl = 1 , kproma
            if ( ldcum(jl) ) then
              zdxtdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                       *(pmfuxt(jl,jk+1,jt)-pmfuxt(jl,jk,jt)      &
                       +pmfdxt(jl,jk+1,jt)-pmfdxt(jl,jk,jt))
              pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
            end if
          end do
        end do
      else
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            llo1 = (pten(jl,jk)-tzero) > d_zero
            zalv = merge(wlhv,wlhs,llo1)
            zrcpm = d_one/pcpen(jl,jk)
            zdtdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))          &
                    *zrcpm*(pmfus(jl,jk)+pmfds(jl,jk)                 &
                    +wlhf*pdpmel(jl,jk)                               &
                    -zalv*(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)   &
                    +plude(jl,jk)))
            ptte(jl,jk) = ptte(jl,jk) + zdtdt
            zdqdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))          &
                    *(pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)          &
                    +(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
            pqte(jl,jk) = pqte(jl,jk) + zdqdt
            pxtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                           *plude(jl,jk)
            pqtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))    &
                           *pqude(jl,jk)
            zsheat(jl) = zsheat(jl)                                   &
                         + zalv*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zmelt(jl) = zmelt(jl) + pdpmel(jl,jk)
          end if
        end do
        do jt = 1 , ktrac
          do jl = 1 , kproma
            if ( ldcum(jl) ) then
              zdxtdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))     &
                       *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
              pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
            end if
          end do
        end do
      end if
    end do
    !
    !-------------------------
    ! 3. UPDATE SURFACE FIELDS
    ! ------------------------
    !
    do jl = 1 , kproma
      prsfc(jl) = prfl(jl)
      pssfc(jl) = psfl(jl)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cudtdq
  !
  ! UPDATES U AND V TENDENCIES,
  ! DOES GLOBAL DIAGNOSTIC OF DISSIPATION
  !
  subroutine cududv(kproma,kbdim,klev,klevp1,ktopm2,ktype,kcbot,    &
                    paphp1,ldcum,puen,pven,pvom,pvol,puu,pud,pvu,   &
                    pvd,pmfu,pmfd)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktopm2
    integer(ik4) , dimension(kbdim) :: kcbot , ktype
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pmfd , pmfu , pud , puen ,     &
           puu , pvd , pven , pvol , pvom , pvu
    intent (in) kcbot , ktype , ldcum , paphp1 , pmfd , pmfu , pud , &
                puen , puu , pvd , pven , pvu
    intent (inout) pvol , pvom
    integer(ik4) :: ik , ikb , jk , jl
    real(rkx) :: zdudt , zdvdt , zp
    real(rkx) , dimension(kbdim,klev) :: zmfdu , zmfdv , zmfuu , zmfuv
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cududv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !---------------------------------------------------
    ! 1.0 CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
    ! --------------------------------------------------
    !
    if ( ktopm2 == 1 ) then
      do jk = 2 , klev
        ik = jk - 1
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          end if
        end do
      end do
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zmfuu(jl,1) = zmfuu(jl,2)
          zmfuv(jl,1) = zmfuv(jl,2)
          zmfdu(jl,1) = zmfdu(jl,2)
          zmfdv(jl,1) = zmfdv(jl,2)
        end if
      end do
    else
      do jk = ktopm2 , klev
        ik = jk - 1
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-puen(jl,ik))
            zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-pven(jl,ik))
            zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-puen(jl,ik))
            zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-pven(jl,ik))
          end if
        end do
      end do
    end if
    do jk = ktopm2 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk > kcbot(jl) ) then
          ikb = kcbot(jl)
          zp = ((paphp1(jl,klevp1)-paphp1(jl,jk))                    &
                /(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zp = merge(zp**2,zp,ktype(jl) == 3)
          zmfuu(jl,jk) = zmfuu(jl,ikb)*zp
          zmfuv(jl,jk) = zmfuv(jl,ikb)*zp
          zmfdu(jl,jk) = zmfdu(jl,ikb)*zp
          zmfdv(jl,jk) = zmfdv(jl,ikb)*zp
        end if
      end do
    end do

    do jk = ktopm2 , klev
      if ( jk < klev ) then
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zdudt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))           &
                    *(zmfuu(jl,jk+1)-zmfuu(jl,jk)+zmfdu(jl,jk+1)      &
                    -zmfdu(jl,jk))
            zdvdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))           &
                    *(zmfuv(jl,jk+1)-zmfuv(jl,jk)+zmfdv(jl,jk+1)      &
                    -zmfdv(jl,jk))
            pvom(jl,jk) = pvom(jl,jk) + zdudt
            pvol(jl,jk) = pvol(jl,jk) + zdvdt
          end if
        end do
      else
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zdudt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))            &
                    *(zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))            &
                    *(zmfuv(jl,jk)+zmfdv(jl,jk))
            pvom(jl,jk) = pvom(jl,jk) + zdudt
            pvol(jl,jk) = pvol(jl,jk) + zdvdt
          end if
        end do
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cududv
  !
  ! THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
  ! FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
  !
  ! INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
  ! AND UPDRAFT VALUES T,Q ETC
  ! IT RETURNS ENTRAINMENT/DETRAINMENT RATES
  !
  ! METHOD.
  ! --------
  !   S. TIEDTKE (1989)
  !
  subroutine cuentr(kproma,kbdim,klev,klevp1,kk,ptenh,pqenh,pqte,   &
                    paphp1,klwmin,ldcum,ktype,kcbot,kctop0,ppbase,  &
                    pmfu,pentr,podetr,khmin,pgeoh,pdmfen,pdmfde)
    implicit none
    integer(ik4) , intent(in) :: kbdim , kk , klev , klevp1 , kproma
    integer(ik4) , dimension(kbdim) :: kcbot , kctop0 , khmin , klwmin ,   &
                                  ktype
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
    real(rkx) , dimension(kbdim,klev) :: pgeoh , pmfu , podetr , pqenh ,&
           pqte , ptenh
    intent (in) kcbot , kctop0 , khmin , klwmin , ktype , ldcum , &
                paphp1 , pentr , pgeoh , pqenh , pqte , ptenh
    intent (out) pdmfde , pdmfen , podetr
    intent (inout) pmfu , ppbase
    integer(ik4) :: ikb , ikh , iklwmin , ikt , jl
    logical :: llo1 , llo2
    real(rkx) :: zarg , zdprho , zentest , zentr , zorgde , zpmid ,     &
               zrrho , ztmzk , zzmzk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuentr'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-----------------------------------------------
    ! 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! ----------------------------------------------
    !
    ! 1.1 SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
    ! ------------------------------------------------
    !
    ! 1.2 SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
    ! ---------------------------------------------
    !
    do jl = 1 , kproma
      ppbase(jl) = paphp1(jl,kcbot(jl))
      zrrho = (rgas*ptenh(jl,kk+1)*(d_one+ep1*pqenh(jl,kk+1)))/paphp1(jl,kk+1)
      zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
      zpmid = d_half*(ppbase(jl)+paphp1(jl,kctop0(jl)))
      zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho

      llo1 = kk < kcbot(jl) .and. ldcum(jl)
      pdmfde(jl) = merge(zentr,d_zero,llo1)
      llo2 = llo1 .and. ktype(jl) == 2 .and.                            &
             (ppbase(jl)-paphp1(jl,kk) < 0.2e5_rkx .or. paphp1(jl,kk) > zpmid)
      pdmfen(jl) = merge(zentr,d_zero,llo2)
      iklwmin = max(klwmin(jl),kctop0(jl)+2)
      llo2 = llo1 .and. ktype(jl) == 3 .and. kk >= iklwmin
      if ( llo2 ) pdmfen(jl) = zentr

      if ( llo2 .and. pqenh(jl,kk+1) > 1.e-5_rkx ) then
        pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
        zentest = max(pqte(jl,kk),d_zero)/pqenh(jl,kk+1)
        zentest = min(entrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
      end if

      llo2 = llo1 .and. ktype(jl) == 1 .and.                            &
             (kk >= iklwmin .or. paphp1(jl,kk) > zpmid)
      if ( llo2 ) pdmfen(jl) = zentr
      !
      ! organized detrainment, detrainment starts at khmin
      !
      llo2 = llo1 .and. ktype(jl) == 1
      ikb = kcbot(jl)
      podetr(jl,kk) = d_zero

      if ( llo2 .and. kk <= khmin(jl) .and. kk >= kctop0(jl) ) then
        ikt = kctop0(jl)
        ikh = khmin(jl)
        if ( ikh > ikt ) then
          zzmzk = -(pgeoh(jl,ikh)-pgeoh(jl,kk))*regrav
          ztmzk = -(pgeoh(jl,ikh)-pgeoh(jl,ikt))*regrav
          zarg = mathpi*(zzmzk/ztmzk)*d_half
          zorgde = tan(zarg)*mathpi*d_half/ztmzk
          zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*(regrav*zrrho)
          podetr(jl,kk) = min(zorgde,entrmax)*pmfu(jl,kk+1)*zdprho
        end if
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuentr
  !
  ! THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
  ! FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
  !
  ! INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
  ! AND UPDRAFT VALUES T,Q ETC
  ! IT RETURNS ENTRAINMENT/DETRAINMENT RATES
  !
  ! METHOD.
  ! --------
  !   S. TIEDTKE (1989)
  !
  subroutine cuentrt(kproma,kbdim,klev,klevp1,kk,ptenh,pqenh,pqte,  &
                     paphp1,klwmin,ldcum,ktype,kcbot,kctop0,ppbase, &
                     pmfu,pentr,pdmfen,pdmfde)
    implicit none
    integer(ik4) , intent(in) :: kbdim , kk , klev , klevp1 , kproma
    integer(ik4) , dimension(kbdim) :: kcbot , kctop0 , klwmin , ktype
    logical , dimension(kbdim) :: ldcum
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
    real(rkx) , dimension(kbdim,klev) :: pmfu , pqenh , pqte , ptenh
    intent (in) kcbot , kctop0 , klwmin , ktype , ldcum , paphp1 , &
                pentr , pqenh , pqte , ptenh
    intent (out) pdmfde , pdmfen
    intent (inout) pmfu , ppbase
    integer(ik4) :: iklwmin , jl
    logical :: llo1 , llo2
    real(rkx) :: zdprho , zentest , zentr , zpmid , zrrho
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuentrt'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !-----------------------------------------------
    ! 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! ----------------------------------------------
    !
    ! 1.1 SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
    ! ------------------------------------------------
    !
    ! 1.2 SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
    ! ---------------------------------------------
    !
    do jl = 1 , kproma
      ppbase(jl) = paphp1(jl,kcbot(jl))
      zrrho = (rgas*ptenh(jl,kk+1)*(d_one+ep1*pqenh(jl,kk+1)))/paphp1(jl,kk+1)
      zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
      zpmid = d_half*(ppbase(jl)+paphp1(jl,kctop0(jl)))
      zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
      llo1 = kk < kcbot(jl) .and. ldcum(jl)
      pdmfde(jl) = merge(zentr,d_zero,llo1)
      llo2 = llo1 .and. ktype(jl) == 2 .and.                            &
             (ppbase(jl)-paphp1(jl,kk) < 0.2e5_rkx .or. paphp1(jl,kk) > zpmid)
      pdmfen(jl) = merge(zentr,d_zero,llo2)
      iklwmin = max(klwmin(jl),kctop0(jl)+2)
      llo2 = llo1 .and. (ktype(jl) == 1 .or. ktype(jl) == 3) .and.        &
             (kk >= iklwmin .or. paphp1(jl,kk) > zpmid)
      if ( llo2 ) pdmfen(jl) = zentr
      if ( llo2 .and. pqenh(jl,kk+1) > 1.e-5_rkx ) then
        pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
        zentest = max(pqte(jl,kk),d_zero)/pqenh(jl,kk+1)
        zentest = min(entrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuentrt
  !
  ! THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
  ! FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
  !
  subroutine cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ptenh,pqenh, &
                   ktrac,pxtenh,pmfuxt,pmfdxt,paphp1,pgeoh,kcbot,   &
                   kctop,kdtop,ktype,lddraf,ldcum,pmfu,pmfd,pmfus,  &
                   pmfds,pmfuq,pmfdq,pmful,pdmfup,pdmfdp,prfl,prain,&
                   pcpcu,pten,psfl,pdpmel,ktopm2)
    implicit none
    integer(ik4) , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
    integer(ik4) , intent(inout):: ktopm2
    integer(ik4) , dimension(kbdim) :: kcbot , kctop , kdtop , ktype
    logical , dimension(kbdim) :: ldcum , lddraf
    real(rkx) , dimension(kbdim,klevp1) :: paphp1
    real(rkx) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pdmfup ,      &
           pdpmel , pgeoh , pmfd , pmfdq , pmfds , pmfu , pmful ,     &
           pmfuq , pmfus , pqen , pqenh , pqsen , pten , ptenh
    real(rkx) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtenh
    real(rkx) , dimension(kbdim) :: prain , prfl , psfl
    intent (in) kcbot , kctop , kdtop , ldcum , paphp1 , pcpcu , pgeoh , &
                pqen , pqenh , pqsen , pten , ptenh , pxtenh
    intent (out) pdpmel
    intent (inout) ktype , lddraf , pdmfdp , pdmfup , pmfd , pmfdq , &
                   pmfds , pmfdxt , pmfu , pmful , pmfuq , pmfus ,   &
                   pmfuxt , prain , prfl , psfl
    integer(ik4) :: ikb , jk , jl , jt
    real(rkx) :: zcons1 , zcons2 , zcucov , zdpevap , zdrfl , zfac ,    &
               zrfl , zrfln , zrmin , zrnew , zrsum , zsnmlt ,        &
               ztmelp2 , zzp
    real(rkx) , dimension(kbdim) :: zpsubcl
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuflx'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! SPECIFY CONSTANTS
    !
    zcons1 = cpd/(wlhf*egrav*dtcum)
    zcons2 = d_one/(egrav*dtcum)
    zcucov = 0.050_rkx
    ztmelp2 = tzero + 2.0_rkx
    !
    ! 1.0 DETERMINE FINAL CONVECTIVE FLUXES
    ! -------------------------------------
    !
    ! itop=klev
    do jl = 1 , kproma
      if ( .not.ldcum(jl) .or. kdtop(jl) < kctop(jl) ) lddraf(jl) = .false.
      if ( .not.ldcum(jl) ) ktype(jl) = 0
    end do
    ktopm2 = 1
    do jk = ktopm2 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          pmfus(jl,jk) = pmfus(jl,jk) - pmfu(jl,jk)                   &
                         *(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk) = pmfuq(jl,jk) - pmfu(jl,jk)*pqenh(jl,jk)
          if ( lddraf(jl) .and. jk >= kdtop(jl) ) then
            pmfds(jl,jk) = pmfds(jl,jk) - pmfd(jl,jk)                 &
                           *(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk) = pmfdq(jl,jk) - pmfd(jl,jk)*pqenh(jl,jk)
          else
            if ( jk == 1 ) then
              call fatal(__FILE__,__LINE__, &
                      'Tiedtke Cumulus: overflow in computing fluxes')
            end if
            pmfd(jl,jk) = d_zero
            pmfds(jl,jk) = d_zero
            pmfdq(jl,jk) = d_zero
            pdmfdp(jl,jk-1) = d_zero
          end if
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            pmfuxt(jl,jk,jt) = pmfuxt(jl,jk,jt) - pmfu(jl,jk)         &
                               *pxtenh(jl,jk,jt)
            if ( lddraf(jl) .and. jk >= kdtop(jl) ) then
              pmfdxt(jl,jk,jt) = pmfdxt(jl,jk,jt) - pmfd(jl,jk)       &
                                 *pxtenh(jl,jk,jt)
            else
              pmfdxt(jl,jk,jt) = d_zero
            end if
          else
            pmfuxt(jl,jk,jt) = d_zero
            pmfdxt(jl,jk,jt) = d_zero
          end if
        end do
      end do

    end do
    do jk = ktopm2 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk > kcbot(jl) ) then
          ikb = kcbot(jl)
          zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))                    &
                /(paphp1(jl,klevp1)-paphp1(jl,ikb)))
          zzp = merge(zzp**2,zzp,ktype(jl) == 3)
          pmfu(jl,jk) = pmfu(jl,ikb)*zzp
          pmfus(jl,jk) = pmfus(jl,ikb)*zzp
          pmfuq(jl,jk) = pmfuq(jl,ikb)*zzp
          pmful(jl,jk) = pmful(jl,ikb)*zzp
        end if
      end do

      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ldcum(jl) .and. jk > kcbot(jl) ) then
            ikb = kcbot(jl)
            zzp = (paphp1(jl,klevp1)-paphp1(jl,jk))                   &
                  /(paphp1(jl,klevp1)-paphp1(jl,ikb))
            zzp = merge(zzp**2,zzp,ktype(jl) == 3)
            pmfuxt(jl,jk,jt) = pmfuxt(jl,ikb,jt)*zzp
          end if
        end do
      end do
    end do
    !
    ! 2. CALCULATE RAIN/SNOW FALL RATES
    !    CALCULATE MELTING OF SNOW
    !    CALCULATE EVAPORATION OF PRECIP
    !--- -------------------------------
    !
    do jl = 1 , kproma
      prfl(jl) = d_zero
      psfl(jl) = d_zero
      prain(jl) = d_zero
    end do
    do jk = ktopm2 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          prain(jl) = prain(jl) + pdmfup(jl,jk)
          if ( pten(jl,jk) > tzero ) then
            prfl(jl) = prfl(jl) + pdmfup(jl,jk) + pdmfdp(jl,jk)
            if ( psfl(jl) > d_zero .and. pten(jl,jk) > ztmelp2 ) then
              zfac = zcons1*(d_one+hcrm1*pqen(jl,jk))  &
                     *(paphp1(jl,jk+1)-paphp1(jl,jk))
              zsnmlt = min(psfl(jl),zfac*(pten(jl,jk)-ztmelp2))
              pdpmel(jl,jk) = zsnmlt
              psfl(jl) = psfl(jl) - zsnmlt
              prfl(jl) = prfl(jl) + zsnmlt
            end if
          else
            psfl(jl) = psfl(jl) + pdmfup(jl,jk) + pdmfdp(jl,jk)
          end if
        end if
      end do
    end do
    do jl = 1 , kproma
      prfl(jl) = max(prfl(jl),d_zero)
      psfl(jl) = max(psfl(jl),d_zero)
      zpsubcl(jl) = prfl(jl) + psfl(jl)
    end do
    do jk = ktopm2 , klev
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk >= kcbot(jl) .and. zpsubcl(jl) > almostzero ) then
          zrfl = zpsubcl(jl)
          if ( ldland(jl) ) then
            zrnew = (max(d_zero,sqrt(zrfl/zcucov)-revap_lnd* &
                     (paphp1(jl,jk+1)-paphp1(jl,jk))           &
                    *max(d_zero,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
          else
            zrnew = (max(d_zero,sqrt(zrfl/zcucov)-revap_ocn* &
                     (paphp1(jl,jk+1)-paphp1(jl,jk))           &
                    *max(d_zero,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
          end if
          zrmin = zrfl - zcucov*max(d_zero,0.80_rkx*pqsen(jl,jk)         &
                  -pqen(jl,jk))*zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
          zrnew = max(zrnew,zrmin)
          zrfln = max(zrnew,d_zero)
          zdrfl = min(d_zero,zrfln-zrfl)
          pdmfup(jl,jk) = pdmfup(jl,jk) + zdrfl
          zpsubcl(jl) = zrfln
        end if
      end do
    end do
    do jl = 1 , kproma
      zrsum = prfl(jl) + psfl(jl)
      zdpevap = zpsubcl(jl) - zrsum
      prfl(jl) = prfl(jl) + zdpevap*prfl(jl)*(d_one/max(almostzero,zrsum))
      psfl(jl) = psfl(jl) + zdpevap*psfl(jl)*(d_one/max(almostzero,zrsum))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuflx
  !
  ! TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
  !
  ! THIS ROUTINE IS CALLED FROM SUBROUTINES:
  !   *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
  !   *CUASC*    (T AND Q AT CLOUD LEVELS)
  !   *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
  ! INPUT ARE UNADJUSTED T AND Q VALUES,
  ! IT RETURNS ADJUSTED VALUES OF T AND Q
  ! NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
  !    KCALL=0    ENV. T AND QS IN*CUINI*
  !    KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
  !    KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
  !
  ! EXTERNALS
  ! ---------
  !   3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
  !     FOR CONDENSATION CALCULATIONS.
  !
  subroutine cuadjtq(kproma,kbdim,klev,kk,pp,pt,pq,ldflag,kcall)
    implicit none
    integer, intent (in) :: kcall, kk, klev, kproma, kbdim
    real(rkx), intent (in) :: pp(kbdim)
    logical, intent (in) :: ldflag(kbdim)
    real(rkx), intent (inout) :: pq(kbdim,klev), pt(kbdim,klev)
    real(rkx):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
    integer(ik4) :: isum, jl, it, it1, iv
    logical :: lo
    real(rkx):: zcond(kbdim)
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cuadjtq'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lookupoverflow = .false.

    zcond = d_zero
    !
    !---------------------------------------------------------
    ! 2. calculate condensation and adjust t and q accordingly
    ! --------------------------------------------------------
    !
    if ( kcall == 1 ) then

      do iv = 1 , 1

        isum = 0
        do jl = 1 , kproma
          if ( ldflag(jl) ) then
            it = nint(pt(jl,kk)*d_1000)
            if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
            it = max(min(it,jptlucu2),jptlucu1)
            zes = tlucua(it)/pp(jl)
            zes = min(qsmax,zes)
            lo = zes < 0.4_rkx
            zcor = d_one/(d_one-ep1*zes)
            zqsat = zes*zcor
            it1 = it+1
            it1 = max(min(it1,jptlucu2),jptlucu1)
            zqst1 = tlucua(it1)/pp(jl)
            zqst1 = min(qsmax,zqst1)
            zqst1 = zqst1/(d_one-ep1*zqst1)
            zdqsdt = (zqst1-zqsat)*d_1000
            zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
            zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
            zcond(jl) = max(zcond(jl),d_zero)
            pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
            pq(jl,kk) = pq(jl,kk)-zcond(jl)
            if ( abs(zcond(jl)) > d_zero) isum = isum+1
          end if
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif

        if ( isum == 0 ) exit

        do jl = 1 , kproma
          if ( ldflag(jl) ) then
            if ( abs(zcond(jl)) > d_zero ) then
              it = nint(pt(jl,kk)*d_1000)
              if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
              it = max(min(it,jptlucu2),jptlucu1)
              zes = tlucua(it)/pp(jl)
              zes = min(qsmax,zes)
              lo = zes < 0.4_rkx
              zcor = d_one/(d_one-ep1*zes)
              zqsat = zes*zcor
              it1 = it+1
              it1 = max(min(it1,jptlucu2),jptlucu1)
              zqst1 = tlucua(it1)/pp(jl)
              zqst1 = min(qsmax,zqst1)
              zqst1 = zqst1/(d_one-ep1*zqst1)
              zdqsdt = (zqst1-zqsat)*d_1000
              zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
              zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
              pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif
      end do

    else if ( kcall == 2 ) then

      do iv = 1 , 1

        isum = 0
        do jl = 1 , kproma
          if ( ldflag(jl) ) then
            it = nint(pt(jl,kk)*d_1000)
            if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
            it = max(min(it,jptlucu2),jptlucu1)
            zes = tlucua(it)/pp(jl)
            zes = min(qsmax,zes)
            lo = zes < 0.4_rkx
            zcor = d_one/(d_one-ep1*zes)
            zqsat = zes*zcor
            it1 = it+1
            it1 = max(min(it1,jptlucu2),jptlucu1)
            zqst1 = tlucua(it1)/pp(jl)
            zqst1 = min(qsmax,zqst1)
            zqst1 = zqst1/(d_one-ep1*zqst1)
            zdqsdt = (zqst1-zqsat)*d_1000
            zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
            zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
            zcond(jl) = min(zcond(jl),d_zero)
            pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
            pq(jl,kk) = pq(jl,kk)-zcond(jl)
            if ( abs(zcond(jl)) > d_zero ) isum = isum+1
          end if
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif

        if ( isum == 0 ) exit

        do jl = 1 , kproma
          if ( ldflag(jl) .and. abs(zcond(jl) ) > d_zero) then
            it = nint(pt(jl,kk)*d_1000)
            if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
            it = max(min(it,jptlucu2),jptlucu1)
            zes = tlucua(it)/pp(jl)
            zes = min(qsmax,zes)
            lo = zes < 0.4_rkx
            zcor = d_one/(d_one-ep1*zes)
            zqsat = zes*zcor
            it1 = it+1
            it1 = max(min(it1,jptlucu2),jptlucu1)
            zqst1 = tlucua(it1)/pp(jl)
            zqst1 = min(qsmax,zqst1)
            zqst1 = zqst1/(d_one-ep1*zqst1)
            zdqsdt = (zqst1-zqsat)*d_1000
            zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
            pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end if
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif
      end do

    else if ( kcall == 0 ) then

      do iv = 1 , 1

        isum = 0
        do jl = 1 , kproma
          it = nint(pt(jl,kk)*d_1000)
          if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
          it = max(min(it,jptlucu2),jptlucu1)
          zes = tlucua(it)/pp(jl)
          zes = min(qsmax,zes)
          lo = zes < 0.4_rkx
          zcor = d_one/(d_one-ep1*zes)
          zqsat = zes*zcor
          it1 = it+1
          it1 = max(min(it1,jptlucu2),jptlucu1)
          zqst1 = tlucua(it1)/pp(jl)
          zqst1 = min(qsmax,zqst1)
          zqst1 = zqst1/(d_one-ep1*zqst1)
          zdqsdt = (zqst1-zqsat)*d_1000
          zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
          zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
          pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
          pq(jl,kk) = pq(jl,kk)-zcond(jl)
          if ( abs(zcond(jl)) > d_zero ) isum = isum+1
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif

        ! isum = 0 !AMT fudge to only make one iteration for temporary speed fix

        if ( isum == 0 ) exit

        do jl = 1 , kproma
          it = nint(pt(jl,kk)*d_1000)
          if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
          it = max(min(it,jptlucu2),jptlucu1)
          zes = tlucua(it)/pp(jl)
          zes = min(qsmax,zes)
          lo = zes < 0.4_rkx
          zcor = d_one/(d_one-ep1*zes)
          zqsat = zes*zcor
          it1 = it+1
          it1 = max(min(it1,jptlucu2),jptlucu1)
          zqst1 = tlucua(it1)/pp(jl)
          zqst1 = min(qsmax,zqst1)
          zqst1 = zqst1/(d_one-ep1*zqst1)
          zdqsdt = (zqst1-zqsat)*d_1000
          zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
          pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do

        if (lookupoverflow) then
          call fatal(__FILE__,__LINE__, &
                     'cumulus tables lookup error: overflow')
        endif
      end do

    else if ( kcall == 4 ) then

      do jl = 1 , kproma
        it = nint(pt(jl,kk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/pp(jl)
        zes = min(qsmax,zes)
        lo = zes < 0.4_rkx
        zcor = d_one/(d_one-ep1*zes)
        zqsat = zes*zcor
        it1 = it+1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/pp(jl)
        zqst1 = min(qsmax,zqst1)
        zqst1 = zqst1/(d_one-ep1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
        zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
        pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
        pq(jl,kk) = pq(jl,kk)-zcond(jl)
      end do

      if (lookupoverflow) then
        call fatal(__FILE__,__LINE__, &
                   'cumulus tables lookup error: overflow')
      endif

      do jl = 1 , kproma
        it = nint(pt(jl,kk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/pp(jl)
        zes = min(qsmax,zes)
        lo = zes < 0.4_rkx
        zcor = d_one/(d_one-ep1*zes)
        zqsat = zes*zcor
        it1 = it+1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/pp(jl)
        zqst1 = min(qsmax,zqst1)
        zqst1 = zqst1/(d_one-ep1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
        pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk) = pq(jl,kk)-zcond1
      end do

      if (lookupoverflow) then
        call fatal(__FILE__,__LINE__, &
                   'cumulus tables lookup error: overflow')
      endif

    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cuadjtq
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
                      ntrac,qtrac,tenc,ccn)
    implicit none
    integer(ik4) , intent(in) :: np    ! number of points
    integer(ik4) , intent(in) :: nk    ! number of levels
    integer(ik4) , intent(in) :: n1    ! start point
    integer(ik4) , intent(in) :: n2    ! end point
    integer(ik4) , intent(in) :: ntrac ! number of tracers
    ! Land/sea mask flag
    logical , dimension(np) , intent(in) :: ldland
    !-----------------------------------------------------
    ! ::::::::::: Input atmospheric condition ::::::::::::
    !-----------------------------------------------------
    ! Temperature K
    real(rkx) , dimension(np,nk) , intent(in) :: t
    ! Water Vapor Specific Humidity kg/kg
    real(rkx) , dimension(np,nk) , intent(in) :: q
    ! U wind m/s
    real(rkx) , dimension(np,nk) , intent(in) :: u
    ! V wind m/s
    real(rkx) , dimension(np,nk) , intent(in) :: v
    ! Total (Liquid + Ice) Mixing ratio kg/kg
    real(rkx) , dimension(np,nk) , intent(in) :: qctot
    ! Vertical Pressure Pa/s
    real(rkx) , dimension(np,nk) , intent(in) :: omega
    ! Turbolent moisture flux kg/(m^2*s)
    real(rkx) , dimension(np,nk+1) , intent(in) :: qhfl
    ! Turbolent sensible heat flux W/(m^2*s)
    real(rkx) , dimension(np,nk+1) , intent(in) :: ahfs
    ! Pressure Pa
    real(rkx) , dimension(np,nk) , intent(in) :: ph
    ! Pressure on sigma levels Pa
    real(rkx) , dimension(np,nk+1) , intent(in) :: pf
    ! Geopotential m^2/s^2
    real(rkx) , dimension(np,nk) , intent(in) :: geo
    ! Geopotential m^2/s^2 on sigma levels
    real(rkx) , dimension(np,nk+1) , intent(in) :: geof
    ! Tracer Mixing ratio kg/kg
    real(rkx) , dimension(np,nk,ntrac) , intent(in) :: qtrac
    ! Cloud condensation nuclei #
    real(rkx) , dimension(np,nk) , intent(in) :: ccn
    !-----------------------------------------------------
    ! ::::::::::::::: Output Tendencies ::::::::::::::::::
    !-----------------------------------------------------
    ! Temperature  tendency K/s
    real(rkx) , dimension(np,nk) , intent(inout) :: tent
    ! specific humidity tendency kg/(kg*s)
    real(rkx) , dimension(np,nk) , intent(inout) :: tenq
    ! U wind m/s^2
    real(rkx) , dimension(np,nk) , intent(inout) :: tenu
    ! V wind m/s^2
    real(rkx) , dimension(np,nk) , intent(inout) :: tenv
    ! Tracer mixing ratio tendency kg/(kg*s)
    real(rkx) , dimension(np,nk,ntrac) , intent(inout) :: tenc
    ! Cloud liquid water mixing ratio tendency kg/(kg*s)
    real(rkx) , dimension(np,nk) , intent(out) :: tenl
    ! Ice water mixing ratio tendency kg/(kg*s)
    real(rkx) , dimension(np,nk) , intent(out) :: teni
    ! Detrained liquid water kg/(m^2*s)
    real(rkx) , dimension(np,nk) , intent(out) :: lude
    ! Cumulus activation flag
    logical , dimension(np) , intent(out) :: ldcum
    ! Type of cumulus convection
    integer(ik4) , dimension(np) , intent(out) :: ktype
    ! Bottom level of cumulus cloud
    integer(ik4) , dimension(np) , intent(out) :: kcbot
    ! Top level of cumulus cloud
    integer(ik4) , dimension(np) , intent(out) :: kctop
    ! Bottom level of shallow cumulus cloud
    integer(ik4) , dimension(np) , intent(out) :: kbotsc
    ! Shallow cumulus flag
    logical , dimension(np) , intent(out) :: ldsc
    ! Updraft temperature K
    real(rkx) , dimension(np,nk) , intent(out) :: tu
    ! Updraft Water Vapor specific humidity kg/kg
    real(rkx) , dimension(np,nk) , intent(out) :: qu
    ! Updraft Cloud Liquid Water mixing ratio kg/kg
    real(rkx) , dimension(np,nk) , intent(out) :: lu
    ! Rain Mass Flux kg/(m^2*s)
    real(rkx) , dimension(np,nk+1) , intent(out) :: mflxr
    ! Snow Mass Flux kg/(m^2*s)
    real(rkx) , dimension(np,nk+1) , intent(out) :: mflxs
    ! Total precipitation produced in convective updrafts
    real(rkx) , dimension(np) , intent(out) :: rain
    ! Massflux updrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) , intent(out) :: mfu
    ! Massflux downdrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) , intent(out) :: mfd
    ! Massflux updrafts rate kg/(m^3*s)
    real(rkx) , dimension(np,nk) , intent(out) :: mfude_rate
    ! Massflux downdrafts rate kg/(m^3*s)
    real(rkx) , dimension(np,nk) , intent(out) :: mfdde_rate
    ! Convectve available potential energy J/kg
    real(rkx) , dimension(np) , intent(out) :: cape
    !-----------------------------------------------------
    ! ::::::::::::::::: Local Variables ::::::::::::::::::
    !-----------------------------------------------------
    ! Vertically averaged updraught velocity m/s
    real(rkx) , dimension(np) :: wmean
    ! Increment of dry static energy J/(kg*s)
    real(rkx) , dimension(np,nk) :: penth
    ! Saturation mixing ratio kg/kg
    real(rkx) , dimension(np,nk) :: qs
    ! Temperature on full sigma levels (where massfluxes are computed) K
    real(rkx) , dimension(np,nk) :: tf
    ! Mixing ratio on full sigma levels kg/kg
    real(rkx) , dimension(np,nk) :: qf
    ! Saturation mixing ratio on full sigma levels kg/kg
    real(rkx) , dimension(np,nk) :: qsf
    ! Temperature in Downdrafts K
    real(rkx) , dimension(np,nk) :: td
    ! Mixing ratio in downdrafts kg/kg
    real(rkx) , dimension(np,nk) :: qd
    ! Flux of dry static energy in updrafts J/(m^2*s)
    real(rkx) , dimension(np,nk) :: mfus
    ! Flux of dry static energy in downdrafts J/(m^2*s)
    real(rkx) , dimension(np,nk) :: mfds
    ! Flux of humidity in updrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: mfuq
    ! Flux of humidity in downdrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: mfdq
    ! Flux difference of precipitation in updrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: dmfup
    ! Flux difference of precipitation in downdrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: dmfdp
    ! Flux of liquid water in updrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: mful
    ! Updraft winds m/s
    real(rkx) , dimension(np,nk) :: uu , vu
    ! Downdraft winds m/s
    real(rkx) , dimension(np,nk) :: ud , vd
    ! Updraft kinetic energy m^2/s^2
    real(rkx) , dimension(np,nk) :: kineu
    ! Downdraft kinetic energy m^2/s^2
    real(rkx) , dimension(np,nk) :: kined
    ! Precipitation Flux kg/(m^2*s)
    real(rkx) , dimension(np) :: rfl
    ! Mass Flux Updraft base of the cloud
    real(rkx) , dimension(np) :: mfub , mfub1
    ! Change in precipitation fluxes due to melting kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: dpmel
    ! Precipitation flux
    real(rkx) , dimension(np,nk) :: xrain
    ! Flux of frozen cloudwater in updrafts kg/(m^2*s)
    real(rkx) , dimension(np,nk) :: lglac
    real(rkx) , dimension(np) :: dqcv
    real(rkx) , dimension(np) :: dhpbl
    real(rkx) , dimension(np) :: wubase
    ! Entrainment , detrainment
    real(rkx) , dimension(np,nk) :: dmfen , dmfde
    ! Flag
    integer(ik4) , dimension(np,nk) :: ilab
    ! Level of max wind
    integer(ik4) , dimension(np) :: ilwmin
    integer(ik4) , dimension(np) :: idtop , ictop0
    ! Departure level for convection
    integer(ik4) , dimension(np) :: idpl
    real(rkx) , dimension(np) :: xcape , rheat
    logical , dimension(np) :: lddraf , llddraf3 , lldcum , llo2
    logical :: llo1
    integer(ik4) :: ikb , itopm2 , k , ik , n
    real(rkx) :: cons2 , cons , dh , dqmin , dz , eps , fac , &
                 mfmax , pbmpt , qumqe , xro , mfa , erate ,  &
                 derate , duten , dvten , tdis , xalv
    real(rkx) , dimension(np) :: sfl
    real(rkx) , dimension(np) :: tau ! adjustment time
    ! scaling factor for momentum and tracer massflux
    real(rkx) , dimension(np) :: mfs , mfuub , mfuvb , xsum12 , xsum22
    real(rkx) , dimension(np) :: mf_shal
    real(rkx) , dimension(np,nk) :: mfuus , mfdus , mfudr , xtent , &
                                    mfddr , xtenu , xtenv , xtenq , &
                                    uv2
    logical :: llconscheck = .false.
    integer(ik4) :: nt
    real(rkx) , dimension(:,:) , allocatable :: xsumc
    real(rkx) , dimension(:,:,:) , allocatable :: xtenc

    ! ---------------------------------------
    ! 0. Compute Saturation specific humidity
    ! ---------------------------------------
    ldcum(:) = .false.
    qs(:,:) = q(:,:)
    call satur(nk060,ph,t,qs)
    !------------------------------------
    ! 1. Specify constants and parameters
    ! -----------------------------------
    cons2 = rmfcfl/(egrav*dtcum)
    cons = d_one/(egrav*dtcum)
    !---------------------------------------------
    ! 2. Initialize values at vertical grid points
    ! --------------------------------------------
    do k = 1 , nk
      do n = n1 , n2
        xtent(n,k) = tent(n,k)
        xtenq(n,k) = tenq(n,k)
        xtenu(n,k) = tenu(n,k)
        xtenv(n,k) = tenv(n,k)
      end do
    end do
    call initcum
    !---------------------------
    ! 3. Cloud base calculations
    ! --------------------------
    ! ---------------------------------
    !   (A) Determine cloud base values
    ! ---------------------------------
    call cloudbase
    !-----------------------------------------------------------
    !   (B) Determine total moisture convergence and
    !       decide on type of cumulus convection on the basis of
    !       the depth of the convection
    !           deep    if cloud depth > 200mb
    !           shallow if cloud depth < 200mb
    ! ----------------------------------------------------------
    ! Calculate column and sub cloud layer moisture convergence
    ! and sub cloud layer moist static energy convergence
    do n = n1 , n2
      dqcv(n) = d_zero
      dhpbl(n) = d_zero
      idtop(n) = 0
    end do
    do k = nk060 , nk
      do n = n1 , n2
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
    do n = n1 , n2
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
      do n = n1 , n2
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
    do n = n1 , n2
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
    call ascent
    !  (C) Check cloud depth and change entrainment rate accordingly
    !      Calculate precipitation rate (for downdraft calculation)
    ! -----------------------------------------------------------------------
    do n = n1 , n2
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
    do k = 2 , nk
      do n = n1 , n2
        rfl(n) = rfl(n) + dmfup(n,k)
      end do
    end do
    do k = 1 , nk
      do n = n1 , n2
        mfd(n,k) = d_zero
        mfds(n,k) = d_zero
        mfdq(n,k) = d_zero
        dmfdp(n,k) = d_zero
        dpmel(n,k) = d_zero
      end do
    end do
    !----------------------------------
    ! 5. Cumulus downdraft calculations
    ! ---------------------------------
    if ( lmfdd ) then
      !-------------------
      ! (A) Determine lfs
      ! --=---------------
      call lfs
      !---------------------------------------
      ! (B) Determine downdraft t,q and fluxes
      ! --------------------------------------
      call ddrafdsc
    end if
    !-----------
    ! 6. Closure
    ! ----------
    !   Recalculate cloud base massflux from a cape closure for deep
    !   convection (ktype=1) and by PBL equilibrum taking downdrafts
    !   into account for shallow convection (ktype=2)
    ! --------------------------------------------------------------
    ! Deep convection
    do n = n1 , n2
      rheat(n) = d_zero
      xcape(n) = d_zero
      mfub1(n) = mfub(n)
    end do
    do k = 1 , nk
      do n = n1 , n2
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
    do n = n1 , n2
      if ( ldcum(n) .and. ktype(n) == 1 ) then
        ikb = kcbot(n)
        ik = kctop(n)
        xcape(n) = max(d_zero,min(xcape(n),5000.0_rkx))
        rheat(n) = max(1.e-4_rkx,rheat(n))
        tau(n) = (geof(n,ik)-geof(n,ikb)) / &
                   ((d_two+min(15.0_rkx,wmean(n)))*egrav)*rtau
        tau(n) = max(dtcum,min(10800.0_rkx,tau(n)))
        tau(n) = max(720.0_rkx,tau(n))
        mfub1(n) = (xcape(n)*mfub(n))/(rheat(n)*tau(n))
        mfub1(n) = max(mfub1(n),0.001_rkx)
        mfmax = (pf(n,ikb)-pf(n,ikb-1))*cons2*rmflic + rmflia
        mfub1(n) = min(mfub1(n),mfmax)
      end if
    end do
    ! Shallow convection and mid_level
    do n = n1 , n2
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
    do k = 1 , nk
      do n = n1 , n2
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
    end do
    !--------------------------
    ! 6.1 Final closure=scaling
    ! -------------------------
    do n = n1 , n2
      if ( ldcum(n) ) mfs(n) = mfub1(n)/max(cmfcmin,mfub(n))
    end do
    do k = 2 , nk
      do n = n1 , n2
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
    do k = 2 , nk
      do n = n1 , n2
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
    end do
    !--------------------------------------------------------
    ! 6.2 In case that either deep or shallow is switched off
    !      reset ldcum to false-> fluxes set to zero
    ! -------------------------------------------------------
    ! exclude pathological ktype=2 kcbot=kctop=nk-1
    do n = n1 , n2
      if ( ktype(n) == 2 .and.        &
           kcbot(n) == kctop(n) .and. &
           kcbot(n) >= nk-1 ) then
        ldcum(n) = .false.
        ktype(n) = 0
      end if
    end do
    if ( .not. lmfscv .or. .not. lmfpen ) then
      do n = n1 , n2
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
    do n = n1 , n2
      if ( lddraf(n) .and. idtop(n) <= kctop(n) ) then
        idtop(n) = kctop(n) + 1
      end if
    end do
    do k = 2 , nk
      do n = n1 , n2
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
    end do
    call cfluxes
    !- rescale DD fluxes if total mass flux becomes negative
    !- correct DD detrainment rates if entrainment becomes negative
    !- correct UD detrainment rates if entrainment becomes negative
    !- conservation correction for precip
    mfs(:) = d_one
    do k = 2 , nk
      ! change for stability
      do n = n1 , n2
        if ( lddraf(n) .and. k >= idtop(n)-1 ) then
          mfmax = mfu(n,k)*0.98_rkx
          if ( mfd(n,k)+mfmax+1.e-15_rkx < d_zero ) then
            mfs(n) = min(mfs(n),-mfmax/mfd(n,k))
          end if
        end if
      end do
    end do
    mfuub(:) = d_zero
    do k = 2 , nk
      do n = n1 , n2
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
    do k = 2 , nk - 1
      do n = n1 , n2
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
    end do
    ! Avoid negative humidities at ddraught top
    do n = n1 , n2
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
    do k = 2 , nk
      do n = n1 , n2
        if ( ldcum(n) .and. k >= kctop(n)-1 .and. k < kcbot(n) ) then
          dz = dtcum*egrav/(pf(n,k+1)-pf(n,k))
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
    end do
    if ( llconscheck ) then
      if ( lmftrac .and. ntrac > 0 ) then
        allocate (xtenc(np,nk,ntrac))
        allocate (xsumc(np,4+ntrac))
        do nt = 1 , ntrac
          do k = 2 , nk
            do n = n1 , n2
              if ( ldcum(n) ) xtenc(n,k,nt) = tenc(n,k,nt)
            end do
          end do
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
      do k = nk , 2 , -1
        do n = n1 , n2
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
      end do
    end if
    call dtdqc
    !-------------------------------------------------
    ! 9. Compute momentum in updraught and downdraught
    ! ------------------------------------------------
    if ( lmfdudv ) then
      do k = nk - 1 , 2 , -1
        ik = k + 1
        do n = n1 , n2
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
      do k = 3 , nk
        ik = k - 1
        do n = n1 , n2
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
      mfs(:) = d_one
      if ( rmfsoluv <= d_one ) then
        do k = 2 , nk
          do n = n1 , n2
            if ( ldcum(n) .and. k >= kctop(n)-1 ) then
              mfmax = (pf(n,k)-pf(n,k-1))*cons
              if ( mfu(n,k) > mfmax .and. k >= kctop(n) ) then
                mfs(n) = min(mfs(n),mfmax/mfu(n,k))
              end if
            end if
          end do
        end do
      end if
      do k = 1 , nk
        do n = n1 , n2
          mfuus(n,k) = mfu(n,k)
          mfdus(n,k) = mfd(n,k)
          if ( ldcum(n) .and. k >= kctop(n)-1 ) then
            mfuus(n,k) = mfu(n,k)*mfs(n)
            mfdus(n,k) = mfd(n,k)*mfs(n)
          end if
        end do
      end do
      ! Recompute Draught properties below for Implicit
      ! based on linear flux profiles
      if ( rmfsoluv > d_zero ) then
        do n = n1 , n2
          if ( ldcum(n) ) then
            k = kcbot(n)
            ik = k - 1
            mfuub(n) = mfuus(n,k)*(uu(n,k)-u(n,ik))
            mfuvb(n) = mfuus(n,k)*(vu(n,k)-v(n,ik))
          end if
        end do
        do k = 2 , nk
          ik = k - 1
          do n = n1 , n2
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
        do n = n1 , n2
          xsum12(n) = d_zero
          xsum22(n) = d_zero
        end do
        do k = 1 , nk
          do n = n1 , n2
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
        do k = 1 , nk
          do n = n1 , n2
            if ( ldcum(n) .and. k>=kctop(n)-1 ) then
              dz = (pf(n,k+1)-pf(n,k))
              tdis = rcpd*xsum12(n)*uv2(n,k)/max(1.e-15_rkx,xsum22(n))
              tent(n,k) = tent(n,k) + tdis
            end if
          end do
        end do
      end if
    end if
    !--------------------------------------------------------
    ! 10. In case that either deep or shallow is switched off
    !     need to set some variables a posteriori to zero
    !--------------------------------------------------------
    if ( .not. lmfscv .or. .not. lmfpen ) then
      do k = 2 , nk
        do n = n1 , n2
          if ( llo2(n) .and. k >= kctop(n)-1 ) then
            tu(n,k) = t(n,k)
            qu(n,k) = q(n,k)
            lu(n,k) = d_zero
            penth(n,k) = d_zero
            mfude_rate(n,k) = d_zero
            mfdde_rate(n,k) = d_zero
          end if
        end do
      end do
      do n = n1 , n2
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
      do n = n1 , n2
        if ( ldcum(n) .and. ktype(n) /= 3 .and. kcbot(n)-kctop(n) >= 1 ) then
          lldcum(n) = .true.
          llddraf3(n) = lddraf(n)
        else
          lldcum(n) = .false.
          llddraf3(n) = .false.
        end if
      end do
      ! check and correct mass fluxes for CFL criterium
      mfs(:) = d_one
      if ( rmfsolct <= d_three ) then
        do k = 2 , nk
          do n = n1 , n2
            if ( lldcum(n) .and. k >= kctop(n) ) then
              mfmax = (pf(n,k)-pf(n,k-1))*0.8_rkx*cons
              if ( mfu(n,k) > mfmax ) then
                mfs(n) = min(mfs(n),mfmax/mfu(n,k))
              end if
            end if
          end do
        end do
      end if
      do k = 1 , nk
        do n = n1 , n2
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
      end do
      if ( lmfsmooth ) then
        ! smoothing of mass fluxes (gradients) at top and bottom of draughts
        do k = 2 , nk - 1
          do n = n1 , n2
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
        do k = nk - 1 , 2 , -1
          do n = n1 , n2
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
    do k = 2 , nk
      do n = n1 , n2
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
    end do

    if ( llconscheck ) then
      !--------------------------------------
      ! 13. Conservation check AND correction
      ! -------------------------------------
      do n = n1 , n2
        xsumc(n,:) = d_zero
      end do
      do k = nk , 2 , -1
        do n = n1 , n2
          if ( ldcum(n) .and. k >= kctop(n)-1 ) then
            dz = (pf(n,k+1)-pf(n,k))*regrav
            xsumc(n,1) = xsumc(n,1) + &
                         (tenq(n,k)-xtenq(n,k))*dz + lude(n,k)
            xalv = mlw(t(n,k))
            xsumc(n,2) = xsumc(n,2) + &
                         cpd*(tent(n,k)-xtent(n,k))*dz - xalv*lude(n,k)
            xsumc(n,3) = xsumc(n,3) + (tenu(n,k)-xtenu(n,k))*dz
            xsumc(n,4) = xsumc(n,4) + (tenv(n,k)-xtenv(n,k))*dz
          end if
        end do
      end do
      if ( lmftrac .and. ntrac > 0 ) then
        do nt = 1 , ntrac
          do k = nk , 2 , -1
            do n = n1 , n2
              if ( ldcum(n) .and. k >= kctop(n)-1 ) then
                dz = (pf(n,k+1)-pf(n,k))*regrav
                xsumc(n,4+nt) = xsumc(n,4+nt) + (tenc(n,k,nt)-xtenc(n,k,nt))*dz
              end if
            end do
          end do
        end do
      end if
      do n = n1 , n2
        if ( ldcum(n) ) then
          xalv = mlw(t(n,nk))
          sfl(n) = mflxr(n,nk+1) + mflxs(n,nk+1)
          write (61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')     &
            n , ' CONS q: ' , -xsumc(n,1)*xalv , sfl(n)*xalv ,   &
            ktype(n) , ' CONS h: ' , xsumc(n,2) , ' CONS uv: ' , &
            xsumc(n,3) , xsumc(n,4)
          if ( lmftrac .and. ntrac > 0 ) then
            write (61,*) ' Conserv Error Tracers 1-' , ntrac , ' :'
            do nt = 1 , ntrac
              write (61,'(i4,e12.4)') nt , xsumc(n,4+nt)
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
    do k = 1 , nk
      do n = n1 , n2
        tenl(n,k) = lude(n,k)*egrav/(pf(n,k+1)-pf(n,k))
        teni(n,k) = (d_one-xalpha(t(n,k)))*tenl(n,k)
        tenl(n,k) = tenl(n,k) - teni(n,k)
        mflxr(n,k) = mflxr(n,k)*1.e-3_rkx
        mflxs(n,k) = mflxs(n,k)*1.e-3_rkx
      end do
    end do
    do n = n1 , n2
      mflxr(n,nk+1) = mflxr(n,nk+1)*1.e-3_rkx
      mflxs(n,nk+1) = mflxs(n,nk+1)*1.e-3_rkx
    end do

    contains
    !
    ! This routine interpolates large-scale fields of t,q etc.
    ! to sigma levels (i.e. grid for massflux scheme), determines level
    ! of maximum vertical velocity and initializes values for updrafts
    ! and downdrafts interface
    !
    subroutine initcum
      implicit none
      real(rkx) , dimension(np) :: wmax
      real(rkx) , dimension(np) :: xph
      logical , dimension(np) :: llflag
      integer(ik4) :: icall , ik , k , n
      real(rkx) :: zs
      !--------------------------------------------------
      ! 1. Specify large scale parameters at half levels
      !    adjust temperature fields if staticly unstable
      !    find level of maximum vertical velocity
      ! -------------------------------------------------
      do k = 2 , nk
        do n = n1 , n2
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
        do n = n1 , n2
          qf(n,k) = min(q(n,k-1),qs(n,k-1))+(qsf(n,k)-qs(n,k-1))
          qf(n,k) = max(qf(n,k),d_zero)
        end do
      end do
      do n = n1 , n2
        tf(n,nk) = (cpd*t(n,nk)+geo(n,nk)-geof(n,nk))*rcpd
        qf(n,nk) = q(n,nk)
        tf(n,1) = t(n,1)
        qf(n,1) = q(n,1)
        ilwmin(n) = nk
        wmax(n) = d_zero
      end do
      do k = nk - 1 , 2 , -1
        do n = n1 , n2
          zs = max(cpd*tf(n,k)+geof(n,k),cpd*tf(n,k+1)+geof(n,k+1))
          tf(n,k) = (zs-geof(n,k))*rcpd
        end do
      end do
      do k = nk , 3 , -1
        do n = n1 , n2
          if ( omega(n,k) < wmax(n) ) then
            wmax(n) = omega(n,k)
            ilwmin(n) = k
          end if
        end do
      end do
      !-------------------------------------------------
      ! 2. Initialize values for updrafts and downdrafts
      ! ------------------------------------------------
      do k = 1 , nk
        ik = k - 1
        if ( k==1 ) ik = 1
        do n = n1 , n2
          tu(n,k) = tf(n,k)
          td(n,k) = tf(n,k)
          qu(n,k) = qf(n,k)
          qd(n,k) = qf(n,k)
          lu(n,k) = d_zero
          uu(n,k) = u(n,ik)
          ud(n,k) = u(n,ik)
          vu(n,k) = v(n,ik)
          vd(n,k) = v(n,ik)
          ilab(n,k) = 0
        end do
      end do
    end subroutine initcum
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
    subroutine solver(ktop,wmask,a,b,r,u)
      implicit none
      integer(ik4) , dimension(np) , intent(in) :: ktop ! cloud top
      logical , dimension(np,nk) , intent(in) :: wmask  ! mask
      real(rkx) , dimension(np,nk) , intent(in) :: a    ! lower diagonal
      real(rkx) , dimension(np,nk) , intent(in) :: b    ! diagonal
      real(rkx) , dimension(np,nk) , intent(in) :: r    ! right hand side
      real(rkx) , dimension(np,nk) , intent(out) :: u   ! solution
      integer(ik4) :: k , n
      u(:,:) = d_zero
      do k = 2 , nk
        do n = n1 , n2
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
    !
    ! Computes specific humidity at saturation used by the diagnostic
    ! cloud scheme to compute relative humidity and liquid water content
    !
    subroutine satur(kt,pr,t,qsat)
      implicit none
      integer(ik4) , intent(in) :: kt                    ! top level
      real(rkx) , dimension(np,nk) , intent(in) :: pr    ! Pressure Pa
      real(rkx) , dimension(np,nk) , intent(in) :: t     ! Temperature K
      real(rkx) , dimension(np,nk) , intent(out) :: qsat ! Satur. MXR kg/kg
      integer(ik4) :: k , n
      real(rkx) :: qs
      do k = kt , nk
        do n = n1 , n2
          qs = min(qsmax,fesat(t(n,k))/pr(n,k))
          qsat(n,k) = qs/(d_one-ep1*qs)
        end do
      end do
    end subroutine satur
    !
    ! Produce adjusted t,q and l values
    !
    subroutine moistadj(kk,sp,t,q,ldflag,jcall)
      implicit none
      integer(ik4) , intent(in) :: kk ! actual level
      real(rkx) , dimension(np) , intent(in) :: sp ! Surface pressure Pa
      real(rkx) , dimension(np,nk) , intent(inout) :: t
      real(rkx) , dimension(np,nk) , intent(inout) :: q
      logical , dimension(np) , intent(in) :: ldflag ! Active convection flag
      integer(ik4) , intent(in) :: jcall ! Calling method
      integer(ik4) :: jl
      real(rkx) :: cond , cond1 , cor , qs , rp
      real(rkx) :: zl , zi , zf
      !---------------------------------------------------------
      ! 1. Calculate condensation and adjust t and q accordingly
      ! --------------------------------------------------------
      if ( jcall == 1 ) then
        do jl = n1 , n2
          if ( ldflag(jl) ) then
            rp = d_one/sp(jl)
            zl = d_one/(t(jl,kk)-c4les)
            zi = d_one/(t(jl,kk)-c4ies)
            qs = c2es*(xalpha(t(jl,kk))*exp(c3les*(t(jl,kk)-tzero)*zl) + &
                  (d_one-xalpha(t(jl,kk)))*exp(c3ies*(t(jl,kk)-tzero)*zi))
            qs = qs*rp
            qs = min(qsmax,qs)
            cor = d_one - ep1*qs
            zf = xalpha(t(jl,kk))*c5alvcp*zl**2 + &
                 (d_one-xalpha(t(jl,kk)))*c5alscp*zi**2
            cond = (q(jl,kk)*cor**2-qs*cor)/(cor**2+qs*zf)
            if ( cond > d_zero ) then
              t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond
              q(jl,kk) = q(jl,kk) - cond
              zl = d_one/(t(jl,kk)-c4les)
              zi = d_one/(t(jl,kk)-c4ies)
              qs = c2es*(xalpha(t(jl,kk)) * &
                exp(c3les*(t(jl,kk)-tzero)*zl)+(d_one-xalpha(t(jl,kk))) * &
                exp(c3ies*(t(jl,kk)-tzero)*zi))
              qs = qs*rp
              qs = xmin(qsmax,qs)
              cor = d_one - ep1*qs
              zf = xalpha(t(jl,kk))*c5alvcp*zl**2 + &
                   (d_one-xalpha(t(jl,kk)))*c5alscp*zi**2
              cond1 = (q(jl,kk)*cor**2-qs*cor)/(cor**2+qs*zf)
              if ( abs(cond) < almostzero ) cond1 = d_zero
              t(jl,kk) = t(jl,kk) + mlwocp(t(jl,kk))*cond1
              q(jl,kk) = q(jl,kk) - cond1
            end if
          end if
        end do
      else if ( jcall == 2 ) then
        do jl = n1 , n2
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
        do jl = n1 , n2
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
        do jl = n1 , n2
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
        do jl = n1 , n2
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
        do jl = n1 , n2
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
    ! This routine does the calculations for cloud ascents, i.e. produce
    ! cloud ascents vertical profiles of T,Q,L,U AND V and corresponding
    ! fluxes as well as precipitation rates.
    ! Lift surface air dry-adiabatically to cloud base and then calculate
    ! moist ascent for entraining/detraining plume.
    ! Entrainment and detrainment rates differ for shallow and deep cumulus
    ! convection. In case there is no penetrative or shallow convection
    ! check for possibility of mid level convection.
    !
    subroutine ascent
      implicit none
      real(rkx) , dimension(np,nk) :: buo
      real(rkx) , dimension(np) :: ndmfen , ndmfde , qold , luold , precip
      real(rkx) , dimension(np) :: dpmean
      real(rkx) , dimension(np) :: zoentr , xph
      logical , dimension(np) :: llflag , llflaguv , llo1
      logical :: llo3
      integer(ik4) :: icall , ik , is , k , n , ikb
      integer(ik4) :: nll , nlm
      integer(ik4) , dimension(np) :: nlx
      real(rkx) :: cldmax , cprc2 , cwifrac , alfaw , bc ,   &
                   be , buoc , zc , cbf , cons2 , zd , zdfi , dkbuo , &
                   dken , dnoprc , zdt , fac , facbuo , zint ,       &
                   kedke , lcrit , leen , lnew , mfmax , mftest ,     &
                   mfulk , mfun , mfuqk , mfusk , prcdgw , prcon ,    &
                   qeen , qude , rnew , rold , scde , seen , &
                   zvi , zvv , zvw , zwu , zco , arg , entrpen
      real(rkx) :: change , zxs , zxe
      logical , dimension(np) :: llklab
      !----------------------
      ! 1. Specify parameters
      ! ---------------------
      cons2 = rmfcfl/(egrav*dtcum)
      facbuo = d_half/(d_one+d_half)
      cldmax = 5.e-3_rkx
      cwifrac = d_half
      cprc2 = d_half
      !----------------------
      ! 2. Set default values
      ! ---------------------
      llo3 = .false.
      do n = n1 , n2
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
      do n = n1 , n2
        llklab(n) = .false.
        if ( .not. ldcum(n) .or. ktype(n) == 3 ) llklab(n) = .true.
      end do
      do k = 1 , nk
        do n = n1 , n2
          if ( k /= kcbot(n) ) lu(n,k) = d_zero
          kineu(n,k) = d_zero
        end do
        do n = n1 , n2
          mfu(n,k) = d_zero
          mfus(n,k) = d_zero
          mfuq(n,k) = d_zero
          mful(n,k) = d_zero
        end do
        do n = n1 , n2
          lude(n,k) = d_zero
          lglac(n,k) = d_zero
          dmfup(n,k) = d_zero
          xrain(n,k) = d_zero
        end do
        do n = n1 , n2
          buo(n,k) = d_zero
          if ( llklab(n) ) ilab(n,k) = 0
          if ( .not. ldcum(n) .and. ph(n,k) < 4.0e4_rkx ) ictop0(n) = k
          dmfen(n,k) = d_zero
          mfude_rate(n,k) = d_zero
        end do
      end do
      do n = n1 , n2
        if ( ktype(n) == 3 ) ldcum(n) = .false.
      end do
      !-----------------------------------------
      ! 3. Initialize values at cloud base level
      ! ----------------------------------------
      do n = n1 , n2
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
      do k = nk - 1 , 3 , -1
        ! Specify cloud base values for midlevel convection
        ! in case there is not already convection
        ! ----------------------------------------------------
        ik = k
        call mcbase(ik)
        is = 0
        nlm = 0
        do n = n1 , n2
          llflag(n) = .false.
          precip(n) = d_zero
          llo1(n) = .false.
          is = is + ilab(n,k+1)
          if ( ilab(n,k+1) == 0 ) ilab(n,k) = 0
          if ( (ldcum(n) .and. ilab(n,k+1) == 2) .or.   &
               (ktype(n) == 3 .and. ilab(n,k+1) == 1) ) then
            llflag(n) = .true.
            nlm = nlm + 1
            nlx(nlm) = n
          end if
          if ( ilab(n,k+1) > 0 ) then
            llflaguv(n) = .true.
          else
            llflaguv(n) = .false.
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
          do n = n1 , n2
            qold(n) = d_zero
          end do
          do nll = 1 , nlm
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
          do n = n1 , n2
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
          do nll = 1 , nlm
            n = nlx(nll)
            if ( qu(n,k) /= qold(n) ) then
              lglac(n,k) = lu(n,k) * ((d_one-xalpha(tu(n,k)))- &
                                     (d_one-xalpha(tu(n,k+1))))
              tu(n,k) = tu(n,k) + wlhfocp*lglac(n,k)
            end if
          end do
          do nll = 1 , nlm
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
          do n = n1 , n2
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
              end if
            end if
          end do
          do n = n1 , n2
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
          do nll = 1 , nlm
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
      do n = n1 , n2
        if ( kctop(n) == -1 ) ldcum(n) = .false.
        kcbot(n) = max(kcbot(n),kctop(n))
        if ( ldcum(n) ) then
          wmean(n) = max(1.e-2_rkx,wmean(n)/max(d_one,dpmean(n)))
          wmean(n) = sqrt(d_two*wmean(n))
        end if
      end do
    end subroutine ascent
    !
    ! This routine calculates entrainment/detrainment rates for updrafts
    ! in cumulus parameterization. Input are environmental values t,q etc
    ! and updraft values t,q etc. It returns entrainment/detrainment rates
    ! Turbulent entrainment is simulated by a constant multiplied by a
    ! vertical scaling function
    !
    subroutine entrainm(kk,ldwork,dmfen,dmfde)
      implicit none
      integer(ik4) , intent(in) :: kk ! input level
      logical , intent(in) :: ldwork                      ! enable/disable
      real(rkx) , dimension(np) , intent(out) :: dmfen    ! Entrainment
      real(rkx) , dimension(np) , intent(out) :: dmfde    ! Detrainment
      integer(ik4) :: n
      real(rkx) :: mf , entrpen , detrpen
      !-----------------------------------------------
      ! 1. Calculate entrainment and detrainment rates
      ! ----------------------------------------------
      if ( ldwork ) then
        do n = n1 , n2
          dmfen(n) = d_zero
          dmfde(n) = d_zero
        end do
        do n = n1 , n2
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
      integer(ik4) , intent(in) :: kk                     ! current level
      integer(ik4) :: n
      !-----------------------------------------------
      ! 1. Calculate entrainment and detrainment rates
      ! ----------------------------------------------
      do n = n1 , n2
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
    ! This routine calculates level of free sinking for cumulus downdrafts
    ! and specifies t,q,u and v values to produce lfs-values for cumulus
    ! downdrafts for massflux cumulus parameterization interface
    ! Check for negative buoyancy of air of equal parts of moist
    ! environmental air and cloud air.
    !
    subroutine lfs
      implicit none
      integer(ik4) , dimension(np) :: ikhsmin
      real(rkx) , dimension(np,nk) :: tenwb , qenwb
      real(rkx) , dimension(np) :: cond , xph , hsmin
      logical , dimension(np) :: llo2
      integer(ik4) :: icall , ik , ike , is , k , n
      real(rkx) :: buo , hsk , mftop , qtest , ttest
      !-------------------------------------
      ! 1. Set default values for downdrafts
      ! ------------------------------------
      do n = n1 , n2
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
        do k = 3 , nk - 2
          do n = n1 , n2
            hsk = cpd*t(n,k) + geo(n,k) + mlw(t(n,k))*qs(n,k)
            if ( hsk < hsmin(n) ) then
              hsmin(n) = hsk
              ikhsmin(n) = k
            end if
          end do
        end do
        ike = nk - 3
        do k = 3 , ike
          ! -----------------------------------------------
          ! 2.1 Calculate wet-bulb temperature and moisture
          !     for environmental air
          ! -----------------------------------------------
          is = 0
          do n = n1 , n2
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
          do n = n1 , n2
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
    end subroutine lfs
    !
    ! This routine calculates cumulus downdraft descent to produce the
    ! vertical profiles for cumulus downdrafts (i.e. t,q,u and v and fluxes)
    ! It calculates moist descent for entraining/detraining plume by
    !   A) Moving air dry-adiabatically to next level below and
    !   B) Correcting for evaporation to obtain saturated state.
    !
    subroutine ddrafdsc
      implicit none
      real(rkx) , dimension(np) :: dmfen , xdmfde , cond , oentr , buoy
      real(rkx) , dimension(np) :: xph
      logical , dimension(np) :: llo2
      integer(ik4) :: icall , ik , is , itopde , k , n
      real(rkx) :: buo , buoyz , buoyv , xdmfdp , dz , entr , mfdqk ,  &
                   mfdsk , qdde , qeen , rain , sdde , seen , zentr , &
                   facbuo , dkbuo , dken
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
      do n = n1 , n2
        oentr(n) = d_zero
        buoy(n) = d_zero
        dmfen(n) = d_zero
        xdmfde(n) = d_zero
        dmfde(n,:) = d_zero
        mfdde_rate(n,:) = d_zero
        kined(n,:) = d_zero
      end do
      do k = 3 , nk
        is = 0
        do n = n1 , n2
          xph(n) = pf(n,k)
          llo2(n) = lddraf(n) .and. mfd(n,k-1) < d_zero
          if ( llo2(n) ) is = is + 1
        end do
        if ( is == 0 ) cycle
        do n = n1 , n2
          if ( llo2(n) ) then
            entr = entrdd*mfd(n,k-1)*(geof(n,k-1)-geof(n,k))*regrav
            dmfen(n) = entr
            xdmfde(n) = entr
          end if
        end do
        if ( k > itopde ) then
          do n = n1 , n2
            if ( llo2(n) ) then
              dmfen(n) = d_zero
              xdmfde(n) = mfd(n,itopde)*(pf(n,k)-pf(n,k-1)) / &
                (pf(n,nk+1)-pf(n,itopde))
            end if
          end do
        end if
        if ( k <= itopde ) then
          do n = n1 , n2
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
        do n = n1 , n2
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
        do n = n1 , n2
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
    end subroutine ddrafdsc
    !
    ! Updates t and q tendencies, precipitation rates, global diagnostics
    !
    subroutine dtdqc
      implicit none
      logical :: lltest
      integer(ik4) :: k , ik , n
      real(rkx) :: ximp , xalv , zp , gq , gs , gh , zs , zq
      real(rkx) , dimension(np,nk) :: xmfus , xmfuq , xmfds , xmfdq
      real(rkx) , dimension(:,:) , allocatable :: dtdt , dqdt , dp
      real(rkx) , dimension(:,:) , allocatable :: bb , r1 , r2
      logical , dimension(:,:) , allocatable :: llcumbas
      !-----------------------------
      ! 1. Setup and initializations
      ! ----------------------------
      ximp = d_one - rmfsoltq
      allocate (dtdt(np,nk))
      allocate (dqdt(np,nk))
      allocate (dp(np,nk))
      do k = 1 , nk
        do n = n1 , n2
          penth(n,k) = d_zero
        end do
      end do
      ! Zero detrained liquid water if diagnostic cloud scheme to be used
      ! This means that detrained liquid water will be evaporated in the
      ! cloud environment and not fed directly into a cloud liquid water
      ! variable
      lltest = .not. lepcld
      if ( lltest ) then
        do k = 1 , nk
          do n = n1 , n2
            lude(n,k) = d_zero
          end do
        end do
      end if
      do k = 1 , nk
        do n = n1 , n2
          if ( ldcum(n) ) then
            dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
            xmfus(n,k) = mfus(n,k)
            xmfds(n,k) = mfds(n,k)
            xmfuq(n,k) = mfuq(n,k)
            xmfdq(n,k) = mfdq(n,k)
          end if
        end do
      end do
      if ( rmfsoltq > d_zero ) then
        !-------------------------------------------
        ! 2. Recompute convective fluxes if implicit
        ! ------------------------------------------
        do k = itopm2 , nk
          ik = k - 1
          do n = n1 , n2
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
        end do
      end if
      !----------------------
      ! 3. Compute tendencies
      ! ---------------------
      do k = itopm2 , nk
        if ( k < nk ) then
          do n = n1 , n2
            if ( ldcum(n) ) then
              xalv = mlw(t(n,k))
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
          do n = n1 , n2
            if ( ldcum(n) ) then
              xalv = mlw(t(n,k))
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
        do k = itopm2 , nk
          do n = n1 , n2
            if ( ldcum(n) ) then
              tent(n,k) = tent(n,k) + dtdt(n,k)
              tenq(n,k) = tenq(n,k) + dqdt(n,k)
              penth(n,k) = dtdt(n,k)*cpd
            end if
          end do
        end do
      else
        !----------------------
        ! 3.2 Implicit solution
        ! ---------------------
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
        ! DTDT and DQDT correspond to the RHS ("constants") of the equation
        ! The solution is in R1 and R2
        allocate (bb(np,nk))
        allocate (r1(np,nk))
        allocate (r2(np,nk))
        allocate (llcumbas(np,nk))
        llcumbas(:,:) = .false.
        bb(:,:) = d_one
        xmfus(:,:) = d_zero
        ! Fill vectors A, B and RHS
        do k = itopm2 , nk
          ik = k + 1
          do n = n1 , n2
            llcumbas(n,k) = ldcum(n) .and. k >= kctop(n) - 1
            if ( llcumbas(n,k) ) then
              zp = rmfsoltq*dp(n,k)*dtcum
              xmfus(n,k) = -zp*(mfu(n,k)+mfd(n,k))
              dtdt(n,k) = dtdt(n,k)*dtcum + t(n,k)
              dqdt(n,k) = dqdt(n,k)*dtcum + q(n,k)
              if ( k < nk ) then
                bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
              else
                bb(n,k) = d_one
              end if
            end if
          end do
        end do
        call solver(kctop,llcumbas,xmfus,bb,dtdt,r1)
        call solver(kctop,llcumbas,xmfus,bb,dqdt,r2)
        ! Compute tendencies
        do k = itopm2 , nk
          do n = n1 , n2
            if ( llcumbas(n,k) ) then
              tent(n,k) = tent(n,k) + (r1(n,k)-t(n,k))/dtcum
              tenq(n,k) = tenq(n,k) + (r2(n,k)-q(n,k))/dtcum
              penth(n,k) = (r1(n,k)-t(n,k))/dtcum
            end if
          end do
        end do
        deallocate (llcumbas)
        deallocate (r2)
        deallocate (r1)
        deallocate (bb)
      end if
      deallocate (dp)
      deallocate (dqdt)
      deallocate (dtdt)
    end subroutine dtdqc
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
      real(rkx) , dimension(np,nk) , intent(in) :: mfu   ! Mass Flux Updraft
      real(rkx) , dimension(np,nk) , intent(in) :: mfd   ! Mass Flux Downdraft
      real(rkx) , dimension(np,nk) :: uen , ven , mfuu , mfdu , mfuv , mfdv
      integer(ik4) :: ik , ikb , k , n
      real(rkx) :: zp , ximp
      real(rkx) , dimension(:,:) , allocatable :: dudt , dvdt , dp
      real(rkx) , dimension(:,:) , allocatable :: bb , r1 , r2
      logical , dimension(:,:) , allocatable :: llcumbas
      ximp = d_one - rmfsoluv
      allocate (dudt(np,nk))
      allocate (dvdt(np,nk))
      allocate (dp(np,nk))
      do k = 1 , nk
        do n = n1 , n2
          if ( ldcum(n) ) then
            uen(n,k) = u(n,k)
            ven(n,k) = v(n,k)
            dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
          end if
        end do
      end do
      !--------------------------------------------------
      ! 1. Calculate fluxes and update u and v tendencies
      ! -------------------------------------------------
      do k = itopm2 , nk
        ik = k - 1
        do n = n1 , n2
          if ( ldcum(n) ) then
            mfuu(n,k) = mfu(n,k) * (uu(n,k) - ximp*uen(n,ik))
            mfuv(n,k) = mfu(n,k) * (vu(n,k) - ximp*ven(n,ik))
            mfdu(n,k) = mfd(n,k) * (ud(n,k) - ximp*uen(n,ik))
            mfdv(n,k) = mfd(n,k) * (vd(n,k) - ximp*ven(n,ik))
          end if
        end do
      end do
      ! linear fluxes below cloud
      if ( abs(rmfsoluv) < almostzero ) then
        do k = itopm2 , nk
          do n = n1 , n2
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
        end do
      end if
      ! ----------------------
      ! 1.2 Compute tendencies
      ! ----------------------
      do k = itopm2 , nk
        if ( k < nk ) then
          ik = k + 1
          do n = n1 , n2
            if ( ldcum(n) ) then
              dudt(n,k) = dp(n,k)*(mfuu(n,ik)-mfuu(n,k)+mfdu(n,ik)-mfdu(n,k))
              dvdt(n,k) = dp(n,k)*(mfuv(n,ik)-mfuv(n,k)+mfdv(n,ik)-mfdv(n,k))
            end if
          end do
        else
          do n = n1 , n2
            if ( ldcum(n) ) then
              dudt(n,k) = -dp(n,k)*(mfuu(n,k)+mfdu(n,k))
              dvdt(n,k) = -dp(n,k)*(mfuv(n,k)+mfdv(n,k))
            end if
          end do
        end if
      end do
      if ( abs(rmfsoluv) < almostzero ) then
        !----------------------
        ! 1.3 Update tendencies
        ! ---------------------
        do k = itopm2 , nk
          do n = n1 , n2
            if ( ldcum(n) ) then
              tenu(n,k) = tenu(n,k) + dudt(n,k)
              tenv(n,k) = tenv(n,k) + dvdt(n,k)
            end if
          end do
        end do
      else
        !------------------------------------------------------------------
        ! 1.4 Implicit solution
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k; reuse A and B
        ! DUDT and DVDT correspond to the RHS ("constants") of the equation
        ! The solution is in R1 and R2
        ! -----------------------------------------------------------------
        allocate (bb(np,nk))
        allocate (r1(np,nk))
        allocate (r2(np,nk))
        allocate (llcumbas(np,nk))
        llcumbas(:,:) = .false.
        bb(:,:) = d_one
        mfuu(:,:) = d_zero
        ! Fill vectors A, B and RHS
        do k = itopm2 , nk
          ik = k + 1
          do n = n1 , n2
            llcumbas(n,k) = ldcum(n) .and. k >= kctop(n) - 1
            if ( llcumbas(n,k) ) then
              zp = rmfsoluv*dp(n,k)*dtcum
              mfuu(n,k) = -zp*(mfu(n,k)+mfd(n,k))
              dudt(n,k) = dudt(n,k)*dtcum + uen(n,k)
              dvdt(n,k) = dvdt(n,k)*dtcum + ven(n,k)
              if ( k < nk ) then
                bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
              else
                bb(n,k) = d_one
              end if
            end if
          end do
        end do
        call solver(kctop,llcumbas,mfuu,bb,dudt,r1)
        call solver(kctop,llcumbas,mfuu,bb,dvdt,r2)
        do k = itopm2 , nk
          do n = n1 , n2
            if ( llcumbas(n,k) ) then
              tenu(n,k) = tenu(n,k) + (r1(n,k)-uen(n,k))/dtcum
              tenv(n,k) = tenv(n,k) + (r2(n,k)-ven(n,k))/dtcum
            end if
          end do
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
    ! This routine does the final calculation of convective fluxes in the
    ! cloud layer and in the subcloud layer
    !
    subroutine cfluxes
      implicit none
      real(rkx) , dimension(np) :: rhebc
      integer(ik4) :: ik , ikb , k , n
      integer(ik4) , dimension(np) :: idbas
      logical :: llddraf
      real(rkx) :: alfaw , cons1 , cons1a , cons2 , denom , drfl , &
                   drfl1 , fac , pdr , pds , xrfl , xrfln , rmin ,  &
                   rnew , snmlt , zp , rcpecons , rcucov
      !-------------------
      ! 0. Setup constants
      ! ------------------
      cons1a = cpd/(wlhf*egrav*rtaumel)
      cons2 = rmfcfl/(egrav*dtcum)
      !-------------------------------------
      ! 1. Determine final convective fluxes
      ! ------------------------------------
      do n = n1 , n2
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
      do k = itopm2 , nk
        ikb = min(k+1,nk)
        do n = n1 , n2
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
      mflxr(:,nk+1) = d_zero
      mflxs(:,nk+1) = d_zero
      !----------------------------------
      ! 1.1 Scale fluxes below cloud base
      !      linear dcrease
      ! ---------------------------------
      do n = n1 , n2
        if ( ldcum(n) ) then
          ikb = kcbot(n)
          ik = ikb + 1
          zp = ((pf(n,nk+1)-pf(n,ik))/(pf(n,nk+1)-pf(n,ikb)))
          if ( ktype(n) == 3 ) zp = zp*zp
          mfu(n,ik) = mfu(n,ikb)*zp
          mfus(n,ik) = (mfus(n,ikb) - mlw(tf(n,ikb)) * mful(n,ikb))*zp
          mfuq(n,ik) = (mfuq(n,ikb)+mful(n,ikb))*zp
          mful(n,ik) = d_zero
        end if
      end do
      do k = itopm2 , nk
        do n = n1 , n2
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
      do k = itopm2 , nk
        do n = n1 , n2
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
      do k = itopm2 , nk
        do n = n1 , n2
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
    end subroutine cfluxes
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
    subroutine cloudbase
      implicit none
      integer(ik4) , dimension(np) :: ictop , icbot , ibotsc , iidpl
      integer(ik4) , dimension(np,nk) :: iilab
      logical , dimension(np) :: ll_ldbase , llgo_on , lldeep , lldcum , &
                lldsc , llfirst , llresetn
      logical :: llreset
      integer(ik4) :: icall , ik , is , k , n , kk , kt1 , kt2 , kt , kb
      real(rkx) , dimension(np,nk) :: xs , suh , wu2h , buoh
      real(rkx) , dimension(np,nk+1) :: xsenh , xqenh
      real(rkx) , dimension(np) :: qold , xph
      real(rkx) , dimension(np) :: zmix
      real(rkx) , dimension(np) :: dz , cbase
      real(rkx) , dimension(np,nk) :: xlu , xqu , xtu , xuu , xvu
      ! local for CAPE at every departure level
      real(rkx) , dimension(np,nk) :: xcape
      real(rkx) :: buof , c2 , epsadd
      real(rkx) :: rho    ! Density at surface (kg/m^3)
      real(rkx) :: khvfl  ! Surface buoyancy flux (K m/s)
      real(rkx) :: ws     ! Sigma_w at lowest model level (m/s)
      real(rkx) :: qexc   ! Humidity excess at lowest model level (kg/kg)
      real(rkx) :: texc   ! Temperature excess at lowest model level (K)
      real(rkx) :: eps    ! Fractional entrainment rate [m^-1]
      real(rkx) :: tvenh  ! Environment virtual temperature (K)
      real(rkx) :: tvuh   ! Updraft virtual temperature (K)
      real(rkx) :: xlglac ! Updraft liquid water frozen in one layer
      real(rkx) :: qsu , cor , dq , alfaw , facw , faci , fac ,     &
                   esdp , dqsdt , dtdp , dp , pdifftop , pdiffbot , &
                   sf , xqf , xaw , xbw
      real(rkx) , dimension(np) :: tven1 , tvu1
      real(rkx) :: tven2 , tvu2 ! pseudoadiabatic T_v
      real(rkx) , dimension(np) :: dtvtrig ! virtual temperatures
      real(rkx) :: work1 , work2 ! work for T and w perturbations
      real(rkx) :: xtmp , entrpen
      !-----------------------------------
      ! 0. Initialize constants and fields
      ! ----------------------------------
      c2 = 0.55_rkx
      xaw = d_one
      xbw = d_one
      epsadd = 1.e-4_rkx
      do n = n1 , n2
        wubase(n) = d_zero
        llgo_on(n) = .true.
        llfirst(n) = .true.
        idpl(n) = nk
      end do
      kt1 = nk350
      kt2 = nk060
      do k = 1 , nk
        do n = n1 , n2
          xtu(n,k) = tu(n,k)
          xqu(n,k) = qu(n,k)
          xlu(n,k) = lu(n,k)
          xuu(n,k) = uu(n,k)
          xvu(n,k) = vu(n,k)
          iilab(n,k) = ilab(n,k)
          xcape(n,k) = d_zero
        end do
      end do
      !-------------------------------------------
      ! 1. Prepare fields by linear interpolation
      !    of specific humidity and static energy
      ! ------------------------------------------
      do k = 1 , nk
        do n = n1 , n2
          wu2h(n,k) = d_zero
          xs(n,k) = cpd*t(n,k) + geo(n,k)
          xqenh(n,k) = qf(n,k)
          xsenh(n,k) = cpd*tf(n,k) + geof(n,k)
        end do
      end do
      do kk = nk , kt1 , -1
        !
        ! Find first departure level that produces deepest cloud top
        ! or take surface level for shallow convection and Sc
        !
        ! ----------------------------------------------
        ! 1.2 Initialise fields at departure model level
        ! ----------------------------------------------
        !
        is = 0
        do n = n1 , n2
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
            dtvtrig(n) = d_zero
            xuu(n,kk) = u(n,kk)*(pf(n,kk+1)-pf(n,kk))
            xvu(n,kk) = v(n,kk)*(pf(n,kk+1)-pf(n,kk))
          end if
        end do
        if ( is /= 0 ) then
          if ( kk == nk ) then
            do n = n1 , n2
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
            do n = n1 , n2
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
                  do k = kk + 1 , kk - 1 , -1
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
        do k = kk - 1 , kt2 , -1
          is = 0
          if ( kk == nk ) then
            ! 1/z mixing for shallow
            do n = n1 , n2
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
            do n = n1 , n2
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
          do n = n1 , n2
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
            do n = n1 , n2
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
          do n = n1 , n2
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
          do k = nk , 1 , -1
            do n = n1 , n2
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
          do n = n1 , n2
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
            do k = nk , 1 , -1
              do n = n1 , n2
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
            end do
          end if
          do n = n1 , n2
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
      do n = n1 , n2
        cape(n) = maxval(xcape(n,:))
      end do
    end subroutine cloudbase
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
      logical , dimension(np) , intent(in) :: ldcum
      logical , dimension(np) , intent(in) :: lddraf
      real(rkx) , dimension(np,nk) , intent(in) :: mfu
      real(rkx) , dimension(np,nk) , intent(in) :: mfd
      real(rkx) , dimension(np,nk) , intent(in) :: pudrate
      real(rkx) , dimension(np,nk) , intent(in) :: pddrate
      integer(ik4) :: ik , k , n , nt
      real(rkx) :: zp , mfa , ximp , erate , posi
      real(rkx) , dimension(:,:,:) , allocatable :: cen , cu , cd , xtenc , mfc
      real(rkx) , dimension(:,:) , allocatable :: dp , bb , r1
      logical , dimension(:,:) , allocatable :: llcumask , llcumbas

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
      do k = 2 , nk
        do n = n1 , n2
          llcumask(n,k) = .false.
          if ( ldcum(n) ) then
            dp(n,k) = egrav/(pf(n,k+1)-pf(n,k))
            if ( k >= kctop(n)-1 ) llcumask(n,k) = .true.
          end if
        end do
      end do
      do nt = 1 , ntrac
        !---------------------------------------
        ! 1. Define tracers at full sigma levels
        ! --------------------------------------
        do k = 2 , nk
          ik = k - 1
          do n = n1 , n2
            cen(n,k,nt) = qtrac(n,k,nt)
            cd(n,k,nt) = qtrac(n,ik,nt)
            cu(n,k,nt) = qtrac(n,ik,nt)
            mfc(n,k,nt) = d_zero
            xtenc(n,k,nt) = d_zero
          end do
        end do
        do n = n1 , n2
          cu(n,nk,nt) = qtrac(n,nk,nt)
        end do
        !--------------------------
        ! 2. Compute updraft values
        ! -------------------------
        do k = nk - 1 , 3 , -1
          ik = k + 1
          do n = n1 , n2
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
        end do
        !----------------------------
        ! 3. Compute downdraft values
        ! ---------------------------
        do k = 3 , nk
          ik = k - 1
          do n = n1 , n2
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
        end do
        ! In order to avoid negative Tracer at nk adjust cd
        k = nk
        ik = k - 1
        do n = n1 , n2
          if ( lddraf(n) ) then
            posi = -dp(n,k) *(mfu(n,k)*cu(n,k,nt) + &
              mfd(n,k)*cd(n,k,nt)-(mfu(n,k)+mfd(n,k))*qtrac(n,ik,nt))
            if ( qtrac(n,k,nt)+posi*dtcum < d_zero ) then
              mfa = d_one/min(-cmfcmin,mfd(n,k))
              cd(n,k,nt) = ((mfu(n,k)+mfd(n,k))*qtrac(n,ik,nt) - &
                mfu(n,k)*cu(n,k,nt)+qtrac(n,k,nt) / &
                (dtcum*dp(n,k)))*mfa
            end if
          end if
        end do
      end do

      do nt = 1 , ntrac
        !------------------
        ! 4. Compute fluxes
        ! -----------------
        do k = 2 , nk
          ik = k - 1
          do n = n1 , n2
            if ( llcumask(n,k) ) then
              mfa = mfu(n,k) + mfd(n,k)
              mfc(n,k,nt) = mfu(n,k)*cu(n,k,nt) + &
                mfd(n,k)*cd(n,k,nt) - ximp*mfa*cen(n,ik,nt)
            end if
          end do
        end do
        !----------------------------
        ! 5. Compute tendencies = rhs
        ! ---------------------------
        do k = 2 , nk - 1
          ik = k + 1
          do n = n1 , n2
            if ( llcumask(n,k) ) then
              xtenc(n,k,nt) = dp(n,k)*(mfc(n,ik,nt)-mfc(n,k,nt))
            end if
          end do
        end do
        k = nk
        do n = n1 , n2
          if ( ldcum(n) ) xtenc(n,k,nt) = -dp(n,k)*mfc(n,k,nt)
        end do
      end do
      if ( abs(rmfsolct) < almostzero ) then
        !---------------------
        ! 6. Update tendencies
        ! --------------------
        do nt = 1 , ntrac
          do k = 2 , nk
            do n = n1 , n2
              if ( llcumask(n,k) ) then
                tenc(n,k,nt) = tenc(n,k,nt)+xtenc(n,k,nt)
              end if
            end do
          end do
        end do
      else
        !---------------------
        ! 7. Implicit solution
        ! ---------------------
        ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
        ! tenc corresponds to the RHS ("constants") of the equation
        ! The solution is in R1
        allocate (bb(np,nk))
        allocate (r1(np,nk))
        allocate (llcumbas(np,nk))
        llcumbas(:,:) = .false.
        bb(:,:) = d_one
        do nt = 1 , ntrac
          ! Fill vectors A, B and RHS
          do k = 2 , nk
            ik = k + 1
            do n = n1 , n2
              llcumbas(n,k) = llcumask(n,k)
              if ( llcumbas(n,k) ) then
                zp = rmfsolct*dp(n,k)*dtcum
                mfc(n,k,nt) = -zp*(mfu(n,k)+mfd(n,k))
                xtenc(n,k,nt) = xtenc(n,k,nt)*dtcum + qtrac(n,k,nt)
                ! for implicit solution including tendency source term
                if ( k < nk ) then
                  bb(n,k) = d_one + zp*(mfu(n,ik)+mfd(n,ik))
                else
                  bb(n,k) = d_one
                end if
              end if
            end do
          end do
          call solver(kctop,llcumbas,mfc(:,:,nt),bb,xtenc(:,:,nt),r1)
          ! Compute tendencies
          do k = 2 , nk
            do n = n1 , n2
              !  for implicit solution including tendency source term
              !  tenc(n,k,nt) = (r1(n,k)-qtrac(n,k,nt))/dtcum
              if ( llcumbas(n,k) ) then
                tenc(n,k,nt) = tenc(n,k,nt) + (r1(n,k)-qtrac(n,k,nt))/dtcum
              end if
            end do
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

    pure real(rkx) function xmin(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      xmin = y - d_half*(abs(x-y)-(x-y))
    end function xmin
    pure real(rkx) function lwocp(t)
      implicit none
      real(rkx) , intent(in) :: t
      real(rkx) :: gtzero
      gtzero = max(d_zero,sign(d_one,t-tzero))
      lwocp = gtzero*wlhvocp + (d_one-gtzero)*wlhsocp
    end function lwocp
    pure real(rkx) function xalpha(t)
      implicit none
      real(rkx) , intent(in) :: t
      xalpha = min(d_one,((max(rtice,min(rtwat,t))-rtice)*rtwat_rtice_r)**2)
    end function xalpha
    !
    ! Magnus Tetens formula
    !
    pure real(rkx) function fesat(t)
      implicit none
      real(rkx) , intent(in) :: t
      fesat = c2es*(xalpha(t)*exp((c3les*((t-tzero)/(t-c4les)))) + &
            (d_one-xalpha(t))*exp((c3ies*((t-tzero)/(t-c4ies)))))
    end function fesat
    pure real(rkx) function fdqsat(t)
      implicit none
      real(rkx) , intent(in) :: t
      fdqsat = xalpha(t)*c5alvcp*(d_one/(t-c4les)**2) + &
              (d_one-xalpha(t))*c5alscp*(d_one/(t-c4ies)**2)
    end function fdqsat
    pure real(rkx) function mlwocp(t)
      implicit none
      real(rkx) , intent(in) :: t
      mlwocp = xalpha(t)*wlhvocp+(d_one-xalpha(t))*wlhsocp
    end function mlwocp
    pure real(rkx) function mlw(t)
      implicit none
      real(rkx) , intent(in) :: t
      mlw = xalpha(t)*wlhv+(d_one-xalpha(t))*wlhs
    end function mlw
    pure real(rkx) function esw(t)
      implicit none
      real(rkx) , intent(in) :: t
      esw = c3les*(t-tzero)/(t-c4les)
    end function esw
    pure real(rkx) function esi(t)
      implicit none
      real(rkx) , intent(in) :: t
      esi = c3ies*(t-tzero)/(t-c4ies)
    end function esi
  end subroutine ntiedtke

end module mod_cu_tiedtke
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
