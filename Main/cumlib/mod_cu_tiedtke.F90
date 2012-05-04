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
!
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_constants
  use mod_cu_common
  use mod_cu_tables
  use mod_service
!
  private
!
  public :: allocate_mod_cu_tiedtke , tiedtkedrv
!
  public :: entrpen , entrscv , entrmid , entrdd , cmfctop , cmfcmax , &
            cmfcmin , cmfdeps , rhcdd , cprcon , iconv , nmctop ,      &
            lmfpen , lmfscv , lmfmid , lmfdd , lmfdudv , cmtcape , zdlev
!
  ! evaporation coefficient for kuo0
  real(dp) , pointer , dimension(:) :: cevapcu

  real(dp) , parameter :: centrmax = 3.0D-4

  real(dp) :: entrpen      ! entrainment rate for penetrative convection
  real(dp) :: entrscv      ! entrainment rate for shallow convection
  real(dp) :: entrmid      ! entrainment rate for midlevel convection
  real(dp) :: entrdd       ! entrainment rate for cumulus downdrafts
  real(dp) :: cmfctop      ! relat. cloud massflux at level above nonbuoyanc
  real(dp) :: cmtcape      ! CAPE adjustment timescale parameter
  real(dp) :: zdlev        ! Restrict rainfall up to this elevation
  real(dp) :: cmfcmax      ! maximum massflux value allowed for
  real(dp) :: cmfcmin      ! minimum massflux value (for safety)
  real(dp) :: cmfdeps      ! fractional massflux for downdrafts at lfs
  real(dp) :: rhcdd        ! relative saturation in downdrafts
  real(dp) :: cprcon       ! coefficients for determining conversion
                           ! from cloud water to rain
  integer :: iconv
  integer :: nmctop    !  max. level for cloud base of mid level conv.
  logical :: lmfpen    !  true if penetrative convection is switched on
  logical :: lmfscv    !  true if shallow convection is switched on
  logical :: lmfmid    !  true if midlevel convection is switched on
  logical :: lmfdd     !  true if cumulus downdraft is switched on
  logical :: lmfdudv   !  true if cumulus friction is switched on

  integer , pointer , dimension(:,:) :: ilab
  integer , pointer , dimension(:) :: ktype

  logical , pointer , dimension(:) :: ldland

  real(dp) :: ztmx
  real(dp) , pointer , dimension(:,:,:) :: pxtm1 , pxtte 

  real(dp) , pointer , dimension(:,:) :: ptm1 , pqm1 , pum1 , pvm1 , &
        pxlm1 , pxim1 , pxite , papp1 , paphp1 , pxtec , pqtec

  real(dp) , pointer , dimension(:) :: prsfc , pssfc , paprc , &
        paprs , ptopmax , xphfx

  real(dp) , pointer , dimension(:,:) :: ptte , pvom , pvol , pqte , &
        pxlte , pverv , xpg
  integer , pointer , dimension(:) :: kctop , kcbot

  integer :: nipoi

  contains
!
! This subroutines allocates space
!
  subroutine allocate_mod_cu_tiedtke
    implicit none
    nipoi = (ici2-ici1+1)*(jci2-jci1+1)
    call getmem1d(cevapcu,1,kz,'mod_cu_tiedtke:cevapcu')
    call getmem2d(ptte,1,nipoi,1,kz,'mod_cu_tiedtke:ptte')
    call getmem2d(pvom,1,nipoi,1,kz,'mod_cu_tiedtke:pvom')
    call getmem2d(pvol,1,nipoi,1,kz,'mod_cu_tiedtke:pvol')
    call getmem2d(pqte,1,nipoi,1,kz,'mod_cu_tiedtke:pqte')
    call getmem2d(pxlte,1,nipoi,1,kz,'mod_cu_tiedtke:pxlte')
    call getmem2d(pverv,1,nipoi,1,kz,'mod_cu_tiedtke:pverv')
    call getmem2d(xpg,1,nipoi,1,kz,'mod_cu_tiedtke:xpg')
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
    call getmem2d(pxtec,1,nipoi,1,kz,'mod_cu_tiedtke:pxtec')
    call getmem2d(pqtec,1,nipoi,1,kz,'mod_cu_tiedtke:pqtec')
    call getmem1d(kctop,1,nipoi,'mod_cu_tiedtke:kctop')
    call getmem1d(kcbot,1,nipoi,'mod_cu_tiedtke:kcbot')
    call getmem2d(paphp1,1,nipoi,1,kzp1,'mod_cu_tiedtke:paphp1')
    call getmem1d(prsfc,1,nipoi,'mod_cu_tiedtke:prsfc')
    call getmem1d(pssfc,1,nipoi,'mod_cu_tiedtke:pssfc')
    call getmem1d(paprc,1,nipoi,'mod_cu_tiedtke:paprc')
    call getmem1d(paprs,1,nipoi,'mod_cu_tiedtke:paprs')
    call getmem1d(ptopmax,1,nipoi,'mod_cu_tiedtke:ptopmax')
    call getmem1d(xphfx,1,nipoi,'mod_cu_tiedtke:xphfx')
    call getmem1d(ldland,1,nipoi,'mod_cu_tiedtke:ldland')

  end subroutine allocate_mod_cu_tiedtke
!
! This subroutines calls cucall
!
  subroutine tiedtkedrv(ktau)
    implicit none
!
    integer(8) , intent(in) :: ktau

!   local variables
    integer :: i , j , k , ii , kclth
    real(dp) :: akclth
    character (len=64) :: subroutine_name='tiedtkedrv'
    integer :: idindx=0
    !
    call time_begin(subroutine_name,idindx)

!   need to translate REGCM to TIEDTKE vars...

    ztmx = dtmdl

    total_precip_points = 0
    ilab(:,:) = 2
    cevapcu(:) = cevapu

!   evaporation coefficient for kuo0 
    if (lchem) then    
      do k = 1 , kz
        ii = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! tracers input profile : implicit loop on tracer
            pxtm1(ii,k,:) = chias(j,i,k,:)
            pxtte(ii,k,:) = tchiten(j,i,k,:)/sfcps(j,i)
            ii = ii + 1
          end do
        end do 
      end do
    else
      pxtm1(:,:,:) = d_zero ! tracers input profiles
      pxtte(:,:,:) = d_zero ! tracer tendencies
    end if

    ii = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
!     AMT NOTE: This is used in the switch between deep and shallow convectio
!     The simpler switch on pressure difference still used in ECMWF 
!     is commented out - possibly tests should be made to reinstate 
!     the ECMWF version of this - this array will then be obsolete 
        xphfx(ii) = qfx(j,i)
        ! Land/water flag - is correctly set?
        ldland(ii) = (lmask(j,i) == 0)
        ii = ii + 1
      end do
    end do

    ii = 1
    do k = 1 , kz
      ii = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Pascal
          papp1(ii,k) = (hlev(k)*sfcps(j,i)+ptop)*d_1000

          ptm1(ii,k)  = tas(j,i,k)  ! temperature
          pum1(ii,k)  = uas(j,i,k)  ! u (guessing!)
          pvm1(ii,k)  = vas(j,i,k)  ! v     "
          pqm1(ii,k)  = qvas(j,i,k) ! humidity
          pxlm1(ii,k) = qcas(j,i,k) ! cloud liquid water

          ptte(ii,k)  = tten(j,i,k)/sfcps(j,i)
          pvom(ii,k)  = uten(j,i,k)/sfcps(j,i)
          pvol(ii,k)  = vten(j,i,k)/sfcps(j,i)
          pqte(ii,k)  = qvten(j,i,k)/sfcps(j,i)
          pxlte(ii,k) = qcten(j,i,k)/sfcps(j,i)

          ! IS vertical velocity in Pa/s or in m/s?
          pverv(ii,k)  = d_half*(svv(j,i,k)+svv(j,i,k+1))

          pxim1(ii,k) = d_zero ! cloud ice water
          pxite(ii,k) = d_zero ! ice tend

          ! scheme diagnostic output - tendencies due to convection
          pxtec(ii,k) = d_zero ! detrained cloud water tendancy
          pqtec(ii,k) = d_zero ! detrained humidity tendancy

          xpg(ii,k) = hgt(j,i,k)*egrav  !   geopotential
          ii = ii + 1
        end do
      end do
    end do

    do k = 1 , kzp1
      ii = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! 1st guess pressure at full levels
          paphp1(ii,k) = (flev(k)*sfcps(j,i)+ptop)*d_1000
          ii = ii + 1
        end do
      end do
    end do

    ! Output variables (1d)
    prsfc(:) = d_zero ! CHECK - surface rain flux
    pssfc(:) = d_zero ! CHECK - surface snow flux      
    paprc(:) = d_zero ! total precip cumulative 
    paprs(:) = d_zero ! total snow cumulative 

    ptopmax(:) = ptop*d_1000  ! pressure top limit for convection 
    kctop(:) = 0
    kcbot(:) = 0

    call cucall(nipoi,nipoi,kz,kzp1,kzm1,ilab,ntr,pxtm1,pxtte,ptm1,     &
                pqm1,pum1,pvm1,pxlm1,pxim1,ptte,pqte,pvom,pvol,pxlte, &
                pxite,pverv,pxtec,pqtec,xphfx,papp1,paphp1,xpg,prsfc, &
                pssfc,paprc,paprs,ktype,ldland,kctop,kcbot,ptopmax)
    !
    ! postprocess some fields including precipitation fluxes
    !
    ii = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        if (ktype(ii) > 0) then
          if ( paprc(ii)+paprs(ii) > dlowval ) then
            total_precip_points = total_precip_points + 1
            ! total precip cumulative 
            rainc(j,i) = rainc(j,i) + paprc(ii)+paprs(ii)

            ! rainfall for bats
            if ( ktau == 0 .and. debug_level > 2 ) then
              lmpcpc(j,i)= lmpcpc(j,i) + (prsfc(ii)+pssfc(ii))
            else
              lmpcpc(j,i)= lmpcpc(j,i) + (prsfc(ii)+pssfc(ii))*aprdiv
            end if
          end if
          ii = ii + 1
        end if
      end do
    end do

    ! update tendencies - note that rate were ADDED in cudtdq
    !                     thus here we must reset the rates. 
    do k = 1 , kz
      ii = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if (ktype(ii) > 0) then
            tten(j,i,k)  = ptte(ii,k)  * sfcps(j,i)
            uten(j,i,k)  = pvom(ii,k)  * sfcps(j,i)
            vten(j,i,k)  = pvol(ii,k)  * sfcps(j,i)
            qvten(j,i,k) = pqte(ii,k)  * sfcps(j,i)
            qcten(j,i,k) = pxlte(ii,k) * sfcps(j,i)
            if ( lchem ) then
              tchiten(j,i,k,:) = pxtte(ii,k,:) * sfcps(j,i)
            end if
          end if
          ii = ii + 1
        end do
      end do
    end do

    ii = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        if (ktype(ii) > 0) then
          kclth = kcbot(ii) - kctop(ii) + 1
          akclth = d_one/dble(kclth)
          do k = kctop(ii) , kcbot(ii)
            rcldlwc(j,i,k) = cllwcv
            rcldfra(j,i,k) = d_one - (d_one-clfrcv)**akclth
          end do
        end if
        ii = ii + 1
      end do
    end do

    call time_end(subroutine_name,idindx)
  end subroutine tiedtkedrv
!
  subroutine cucall(kproma,kbdim,klev,klevp1,klevm1,ilab,ktrac,     &
                    pxtm1,pxtte,ptm1,pqm1,pum1,pvm1,pxlm1,pxim1,    &
                    ptte,pqte,pvom,pvol,pxlte,pxite,pverv,pxtec,    &
                    pqtec,pqhfla,papp1,paphp1,pgeo,prsfc,pssfc,     &
                    paprc,paprs,ktype,ldland,kctop,kcbot,ptopmax)
!
!
!     *CUCALL* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                     *CUMASTR* (CUMULUS PARAMETERIZATION)
!
!     M.Tiedtke      E.C.M.W.F.     12/1989
!
!**   PURPOSE.
!     --------
!
!          *CUCALL* - INTERFACE FOR *CUMASTR*:
!                     PROVIDES INPUT FOR CUMASTR
!                     RECEIVES UPDATED TENDENCIES, PRECIPITATION.
!
!**   INTERFACE.
!     ----------
!
!          *CUCALL* IS CALLED FROM *PHYSC*
!
!     EXTERNALS.
!     ----------
!
!          CUMASTR, CUMASTRT OR CUMASTRH
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
  integer , dimension(kbdim,klev) :: ilab
  integer , dimension(kbdim) :: ktype
  logical , dimension(kbdim) :: ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: papp1 , pgeo , pqm1 , pqte ,   &
         pqtec , ptm1 , ptte , pum1 , pverv , pvm1 , pvol , pvom ,  &
         pxim1 , pxite , pxlm1 , pxlte , pxtec
  real(dp) , dimension(kbdim) :: paprc , paprs , pqhfla , prsfc ,    &
                                pssfc , ptopmax
  integer , dimension(kbdim) , intent(out) :: kctop , kcbot
  real(dp) , dimension(kbdim,klev,ktrac) :: pxtm1 , pxtte
  intent (in) papp1 , pqm1 , ptm1 , pum1 , pvm1 , pxim1 , pxite ,   &
              pxlm1 , pxlte , pxtm1
  intent (inout) ptopmax
!
  integer , dimension(kbdim) :: itopec2
  integer :: ilevmin , it , jk , jl , jt
  logical , dimension(kbdim) :: locum
  real(dp) , dimension(kbdim,klev) :: zlu , zlude , zmfd , zmfu ,    &
         zqp1 , zqsat , zqu , zqude , ztp1 , ztu , ztvp1 , zup1 ,   &
         zvp1 , zxp1
  real(dp) , dimension(kbdim) :: zrain , ztopmax
  real(dp) :: zxip1 , zxlp1
  real(dp) , dimension(kbdim,klev,ktrac) :: zxtp1 , zxtu
  character (len=64) :: subroutine_name='cucall'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!  Executable statements
!
  lookupoverflow = .false.
!
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*    -----------------------------------
!
  do jk = 1 , klev
    do jl = 1 , kproma
      ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*ztmx
      zqp1(jl,jk) = max(d_zero,pqm1(jl,jk)+pqte(jl,jk)*ztmx)
      zxlp1 = pxlm1(jl,jk) + pxlte(jl,jk)*ztmx
      zxip1 = pxim1(jl,jk) + pxite(jl,jk)*ztmx
      zxp1(jl,jk) = max(d_zero,zxlp1+zxip1)
      ztvp1(jl,jk) = ztp1(jl,jk)*(d_one+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
      zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*ztmx
      zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*ztmx
      it = int(ztp1(jl,jk)*d_1000)
      if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
      zqsat(jl,jk) = min(d_half,zqsat(jl,jk))
      zqsat(jl,jk) = zqsat(jl,jk)/(d_one-vtmpc1*zqsat(jl,jk))
    end do
 
    if ( lookupoverflow ) then
      call fatal(__FILE__,__LINE__,'Cumulus Tables lookup error: OVERFLOW')
    end if
 
    do jt = 1 , ktrac
      do jl = 1 , kproma
        zxtp1(jl,jk,jt) = pxtm1(jl,jk,jt) + pxtte(jl,jk,jt)*ztmx
      end do
    end do
 
  end do
  do jl = 1 , kproma
    zrain(jl) = d_zero
    locum(jl) = .false.
  end do
!
!-----------------------------------------------------------------------
!*    2.     CALL 'CUMASTR'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION)
!*    -----------------------------------------------------------
!
  select case (iconv)
  case (1)
    call cumastr(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,    &
                 zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,pxtte,&
                 pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,pvom,pvol,&
                 prsfc,pssfc,paprc,paprs,pxtec,pqtec,zqude,locum,   &
                 ktype,kcbot,kctop,ztu,zqu,zlu,zlude,zmfu,zmfd,     &
                 zrain)
  case (2)
    call cumastrt(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                  zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                  pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                  pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,pqtec,    &
                  zqude,locum,ktype,kcbot,kctop,ztu,zqu,zlu,zlude,  &
                  zmfu,zmfd,zrain)
  case (3)
    call cumastrh(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                  zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                  pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                  pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,pqtec,    &
                  zqude,locum,ktype,kcbot,kctop,ztu,zqu,zlu,zlude,  &
                  zmfu,zmfd,zrain)
 
  case default
  end select
!
!     ------------------------------------------------------------------
!*    3.     PRESSURE ALTITUDE OF CONVECTIVE CLOUD TOPS.
!     -------- -------- -- ---------- ----- -----
!
  ilevmin = klev - 4
!
  do jl = 1 , kproma
    itopec2(jl) = klevp1
  end do
!
  do jk = 1 , ilevmin
    do jl = 1 , kproma
      if ( ilab(jl,jk) == 2 .and. itopec2(jl) == klevp1 ) itopec2(jl) = jk
    end do
  end do

!
  ztopmax(1:kproma) = ptopmax(1:kproma)
  do jl = 1 , kproma
    if ( itopec2(jl) == 1 ) then
      ptopmax(jl) = papp1(jl,1)
    else if ( itopec2(jl) /= klevp1 ) then
      ptopmax(jl) = paphp1(jl,itopec2(jl))
    else
      ptopmax(jl) = 99999.0D0
    end if
    ptopmax(jl) = min(ptopmax(jl),ztopmax(jl))
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cucall
!
!---------------------------------------------------------------------
!
  subroutine cumastr(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,pqen,&
                     pxen,puen,pven,ptven,ktrac,ldland,pxten,pxtu,  &
                     pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,ptte,pqte,&
                     pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,pqtec, &
                     pqude,ldcum,ktype,kcbot,kctop,ptu,pqu,plu,     &
                     plude,pmfu,pmfd,prain)
!
!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTR* IS CALLED FROM *CUCALL*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASC:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN SUBROUTINE CUPARAM)
!     ------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
  integer , dimension(kbdim,klev) :: ilab
  integer , dimension(kbdim) :: kcbot , kctop , ktype
  logical , dimension(kbdim) :: ldcum , ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(dp) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(dp) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ihmin , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(dp) :: zalvdcp , zalvs , zb , zbi , zcons2 , zcor , zdepth ,  &
             zdhdz , zdqmin , zdqsdt , zdz , zeps , zes , zfac ,    &
             zgam , zhhat , zhsat , zmfmax , zpbmpt , zqalv ,       &
             zqsat , zqst1 , zqumqe , zrh , zro , ztau , zzz
  real(dp) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                zhcbase , zheat , zhmin , zmfub ,   &
                                zmfub1 , zrfl , zsfl
  real(dp) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zhhatt , zmfdq , zmfds , zmful , &
         zmfuq , zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , &
         zuu , zvd , zvu , zxenh
  real(dp) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
  character (len=64) :: subroutine_name='cumastr'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!     Executable statements
 
  lookupoverflow = .false.
!
!-----------------------------------------------------------------------
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!     --------------------------------
!
  zcons2 = d_one/(egrav*dtcum)

! *AMT* NOTE!
! this paramter is the CAPE adjustment timescale which in the global model
! was a function of horizontal resolution (nn wavenumber of a spectral model)
! this is translated roughly into horizontal resolution in meters
!
!      ztau = min(3.0D0*3600.0D0,7200.0D0*63.0D0/nn)
       ztau=min(3.0D0*3600.0D0,cmtcape*ds)
!
!----------------------------------------------------------------------
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!     ---------------------------------------------------
!
  call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
             puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
             zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
             zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
             zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
             zdpmel,plu,plude,pqude,ilab)
!
!-----------------------------------------------------------------------
!*    3.0          CLOUD BASE CALCULATIONS
!     -----------------------
!
!*    (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!     ---------------------------------------
!
  call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
              paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
              kcbot,ilab)
!
!*    (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*    THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!     -----------------------------------------
!
  jk = 1
  do jl = 1 , kproma
    zdqpbl(jl) = 0.00D0
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
!*    (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*    AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*    THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!     ------------------------------------------
!
  do jl = 1 , kproma
    ikb = kcbot(jl)
    zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
    zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
    llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010D0,llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,-1.10D0*pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
  end do
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!     -------------------------------------------
!
!*    (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*    CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*    FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!     -------------------------------------------------
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
      zes = min(d_half,zes)
      lo = zes < 0.40D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
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
    ktype(jl) = merge(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.0D4)
    zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
  end do
!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__,'Cumulus Tables lookup error: OVERFLOW')
  end if
!
!     FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!     -------------------------------------------
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
!
  zb = 25.0D0
  zbi = d_one/(zb*egrav)
  do jk = klev , 1 , -1
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl) == 1 .and. ihmin(jl) == kcbot(jl)
      if ( llo1 .and. jk < kcbot(jl) .and. jk >= ictop0(jl) ) then
        zalvs = merge(wlhv,wlhs,ztenh(jl,jk) > tzero)
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zdhdz = (zcpen(jl,jk-1)*pten(jl,jk-1)-zcpen(jl,jk)          &
                *pten(jl,jk)+zalvs*(pqen(jl,jk-1)-pqen(jl,jk))      &
                +(pgeo(jl,jk-1)-pgeo(jl,jk)))                       &
                *egrav/(pgeo(jl,jk-1)-pgeo(jl,jk))
        zdepth = zgeoh(jl,jk) - zgeoh(jl,ikb)
        zfac = sqrt(d_one+zdepth*zbi)
        zhmin(jl) = zhmin(jl) + zdhdz*zfac*zdz
        zrh = -zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
        if ( zhmin(jl) > zrh ) ihmin(jl) = jk
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl) == 1 ) then
      if ( ihmin(jl) < ictop0(jl) ) ihmin(jl) = ictop0(jl)
    end if
  end do
!
!*    (B) DO ASCENT IN 'CUASC'IN ABSENCE OF DOWNDRAFTS
!     --------------------------------------------
!
  call cuasc(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven, &
             ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,   &
             zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,ilab,&
             ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,zmful,&
             plude,pqude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,zcpen,  &
             zcpcu,kcbot,kctop,ictop0)
!
!*    (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!     CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!     ---------------------------------------------------------
!
  do jl = 1 , kproma
    zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
    if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.D4 ) then
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
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!     ------------------------------
!
!
  if ( lmfdd ) then
!
!*      (A) DETERMINE LFS IN 'CUDLFS'
!       -------------------------
!
    call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
!
!*      (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!       -----------------------------------------------
!
    call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                 ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                 zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)
!
  end if
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*    CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*    AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*    ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!     -------------------------------------------
!
  do jl = 1 , kproma
    zheat(jl) = d_zero
    zcape(jl) = d_zero
    zmfub1(jl) = zmfub(jl)
  end do
!
  do jk = 1 , klev
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl) == 1
      if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zheat(jl) = zheat(jl)                                       &
                    + ((pten(jl,jk-1)-pten(jl,jk)                   &
                    +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)           &
                    +vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                    *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
        zcape(jl) = zcape(jl)                                        &
                    + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)  &
                    +egrav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))          &
                    -egrav*plu(jl,jk))*zdz
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl) == 1 ) then
      ikb = kcbot(jl)
      zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
      zmfub1(jl) = max(zmfub1(jl),0.0010D0)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub1(jl) = min(zmfub1(jl),zmfmax)
    end if
  end do
!
!*    RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*    DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*    FOR SHALLOW CONVECTION (KTYPE=2)
!     --------------------------------------------
!
  do jl = 1 , kproma
    if ( ktype(jl) == 2 ) then
      ikb = kcbot(jl)
      llo1 = pmfd(jl,ikb) < d_zero .and. loddraf(jl)
      zeps = merge(cmfdeps,d_zero,llo1)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb) &
               - (d_one-zeps)*zqenh(jl,ikb)
      zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. &
             ldcum(jl) .and. zmfub(jl) < zmfmax
      zmfub1(jl) = merge(zdqpbl(jl)/(egrav                          &
                   *max(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl)        &
                   -zmfub(jl)) < 0.20D0*zmfub(jl))
    end if
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
        pmfd(jl,jk) = pmfd(jl,jk)*zfac
        zmfds(jl,jk) = zmfds(jl,jk)*zfac
        zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
      end if
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
          zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
        end if
      end do
    end do
!
  end do
!
!*    NEW VALUES OF CLOUD BASE MASS FLUX
!     ----------------------------------
!
  do jl = 1 , kproma
    if ( ldcum(jl) ) zmfub(jl) = zmfub1(jl)
  end do


!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*    FOR PENETRATIVE CONVECTION (TYPE=1),
!*    FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*    AND FOR MID-LEVEL CONVECTION (TYPE=3).
!     --------------------------------------------------
!
  call cuasc(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven, &
             ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,   &
             zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,ilab,&
             ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,zmful,&
             plude,pqude,zdmfup,ihmin,zhhatt,zhcbase,zqsenh,zcpen,  &
             zcpcu,kcbot,kctop,ictop0)
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!     ------------------------------------------
!
  call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
             zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
             ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
             zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
             itopm2)
!
!-----------------------------------------------------------------------
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!     --------------------------------------------------
!
  call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
              pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
              zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
              zsfl,zcpen,pqtec,pqude,prsfc,pssfc,paprc,paprs)

!*AMT*
!  do jl = 1 , kproma
!    if (ktype(jl) > 0) write (50,*) zcape(jl),ktype(jl),prsfc(jl)
!  end do
!*AMT*
!
!-----------------------------------------------------------------------
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!     --------------------------------------------------
!
!
  if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                             kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                             zuu,zud,zvu,zvd,pmfu,pmfd)
!
  call time_end(subroutine_name,idindx)
  end subroutine cumastr
!
  subroutine cumastrh(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,    &
                      pqen,pxen,puen,pven,ptven,ktrac,ldland,pxten, &
                      pxtu,pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,    &
                      ptte,pqte,pvom,pvol,prsfc,pssfc,paprc,paprs,  &
                      pxtec,pqtec,pqude,ldcum,ktype,kcbot,kctop,ptu,&
                      pqu,plu,plude,pmfu,pmfd,prain)
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTRH* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 3
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASCT:  CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
!     ----------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
  integer , dimension(kbdim,klev) :: ilab
  integer , dimension(kbdim) :: kcbot , kctop , ktype
  logical , dimension(kbdim) :: ldcum , ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(dp) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(dp) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(dp) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
             zdz , zeps , zes , zfac , zgam , zhhat , zhsat ,       &
             zmfmax , zpbmpt , zqalv , zqsat , zqst1 , zqumqe ,     &
             zro , ztau , zzz
  real(dp) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                zhcbase , zheat , zmfub , zmfub1 ,  &
                                zrfl , zsfl
  real(dp) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
         zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
         zvd , zvu , zxenh
  real(dp) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
  character (len=64) :: subroutine_name='cumastrh'
  integer :: idindx=0
  !
  call time_begin(subroutine_name,idindx)
!
!     Executable statements
 
  lookupoverflow = .false.
!-----------------------------------------------------------------------
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!     --------------------------------
!
  zcons2 = d_one/(egrav*dtcum)
! *AMT* NOTE!
! this paramter is the CAPE adjustment timescale which in the global model
! was a function of horizontal resolution (nn wavenumber of a spectral model)
! this is translated roughly into horizontal resolution in meters
!      ztau = min(3.0D0*3600.0D0,7200.0D0*63.0D0/nn)
  ztau=min(3.0D0*3600.0D0,22.7D0*ds)
!
!-----------------------------------------------------------------------
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!     ---------------------------------------------------
!
  call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
             puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
             zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
             zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
             zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
             zdpmel,plu,plude,pqude,ilab)
!
!-----------------------------------------------------------------------
!*    3.0          CLOUD BASE CALCULATIONS
!     -----------------------
!
!*    (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!     ---------------------------------------
!
  call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
              paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
              kcbot,ilab)
!
!*    (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*    THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!     -----------------------------------------
!
  jk = 1
  do jl = 1 , kproma
    zdqpbl(jl) = 0.00D0
    zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
    idtop(jl) = 0
  end do
  do jk = 2 , klev
    do jl = 1 , kproma
      zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)                           &
                  *(paphp1(jl,jk+1)-paphp1(jl,jk))
      if ( jk >= kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)    &
           *(paphp1(jl,jk+1)-paphp1(jl,jk))
    end do
  end do
!
!*    (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*    AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*    THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!     ------------------------------------------
!
  do jl = 1 , kproma
    ikb = kcbot(jl)
    zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
    zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
    llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010D0,llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,-1.10D0*pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
  end do
!
!-----------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!     -------------------------------------------
!
!*    (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*    CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*    FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!     -------------------------------------------------
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
      zes = min(d_half,zes)
      lo = zes < 0.40D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
      zhhat = zhsat - (zzz+zgam*zzz)/(d_one+zgam*zzz*zqalv)         &
              *max(zqsenh(jl,jk)-zqenh(jl,jk),d_zero)
      if ( jk < ictop0(jl) .and. zhcbase(jl) > zhhat ) ictop0(jl) = jk
    end do
  end do
!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__,'Cumulus Tables lookup error: OVERFLOW')
  end if
!!
!!    DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!    (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!     DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.D4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl) == 1)
!     ENDDO
!!
!*    (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!     ----------------------------------------------
!
  call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
              ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
              zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
              ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
              zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
              ictop0)
!
!*    (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!     CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!     ---------------------------------------------------------
!
  do jl = 1 , kproma
    zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
    if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.D4 ) then
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
!-----------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!     ------------------------------
!
  if ( lmfdd ) then
!
!*      (A) DETERMINE LFS IN 'CUDLFS'
!       -------------------------
!
    call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
!
!*      (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!       -----------------------------------------------
!
    call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                 ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                 zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)
!
  end if
!
!*    5.5          RECALCULATE CLOUD BASE MASSFLUX FROM A
!*    CAPE CLOSURE FOR DEEP CONVECTION (KTYPE=1)
!*    AND BY PBL EQUILIBRUM TAKING DOWNDRAFTS INTO
!*    ACCOUNT FOR SHALLOW CONVECTION (KTYPE=2)
!     -------------------------------------------
!
  do jl = 1 , kproma
    zheat(jl) = d_zero
    zcape(jl) = d_zero
    zmfub1(jl) = zmfub(jl)
  end do
!
  do jk = 1 , klev
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl) == 1
      if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zheat(jl) = zheat(jl)                                       &
                    + ((pten(jl,jk-1)-pten(jl,jk)                   &
                    +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)           &
                    +vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                    *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
        zcape(jl) = zcape(jl)                                       &
                    + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) &
                    +egrav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))         &
                    -egrav*plu(jl,jk))*zdz
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl) == 1 ) then
      ikb = kcbot(jl)
      zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
      zmfub1(jl) = max(zmfub1(jl),0.0010D0)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      zmfub1(jl) = min(zmfub1(jl),zmfmax)
    end if
  end do
!
!*    RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*    DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!*    FOR SHALLOW CONVECTION (KTYPE=2)
!     --------------------------------------------
!
  do jl = 1 , kproma
    if ( ktype(jl) == 2 ) then
      ikb = kcbot(jl)
      llo1 = pmfd(jl,ikb) < d_zero .and. loddraf(jl)
      zeps = merge(cmfdeps,d_zero,llo1)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)         &
               - (d_one-zeps)*zqenh(jl,ikb)
      zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. &
             ldcum(jl) .and. zmfub(jl) < zmfmax
      zmfub1(jl) = merge(zdqpbl(jl)/(egrav                          &
             *max(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl)        &
                   -zmfub(jl)) < 0.20D0*zmfub(jl))
    end if
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
        pmfd(jl,jk) = pmfd(jl,jk)*zfac
        zmfds(jl,jk) = zmfds(jl,jk)*zfac
        zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
        zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
      end if
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
          zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
        end if
      end do
    end do
!
  end do
!
!*    NEW VALUES OF CLOUD BASE MASS FLUX
!     ----------------------------------
!
  do jl = 1 , kproma
    if ( ldcum(jl) ) zmfub(jl) = zmfub1(jl)
  end do
!
!-----------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*    FOR PENETRATIVE CONVECTION (TYPE=1),
!*    FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*    AND FOR MID-LEVEL CONVECTION (TYPE=3).
!     --------------------------------------------------
!
  call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
              ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
              zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
              ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
              zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
              ictop0)
!
!-----------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!     ------------------------------------------
!
  call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
             zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
             ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
             zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
             itopm2)
!
!-----------------------------------------------------------------------
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!     --------------------------------------------------
!
  call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
              pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
              zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
              zsfl,zcpen,pqtec,pqude,prsfc,pssfc,paprc,paprs)
!
!-----------------------------------------------------------------------
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!     --------------------------------------------------
!
  if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                             kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                             zuu,zud,zvu,zvd,pmfu,pmfd)
!
  call time_end(subroutine_name,idindx)
  end subroutine cumastrh
!
  subroutine cumastrt(kproma,kbdim,klev,klevp1,klevm1,ilab,pten,    &
                      pqen,pxen,puen,pven,ptven,ktrac,ldland,pxten, &
                      pxtu,pxtte,pverv,pqsen,pqhfla,paphp1,pgeo,    &
                      ptte,pqte,pvom,pvol,prsfc,pssfc,paprc,paprs,  &
                      pxtec,pqtec,pqude,ldcum,ktype,kcbot,kctop,ptu,&
                      pqu,plu,plude,pmfu,pmfd,prain)
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U AND V DUE TO CONVECTIVE PROCESSES.
!     PROCESSES CONSIDERED ARE: CONVECTIVE FLUXES, FORMATION OF
!     PRECIPITATION, EVAPORATION OF FALLING RAIN BELOW CLOUD BASE,
!     SATURATED CUMULUS DOWNDRAFTS.
!
!**   INTERFACE.
!     ----------
!
!          *CUMASTRT* IS CALLED FROM *CUCALL* IN CASE ICONV .EQ. 2
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE
!     T,Q,U,V,PHI AND P AND MOISTURE TENDENCIES.
!     IT RETURNS ITS OUTPUT TO THE SAME SPACE
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES
!      2.RATES OF CONVECTIVE PRECIPITATION
!        (USED IN SUBROUTINE SURF)
!
!     METHOD
!     -------
!
!     PARAMETERIZATION IS DONE USING A MASSFLUX-SCHEME.
!        (1) DEFINE CONSTANTS AND PARAMETERS
!        (2) SPECIFY VALUES (T,Q,QS...) AT HALF LEVELS AND
!            INITIALIZE UPDRAFT- AND DOWNDRAFT-VALUES IN 'CUINI'
!        (3) CALCULATE CLOUD BASE IN 'CUBASE'
!            AND SPECIFY CLOUD BASE MASSFLUX FROM PBL MOISTURE BUDGET
!        (4) DO CLOUD ASCENT IN 'CUASC' IN ABSENCE OF DOWNDRAFTS
!        (5) DO DOWNDRAFT CALCULATIONS:
!              (A) DETERMINE VALUES AT LFS IN 'CUDLFS'
!              (B) DETERMINE MOIST DESCENT IN 'CUDDRAF'
!              (C) RECALCULATE CLOUD BASE MASSFLUX CONSIDERING THE
!                  EFFECT OF CU-DOWNDRAFTS
!        (6) DO FINAL CLOUD ASCENT IN 'CUASC'
!        (7) DO FINAL ADJUSMENTS TO CONVECTIVE FLUXES IN 'CUFLX',
!            DO EVAPORATION IN SUBCLOUD LAYER
!        (8) CALCULATE INCREMENTS OF T AND Q IN 'CUDTDQ'
!        (9) CALCULATE INCREMENTS OF U AND V IN 'CUDUDV'
!
!     EXTERNALS.
!     ----------
!
!       CUINI:  INITIALIZES VALUES AT VERTICAL GRID USED IN CU-PARAMETR.
!       CUBASE: CLOUD BASE CALCULATION FOR PENETR.AND SHALLOW CONVECTION
!       CUASCT: CLOUD ASCENT FOR ENTRAINING PLUME
!       CUDLFS: DETERMINES VALUES AT LFS FOR DOWNDRAFTS
!       CUDDRAF:DOES MOIST DESCENT FOR CUMULUS DOWNDRAFTS
!       CUFLX:  FINAL ADJUSTMENTS TO CONVECTIVE FLUXES (ALSO IN PBL)
!       CUDQDT: UPDATES TENDENCIES FOR T AND Q
!       CUDUDV: UPDATES TENDENCIES FOR U AND V
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.T.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.T.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.T.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.T.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.T.  CUMULUS FRICTION SWITCHED ON
!
!     MODEL PARAMETERS (DEFINED IN MODULE mo_cumulus_flux)
!     ----------------------------------------------------
!     ENTRPEN    ENTRAINMENT RATE FOR PENETRATIVE CONVECTION
!     ENTRSCV    ENTRAINMENT RATE FOR SHALLOW CONVECTION
!     ENTRMID    ENTRAINMENT RATE FOR MIDLEVEL CONVECTION
!     ENTRDD     ENTRAINMENT RATE FOR CUMULUS DOWNDRAFTS
!     CMFCTOP    RELATIVE CLOUD MASSFLUX AT LEVEL ABOVE NONBUOYANCY LEVEL
!     CMFCMAX    MAXIMUM MASSFLUX VALUE ALLOWED FOR
!     CMFCMIN    MINIMUM MASSFLUX VALUE (FOR SAFETY)
!     CMFDEPS    FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS
!     CPRCON     COEFFICIENT FOR CONVERSION FROM CLOUD WATER TO RAIN
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
  integer , dimension(kbdim,klev) :: ilab
  integer , dimension(kbdim) :: kcbot , kctop , ktype
  logical , dimension(kbdim) :: ldcum , ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(dp) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(dp) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(dp) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
             zeps , zes , zfac , zgam , zhhat , zhsat , zmfmax ,    &
             zpbmpt , zqalv , zqsat , zqst1 , zqumqe , zzz
  real(dp) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
         zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
         zvd , zvu , zxenh
  real(dp) , dimension(kbdim) :: zdqcv , zdqpbl , zentr , zhcbase ,  &
                                zmfub , zmfub1 , zrfl , zsfl
  real(dp) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
  character (len=64) :: subroutine_name='cumastrt'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!     Executable statements
!
  lookupoverflow = .false.
!
!---------------------------------------------------------------------
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!     --------------------------------
!
!
  zcons2 = d_one/(egrav*dtcum)
!
!---------------------------------------------------------------------
!*    2.           INITIALIZE VALUES AT VERTICAL GRID POINTS IN 'CUINI'
!     ---------------------------------------------------
!
  call cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen,pxen,  &
             puen,pven,ptven,ktrac,pxten,zxtenh,pxtu,zxtd,zmfuxt,   &
             zmfdxt,pverv,pgeo,paphp1,zgeoh,ztenh,zqenh,zqsenh,     &
             zxenh,ilwmin,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,pmfu,pmfd,&
             zmfus,zmfds,zmfuq,zmfdq,zdmfup,zdmfdp,zcpen,zcpcu,     &
             zdpmel,plu,plude,pqude,ilab)
!
!---------------------------------------------------------------------
!*    3.0          CLOUD BASE CALCULATIONS
!     -----------------------
!
!
!*    (A) DETERMINE CLOUD BASE VALUES IN 'CUBASE'
!     ---------------------------------------
!
  call cubase(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,zgeoh,    &
              paphp1,ptu,pqu,plu,puen,pven,zuu,zvu,zcpcu,ldcum,     &
              kcbot,ilab)
!
!*    (B) DETERMINE TOTAL MOISTURE CONVERGENCE AND
!*    THEN DECIDE ON TYPE OF CUMULUS CONVECTION
!     -----------------------------------------
!
  jk = 1
  do jl = 1 , kproma
    zdqpbl(jl) = 0.00D0
    zdqcv(jl) = pqte(jl,jk)*(paphp1(jl,jk+1)-paphp1(jl,jk))
    idtop(jl) = 0
  end do
  do jk = 2 , klev
    do jl = 1 , kproma
      zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)                           &
                  *(paphp1(jl,jk+1)-paphp1(jl,jk))
      if ( jk >= kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)    &
           *(paphp1(jl,jk+1)-paphp1(jl,jk))
    end do
  end do
!
!*    (C) DETERMINE MOISTURE SUPPLY FOR BOUNDARY LAYER
!*    AND DETERMINE CLOUD BASE MASSFLUX IGNORING
!*    THE EFFECTS OF DOWNDRAFTS AT THIS STAGE
!     ------------------------------------------
!
  do jl = 1 , kproma
    ikb = kcbot(jl)
    zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
    zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
    llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010D0,llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl) > max(d_zero,-1.10D0 *pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl) == 1)
  end do
!
!---------------------------------------------------------------------
!*    4.0          DETERMINE CLOUD ASCENT FOR ENTRAINING PLUME
!     -------------------------------------------
!
!*    (A) ESTIMATE CLOUD HEIGHT FOR ENTRAINMENT/DETRAINMENT
!*    CALCULATIONS IN CUASC (MAX.POSSIBLE CLOUD HEIGHT
!*    FOR NON-ENTRAINING PLUME, FOLLOWING A.-S.,1974)
!     -------------------------------------------------
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
      zes = min(d_half,zes)
      lo = zes < 0.40D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
      zhhat = zhsat - (zzz+zgam*zzz)/(d_one+zgam*zzz*zqalv)         &
              *max(zqsenh(jl,jk)-zqenh(jl,jk),d_zero)
      if ( jk < ictop0(jl) .and. zhcbase(jl) > zhhat ) ictop0(jl) = jk
    end do
  end do
!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__,'Cumulus Tables lookup error: OVERFLOW')
  end if
!!
!!    DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!    (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!     DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.D4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl) == 1)
!     ENDDO
!!
!*    (B) DO ASCENT IN 'CUASCT' IN ABSENCE OF DOWNDRAFTS
!     ----------------------------------------------
!
  call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
              ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
              zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
              ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
              zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
              ictop0)
!
!*    (C) CHECK CLOUD DEPTH AND CHANGE ENTRAINMENT RATE ACCORDINGLY
!     CALCULATE PRECIPITATION RATE (FOR DOWNDRAFT CALCULATION)
!     ---------------------------------------------------------
!
  do jl = 1 , kproma
    zpbmpt = paphp1(jl,kcbot(jl)) - paphp1(jl,kctop(jl))
    if ( ldcum(jl) .and. ktype(jl) == 1 .and. zpbmpt < 2.D4 ) then
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
!---------------------------------------------------------------------
!*    5.0          CUMULUS DOWNDRAFT CALCULATIONS
!     ------------------------------
!
!
  if ( lmfdd ) then
!
!*      (A) DETERMINE LFS IN 'CUDLFS'
!       -------------------------
!
    call cudlfs(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,     &
                ktrac,zxtenh,pxtu,zxtd,zmfdxt,zgeoh,paphp1,ptu,pqu, &
                zuu,zvu,ldcum,kcbot,kctop,zmfub,zrfl,ztd,zqd,zud,   &
                zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,idtop,loddraf)
!
!*      (B)  DETERMINE DOWNDRAFT T,Q AND FLUXES IN 'CUDDRAF'
!       -----------------------------------------------
!
    call cuddraf(kproma,kbdim,klev,klevp1,ztenh,zqenh,puen,pven,    &
                 ktrac,zxtenh,zxtd,zmfdxt,zgeoh,paphp1,zrfl,ztd,zqd,&
                 zud,zvd,pmfd,zmfds,zmfdq,zdmfdp,zcpcu,loddraf)
!
!*      (C)  RECALCULATE CONVECTIVE FLUXES DUE TO EFFECT OF
!*      DOWNDRAFTS ON BOUNDARY LAYER MOISTURE BUDGET
!       --------------------------------------------
!
    do jl = 1 , kproma
      if ( loddraf(jl) ) then
        ikb = kcbot(jl)
        llo1 = pmfd(jl,ikb) < d_zero
        zeps = merge(cmfdeps,d_zero,llo1)
        zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)       &
                 - (d_one-zeps)*zqenh(jl,ikb)
        zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1 = zdqpbl(jl) > d_zero .and. zqumqe > zdqmin .and. ldcum(jl) &
               .and. zmfub(jl) < zmfmax
        zmfub1(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),zmfub(jl),llo1)
        zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),(ktype(jl) == 1 .or.  &
                     ktype(jl) == 2) .and. &
                     abs(zmfub1(jl)-zmfub(jl)) < 0.20D0*zmfub(jl))
      end if
    end do
    do jk = 1 , klev
      do jl = 1 , kproma
        if ( loddraf(jl) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
          pmfd(jl,jk) = pmfd(jl,jk)*zfac
          zmfds(jl,jk) = zmfds(jl,jk)*zfac
          zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
        end if
      end do
!
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( loddraf(jl) ) then
            zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
            zmfdxt(jl,jk,jt) = zmfdxt(jl,jk,jt)*zfac
          end if
        end do
      end do
!
    end do
!
!*      NEW VALUES OF CLOUD BASE MASS FLUX
!       ----------------------------------
    do jl = 1 , kproma
      if ( loddraf(jl) ) zmfub(jl) = zmfub1(jl)
    end do
!
  end if
!
!---------------------------------------------------------------------
!*    6.0          DETERMINE FINAL CLOUD ASCENT FOR ENTRAINING PLUME
!*    FOR PENETRATIVE CONVECTION (TYPE=1),
!*    FOR SHALLOW TO MEDIUM CONVECTION (TYPE=2)
!*    AND FOR MID-LEVEL CONVECTION (TYPE=3).
!     ----------------------------------------------------
!
!
  call cuasct(kproma,kbdim,klev,klevp1,klevm1,ztenh,zqenh,puen,pven,&
              ktrac,zxtenh,pxten,pxtu,zmfuxt,pten,pqen,pqsen,pgeo,  &
              zgeoh,paphp1,pqte,pverv,ilwmin,ldcum,ldland,ktype,    &
              ilab,ptu,pqu,plu,zuu,zvu,pmfu,zmfub,zentr,zmfus,zmfuq,&
              zmful,plude,pqude,zdmfup,zcpen,zcpcu,kcbot,kctop,     &
              ictop0)
!
!---------------------------------------------------------------------
!*    7.0          DETERMINE FINAL CONVECTIVE FLUXES IN 'CUFLX'
!     ------------------------------------------
!
  call cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ztenh,zqenh,ktrac, &
             zxtenh,zmfuxt,zmfdxt,paphp1,zgeoh,kcbot,kctop,idtop,   &
             ktype,loddraf,ldcum,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq, &
             zmful,zdmfup,zdmfdp,zrfl,prain,zcpcu,pten,zsfl,zdpmel, &
             itopm2)
!
!
!---------------------------------------------------------------------
!*    8.0          UPDATE TENDENCIES FOR T AND Q IN SUBROUTINE CUDTDQ
!     --------------------------------------------------
!
  call cudtdq(kproma,kbdim,klev,klevp1,itopm2,ldcum,ktrac,paphp1,   &
              pten,ptte,pqte,pxtte,pxtec,zmfuxt,zmfdxt,zmfus,zmfds, &
              zmfuq,zmfdq,zmful,zdmfup,zdmfdp,plude,zdpmel,zrfl,    &
              zsfl,zcpen,pqtec,pqude,prsfc,pssfc,paprc,paprs)
!
!---------------------------------------------------------------------
!*    9.0          UPDATE TENDENCIES FOR U AND U IN SUBROUTINE CUDUDV
!     --------------------------------------------------
!
  if ( lmfdudv ) call cududv(kproma,kbdim,klev,klevp1,itopm2,ktype, &
                             kcbot,paphp1,ldcum,puen,pven,pvom,pvol,&
                             zuu,zud,zvu,zvd,pmfu,pmfd)
!
  call time_end(subroutine_name,idindx)
  end subroutine cumastrt
!
  subroutine cuini(kproma,kbdim,klev,klevp1,klevm1,pten,pqen,pqsen, &
                   pxen,puen,pven,ptven,ktrac,pxten,pxtenh,pxtu,    &
                   pxtd,pmfuxt,pmfdxt,pverv,pgeo,paphp1,pgeoh,ptenh,&
                   pqenh,pqsenh,pxenh,klwmin,ptu,pqu,ptd,pqd,puu,   &
                   pvu,pud,pvd,pmfu,pmfd,pmfus,pmfds,pmfuq,pmfdq,   &
                   pdmfup,pdmfdp,pcpen,pcpcu,pdpmel,plu,plude,pqude,&
                   klab)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , kproma , ktrac , klevm1 , klevp1
  integer , dimension(kbdim,klev) :: klab
  integer , dimension(kbdim) :: klwmin
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfdp ,       &
         pdmfup , pdpmel , pgeo , pgeoh , plu , plude , pmfd ,      &
         pmfdq , pmfds , pmfu , pmfuq , pmfus , pqd , pqen , pqenh ,&
         pqsen , pqsenh , pqu , pqude , ptd , pten , ptenh , ptu ,  &
         ptven , pud , puen , puu , pvd , pven , pverv , pvu ,      &
         pxen , pxenh
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtd , &
         pxten , pxtenh , pxtu
  intent (in) paphp1 , pgeo , pqen , pqsen , pten , ptven , puen , &
              pven , pverv , pxen , pxten
  intent (out) klab , klwmin , pdmfdp , pdmfup , pdpmel , plu ,     &
               plude , pmfd , pmfdq , pmfds , pmfdxt , pmfu ,       &
               pmfuq , pmfus , pmfuxt , pqd , pqu , pqude , ptd ,   &
               ptu , pud , puu , pvd , pvu , pxenh , pxtd , pxtu
  intent (inout) pcpcu , pcpen , pgeoh , pqenh , pqsenh , ptenh ,   &
                 pxtenh
!
  integer :: icall , ik , jk , jl , jt
  logical , dimension(kbdim) :: loflag
  real(dp) :: zarg , zcpm , zzs
  real(dp) , dimension(kbdim) :: zph , zwmax
  character (len=64) :: subroutine_name='cumini'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*    ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*    FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!     ----------------------------------------------
!
  do jk = 1 , klev
    do jl = 1 , kproma
!         pcpen(jl,jk) = cpd*(1.+vtmpc2*max(pqen(jl,jk),0.00D0))
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
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pxtenh(jl,jk,jt) = (pxten(jl,jk,jt)+pxten(jl,jk-1,jt))*d_half
      end do
    end do
!
!
    ik = jk
    icall = 0
    call cuadjtq(kproma,kbdim,klev,ik,zph,ptenh,pqsenh,loflag,icall)
!
    do jl = 1 , kproma
      pxenh(jl,jk) = (pxen(jl,jk)+pxen(jl,jk-1))*d_half
      pqenh(jl,jk) = min(pqen(jl,jk-1),pqsen(jl,jk-1))              &
                     + (pqsenh(jl,jk)-pqsen(jl,jk-1))
      pqenh(jl,jk) = max(pqenh(jl,jk),d_zero)
!         pcpcu(jl,jk) = cpd*(1.+vtmpc2*pqenh(jl,jk))
      pcpcu(jl,jk) = cpd
    end do
  end do
!
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
!
  do jt = 1 , ktrac
    do jl = 1 , kproma
      pxtenh(jl,klev,jt) = pxten(jl,klev,jt)
      pxtenh(jl,1,jt) = pxten(jl,1,jt)
    end do
  end do
!
  do jk = klevm1 , 2 , -1
    do jl = 1 , kproma
      zzs = max(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),pcpcu(jl,jk+1)* &
             ptenh(jl,jk+1)+pgeoh(jl,jk+1))
      ptenh(jl,jk) = (zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
    end do
  end do
!
  do jk = klev , 3 , -1
    do jl = 1 , kproma
      if ( pverv(jl,jk) < zwmax(jl) ) then
        zwmax(jl) = pverv(jl,jk)
        klwmin(jl) = jk
      end if
    end do
  end do
!
!-----------------------------------------------------------------------
!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*    ---------------------------------------------
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
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pxtu(jl,jk,jt) = pxtenh(jl,jk,jt)
        pxtd(jl,jk,jt) = pxtenh(jl,jk,jt)
        pmfuxt(jl,jk,jt) = d_zero
        pmfdxt(jl,jk,jt) = d_zero
      end do
    end do
!
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cuini
!
  subroutine cuasc(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,puen,&
                   pven,ktrac,pxtenh,pxten,pxtu,pmfuxt,pten,pqen,   &
                   pqsen,pgeo,pgeoh,paphp1,pqte,pverv,klwmin,ldcum, &
                   ldland,ktype,klab,ptu,pqu,plu,puu,pvu,pmfu,pmfub,&
                   pentr,pmfus,pmfuq,pmful,plude,pqude,pdmfup,khmin,&
                   phhatt,phcbase,pqsenh,pcpen,pcpcu,kcbot,kctop,   &
                   kctop0)
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac , &
             klevm1
  integer , dimension(kbdim) :: kcbot , kctop , kctop0 , khmin ,    &
                                klwmin , ktype
  integer , dimension(kbdim,klev) :: klab
  logical , dimension(kbdim) :: ldcum , ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
         pgeoh , phhatt , plu , plude , pmfu , pmful , pmfuq ,      &
         pmfus , pqen , pqenh , pqsen , pqsenh , pqte , pqu ,       &
         pqude , pten , ptenh , ptu , puen , puu , pven , pverv ,   &
         pvu
  real(dp) , dimension(kbdim) :: pentr , phcbase , pmfub
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
         pxtu
  intent (in) ldland , pcpcu , phcbase , phhatt , pqsenh , pxtenh
  intent (inout) kcbot , kctop , kctop0 , klab , ktype , ldcum ,    &
                 plu , plude , pmfu , pmfub , pmful , pmfuq ,       &
                 pmfus , pmfuxt , pqu , pqude , ptu , puu , pvu ,   &
                 pxtu
!
  integer :: icall , ik , ikb , ikt , jk , jl , jt
  logical , dimension(kbdim) :: loflag
  real(dp) :: zalvs , zbuo , zbuoyz , zcons2 , zdmfdu ,      &
             zdmfeu , zdnoprc , zdprho , zdrodz , zdt , zdz , zfac ,&
             zga , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,     &
             zmfusk , zmfuxtk , zmse , znevn , zodmax , zprcon ,    &
             zqcod , zqeen , zqude , zscde , zscod , zseen ,        &
             ztglace , zxteen , zxtude , zz , zzdmf
  real(dp) , dimension(kbdim) :: zbuoy , zdmfde , zdmfen , zmfuu ,   &
                                zmfuv , zpbase , zph , zqold
  real(dp) , dimension(kbdim,klev) :: zodetr , zoentr
!
  character (len=64) :: subroutine_name='cuasc'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!----------------------------------------------------------------------
!*    1.           SPECIFY PARAMETERS
!     ------------------
!
  zcons2 = d_one/(egrav*dtcum)
  ztglace = tzero - 13.0D0
  zqold(1:kproma) = d_zero

! AMT NOTE!!! in the original scheme, this level which restricts rainfall 
! below a certain pressure (from the surface) is hard wired according to the 
! vertical resolution of the model
!      if ( klev /= 11 ) then
!        zdlev = 3.0D4
!      else if ( nn == 21 ) then
!        zdlev = 1.5D4
!      else if ( nn == 31 ) then
!        zdlev = 2.0D4 
!      else
!        zdlev = 3.0D4
!      end if
!
!----------------------------------------------------------------------
!     2.           SET DEFAULT VALUES
!     ------------------
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
      if ( .not.ldcum(jl) .and. paphp1(jl,jk) < 4.D4 ) kctop0(jl) = jk
      if ( jk < kcbot(jl) ) klab(jl,jk) = 0
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pmfuxt(jl,jk,jt) = d_zero
      end do
    end do
!
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      zoentr(jl,jk) = d_zero
      zodetr(jl,jk) = d_zero
    end do
  end do
!
!----------------------------------------------------------------------
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!     ----------------------------------
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
!
  do jt = 1 , ktrac
    do jl = 1 , kproma
      if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = d_zero
      pmfuxt(jl,klev,jt) = pmfub(jl)*pxtu(jl,klev,jt)
    end do
  end do
!
  do jl = 1 , kproma
    ldcum(jl) = .false.
  end do
!
!----------------------------------------------------------------------
!     3.5          FIND ORGANIZED ENTRAINMENT AT CLOUD BASE
!     ----------------------------------------
!
  do jl = 1 , kproma
    if ( ktype(jl) == 1 ) then
      ikb = kcbot(jl)
      zbuoy(jl) = egrav*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb)     &
                  + egrav*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
      if ( zbuoy(jl) > d_zero ) then
        zdz = (pgeo(jl,ikb-1)-pgeo(jl,ikb))*regrav
        zdrodz = -log(pten(jl,ikb-1)/pten(jl,ikb))                  &
                 /zdz - egrav/(rgas*ptenh(jl,ikb)                   &
                 *(d_one+vtmpc1*pqenh(jl,ikb)))
!           nb zoentr is here a fractional value
        zoentr(jl,ikb-1) = zbuoy(jl)*d_half/(d_one+zbuoy(jl)*zdz)   &
                           + zdrodz
        zoentr(jl,ikb-1) = min(zoentr(jl,ikb-1),centrmax)
        zoentr(jl,ikb-1) = max(zoentr(jl,ikb-1),d_zero)
      end if
    end if
  end do
!
!
!----------------------------------------------------------------------
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!     BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!     BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!     THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!     -------------------------------------------------
!
  do jk = klevm1 , 2 , -1
!
!       SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!       IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!       ----------------------------------------------------
!
    ik = jk
    if ( lmfmid .and. ik < klevm1 .and. ik > nmctop )                   &
         call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,    &
         puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
         ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
         pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)
!
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
!       RESET PMFUB IF NECESSARY
!
    do jl = 1 , kproma
      if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        pmfub(jl) = min(pmfub(jl),zmfmax)
      end if
    end do
!
!
!*      SPECIFY TURBULENT ENTRAINMENT AND DETRAINMENTS
!       RATES PLUS ORGANIZED DETRAINMENT RATES IN *CUENTR*
!       -------------------------------------
!
    ik = jk
    call cuentr(kproma,kbdim,klev,klevp1,ik,ptenh,pqenh,pqte,paphp1,&
                klwmin,ldcum,ktype,kcbot,kctop0,zpbase,pmfu,pentr,  &
                zodetr,khmin,pgeoh,zdmfen,zdmfde)
!
!       DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!       THE CLOUD ENSEMBLE ENTRAINS ENVIRONMENTAL VALUES
!       IN TURBULENT DETRAINMENT CLOUD ENSEMBLE VALUES
!       ARE DETRAINED
!       IN ORGANIZED DETRAINMENT THE DRY STATIC ENERGY AND
!       MOISTURE THAT ARE NEUTRAL COMPARED TO THE
!       ENVIRONMENTAL AIR ARE DETRAINED
!       ---------------------------------------------------
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
        zdmfde(jl) = min(zdmfde(jl),0.750D0*pmfu(jl,jk+1))
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
!             limit organized detrainment to not allowing for too
!             deep clouds
          zalvs = merge(wlhv,wlhs,ptu(jl,jk+1) > tzero)
          zmse = pcpcu(jl,jk+1)*ptu(jl,jk+1) + zalvs*pqu(jl,jk+1)   &
                 + pgeoh(jl,jk+1)
          ikt = kctop0(jl)
          znevn = (pgeoh(jl,ikt)-pgeoh(jl,jk+1))                    &
                  *(zmse-phhatt(jl,jk+1))*regrav
          if ( znevn <= d_zero ) znevn = 1.
          zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
          zodmax = ((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
          zodmax = max(zodmax,d_zero)
          zodetr(jl,jk) = min(zodetr(jl,jk),zodmax)
        end if
        zodetr(jl,jk) = min(zodetr(jl,jk),0.750D0*pmfu(jl,jk))
        pmfu(jl,jk) = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
        zqeen = pqenh(jl,jk+1)*zdmfen(jl)
        zqeen = zqeen + pqenh(jl,jk+1)*zoentr(jl,jk)
        zseen = (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
        zseen = zseen +                                             &
                (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zoentr(jl,jk)
        zscde = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
!           find moist static energy that give nonbuoyant air
        zalvs = merge(wlhv,wlhs,ptenh(jl,jk+1) > tzero)
        zga = zalvs*pqsenh(jl,jk+1)/(rwat*(ptenh(jl,jk+1)**2))
        zdt = (plu(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))&
              /(d_one/ptenh(jl,jk+1)+vtmpc1*zga)
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
        ptu(jl,jk) = max(100.0D0,ptu(jl,jk))
        ptu(jl,jk) = min(400.0D0,ptu(jl,jk))
        zqold(jl) = pqu(jl,jk)
      end if
    end do
!
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
!       DO CORRECTIONS FOR MOIST ASCENT
!       BY ADJUSTING T,Q AND L IN *CUADJTQ*
!       -----------------------------------
!
    ik = jk
    icall = 1
    call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)
!
    do jl = 1 , kproma
      if ( loflag(jl) ) then
        if ( pqu(jl,jk) < zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+vtmpc1*pqu(jl,jk)-plu(jl,jk))    &
                 - ptenh(jl,jk)*(d_one+vtmpc1*pqenh(jl,jk))
          if ( klab(jl,jk+1) == 1 ) zbuo = zbuo + d_half
          if ( zbuo > d_zero .and. pmfu(jl,jk) >= 0.010D0*pmfub(jl) .and.&
               jk >= kctop0(jl) ) then
            kctop(jl) = jk
            ldcum(jl) = .true.
            zdnoprc = merge(zdlev,1.5D4,ldland(jl))
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
!
    if ( lmfdudv ) then
      do jl = 1 , kproma
        zdmfen(jl) = zdmfen(jl) + zoentr(jl,jk)
        zdmfde(jl) = zdmfde(jl) + zodetr(jl,jk)
      end do
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) then
            zz = merge(3.0D0,2.0D0,zdmfen(jl) == d_zero)
          else
            zz = merge(d_one,d_zero,zdmfen(jl) == d_zero)
          end if
          zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
          zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
          zdmfdu = min(zdmfdu,0.750D0*pmfu(jl,jk+1))
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
!       COMPUTE ORGANIZED ENTRAINMENT
!       FOR USE AT NEXT LEVEL
!       ------------------------------
!
    do jl = 1 , kproma
      if ( loflag(jl) .and. ktype(jl) == 1 ) then
        zbuoyz = egrav*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk)         &
                 + egrav*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))             &
                 - egrav*plu(jl,jk)
        zbuoyz = max(zbuoyz,0.00D0)
        zdz = (pgeo(jl,jk-1)-pgeo(jl,jk))*regrav
        zdrodz = -log(pten(jl,jk-1)/pten(jl,jk))                    &
                 /zdz - egrav/(rgas*ptenh(jl,jk)                    &
                 *(d_one+vtmpc1*pqenh(jl,jk)))
        zbuoy(jl) = zbuoy(jl) + zbuoyz*zdz
        zoentr(jl,jk-1) = zbuoyz*d_half/(d_one+zbuoy(jl)) + zdrodz
        zoentr(jl,jk-1) = min(zoentr(jl,jk-1),centrmax)
        zoentr(jl,jk-1) = max(zoentr(jl,jk-1),d_zero)
!
      end if
    end do
!
  end do
!
!----------------------------------------------------------------------
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!     ----------------------------------------------------
!
!     (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!     AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!     FROM PREVIOUS CALCULATIONS ABOVE)
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
!
  if ( lmfdudv ) then
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        jk = kctop(jl) - 1
        puu(jl,jk) = puu(jl,jk+1)
        pvu(jl,jk) = pvu(jl,jk+1)
      end if
    end do
  end if
!
  call time_end(subroutine_name,idindx)
  end subroutine cuasc
!
  subroutine cuasct(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,    &
                    puen,pven,ktrac,pxtenh,pxten,pxtu,pmfuxt,pten,  &
                    pqen,pqsen,pgeo,pgeoh,paphp1,pqte,pverv,klwmin, &
                    ldcum,ldland,ktype,klab,ptu,pqu,plu,puu,pvu,    &
                    pmfu,pmfub,pentr,pmfus,pmfuq,pmful,plude,pqude, &
                    pdmfup,pcpen,pcpcu,kcbot,kctop,kctop0)
!
!          THIS ROUTINE DOES THE CALCULATIONS FOR CLOUD ASCENTS
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD ASCENTS FOR CU-PARAMETRIZATION
!          (VERTICAL PROFILES OF T,Q,L,U AND V AND CORRESPONDING
!           FLUXES AS WELL AS PRECIPITATION RATES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTRT*.
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          AND THEN CALCULATE MOIST ASCENT FOR
!          ENTRAINING/DETRAINING PLUME.
!          ENTRAINMENT AND DETRAINMENT RATES DIFFER FOR
!          SHALLOW AND DEEP CUMULUS CONVECTION.
!          IN CASE THERE IS NO PENETRATIVE OR SHALLOW CONVECTION
!          CHECK FOR POSSIBILITY OF MID LEVEL CONVECTION
!          (CLOUD BASE VALUES CALCULATED IN *CUBASMC*)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!          *CUENTR*  CALCULATE ENTRAINMENT/DETRAINMENT RATES
!          *CUBASMC* CALCULATE CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma , ktrac
  integer , dimension(kbdim) :: kcbot , kctop , kctop0 , klwmin ,   &
                                ktype
  integer , dimension(kbdim,klev) :: klab
  logical , dimension(kbdim) :: ldcum , ldland
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
         pgeoh , plu , plude , pmfu , pmful , pmfuq , pmfus , pqen ,&
         pqenh , pqsen , pqte , pqu , pqude , pten , ptenh , ptu ,  &
         puen , puu , pven , pverv , pvu
  real(dp) , dimension(kbdim) :: pentr , pmfub
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
         pxtu
  intent (in) ldland , pcpcu , pxtenh
  intent (out) pqude
  intent (inout) kcbot , kctop , klab , ktype , ldcum , plu ,       &
                 plude , pmfu , pmfub , pmful , pmfuq , pmfus ,     &
                 pmfuxt , pqu , ptu , puu , pvu , pxtu
!
  integer :: icall , ik , jk , jl , jt
  logical , dimension(kbdim) :: loflag
  real(dp) :: zbuo , zcons2 , zdmfdu , zdmfeu , zdnoprc ,    &
             zfac , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,    &
             zmfusk , zmfuxtk , zprcon , zqeen , zqude , zscde ,    &
             zseen , ztglace , zxteen , zxtude , zz , zzdmf
  real(dp) , dimension(kbdim) :: zdmfde , zdmfen , zmfuu , zmfuv ,   &
                                zpbase , zph , zqold
  character (len=64) :: subroutine_name='cuasct'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.           SPECIFY PARAMETERS
!     ------------------
!
  zcons2 = d_one/(egrav*dtcum)
  ztglace = tzero - 13.0D0

! AMT NOTE!!! in the original scheme, this level which restricts rainfall 
! below a certain pressure (from the surface) is hard wired according to the 
! vertical resolution of the model
!      if ( klev /= 11 ) then
!        zdlev = 3.0D4
!      else if ( nn == 21 ) then
!        zdlev = 1.5D4
!      else if ( nn == 31 ) then
!        zdlev = 2.0D4 
!      else
!        zdlev = 3.0D4
!      end if
!
!----------------------------------------------------------------------
!     2.           SET DEFAULT VALUES
!     ------------------
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
      if ( .not.ldcum(jl) .and. paphp1(jl,jk) < 4.D4 ) kctop0(jl) = jk
      if ( jk < kcbot(jl) ) klab(jl,jk) = 0
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pmfuxt(jl,jk,jt) = d_zero
      end do
    end do
!
  end do
!
!----------------------------------------------------------------------
!     3.0          INITIALIZE VALUES AT LIFTING LEVEL
!     ----------------------------------
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
!
  do jt = 1 , ktrac
    do jl = 1 , kproma
      if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = d_zero
      pmfuxt(jl,klev,jt) = pmfub(jl)*pxtu(jl,klev,jt)
    end do
  end do
!
  do jl = 1 , kproma
    ldcum(jl) = .false.
  end do
!
!
!----------------------------------------------------------------------
!     4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
!     BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
!     BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!     THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
!     -------------------------------------------------
!
  do jk = klevm1 , 2 , -1
!
!       SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
!       IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
!       ----------------------------------------------------
!
    ik = jk
    if ( lmfmid .and. ik < klevm1 .and. ik > nmctop )                   &
         call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,    &
         puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
         ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
         pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)
!
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
!       RESET PMFUB IF NECESSARY
!
    do jl = 1 , kproma
      if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        pmfub(jl) = min(pmfub(jl),zmfmax)
      end if
    end do
!
!
!*      SPECIFY ENTRAINMENT RATES IN *CUENTRT*
!       --------------------------------------
!
    ik = jk
    call cuentrt(kproma,kbdim,klev,klevp1,ik,ptenh,pqenh,pqte,      &
                 paphp1,klwmin,ldcum,ktype,kcbot,kctop0,zpbase,pmfu,&
                 pentr,zdmfen,zdmfde)
!
!       DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
!       ---------------------------------------------------
!
    do jl = 1 , kproma
      if ( loflag(jl) ) then
        if ( jk < kcbot(jl) ) then
          zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2)
          zdmfen(jl) = max(zdmfen(jl)-max(zmftest-zmfmax,d_zero),d_zero)
        end if
        zdmfde(jl) = min(zdmfde(jl),0.750D0*pmfu(jl,jk+1))
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
        ptu(jl,jk) = max(100.0D0,ptu(jl,jk))
        ptu(jl,jk) = min(400.0D0,ptu(jl,jk))
        zqold(jl) = pqu(jl,jk)
      end if
    end do
!
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
!       DO CORRECTIONS FOR MOIST ASCENT
!       BY ADJUSTING T,Q AND L IN *CUADJTQ*
!       -----------------------------------
!
    ik = jk
    icall = 1
    call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)
!
    do jl = 1 , kproma
      if ( loflag(jl) .and. pqu(jl,jk) < zqold(jl) ) then
        klab(jl,jk) = 2
        plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
        zbuo = ptu(jl,jk)*(d_one+vtmpc1*pqu(jl,jk)-plu(jl,jk))      &
               - ptenh(jl,jk)*(d_one+vtmpc1*pqenh(jl,jk))
        if ( klab(jl,jk+1) == 1 ) zbuo = zbuo + d_half
        if ( zbuo > d_zero .and. pmfu(jl,jk) >= 0.10D0*pmfub(jl) ) then
          kctop(jl) = jk
          ldcum(jl) = .true.
          zdnoprc = merge(zdlev,1.5D4,ldland(jl))
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
!
    if ( lmfdudv ) then
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) then
            zz = merge(3.0D0,2.0D0,zdmfen(jl) == d_zero)
          else
            zz = merge(d_one,d_zero,zdmfen(jl) == d_zero)
          end if
          zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
          zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
          zdmfdu = min(zdmfdu,0.750D0*pmfu(jl,jk+1))
          zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk) - zdmfdu*puu(jl,jk+1)
          zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk) - zdmfdu*pvu(jl,jk+1)
          if ( pmfu(jl,jk) > d_zero ) then
            puu(jl,jk) = zmfuu(jl)*(d_one/pmfu(jl,jk))
            pvu(jl,jk) = zmfuv(jl)*(d_one/pmfu(jl,jk))
          end if
        end if
      end do
    end if
!
  end do
!
!
!----------------------------------------------------------------------
!     5.           DETERMINE CONVECTIVE FLUXES ABOVE NON-BUOYANCY LEVEL
!     ----------------------------------------------------
!     (NOTE: CLOUD VARIABLES LIKE T,Q AND L ARE NOT
!     AFFECTED BY DETRAINMENT AND ARE ALREADY KNOWN
!     FROM PREVIOUS CALCULATIONS ABOVE)
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
!
  if ( lmfdudv ) then
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        jk = kctop(jl) - 1
        puu(jl,jk) = puu(jl,jk+1)
        pvu(jl,jk) = pvu(jl,jk+1)
      end if
    end do
  end if
!
  call time_end(subroutine_name,idindx)
  end subroutine cuasct
!
  subroutine cubase(kproma,kbdim,klev,klevp1,klevm1,ptenh,pqenh,    &
                    pgeoh,paph,ptu,pqu,plu,puen,pven,puu,pvu,pcpcu, &
                    ldcum,kcbot,klab)
!
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES (T AND Q)
!          FOR CUMULUS PARAMETERIZATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE VALUES FOR CU-PARAMETRIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD BASE VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CONDENSATION LEVEL
!
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD BASE
!          (NON ENTRAINING PLUME,I.E.CONSTANT MASSFLUX)
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevm1 , klevp1 , kproma
  integer , dimension(kbdim) :: kcbot
  integer , dimension(kbdim,klev) :: klab
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klevp1) :: paph
  real(dp) , dimension(kbdim,klev) :: pcpcu , pgeoh , plu , pqenh ,  &
         pqu , ptenh , ptu , puen , puu , pven , pvu
  intent (in) paph , pcpcu , pgeoh , pqenh , ptenh , puen , pven
  intent (inout) kcbot , klab , ldcum , plu , pqu , ptu , puu , pvu
!
! Local variables
!
  integer :: icall , ik , ikb , is , jk , jl
  logical , dimension(kbdim) :: loflag
  real(dp) :: zbuo , zz
  real(dp) , dimension(kbdim) :: zph , zqold
  character (len=64) :: subroutine_name='cubase'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!     1.           INITIALIZE VALUES AT LIFTING LEVEL
!     ----------------------------------
!
  do jl = 1 , kproma
    klab(jl,klev) = 1
    kcbot(jl) = klevm1
    ldcum(jl) = .false.
    puu(jl,klev) = puen(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
    pvu(jl,klev) = pven(jl,klev)*(paph(jl,klevp1)-paph(jl,klev))
  end do
!
!
!----------------------------------------------------------------------
!     2.0          DO ASCENT IN SUBCLOUD LAYER,
!     CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!     ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!     CHECK FOR BUOYANCY AND SET FLAGS
!     -------------------------------------
!
  do jk = klevm1 , 2 , -1
    is = 0
    do jl = 1 , kproma
      is = is + merge(1,0,klab(jl,jk+1) == 1)
      loflag(jl) = klab(jl,jk+1) == 1
      zph(jl) = paph(jl,jk)
    end do
    if ( is /= 0 ) then
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          pqu(jl,jk) = pqu(jl,jk+1)
          ptu(jl,jk) = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)  &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+vtmpc1*pqu(jl,jk)) - ptenh(jl,jk)&
                 *(d_one+vtmpc1*pqenh(jl,jk)) + d_half
          if ( zbuo > d_zero ) klab(jl,jk) = 1
          zqold(jl) = pqu(jl,jk)
        end if
      end do
!
      ik = jk
      icall = 1
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)
!
      do jl = 1 , kproma
        if ( loflag(jl) .and. pqu(jl,jk) < zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(d_one+vtmpc1*pqu(jl,jk)-plu(jl,jk))    &
                 - ptenh(jl,jk)*(d_one+vtmpc1*pqenh(jl,jk)) + d_half
          if ( zbuo > 0. ) then
            kcbot(jl) = jk
            ldcum(jl) = .true.
          end if
        end if
      end do
!
!         CALCULATE AVERAGES OF U AND V FOR SUBCLOUD ARA,.
!         THE VALUES WILL BE USED TO DEFINE CLOUD BASE VALUES.
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
    end if
!
  end do
!
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
!
  call time_end(subroutine_name,idindx)
  end subroutine cubase
!
  subroutine cubasmc(kproma,kbdim,klev,kk,klab,pten,pqen,pqsen,puen,&
                     pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh, &
                     ldcum,ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,&
                     puu,pvu,pmfus,pmfuq,pmful,pdmfup,pmfuu,pcpen,  &
                     pmfuv)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
  implicit none
!
  integer , intent(in) :: kbdim , kk , klev , kproma , ktrac
  integer , dimension(kbdim) :: kcbot , ktype
  integer , dimension(kbdim,klev) :: klab
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klev) :: pcpen , pdmfup , pgeo , pgeoh ,&
         plu , pmfu , pmful , pmfuq , pmfus , pqen , pqsen , pqu ,  &
         pten , ptu , puen , puu , pven , pverv , pvu
  real(dp) , dimension(kbdim) :: pentr , pmfub , pmfuu , pmfuv
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtu
  intent (in) ldcum , pcpen , pgeo , pgeoh , pqen , pqsen , pten , &
              puen , pven , pverv , pxten
  intent (out) kcbot , ktype , pdmfup , pentr , plu , pmfu , pmful ,&
               pmfuq , pmfus , pmfuu , pmfuv , pmfuxt
  intent (inout) klab , pmfub , pqu , ptu , puu , pvu , pxtu
!
  integer :: jl , jt
  logical , dimension(kbdim) :: llo3
  real(dp) :: zzzmb
  character (len=64) :: subroutine_name='cubasmc'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!     -------------------------------------------
!
  do jl = 1 , kproma
    llo3(jl) = .false.
    if ( .not.ldcum(jl) .and. klab(jl,kk+1) == 0 .and. &
          pqen(jl,kk) > 0.900D0*pqsen(jl,kk) ) then
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
!
  call time_end(subroutine_name,idindx)
  end subroutine cubasmc
!
  subroutine cuddraf(kproma,kbdim,klev,klevp1,ptenh,pqenh,puen,pven,&
                     ktrac,pxtenh,pxtd,pmfdxt,pgeoh,paphp1,prfl,ptd,&
                     pqd,pud,pvd,pmfd,pmfds,pmfdq,pdmfdp,pcpcu,     &
                     lddraf)
!
!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!
!          METHOD.
!          --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT
!
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
  logical , dimension(kbdim) :: lddraf
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
         pmfdq , pmfds , pqd , pqenh , ptd , ptenh , pud , puen ,   &
         pvd , pven
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh
  real(dp) , dimension(kbdim) :: prfl
  intent (in) lddraf , paphp1 , pcpcu , pgeoh , pqenh , ptenh , &
              puen , pven , pxtenh
  intent (out) pdmfdp
  intent (inout) pmfd , pmfdq , pmfds , pmfdxt , pqd , prfl , ptd , &
                 pud , pvd , pxtd
!
  integer :: icall , ik , is , itopde , jk , jl , jt
  logical :: llo1
  logical , dimension(kbdim) :: llo2
  real(dp) :: zbuo , zdmfdp , zentr , zmfdqk , zmfdsk , zmfduk ,     &
             zmfdvk , zmfdxtk , zqdde , zqeen , zsdde , zseen ,     &
             zxtdde , zxteen
  real(dp) , dimension(kbdim) :: zcond , zdmfde , zdmfen , zph
  character (len=64) :: subroutine_name='cuddraf'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!     1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
!     (A) CALCULATING ENTRAINMENT RATES, ASSUMING
!     LINEAR DECREASE OF MASSFLUX IN PBL
!     (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
!     AND MOISTENING IS CALCULATED IN *CUADJTQ*
!     (C) CHECKING FOR NEGATIVE BUOYANCY AND
!     SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
!     -------------------------------------------------
!
  do jk = 3 , klev
    is = 0
    do jl = 1 , kproma
      zph(jl) = paphp1(jl,jk)
      llo2(jl) = lddraf(jl) .and. pmfd(jl,jk-1) < d_zero
      is = is + merge(1,0,llo2(jl))
    end do
    if ( is /= 0 ) then
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
!
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
          ptd(jl,jk) = min(400.0D0,ptd(jl,jk))
          ptd(jl,jk) = max(100.0D0,ptd(jl,jk))
          zcond(jl) = pqd(jl,jk)
        end if
      end do
!
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
!
      ik = jk
      icall = 2
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptd,pqd,llo2,icall)
!
      do jl = 1 , kproma
        if ( llo2(jl) ) then
          zcond(jl) = zcond(jl) - pqd(jl,jk)
          zbuo = ptd(jl,jk)*(d_one+vtmpc1*pqd(jl,jk)) - &
                 ptenh(jl,jk)*(d_one+vtmpc1*pqenh(jl,jk))
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
!
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( llo2(jl) ) pmfdxt(jl,jk,jt) = pxtd(jl,jk,jt)*pmfd(jl,jk)
        end do
      end do
!
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
    end if
!
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cuddraf
!
  subroutine cudlfs(kproma,kbdim,klev,klevp1,ptenh,pqenh,puen,pven, &
                    ktrac,pxtenh,pxtu,pxtd,pmfdxt,pgeoh,paphp1,ptu, &
                    pqu,puu,pvu,ldcum,kcbot,kctop,pmfub,prfl,ptd,   &
                    pqd,pud,pvd,pmfd,pmfds,pmfdq,pdmfdp,pcpcu,kdtop,&
                    lddraf)
!
!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!
!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!
!          METHOD.
!          --------
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
  integer , dimension(kbdim) :: kcbot , kctop , kdtop
  logical , dimension(kbdim) :: ldcum , lddraf
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
         pmfdq , pmfds , pqd , pqenh , pqu , ptd , ptenh , ptu ,    &
         pud , puen , puu , pvd , pven , pvu
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh , &
         pxtu
  real(dp) , dimension(kbdim) :: pmfub , prfl
  intent (in) kcbot , kctop , ldcum , paphp1 , pcpcu , pgeoh , &
              pmfub , pqenh , pqu , ptenh , ptu , puen , puu , &
              pven , pvu , pxtenh , pxtu
  intent (out) kdtop , pmfdq , pmfds , pmfdxt , pud , pvd
  intent (inout) lddraf , pdmfdp , pmfd , pqd , prfl , ptd , pxtd
!
  integer :: icall , ik , is , jk , jl , jt , ke
  logical , dimension(kbdim) :: llo2 , llo3
  real(dp) :: zbuo , zmftop , zqtest , zttest
  real(dp) , dimension(kbdim) :: zcond , zph
  real(dp) , dimension(kbdim,klev) :: zqenwb , ztenwb
  character (len=64) :: subroutine_name='cudlf'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!     1.           SET DEFAULT VALUES FOR DOWNDRAFTS
!     ---------------------------------
!
  do jl = 1 , kproma
    lddraf(jl) = .false.
    kdtop(jl) = klevp1
  end do
!
  if ( lmfdd ) then
!
!----------------------------------------------------------------------
!       2.           DETERMINE LEVEL OF FREE SINKING BY
!       DOING A SCAN FROM TOP TO BASE OF CUMULUS CLOUDS
!
!       FOR EVERY POINT AND PROCEED AS FOLLOWS:
!
!       (1) DETEMINE WET BULB ENVIRONMENTAL T AND Q
!       (2) DO MIXING WITH CUMULUS CLOUD AIR
!       (3) CHECK FOR NEGATIVE BUOYANCY
!
!       THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
!       OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
!       TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
!       EVAPORATION OF RAIN AND CLOUD WATER)
!       ----------------------------------------------------
!
!
    ke = klev - 3
    do jk = 3 , ke
!
!
!         2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
!         FOR ENVIRONMENTAL AIR IN *CUADJTQ*
!         -------------------------------------------
!
      is = 0
      do jl = 1 , kproma
        ztenwb(jl,jk) = ptenh(jl,jk)
        zqenwb(jl,jk) = pqenh(jl,jk)
        zph(jl) = paphp1(jl,jk)
        llo2(jl) = ldcum(jl) .and. prfl(jl) > d_zero .and.            &
                   .not.lddraf(jl) .and.                            &
                   (jk < kcbot(jl) .and. jk > kctop(jl))
        is = is + merge(1,0,llo2(jl))
      end do
      if ( is /= 0 ) then
!
        ik = jk
        icall = 2
        call cuadjtq(kproma,kbdim,klev,ik,zph,ztenwb,zqenwb,llo2,icall)
!
!
!           2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
!           AND CHECK FOR NEGATIVE BUOYANCY.
!           THEN SET VALUES FOR DOWNDRAFT AT LFS.
!           ----------------------------------------
!
        do jl = 1 , kproma
          llo3(jl) = .false.
          if ( llo2(jl) ) then
            zttest = d_half*(ptu(jl,jk)+ztenwb(jl,jk))
            zqtest = d_half*(pqu(jl,jk)+zqenwb(jl,jk))
            zbuo = zttest*(d_one+vtmpc1*zqtest) - ptenh(jl,jk)      &
                   *(d_one+vtmpc1*pqenh(jl,jk))
            zcond(jl) = pqenh(jl,jk) - zqenwb(jl,jk)
            zmftop = -cmfdeps*pmfub(jl)
            if ( zbuo < d_zero .and. prfl(jl) > 10.0D0*zmftop*zcond(jl) ) then
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
!
        do jt = 1 , ktrac
          do jl = 1 , kproma
            if ( llo3(jl) ) then
              pxtd(jl,jk,jt) = d_half*(pxtu(jl,jk,jt)+pxtenh(jl,jk,jt))
              pmfdxt(jl,jk,jt) = pmfd(jl,jk)*pxtd(jl,jk,jt)
            end if
          end do
        end do
!
        if ( lmfdudv ) then
          do jl = 1 , kproma
            if ( pmfd(jl,jk) < d_zero ) then
              pud(jl,jk) = d_half*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk) = d_half*(pvu(jl,jk)+pven(jl,jk-1))
            end if
          end do
        end if
      end if
!
    end do
  end if
!
  call time_end(subroutine_name,idindx)
  end subroutine cudlfs
!
  subroutine cudtdq(kproma,kbdim,klev,klevp1,ktopm2,ldcum,ktrac,    &
                    paphp1,pten,ptte,pqte,pxtte,pxtec,pmfuxt,pmfdxt,&
                    pmfus,pmfds,pmfuq,pmfdq,pmful,pdmfup,pdmfdp,    &
                    plude,pdpmel,prfl,psfl,pcpen,pqtec,pqude,prsfc, &
                    pssfc,paprc,paprs)
!
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktopm2 , ktrac
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: paprc , paprs , prfl , prsfc ,      &
                                psfl , pssfc
  real(dp) , dimension(kbdim,klev) :: pcpen , pdmfdp , pdmfup ,      &
         pdpmel , plude , pmfdq , pmfds , pmful , pmfuq , pmfus ,   &
         pqte , pqtec , pqude , pten , ptte , pxtec
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtte
  intent (in) ldcum , paphp1 , pcpen , pdmfdp , pdmfup , pdpmel ,   &
              plude , pmfdq , pmfds , pmfdxt , pmful , pmfuq ,      &
              pmfus , pmfuxt , pqude , prfl , psfl , pten
  intent (out) pqtec , prsfc , pssfc , pxtec
  intent (inout) paprc , paprs , pqte , ptte , pxtte
!
  integer :: jk , jl , jt
  logical :: llo1
  real(dp) :: zalv , zdqdt , zdtdt , zdxtdt , zrcpm
  real(dp) , dimension(kbdim) :: zmelt , zsheat
  character (len=64) :: subroutine_name='cudtdq'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!     ------------------------------------
!
  do jl = 1 , kproma
    zmelt(jl) = d_zero
    zsheat(jl) = d_zero
  end do
!
  do jk = ktopm2 , klev
!
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
!
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
!
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
!
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( ldcum(jl) ) then
            zdxtdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))     &
                     *(pmfuxt(jl,jk,jt)+pmfdxt(jl,jk,jt))
            pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + zdxtdt
          end if
        end do
      end do
!
    end if
!
  end do
!
!
!---------------------------------------------------------------------
!     3.          UPDATE SURFACE FIELDS
!     ---------------------
!
  do jl = 1 , kproma
    prsfc(jl) = prfl(jl)
    pssfc(jl) = psfl(jl)
    paprc(jl) = paprc(jl) + ztmx*(prfl(jl)+psfl(jl))
    paprs(jl) = paprs(jl) + ztmx*psfl(jl)
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cudtdq
!
  subroutine cududv(kproma,kbdim,klev,klevp1,ktopm2,ktype,kcbot,    &
                    paphp1,ldcum,puen,pven,pvom,pvol,puu,pud,pvu,   &
                    pvd,pmfu,pmfd)
!
!
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!**   INTERFACE.
!     ----------
!
!          *CUDUDV* IS CALLED FROM *CUMASTR*
!
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktopm2
  integer , dimension(kbdim) :: kcbot , ktype
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pmfd , pmfu , pud , puen ,     &
         puu , pvd , pven , pvol , pvom , pvu
  intent (in) kcbot , ktype , ldcum , paphp1 , pmfd , pmfu , pud , &
              puen , puu , pvd , pven , pvu
  intent (inout) pvol , pvom
!
  integer :: ik , ikb , jk , jl
  real(dp) :: zdudt , zdvdt , zzp
  real(dp) , dimension(kbdim,klev) :: zmfdu , zmfdv , zmfuu , zmfuv
  character (len=64) :: subroutine_name='cududv'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!     ----------------------------------------------
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
        zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))                    &
              /(paphp1(jl,klevp1)-paphp1(jl,ikb)))
        zzp = merge(zzp**2,zzp,ktype(jl) == 3)
        zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
        zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
        zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
        zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
      end if
    end do
  end do
!
  do jk = ktopm2 , klev
!
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
!
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
!
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cududv
!
  subroutine cuentr(kproma,kbdim,klev,klevp1,kk,ptenh,pqenh,pqte,   &
                    paphp1,klwmin,ldcum,ktype,kcbot,kctop0,ppbase,  &
                    pmfu,pentr,podetr,khmin,pgeoh,pdmfen,pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
  implicit none
!
  integer , intent(in) :: kbdim , kk , klev , klevp1 , kproma
  integer , dimension(kbdim) :: kcbot , kctop0 , khmin , klwmin ,   &
                                ktype
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
  real(dp) , dimension(kbdim,klev) :: pgeoh , pmfu , podetr , pqenh ,&
         pqte , ptenh
  intent (in) kcbot , kctop0 , khmin , klwmin , ktype , ldcum , &
              paphp1 , pentr , pgeoh , pqenh , pqte , ptenh
  intent (out) pdmfde , pdmfen , podetr
  intent (inout) pmfu , ppbase
!
! Local variables
!
  integer :: ikb , ikh , iklwmin , ikt , jl
  logical :: llo1 , llo2
  real(dp) :: zarg , zdprho , zentest , zentr , zorgde , zpmid ,     &
             zrrho , ztmzk , zzmzk
  character (len=64) :: subroutine_name='cuentr'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!     -------------------------------------------
!
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!     --------------------------------------------
!
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!     -----------------------------------------
!
  pdmfde(:) = d_zero
  pdmfen(:) = d_zero
  podetr(:,:) = d_zero

  

  if (.false.) then !AMT fudge timing test

  do jl = 1 , kproma

    ppbase(jl) = paphp1(jl,kcbot(jl))
    if ( .not. ldcum(jl) ) cycle
    zrrho = (rgas*ptenh(jl,kk+1)*(d_one+vtmpc1*pqenh(jl,kk+1)))/paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
    zpmid = d_half*(ppbase(jl)+paphp1(jl,kctop0(jl)))
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho

    llo1 = kk < kcbot(jl) .and. ldcum(jl)
    pdmfde(jl) = merge(zentr,d_zero,llo1)
    llo2 = llo1 .and. ktype(jl) == 2 .and.                            &
           (ppbase(jl)-paphp1(jl,kk) < 0.2D5 .or. paphp1(jl,kk) > zpmid)
    pdmfen(jl) = merge(zentr,d_zero,llo2)
    iklwmin = max(klwmin(jl),kctop0(jl)+2)
    llo2 = llo1 .and. ktype(jl) == 3 .and. kk >= iklwmin
    if ( llo2 ) pdmfen(jl) = zentr

    if ( llo2 .and. pqenh(jl,kk+1) > 1.D-5 ) then
      pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
      zentest = max(pqte(jl,kk),d_zero)/pqenh(jl,kk+1)
      zentest = min(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
      pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
    end if

    llo2 = llo1 .and. ktype(jl) == 1 .and.                            &
           (kk >= iklwmin .or. paphp1(jl,kk) > zpmid)
    if ( llo2 ) pdmfen(jl) = zentr
!
!       organized detrainment, detrainment starts at khmin
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
        podetr(jl,kk) = min(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho
      end if
    end if
  end do

  else
  
  do jl = 1 , kproma
    ppbase(jl) = paphp1(jl,kcbot(jl))

    if ( .not. ldcum(jl) ) cycle

    zrrho = (rgas*ptenh(jl,kk+1)*(d_one+vtmpc1*pqenh(jl,kk+1)))/paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
    zpmid = d_half*(ppbase(jl)+paphp1(jl,kctop0(jl)))
  
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
    llo1 = kk < kcbot(jl) .and. ldcum(jl)
    pdmfde(jl) = merge(zentr,d_zero,llo1)

    select case(ktype(jl))

    case(1) 
      llo2 = llo1 .and.                            &
           (kk >= iklwmin .or. paphp1(jl,kk) > zpmid)
      if ( llo2 ) pdmfen(jl) = zentr
  !
  !       organized detrainment, detrainment starts at khmin
  !
      ikb = kcbot(jl)
      podetr(jl,kk) = d_zero
      if ( llo1 .and. kk <= khmin(jl) .and. kk >= kctop0(jl) ) then
        ! AMT: changed detrainment profile to linear instead of TAN
        !      to speed up this section of code by factor of 5
        ikt = kctop0(jl)
        ikh = khmin(jl)
        zarg=(paphp1(jl,ikh)-paphp1(jl,kk))/max(d_one,(paphp1(jl,ikh)-paphp1(jl,ikt)))
        podetr(jl,kk)=zarg*centrmax*pmfu(jl,kk+1)*zdprho
      end if

    case(2)
      llo2 = llo1 .and.                            &
           (ppbase(jl)-paphp1(jl,kk) < 0.2D5 .or. paphp1(jl,kk) > zpmid)
      pdmfen(jl) = merge(zentr,d_zero,llo2)
      iklwmin = max(klwmin(jl),kctop0(jl)+2)

    case(3)
      llo2 = llo1 .and. kk >= iklwmin
      if ( llo2 ) pdmfen(jl) = zentr
      if ( llo2 .and. pqenh(jl,kk+1) > 1.D-5 ) then
        pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
        zentest = max(pqte(jl,kk),d_zero)/pqenh(jl,kk+1)
        zentest = min(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
        pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
      end if
    case default
      cycle
    end select

  end do

  end if !AMT temporary if statement for optimisation
!
  call time_end(subroutine_name,idindx)
  end subroutine cuentr
!
  subroutine cuentrt(kproma,kbdim,klev,klevp1,kk,ptenh,pqenh,pqte,  &
                     paphp1,klwmin,ldcum,ktype,kcbot,kctop0,ppbase, &
                     pmfu,pentr,pdmfen,pdmfde)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!
!          EXTERNALS
!          ---------
!          NONE
!
  implicit none
!
  integer , intent(in) :: kbdim , kk , klev , klevp1 , kproma
  integer , dimension(kbdim) :: kcbot , kctop0 , klwmin , ktype
  logical , dimension(kbdim) :: ldcum
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
  real(dp) , dimension(kbdim,klev) :: pmfu , pqenh , pqte , ptenh
  intent (in) kcbot , kctop0 , klwmin , ktype , ldcum , paphp1 , &
              pentr , pqenh , pqte , ptenh
  intent (out) pdmfde , pdmfen
  intent (inout) pmfu , ppbase
!
  integer :: iklwmin , jl
  logical :: llo1 , llo2
  real(dp) :: zdprho , zentest , zentr , zpmid , zrrho
  character (len=64) :: subroutine_name='cuentrt'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!     -------------------------------------------
!
!
!
!*    1.1          SPECIFY ENTRAINMENT RATES FOR SHALLOW CLOUDS
!     --------------------------------------------
!
!
!
!*    1.2          SPECIFY ENTRAINMENT RATES FOR DEEP CLOUDS
!     -----------------------------------------
!
  do jl = 1 , kproma
    ppbase(jl) = paphp1(jl,kcbot(jl))
    zrrho = (rgas*ptenh(jl,kk+1)*(d_one+vtmpc1*pqenh(jl,kk+1)))/paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
    zpmid = d_half*(ppbase(jl)+paphp1(jl,kctop0(jl)))
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
    llo1 = kk < kcbot(jl) .and. ldcum(jl)
    pdmfde(jl) = merge(zentr,d_zero,llo1)
    llo2 = llo1 .and. ktype(jl) == 2 .and.                            &
           (ppbase(jl)-paphp1(jl,kk) < 0.2D5 .or. paphp1(jl,kk) > zpmid)
    pdmfen(jl) = merge(zentr,d_zero,llo2)
    iklwmin = max(klwmin(jl),kctop0(jl)+2)
    llo2 = llo1 .and. (ktype(jl) == 1 .or. ktype(jl) == 3) .and.        &
           (kk >= iklwmin .or. paphp1(jl,kk) > zpmid)
    if ( llo2 ) pdmfen(jl) = zentr
    if ( llo2 .and. pqenh(jl,kk+1) > 1.D-5 ) then
      pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
      zentest = max(pqte(jl,kk),d_zero)/pqenh(jl,kk+1)
      zentest = min(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
      pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
    end if
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cuentrt
!
  subroutine cuflx(kproma,kbdim,klev,klevp1,pqen,pqsen,ptenh,pqenh, &
                   ktrac,pxtenh,pmfuxt,pmfdxt,paphp1,pgeoh,kcbot,   &
                   kctop,kdtop,ktype,lddraf,ldcum,pmfu,pmfd,pmfus,  &
                   pmfds,pmfuq,pmfdq,pmful,pdmfup,pdmfdp,prfl,prain,&
                   pcpcu,pten,psfl,pdpmel,ktopm2)
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!
!          PURPOSE
!          -------
!
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!
!          EXTERNALS
!          ---------
!          NONE
!
  implicit none
!
  integer , intent(in) :: kbdim , klev , klevp1 , kproma , ktrac
  integer , intent(inout):: ktopm2
  integer , dimension(kbdim) :: kcbot , kctop , kdtop , ktype
  logical , dimension(kbdim) :: ldcum , lddraf
  real(dp) , dimension(kbdim,klevp1) :: paphp1
  real(dp) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pdmfup ,      &
         pdpmel , pgeoh , pmfd , pmfdq , pmfds , pmfu , pmful ,     &
         pmfuq , pmfus , pqen , pqenh , pqsen , pten , ptenh
  real(dp) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtenh
  real(dp) , dimension(kbdim) :: prain , prfl , psfl
  intent (in) kcbot , kctop , kdtop , ldcum , paphp1 , pcpcu , pgeoh , &
              pqen , pqenh , pqsen , pten , ptenh , pxtenh
  intent (out) pdpmel
  intent (inout) ktype , lddraf , pdmfdp , pdmfup , pmfd , pmfdq , &
                 pmfds , pmfdxt , pmfu , pmful , pmfuq , pmfus ,   &
                 pmfuxt , prain , prfl , psfl
!
  integer :: ikb , jk , jl , jt
  real(dp) :: zcons1 , zcons2 , zcucov , zdpevap , zdrfl , zfac ,    &
             zrfl , zrfln , zrmin , zrnew , zrsum , zsnmlt ,        &
             ztmelp2 , zzp
  real(dp) , dimension(kbdim) :: zpsubcl
  character (len=64) :: subroutine_name='cuflx'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)
!
!*    SPECIFY CONSTANTS
!
  zcons1 = cpd/(wlhf*egrav*dtcum)
  zcons2 = d_one/(egrav*dtcum)
  zcucov = 0.050D0
  ztmelp2 = tzero + 2.0D0
!
!
!*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
!     ---------------------------------
!
!     itop=klev
  do jl = 1 , kproma
!       itop=MIN(itop,kctop(jl))
    if ( .not.ldcum(jl) .or. kdtop(jl) < kctop(jl) ) lddraf(jl)       &
         = .false.
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
          pmfd(jl,jk) = d_zero
          pmfds(jl,jk) = d_zero
          pmfdq(jl,jk) = d_zero
          pdmfdp(jl,jk-1) = d_zero
        end if
      end if
    end do
!
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
!
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
!
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
!
  end do
!
!
!*    2.            CALCULATE RAIN/SNOW FALL RATES
!*    CALCULATE MELTING OF SNOW
!*    CALCULATE EVAPORATION OF PRECIP
!     -------------------------------
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
            zfac = zcons1*(d_one+vtmpc2*pqen(jl,jk))                &
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
      if ( ldcum(jl) .and. jk >= kcbot(jl) .and. zpsubcl(jl) > 1.D-20 ) then
        zrfl = zpsubcl(jl)
        zrnew = (max(d_zero,sqrt(zrfl/zcucov)-cevapcu(jk)* &
                 (paphp1(jl,jk+1)-paphp1(jl,jk))           &
                *max(d_zero,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
        zrmin = zrfl - zcucov*max(d_zero,0.80D0*pqsen(jl,jk)         &
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
    prfl(jl) = prfl(jl) + zdpevap*prfl(jl)*(d_one/max(1.D-20,zrsum))
    psfl(jl) = psfl(jl) + zdpevap*psfl(jl)*(d_one/max(1.D-20,zrsum))
  end do
!
  call time_end(subroutine_name,idindx)
  end subroutine cuflx
!
  subroutine cuadjtq(kproma,kbdim,klev,kk,pp,pt,pq,ldflag,kcall)
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          D.SALMOND         CRAY(UK))      12/8/91
!
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *CUBASE*   (T AND Q AT CONDENSTION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!          NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
!               KCALL=0    ENV. T AND QS IN*CUINI*
!               KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!               KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!
!          EXTERNALS
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SETPHYS*.
!
  implicit none

  !  scalar arguments with intent(in):
  integer, intent (in) :: kcall, kk, klev, kproma, kbdim

  !  array arguments with intent(in):
  real(dp), intent (in) :: pp(kbdim)
  logical, intent (in) :: ldflag(kbdim)

  !  array arguments with intent(inout):
  real(dp), intent (inout) :: pq(kbdim,klev), pt(kbdim,klev)

  !  local scalars: 
  real(dp):: zcond1, zqst1, zdqsdt, zqsat, zes, zcor, zlcdqsdt
  integer :: isum, jl, it, it1
  logical :: lo

  !  local arrays: 
  real(dp):: zcond(kbdim)

  !  Executable statements 
  character (len=64) :: subroutine_name='cuadjtq'
  integer :: idindx=0

  call time_begin(subroutine_name,idindx)

  lookupoverflow = .false.

  zcond = d_zero
!
!----------------------------------------------------------------------
!
!     2.           calculate condensation and adjust t and q accordingly
!                  -----------------------------------------------------
!
  if ( kcall == 1 ) then

    isum = 0
    do jl = 1 , kproma
      if ( ldflag(jl) ) then
        it = nint(pt(jl,kk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/pp(jl)
        zes = min(d_half,zes)
        lo = zes < 0.4D0
        zcor = d_one/(d_one-vtmpc1*zes)
        zqsat = zes*zcor
        it1 = it+1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/pp(jl)
        zqst1 = min(d_half,zqst1)
        zqst1 = zqst1/(d_one-vtmpc1*zqst1)
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
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

    isum = 0 !AMT fudge to only make one iteration for temporary speed fix

    if ( isum == 0 ) go to 230

    do jl = 1 , kproma
      if ( ldflag(jl) ) then
        if ( abs(zcond(jl)) > d_zero ) then
          it = nint(pt(jl,kk)*d_1000)
          if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
          it = max(min(it,jptlucu2),jptlucu1)
          zes = tlucua(it)/pp(jl)
          zes = min(d_half,zes)
          lo = zes < 0.4D0
          zcor = d_one/(d_one-vtmpc1*zes)
          zqsat = zes*zcor
          it1 = it+1
          it1 = max(min(it1,jptlucu2),jptlucu1)
          zqst1 = tlucua(it1)/pp(jl)
          zqst1 = min(d_half,zqst1)
          zqst1 = zqst1/(d_one-vtmpc1*zqst1)
          zdqsdt = (zqst1-zqsat)*d_1000
          zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
          pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end if
      end if
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

230 continue

  end if

  if ( kcall == 2 ) then

    isum = 0
    do jl = 1 , kproma
      if ( ldflag(jl) ) then
        it = nint(pt(jl,kk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/pp(jl)
        zes = min(d_half,zes)
        lo = zes < 0.4D0
        zcor = d_one/(d_one-vtmpc1*zes)
        zqsat = zes*zcor
        it1 = it+1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/pp(jl)
        zqst1 = min(d_half,zqst1)
        zqst1 = zqst1/(d_one-vtmpc1*zqst1)
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
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

    isum = 0 !AMT fudge to only make one iteration for temporary speed fix

    if(isum == 0) go to 330

    do jl = 1 , kproma
      if ( ldflag(jl) .and. abs(zcond(jl) ) > d_zero) then
        it = nint(pt(jl,kk)*d_1000)
        if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
        it = max(min(it,jptlucu2),jptlucu1)
        zes = tlucua(it)/pp(jl)
        zes = min(d_half,zes)
        lo = zes < 0.4D0
        zcor = d_one/(d_one-vtmpc1*zes)
        zqsat = zes*zcor
        it1 = it+1
        it1 = max(min(it1,jptlucu2),jptlucu1)
        zqst1 = tlucua(it1)/pp(jl)
        zqst1 = min(d_half,zqst1)
        zqst1 = zqst1/(d_one-vtmpc1*zqst1)
        zdqsdt = (zqst1-zqsat)*d_1000
        zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
        pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
        pq(jl,kk) = pq(jl,kk)-zcond1
      end if
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

330 continue

  end if

  if ( kcall == 0 ) then

    isum = 0
    do jl = 1 , kproma
      it = nint(pt(jl,kk)*d_1000)
      if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/pp(jl)
      zes = min(d_half,zes)
      lo = zes < 0.4D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it+1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/pp(jl)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
      zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
      pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
      pq(jl,kk) = pq(jl,kk)-zcond(jl)
      if ( abs(zcond(jl)) > d_zero ) isum = isum+1
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

    isum = 0 !AMT fudge to only make one iteration for temporary speed fix

    if ( isum == 0 ) go to 430

    do jl = 1 , kproma
      it = nint(pt(jl,kk)*d_1000)
      if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/pp(jl)
      zes = min(d_half,zes)
      lo = zes < 0.4D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it+1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/pp(jl)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
      zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
      pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
      pq(jl,kk) = pq(jl,kk)-zcond1
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

430 continue

  end if

  if ( kcall == 4 ) then

    do jl = 1 , kproma
      it = nint(pt(jl,kk)*d_1000)
      if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/pp(jl)
      zes = min(d_half,zes)
      lo = zes < 0.4D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it+1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/pp(jl)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
      zcond(jl) = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
      pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond(jl)
      pq(jl,kk) = pq(jl,kk)-zcond(jl)
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

    do jl = 1 , kproma
      it = nint(pt(jl,kk)*d_1000)
      if ( it < jptlucu1 .or. it > jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/pp(jl)
      zes = min(d_half,zes)
      lo = zes < 0.4D0
      zcor = d_one/(d_one-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it+1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/pp(jl)
      zqst1 = min(d_half,zqst1)
      zqst1 = zqst1/(d_one-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*d_1000
      zlcdqsdt = merge(zdqsdt*tlucuc(it),zqsat*zcor*tlucub(it),lo)
      zcond1 = (pq(jl,kk)-zqsat)/(d_one+zlcdqsdt)
      pt(jl,kk) = pt(jl,kk)+tlucuc(it)*zcond1
      pq(jl,kk) = pq(jl,kk)-zcond1
    end do

    if (lookupoverflow) then 
      call fatal(__FILE__,__LINE__,'cumulus tables lookup error: overflow')
    endif

  end if
  call time_end(subroutine_name,idindx)
  end subroutine cuadjtq

end module mod_cu_tiedtke
