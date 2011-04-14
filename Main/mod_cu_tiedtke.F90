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
  use m_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_constants
  use mod_cu_tables
  use mod_message
  use mod_date
!
  private
!
  public :: tiedtkedrv
!
  public :: entrpen , entrscv , entrmid , entrdd , cmfctop , cmfcmax , &
            cmfcmin , cmfdeps , rhcdd , cprcon , iconv , nmctop ,      &
            lmfpen , lmfscv , lmfmid , lmfdd , lmfdudv
!
  ! evaporation coefficient for kuo0
  real(dp) , allocatable , dimension(:) :: cevapcu

  real(dp) , parameter :: centrmax = 3.0D-4

  real(dp) :: entrpen      !    entrainment rate for penetrative convection
  real(dp) :: entrscv      !    entrainment rate for shallow convection
  real(dp) :: entrmid      !    entrainment rate for midlevel convection
  real(dp) :: entrdd       !    entrainment rate for cumulus downdrafts
  real(dp) :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  real(dp) :: cmfcmax      !    maximum massflux value allowed for
  real(dp) :: cmfcmin      !    minimum massflux value (for safety)
  real(dp) :: cmfdeps      !    fractional massflux for downdrafts at lfs
  real(dp) :: rhcdd        !    relative saturation in downdrafts
  real(dp) :: cprcon       !    coefficients for determining conversion
                           !    from cloud water to rain
  integer :: iconv
  integer :: nmctop    !  max. level for cloud base of mid level conv.
  logical :: lmfpen    !  true if penetrative convection is switched on
  logical :: lmfscv    !  true if shallow convection is switched on
  logical :: lmfmid    !  true if midlevel convection is switched on
  logical :: lmfdd     !  true if cumulus downdraft is switched on
  logical :: lmfdudv   !  true if cumulus friction is switched on
!
  contains
!
! This subroutines calls cucall
!
  subroutine tiedtkedrv(j)
    implicit none
    integer , intent(in) :: j
    !
    ! CALL CUCALL
    !
  end subroutine tiedtkedrv
!
  subroutine cucall(kproma,kbdim,klev,klevp1,klevm1,ilab,ktrac,     &
                    pxtm1,pxtte,ptm1,pqm1,pum1,pvm1,pxlm1,pxim1,    &
                    ptte,pqte,pvom,pvol,pxlte,pxite,pverv,pxtec,    &
                    pqtec,pqhfla,papp1,paphp1,pgeo,prsfc,pssfc,     &
                    paprc,paprs,ktype,ldland,ptopmax)
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: papp1 , pgeo , pqm1 , pqte ,   &
         pqtec , ptm1 , ptte , pum1 , pverv , pvm1 , pvol , pvom ,  &
         pxim1 , pxite , pxlm1 , pxlte , pxtec
  real(8) , dimension(kbdim) :: paprc , paprs , pqhfla , prsfc ,    &
                                pssfc , ptopmax
  real(8) , dimension(kbdim,klev,ktrac) :: pxtm1 , pxtte
  intent (in) papp1 , pqm1 , ptm1 , pum1 , pvm1 , pxim1 , pxite ,   &
              pxlm1 , pxlte , pxtm1
  intent (inout) ptopmax
!
  integer , dimension(kbdim) :: icbot , ictop , itopec2
  integer :: ilevmin , it , jk , jl , jt
  logical , dimension(kbdim) :: locum
  real(8) , dimension(kbdim,klev) :: zlu , zlude , zmfd , zmfu ,    &
         zqp1 , zqsat , zqu , zqude , ztp1 , ztu , ztvp1 , zup1 ,   &
         zvp1 , zxp1
  real(8) , dimension(kbdim) :: zrain , ztopmax
  real(8) :: ztmst , zxip1 , zxlp1
  real(8) , dimension(kbdim,klev,ktrac) :: zxtp1 , zxtu
!
!  Executable statements
!
  lookupoverflow = .false.
!
!-----------------------------------------------------------------------
!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*    -----------------------------------
!
  ztmst = dt
  do jk = 1 , klev
    do jl = 1 , kproma
      ztp1(jl,jk) = ptm1(jl,jk) + ptte(jl,jk)*ztmst
      zqp1(jl,jk) = max(0.0D0,pqm1(jl,jk)+pqte(jl,jk)*ztmst)
      zxlp1 = pxlm1(jl,jk) + pxlte(jl,jk)*ztmst
      zxip1 = pxim1(jl,jk) + pxite(jl,jk)*ztmst
      zxp1(jl,jk) = max(0.0D0,zxlp1+zxip1)
      ztvp1(jl,jk) = ztp1(jl,jk)                                    &
                     *(1.0D0+vtmpc1*zqp1(jl,jk)-zxp1(jl,jk))
      zup1(jl,jk) = pum1(jl,jk) + pvom(jl,jk)*ztmst
      zvp1(jl,jk) = pvm1(jl,jk) + pvol(jl,jk)*ztmst
      it = int(ztp1(jl,jk)*1000.0D0)
      if ( it<jptlucu1 .or. it>jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zqsat(jl,jk) = tlucua(it)/papp1(jl,jk)
      zqsat(jl,jk) = min(0.50D0,zqsat(jl,jk))
      zqsat(jl,jk) = zqsat(jl,jk)/(1.0D0-vtmpc1*zqsat(jl,jk))
    end do
 
    if ( lookupoverflow ) then
      call fatal(__FILE__,__LINE__, &
           'Cumulus Tables lookup error: OVERFLOW')
    end if
 
    do jt = 1 , ktrac
      do jl = 1 , kproma
        zxtp1(jl,jk,jt) = pxtm1(jl,jk,jt) + pxtte(jl,jk,jt)*ztmst
      end do
    end do
 
  end do
  do jl = 1 , kproma
    zrain(jl) = 0.0D0
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
                 ktype,icbot,ictop,ztu,zqu,zlu,zlude,zmfu,zmfd,     &
                 zrain)
  case (2)
    call cumastrt(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                  zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                  pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                  pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,pqtec,    &
                  zqude,locum,ktype,icbot,ictop,ztu,zqu,zlu,zlude,  &
                  zmfu,zmfd,zrain)
  case (3)
    call cumastrh(kproma,kbdim,klev,klevp1,klevm1,ilab,ztp1,zqp1,   &
                  zxp1,zup1,zvp1,ztvp1,ktrac,ldland,zxtp1,zxtu,     &
                  pxtte,pverv,zqsat,pqhfla,paphp1,pgeo,ptte,pqte,   &
                  pvom,pvol,prsfc,pssfc,paprc,paprs,pxtec,pqtec,    &
                  zqude,locum,ktype,icbot,ictop,ztu,zqu,zlu,zlude,  &
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
      if ( ilab(jl,jk)==2 .and. itopec2(jl)==klevp1 ) itopec2(jl)   &
           = jk
    end do
  end do
!
  ztopmax(1:kproma) = ptopmax(1:kproma)
 
  do jl = 1 , kproma
    if ( itopec2(jl)==1 ) then
      ptopmax(jl) = papp1(jl,1)
    else if ( itopec2(jl)/=klevp1 ) then
      ptopmax(jl) = paphp1(jl,itopec2(jl))
    else
      ptopmax(jl) = 99999.0D0
    end if
    ptopmax(jl) = min(ptopmax(jl),ztopmax(jl))
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(8) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(8) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ihmin , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(8) :: zalvdcp , zalvs , zb , zbi , zcons2 , zcor , zdepth ,  &
             zdhdz , zdqmin , zdqsdt , zdz , zeps , zes , zfac ,    &
             zgam , zhhat , zhsat , zmfmax , zpbmpt , zqalv ,       &
             zqsat , zqst1 , zqumqe , zrh , zro , ztau , zzz
  real(8) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                zhcbase , zheat , zhmin , zmfub ,   &
                                zmfub1 , zrfl , zsfl
  real(8) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zhhatt , zmfdq , zmfds , zmful , &
         zmfuq , zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , &
         zuu , zvd , zvu , zxenh
  real(8) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
!
!     Executable statements
 
  lookupoverflow = .false.
!
!-----------------------------------------------------------------------
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!     --------------------------------
!
  zcons2 = 1.0D0/(egrav*dt)

! *AMT* NOTE!
! this paramter is the CAPE adjustment timescale which in the global model
! was a function of horizontal resolution (nn wavenumber of a spectral model)
! this is translated roughly into horizontal resolution in meters
! !!!WARNING: this is defined twice and should be removed into a module!!!!
!
!      ztau = min(3.0D0*3600.0D0,7200.0D0*63.0D0/nn)
  ztau=min(3.0D0*3600.0D0,22.7D0*ds)

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
      zdqcv(jl) = zdqcv(jl) + pqte(jl,jk)                           &
                  *(paphp1(jl,jk+1)-paphp1(jl,jk))
      if ( jk>=kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)    &
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
    llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010D0,  &
                llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl)>max(0.0D0,-1.10D0               &
                *pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl)==1)
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
    zalvs = merge(wlhv,wlhs,ptu(jl,ikb)>tzero)
    zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb)         &
                  + zalvs*pqu(jl,ikb)
    ictop0(jl) = kcbot(jl) - 1
  end do
  do jk = klevm1 , 3 , -1
    do jl = 1 , kproma
      zalvs = merge(wlhv,wlhs,ztenh(jl,jk)>tzero)
      zalvdcp = zalvs/zcpcu(jl,jk)
      zqalv = 1.0D0/zalvs
      zhsat = zcpcu(jl,jk)*ztenh(jl,jk) + zgeoh(jl,jk)              &
              + zalvs*zqsenh(jl,jk)
      it = nint(ztenh(jl,jk)*1000.0D0)
      if ( it<jptlucu1 .or. it>jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/paphp1(jl,jk)
      zes = min(0.50D0,zes)
      lo = zes<0.40D0
      zcor = 1.0D0/(1.0D0-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(0.50D0,zqst1)
      zqst1 = zqst1/(1.0D0-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*1000.0D0
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
      zhhat = zhsat - (zzz+zgam*zzz)/(1.0D0+zgam*zzz*zqalv)         &
              *max(zqsenh(jl,jk)-zqenh(jl,jk),0.0D0)
      zhhatt(jl,jk) = zhhat
      if ( jk<ictop0(jl) .and. zhcbase(jl)>zhhat ) ictop0(jl) = jk
    end do
  end do
!!
!!    DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!    (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!     DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             
!     & paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.D4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
!     ENDDO
!!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__, &
         'Cumulus Tables lookup error: OVERFLOW')
  end if
!
!     FIND LOWEST POSSIBLE ORG. DETRAINMENT LEVEL
!     -------------------------------------------
!
  do jl = 1 , kproma
    zhmin(jl) = 0.0D0
    ihmin(jl) = 0
    llo1 = ldcum(jl) .and. ktype(jl)==1
    if ( llo1 ) then
      ikb = kcbot(jl)
      ihmin(jl) = ikb
    end if
  end do
!
  zb = 25.0D0
  zbi = 1.0D0/(zb*egrav)
  do jk = klev , 1 , -1
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl)==1 .and. ihmin(jl)==kcbot(jl)
      if ( llo1 .and. jk<kcbot(jl) .and. jk>=ictop0(jl) ) then
        zalvs = merge(wlhv,wlhs,ztenh(jl,jk)>tzero)
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)                                         &
              /(rgas*ztenh(jl,jk)*(1.0D0+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zdhdz = (zcpen(jl,jk-1)*pten(jl,jk-1)-zcpen(jl,jk)          &
                *pten(jl,jk)+zalvs*(pqen(jl,jk-1)-pqen(jl,jk))      &
                +(pgeo(jl,jk-1)-pgeo(jl,jk)))                       &
                *egrav/(pgeo(jl,jk-1)-pgeo(jl,jk))
        zdepth = zgeoh(jl,jk) - zgeoh(jl,ikb)
        zfac = sqrt(1.0D0+zdepth*zbi)
        zhmin(jl) = zhmin(jl) + zdhdz*zfac*zdz
        zrh = -zalvs*(zqsenh(jl,jk)-zqenh(jl,jk))*zfac
        if ( zhmin(jl)>zrh ) ihmin(jl) = jk
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl)==1 ) then
      if ( ihmin(jl)<ictop0(jl) ) ihmin(jl) = ictop0(jl)
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
    if ( ldcum(jl) .and. ktype(jl)==1 .and. zpbmpt<2.D4 ) ktype(jl) &
         = 2
    if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
    if ( ktype(jl)==2 ) zentr(jl) = entrscv
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
    zheat(jl) = 0.0D0
    zcape(jl) = 0.0D0
    zmfub1(jl) = zmfub(jl)
  end do
!
  do jk = 1 , klev
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl)==1
      if ( llo1 .and. jk<=kcbot(jl) .and. jk>kctop(jl) ) then
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)                                         &
              /(rgas*ztenh(jl,jk)*(1.0D0+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zheat(jl) = zheat(jl)                                       &
                    + ((pten(jl,jk-1)-pten(jl,jk)                   &
                    +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)             &
                    +vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                    *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
        zcape(jl) = zcape(jl)                                       &
                    + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)   &
                    +egrav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))           &
                    -egrav*plu(jl,jk))*zdz
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl)==1 ) then
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
    if ( ktype(jl)==2 ) then
      ikb = kcbot(jl)
      llo1 = pmfd(jl,ikb)<0.0D0 .and. loddraf(jl)
      zeps = merge(cmfdeps,0.0D0,llo1)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)         &
               - (1.0D0-zeps)*zqenh(jl,ikb)
      zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl)   &
             .and. zmfub(jl)<zmfmax
      zmfub1(jl) = merge(zdqpbl(jl)/(egrav                            &
                   *max(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl)        &
                   -zmfub(jl))<0.20D0*zmfub(jl))
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(8) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(8) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(8) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
             zdz , zeps , zes , zfac , zgam , zhhat , zhsat ,       &
             zmfmax , zpbmpt , zqalv , zqsat , zqst1 , zqumqe ,     &
             zro , ztau , zzz
  real(8) , dimension(kbdim) :: zcape , zdqcv , zdqpbl , zentr ,    &
                                zhcbase , zheat , zmfub , zmfub1 ,  &
                                zrfl , zsfl
  real(8) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
         zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
         zvd , zvu , zxenh
  real(8) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
!
!     Executable statements
 
  lookupoverflow = .false.
!-----------------------------------------------------------------------
!     1.           SPECIFY CONSTANTS AND PARAMETERS
!     --------------------------------
!
  zcons2 = 1.0D0/(egrav*dt)
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
      if ( jk>=kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)    &
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
    llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),0.010D0,  &
                llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl)>max(0.0D0,-1.10D0               &
                *pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl)==1)
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
    zalvs = merge(wlhv,wlhs,ptu(jl,ikb)>tzero)
    zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb)         &
                  + zalvs*pqu(jl,ikb)
    ictop0(jl) = kcbot(jl) - 1
  end do
  do jk = klevm1 , 3 , -1
    do jl = 1 , kproma
      zalvs = merge(wlhv,wlhs,ztenh(jl,jk)>tzero)
      zalvdcp = zalvs/zcpcu(jl,jk)
      zqalv = 1.0D0/zalvs
      zhsat = zcpcu(jl,jk)*ztenh(jl,jk) + zgeoh(jl,jk)              &
              + zalvs*zqsenh(jl,jk)
      it = nint(ztenh(jl,jk)*1000.0D0)
      if ( it<jptlucu1 .or. it>jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/paphp1(jl,jk)
      zes = min(0.50D0,zes)
      lo = zes<0.40D0
      zcor = 1.0D0/(1.0D0-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(0.50D0,zqst1)
      zqst1 = zqst1/(1.0D0-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*1000.0D0
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
      zhhat = zhsat - (zzz+zgam*zzz)/(1.0D0+zgam*zzz*zqalv)         &
              *max(zqsenh(jl,jk)-zqenh(jl,jk),0.0D0)
      if ( jk<ictop0(jl) .and. zhcbase(jl)>zhhat ) ictop0(jl) = jk
    end do
  end do
!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__, &
         'Cumulus Tables lookup error: OVERFLOW')
  end if
!!
!!    DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!    (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!     DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             
!     & paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.D4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
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
    if ( ldcum(jl) .and. ktype(jl)==1 .and. zpbmpt<2.D4 ) ktype(jl) &
         = 2
    if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
    if ( ktype(jl)==2 ) zentr(jl) = entrscv
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
    zheat(jl) = 0.0D0
    zcape(jl) = 0.0D0
    zmfub1(jl) = zmfub(jl)
  end do
!
  do jk = 1 , klev
    do jl = 1 , kproma
      llo1 = ldcum(jl) .and. ktype(jl)==1
      if ( llo1 .and. jk<=kcbot(jl) .and. jk>kctop(jl) ) then
        ikb = kcbot(jl)
        zro = paphp1(jl,jk)                                         &
              /(rgas*ztenh(jl,jk)*(1.0D0+vtmpc1*zqenh(jl,jk)))
        zdz = (paphp1(jl,jk)-paphp1(jl,jk-1))/(egrav*zro)
        zheat(jl) = zheat(jl)                                       &
                    + ((pten(jl,jk-1)-pten(jl,jk)                   &
                    +egrav*zdz/zcpcu(jl,jk))/ztenh(jl,jk)             &
                    +vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk)))            &
                    *(egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
        zcape(jl) = zcape(jl)                                       &
                    + (egrav*(ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk)   &
                    +egrav*vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))           &
                    -egrav*plu(jl,jk))*zdz
      end if
    end do
  end do
!
  do jl = 1 , kproma
    if ( ldcum(jl) .and. ktype(jl)==1 ) then
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
    if ( ktype(jl)==2 ) then
      ikb = kcbot(jl)
      llo1 = pmfd(jl,ikb)<0.0D0 .and. loddraf(jl)
      zeps = merge(cmfdeps,0.0D0,llo1)
      zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)         &
               - (1.0D0-zeps)*zqenh(jl,ikb)
      zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
      zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
      llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl)   &
             .and. zmfub(jl)<zmfmax
      zmfub1(jl) = merge(zdqpbl(jl)/(egrav                            &
             *max(zqumqe,zdqmin)),zmfub(jl),llo1)
      zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),abs(zmfub1(jl)        &
                   -zmfub(jl))<0.20D0*zmfub(jl))
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: paprc , paprs , pqhfla , prain ,    &
                                prsfc , pssfc
  real(8) , dimension(kbdim,klev) :: pgeo , plu , plude , pmfd ,    &
         pmfu , pqen , pqsen , pqte , pqtec , pqu , pqude , pten ,  &
         ptte , ptu , ptven , puen , pven , pverv , pvol , pvom ,   &
         pxen , pxtec
  real(8) , dimension(kbdim,klev,ktrac) :: pxten , pxtte , pxtu
  intent (in) pqhfla
  intent (inout) ktype , ldcum , pmfd
!
  integer , dimension(kbdim) :: ictop0 , idtop , ilwmin
  integer :: ikb , it , it1 , itopm2 , jk , jl , jt
  logical :: llo1 , lo
  logical , dimension(kbdim) :: loddraf
  real(8) :: zalvdcp , zalvs , zcons2 , zcor , zdqmin , zdqsdt ,    &
             zeps , zes , zfac , zgam , zhhat , zhsat , zmfmax ,    &
             zpbmpt , zqalv , zqsat , zqst1 , zqumqe , zzz
  real(8) , dimension(kbdim,klev) :: zcpcu , zcpen , zdmfdp ,       &
         zdmfup , zdpmel , zgeoh , zmfdq , zmfds , zmful , zmfuq ,  &
         zmfus , zqd , zqenh , zqsenh , ztd , ztenh , zud , zuu ,   &
         zvd , zvu , zxenh
  real(8) , dimension(kbdim) :: zdqcv , zdqpbl , zentr , zhcbase ,  &
                                zmfub , zmfub1 , zrfl , zsfl
  real(8) , dimension(kbdim,klev,ktrac) :: zmfdxt , zmfuxt , zxtd , &
         zxtenh
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
  zcons2 = 1.0D0/(egrav*dt)
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
      if ( jk>=kcbot(jl) ) zdqpbl(jl) = zdqpbl(jl) + pqte(jl,jk)    &
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
    llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl)
    zmfub(jl) = merge(zdqpbl(jl)                                    &
                /(egrav*max(zqumqe,zdqmin)),0.010D0,llo1)
    zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
    zmfub(jl) = min(zmfub(jl),zmfmax)
    if ( .not.llo1 ) ldcum(jl) = .false.
    ktype(jl) = merge(1,2,zdqcv(jl)>max(0.0D0,-1.10D0               &
                *pqhfla(jl)*egrav))
    zentr(jl) = merge(entrpen,entrscv,ktype(jl)==1)
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
    zalvs = merge(wlhv,wlhs,ptu(jl,ikb)>tzero)
    zhcbase(jl) = zcpcu(jl,ikb)*ptu(jl,ikb) + zgeoh(jl,ikb)         &
                  + zalvs*pqu(jl,ikb)
    ictop0(jl) = kcbot(jl) - 1
  end do
  do jk = klevm1 , 3 , -1
    do jl = 1 , kproma
      zalvs = merge(wlhv,wlhs,ztenh(jl,jk)>tzero)
      zalvdcp = zalvs/zcpcu(jl,jk)
      zqalv = 1.0D0/zalvs
      zhsat = zcpcu(jl,jk)*ztenh(jl,jk) + zgeoh(jl,jk)              &
              + zalvs*zqsenh(jl,jk)
      it = nint(ztenh(jl,jk)*1000.0D0)
      if ( it<jptlucu1 .or. it>jptlucu2 ) lookupoverflow = .true.
      it = max(min(it,jptlucu2),jptlucu1)
      zes = tlucua(it)/paphp1(jl,jk)
      zes = min(0.50D0,zes)
      lo = zes<0.40D0
      zcor = 1.0D0/(1.0D0-vtmpc1*zes)
      zqsat = zes*zcor
      it1 = it + 1
      it1 = max(min(it1,jptlucu2),jptlucu1)
      zqst1 = tlucua(it1)/paphp1(jl,jk)
      zqst1 = min(0.50D0,zqst1)
      zqst1 = zqst1/(1.0D0-vtmpc1*zqst1)
      zdqsdt = (zqst1-zqsat)*1000.0D0
      zgam = merge(zalvdcp*zdqsdt,zqsat*zcor*tlucub(it),lo)
      zzz = zcpcu(jl,jk)*ztenh(jl,jk)*vtmpc1
      zhhat = zhsat - (zzz+zgam*zzz)/(1.0D0+zgam*zzz*zqalv)         &
              *max(zqsenh(jl,jk)-zqenh(jl,jk),0.0D0)
      if ( jk<ictop0(jl) .and. zhcbase(jl)>zhhat ) ictop0(jl) = jk
    end do
  end do
!
  if ( lookupoverflow ) then
    call fatal(__FILE__,__LINE__, &
         'Cumulus Tables lookup error: OVERFLOW')
  end if
!!
!!    DEEP CONVECTION IF CLOUD DEPTH > 200 HPA, ELSE SHALLOW
!!    (CLOUD DEPTH FROM NON-ENTRAINIG PLUME)
!!
!     DO jl=1,kproma
!     ktype(jl)=MERGE(1,2,                                             &
!     paphp1(jl,kcbot(jl))-paphp1(jl,ictop0(jl)).gt.2.D4)
!     zentr(jl)=MERGE(entrpen,entrscv,ktype(jl).eq.1)
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
    if ( ldcum(jl) .and. ktype(jl)==1 .and. zpbmpt<2.D4 ) ktype(jl) &
         = 2
    if ( ldcum(jl) ) ictop0(jl) = kctop(jl)
    if ( ktype(jl)==2 ) zentr(jl) = entrscv
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
        llo1 = pmfd(jl,ikb)<0.0D0
        zeps = merge(cmfdeps,0.0D0,llo1)
        zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zeps*zqd(jl,ikb)       &
                 - (1.0D0-zeps)*zqenh(jl,ikb)
        zdqmin = max(0.010D0*zqenh(jl,ikb),1.D-10)
        zmfmax = (paphp1(jl,ikb)-paphp1(jl,ikb-1))*zcons2
        llo1 = zdqpbl(jl)>0.0D0 .and. zqumqe>zdqmin .and. ldcum(jl) &
               .and. zmfub(jl)<zmfmax
        zmfub1(jl) = merge(zdqpbl(jl)/(egrav*max(zqumqe,zdqmin)),     &
                     zmfub(jl),llo1)
        zmfub1(jl) = merge(zmfub1(jl),zmfub(jl),(ktype(jl)==1 .or.  &
                     ktype(jl)==2) .and. abs(zmfub1(jl)-zmfub(jl))  &
                     <0.20D0*zmfub(jl))
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfdp ,       &
         pdmfup , pdpmel , pgeo , pgeoh , plu , plude , pmfd ,      &
         pmfdq , pmfds , pmfu , pmfuq , pmfus , pqd , pqen , pqenh ,&
         pqsen , pqsenh , pqu , pqude , ptd , pten , ptenh , ptu ,  &
         ptven , pud , puen , puu , pvd , pven , pverv , pvu ,      &
         pxen , pxenh
  real(8) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtd , &
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
  real(8) :: zarg , zcpm , zzs
  real(8) , dimension(kbdim) :: zph , zwmax
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
      zcpm = (pcpen(jl,jk)+pcpen(jl,jk-1))*0.50D0
      ptenh(jl,jk) = (max(pcpen(jl,jk-1)*pten(jl,jk-1)+pgeo(jl,jk-1)&
                     ,pcpen(jl,jk)*pten(jl,jk)+pgeo(jl,jk))         &
                     -pgeoh(jl,jk))/zcpm
      pqsenh(jl,jk) = pqsen(jl,jk-1)
      zph(jl) = paphp1(jl,jk)
      loflag(jl) = .true.
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pxtenh(jl,jk,jt) = (pxten(jl,jk,jt)+pxten(jl,jk-1,jt))      &
                           *0.50D0
      end do
    end do
!
!
    ik = jk
    icall = 0
    call cuadjtq(kproma,kbdim,klev,ik,zph,ptenh,pqsenh,loflag,icall)
!
    do jl = 1 , kproma
      pxenh(jl,jk) = (pxen(jl,jk)+pxen(jl,jk-1))*0.50D0
      pqenh(jl,jk) = min(pqen(jl,jk-1),pqsen(jl,jk-1))              &
                     + (pqsenh(jl,jk)-pqsen(jl,jk-1))
      pqenh(jl,jk) = max(pqenh(jl,jk),0.0D0)
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
    zwmax(jl) = 0.0D0
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
      zzs = max(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk),             &
            pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
      ptenh(jl,jk) = (zzs-pgeoh(jl,jk))/pcpcu(jl,jk)
    end do
  end do
!
  do jk = klev , 3 , -1
    do jl = 1 , kproma
      if ( pverv(jl,jk)<zwmax(jl) ) then
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
    if ( jk==1 ) ik = 1
    do jl = 1 , kproma
      ptu(jl,jk) = ptenh(jl,jk)
      ptd(jl,jk) = ptenh(jl,jk)
      pqu(jl,jk) = pqenh(jl,jk)
      pqd(jl,jk) = pqenh(jl,jk)
      plu(jl,jk) = 0.0D0
      puu(jl,jk) = puen(jl,ik)
      pud(jl,jk) = puen(jl,ik)
      pvu(jl,jk) = pven(jl,ik)
      pvd(jl,jk) = pven(jl,ik)
      pmfu(jl,jk) = 0.0D0
      pmfd(jl,jk) = 0.0D0
      pmfus(jl,jk) = 0.0D0
      pmfds(jl,jk) = 0.0D0
      pmfuq(jl,jk) = 0.0D0
      pmfdq(jl,jk) = 0.0D0
      pdmfup(jl,jk) = 0.0D0
      pdmfdp(jl,jk) = 0.0D0
      pdpmel(jl,jk) = 0.0D0
      plude(jl,jk) = 0.0D0
      pqude(jl,jk) = 0.0D0
      klab(jl,jk) = 0
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pxtu(jl,jk,jt) = pxtenh(jl,jk,jt)
        pxtd(jl,jk,jt) = pxtenh(jl,jk,jt)
        pmfuxt(jl,jk,jt) = 0.0D0
        pmfdxt(jl,jk,jt) = 0.0D0
      end do
    end do
!
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
         pgeoh , phhatt , plu , plude , pmfu , pmful , pmfuq ,      &
         pmfus , pqen , pqenh , pqsen , pqsenh , pqte , pqu ,       &
         pqude , pten , ptenh , ptu , puen , puu , pven , pverv ,   &
         pvu
  real(8) , dimension(kbdim) :: pentr , phcbase , pmfub
  real(8) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
         pxtu
  intent (in) ldland , pcpcu , phcbase , phhatt , pqsenh , pxtenh
  intent (inout) kcbot , kctop , kctop0 , klab , ktype , ldcum ,    &
                 plu , plude , pmfu , pmfub , pmful , pmfuq ,       &
                 pmfus , pmfuxt , pqu , pqude , ptu , puu , pvu ,   &
                 pxtu
!
  integer :: icall , ik , ikb , ikt , jk , jl , jt
  logical , dimension(kbdim) :: loflag
  real(8) :: zalvs , zbuo , zbuoyz , zcons2 , zdlev , zdmfdu ,      &
             zdmfeu , zdnoprc , zdprho , zdrodz , zdt , zdz , zfac ,&
             zga , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,     &
             zmfusk , zmfuxtk , zmse , znevn , zodmax , zprcon ,    &
             zqcod , zqeen , zqude , zscde , zscod , zseen ,        &
             ztglace , zxteen , zxtude , zz , zzdmf
  real(8) , dimension(kbdim) :: zbuoy , zdmfde , zdmfen , zmfuu ,   &
                                zmfuv , zpbase , zph , zqold
  real(8) , dimension(kbdim,klev) :: zodetr , zoentr
!
!----------------------------------------------------------------------
!*    1.           SPECIFY PARAMETERS
!     ------------------
!
  zcons2 = 1.0D0/(egrav*dt)
  ztglace = tzero - 13.0D0
  zqold(1:kproma) = 0.0D0

! AMT NOTE!!! in the original scheme, this level which restricts rainfall 
! below a certain pressure (from the surface) is hard wired according to the 
! vertical resolution of the model - This has been fixed to 150 hPa.
! WARNING - this should be set outside of the routines.
!      if ( klev/=11 ) then
!        zdlev = 3.0D4
!      else if ( nn==21 ) then
!        zdlev = 1.5D4
!      else if ( nn==31 ) then
!        zdlev = 2.0D4 
!      else
!        zdlev = 3.0D4
!      end if
   zdlev = 1.5D4

!
!----------------------------------------------------------------------
!     2.           SET DEFAULT VALUES
!     ------------------
!
  do jl = 1 , kproma
    zmfuu(jl) = 0.0D0
    zmfuv(jl) = 0.0D0
    if ( .not.ldcum(jl) ) ktype(jl) = 0
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      plu(jl,jk) = 0.0D0
      pmfu(jl,jk) = 0.0D0
      pmfus(jl,jk) = 0.0D0
      pmfuq(jl,jk) = 0.0D0
      pmful(jl,jk) = 0.0D0
      plude(jl,jk) = 0.0D0
      pqude(jl,jk) = 0.0D0
      pdmfup(jl,jk) = 0.0D0
      if ( .not.ldcum(jl) .or. ktype(jl)==3 ) klab(jl,jk) = 0
      if ( .not.ldcum(jl) .and. paphp1(jl,jk)<4.D4 ) kctop0(jl) = jk
      if ( jk<kcbot(jl) ) klab(jl,jk) = 0
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pmfuxt(jl,jk,jt) = 0.0D0
      end do
    end do
!
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      zoentr(jl,jk) = 0.0D0
      zodetr(jl,jk) = 0.0D0
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
      pmfub(jl) = 0.0D0
      pqu(jl,klev) = 0.0D0
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
      if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = 0.0D0
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
    if ( ktype(jl)==1 ) then
      ikb = kcbot(jl)
      zbuoy(jl) = egrav*(ptu(jl,ikb)-ptenh(jl,ikb))/ptenh(jl,ikb)     &
                  + egrav*vtmpc1*(pqu(jl,ikb)-pqenh(jl,ikb))
      if ( zbuoy(jl)>0.0D0 ) then
        zdz = (pgeo(jl,ikb-1)-pgeo(jl,ikb))*regrav
        zdrodz = -log(pten(jl,ikb-1)/pten(jl,ikb))                  &
                 /zdz - egrav/(rgas*ptenh(jl,ikb)                     &
                 *(1.0D0+vtmpc1*pqenh(jl,ikb)))
!           nb zoentr is here a fractional value
        zoentr(jl,ikb-1) = zbuoy(jl)*0.50D0/(1.0D0+zbuoy(jl)*zdz)   &
                           + zdrodz
        zoentr(jl,ikb-1) = min(zoentr(jl,ikb-1),centrmax)
        zoentr(jl,ikb-1) = max(zoentr(jl,ikb-1),0.0D0)
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
    if ( lmfmid .and. ik<klevm1 .and. ik>nmctop )                   &
         call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,    &
         puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
         ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
         pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)
!
    do jl = 1 , kproma
      if ( klab(jl,jk+1)==0 ) klab(jl,jk) = 0
      loflag(jl) = klab(jl,jk+1)>0
      zph(jl) = paphp1(jl,jk)
      if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        if ( pmfub(jl)>zmfmax ) then
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
        if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          if ( pmfub(jl)>zmfmax ) then
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
      if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
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
        if ( jk<kcbot(jl) ) then
          zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))      &
                   *zcons2)
          zdmfen(jl) = max(zdmfen(jl)-max(zmftest-zmfmax,0.0D0),    &
                       0.0D0)
        end if
        zdmfde(jl) = min(zdmfde(jl),0.750D0*pmfu(jl,jk+1))
        pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
        if ( ktype(jl)==1 .and. jk<kcbot(jl) ) then
          zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
          zoentr(jl,jk) = zoentr(jl,jk)*zdprho*pmfu(jl,jk+1)
          zmftest = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
          zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))      &
                   *zcons2)
          zoentr(jl,jk) = max(zoentr(jl,jk)-max(zmftest-zmfmax,0.0D0&
                          ),0.0D0)
        else
          zoentr(jl,jk) = 0.0D0
        end if
        if ( ktype(jl)==1 .and. jk<kcbot(jl) .and. jk<=khmin(jl) )  &
             then
!             limit organized detrainment to not allowing for too
!             deep clouds
          zalvs = merge(wlhv,wlhs,ptu(jl,jk+1)>tzero)
          zmse = pcpcu(jl,jk+1)*ptu(jl,jk+1) + zalvs*pqu(jl,jk+1)   &
                 + pgeoh(jl,jk+1)
          ikt = kctop0(jl)
          znevn = (pgeoh(jl,ikt)-pgeoh(jl,jk+1))                    &
                  *(zmse-phhatt(jl,jk+1))*regrav
          if ( znevn<=0.0D0 ) znevn = 1.0D0
          zdprho = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
          zodmax = ((phcbase(jl)-zmse)/znevn)*zdprho*pmfu(jl,jk+1)
          zodmax = max(zodmax,0.0D0)
          zodetr(jl,jk) = min(zodetr(jl,jk),zodmax)
        end if
        zodetr(jl,jk) = min(zodetr(jl,jk),0.750D0*pmfu(jl,jk))
        pmfu(jl,jk) = pmfu(jl,jk) + zoentr(jl,jk) - zodetr(jl,jk)
        zqeen = pqenh(jl,jk+1)*zdmfen(jl)
        zqeen = zqeen + pqenh(jl,jk+1)*zoentr(jl,jk)
        zseen = (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))      &
                *zdmfen(jl)
        zseen = zseen +                                             &
                (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))      &
                *zoentr(jl,jk)
        zscde = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))        &
                *zdmfde(jl)
!           find moist static energy that give nonbuoyant air
        zalvs = merge(wlhv,wlhs,ptenh(jl,jk+1)>tzero)
        zga = zalvs*pqsenh(jl,jk+1)/(rwat*(ptenh(jl,jk+1)**2))
        zdt = (plu(jl,jk+1)-vtmpc1*(pqsenh(jl,jk+1)-pqenh(jl,jk+1)))&
              /(1.0D0/ptenh(jl,jk+1)+vtmpc1*zga)
        zscod = pcpcu(jl,jk+1)*ptenh(jl,jk+1) + pgeoh(jl,jk+1)      &
                + pcpcu(jl,jk+1)*zdt
        zscod = max(zscod,0.0D0)
        zscde = zscde + zodetr(jl,jk)*zscod
        zqude = pqu(jl,jk+1)*zdmfde(jl)
        zqcod = pqsenh(jl,jk+1) + zga*zdt
        zqcod = max(zqcod,0.0D0)
        zqude = zqude + zodetr(jl,jk)*zqcod
        pqude(jl,jk) = zqude
        plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
        plude(jl,jk) = plude(jl,jk) + plu(jl,jk+1)*zodetr(jl,jk)
        zmfusk = pmfus(jl,jk+1) + zseen - zscde
        zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
        zmfulk = pmful(jl,jk+1) - plude(jl,jk)
        plu(jl,jk) = zmfulk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
        pqu(jl,jk) = zmfuqk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
        ptu(jl,jk) = (zmfusk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))       &
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
          pxtu(jl,jk,jt) = zmfuxtk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
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
        if ( pqu(jl,jk)<zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(1.0D0+vtmpc1*pqu(jl,jk)-plu(jl,jk))    &
                 - ptenh(jl,jk)*(1.0D0+vtmpc1*pqenh(jl,jk))
          if ( klab(jl,jk+1)==1 ) zbuo = zbuo + 0.50D0
          if ( zbuo>0.0D0 .and. pmfu(jl,jk)>=0.010D0*pmfub(jl) .and.&
               jk>=kctop0(jl) ) then
            kctop(jl) = jk
            ldcum(jl) = .true.
            zdnoprc = merge(zdlev,1.5D4,ldland(jl))
            zprcon = merge(0.0D0,cprcon,zpbase(jl)-paphp1(jl,jk)    &
                     <zdnoprc)
            zlnew = plu(jl,jk)                                      &
                    /(1.0D0+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
            pdmfup(jl,jk) = max(0.0D0,(plu(jl,jk)-zlnew)*pmfu(jl,jk)&
                            )
            plu(jl,jk) = zlnew
          else
            klab(jl,jk) = 0
            pmfu(jl,jk) = 0.0D0
          end if
        end if
      end if
    end do
    do jl = 1 , kproma
      if ( loflag(jl) ) then
        pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
        pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))       &
                       *pmfu(jl,jk)
        pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
      end if
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( loflag(jl) ) pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)         &
             *pmfu(jl,jk)
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
          if ( ktype(jl)==1 .or. ktype(jl)==3 ) then
            zz = merge(3.0D0,2.0D0,zdmfen(jl)==0.0D0)
          else
            zz = merge(1.0D0,0.0D0,zdmfen(jl)==0.0D0)
          end if
          zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
          zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
          zdmfdu = min(zdmfdu,0.750D0*pmfu(jl,jk+1))
          zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk)                &
                      - zdmfdu*puu(jl,jk+1)
          zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk)                &
                      - zdmfdu*pvu(jl,jk+1)
          if ( pmfu(jl,jk)>0.0D0 ) then
            puu(jl,jk) = zmfuu(jl)*(1.0D0/pmfu(jl,jk))
            pvu(jl,jk) = zmfuv(jl)*(1.0D0/pmfu(jl,jk))
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
      if ( loflag(jl) .and. ktype(jl)==1 ) then
        zbuoyz = egrav*(ptu(jl,jk)-ptenh(jl,jk))/ptenh(jl,jk)         &
                 + egrav*vtmpc1*(pqu(jl,jk)-pqenh(jl,jk))             &
                 - egrav*plu(jl,jk)
        zbuoyz = max(zbuoyz,0.00D0)
        zdz = (pgeo(jl,jk-1)-pgeo(jl,jk))*regrav
        zdrodz = -log(pten(jl,jk-1)/pten(jl,jk))                    &
                 /zdz - egrav/(rgas*ptenh(jl,jk)                      &
                 *(1.0D0+vtmpc1*pqenh(jl,jk)))
        zbuoy(jl) = zbuoy(jl) + zbuoyz*zdz
        zoentr(jl,jk-1) = zbuoyz*0.50D0/(1.0D0+zbuoy(jl)) + zdrodz
        zoentr(jl,jk-1) = min(zoentr(jl,jk-1),centrmax)
        zoentr(jl,jk-1) = max(zoentr(jl,jk-1),0.0D0)
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
    if ( kctop(jl)==klevm1 ) ldcum(jl) = .false.
    kcbot(jl) = max(kcbot(jl),kctop(jl))
  end do
  do jl = 1 , kproma
    if ( ldcum(jl) ) then
      jk = kctop(jl) - 1
      zzdmf = cmfctop
      zdmfde(jl) = (1.0D0-zzdmf)*pmfu(jl,jk+1)
      plude(jl,jk) = zdmfde(jl)*plu(jl,jk+1)
      pqude(jl,jk) = zdmfde(jl)*pqu(jl,jk+1)
      pmfu(jl,jk) = pmfu(jl,jk+1) - zdmfde(jl)
      pdmfup(jl,jk) = 0.0D0
      pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                     *pmfu(jl,jk)
      pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
      pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
      if ( jk>=2 ) then
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pcpen , pdmfup , pgeo ,&
         pgeoh , plu , plude , pmfu , pmful , pmfuq , pmfus , pqen ,&
         pqenh , pqsen , pqte , pqu , pqude , pten , ptenh , ptu ,  &
         puen , puu , pven , pverv , pvu
  real(8) , dimension(kbdim) :: pentr , pmfub
  real(8) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtenh ,&
         pxtu
  intent (in) ldland , pcpcu , pxtenh
  intent (out) pqude
  intent (inout) kcbot , kctop , klab , ktype , ldcum , plu ,       &
                 plude , pmfu , pmfub , pmful , pmfuq , pmfus ,     &
                 pmfuxt , pqu , ptu , puu , pvu , pxtu
!
  integer :: icall , ik , jk , jl , jt
  logical , dimension(kbdim) :: loflag
  real(8) :: zbuo , zcons2 , zdlev , zdmfdu , zdmfeu , zdnoprc ,    &
             zfac , zlnew , zmfmax , zmftest , zmfulk , zmfuqk ,    &
             zmfusk , zmfuxtk , zprcon , zqeen , zqude , zscde ,    &
             zseen , ztglace , zxteen , zxtude , zz , zzdmf
  real(8) , dimension(kbdim) :: zdmfde , zdmfen , zmfuu , zmfuv ,   &
                                zpbase , zph , zqold
!
!----------------------------------------------------------------------
!*    1.           SPECIFY PARAMETERS
!     ------------------
!
  zcons2 = 1.0D0/(egrav*dt)
  ztglace = tzero - 13.0D0

! AMT NOTE!!! in the original scheme, this level which restricts rainfall 
! below a certain pressure (from the surface) is hard wired according to the 
! vertical resolution of the model - This has been fixed to 150 hPa.
! WARNING - this should be set outside of the routines.
!      if ( klev/=11 ) then
!        zdlev = 3.0D4
!      else if ( nn==21 ) then
!        zdlev = 1.5D4
!      else if ( nn==31 ) then
!        zdlev = 2.0D4 
!      else
!        zdlev = 3.0D4
!      end if
   zdlev = 1.5D4

!
!----------------------------------------------------------------------
!     2.           SET DEFAULT VALUES
!     ------------------
!
  do jl = 1 , kproma
    zmfuu(jl) = 0.0D0
    zmfuv(jl) = 0.0D0
    if ( .not.ldcum(jl) ) ktype(jl) = 0
  end do
  do jk = 1 , klev
    do jl = 1 , kproma
      plu(jl,jk) = 0.0D0
      pmfu(jl,jk) = 0.0D0
      pmfus(jl,jk) = 0.0D0
      pmfuq(jl,jk) = 0.0D0
      pmful(jl,jk) = 0.0D0
      plude(jl,jk) = 0.0D0
      pqude(jl,jk) = 0.0D0
      pdmfup(jl,jk) = 0.0D0
      if ( .not.ldcum(jl) .or. ktype(jl)==3 ) klab(jl,jk) = 0
      if ( .not.ldcum(jl) .and. paphp1(jl,jk)<4.D4 ) kctop0(jl) = jk
      if ( jk<kcbot(jl) ) klab(jl,jk) = 0
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        pmfuxt(jl,jk,jt) = 0.0D0
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
      pmfub(jl) = 0.0D0
      pqu(jl,klev) = 0.0D0
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
      if ( .not.ldcum(jl) ) pxtu(jl,klev,jt) = 0.0D0
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
    if ( lmfmid .and. ik<klevm1 .and. ik>nmctop )                   &
         call cubasmc(kproma,kbdim,klev,ik,klab,pten,pqen,pqsen,    &
         puen,pven,ktrac,pxten,pxtu,pmfuxt,pverv,pgeo,pgeoh,ldcum,  &
         ktype,pmfu,pmfub,pentr,kcbot,ptu,pqu,plu,puu,pvu,pmfus,    &
         pmfuq,pmful,pdmfup,zmfuu,pcpen,zmfuv)
!
    do jl = 1 , kproma
      if ( klab(jl,jk+1)==0 ) klab(jl,jk) = 0
      loflag(jl) = klab(jl,jk+1)>0
      zph(jl) = paphp1(jl,jk)
      if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
        zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
        if ( pmfub(jl)>zmfmax ) then
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
        if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
          zmfmax = (paphp1(jl,jk)-paphp1(jl,jk-1))*zcons2
          if ( pmfub(jl)>zmfmax ) then
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
      if ( ktype(jl)==3 .and. jk==kcbot(jl) ) then
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
        if ( jk<kcbot(jl) ) then
          zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zmfmax = min(zmftest,(paphp1(jl,jk)-paphp1(jl,jk-1))      &
                   *zcons2)
          zdmfen(jl) = max(zdmfen(jl)-max(zmftest-zmfmax,0.0D0),    &
                       0.0D0)
        end if
        zdmfde(jl) = min(zdmfde(jl),0.750D0*pmfu(jl,jk+1))
        pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
        zqeen = pqenh(jl,jk+1)*zdmfen(jl)
        zseen = (pcpcu(jl,jk+1)*ptenh(jl,jk+1)+pgeoh(jl,jk+1))      &
                *zdmfen(jl)
        zscde = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1))        &
                *zdmfde(jl)
        zqude = pqu(jl,jk+1)*zdmfde(jl)
        pqude(jl,jk) = zqude
        plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
        zmfusk = pmfus(jl,jk+1) + zseen - zscde
        zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
        zmfulk = pmful(jl,jk+1) - plude(jl,jk)
        plu(jl,jk) = zmfulk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
        pqu(jl,jk) = zmfuqk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
        ptu(jl,jk) = (zmfusk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))       &
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
          pxtu(jl,jk,jt) = zmfuxtk*(1.0D0/max(cmfcmin,pmfu(jl,jk)))
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
      if ( loflag(jl) .and. pqu(jl,jk)<zqold(jl) ) then
        klab(jl,jk) = 2
        plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
        zbuo = ptu(jl,jk)*(1.0D0+vtmpc1*pqu(jl,jk)-plu(jl,jk))      &
               - ptenh(jl,jk)*(1.0D0+vtmpc1*pqenh(jl,jk))
        if ( klab(jl,jk+1)==1 ) zbuo = zbuo + 0.50D0
        if ( zbuo>0.0D0 .and. pmfu(jl,jk)>=0.10D0*pmfub(jl) ) then
          kctop(jl) = jk
          ldcum(jl) = .true.
          zdnoprc = merge(zdlev,1.5D4,ldland(jl))
          zprcon = merge(0.0D0,cprcon,zpbase(jl)-paphp1(jl,jk)      &
                   <zdnoprc)
          zlnew = plu(jl,jk)                                        &
                  /(1.0D0+zprcon*(pgeoh(jl,jk)-pgeoh(jl,jk+1)))
          pdmfup(jl,jk) = max(0.0D0,(plu(jl,jk)-zlnew)*pmfu(jl,jk))
          plu(jl,jk) = zlnew
        else
          klab(jl,jk) = 0
          pmfu(jl,jk) = 0.0D0
        end if
      end if
    end do
    do jl = 1 , kproma
      if ( loflag(jl) ) then
        pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
        pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))       &
                       *pmfu(jl,jk)
        pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
      end if
    end do
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( loflag(jl) ) pmfuxt(jl,jk,jt) = pxtu(jl,jk,jt)         &
             *pmfu(jl,jk)
      end do
    end do
!
    if ( lmfdudv ) then
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          if ( ktype(jl)==1 .or. ktype(jl)==3 ) then
            zz = merge(3.0D0,2.0D0,zdmfen(jl)==0.0D0)
          else
            zz = merge(1.0D0,0.0D0,zdmfen(jl)==0.0D0)
          end if
          zdmfeu = zdmfen(jl) + zz*zdmfde(jl)
          zdmfdu = zdmfde(jl) + zz*zdmfde(jl)
          zdmfdu = min(zdmfdu,0.750D0*pmfu(jl,jk+1))
          zmfuu(jl) = zmfuu(jl) + zdmfeu*puen(jl,jk)                &
                      - zdmfdu*puu(jl,jk+1)
          zmfuv(jl) = zmfuv(jl) + zdmfeu*pven(jl,jk)                &
                      - zdmfdu*pvu(jl,jk+1)
          if ( pmfu(jl,jk)>0.0D0 ) then
            puu(jl,jk) = zmfuu(jl)*(1.0D0/pmfu(jl,jk))
            pvu(jl,jk) = zmfuv(jl)*(1.0D0/pmfu(jl,jk))
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
    if ( kctop(jl)==klevm1 ) ldcum(jl) = .false.
    kcbot(jl) = max(kcbot(jl),kctop(jl))
  end do
  do jl = 1 , kproma
    if ( ldcum(jl) ) then
      jk = kctop(jl) - 1
      zzdmf = cmfctop
      zdmfde(jl) = (1.0D0-zzdmf)*pmfu(jl,jk+1)
      plude(jl,jk) = zdmfde(jl)*plu(jl,jk+1)
      pqude(jl,jk) = zdmfde(jl)*pqu(jl,jk+1)
      pmfu(jl,jk) = pmfu(jl,jk+1) - zdmfde(jl)
      pdmfup(jl,jk) = 0.0D0
      pmfus(jl,jk) = (pcpcu(jl,jk)*ptu(jl,jk)+pgeoh(jl,jk))         &
                     *pmfu(jl,jk)
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
  real(8) , dimension(kbdim,klevp1) :: paph
  real(8) , dimension(kbdim,klev) :: pcpcu , pgeoh , plu , pqenh ,  &
         pqu , ptenh , ptu , puen , puu , pven , pvu
  intent (in) paph , pcpcu , pgeoh , pqenh , ptenh , puen , pven
  intent (inout) kcbot , klab , ldcum , plu , pqu , ptu , puu , pvu
!
! Local variables
!
  integer :: icall , ik , ikb , is , jk , jl
  logical , dimension(kbdim) :: loflag
  real(8) :: zbuo , zz
  real(8) , dimension(kbdim) :: zph , zqold
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
      is = is + merge(1,0,klab(jl,jk+1)==1)
      loflag(jl) = klab(jl,jk+1)==1
      zph(jl) = paph(jl,jk)
    end do
    if ( is/=0 ) then
      do jl = 1 , kproma
        if ( loflag(jl) ) then
          pqu(jl,jk) = pqu(jl,jk+1)
          ptu(jl,jk) = (pcpcu(jl,jk+1)*ptu(jl,jk+1)+pgeoh(jl,jk+1)  &
                       -pgeoh(jl,jk))/pcpcu(jl,jk)
          zbuo = ptu(jl,jk)*(1.0D0+vtmpc1*pqu(jl,jk)) - ptenh(jl,jk)&
                 *(1.0D0+vtmpc1*pqenh(jl,jk)) + 0.50D0
          if ( zbuo>0.0D0 ) klab(jl,jk) = 1
          zqold(jl) = pqu(jl,jk)
        end if
      end do
!
      ik = jk
      icall = 1
      call cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)
!
      do jl = 1 , kproma
        if ( loflag(jl) .and. pqu(jl,jk)<zqold(jl) ) then
          klab(jl,jk) = 2
          plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
          zbuo = ptu(jl,jk)*(1.0D0+vtmpc1*pqu(jl,jk)-plu(jl,jk))    &
                 - ptenh(jl,jk)*(1.0D0+vtmpc1*pqenh(jl,jk)) + 0.50D0
          if ( zbuo>0. ) then
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
          if ( jk>=kcbot(jl) ) then
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
        zz = 1.0D0/(paph(jl,klevp1)-paph(jl,ikb))
        puu(jl,klev) = puu(jl,klev)*zz
        pvu(jl,klev) = pvu(jl,klev)*zz
      else
        puu(jl,klev) = puen(jl,klevm1)
        pvu(jl,klev) = pven(jl,klevm1)
      end if
    end do
  end if
!
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
  real(8) , dimension(kbdim,klev) :: pcpen , pdmfup , pgeo , pgeoh ,&
         plu , pmfu , pmful , pmfuq , pmfus , pqen , pqsen , pqu ,  &
         pten , ptu , puen , puu , pven , pverv , pvu
  real(8) , dimension(kbdim) :: pentr , pmfub , pmfuu , pmfuv
  real(8) , dimension(kbdim,klev,ktrac) :: pmfuxt , pxten , pxtu
  intent (in) ldcum , pcpen , pgeo , pgeoh , pqen , pqsen , pten , &
              puen , pven , pverv , pxten
  intent (out) kcbot , ktype , pdmfup , pentr , plu , pmfu , pmful ,&
               pmfuq , pmfus , pmfuu , pmfuv , pmfuxt
  intent (inout) klab , pmfub , pqu , ptu , puu , pvu , pxtu
!
  integer :: jl , jt
  logical , dimension(kbdim) :: llo3
  real(8) :: zzzmb
!
!----------------------------------------------------------------------
!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!     -------------------------------------------
!
  do jl = 1 , kproma
    llo3(jl) = .false.
    if ( .not.ldcum(jl) .and. klab(jl,kk+1)==0 .and. pqen(jl,kk)    &
         >0.900D0*pqsen(jl,kk) ) then
      llo3(jl) = .true.
      ptu(jl,kk+1) = (pcpen(jl,kk)*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,&
                     kk+1))/pcpen(jl,kk)
      pqu(jl,kk+1) = pqen(jl,kk)
      plu(jl,kk+1) = 0.0D0
      zzzmb = max(cmfcmin,-pverv(jl,kk)*regrav)
      zzzmb = min(zzzmb,cmfcmax)
      pmfub(jl) = zzzmb
      pmfu(jl,kk+1) = pmfub(jl)
      pmfus(jl,kk+1) = pmfub(jl)                                    &
                       *(pcpen(jl,kk)*ptu(jl,kk+1)+pgeoh(jl,kk+1))
      pmfuq(jl,kk+1) = pmfub(jl)*pqu(jl,kk+1)
      pmful(jl,kk+1) = 0.0D0
      pdmfup(jl,kk+1) = 0.0D0
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
         pmfdq , pmfds , pqd , pqenh , ptd , ptenh , pud , puen ,   &
         pvd , pven
  real(8) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh
  real(8) , dimension(kbdim) :: prfl
  intent (in) lddraf , paphp1 , pcpcu , pgeoh , pqenh , ptenh , &
              puen , pven , pxtenh
  intent (out) pdmfdp
  intent (inout) pmfd , pmfdq , pmfds , pmfdxt , pqd , prfl , ptd , &
                 pud , pvd , pxtd
!
  integer :: icall , ik , is , itopde , jk , jl , jt
  logical :: llo1
  logical , dimension(kbdim) :: llo2
  real(8) :: zbuo , zdmfdp , zentr , zmfdqk , zmfdsk , zmfduk ,     &
             zmfdvk , zmfdxtk , zqdde , zqeen , zsdde , zseen ,     &
             zxtdde , zxteen
  real(8) , dimension(kbdim) :: zcond , zdmfde , zdmfen , zph
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
      llo2(jl) = lddraf(jl) .and. pmfd(jl,jk-1)<0.0D0
      is = is + merge(1,0,llo2(jl))
    end do
    if ( is/=0 ) then
      do jl = 1 , kproma
        if ( llo2(jl) ) then
          zentr = entrdd*pmfd(jl,jk-1)*rgas*ptenh(jl,jk-1)          &
                  /(egrav*paphp1(jl,jk-1))                            &
                  *(paphp1(jl,jk)-paphp1(jl,jk-1))
          zdmfen(jl) = zentr
          zdmfde(jl) = zentr
        end if
      end do
      itopde = klev - 2
      if ( jk>itopde ) then
        do jl = 1 , kproma
          if ( llo2(jl) ) then
            zdmfen(jl) = 0.0D0
            zdmfde(jl) = pmfd(jl,itopde)                            &
                         *(paphp1(jl,jk)-paphp1(jl,jk-1))           &
                         /(paphp1(jl,klevp1)-paphp1(jl,itopde))
          end if
        end do
      end if
!
      do jl = 1 , kproma
        if ( llo2(jl) ) then
          pmfd(jl,jk) = pmfd(jl,jk-1) + zdmfen(jl) - zdmfde(jl)
          zseen = (pcpcu(jl,jk-1)*ptenh(jl,jk-1)+pgeoh(jl,jk-1))    &
                  *zdmfen(jl)
          zqeen = pqenh(jl,jk-1)*zdmfen(jl)
          zsdde = (pcpcu(jl,jk-1)*ptd(jl,jk-1)+pgeoh(jl,jk-1))      &
                  *zdmfde(jl)
          zqdde = pqd(jl,jk-1)*zdmfde(jl)
          zmfdsk = pmfds(jl,jk-1) + zseen - zsdde
          zmfdqk = pmfdq(jl,jk-1) + zqeen - zqdde
          pqd(jl,jk) = zmfdqk*(1.0D0/min(-cmfcmin,pmfd(jl,jk)))
          ptd(jl,jk) = (zmfdsk*(1.0D0/min(-cmfcmin,pmfd(jl,jk)))    &
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
            pxtd(jl,jk,jt) = zmfdxtk*(1.0D0/min(-cmfcmin,pmfd(jl,jk)&
                             ))
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
          zbuo = ptd(jl,jk)*(1.0D0+vtmpc1*pqd(jl,jk)) - ptenh(jl,jk)&
                 *(1.0D0+vtmpc1*pqenh(jl,jk))
          llo1 = zbuo<0.0D0 .and.                                   &
                 (prfl(jl)-pmfd(jl,jk)*zcond(jl)>0.0D0)
          pmfd(jl,jk) = merge(pmfd(jl,jk),0.0D0,llo1)
          pmfds(jl,jk) = (pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))     &
                         *pmfd(jl,jk)
          pmfdq(jl,jk) = pqd(jl,jk)*pmfd(jl,jk)
          zdmfdp = -pmfd(jl,jk)*zcond(jl)
          pdmfdp(jl,jk-1) = zdmfdp
          prfl(jl) = prfl(jl) + zdmfdp
        end if
      end do
!
      do jt = 1 , ktrac
        do jl = 1 , kproma
          if ( llo2(jl) ) pmfdxt(jl,jk,jt) = pxtd(jl,jk,jt)         &
               *pmfd(jl,jk)
        end do
      end do
!
      if ( lmfdudv ) then
        do jl = 1 , kproma
          if ( llo2(jl) .and. pmfd(jl,jk)<0.0D0 ) then
            zmfduk = pmfd(jl,jk-1)*pud(jl,jk-1) + zdmfen(jl)        &
                     *puen(jl,jk-1) - zdmfde(jl)*pud(jl,jk-1)
            zmfdvk = pmfd(jl,jk-1)*pvd(jl,jk-1) + zdmfen(jl)        &
                     *pven(jl,jk-1) - zdmfde(jl)*pvd(jl,jk-1)
            pud(jl,jk) = zmfduk*(1.0D0/min(-cmfcmin,pmfd(jl,jk)))
            pvd(jl,jk) = zmfdvk*(1.0D0/min(-cmfcmin,pmfd(jl,jk)))
          end if
        end do
      end if
    end if
!
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pgeoh , pmfd ,&
         pmfdq , pmfds , pqd , pqenh , pqu , ptd , ptenh , ptu ,    &
         pud , puen , puu , pvd , pven , pvu
  real(8) , dimension(kbdim,klev,ktrac) :: pmfdxt , pxtd , pxtenh , &
         pxtu
  real(8) , dimension(kbdim) :: pmfub , prfl
  intent (in) kcbot , kctop , ldcum , paphp1 , pcpcu , pgeoh , &
              pmfub , pqenh , pqu , ptenh , ptu , puen , puu , &
              pven , pvu , pxtenh , pxtu
  intent (out) kdtop , pmfdq , pmfds , pmfdxt , pud , pvd
  intent (inout) lddraf , pdmfdp , pmfd , pqd , prfl , ptd , pxtd
!
  integer :: icall , ik , is , jk , jl , jt , ke
  logical , dimension(kbdim) :: llo2 , llo3
  real(dp) :: zbuo , zmftop , zqtest , zttest
  real(8) , dimension(kbdim) :: zcond , zph
  real(8) , dimension(kbdim,klev) :: zqenwb , ztenwb
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
        llo2(jl) = ldcum(jl) .and. prfl(jl)>0.0D0 .and.             &
                   .not.lddraf(jl) .and.                            &
                   (jk<kcbot(jl) .and. jk>kctop(jl))
        is = is + merge(1,0,llo2(jl))
      end do
      if ( is/=0 ) then
!
        ik = jk
        icall = 2
        call cuadjtq(kproma,kbdim,klev,ik,zph,ztenwb,zqenwb,llo2,   &
                     icall)
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
            zttest = 0.50D0*(ptu(jl,jk)+ztenwb(jl,jk))
            zqtest = 0.50D0*(pqu(jl,jk)+zqenwb(jl,jk))
            zbuo = zttest*(1.0D0+vtmpc1*zqtest) - ptenh(jl,jk)      &
                   *(1.0D0+vtmpc1*pqenh(jl,jk))
            zcond(jl) = pqenh(jl,jk) - zqenwb(jl,jk)
            zmftop = -cmfdeps*pmfub(jl)
            if ( zbuo<0.0D0 .and. prfl(jl)>10.0D0*zmftop*zcond(jl) )&
                 then
              llo3(jl) = .true.
              kdtop(jl) = jk
              lddraf(jl) = .true.
              ptd(jl,jk) = zttest
              pqd(jl,jk) = zqtest
              pmfd(jl,jk) = zmftop
              pmfds(jl,jk) = pmfd(jl,jk)                            &
                             *(pcpcu(jl,jk)*ptd(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk) = pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1) = -0.50D0*pmfd(jl,jk)*zcond(jl)
              prfl(jl) = prfl(jl) + pdmfdp(jl,jk-1)
            end if
          end if
        end do
!
        do jt = 1 , ktrac
          do jl = 1 , kproma
            if ( llo3(jl) ) then
              pxtd(jl,jk,jt) = 0.50D0*(pxtu(jl,jk,jt)+pxtenh(jl,jk, &
                               jt))
              pmfdxt(jl,jk,jt) = pmfd(jl,jk)*pxtd(jl,jk,jt)
            end if
          end do
        end do
!
        if ( lmfdudv ) then
          do jl = 1 , kproma
            if ( pmfd(jl,jk)<0.0D0 ) then
              pud(jl,jk) = 0.50D0*(puu(jl,jk)+puen(jl,jk-1))
              pvd(jl,jk) = 0.50D0*(pvu(jl,jk)+pven(jl,jk-1))
            end if
          end do
        end if
      end if
!
    end do
  end if
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: paprc , paprs , prfl , prsfc ,      &
                                psfl , pssfc
  real(8) , dimension(kbdim,klev) :: pcpen , pdmfdp , pdmfup ,      &
         pdpmel , plude , pmfdq , pmfds , pmful , pmfuq , pmfus ,   &
         pqte , pqtec , pqude , pten , ptte , pxtec
  real(8) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtte
  intent (in) ldcum , paphp1 , pcpen , pdmfdp , pdmfup , pdpmel ,   &
              plude , pmfdq , pmfds , pmfdxt , pmful , pmfuq ,      &
              pmfus , pmfuxt , pqude , prfl , psfl , pten
  intent (out) pqtec , prsfc , pssfc , pxtec
  intent (inout) paprc , paprs , pqte , ptte , pxtte
!
  integer :: jk , jl , jt
  logical :: llo1
  real(8) :: zalv , zdiagt , zdqdt , zdtdt , zdxtdt , zrcpm
  real(8) , dimension(kbdim) :: zmelt , zsheat
!
!----------------------------------------------------------------------
!*    1.0          SPECIFY PARAMETERS
!     ------------------
!
  zdiagt = dt
             ! reciprocal value of specific heat of moist air
!
!
!----------------------------------------------------------------------
!*    2.0          INCREMENTATION OF T AND Q TENDENCIES
!     ------------------------------------
!
  do jl = 1 , kproma
    zmelt(jl) = 0.0D0
    zsheat(jl) = 0.0D0
  end do
!
  do jk = ktopm2 , klev
!
    if ( jk<klev ) then
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          llo1 = (pten(jl,jk)-tzero)>0.0D0
          zalv = merge(wlhv,wlhs,llo1)
          zrcpm = 1.0D0/pcpen(jl,jk)
          zdtdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))             &
                  *zrcpm*(pmfus(jl,jk+1)-pmfus(jl,jk)+pmfds(jl,jk+1)&
                  -pmfds(jl,jk)-wlhf*pdpmel(jl,jk)                   &
                  -zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl,jk)   &
                  -(pdmfup(jl,jk)+pdmfdp(jl,jk))))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))             &
                  *(pmfuq(jl,jk+1)-pmfuq(jl,jk)+pmfdq(jl,jk+1)      &
                  -pmfdq(jl,jk)+pmful(jl,jk+1)-pmful(jl,jk)         &
                  -plude(jl,jk)-(pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                         *plude(jl,jk)
          pqtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
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
          llo1 = (pten(jl,jk)-tzero)>0.0D0
          zalv = merge(wlhv,wlhs,llo1)
          zrcpm = 1.0D0/pcpen(jl,jk)
          zdtdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))            &
                  *zrcpm*(pmfus(jl,jk)+pmfds(jl,jk)                 &
                  +wlhf*pdpmel(jl,jk)                                &
                  -zalv*(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)   &
                  +plude(jl,jk)))
          ptte(jl,jk) = ptte(jl,jk) + zdtdt
          zdqdt = -(egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))            &
                  *(pmfuq(jl,jk)+pmfdq(jl,jk)+plude(jl,jk)          &
                  +(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
          pqte(jl,jk) = pqte(jl,jk) + zdqdt
          pxtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
                         *plude(jl,jk)
          pqtec(jl,jk) = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))      &
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
    paprc(jl) = paprc(jl) + zdiagt*(prfl(jl)+psfl(jl))
    paprs(jl) = paprs(jl) + zdiagt*psfl(jl)
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pmfd , pmfu , pud , puen ,     &
         puu , pvd , pven , pvol , pvom , pvu
  intent (in) kcbot , ktype , ldcum , paphp1 , pmfd , pmfu , pud , &
              puen , puu , pvd , pven , pvu
  intent (inout) pvol , pvom
!
  integer :: ik , ikb , jk , jl
  real(8) :: zdudt , zdvdt , zzp
  real(8) , dimension(kbdim,klev) :: zmfdu , zmfdv , zmfuu , zmfuv
!
!----------------------------------------------------------------------
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!     ----------------------------------------------
!
  if ( ktopm2==1 ) then
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
      if ( ldcum(jl) .and. jk>kcbot(jl) ) then
        ikb = kcbot(jl)
        zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))                    &
              /(paphp1(jl,klevp1)-paphp1(jl,ikb)))
        zzp = merge(zzp**2,zzp,ktype(jl)==3)
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
    if ( jk<klev ) then
      do jl = 1 , kproma
        if ( ldcum(jl) ) then
          zdudt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))             &
                  *(zmfuu(jl,jk+1)-zmfuu(jl,jk)+zmfdu(jl,jk+1)      &
                  -zmfdu(jl,jk))
          zdvdt = (egrav/(paphp1(jl,jk+1)-paphp1(jl,jk)))             &
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
  real(8) , dimension(kbdim,klev) :: pgeoh , pmfu , podetr , pqenh ,&
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
  real(8) :: zarg , zdprho , zentest , zentr , zorgde , zpmid ,     &
             zrrho , ztmzk , zzmzk
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
    zrrho = (rgas*ptenh(jl,kk+1)*(1.0D0+vtmpc1*pqenh(jl,kk+1)))     &
            /paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
    zpmid = 0.50D0*(ppbase(jl)+paphp1(jl,kctop0(jl)))
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
    llo1 = kk<kcbot(jl) .and. ldcum(jl)
    pdmfde(jl) = merge(zentr,0.0D0,llo1)
    llo2 = llo1 .and. ktype(jl)==2 .and.                            &
           (ppbase(jl)-paphp1(jl,kk)<0.2D5 .or. paphp1(jl,kk)>zpmid)
    pdmfen(jl) = merge(zentr,0.0D0,llo2)
    iklwmin = max(klwmin(jl),kctop0(jl)+2)
    llo2 = llo1 .and. ktype(jl)==3 .and. kk>=iklwmin
    if ( llo2 ) pdmfen(jl) = zentr
    if ( llo2 .and. pqenh(jl,kk+1)>1.D-5 ) then
      pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
      zentest = max(pqte(jl,kk),0.0D0)/pqenh(jl,kk+1)
      zentest = min(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
      pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
    end if
    llo2 = llo1 .and. ktype(jl)==1 .and.                            &
           (kk>=iklwmin .or. paphp1(jl,kk)>zpmid)
    if ( llo2 ) pdmfen(jl) = zentr
!
!       organized detrainment, detrainment starts at khmin
!
    llo2 = llo1 .and. ktype(jl)==1
    ikb = kcbot(jl)
    podetr(jl,kk) = 0.0D0
    if ( llo2 .and. kk<=khmin(jl) .and. kk>=kctop0(jl) ) then
      ikt = kctop0(jl)
      ikh = khmin(jl)
      if ( ikh>ikt ) then
        zzmzk = -(pgeoh(jl,ikh)-pgeoh(jl,kk))*regrav
        ztmzk = -(pgeoh(jl,ikh)-pgeoh(jl,ikt))*regrav
        zarg = 3.14150D0*(zzmzk/ztmzk)*0.50D0
        zorgde = tan(zarg)*3.14150D0*0.50D0/ztmzk
        zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*(regrav*zrrho)
        podetr(jl,kk) = min(zorgde,centrmax)*pmfu(jl,kk+1)*zdprho
      end if
    end if
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim) :: pdmfde , pdmfen , pentr , ppbase
  real(8) , dimension(kbdim,klev) :: pmfu , pqenh , pqte , ptenh
  intent (in) kcbot , kctop0 , klwmin , ktype , ldcum , paphp1 , &
              pentr , pqenh , pqte , ptenh
  intent (out) pdmfde , pdmfen
  intent (inout) pmfu , ppbase
!
  integer :: iklwmin , jl
  logical :: llo1 , llo2
  real(8) :: zdprho , zentest , zentr , zpmid , zrrho
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
    zrrho = (rgas*ptenh(jl,kk+1)*(1.0D0+vtmpc1*pqenh(jl,kk+1)))     &
            /paphp1(jl,kk+1)
    zdprho = (paphp1(jl,kk+1)-paphp1(jl,kk))*regrav
    zpmid = 0.50D0*(ppbase(jl)+paphp1(jl,kctop0(jl)))
    zentr = pentr(jl)*pmfu(jl,kk+1)*zdprho*zrrho
    llo1 = kk<kcbot(jl) .and. ldcum(jl)
    pdmfde(jl) = merge(zentr,0.0D0,llo1)
    llo2 = llo1 .and. ktype(jl)==2 .and.                            &
           (ppbase(jl)-paphp1(jl,kk)<0.2D5 .or. paphp1(jl,kk)>zpmid)
    pdmfen(jl) = merge(zentr,0.0D0,llo2)
    iklwmin = max(klwmin(jl),kctop0(jl)+2)
    llo2 = llo1 .and. (ktype(jl)==1 .or. ktype(jl)==3) .and.        &
           (kk>=iklwmin .or. paphp1(jl,kk)>zpmid)
    if ( llo2 ) pdmfen(jl) = zentr
    if ( llo2 .and. pqenh(jl,kk+1)>1.D-5 ) then
      pmfu(jl,kk+1) = max(pmfu(jl,kk+1),cmfcmin)
      zentest = max(pqte(jl,kk),0.0D0)/pqenh(jl,kk+1)
      zentest = min(centrmax,zentest/(pmfu(jl,kk+1)*zrrho))
      pdmfen(jl) = zentr + zentest*pmfu(jl,kk+1)*zrrho*zdprho
    end if
  end do
!
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
  real(8) , dimension(kbdim,klevp1) :: paphp1
  real(8) , dimension(kbdim,klev) :: pcpcu , pdmfdp , pdmfup ,      &
         pdpmel , pgeoh , pmfd , pmfdq , pmfds , pmfu , pmful ,     &
         pmfuq , pmfus , pqen , pqenh , pqsen , pten , ptenh
  real(8) , dimension(kbdim,klev,ktrac) :: pmfdxt , pmfuxt , pxtenh
  real(8) , dimension(kbdim) :: prain , prfl , psfl
  intent (in) kcbot , kctop , kdtop , ldcum , paphp1 , pcpcu , pgeoh , &
              pqen , pqenh , pqsen , pten , ptenh , pxtenh
  intent (out) pdpmel
  intent (inout) ktype , lddraf , pdmfdp , pdmfup , pmfd , pmfdq , &
                 pmfds , pmfdxt , pmfu , pmful , pmfuq , pmfus ,   &
                 pmfuxt , prain , prfl , psfl
!
  integer :: ikb , jk , jl , jt
  real(8) :: zcons1 , zcons2 , zcucov , zdpevap , zdrfl , zfac ,    &
             zrfl , zrfln , zrmin , zrnew , zrsum , zsnmlt ,        &
             ztmelp2 , zzp
  real(8) , dimension(kbdim) :: zpsubcl
!
!*    SPECIFY CONSTANTS
!
  zcons1 = cpd/(wlhf*egrav*dt)
  zcons2 = 1.0D0/(egrav*dt)
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
    if ( .not.ldcum(jl) .or. kdtop(jl)<kctop(jl) ) lddraf(jl)       &
         = .false.
    if ( .not.ldcum(jl) ) ktype(jl) = 0
  end do
  ktopm2 = 1
  do jk = ktopm2 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) .and. jk>=kctop(jl)-1 ) then
        pmfus(jl,jk) = pmfus(jl,jk) - pmfu(jl,jk)                   &
                       *(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
        pmfuq(jl,jk) = pmfuq(jl,jk) - pmfu(jl,jk)*pqenh(jl,jk)
        if ( lddraf(jl) .and. jk>=kdtop(jl) ) then
          pmfds(jl,jk) = pmfds(jl,jk) - pmfd(jl,jk)                 &
                         *(pcpcu(jl,jk)*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfdq(jl,jk) = pmfdq(jl,jk) - pmfd(jl,jk)*pqenh(jl,jk)
        else
          pmfd(jl,jk) = 0.0D0
          pmfds(jl,jk) = 0.0D0
          pmfdq(jl,jk) = 0.0D0
          pdmfdp(jl,jk-1) = 0.0D0
        end if
      end if
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk>=kctop(jl)-1 ) then
          pmfuxt(jl,jk,jt) = pmfuxt(jl,jk,jt) - pmfu(jl,jk)         &
                             *pxtenh(jl,jk,jt)
          if ( lddraf(jl) .and. jk>=kdtop(jl) ) then
            pmfdxt(jl,jk,jt) = pmfdxt(jl,jk,jt) - pmfd(jl,jk)       &
                               *pxtenh(jl,jk,jt)
          else
            pmfdxt(jl,jk,jt) = 0.0D0
          end if
        else
          pmfuxt(jl,jk,jt) = 0.0D0
          pmfdxt(jl,jk,jt) = 0.0D0
        end if
      end do
    end do
!
  end do
  do jk = ktopm2 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) .and. jk>kcbot(jl) ) then
        ikb = kcbot(jl)
        zzp = ((paphp1(jl,klevp1)-paphp1(jl,jk))                    &
              /(paphp1(jl,klevp1)-paphp1(jl,ikb)))
        zzp = merge(zzp**2,zzp,ktype(jl)==3)
        pmfu(jl,jk) = pmfu(jl,ikb)*zzp
        pmfus(jl,jk) = pmfus(jl,ikb)*zzp
        pmfuq(jl,jk) = pmfuq(jl,ikb)*zzp
        pmful(jl,jk) = pmful(jl,ikb)*zzp
      end if
    end do
!
    do jt = 1 , ktrac
      do jl = 1 , kproma
        if ( ldcum(jl) .and. jk>kcbot(jl) ) then
          ikb = kcbot(jl)
          zzp = (paphp1(jl,klevp1)-paphp1(jl,jk))                   &
                /(paphp1(jl,klevp1)-paphp1(jl,ikb))
          zzp = merge(zzp**2,zzp,ktype(jl)==3)
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
    prfl(jl) = 0.0D0
    psfl(jl) = 0.0D0
    prain(jl) = 0.0D0
  end do
  do jk = ktopm2 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) ) then
        prain(jl) = prain(jl) + pdmfup(jl,jk)
        if ( pten(jl,jk)>tzero ) then
          prfl(jl) = prfl(jl) + pdmfup(jl,jk) + pdmfdp(jl,jk)
          if ( psfl(jl)>0.0D0 .and. pten(jl,jk)>ztmelp2 ) then
            zfac = zcons1*(1.0D0+vtmpc2*pqen(jl,jk))                &
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
    prfl(jl) = max(prfl(jl),0.0D0)
    psfl(jl) = max(psfl(jl),0.0D0)
    zpsubcl(jl) = prfl(jl) + psfl(jl)
  end do
  do jk = ktopm2 , klev
    do jl = 1 , kproma
      if ( ldcum(jl) .and. jk>=kcbot(jl) .and. zpsubcl(jl)>1.D-20 ) &
           then
        zrfl = zpsubcl(jl)
        zrnew = (max(0.0D0,sqrt(zrfl/zcucov)-cevapcu(jk)*(paphp1(jl,&
                jk+1)-paphp1(jl,jk))                                &
                *max(0.0D0,pqsen(jl,jk)-pqen(jl,jk))))**2*zcucov
        zrmin = zrfl - zcucov*max(0.0D0,0.80D0*pqsen(jl,jk)         &
                -pqen(jl,jk))*zcons2*(paphp1(jl,jk+1)-paphp1(jl,jk))
        zrnew = max(zrnew,zrmin)
        zrfln = max(zrnew,0.0D0)
        zdrfl = min(0.0D0,zrfln-zrfl)
        pdmfup(jl,jk) = pdmfup(jl,jk) + zdrfl
        zpsubcl(jl) = zrfln
      end if
    end do
  end do
  do jl = 1 , kproma
    zrsum = prfl(jl) + psfl(jl)
    zdpevap = zpsubcl(jl) - zrsum
    prfl(jl) = prfl(jl) + zdpevap*prfl(jl)*(1.0D0/max(1.D-20,zrsum))
    psfl(jl) = psfl(jl) + zdpevap*psfl(jl)*(1.0D0/max(1.D-20,zrsum))
  end do
!
  end subroutine cuflx
!
  subroutine cuadjtq(kproma,kbdim,klev,ik,zph,ptu,pqu,loflag,icall)
    implicit none
    integer , intent(in) :: kbdim , klev , kproma , ik , icall
    logical , intent(in) , dimension(kbdim) :: loflag
    real(8) , intent(in) , dimension(kbdim) :: zph
    real(8) , intent(in) , dimension(kbdim,klev) :: ptu , pqu
  end subroutine cuadjtq
!
end module mod_cu_tiedtke
