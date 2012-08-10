module mod_cu_tiedtke_38r2

  use mod_constants
  use mod_dynparam

  private

  integer , parameter :: n_vmass = 0        ! Using or not vector mass
  logical , parameter :: lphylin = .false.  ! Linearized physics is activated ?
  real(dp) , parameter :: rlpal1 = 0.15D0   ! Smoothing coefficient
  real(dp) , parameter :: rlpal2 = 20.0D0   ! Smoothing coefficient
  real(dp) , parameter :: rmfsoluv = 1.0D0  ! Mass flux solver for momentum
  real(dp) , parameter :: detrpen = 0.75D-4 ! Detrainment rate for penetrative
                                            !  convection
  real(dp) , parameter :: entrorg = 1.75D-3 ! Entrainment for positively
                                            !  buoyant convection 1/(m)
  real(dp) , parameter :: entrdd = 3.0D-4   ! Average entrainment rate
                                            !  for downdrafts
  real(dp) , parameter :: rmfcmax = 1.0D0   ! Maximum massflux value allowed
                                            !  for updrafts etc
  real(dp) , parameter :: rmfcmin = 1.0D-10 ! Minimum massflux value (safety)
  real(dp) , parameter :: rmfdeps = 0.30D0  ! Fractional massflux for
                                            !   downdrafts at lfs
  real(dp) , parameter :: zqmax = 0.5D0

  logical , public :: lmfmid    ! True if midlevel convection is switched on
  logical , public :: lmfdd     ! True if cumulus downdraft is switched on

  integer , parameter :: njkt2 = 2

  contains
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          P.BECHTOLD        E.C.M.W.F.     06/07
!          P.BECHTOLD        E.C.M.W.F.     03/12
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES
!          METHOD.
!          --------
!          TURBULENT ENTRAINMENT IS SIMULATED BY A CONSTANT
!          MULTIPLIED BY A VERTICAL SCALING FUNCTION
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           CURRENT LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    INPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    INPUT PARAMETERS (REAL):
!    *PQSEN*        SATURATION SPEC. HUMIDITY                   KG/KG
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA
!    *PGEOH*        PROVISIONAL GEOPOTENTIAL ON HALF LEVELS      PA
!    *PMFU*         MASSFLUX IN UPDRAFTS                        KG/(M2*S)
!    OUTPUT PARAMETERS (REAL):
!    *PDMFEN*       ENTRAINMENT RATE                            KG/(M2*S)
!    *PDMFDE*       DETRAINMENT RATE                            KG/(M2*S)
!          EXTERNALS
!          ---------
!          NONE
!----------------------------------------------------------------------
!
  subroutine cuentr(kidia,kfdia,klon,klev,kk,kcbot,ktype,ldcum,ldwork, &
                    pqsen,paph,pgeoh,pmfu,pdmfen,pdmfde)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer , intent(in) :: kk
    integer , dimension(klon) , intent(in) :: kcbot , ktype
    logical , dimension(klon) , intent(in) :: ldcum
    logical , intent(in) :: ldwork
    real(dp) , dimension(klon,klev) , intent(in) :: pqsen , pmfu
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph , pgeoh
    real(dp) , dimension(klon) , intent(out) :: pdmfen , pdmfde
    logical :: llo1
    integer :: jl
    real(dp) , dimension(klon) :: zentr
    real(dp) :: zdz , zmf
    !
    ! 1. Calculate entrainment and detrainment rates
    !
    if ( ldwork ) then
      do jl = kidia , kfdia
        pdmfen(jl) = d_zero
        pdmfde(jl) = d_zero
        zentr(jl) = d_zero
      end do
      !
      ! 1.1 Specify entrainment rates
      !
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zdz = (pgeoh(jl,kk)-pgeoh(jl,kk+1))*regrav
          zmf = pmfu(jl,kk+1)*zdz
          llo1 = kk < kcbot(jl)
          if ( llo1 ) then
            pdmfen(jl) = zentr(jl)*zmf
            pdmfde(jl) = detrpen*zmf
          end if
        end if
      end do
    end if
  end subroutine cuentr
!
!----------------------------------------------------------------------
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          PURPOSE
!          -------
!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)
!     PARAMETER     DESCRIPTION                                   UNITS 
!     ---------     -----------                                   ----- 
!     INPUT PARAMETERS (INTEGER): 
!    *KIDIA*        START POINT 
!    *KFDIA*        END POINT 
!    *KLON*         NUMBER OF GRID POINTS PER PACKET 
!    *KTDIA*        START OF THE VERTICAL LOOP 
!    *KLEV*         NUMBER OF LEVELS 
!    INPUT PARAMETERS (REAL): 
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K 
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG 
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG 
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    OUTPUT PARAMETERS (INTEGER):
!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY 
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CONDENSATION LEVEL
!    OUTPUT PARAMETERS (REAL):
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
!    *PQSENH*       ENV. SPEC. SATURATION HUMIDITY (T+1)
!                   ON HALF LEVELS                              KG/KG
!    *PTU*          TEMPERATURE IN UPDRAFTS                       K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                  KG/KG
!    *PTD*          TEMPERATURE IN DOWNDRAFTS                     K
!    *PQU*          SPEC. HUMIDITY IN DOWNDRAFTS                KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                       M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                     M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                     M/S
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!             05-02-11 : Optimisation (NJKT2) P. BECHTOLD
!----------------------------------------------------------------------
!
  subroutine cuinin(kidia,kfdia,klon,ktdia,klev,pten,pqen,pqsen,puen,pven, &
                    pvervel,pgeo,paph,pap,klwmin,klab,ptenh,pqenh,pqsenh,  &
                    pgeoh,ptu,pqu,ptd,pqd,puu,pvu,pud,pvd,plu)  
    implicit none

    integer , intent(in) :: klon 
    integer , intent(in) :: klev 
    integer , intent(in) :: kidia 
    integer , intent(in) :: kfdia 
    integer :: ktdia ! argument not used
    real(dp) , dimension(klon,klev) , intent(in) :: pten , pqen , pqsen , &
                     puen , pven , pvervel , pgeo , pap
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph , pgeoh
    integer , dimension(klon) , intent(out) :: klwmin
    integer , dimension(klon,klev) , intent(out) :: klab
    real(dp) , dimension(klon,klev) , intent(out) :: pqenh , ptu , pqu , &
                     ptd , pqd , puu , pvu , pud , pvd , plu
    real(dp) , dimension(klon,klev) , intent(inout) :: ptenh , pqsenh

    real(dp) , dimension(klon) :: zwmax , zph
    logical , dimension(klon) :: llflag
    integer :: icall , ik , jk , jl
    real(dp) :: zalfa , zzs
    !
    ! 1. specify large scale parameters at half levels
    !    adjust temperature fields if staticly unstable
    !    find level of maximum vertical velocity
    !
    zalfa = log(d_two)
    do jk = 2 , klev
      do jl = kidia , kfdia
        ptenh(jl,jk) = (max(cpd*pten(jl,jk-1)+pgeo(jl,jk-1), &
                        cpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))*rcpd  
        pqenh(jl,jk) = pqen(jl,jk-1)
        pqsenh(jl,jk) = pqsen(jl,jk-1)
        zph(jl) = paph(jl,jk)
        llflag(jl) = .true.
      end do
      if ( jk >= klev-1 .or. jk < njkt2 ) cycle
      ik = jk
      if ( lphylin ) then
        icall = 0
        call cuadjtqs(kidia,kfdia,klon,ktdia,klev,ik, &
                      zph,ptenh,pqsenh,llflag,icall)  
      else
        icall = 3
        call cuadjtq(kidia,kfdia,klon,ktdia,klev,ik, &
                     zph,ptenh,pqsenh,llflag,icall)  
      end if
      do jl = kidia , kfdia
        pqenh(jl,jk) = min(pqen(jl,jk-1),pqsen(jl,jk-1)) + &
                          (pqsenh(jl,jk)-pqsen(jl,jk-1))  
        pqenh(jl,jk) = max(pqenh(jl,jk),d_zero)
      end do
    end do
    do jl = kidia , kfdia
      ptenh(jl,klev) = (cpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))*rcpd
      pqenh(jl,klev) = pqen(jl,klev)
      ptenh(jl,1) = pten(jl,1)
      pqenh(jl,1) = pqen(jl,1)
      klwmin(jl) = klev
      zwmax(jl) = d_zero
    end do
    do jk = klev-1 , 2 , -1
      do jl = kidia , kfdia
        zzs = max(cpd*ptenh(jl,jk)  +pgeoh(jl,jk) , &
                  cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))  
        ptenh(jl,jk) = (zzs-pgeoh(jl,jk))*rcpd
      end do
    end do
    do jk = klev , 3 , -1
!dir$ ivdep
!ocl novrec
      do jl = kidia , kfdia
        if ( pvervel(jl,jk) < zwmax(jl) ) then
          zwmax(jl) = pvervel(jl,jk)
          klwmin(jl) = jk
        end if
      end do
    end do
    !
    ! 2.0 initialize values for updrafts and downdrafts
    !
    do jk = 1 , klev
      ik = jk-1
      if ( jk == 1 ) ik = 1
      do jl = kidia , kfdia
        ptu(jl,jk) = ptenh(jl,jk)
        ptd(jl,jk) = ptenh(jl,jk)
        pqu(jl,jk) = pqenh(jl,jk)
        pqd(jl,jk) = pqenh(jl,jk)
        plu(jl,jk) = d_zero
        puu(jl,jk) = puen(jl,ik)
        pud(jl,jk) = puen(jl,ik)
        pvu(jl,jk) = pven(jl,ik)
        pvd(jl,jk) = pven(jl,ik)
        klab(jl,jk) = 0
      end do
    end do
  end subroutine cuinin
!
!**   *CUADJTQS* - SIMPLIFIED VERSION OF MOIST ADJUSTMENT
!     J.F. MAHFOUF      ECMWF
!     PURPOSE.
!     --------
!     TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!     INTERFACE
!     ---------
!     THIS ROUTINE IS CALLED FROM SUBROUTINES:
!       *COND*
!       *CUBMADJ*
!       *CUBMD*
!       *CONDAD*
!       *CUBMADJAD*
!       *CUBMDAD*
!     INPUT ARE UNADJUSTED T AND Q VALUES,
!     IT RETURNS ADJUSTED VALUES OF T AND Q
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
!     INPUT PARAMETERS (REAL):
!    *PSP*          PRESSURE                                        PA
!     UPDATED PARAMETERS (REAL):
!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
!          MODIFICATIONS
!          -------------
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------
!
  subroutine cuadjtqs(kidia,kfdia,klon,ktdia,klev,kk,psp,pt,pq,ldflag,kcall)
    implicit none
    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer , intent(in) :: kk
    real(dp) , dimension(klon) , intent(in) :: psp
    real(dp) , dimension(klon,klev) , intent(inout) :: pt , pq
    logical , dimension(klon) , intent(in) :: ldflag
    integer , intent(in) :: kcall

    real(dp) , dimension(klon) :: z3es , z4es , z5alcp , zaldcp
    integer :: jl
    real(dp) :: zqmax , zqp , zcond , zcond1 , ztarg , zcor , &
                zqsat , zfoeew , z2s
    !
    ! calculate condensation and adjust t and q accordingly
    !
    ! ice-water thermodynamical functions
    !
    do jl = kidia , kfdia
      if ( pt(jl,kk) > tzero ) then
        z3es(jl) = c3les
        z4es(jl) = c4les
        z5alcp(jl) = c5alvcp
        zaldcp(jl) = wlhvocp
      else
        z3es(jl) = c3ies
        z4es(jl) = c4ies
        z5alcp(jl) = c5alscp
        zaldcp(jl) = wlhsocp
      end if
    end do

    if ( kcall == 1 ) then
!dir$    ivdep
!ocl novrec
      do jl = kidia , kfdia
        if ( ldflag(jl) ) then
          zqp = d_one/psp(jl)
          ztarg = pt(jl,kk)
          zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
          zqsat = zqp*zfoeew
          if ( zqsat > zqmax ) then
            zqsat = zqmax
          end if
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
          zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          zcond = max(zcond,d_zero)
          ! if ( dabs(zcond) > dlowval ) then
          pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond
          pq(jl,kk) = pq(jl,kk)-zcond
          ztarg = pt(jl,kk)
          zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
          zqsat = zqp*zfoeew
          if ( zqsat > zqmax ) then
            zqsat = zqmax
          end if
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          if ( dabs(zcond) < dlowval ) zcond1 = d_zero
          pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          ! end if
        end if
      end do
    end if

    if ( kcall == 2 ) then
!dir$    ivdep
!ocl novrec
      do jl = kidia , kfdia
        if ( ldflag(jl) ) then
          zqp = d_one/psp(jl)
          ztarg = pt(jl,kk)
          zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
          zqsat = zqp*zfoeew
          if ( zqsat > zqmax ) then
            zqsat = zqmax
          end if
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
          zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          zcond = min(zcond,d_zero)
          ! if ( dabs(zcond) > dlowval ) then
          pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond
          pq(jl,kk) = pq(jl,kk)-zcond
          ztarg = pt(jl,kk)
          zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
          zqsat = zqp*zfoeew
          if ( zqsat > zqmax ) then
            zqsat = zqmax
          end if
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          if ( dabs(zcond) < dlowval ) zcond1 = d_zero
          pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          ! end if
        end if
      end do
    end if

    if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
      do jl = kidia , kfdia
        zqp = d_one/psp(jl)
        ztarg = pt(jl,kk)
        zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
        zqsat = zqp*zfoeew
        if ( zqsat > zqmax ) then
          zqsat = zqmax
        end if
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
        pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk) = pq(jl,kk)-zcond1
        ztarg = pt(jl,kk)
        zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
        zqsat = zqp*zfoeew
        if ( zqsat > zqmax ) then
          zqsat = zqmax
        end if
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
        pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk) = pq(jl,kk)-zcond1
      end do
    end if

    if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
      do jl = kidia , kfdia
        zqp = d_one/psp(jl)
        ztarg = pt(jl,kk)
        zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
        zqsat = zqp*zfoeew
        if ( zqsat > zqmax ) then
          zqsat = zqmax
        end if
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
        zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
        pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond
        pq(jl,kk) = pq(jl,kk)-zcond
        ztarg = pt(jl,kk)
        zfoeew = c2es*exp(z3es(jl)*(ztarg-tzero)/(ztarg-z4es(jl)))
        zqsat = zqp*zfoeew
        if ( zqsat > zqmax ) then
          zqsat = zqmax
        end if
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        z2s = z5alcp(jl)/(ztarg-z4es(jl))**d_two
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
        pt(jl,kk) = pt(jl,kk)+zaldcp(jl)*zcond1
        pq(jl,kk) = pq(jl,kk)-zcond1
      end do
    end if
  end subroutine cuadjtqs
!
!-------------------------------------------------------------------------
!**   *CUADJTQ* - MOIST ADJUSTMENT
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          MODIFICATIONS
!          -------------
!          D.SALMOND         CRAY(UK))      12/8/91
!          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
!          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
!          D.SALMOND & M.HAMRUD ECMWF       99-06-04   Optimisation
!          J.HAGUE                          03-01-13   MASS Vector Functions
!          J.HAGUE                          03-07-07   More MASS V.F.
!        M.Hamrud              01-Oct-2003 CY28 Cleaning
!        J.Hague & D.Salmond   22-Nov-2005 Optimisations 
!          PURPOSE.
!          --------
!          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM SUBROUTINES:
!              *COND*     (T AND Q AT CONDENSATION LEVEL)
!              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
!              *CUASC*    (T AND Q AT CLOUD LEVELS)
!              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
!              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
!          INPUT ARE UNADJUSTED T AND Q VALUES,
!          IT RETURNS ADJUSTED VALUES OF T AND Q
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
!     INPUT PARAMETERS (REAL):
!    *PSP*          PRESSURE                                        PA
!     UPDATED PARAMETERS (REAL):
!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
!          EXTERNALS   
!          ---------
!          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
!          FOR CONDENSATION CALCULATIONS.
!          THE TABLES ARE INITIALISED IN *SUPHEC*.
!-------------------------------------------------------------------------
!
  subroutine cuadjtq(kidia,kfdia,klon,ktdia,klev,kk,psp,pt,pq,ldflag,kcall)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer , intent(in) :: kk
    real(dp) , dimension(klon) , intent(in) :: psp
    real(dp) , dimension(klon,klev) , intent(inout) :: pt , pq
    logical , dimension(klon) , intent(in) :: ldflag
    integer , intent(in) :: kcall
    integer :: jl , jlen

    real(dp) , dimension(kfdia-kidia+1) :: ztmp0
    real(dp) , dimension(kfdia-kidia+1) :: ztmp1
    real(dp) , dimension(kfdia-kidia+1) :: ztmp2
    real(dp) , dimension(kfdia-kidia+1) :: ztmp3
    real(dp) , dimension(kfdia-kidia+1) :: ztmp4
    real(dp) , dimension(kfdia-kidia+1) :: ztmp5
    real(dp) , dimension(kfdia-kidia+1) :: ztmp6

    real(dp) :: z1s , z2s , zcond , zcond1 , zcor , zfoeewi , zfoeewl , &
                zoealfa , zqsat , ztarg , zqp
    real(dp) :: zl , zi , zf

!dir$ vfunction exphf

    if ( .not. lphylin ) then
      if ( n_vmass > 0 ) then
        jlen = kfdia-kidia+1
      end if
      !
      ! calculate condensation and adjust t and q accordingly
      !
      if ( kcall == 1 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ! zqsat = foeewmcu(pt(jl,kk))*zqp
            ! foeewmcu(ptare) = c2es*(foealfcu(ptare)* &
            !                   exp(c3les*(ptare-tzero)/(ptare-c4les))+&
            !                 (d_one-foealfcu(ptare))* &
            !                   exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
            zl = d_one/(pt(jl,kk)-c4les)
            zi = d_one/(pt(jl,kk)-c4ies)
            zqsat = c2es *(foealfcu(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tzero)*zl)+&
                    (d_one-foealfcu(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tzero)*zi))
            zqsat = zqsat*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            ! foedemcu(ptare) = foealfcu(ptare)*c5alvcp*  &
            !                   (d_one/(ptare-c4les)**d_two)+ &
            !                   (d_one-foealfcu(ptare))*  &
            !                   c5alscp*(d_one/(ptare-c4ies)**d_two)
            zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                 (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
            zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            ! zcond = max(zcond,d_zero)
            if ( zcond > d_zero ) then
              pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ! zqsat = foeewmcu(pt(jl,kk))*zqp
              zl = d_one/(pt(jl,kk)-c4les)
              zi = d_one/(pt(jl,kk)-c4ies)
              zqsat = c2es*(foealfcu(pt(jl,kk)) * &
                      exp(c3les*(pt(jl,kk)-tzero)*zl) + &
                      (d_one-foealfcu(pt(jl,kk))) * &
                      exp(c3ies*(pt(jl,kk)-tzero)*zi))
              zqsat = zqsat*zqp
              zqsat = minj(zqmax,zqsat)
              zcor = d_one-retv*zqsat
              ! zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
              !          (zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
              zf = foealfcu(pt(jl,kk))*c5alvcp*zl**d_two + &
                   (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**d_two
              zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
              if ( dabs(zcond) < dlowval ) zcond1 = d_zero
              pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 2 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            zcond = min(zcond,d_zero)
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            if ( dabs(zcond) < dlowval ) zcond1 = min(zcond1,d_zero)
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end if
        end do
      end if

      if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(zqmax,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(zqmax,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end if
        end do
      end if

      if ( kcall == 5 ) then  ! same as 4 but with ldflag all true
!dir$    ivdep
!ocl novrec
        if ( n_vmass <= 0 )  then ! not using vector mass
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            zqsat = foeewm(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        else ! using vector vmass
          do jl = kidia , kfdia
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          call vrec(ztmp5,ztmp3,jlen)
          call vrec(ztmp6,ztmp4,jlen)
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = c2es*(foealfa(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(zqmax,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !         (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (d_one/(ptare-c4les)**d_two)+&
            !                 (d_one-foealfa(ptare))*c5alscp*&
            !                 (d_one/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
            zcond = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            ztmp0(jl-kidia+1) = zqp
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          call vrec(ztmp5,ztmp3,jlen)
          call vrec(ztmp6,ztmp4,jlen)
          do jl = kidia , kfdia
            zqp = ztmp0(jl-kidia+1)
            zqsat = c2es*(foealfa(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfa(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(zqmax,zqsat)
            zcor = d_one-retv*zqsat
            ! zcond1=(pq(jl,kk)*zcor**d_two-zqsat*zcor) / &
            !        (zcor**d_two+zqsat*foedem(pt(jl,kk)))
            ! foedem(ptare) = foealfa(ptare)*c5alvcp*&
            !                 (d_one/(ptare-c4les)**d_two)+&
            !                 (d_one-foealfa(ptare))*c5alscp*&
            !                 (d_one/(ptare-c4ies)**d_two)
            zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**d_two) + &
                 (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**d_two)
            zcond1 = (pq(jl,kk)*zcor**d_two-zqsat*zcor)/(zcor**d_two+zqsat*zf)
            pt(jl,kk) = pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        end if
      end if

      if ( kcall == 3 ) then
!dir$    ivdep !ocl novrec
        if ( n_vmass <= 0 ) then ! not using vector mass
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
            zqsat = foeewmcu(pt(jl,kk))*zqp
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv  *zqsat)
            zqsat = zqsat*zcor
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        else
          do jl = kidia , kfdia
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          do jl = kidia , kfdia
            zqp = d_one/psp(jl)
            zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(zqmax,zqsat)
            zcor = d_one-retv*zqsat
            zcond1 = (pq(jl,kk)*zcor**d_two - &
                     zqsat*zcor)/(zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
            ztmp0(jl-kidia+1) = zqp
            ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
            ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
            ztmp3(jl-kidia+1) = pt(jl,kk)-c4les
            ztmp4(jl-kidia+1) = pt(jl,kk)-c4ies
          end do
          call vdiv(ztmp5,ztmp1,ztmp3,jlen)
          call vdiv(ztmp6,ztmp2,ztmp4,jlen)
          call vexp(ztmp1,ztmp5,jlen)
          call vexp(ztmp2,ztmp6,jlen)
          do jl = kidia,kfdia
            zqp = ztmp0(jl-kidia+1)
            zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                    (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
            zqsat = minj(zqmax,zqsat)
            zcor = d_one-retv*zqsat
            zcond1 = (pq(jl,kk)*zcor**d_two - &
                     zqsat*zcor)/(zcor**d_two+zqsat*foedemcu(pt(jl,kk)))
            pt(jl,kk) = pt(jl,kk)+foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk)-zcond1
          end do
        end if
      end if
    else
    !
    !  calculate condensation and adjust t and q accordingly
    !
      if ( kcall == 1 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            zcond = max(zcond,d_zero)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                      ( zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
              zcor = d_one/(d_one-retv*zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                    (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 2 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            zcond = min(zcond,d_zero)
            if ( dabs(zcond) > dlowval ) then
              pt(jl,kk) = pt(jl,kk) + &
                         (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
              pq(jl,kk) = pq(jl,kk)-zcond
              ztarg = pt(jl,kk)
              zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
              z1s = tanh(rlpal2*(zqsat-zqmax))
              zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
              zcor = d_one/(d_one-retv  *zqsat)
              zqsat = zqsat*zcor
              z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                    (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
              zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
              pt(jl,kk) = pt(jl,kk) + &
                          (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
              pq(jl,kk) = pq(jl,kk)-zcond1
            end if
          end if
        end do
      end if

      if ( kcall == 0 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          ztarg = pt(jl,kk)
          zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
          zcor = d_one/(d_one-retv  *zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
          zcond1=(pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
          ztarg = pt(jl,kk)
          zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
          z1s = tanh(rlpal2*(zqsat-zqmax))
          zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
          pt(jl,kk) = pt(jl,kk)+(zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
          pq(jl,kk) = pq(jl,kk)-zcond1
        end do
      end if

      if ( kcall == 4 ) then
!dir$    ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldflag(jl) ) then
            zqp = d_one/psp(jl)
            ztarg = pt(jl,kk)
            zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa *c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                       (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond
            pq(jl,kk) = pq(jl,kk)-zcond
            ztarg = pt(jl,kk)
            zoealfa = d_half*(tanh(rlpal1*(ztarg-mpcrt))+d_one)
            zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
            zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
            zqsat = zqp*(zoealfa*zfoeewl+(d_one-zoealfa)*zfoeewi)
            z1s = tanh(rlpal2*(zqsat-zqmax))
            zqsat = d_half*((d_one-z1s)*zqsat+(d_one+z1s)*zqmax)
            zqsat = min(zqmax,zqsat)
            zcor = d_one/(d_one-retv*zqsat)
            zqsat = zqsat*zcor
            z2s = zoealfa*c5alvcp*(d_one/(ztarg-c4les)**d_two) + &
                  (d_one-zoealfa)*c5alscp*(d_one/(ztarg-c4ies)**d_two)
            zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*z2s)
            pt(jl,kk) = pt(jl,kk) + &
                   (zoealfa*wlhvocp+(d_one-zoealfa)*wlhsocp)*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end do
      end if
    end if
  end subroutine cuadjtq
!
!**** *CUDUDV* - UPDATES U AND V TENDENCIES,
!                DOES GLOBAL DIAGNOSTIC OF DISSIPATION
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!          P.BECHTOLD        E.C.M.W.F.    11/02/05 IMPLICIT SOLVER
!**   INTERFACE.
!     ----------
!          *CUDUDV* IS CALLED FROM *CUMASTR*
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    INPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    INPUT PARAMETERS (REAL):
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                       M/S
!    UPDATED PARAMETERS (REAL):
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!            METHOD
!            -------
!       EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!       DEPENDING ON VALUE OF RMFSOLUV:
!       0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT
!
!       FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!       FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                              TO CORRECT TENDENCIES BELOW CLOUD BASE
!
!            EXTERNALS
!            ---------
!            CUBIDIAG
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------
!
  subroutine cududv(kidia,kfdia,klon,ktdia,klev,ktopm2,ktype,kcbot,kctop, &
                    ldcum,ptsphy,paph,puen,pven,pmfu,pmfd,puu,pud,pvu,    &
                    pvd,ptenu,ptenv)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer , intent(in) :: ktdia 
    integer , intent(in) :: ktopm2
    integer , dimension(klon) , intent(in) :: ktype , kcbot , kctop
    logical , dimension(klon) , intent(in) :: ldcum
    real(dp) , intent(in) :: ptsphy
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph
    real(dp) , dimension(klon,klev) , intent(in) :: puen , pven , &
             pmfu , pmfd , puu , pud , pvu , pvd
    real(dp) , dimension(klon,klev) , intent(inout) :: ptenu , ptenv

    real(dp) , dimension(klon,klev) :: zuen , zven , zmfuu , zmfdu , &
             zmfuv , zmfdv , zdudt , zdvdt , zdp , zb ,  zr1 ,  zr2
    logical , dimension(klon,klev) :: llcumbas
    integer :: ik , ikb , jk , jl
    real(dp) :: zzp , ztsphy
    real(dp) , parameter :: zimp = d_one-rmfsoluv

    ztsphy = d_one/ptsphy

    do jk = 1 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zuen(jl,jk) = puen(jl,jk)
          zven(jl,jk) = pven(jl,jk)
          zdp(jl,jk) = egrav/(paph(jl,jk+0)-paph(jl,jk))
        end if
      end do
    end do
    !
    ! 1.0 Calculate fluxes and update u and v tendencies
    !
    do jk = ktopm2 , klev
      ik = jk-1
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-zimp*zuen(jl,ik))
          zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-zimp*zven(jl,ik))
          zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-zimp*zuen(jl,ik))
          zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-zimp*zven(jl,ik))
        end if
      end do
    end do
    !
    ! Linear fluxes below cloud
    !
    if ( dabs(rmfsoluv) < dlowval ) then
      do jk = ktopm2 , klev
!dir$ ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldcum(jl) .and. jk > kcbot(jl) ) then
            ikb = kcbot(jl)
            zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            if ( ktype(jl) == 3 ) then
              zzp = zzp*zzp
            end if
            zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
            zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
            zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
            zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
          end if
        end do
      end do
    end if
    !
    ! 1.2 Compute tendencies
    !
    do jk = ktopm2 , klev
      if ( jk < klev ) then
        ik = jk+1
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = zdp(jl,jk) * &
                      (zmfuu(jl,ik)-zmfuu(jl,jk)+zmfdu(jl,ik)-zmfdu(jl,jk))
            zdvdt(jl,jk) = zdp(jl,jk) * &
                      (zmfuv(jl,ik)-zmfuv(jl,jk)+zmfdv(jl,ik)-zmfdv(jl,jk))
          end if
        end do
      else
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = -zdp(jl,jk)*(zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt(jl,jk) = -zdp(jl,jk)*(zmfuv(jl,jk)+zmfdv(jl,jk))
          end if
        end do
      end if
    end do

    if ( dabs(rmfsoluv) < dlowval ) then
      !
      ! 1.3 Update tendencies
      !
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            ptenu(jl,jk) = ptenu(jl,jk)+zdudt(jl,jk)
            ptenv(jl,jk) = ptenv(jl,jk)+zdvdt(jl,jk)
          end if
        end do
      end do
    else
      !
      ! 1.6 Implicit solution
      !
      ! Fill bi-diagonal matrix vectors a=k-1, b=k;
      ! reuse zmfuu=a and zb=b;
      ! zdudt and zdvdt correspond to the rhs ("constants") of the equation
      ! the solution is in zr1 and zr2
      !
      llcumbas(:,:) = .false.
      zb(:,:) = d_one
      zmfuu(:,:) = d_zero
      !
      ! Fill vectors a, b and rhs
      !
      do jk = ktopm2 , klev
        ik = jk+1
        do jl = kidia , kfdia
          llcumbas(jl,jk) = ldcum(jl) .and. jk >= kctop(jl)-1
          if ( llcumbas(jl,jk) ) then
            zzp = rmfsoluv*zdp(jl,jk)*ptsphy
            zmfuu(jl,jk) = -zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdudt(jl,jk) = zdudt(jl,jk)*ptsphy+zuen(jl,jk)
            zdvdt(jl,jk) = zdvdt(jl,jk)*ptsphy+zven(jl,jk)
            ! zdudt(jl,jk) = (ptenu(jl,jk)+zdudt(jl,jk))*ptsphy+zuen(jl,jk)
            ! zdvdt(jl,jk) = (ptenv(jl,jk)+zdvdt(jl,jk))*ptsphy+zven(jl,jk)
            if ( jk<klev ) then
              zb(jl,jk) = d_one+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            else
              zb(jl,jk) = d_one
            end if
          end if
        end do
      end do
      !
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfuu,zb,zdudt,zr1)
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfuu,zb,zdvdt,zr2)
      !
      ! Compute tendencies
      !
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( llcumbas(jl,jk) ) then
            ptenu(jl,jk) = ptenu(jl,jk)+(zr1(jl,jk)-zuen(jl,jk))*ztsphy
            ptenv(jl,jk) = ptenv(jl,jk)+(zr2(jl,jk)-zven(jl,jk))*ztsphy
            ! ptenu(jl,jk) = (zr1(jl,jk)-zuen(jl,jk))*ztsphy
            ! ptenv(jl,jk) = (zr2(jl,jk)-zven(jl,jk))*ztsphy
          end if
        end do
      end do
    end if
  end subroutine cududv
!
!          P. Bechtold         E.C.M.W.F.     07/03
!          PURPOSE.
!          --------
!          SOLVES BIDIAGONAL SYSTEM
!          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
!          IT RETURNS UPDATED VALUE OF QUANTITY
!          METHOD.
!          --------
!          NUMERICAL RECIPES (Cambridge Press)
!          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
!          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
!          M  x  U  = R
!          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C
!          (  B(kctop-1)    C(kctop-1)    0          0        )
!          (  A(kctop)      B(kctop)      C(kctop)   0        )
!          (  0             A(jk)         B(jk)      C(jk)    )
!          (  0             0             A(klev)    B(klev)  )
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!    INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVELS
!    INPUT PARAMETERS (REAL):
!    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS
!    *PR*           RHS VECTOR CONTAINING "CONSTANTS"
!    OUTPUT PARAMETERS (REAL):
!    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY
!          EXTERNALS
!          ---------
!          NONE
!----------------------------------------------------------------------
!
  subroutine cubidiag(kidia,kfdia,klon,klev,kctop,ld_lcumask,pa,pb,pr,pu)
    implicit none
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , dimension(klon) , intent(in) :: kctop
    logical , dimension(klon,klev) , intent(in) :: ld_lcumask
    real(dp) , dimension(klon,klev) , intent(in) :: pa , pb , pr
    real(dp) , dimension(klon,klev) , intent(out) :: pu
    integer :: jk , jl
    real(dp) :: zbet
    real(dp) , parameter :: eps = 1.0D-35

    pu(:,:) = d_zero
    !
    ! Forward substitution
    !
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ld_lcumask(jl,jk) ) then
          if ( jk == kctop(jl)-1 ) then
            zbet = d_one/(pb(jl,jk)+eps)
            pu(jl,jk) = pr(jl,jk) * zbet
          else if ( jk > kctop(jl)-1 ) then
            zbet = d_one/(pb(jl,jk)+eps)
            pu(jl,jk) = (pr(jl,jk)-pa(jl,jk)*pu(jl,jk-1))*zbet
          end if
        end if
      end do
    end do
  end subroutine cubidiag
!
!-------------------------------------------------------------------------
!
! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION
!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96
!       Modified J. HAGUE          13/01/03 MASS Vector Functions       
!       PURPOSE.
!       --------
!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
!       AND LIQUID WATER CONTENT  
!       INTERFACE
!       ---------
!       THIS ROUTINE IS CALLED FROM *CALLPAR*.
!       PARAMETER     DESCRIPTION                                 UNITS
!       ---------     -----------                                 -----
!       INPUT PARAMETERS (INTEGER):
!      *KIDIA*        START POINT
!      *KFDIA*        END POINT
!      *KLON*         NUMBER OF GRID POINTS PER PACKET
!      *KTDIA*        START OF THE VERTICAL LOOP
!      *KLEV*         NUMBER OF LEVELS
!       INPUT PARAMETERS (REAL):
!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
!      *PT*            TEMPERATURE AT T-DT                          K
!       INPUT PARAMETERS (INTEGER):
!      *KFLAG*         FLAG TO DETECT CALL FROM
!                      CONVECTION  KFLAG=1
!                      OTHER       KFLAG=2
!       OUTPUT PARAMETER (REAL):
!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG
!-------------------------------------------------------------------------
!
  subroutine satur(kidia,kfdia,klon,ktdia,klev,paprsf,pt,pqsat,kflag)  
    implicit none

    integer , intent(in) :: klon 
    integer , intent(in) :: klev 
    integer , intent(in) :: kidia 
    integer , intent(in) :: kfdia 
    integer , intent(in) :: ktdia 
    real(dp) , dimension(klon,klev) , intent(in) :: paprsf , pt
    real(dp) , dimension(klon,klev) , intent(out) :: pqsat
    integer , intent(in) :: kflag 
    integer :: jk , jl , jlen

    real(dp) :: zcor , zew , zfoeew , zqmax , zqs , ztarg
    real(dp) :: zalfa , zfoeewl , zfoeewi
    real(dp) , dimension(kidia:kfdia) :: z_exparg1
    real(dp) , dimension(kidia:kfdia) :: z_exparg2
    real(dp) , dimension(kidia:kfdia) :: z_expout1
    real(dp) , dimension(kidia:kfdia) :: z_expout2

!dir$ vfunction exphf

    !
    ! calculate saturation specific humidity
    !
    if ( lphylin ) then
      do jk = ktdia , klev
        do jl = kidia , kfdia
          ztarg = pt(jl,jk)
          zalfa = foealfa(ztarg)
          zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
          zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
          zfoeew = zalfa*zfoeewl+(d_one-zalfa)*zfoeewi
          zqs = zfoeew/paprsf(jl,jk)
          if ( zqs > zqmax ) then
            zqs = zqmax
          end if
          zcor = d_one/(d_one-retv*zqs)
          pqsat(jl,jk) = zqs*zcor
        end do
      end do
    else
      if ( n_vmass <= 0 ) then ! not using vector mass
        do jk = ktdia , klev
          do jl = kidia , kfdia
            if ( kflag == 1 ) then
              zew = foeewmcu(pt(jl,jk))
            else
              zew = foeewm(pt(jl,jk))
            end if
            zqs = zew/paprsf(jl,jk)
            zqs = min(zqmax,zqs)
            zcor = d_one/(d_one-retv*zqs)
            pqsat(jl,jk) = zqs*zcor
          end do
        end do
      else ! using vector mass
        jlen = kfdia-kidia+1
        do jk = ktdia , klev
          do jl = kidia , kfdia
            z_exparg1(jl) = foeles_v(pt(jl,jk))
            z_exparg2(jl) = foeies_v(pt(jl,jk))
          end do
          call vexp(z_expout1,z_exparg1,jlen)
          call vexp(z_expout2,z_exparg2,jlen)
          do jl = kidia , kfdia
            if ( kflag == 1 ) then
              zew = foeewmcu_v( pt(jl,jk),z_expout1(jl),z_expout2(jl) )
            else
              zew = foeewm_v( pt(jl,jk),z_expout1(jl),z_expout2(jl) )
            end if
            ! zqs = zew/paprsf(jl,jk)
            ! zqs = min(zqmax,zqs)
            !! zcor = _one_/(_one_-retv*zqs)
            ! pqsat(jl,jk) = zqs/(_one_-retv*zqs)
            zqs = min(zqmax*paprsf(jl,jk),zew)
            pqsat(jl,jk) = zqs/(paprsf(jl,jk)-retv*zqs)
          end do
        end do
      end if
    end if
  end subroutine satur
!
!----------------------------------------------------------------------
!
!          M.TIEDTKE         E.C.M.W.F.     12/89
!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES CLOUD BASE VALUES
!          FOR MIDLEVEL CONVECTION
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          IT RETURNS CLOUDBASE VALUES FOR MIDLEVEL CONVECTION
!          METHOD.
!          --------
!          S. TIEDTKE (1989)
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           ACTUAL LEVEL
!    INPUT PARAMETERS (REAL):
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PLRAIN*       RAIN WATER CONTENT IN UPDRAFTS                KG/KG
!    INPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    UPDATED PARAMETERS (INTEGER):
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS
!    *KCBOT*        CLOUD BASE LEVEL
!    OUTPUT PARAMETERS (REAL):
!    *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE            KG/(M2*S)
!    *PENTR*        FRACTIONAL MASS ENTRAINMENT RATE               1/M
!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!          EXTERNALS
!          ---------
!          NONE
!----------------------------------------------------------------------
!
  subroutine cubasmcn(kidia,kfdia,klon,ktdia,klev,kk,pten,pqen,pqsen, &
                      pvervel,pgeo,pgeoh,ldcum,ktype,klab,kcbot,pmfu, &
                      pmfub,plrain,ptu,pqu,plu,pmfus,pmfuq,pmful,pdmfup)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer ,  intent(in) :: kk
    real(dp) , dimension(klon,klev) , intent(in) :: pten , pqen , &
                 pqsen , pvervel , pgeo
    real(dp) , dimension(klon,klev+1) , intent(in) :: pgeoh
    logical , dimension(klon) , intent(in) :: ldcum
    integer , dimension(klon) , intent(out) :: ktype , kcbot
    integer , dimension(klon,klev) , intent(inout) :: klab
    real(dp) , dimension(klon,klev) , intent(out) :: pmfu , plrain , ptu , &
                pqu , plu , pmfus , pmfuq , pmful , pdmfup
    real(dp) , dimension(klon) , intent(out) :: pmfub
    integer :: jl
    real(dp) :: zzzmb

    !
    ! 1. Calculate entrainment and detrainment rates
    !
!dir$ ivdep
!ocl novrec
    do jl = kidia , kfdia
      if ( .not. ldcum(jl) .and. klab(jl,kk+1) == 0 ) then
        if ( lmfmid .and. pgeo(jl,kk) > 5000.0D0 .and. &
                          pgeo(jl,kk) < 1.D5 .and.     &
                          pqen(jl,kk) > 0.80D0*pqsen(jl,kk) ) then  
          ptu(jl,kk+1) = (cpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))/cpd
          pqu(jl,kk+1) = pqen(jl,kk)
          plu(jl,kk+1) = d_zero
          zzzmb = max(rmfcmin,-pvervel(jl,kk)/egrav)
          zzzmb = min(zzzmb,rmfcmax)
          pmfub(jl) = zzzmb
          pmfu(jl,kk+1) = pmfub(jl)
          pmfus(jl,kk+1) = pmfub(jl)*(cpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
          pmfuq(jl,kk+1) = pmfub(jl)*pqu(jl,kk+1)
          pmful(jl,kk+1) = d_zero
          pdmfup(jl,kk+1) = d_zero
          plrain(jl,kk+1) = d_zero
          kcbot(jl) = kk
          klab(jl,kk+1) = 1
          ktype(jl) = 3
        end if
      end if
    end do
  end subroutine cubasmcn
!
!----------------------------------------------------------------------
!
!          THIS ROUTINE CALCULATES LEVEL OF FREE SINKING FOR
!          CUMULUS DOWNDRAFTS AND SPECIFIES T,Q,U AND V VALUES
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!          PURPOSE.
!          --------
!          TO PRODUCE LFS-VALUES FOR CUMULUS DOWNDRAFTS
!          FOR MASSFLUX CUMULUS PARAMETERIZATION
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONMENTAL VALUES OF T,Q,U,V,P,PHI
!          AND UPDRAFT VALUES T,Q,U AND V AND ALSO
!          CLOUD BASE MASSFLUX AND CU-PRECIPITATION RATE.
!          IT RETURNS T,Q,U AND V VALUES AND MASSFLUX AT LFS.
!          METHOD.
!          CHECK FOR NEGATIVE BUOYANCY OF AIR OF EQUAL PARTS OF
!          MOIST ENVIRONMENTAL AIR AND CLOUD AIR.
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    INPUT PARAMETERS (LOGICAL):
!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    INPUT PARAMETERS (REAL):
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                        M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                        M/S
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
!    UPDATED PARAMETERS (REAL):
!    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)
!    OUTPUT PARAMETERS (REAL):
!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
!    OUTPUT PARAMETERS (INTEGER):
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS
!    OUTPUT PARAMETERS (LOGICAL):
!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR CALCULATING WET BULB T AND Q AT LFS
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             99-06-04 : Optimisation        D.SALMOND
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Lopez      20-Jun-2007 CY32R2 Bug correction in latent heat
!                                         when LPHYLIN=T.
!----------------------------------------------------------------------
!
  subroutine cudlfsn(kidia,kfdia,klon,ktdia,klev,kcbot,kctop,ldland,ldcum, &
                     ptenh,pqenh,puen,pven,pten,pqsen,pgeo,pgeoh,paph,ptu, &
                     pqu,plu,puu,pvu,pmfub,prfl,ptd,pqd,pmfd,pmfds,pmfdq,  &
                     pdmfdp,kdtop,lddraf)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer , dimension(klon) :: kcbot , kctop  ! argument not used
    logical , dimension(klon) :: ldland , ldcum ! argument not used
    real(dp) , dimension(klon,klev) , intent(in) :: ptenh , pqenh ,  &
               puen , pven , pten , pqsen , pgeo , ptu , pqu , plu , &
               puu , pvu
    real(dp) , dimension(klon,klev+1) , intent(in) :: pgeoh , paph
    real(dp) , dimension(klon) , intent(in) :: pmfub
    real(dp) , dimension(klon) , intent(inout) :: prfl
    real(dp) , dimension(klon,klev) , intent(inout) :: pmfd
    real(dp) , dimension(klon,klev) , intent(out) :: ptd , pqd , pmfds , &
               pmfdq , pdmfdp
    integer , dimension(klon) , intent(out) :: kdtop
    logical , dimension(klon) , intent(out) :: lddraf
    integer , dimension(klon) :: ikhsmin
    real(dp) , dimension(klon,klev) :: ztenwb , zqenwb
    real(dp) , dimension(klon) :: zcond , zph , zhsmin
    logical , dimension(klon) :: llo2

    integer :: icall , ik , ike , is , jk , jl
    real(dp) :: zbuo , zhsk , zmftop , zoealfa , zoelhm , zqtest , &
                ztarg , zttest
    !
    ! 1. Set default values for downdrafts
    !
    do jl = kidia , kfdia
      lddraf(jl) = .false.
      kdtop(jl) = klev+1
      ikhsmin(jl) = klev+1
      zhsmin(jl) = 1.0D8
    end do

    if ( lmfdd ) then
      !
      ! 2. Determine level of free sinking:
      !    downdrafts shall start at model level of minimum
      !    of saturation moist static energy or below respectively
      !    For every point and proceed as follows:
      !      (1) determine level of minimum of hs
      !      (2) determine wet bulb environmental t and q
      !      (3) do mixing with cumulus cloud air
      !      (4) check for negative buoyancy
      !      (5) if buoyancy>0 repeat (2) to (4) for next level below
      !    The assumption is that air of downdrafts is mixture
      !    of 50% cloud air + 50% environmental air at wet bulb
      !    temperature (i.e. which became saturated due to
      !    evaporation of rain and cloud water)
      !
      do jk = 3 , klev-2
        if ( lphylin ) then
          do jl = kidia , kfdia
            ztarg = pten(jl,jk)
            zoealfa = 0.545D0*(tanh(0.17D0*(ztarg-mpcrt))+d_one)
            zoelhm  = zoealfa*wlhv+(d_one-zoealfa)*wlhs
            zhsk = cpd*pten(jl,jk)+pgeo(jl,jk)+zoelhm*pqsen(jl,jk)
            if ( zhsk < zhsmin(jl) ) then
              zhsmin(jl) = zhsk
              ikhsmin(jl) = jk
            end if
          end do
        else
          do jl = kidia , kfdia
            zhsk=cpd*pten(jl,jk)+pgeo(jl,jk)+foelhmcu(pten(jl,jk))*pqsen(jl,jk)
            if ( zhsk < zhsmin(jl) ) then
              zhsmin(jl) = zhsk
              ikhsmin(jl) = jk
            end if
          end do
        end if
      end do
      ike = klev-3
      do jk = 3 , ike
        !
        ! 2.1 Calculate wet-bulb temperature and moisture
        !     for environmental air in *cuadjtq*
        !
        is = 0
        do jl = kidia , kfdia
          ztenwb(jl,jk) = ptenh(jl,jk)
          zqenwb(jl,jk) = pqenh(jl,jk)
          zph(jl) = paph(jl,jk)
          llo2(jl) = ldcum(jl) .and. prfl(jl) > d_zero .and. &
                     .not. lddraf(jl) .and. &
                     (jk < kcbot(jl) .and. jk > kctop(jl) ) .and. &
                     jk >= ikhsmin(jl)
          if ( llo2(jl) )then
            is = is+1
          end if
        end do
        if ( is == 0 ) cycle
        ik = jk
        icall = 2
        call cuadjtq(kidia,kfdia,klon,ktdia,klev,ik, &
                     zph,ztenwb,zqenwb,llo2,icall)
        !
        ! 2.2 Do mixing of cumulus and environmental air
        !     and check for negative buoyancy.
        !     Then set values for downdraft at lfs.
        !
!dir$ ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( llo2(jl) ) then
            zttest = d_half*(ptu(jl,jk)+ztenwb(jl,jk))
            zqtest = d_half*(pqu(jl,jk)+zqenwb(jl,jk))
            zbuo = zttest*(d_one+retv*zqtest) - &
                   ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
            zcond(jl) = pqenh(jl,jk)-zqenwb(jl,jk)
            zmftop = -rmfdeps*pmfub(jl)
            if ( zbuo < d_zero .and. prfl(jl) > 10.0D0*zmftop*zcond(jl)) then
              kdtop(jl) = jk
              lddraf(jl) = .true.
              ptd(jl,jk) = zttest
              pqd(jl,jk) = zqtest
              pmfd(jl,jk) = zmftop
              pmfds(jl,jk) = pmfd(jl,jk)*(rcpd*ptd(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk) = pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1) = -d_half*pmfd(jl,jk)*zcond(jl)
              prfl(jl) = prfl(jl)+pdmfdp(jl,jk-1)
            end if
          end if
        end do
      end do
    end if
  end subroutine cudlfsn
!
!          THIS ROUTINE CALCULATES CUMULUS DOWNDRAFT DESCENT
!          M.TIEDTKE         E.C.M.W.F.    12/86 MODIF. 12/89
!          PURPOSE.
!          --------
!          TO PRODUCE THE VERTICAL PROFILES FOR CUMULUS DOWNDRAFTS
!          (I.E. T,Q,U AND V AND FLUXES)
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT IS T,Q,P,PHI,U,V AT HALF LEVELS.
!          IT RETURNS FLUXES OF S,Q AND EVAPORATION RATE
!          AND U,V AT LEVELS WHERE DOWNDRAFT OCCURS
!          METHOD.
!          --------
!          CALCULATE MOIST DESCENT FOR ENTRAINING/DETRAINING PLUME BY
!          A) MOVING AIR DRY-ADIABATICALLY TO NEXT LEVEL BELOW AND
!          B) CORRECTING FOR EVAPORATION TO OBTAIN SATURATED STATE.
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    INPUT PARAMETERS (LOGICAL):
!    *LDDRAF*       .TRUE. IF DOWNDRAFTS EXIST
!    INPUT PARAMETERS (REAL):
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)
!    UPDATED PARAMETERS (REAL):
!    *PRFL*         PRECIPITATION RATE                           KG/(M2*S)
!    OUTPUT PARAMETERS (REAL):
!    *PTD*          TEMPERATURE IN DOWNDRAFTS                      K
!    *PQD*          SPEC. HUMIDITY IN DOWNDRAFTS                 KG/KG
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                       KG/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS       J/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS         KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS     KG/(M2*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                   KG/(M2*S)
!    *PKINED*       DOWNDRAFT KINETIC ENERGY                     M2/S2
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO EVAPORATION IN
!          SATURATED DESCENT
!          REFERENCE
!          ---------
!          (TIEDTKE,1989)
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             03-08-28 : Clean-up detrainment rates   P. BECHTOLD
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!
!----------------------------------------------------------------------
!
  subroutine cuddrafn(kidia,kfdia,klon,ktdia,klev,lddraf,ptenh,pqenh, &
                      puen,pven,pgeo,pgeoh,paph,prfl,ptd,pqd,pmfu,    &
                      pmfd,pmfds,pmfdq,pdmfdp,pdmfde,pmfdde_rate,pkined)
    implicit none

    integer , intent(in) :: klon 
    integer , intent(in) :: klev 
    integer , intent(in) :: kidia 
    integer , intent(in) :: kfdia 
    integer :: ktdia ! argument not used
    logical , dimension(klon) , intent(in) :: lddraf
    real(dp) , dimension(klon,klev) , intent(in) :: ptenh , pqenh , &
               puen , pven , pgeo , pmfu
    real(dp) , dimension(klon,klev+1) , intent(in) :: pgeoh , paph
    real(dp) , dimension(klon) , intent(inout) :: prfl
    real(dp) , dimension(klon,klev) , intent(inout) :: ptd , pqd , &
               pmfd , pmfds , pmfdq
    real(dp) , dimension(klon,klev) , intent(out) :: pdmfdp , pdmfde , &
               pmfdde_rate , pkined
    real(dp) , dimension(klon) :: zdmfen , zdmfde , zcond , zoentr , zbuoy
    real(dp) , dimension(klon) :: zph
    logical , dimension(klon) :: llo2
    integer :: icall , ik , is , njkt3 , itopde , jk , jl
    real(dp) :: zbuo , zbuoyz , zbuoyv , zdmfdp , zdz , zentr , zmfdqk ,  &
                zmfdsk , zqdde , zqeen , zrain , zsdde , zseen , zzentr , &
                zdkbuo , zdken
    real(dp) , parameter :: zfacbuo = d_half/(d_one+d_half)
    real(dp) , parameter :: z_cwdrag = (3.0D0/8.0D0)*0.506D0/0.2D0

    njkt3 = kz-2
    itopde = njkt3
    !
    ! 1. Calculate moist descent for cumulus downdraft by
    !      (a) Calculating entrainment/detrainment rates, 
    !          including organized entrainment dependent on
    !          negative buoyancy and assuming
    !          linear decrease of massflux in pbl
    !      (b) Doing moist descent - evaporative cooling
    !          and moistening is calculated in *cuadjtq*
    !      (c) Checking for negative buoyancy and
    !          specifying final t,q,u,v and downward fluxes
    !
    do jl = kidia , kfdia
      zoentr(jl) = d_zero
      zbuoy(jl) = d_zero
      zdmfen(jl) = d_zero
      zdmfde(jl) = d_zero
      pdmfde(jl,:) = d_zero
      pmfdde_rate(jl,:) = d_zero
      pkined(jl,:) = d_zero
    end do

    do jk = 3 , klev
      is = 0
      do jl = kidia , kfdia
        zph(jl) = paph(jl,jk)
        llo2(jl) = lddraf(jl) .and. pmfd(jl,jk-1) < d_zero
        if ( llo2(jl) ) then
          is = is+1
        end if
      end do
      if ( is == 0 ) cycle
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          zentr = entrdd*pmfd(jl,jk-1)*(pgeoh(jl,jk-1)-pgeoh(jl,jk))*regrav
          zdmfen(jl) = zentr
          zdmfde(jl) = zentr
        end if
      end do
      if ( jk > itopde ) then
        do jl = kidia , kfdia
          if ( llo2(jl) ) then
            zdmfen(jl) = d_zero
            zdmfde(jl) = pmfd(jl,itopde)*(paph(jl,jk)-paph(jl,jk-1)) / &
                                         (paph(jl,klev+1)-paph(jl,itopde))
          end if
        end do
      end if
      if ( jk <= itopde ) then
        do jl = kidia , kfdia
          if ( llo2(jl) ) then
            zdz = -(pgeoh(jl,jk-1)-pgeoh(jl,jk))*regrav
            zzentr = zoentr(jl)*zdz*pmfd(jl,jk-1)
            zdmfen(jl) = zdmfen(jl)+zzentr
            zdmfen(jl) = max(zdmfen(jl),0.3D0*pmfd(jl,jk-1))
            zdmfen(jl) = max(zdmfen(jl),-0.75D0*pmfu(jl,jk) - &
                                 (pmfd(jl,jk-1)-zdmfde(jl)))  
            zdmfen(jl) = min(zdmfen(jl),d_zero)
          end if
          pdmfde(jl,jk) = zdmfen(jl)-zdmfde(jl)
        end do
      end if
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          pmfd(jl,jk) = pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
          zseen = (cpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
          zqeen = pqenh(jl,jk-1)*zdmfen(jl)
          zsdde = (cpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
          zqdde = pqd(jl,jk-1)*zdmfde(jl)
          zmfdsk = pmfds(jl,jk-1)+zseen-zsdde
          zmfdqk = pmfdq(jl,jk-1)+zqeen-zqdde
          pqd(jl,jk) = zmfdqk*(d_one/min(-rmfcmin,pmfd(jl,jk)))
          ptd(jl,jk) = (zmfdsk*(d_one/min(-rmfcmin, &
                                           pmfd(jl,jk)))-pgeoh(jl,jk))/cpd
          ptd(jl,jk) = min(400.0D0,ptd(jl,jk))
          ptd(jl,jk) = max(100.0D0,ptd(jl,jk))
          zcond(jl) = pqd(jl,jk)
        end if
      end do
      ik = jk
      icall = 2
      call cuadjtq(kidia,kfdia,klon,ktdia,klev,ik,zph,ptd,pqd,llo2,icall)
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          zcond(jl) = zcond(jl)-pqd(jl,jk)
          zbuo = ptd(jl,jk)*(d_one+retv*pqd(jl,jk)) - &
                             ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
          if ( prfl(jl) > d_zero .and. pmfu(jl,jk) > d_zero ) then
            zrain = prfl(jl)/pmfu(jl,jk)
            zbuo = zbuo-ptd(jl,jk)*zrain
          end if
          if ( zbuo >= d_zero .or. prfl(jl) <= (pmfd(jl,jk)*zcond(jl)) ) then
            pmfd(jl,jk) = d_zero
            zbuo = d_zero
          end if
          pmfds(jl,jk) = (cpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
          pmfdq(jl,jk) = pqd(jl,jk)*pmfd(jl,jk)
          zdmfdp = -pmfd(jl,jk)*zcond(jl)
          pdmfdp(jl,jk-1) = zdmfdp
          prfl(jl) = prfl(jl)+zdmfdp
          !
          ! Compute organized entrainment for use at next level
          !
          zbuoyz = zbuo/ptenh(jl,jk)
          zbuoyv = zbuoyz
          zbuoyz = min(zbuoyz,d_zero)
          zdz = -(pgeo(jl,jk-1)-pgeo(jl,jk))
          zbuoy(jl) = zbuoy(jl)+zbuoyz*zdz
          zoentr(jl) = egrav*zbuoyz*d_half/(d_one+zbuoy(jl))
          !
          ! Store downdraught detrainment rates
          !
          pmfdde_rate(jl,jk) = -zdmfde(jl)
          !
          ! Compute kinetic energy
          !
          zdkbuo = zdz*zbuoyv*zfacbuo
          if ( zdmfen(jl) < d_zero )then
            zdken = min(d_one,(d_one+z_cwdrag) * &
                               zdmfen(jl)/min(-rmfcmin,pmfd(jl,jk-1)))  
          else
            zdken = min(d_one,(d_one+z_cwdrag) * &
                               zdmfde(jl)/min(-rmfcmin,pmfd(jl,jk-1)))  
          end if
          pkined(jl,jk) = max(d_zero,(pkined(jl,jk-1) * &
                                     (d_one-zdken)+zdkbuo)/(d_one+zdken))
        end if
      end do
    end do
  end subroutine cuddrafn
!
!-----------------------------------------------------------------------------
!
  real(dp) function minj(x,y)
    implicit none
    real(dp) , intent(in) :: x , y
    minj = y - d_half*(dabs(x-y)-(x-y))
  end function minj
  real(dp) function maxj(x,y)
    implicit none
    real(dp) , intent(in) :: x , y
    maxj = y + d_half*(dabs(x-y)+(x-y))
  end function maxj
  real(dp) function foedelta(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedelta = max(d_zero,sign(d_one,ptare-tzero))
  end function foedelta
  real(dp) function foeew(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeew = c2es*exp((c3les*foedelta(ptare) + &
            c3ies*(d_one-foedelta(ptare)))*(ptare-tzero) / &
            (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare)))))
  end function foeew
  real(dp) function foede(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foede = (foedelta(ptare)*c5alvcp+(d_one-foedelta(ptare))*c5alscp) / &
       (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**d_two
  end function foede
  real(dp) function foedesu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedesu = (foedelta(ptare)*c5les+(d_one-foedelta(ptare))*c5ies) / &
         (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**d_two
  end function foedesu
  real(dp) function foelh(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelh = foedelta(ptare)*wlhv + (d_one-foedelta(ptare))*wlhs
  end function foelh
  real(dp) function foeldcp(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcp = foedelta(ptare)*wlhvocp + (d_one-foedelta(ptare))*wlhsocp
  end function foeldcp
  real(dp) function foealfa(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfa = min(d_one,((max(rtice,min(rtwat,ptare))-rtice) * &
                  rtwat_rtice_r)**d_two)
  end function foealfa
  real(dp) function foeewm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewm = c2es*(foealfa(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
          (d_one-foealfa(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewm
  real(dp) function foedem(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedem = foealfa(ptare)*c5alvcp*(d_one/(ptare-c4les)**d_two) + &
            (d_one-foealfa(ptare))*c5alscp*(d_one/(ptare-c4ies)**d_two)
  end function foedem
  real(dp) function foeldcpm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpm = foealfa(ptare)*wlhvocp+(d_one-foealfa(ptare))*wlhsocp
  end function foeldcpm
  real(dp) function foelhm(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhm = foealfa(ptare)*wlhv+(d_one-foealfa(ptare))*wlhs
  end function foelhm
  real(dp) function foetb(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foetb = foealfa(ptare)*c3les*(tzero-c4les)*(d_one/(ptare-c4les)**d_two)+ &
      (d_one-foealfa(ptare))*c3ies*(tzero-c4ies)*(d_one/(ptare-c4ies)**d_two)
  end function foetb
  real(dp) function foealfcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foealfcu = min(d_one, &
           ((max(rtice,min(rtwat,ptare))-rtice)*rtwat_rtice_r)**d_two)
  end function foealfcu
  real(dp) function foeewmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewmcu = c2es*(foealfcu(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
            (d_one-foealfcu(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewmcu
  real(dp) function foedemcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foedemcu = foealfcu(ptare)*c5alvcp*(d_one/(ptare-c4les)**d_two) + &
           (d_one-foealfcu(ptare))*c5alscp*(d_one/(ptare-c4ies)**d_two)
  end function foedemcu
  real(dp) function foeldcpmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeldcpmcu = foealfcu(ptare)*wlhvocp+(d_one-foealfcu(ptare))*wlhsocp
  end function foeldcpmcu
  real(dp) function foelhmcu(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foelhmcu = foealfcu(ptare)*wlhv+(d_one-foealfcu(ptare))*wlhs
  end function foelhmcu
  real(dp) function foeewmo(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeewmo = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeewmo
  real(dp) function foeeliq(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeeliq = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeeliq
  real(dp) function foeeice(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeeice = c2es*exp(c3ies*(ptare-tzero)/(ptare-c4ies))
  end function foeeice
  real(dp) function foeles_v(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeles_v = c3les*(ptare-tzero)/(ptare-c4les)
  end function foeles_v
  real(dp) function foeies_v(ptare)
    implicit none
    real(dp) , intent(in) :: ptare
    foeies_v = c3ies*(ptare-tzero)/(ptare-c4ies)
  end function foeies_v
  real(dp) function foeewm_v(ptare,exp1,exp2)
    implicit none
    real(dp) , intent(in) :: ptare , exp1 , exp2
    foeewm_v = c2es*(foealfa(ptare)*exp1+(d_one-foealfa(ptare))*exp2)
  end function foeewm_v
  real(dp) function foeewmcu_v(ptare,exp1,exp2)
    implicit none
    real(dp) , intent(in) :: ptare , exp1 , exp2
    foeewmcu_v = c2es*(foealfcu(ptare)*exp1+(d_one-foealfcu(ptare))*exp2)
  end function foeewmcu_v
  subroutine vdiv(z,x,y,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: z
    real(dp) , dimension(:) , intent(in) :: x , y
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      z(i) = x(i)/y(i)
    end do
  end subroutine vdiv
  subroutine vexp(y,x,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: y
    real(dp) , dimension(:) , intent(in) :: x
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      y(i) = dexp(x(i))
    end do
  end subroutine vexp
  subroutine vrec(y,x,n)
    implicit none
    real(dp) , dimension(:) , intent(out) :: y
    real(dp) , dimension(:) , intent(in) :: x
    integer , intent(in) :: n
    integer :: i
    do i = 1 , n
      y(i) = d_one/x(i)
    end do
  end subroutine vrec

end module mod_cu_tiedtke_38r2
